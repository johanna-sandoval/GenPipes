#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import logging
import os
import tempfile
import textwrap
from uuid import uuid4

from config import config
from utils import utils

# Output comment separator line
separator_line = "#" + "-" * 79

logger = logging.getLogger(__name__)

def create_scheduler(s_type, config_files, container=None):
    if s_type == "pbs":
        return PBSScheduler(config_files, container=container)
    elif s_type == "batch":
        return BatchScheduler(config_files)
    elif s_type == "daemon":
        return DaemonScheduler(config_files)
    elif s_type == "slurm":
        return SlurmScheduler(config_files, container=container)
    else:
        raise Exception("Error: scheduler type \"" + s_type + "\" is invalid!")

class Scheduler(object):
    def __init__(self, config_files, container=None, **kwargs):

        self.name = 'generic'
        self._config_files = config_files
        self._container = container
        self._host_cvmfs_cache = None
        self._cvmfs_cache = None
        self._bind = None

    def submit(self, pipeline):
        # Needs to be defined in scheduler child class
        raise NotImplementedError

    @property
    def container(self):
        return self._container

    @property
    def host_cvmfs_cache(self):

        if self._host_cvmfs_cache is None:

            self._host_cvmfs_cache = config.param("container", 'host_cvmfs_cache'
                                            , required=False, type="string")

            if not self._host_cvmfs_cache:

                tmp_dir = config.param("DEFAULT", 'tmp_dir', required=True)
                tmp_dir = utils.expandvars(tmp_dir)

                if not tmp_dir:
                    tmp_dir = None

                self._host_cvmfs_cache = tempfile.mkdtemp(prefix="genpipes_cvmfs_", dir=tmp_dir)

        return self._host_cvmfs_cache

    @property
    def cvmfs_cache(self):

        if self._cvmfs_cache is None:
            self._cvmfs_cache = config.param("container", 'cvmfs_cache', required=False, type="string")
            if not self._cvmfs_cache:
                self._cvmfs_cache = "/cvmfs-cache"

        return self._cvmfs_cache

    @property
    def bind(self):
        if self._bind is None:
            self._bind = config.param("container", 'bind_list', required=False, type='list')

            if not self._bind:
                bind = ['/tmp', '/home']
        return self._bind

    @property
    def disable_modulercfile(self):
        if self.container:
            return 'unset MODULERCFILE'
        return ""

    @property
    def container_line(self):
        if self.container:
            if self.container.type == 'docker':
                v_opt = ' -v {}:{}'.format(self.host_cvmfs_cache, self.cvmfs_cache)
                for b in self.bind:
                    v_opt += ' -v {0}:{0}'.format(b)
                network = " --network host"
                user = " --user $UID:$GROUPS"

                mount_point=''
                return ('docker run --env-file <( env| cut -f1 -d= ) --privileged '
                        '{network} {user} {v_opt} {name} '
                        .format(network=network, user=user, v_opt=v_opt
                                , name=self.container.name))

            elif self.container.type == 'singularity':
                b_opt = ' -B {}:{}'.format(self.host_cvmfs_cache, self.cvmfs_cache)
                for b in self.bind:
                    b_opt += ' -B {0}:{0}'.format(b)

                return ("singularity run {b_opt} {name}   "
                        .format(b_opt=b_opt, name=self.container.name))
        else:
            return ""


    def print_header(self, pipeline):
        print("""\
#!/bin/bash
# Exit immediately on error
{scheduler.disable_modulercfile}
set -eu -o pipefail

{separator_line}
# {pipeline.__class__.__name__} {scheduler.name} Job Submission Bash script
# Version: {pipeline.version}
# Created on: {pipeline.timestamp}
# Steps:
{steps}
{separator_line}"""
            .format(
                separator_line=separator_line,
                pipeline=pipeline,
                scheduler=self,
                steps="\n".join(["#   " + step.name + ": " + str(len(step.jobs)) + " job" + ("s" if len(step.jobs) > 1 else "" if step.jobs else "... skipping") for step in pipeline.step_range]) + \
                "\n#   TOTAL: " + str(len(pipeline.jobs)) + " job" + ("s" if len(pipeline.jobs) > 1 else "" if pipeline.jobs else "... skipping")
            )
        )

        if pipeline.jobs:
            print(
"""
OUTPUT_DIR={pipeline.output_dir}
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/{pipeline.__class__.__name__}_job_list_$TIMESTAMP
export CONFIG_FILES="{config_files}"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
"""
                .format(
                    pipeline=pipeline,
                    config_files=",".join([c.name for c in self._config_files ])
                )
            )

        if pipeline.json:
            json_files = []
            for step in pipeline.step_range:
                for job in step.jobs:
                    for sample in job.samples:
                        json_files.append(os.path.join(pipeline.output_dir, "json", sample.json_file))
            json_files = list(set(json_files))
            for j_file in json_files:
                print(
"""sed -i "s/\\"submission_date\\": \\"\\",/\\"submission_date\\": \\"$TIMESTAMP\\",/" {file}"""
                    .format(
                        file=j_file
                    )
                )

        ## Print a copy of sample JSONs for the genpipes dashboard
        if pipeline.json and pipeline.portal_output_dir != "":
            copy_commands = []
            test_copy_commands = []
            for i, sample in enumerate(pipeline.sample_list):
                unique_uuid = uuid4().get_hex()
                input_file = pipeline.sample_paths[i]
                output_file = os.path.join(pipeline.portal_output_dir, '$USER.' + sample.name + '.' + unique_uuid + '.json')
                #test_output_file = os.path.join("/lb/project/mugqic/analyste_dev/portal_test_dir/", '$USER.' + sample.name + '.' + unique_uuid + '.json')
                copy_commands.append("cp \"{input_file}\" \"{output_file}\"".format(
                    input_file=input_file, output_file=output_file))
                #test_copy_commands.append("cp \"{input_file}\" \"{output_file}\"".format(
                    #input_file=input_file, output_file=test_output_file))
            print(textwrap.dedent("""
                #------------------------------------------------------------------------------
                # Print a copy of sample JSONs for the genpipes dashboard
                #------------------------------------------------------------------------------
                {copy_commands}
            """).format(copy_commands='\n'.join(copy_commands)))
            #print(textwrap.dedent("""
                ##------------------------------------------------------------------------------
                ## Print a copy of sample JSONs for testing of the dashboard
                ##------------------------------------------------------------------------------
                #{copy_commands}
            #""").format(copy_commands='\n'.join(test_copy_commands)))

    def print_step(self, step):
        print("""
{separator_line}
# STEP: {step.name}
{separator_line}
STEP={step.name}
mkdir -p $JOB_OUTPUT_DIR/$STEP
""".format(separator_line=separator_line, step=step)
        )

    def job2json(self, pipeline, step, job, job_status):
        if not pipeline.json:
            return ""

        json_file_list = ",".join([os.path.join(pipeline.output_dir, "json", sample.json_file) for sample in job.samples])
        return """\
module load {module_python}
{job2json_script} \\
  -u \\"$USER\\" \\
  -c \\"{config_files}\\" \\
  -s \\"{step.name}\\" \\
  -j \\"$JOB_NAME\\" \\
  -d \\"$JOB_DONE\\" \\
  -l \\"$JOB_OUTPUT\\" \\
  -o \\"{jsonfiles}\\" \\
  -f {status}
module unload {module_python} {command_separator}""".format(
            job2json_script=os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "utils", "job2json.py"),
            module_python=config.param('DEFAULT', 'module_python'),
            step=step,
            jsonfiles=json_file_list,
            config_files=",".join([ os.path.relpath(c.name, pipeline.output_dir) for c in self._config_files ]),
            status=job_status,
            command_separator="&&" if (job_status=='\\"running\\"') else ""
        ) if json_file_list else ""



class PBSScheduler(Scheduler):

    def __init__(self, *args, **kwargs):
        super(PBSScheduler, self).__init__(*args, **kwargs)
        self.name = 'PBS/TORQUE'

    def submit(self, pipeline):
        self.print_header(pipeline)
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    if job.dependency_jobs:
                        # Chunk JOB_DEPENDENCIES on multiple lines to avoid lines too long
                        max_dependencies_per_line = 50
                        dependency_chunks = [job.dependency_jobs[i:i + max_dependencies_per_line] for i in range(0, len(job.dependency_jobs), max_dependencies_per_line)]
                        job_dependencies = "JOB_DEPENDENCIES=" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunks[0]])
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += "\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunk])
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    #sleepTime = random.randint(10, 100)
                    print("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << '{limit_string}'
{job.command_with_modules}
{limit_string}
)""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line,
                            limit_string=os.path.basename(job.done)
                        )
                    )

                    cmd = """\
echo "rm -f $JOB_DONE && {job2json_start} {container_line} $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
{job2json_end}
if [ \$MUGQIC_STATE -eq 0 ] ; then
  touch $JOB_DONE ;
fi
exit \$MUGQIC_STATE" | \\
""".format(
                        container_line=self.container_line,
                        job2json_start=self.job2json(pipeline, step, job, '\\"running\\"'),
                        job2json_end=self.job2json(pipeline, step, job, '\\$MUGQIC_STATE')
                    )
                        #sleep_time=sleepTime

                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                    job_name_prefix = job.name.split(".")[0]
                    cmd += \
                        config.param(job_name_prefix, 'cluster_submit_cmd') + " " + \
                        config.param(job_name_prefix, 'cluster_other_arg') + " " + \
                        config.param(job_name_prefix, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        config.param(job_name_prefix, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        config.param(job_name_prefix, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        config.param(job_name_prefix, 'cluster_walltime') + " " + \
                        config.param(job_name_prefix, 'cluster_queue') + " " + \
                        config.param(job_name_prefix, 'cluster_cpu')
                    #cmd += \
                        #config.param(job_name_prefix, 'cluster_submit_cmd') + " " + \
                        #config.param(job_name_prefix, 'cluster_other_arg') + " " + \
                        #config.param(job_name_prefix, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        #config.param(job_name_prefix, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        #config.param(job_name_prefix, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        #config.param(job_name_prefix, 'cluster_queue') + " -l walltime=1:00:0 -l nodes=1:ppn=1 "

                    if job.dependency_jobs:
                        cmd += " " + config.param(job_name_prefix, 'cluster_dependency_arg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(job_name_prefix, 'cluster_submit_cmd_suffix')

                    if config.param(job_name_prefix, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST\n"

                    print cmd

        # Check cluster maximum job submission
        cluster_max_jobs = config.param('DEFAULT', 'cluster_max_jobs', type='posint', required=False)
        if cluster_max_jobs and len(pipeline.jobs) > cluster_max_jobs:
            logging.warning("Number of jobs: " + str(len(pipeline.jobs)) + " > Cluster maximum number of jobs: " + str(
                cluster_max_jobs) + "!")


class BatchScheduler(Scheduler):

    def __init__(self, *args, **kwargs):
        super(BatchScheduler, self ).__init__(*args, **kwargs)
        self.name = 'Batch'

    def submit(self, pipeline):
        self.print_header(pipeline)
        if pipeline.jobs:
            print("SEPARATOR_LINE=`seq -s - 80 | sed 's/[0-9]//g'`")
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    print("""
{separator_line}
# JOB: {job.name}
{separator_line}
JOB_NAME={job.name}
JOB_DONE={job.done}
printf "\\n$SEPARATOR_LINE\\n"
echo "Begin MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`" && \\
rm -f $JOB_DONE && {job2json_start} \\
{job.command_with_modules}
MUGQIC_STATE=$PIPESTATUS
echo "End MUGQIC Job $JOB_NAME at `date +%FT%H:%M:%S`"
echo MUGQICexitStatus:$MUGQIC_STATE
{job2json_end}
if [ $MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; else exit $MUGQIC_STATE ; fi
""".format(
                            job=job,
                            separator_line=separator_line,
                            job2json_start=self.job2json(pipeline, step, job, '\\"running\\"'),
                            job2json_end=self.job2json(pipeline, step, job, '\\$MUGQIC_STATE')
                        )
                    )


class SlurmScheduler(Scheduler):

    def __init__(self, *args, **kwargs):
        super(SlurmScheduler, self).__init__(*args, **kwargs)
        self.name = 'SLURM'

    def submit(self, pipeline):
        self.print_header(pipeline)
        for step in pipeline.step_range:
            if step.jobs:
                self.print_step(step)
                for job in step.jobs:
                    if job.dependency_jobs:
                        # Chunk JOB_DEPENDENCIES on multiple lines to avoid lines too long
                        max_dependencies_per_line = 50
                        dependency_chunks = [job.dependency_jobs[i:i + max_dependencies_per_line] for i in range(0, len(job.dependency_jobs), max_dependencies_per_line)]
                        job_dependencies = "JOB_DEPENDENCIES=" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunks[0]])
                        for dependency_chunk in dependency_chunks[1:]:
                            job_dependencies += "\nJOB_DEPENDENCIES=$JOB_DEPENDENCIES:" + ":".join(["$" + dependency_job.id for dependency_job in dependency_chunk])
                    else:
                        job_dependencies = "JOB_DEPENDENCIES="

                    print("""
{separator_line}
# JOB: {job.id}: {job.name}
{separator_line}
JOB_NAME={job.name}
{job_dependencies}
JOB_DONE={job.done}
JOB_OUTPUT_RELATIVE_PATH=$STEP/${{JOB_NAME}}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$JOB_OUTPUT_DIR/$STEP/${{JOB_NAME}}_$TIMESTAMP.sh
cat << '{limit_string}' > $COMMAND
{job.command_with_modules}
{limit_string}
chmod 755 $COMMAND
""".format(
                            job=job,
                            job_dependencies=job_dependencies,
                            separator_line=separator_line,
                            limit_string=os.path.basename(job.done)
                    )
                        )

                    cmd = """\
echo "#! /bin/bash
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
rm -f $JOB_DONE && {job2json_start} {container_line}  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
{job2json_end}
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch
echo '#######################################'
exit \$MUGQIC_STATE" | \\
""".format(
                        job=job,
                        job2json_start=self.job2json(pipeline, step, job, '\\"running\\"'),
                        job2json_end=self.job2json(pipeline, step, job, '\\$MUGQIC_STATE') ,
                        container_line=self.container_line
)

                    # Cluster settings section must match job name prefix before first "."
                    # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                    job_name_prefix = job.name.split(".")[0]
                    cmd += \
                        config.param(job_name_prefix, 'cluster_submit_cmd') + " " + \
                        config.param(job_name_prefix, 'cluster_other_arg') + " " + \
                        config.param(job_name_prefix, 'cluster_work_dir_arg') + " $OUTPUT_DIR " + \
                        config.param(job_name_prefix, 'cluster_output_dir_arg') + " $JOB_OUTPUT " + \
                        config.param(job_name_prefix, 'cluster_job_name_arg') + " $JOB_NAME " + \
                        config.param(job_name_prefix, 'cluster_walltime') + " " + \
                        config.param(job_name_prefix, 'cluster_queue') + " " + \
                        config.param(job_name_prefix, 'cluster_cpu')
                    if job.dependency_jobs:
                        cmd += " " + config.param(job_name_prefix, 'cluster_dependency_arg') + "$JOB_DEPENDENCIES"
                    cmd += " " + config.param(job_name_prefix, 'cluster_submit_cmd_suffix')

                    if config.param(job_name_prefix, 'cluster_cmd_produces_job_id'):
                        cmd = job.id + "=$(" + cmd + ")"
                    else:
                        cmd += "\n" + job.id + "=" + job.name

                    # Write job parameters in job list file
                    cmd += "\necho \"$" + job.id + "\t$JOB_NAME\t$JOB_DEPENDENCIES\t$JOB_OUTPUT_RELATIVE_PATH\" >> $JOB_LIST\n"

                    #add 0.2s sleep to let slurm submiting the job correctly
                    cmd += "\nsleep 0.1\n"

                    print cmd

        # Check cluster maximum job submission
        cluster_max_jobs = config.param('DEFAULT', 'cluster_max_jobs', type='posint', required=False)
        if cluster_max_jobs and len(pipeline.jobs) > cluster_max_jobs:
            logger.warning("Number of jobs: " + str(len(pipeline.jobs)) + " > Cluster maximum number of jobs: " + str(
                cluster_max_jobs) + "!")


class DaemonScheduler(Scheduler):

    def __init__(self, *args, **kwargs):
        super(DaemonScheduler, self).__init__(*args, **kwargs)
        self.name = 'DAEMON'

    def submit(self, pipeline):
        print self.json(pipeline)

    def json(self, pipeline):
        #with open('sample.json', 'w') as json_file:
        #json.dump(the_dump, json_file)
        return json.dumps(
            {'pipeline': {
                'output_dir': pipeline.output_dir,
                'samples': [{
                    'name': sample.name,
                    'readsets': [{
                        "name": readset.name,
                        "library": readset.library,
                        "runType": readset.run_type,
                        "run": readset.run,
                        "lane": readset.lane,
                        "adapter1": readset.adapter1,
                        "adapter2": readset.adapter2,
                        "qualityoffset": readset.quality_offset,
                        "bed": [bed for bed in readset.beds],
                        "fastq1": readset.fastq1,
                        "fastq2": readset.fastq2,
                        "bam": readset.bam,
                     } for readset in pipeline.readsets if readset.sample.name == sample.name]
                } for sample in pipeline.samples],
                'steps': [{
                    'name': step.name,
                    'jobs': [{
                        "job_name": job.name,
                        "job_id": job.id,
                        "job_command": job.command_with_modules,
                        "job_input_files": job.input_files,
                        "job_output_files": job.output_files,
                        "job_dependencies": [dependency_job.id for dependency_job in job.dependency_jobs],
                        "job_cluster_options": {
                            # Cluster settings section must match job name prefix before first "."
                            # e.g. "[trimmomatic] cluster_cpu=..." for job name "trimmomatic.readset1"
                            'cluster_submit_cmd': config.param(job.name.split(".")[0], 'cluster_submit_cmd'),
                            'cluster_other_arg': config.param(job.name.split(".")[0], 'cluster_other_arg'),
                            'cluster_work_dir_arg': config.param(job.name.split(".")[0], 'cluster_work_dir_arg') + " " + pipeline.output_dir,
                            'cluster_output_dir_arg': config.param(job.name.split(".")[0], 'cluster_output_dir_arg') + " " + os.path.join(pipeline.output_dir, "job_output", step.name, job.name + ".o"),
                            'cluster_job_name_arg': config.param(job.name.split(".")[0], 'cluster_job_name_arg') + " " + job.name,
                            'cluster_walltime': config.param(job.name.split(".")[0], 'cluster_walltime'),
                            'cluster_queue': config.param(job.name.split(".")[0], 'cluster_queue'),
                            'cluster_cpu': config.param(job.name.split(".")[0], 'cluster_cpu')
                        },
                        "job_done": job.done
                    } for job in step.jobs]
                } for step in pipeline.step_range]
            }}, indent=4)
