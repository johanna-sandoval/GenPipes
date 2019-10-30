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

# Python Standard Modules
import os
import re

# MUGQIC Modules
from core.config import config
from core.job import Job

def index(
    input,
    output,
    basecalls_dir,
    mismatches,
    lane,
    mask
    ):

    return Job(
        [input],
        [output],
        [
            ["index", "module_java"],
            ["index", "module_mugqic_tools"]
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} \\
  {java_other_options} \\
  -Xmx{ram} \\
  -jar $JAVA_TOOLS/{jar} \\
  MAX_MISMATCHES={mismatches} \\
  NUM_PROCESSORS={threads} \\
  BARCODE_FILE={barcode_file} \\
  BASECALLS_DIR={basecalls_dir} \\
  LANE={lane_number} \\
  READ_STRUCTURE={read_structure} \\
  METRICS_FILE={output} \\
  TMP_DIR={tmp_dir}""".format(
            tmp_dir=config.param('index', 'tmp_dir'),
            java_other_options=config.param('index', 'java_other_options'),
            ram=config.param('index', 'ram'),
            jar=config.param('index', 'jar'),
            mismatches=mismatches,
            threads=config.param('index', 'threads'),
            barcode_file=config.param('index', 'barcode_file'),
            basecalls_dir=basecalls_dir,
            lane_number=lane,
            read_structure=mask,
            output=output
        )
    )

def bcl2fastq(
    input,
    fastq_outputs,
    output_dir,
    sample_sheet,
    run,
    lane,
    extra_option
    demultiplex=False,
    mismatches=None,
    mask=None,
    ):

    if demultiplex:
        demultiplex_parameters = """\
  --barcode-mismatches {number_of_mismatches} \\
  --use-bases-mask {mask}""".format(
            number_of_mismatches=mismatches,
            mask=mask
        )
    else:
        command_suffix = ""

    return Job(
        [input],
        fastq_outputs,
        [
            ['fastq', 'module_bcl_to_fastq']
        ],
        command="""\
bcl2fastq \\
  --runfolder-dir {run_dir} \\
  --output-dir {output_dir} \\
  --tiles {tiles} \\
  --sample-sheet {sample_sheet} \\
  --create-fastq-for-index-reads \\
  {demultiplex_parameters} \\
  {other_options} \\
  {extra_option}""".format(
            run_dir=run,
            output_dir=output_dir,
            tiles="s_" + str(lane),
            sample_sheet=sample_sheet,
            demultiplex_parameters=demultiplex_parameters,
            other_options=config.param('fastq', 'other_options'),
            extra_option=extra_option
        )
    )

def aggregate_fastqs(
    readset,
    merge,
    mask
    ):

    # For HaloPlex-like masks, do the necessary name swapping
    if ''.join(i for i in mask if not i.isdigit()) == "Y,I,Y,Y":
        rename R2 => I2
        rename R3 => R2
        index2=8

    read_inputs = []
    index_inputs = []
    read_outputs = []
    index_outputs = []
    read1_command = ""
    read2_command = ""

    # If 10X libraries : 4 indexes per sample
    if re.search("tenX", readset.library_type):
        read_inputs.append()
        index_inputs.append()
        read_outputs.append()
        index_outputs.append()
        cat 4 files > 1 file (R1)
        cat 4 files > 1 file (I1)

        # For paired-end sequencing, do not forget the fastq of the reverse reads
        if readset.run_type == "PAIRED_END" :
            read_inputs.append()
            read_outputs.append()
            cat 4 files > 1 file (R2)

        # For dual index multiplexing, do not forget the fastq of the second index
        if index2 != 0 :
            index_inputs.append()
            index_outputs.append()
            cat 4 files > 1 file (I2)

        # If True, then merge the 'Undeternined' reads
        if merge:
            read_inputs.append()
            index_inputs.append()
            cat 5 files > 1 file (R1)
            cat 5 files > 1 file (I1)

            # For paired-end sequencing, do not forget the fastq of the reverse reads
            if readset.run_type == "PAIRED_END" :
                read_inputs.append()
                read_outputs.append()
                cat 5 files > 1 file (R2)

            # For dual index multiplexing, do not forget the fastq of the second index
            if index2 != 0 :
                index_inputs.append()
                index_outputs.append()
                cat 5 files > 1 file (I2)

    # not a 10X library : 1 index per sample
    else:
        # If ask to merge the Undeternined reads
        if merge:
           
            cat 2 files > 1 file (R1)
            cat 2 files > 1 file (I1)
            # For paired-end sequencing, do not forget the fastq of the reverse reads
            if readset.run_type == "PAIRED_END" :
                cat 2 files > 1 file (R2)
            # For dual index multiplexing, do not forget the fastq of the second index
            if index2 != 0 :
                cat 2 files > 1 file (I2)

    inputs = read_inputs + index_inputs
    ouputs = read_outputs + index_outputs

    return Job(
        inputs,
        ouputs,
        [[]],
        command="""\
""".format(
        )
    )