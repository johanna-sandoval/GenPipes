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
import logging
import math
import os
import re
import sys


# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from pipelines import common
from pipelines.chipseq import chipseq


from bfx import picard
from bfx import samtools
from bfx import macs2
from bfx import bowtie2
from bfx import samtools


log = logging.getLogger(__name__)

class AtacSeq(chipseq.ChipSeq):
    """
    ATAC-Seq Pipeline
    ==============

    ATAC-Seq experiments allow researchers to understand chromatin accessibility using transposase Tn5 transposition. 
    Similar techniques include MNase-seq, FAIRE-Seq, and DNAse-seq.
    This pipeline analyzes both ATAC-Seq experimental data based on the [ENCODE ATAC-Seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).
    Tt starts by trimming adaptors and low quality bases.
    It then maps the reads to a reference genome using Bowtie2.
    

    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        super(AtacSeq, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {'alignment_output_directory': 'alignment',
                'trim_output_directory': 'trim',
                'report_output_directory': 'report',
                'metrics_output_directory': 'metrics',
                'homer_output_directory': 'tags',
                'graphs_output_directory': 'graphs',
                'tracks_output_directory': 'tracks',
                'macs_output_directory': 'peak_call',
                'anno_output_directory': 'annotation',
                'ihecA_output_directory': 'ihec_alignment',
                'ihecM_output_directory': 'ihec_metrics'
                }
        return dirs

    def bowtie_align(self):
        """
        Single ended or paired-ended reads are aligned to the genome using bowtie2. For PE mode, -X2000 is used to allow fragments up to 2Kb to align to capture open chromatin regions.
        The reads are then sorted by coordinates using samtools.

        For more detailed information about the Bowtie2 process visit: [Bowtie] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
        """

        jobs = []

        ## ENCODE standard multimapping=4
        multimapping = config.param('bowtie_align', 'multimapping')
        threads = config.param('bowtie_align', 'threads')
        index = config.param('bowtie_align', 'genome_bwt2_index')


        for readset in self.readsets:
            sample_output_dir = os.path.join(self.output_dirs['alignment_output_directory'], readset.sample.name, readset.name)
            trim_file_prefix = os.path.join(self.output_dirs['trim_output_directory'], readset.sample.name, readset.name + ".trim.")
            prefix = os.path.join(sample_output_dir, readset.name)
            
            log = os.path.join(sample_output_dir, readset.name + ".bowtie_align.log")

             # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            job_align = concat_jobs([
                Job(command="mkdir -p " + sample_output_dir),
                bowtie2.align_atac(prefix, fastq1, fastq2, multimapping, index, threads, 'bowtie_align', log)
            ])
            job_align.name = "bowtie_align.align_atac." + readset.name
            job_align.samples = [readset.sample]


            jobs.append(job_align)

        return jobs    



    def samtools_merge_bams(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [samtools](http://samtools.sourceforge.net/).

        This step takes as input files the aligned bams/sams from the hicup_align step
        """

        jobs = []
        for sample in self.samples:
            sample_output = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, sample.name + ".merged.bam")
            readset_bams = [os.path.join(self.output_dirs['alignment_output_directory'], readset.sample.name, readset.name, readset.name + ".bam") for readset in sample.readsets]

            mkdir_job = Job(command="mkdir -p " + self.output_dirs['alignment_output_directory'])

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, os.path.dirname(sample_output))

                job = concat_jobs([
                    mkdir_job,
                    Job(
                        readset_bams,
                        [sample_output],
                        command="ln -s -f " + target_readset_bam + " " + sample_output
                    )
                ], name="symlink_readset_sample_bam." + sample.name)

            elif len(sample.readsets) > 1:

                samtools_merge = samtools.merge(sample_output, readset_bams)

                job = concat_jobs([
                    mkdir_job,
                    samtools_merge
                ])
                job.name = "samtools_merge_bams." + sample.name

            job.samples = [sample]
            jobs.append(job)

        return jobs


    def samtools_sort_index(self):
        jobs = []
        for sample in self.samples:
            sample_prefix = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, sample.name)
            sample_bam = sample_prefix + ".merged.bam"
            sorted_readset_bam = sample_prefix + ".merged.sorted"

            job_rename = Job(command = "mv " + sorted_readset_bam + ".bam " +  sample_bam)

            jobs.append(concat_jobs([
                samtools.sort(sample_bam, sorted_readset_bam),
                job_rename,
                samtools.index(sample_bam)
                ], 
                name = "samtools_sort_index." + sample.name ,
                samples = [sample.name]))

        return jobs




    def samtools_view_filter(self):
        jobs = []
        for sample in self.samples:
            sample_prefix = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, sample.name)
            input_bam = sample_prefix +  ".sorted.dup.bam"
            filtered_bam = sample_prefix + ".filtered.bam"
            ## remove chrM too: chrM, chM, chrMT, chMT, M, MT, NC_012920.1
            if self.run_type == "PAIRED_END":
                options = "-b -F 1804 -f 2 -q " + str(config.param('samtools_view_filter', 'min_mapq', type='int'))
            else:
                options = "-b -F 1804 -q " + str(config.param('samtools_view_filter', 'min_mapq', type='int'))
            job = samtools.view(input_bam, filtered_bam, options)
            job.name = "samtools_view_filter." + sample.name
            job.samples = [sample.name]
            jobs.append(job)
        
        return jobs


    @property
    def steps(self):
        return [
            self.samtools_bam_sort,
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.bowtie_align,
            self.samtools_merge_bams,
            self.samtools_sort_index,
            self.picard_mark_duplicates,
            self.samtools_view_filter
        ]

if __name__ == '__main__':
    AtacSeq()
