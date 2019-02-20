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

# MUGQIC Modules
from core.config import *
from core.job import *


def align_atac(sample, fastq1, fastq2=None,  multimapping=0, idx=None, threads=1, ini_section='bowtie_align', log='align_samtools_sort.log'):
    
    other_options = config.param(ini_section, 'other_options', required=False)
    bam = sample + ".bam"

    SE_reads = "-U " + fastq1
    PE_reads = "--mm -X2000 -1 " + fastq1 + " -2 " + fastq2 if fastq2 else ""

    command_align = """bowtie2 {other_options} \\
    -k {multimapping} --local \\
    -x {bwt2_idx} \\
    --threads {nth_bwt2} \\
    {reads} \\
    | samtools view -bS - > {bam} \\
    2> {log}""".format(
            other_options=" \\\n  " + other_options if other_options else "",
            multimapping = multimapping, 
            bwt2_idx = idx if idx else config.param(ini_section, 'genome_bwt2_index'),
            nth_bwt2 = threads,
            reads = PE_reads if fastq2 else SE_reads,
            bam = bam,
            log = log
            )

    return Job(input_files = [fastq1, fastq2],
            output_files = [bam],
            module_entries = [['align_atac', 'module_bowtie2'], ['align_atac', 'module_samtools']],
            name = "align_atac." + sample,
            command = command_align
            )

