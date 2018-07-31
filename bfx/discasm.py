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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def run(fastqs1, fastqs2, chimeric_jxt, star_fusion_bam, output_dir):
	
    output_file = os.path.join(output_dir, "trinity_out_dir", "Trinity.fasta")
    return Job(
        [chimeric_jxt,star_fusion_bam],
        [output_file],
        [
	        ['run_discasm', 'module_discasm'],
	        ['run_discasm', 'module_perl'],
	        ['run_discasm', 'module_python'],
	        ['run_discasm', 'module_trinity'],
	        ['run_discasm', 'module_bowtie'],
            ['run_discasm', 'module_jellyfish'],
            ['run_discasm', 'module_salmon'],
        ],

        command="""\
DISCASM {options} \\
        --chimeric_junctions {chimeric_jxt} \\
        --aligned_bam {star_fusion_bam} \\
        --left_fq {fastq1} \\
        --right_fq {fastq2} \\
        --denovo_assembler {denovo_assembler} \\       
        --out_dir {output_dir}""".format(
            chimeric_jxt=chimeric_jxt,
	        star_fusion_bam=star_fusion_bam,
            options=config.param('run_discasm_gmap_fusion', 'trinity_options'),
	        denovo_assembler="Trinity",
            fastq1=fastqs1,
            fastq2=fastqs2,
            output_dir=output_dir,
        ),
    )