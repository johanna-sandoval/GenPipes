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

<<<<<<< HEAD
def run(fastqs1, fastqs2, output_dir):
    output_file = os.path.join(output_dir, "fusion_predictions.abridged.coding_effect.tsv")
    jxt_file = os.path.join(output_dir, "Chimeric.out.junction")
    star_file = os.path.join(output_dir, "std.STAR.bam")
    return Job(
        fastqs1,
        [output_file,jxt_file,star_file],
=======
def run(fastq1, fastq2, output_dir):
    return Job(
        [fastq1, fastq2],
        output_dir,
>>>>>>> 7044fc6c (Updates and bug fixes)
        [
            ['run_star_fusion','module_perl'],
            ['run_star_fusion','module_star'],
            ['run_star_fusion','module_samtools'],
<<<<<<< HEAD
            ['run_star_fusion','module_star_fusion'],
            #['run_star_fusion', 'module_gcc']
=======
            ['run_star_fusion','module_star_fusion']
>>>>>>> 7044fc6c (Updates and bug fixes)
        ],

        command="""\
    STAR-Fusion --CPU {threads} {options} \\
        --genome_lib_dir {genome_build} \\
        --left_fq {fastq1} \\
        --right_fq {fastq2} \\
        --output_dir {output_dir}""".format(
            genome_build=config.param('run_star_fusion', 'genome_build'),
            threads=config.param('run_star_fusion', 'threads', type='posint'),
            options=config.param('run_star_fusion', 'options'),
<<<<<<< HEAD
            fastq1=",".join(fastq1 for fastq1 in fastqs1),
            fastq2=",".join(fastq2 for fastq2 in fastqs2),
=======
            fastq1=fastq1,
            fastq2=fastq2,
>>>>>>> 7044fc6c (Updates and bug fixes)
            output_dir=output_dir,
        ),
    )