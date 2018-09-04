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

def run(fastqs1, fastqs2, transcripts, output_dir):
    output_file = os.path.join(output_dir, "GMAP-fusion.fusion_predictions.tsv")
    return Job(
        [transcripts],
        [output_file],
        [
	        ['run_gmap_fusion', 'module_perl'],
	        ['run_gmap_fusion', 'module_bowtie2'],
            ['run_gmap_fusion', 'module_gmap'],
	        ['run_gmap_fusion', 'module_gmap_fusion'],
            ['run_gmap_fusion', 'module_samtools'],
	        ['run_gmap_fusion', 'module_python'],
        ],

        command="""\
$GMAPF_HOME/GMAP-fusion {options} \\
        --transcripts {transcripts} \\
        --genome_lib_dir {genome_build} \\
        --left_fq {fastq1} \\
        --right_fq {fastq2}
        --output {output_dir}""".format(
            options=config.param('run_gmap_fusion', 'options'),
	        genome_build=config.param('run_gmap_fusion', 'genome_build'),
	        transcripts=transcripts,
	        fastq1=",".join(fastq1 for fastq1 in fastqs1),
	        fastq2=",".join(fastq2 for fastq2 in fastqs2),
	        output_dir=output_dir,
        ),
    )
