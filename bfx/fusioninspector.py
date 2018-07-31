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

def run(fusion_lists, fastqs1, fastqs2, sample_name, output_dir):
    output_file = os.path.join(output_dir, sample_name, sample_name + ".fusion_predictions.final.abridged.FFPM.annotated.coding_effect")
    return Job(
        [fusion_lists],
        [output_file],
        [
	        ['fusion_annotation', 'module_perl'],
	        ['fusion_annotation', 'module_python'],
	        ['fusion_annotation', 'module_htslib'],
	        ['fusion_annotation', 'module_gmap'],
	        ['fusion_annotation', 'module_trinity'],
	        ['fusion_annotation', 'module_star'],
            ['fusion_annotation', 'module_samtools'],
            ['fusion_annotation', 'module_star_fusion'],
        ],

        command="""\
$FUSIONINSPECTOR_HOME/FusionInspector \\
        --fusions {fusion_list} \\
        --genome_lib_dir {genome_build} \\
        --left_fq {fastq1} \\
        --right_fq {fastq2} \\
        --CPU {threads} \\
        --out_prefix {sample_name} \\
        --out_dir {output_dir} \\
        {options}""".format(
            fusion_list=fusion_lists,
	        genome_build=config.param('run_star_fusion', 'genome_build'),
            options=config.param('fusion_annotation', 'options'),
	        threads=config.param('fusion_annotation', 'threads'),
            fastq1=",".join(fastq1 for fastq1 in fastqs1),
            fastq2=",".join(fastq2 for fastq2 in fastqs2),
	        sample_name=sample_name,
            output_dir=os.path.join(output_dir, sample_name),
        ),
    )