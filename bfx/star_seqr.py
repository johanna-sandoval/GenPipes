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

def run(fastq1, fastq2, output_dir):
    return Job(
        [fastq1, fastq2],
        output_dir,
        [
            ['run_star_seqr', 'module_conda']
        ],

        command="""\
    starseqr.py -t {threads} {options} \\
      -i {genome_build} \\
      -g {gene_annot} \\
      -r {reference} \\
      -1 {fastq1} \\
      -2 {fastq2} \\
      -p {output_dir}""".format(
            genome_build=config.param('run_star_seqr', 'genome_build'),
            gene_annot=config.param('run_star_seqr', 'gene_annot'),
            reference=config.param('run_star_seqr', 'reference'),
            threads=config.param('run_star_seqr', 'threads', type='posint'),
            options=config.param('run_star_seqr', 'options'),
            fastq1=fastq1,
            fastq2=fastq2,
            output_dir=output_dir,
        ),
    )