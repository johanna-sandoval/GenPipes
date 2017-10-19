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

def find_intersect_snps(input, output):
    return Job(
        [input],
        [output],
        [
            ['allelic_mapping','module_wasp'],
            ['allelic_mapping','module_python']
        ],

    command="""\
python $WASP_MAPPING/find_intersecting_snps.py \\
    --is_paired_end \\
    --is_sorted \\
    --snp_tab $WASP_DATA/phase1_v3.snp_tab.h5 \\
    --snp_index $WASP_DATA/phase1_v3.snp_index.h5 \\
    --haplotype $WASP_DATA/phase1_v3.haplotypes.h5 \\
    --samples $WASP_DATA/samples.txt \\
    --output_dir find_intersecting_snps \\
    {input}""".format(
        input=input,
        output=output,
    ),
)

def filter_remapped_reads(to_remap_input, remap_input, output):
    return Job(
        [to_remap_input, remap_input],
        [output],
        [
            ['allelic_mapping','module_wasp'],
            ['allelic_mapping','module_python']
        ],

    command="""\
python $WASP_MAPPING/filter_remapped_reads.py \\
    {to_remap_input} \\
    {remap_input} \\
    {output}""".format(
        to_remap_input=to_remap_input,
        remap_input=remap_input,
        output=output,
    ),
)
def rmdup_pe(input, output):
    return Job(
        [input],
        [output],
        [
            ['allelic_mapping','module_wasp'],
            ['allelic_mapping','module_python']
        ],

    command="""\
python $WASP_MAPPING/rmdup_pe.py \\
    {input} \\
    {output}""".format(
        input=input,
        output=output,
    ),
)
