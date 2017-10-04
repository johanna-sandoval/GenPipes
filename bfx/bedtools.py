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

from bfx import ucsc

def graph(input_bam, output_bed_graph, output_wiggle, library_type="PAIRED_END"):

    if library_type == "PAIRED_END":
        if "forward" in output_bed_graph:
            samtools_options = "-F 256 -f 81 "
        elif "reverse" in output_bed_graph:
            samtools_options = "-F 256 -f 97 "
        else:
            raise Exception("Error: PAIRED_END library was provided but no strand orientation could be determined from " + output_bed_graph + "...")
    else:
        samtools_options = "-F 256"

    return concat_jobs([
        Job(
            [input_bam],
            [output_bed_graph, output_wiggle],
            [
                ['bedtools', 'module_samtools'],
                ['bedtools', 'module_bedtools']
            ],
            command="""\
nmblines=$(samtools view {samtools_options} {input_bam} | wc -l) && \\
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \\
genomeCoverageBed {other_options} -bg -split -scale $scalefactor \\
  -ibam {input_bam} \\
  -g {chromosome_size} \\
  > {output_bed_graph}""".format(
                samtools_options=samtools_options,
                input_bam=input_bam,
                chromosome_size=config.param('bedtools_graph', 'chromosome_size', type='filepath'),
                other_options=config.param('bedtools_graph', 'other_options', required=False),
                output_bed_graph=output_bed_graph
            )
        ),
        ucsc.bedGraphToBigWig(
            output_bed_graph,
            output_wiggle,
            False
        )
    ])

def intersect(input_bam, output_bam, readset):
    target_bed = config.param('bedtools_intersect', 'target_bed', required=False)

    if target_bed:
        if target_bed == "auto":
            if readset.beds:
                target_bed = os.path.abspath(readset.beds[0])
        else:
            # Add filepath validation
            target_bed = config.param('bedtools_intersect', 'target_bed', type='filepath')

    return Job(
        [input_bam],
        [output_bam],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools intersect \\
  -a {input_bam} \\
  -b {target_bed} \\
  {other_options} > {output_bam}""".format(
            input_bam=input_bam,
            target_bed=target_bed,
            other_options=config.param('bedtools_intersect', 'other_options'),
            output_bam=output_bam
        )
    )


def intersect_beds(bed1, bed2, output_bed, other_options=""):

    return Job(
        [bed1, bed2],
        [output_bed],
        [
            ['bedtools', 'module_bedtools']
        ],
        command="""\
bedtools intersect {other_options}\\
  -a {bed1} \\
  -b {bed2} \\
  > {output_bed}""".format(
            bed1 = bed1,
            bed2 = bed2,
            other_options = other_options,
            output_bed = output_bed
        )
    )
