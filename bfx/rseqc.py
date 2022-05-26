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

def bam_stat(input, output):
    return Job(
        [input],
        [output],
        [
            ['rseqc','module_rseqc']
        ],

    command="""\
    bam_stat.py {options} \\
        -i {input} \\
        > {output}""".format(
            options=config.param('rseqc', 'bam_stat_options'),
            input=input,
            output=output,
        ),
    )

def gene_body_coverage(input, output):
    return Job(
        [input],
        [output],
        [
            ['rseqc', 'module_rseqc']
        ],

        command="""\
        geneBody_coverage.py \\
            -r {geneBody} \\
            -i {input} \\
            -o {output}""".format(
            geneBody=config.param('rseqc', 'housekeeping'),
            input=input,
            output=output,
        ),
    )

def infer_experiment(input, output):
    return Job(
        [input],
        [output],
        [
            ['rseqc', 'module_rseqc']
        ],

        command="""\
        infer_experiment.py \\
            -r {ref_gene_model} \\
            -i {input} \\
            > {output}""".format(
            ref_gene_model=config.param('rseqc', 'ref_gene_model'),
            input=input,
            output=output,
        ),
    )

def inner_distance(input, output):
    return Job(
        [input],
        [output],
        [
            ['rseqc', 'module_rseqc']
        ],

        command="""\
        inner_distance.py \\
            -r {ref_gene_model} \\
            -i {input} \\
            -o {output}""".format(
            ref_gene_model=config.param('rseqc', 'ref_gene_model'),
            input=input,
            output=output,
        ),
    )

def junction_annotation(input, output):
    return Job(
        [input],
        [output],
        [
            ['rseqc', 'module_rseqc']
        ],

        command="""\
        junction_annotation.py \\
            -r {refseq} \\
            -i {input} \\
            -o {output}""".format(
            refseq=config.param('rseqc', 'refseq'),
            input=input,
            output=output,
        ),
    )

def junction_saturation(input, output):
    return Job(
        [input],
        [output],
        [
            ['rseqc', 'module_rseqc']
        ],

        command="""\
        junction_saturation.py \\
            -r {ref_gene_model} \\
            -i {input} \\
            -o {output}""".format(
            ref_gene_model=config.param('rseqc', 'ref_gene_model'),
            input=input,
            output=output,
        ),
    )

def tin(input):
    output = re.sub("\.bam$", ".summary.txt", os.path.basename(input))
    return Job(
        [input],
        [output],
        [
            ['rseqc', 'module_rseqc']
        ],

        command="""\
        tin.py \\
            -r {ref_gene_model} \\
            -i {input}""".format(
            ref_gene_model=config.param('rseqc', 'ref_gene_model'),
            input=input
        ),
    )