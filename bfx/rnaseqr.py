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
            ['rnaseqr','module_rnaseqr'],
            ['rnaseqr','module_python']
        ],

    command="""\
    bam_stat.py \\
        -i {input} \\
        -q 10 \\
        > {output}""".format(
            input=input,
            output=output,
        ),
    )

def gene_body_coverage(input, output):
    return Job(
        [input],
        [output],
        [
            ['rnaseqr', 'module_rnaseqr'],
            ['rnaseqr', 'module_python']
        ],

        command="""\
        geneBody_coverage.py \\
            -r {geneBody} \\
            -i {input} \\
            -o {output}""".format(
            geneBody=config.param('runseqr', 'geneBody'),
            input=input,
            output=output,
        ),
    )

def infer_experiment(input, output):
    return Job(
        [input],
        [output],
        [
            ['rnaseqr', 'module_rnaseqr'],
            ['rnaseqr', 'module_python']
        ],

        command="""\
        infer_experiment.py \\
            -r {refseq} \\
            -i {input} \\
            > {output}""".format(
            refseq=config.param('runseqr', 'refseq'),
            input=input,
            output=output,
        ),
    )

def inner_distance(input, output):
    return Job(
        [input],
        [output],
        [
            ['rnaseqr', 'module_rnaseqr'],
            ['rnaseqr', 'module_python']
        ],

        command="""\
        inner_distance.py \\
            -r {refseq} \\
            -i {input} \\
            -o {output}""".format(
            refseq=config.param('runseqr', 'refseq'),
            input=input,
            output=output,
        ),
    )

def junction_annotation(input, output):
    return Job(
        [input],
        [output],
        [
            ['rnaseqr', 'module_rnaseqr'],
            ['rnaseqr', 'module_python']
        ],

        command="""\
        junction_annotation.py \\
            -r {refseq} \\
            -i {input} \\
            -o {output}""".format(
            refseq=config.param('runseqr', 'refseq'),
            input=input,
            output=output,
        ),
    )

def junction_saturation(input, output):
    return Job(
        [input],
        [output],
        [
            ['rnaseqr', 'module_rnaseqr'],
            ['rnaseqr', 'module_python']
        ],

        command="""\
        junction_saturation.py \\
            -r {refseq} \\
            -i {input} \\
            -o {output}""".format(
            refseq=config.param('runseqr', 'refseq'),
            input=input,
            output=output,
        ),
    )