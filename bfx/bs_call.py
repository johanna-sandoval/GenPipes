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

def call(input, sample_name, output):

    return Job(
        [input],
        [output],
        [
            ['bscall_call', 'module_bscall'],
        ],
        command="""\
bs_call -r {reference} \\
    -n {sample_name} \\
    {output} \\
    {input}""".format(
        reference=config.param('bscall_call', 'reference'),
        sample_name=sample_name,
        output=" -o ".join(output) if output else "",
        input=input
        # output_type=config.param('bscall_call', 'output_type')
        )
    )
