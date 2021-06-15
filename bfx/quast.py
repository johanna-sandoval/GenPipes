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

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def quast(input, output_dir, prefix):
    inputs = [input]
    outputs = [
        os.path.join(output_dir, "report.html"),
        os.path.join(output_dir, "report.pdf"),
        os.path.join(output_dir, "report.tex"),
        os.path.join(output_dir, "report.tsv"),
        os.path.join(output_dir, "report.txt")
        ]

    return Job(
        inputs,
        outputs,
        [
            ['quast', 'module_quast']
        ],

        command="""\
quast.py {reference} \\
  {output_dir} \\
  {features} \\
  {nthread} \\
  {input}""".format(
      reference="-r " + config.param(prefix, 'reference_genome', required=False),
      features="--features " + config.param(prefix, 'genomic_feature', required=False) if config.param(prefix, 'genomic_feature', required=False) else "",
      nthread="--threads " + config.param(prefix, 'threads', required=False),
      output_dir="--output-dir " + output_dir,
      input=input
      ),
    )