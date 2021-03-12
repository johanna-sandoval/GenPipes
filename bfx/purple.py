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

def run(amber, cobalt, normal_name, tumor_name, output_dir, somatic_snv=None):

    return Job(
        [amber, cobalt],
        [output_dir],
        [
            ['purple', 'module_java'],
            ['purple', 'module_R'],
            ['purple', 'module_perl'],
            ['purple', 'module_circos'],
            ['purple', 'module_purple'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $PURPLE_JAR \\
  -ref_genome {reference_sequence} \\
  -reference {reference} \\
  -tumor {tumor} \\
  -cobalt {cobalt} \\
  -amber {amber} \\
  -gc_profile {gc_profile} \\
  -circos circos \\
  {output_dir}{somatic_snv}""".format(
        tmp_dir=config.param('purple', 'tmp_dir'),
        java_other_options=config.param('purple', 'java_other_options'),
        ram=config.param('purple', 'ram'),
        threads=config.param('purple', 'threads'),
        gc_profile=config.param('purple', 'gc_profile'),
        reference_sequence=config.param('purple', 'genome_fasta', type='filepath'),
        reference=normal_name,
        tumor=tumor_name,
        amber=amber,
        cobalt=cobalt,
        output_dir=output_dir,
        somatic_snv=somatic_snv,
        )
    )