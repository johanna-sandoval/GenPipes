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

def kraken2(input1, input2, prefix, other_options=config.param('kraken2', 'other_options', required=False), nthread=config.param('kraken2', 'threads', required=False), database=config.param('kraken2', 'database', required=False)):
    # unclassified_output = prefix + ".unclassified_sequences#.fastq"
    # classified_output = prefix + ".classified_sequences#.fastq"
    output = prefix + ".kraken2_output"
    report = prefix + ".kraken2_report"

    if input2:  # Paired end reads
        inputs = [input1, input2]
        unclassified_output = prefix + ".unclassified_sequences.fastq"
        classified_output = prefix + ".classified_sequences.fastq"
    else:   # Single end reads
        inputs = [input1]
        unclassified_output = prefix + ".unclassified_sequences#.fastq"
        classified_output = prefix + ".classified_sequences#.fastq"

    # inputs = [input]
    outputs = [
        unclassified_output,
        classified_output,
        output,
        report
        ]

    return Job(
        inputs,
        outputs,
        [
            ['kraken2', 'module_kraken2']
        ],

        command="""\
kraken2 \\
  {other_options} \\
  {nthread} \\
  {database} \\
  {paired} \\
  {unclassified_output} \\
  {classified_output} \\
  {output} \\
  {report} \\
  {inputs}""".format(
      other_options=other_options,
      nthread="--threads " + nthread,
      database="--db " + database,
      unclassified_output="--unclassified-out " + unclassified_output,
      classified_output="--classified-out " + classified_output,
      paired="--paired" if input2 else "",
      output="--output " + output,
      report="--report " + report,
      inputs=" \\\n  ".join(inputs)
      ),
    )
