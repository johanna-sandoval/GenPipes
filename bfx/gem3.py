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

def align(input1, input2, output, report):

    if input2:  # Paired end reads
        inputs = [input1, input2]
        input_files = ["-1 " + input1, "-2 " +  input2]
        # outputs = [paired_output1, unpaired_output1, paired_output2, unpaired_output2]
    else:   # Single end reads
        inputs = [input1]
        input_files = ["-i " + input1]
        # outputs = [single_output]

    ## Get param from config file
    nbthreads = config.param('gem3_align', 'nbthreads', type='int')

    return Job(
        inputs,
        [output],
        [
            ['gem3_align', 'module_gem3']
        ],
        command="""\
gem-mapper \\
  -I {index} \\
  {inputs} \\
  {output_file} \\
  -t {nbthreads} \\
  --report-file={report}""".format(
      index=config.param('gem3_align', 'index'),
      inputs=" \\\n  ".join(input_files),
      output_file=" -o ".join(output) if output else "",
      nbthreads=nbthreads if str(nbthreads) != "" and isinstance(nbthreads, int) and nbthreads > 0 else 1,
      report=report
      )
    )
