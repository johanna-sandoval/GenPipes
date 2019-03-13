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

# MUGQIC Modules
from core.job import *



def tagAlign_SE(input_bam, tagAlign_output, sampleName):
    """
    makes tagAlign file using bedools from Single ended bam
    """
    cmd = 'bedtools bamtobed -i {bam} | '.format(bam = input_bam)
    cmd += 'awk \'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}\' | '
    cmd += 'gzip -nc > {TA}'.format(TA = tagAlign_output)

    return Job (input_files = [input_bam], 
    output_files = [tagAlign_output], 
    module_entries = [['tagAlign', 'module_bedtools']], 
    command = cmd,
    name = "tagAlign.SE." + sampleName
    )


def tagAlign_PE(nSortedBam, bedpeName, tagAlign_output, sampleName):
    """
    makes tagAlign file using bedools from  Pair ended bam
    """
    
    cmd_bedpe = 'LC_COLLATE=C bedtools bamtobed -bedpe -mate1 -i {nSortedBam} | '.format(nSortedBam = nSortedBam)
    cmd_bedpe += 'gzip -nc > {bedpe}'.format(bedpe = bedpeName)

    job_bedpe = Job (input_files = [nSortedBam], 
    output_files = [bedpeName], 
    module_entries = [['tagAlign', 'module_bedtools']], 
    command = cmd_bedpe,
    name = "tagAlign.PE.bedpe." + sampleName
    )

    cmd_TA = 'zcat -f {bedpe} | '.format(bedpe = bedpeName)
    cmd_TA += 'awk \'BEGIN{{OFS="\\t"}}'
    cmd_TA += '{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n'
    cmd_TA += '%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",'
    cmd_TA += '$1,$2,$3,$9,$4,$5,$6,$10}}\' | '
    cmd_TA += 'gzip -nc > {TA}'.format(TA = tagAlign_output)

    job_TA = Job (input_files = [bedpeName], 
    output_files = [tagAlign_output], 
    module_entries = [['tagAlign', 'module_bedtools']], 
    command = cmd_TA,
    name = "tagAlign.PE.tagAlign." + sampleName
    )

    return concat_jobs([job_bedpe, job_TA])


def tagAlign_subsample(tagAlign_input, subSampled_output, num_reads, type, sampleName):
    """
    subsamples a specific number of reads from tagAlign file for QC metrics
    if type = SE, input is the tagAlign file. If type = PE, input is the bedpe file.
    """
    if type == "SE":

        cmd = """zcat {tagAlign_input} |
        shuf -n {num_reads} --random-source={tagAlign_input} |
        gzip -nc > {subSampled_output}
        """.format(tagAlign_input = tagAlign_input,
            num_reads = num_reads,
            subSampled_output = subSampled_output)

    elif type == "PE":

        cmd = """zcat {tagAlign_input} |
        shuf -n {num_reads} --random-source={tagAlign_input}  |
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,"N","1000",$9}}' |
        gzip -nc > {subSampled_output}""".format(tagAlign_input = tagAlign_input,
            num_reads = num_reads,
            subSampled_output = subSampled_output)

    return Job (input_files = [tagAlign_input], 
    output_files = [subSampled_output], 
    command = cmd,
    name = "tagAlign.subsample." + sampleName
    )






