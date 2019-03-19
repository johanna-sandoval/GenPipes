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


def cross_correlation(subSampled_TA, scores_output, plot_output, other, sampleName, threads = 1):
    """
    calculates cross correlation analysis using the run_spp.R sctipr from phantompeakqualtools
    """

    spp_cmd = """Rscript --max-ppsize=500000 $R_TOOLS/run_spp.R
                -c={subSampled_TA} 
                -p=${threads} {other} 
                -savp={plot_output} 
                -out={scores_output}""".format(
                subSampled_TA = subSampled_TA, 
                threads = threads, 
                other = other if other else "",
                plot_output = plot_output, 
                scores_output = scores_output
                )


    return Job (input_files = [subSampled_TA], 
    output_files = [scores_output, plot_output], 
    module_entries = [['metrics', 'module_mugqic_tools'], ['metrics', 'module_R']], 
    command = spp_cmd,
    name = "metrics.cross_correlation." + sampleName
    )


def spr_SE(prefix, tagAlign_input, pr1_output, pr2_output, sampleName):
    """
    Generates pseudo replicates from tagAlign files based on the ATAC ENCODE pipeline for Single Ended experiments
    """

    cmd_nlines = """export TA_lines=$(zcat -f {tagAlign_input} | wc -l) && let nlines=\"($TA_lines + 1)/2\"""".format(tagAlign_input = tagAlign_input)
    #int((get_num_lines(ta)+1)/2)


    cmd_tmp = """zcat {tagAlign_input} |
        shuf --random-source={tagAlign_input} |
        split -d -l $nlines - {prefix}""".format(tagAlign_input = tagAlign_input,
            prefix = prefix)


    cmd_zip1 = """gzip -nc {pr1_tmp}> {pr1_output} &&
            rm {pr1_tmp}""".format(
            pr1_tmp = prefix + "00",
            pr1_output = pr1_output)

    cmd_zip2 = """gzip -nc {pr2_tmp}> {pr2_output} &&
            rm {pr2_tmp}""".format(
            pr2_tmp = prefix + "01",
            pr2_output = pr2_output)

    cmd = """{cmd_nlines} && 
        {cmd_tmp} &&
        {cmd_zip1} && 
        {cmd_zip2}
        """.format(
        cmd_nlines = cmd_nlines,
        cmd_tmp = cmd_tmp,
        cmd_zip1 = cmd_zip1,
        cmd_zip2 = cmd_zip2)

    return Job (input_files = [tagAlign_input], 
    output_files = [pr1_output, pr2_output], 
    command = cmd,
    name = "atac.spr_SE." + sampleName
    )




def spr_PE(prefix, tagAlign_input, pr1_output, pr2_output, sampleName):
    """
    Generates pseudo replicates from tagAlign files based on the ATAC ENCODE pipeline for Single Ended experiments
    """

    cmd_nlines = """export TA_lines=$(zcat -f {tagAlign_input} | wc -l) && let nlines=\"($TA_lines/2 + 1)/2\"""".format(tagAlign_input = tagAlign_input)
    #int((get_num_lines(ta)/2+1)/2)


    cmd_tmp = """zcat -f {tagAlign_input} | sed \'N;s/\\n/\\t/\' |
        shuf --random-source={tagAlign_input} |
        split -d -l $nlines - {prefix}""".format(tagAlign_input = tagAlign_input,
            prefix = prefix)


    cmd_zip1 = """zcat -f {pr1_tmp} |
            awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' |
            gzip -nc > {pr1_output}
            rm {pr1_tmp}""".format(
            pr1_tmp = prefix + "00",
            pr1_output = pr1_output)


    cmd_zip2 = """zcat -f {pr2_tmp} |
            awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' |
            gzip -nc > {pr2_output}
            rm {pr2_tmp}""".format(
            pr2_tmp = prefix + "01",
            pr2_output = pr2_output)


    cmd = """{cmd_nlines} && 
        {cmd_tmp} &&
        {cmd_zip1} && 
        {cmd_zip2}""".format(
        cmd_nlines = cmd_nlines,
        cmd_tmp = cmd_tmp,
        cmd_zip1 = cmd_zip1,
        cmd_zip2 = cmd_zip2)

    return Job (input_files = [tagAlign_input], 
    output_files = [pr1_output, pr2_output], 
    command = cmd,
    name = "atac.spr_PE." + sampleName
    )



def tn5_shift_ta(tagAlign_input, output, sampleName):

    cmd = """zcat -f {tagAlign_input} |
        awk 'BEGIN {{OFS = "\\t"}}{{ if ($6 == "+") {{$2 = $2 + 4}} else if ($6 == "-") {{$3 = $3 - 5}} print $0}}' | 
        gzip -nc > {output}""".format(
        tagAlign_input = tagAlign_input,
        output = output)

    return Job (input_files = [tagAlign_input], 
    output_files = [output], 
    command = cmd,
    name = "atac.tn5_shift_ta." + sampleName
    )





