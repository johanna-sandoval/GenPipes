#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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
import argparse
import collections
import csv
import logging
import math
import os
import re
import sys
import subprocess

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
from pipelines import common
from bfx.sequence_dictionary import parse_sequence_dictionary_file, split_by_size
import utils.utils

#Expression tools
from bfx import bedtools
from bfx import bwa
from bfx import cufflinks
from bfx import stringtie
from bfx import ballgown
from bfx import differential_expression
from bfx import gq_seq_utils
from bfx import htseq
from bfx import metrics
from bfx import fastqc
from bfx import samtools
from bfx import bvatools
from bfx import rmarkdown
from bfx import tools
from bfx import ucsc
from bfx import star

#Variant tools
from bfx import picard2
from bfx import gatk4
from bfx import sambamba
from bfx import adapters
from bfx import skewer
from bfx import htslib
from bfx import vt
from bfx import gemini
from bfx import snpeff
from bfx import vcfanno
from bfx import star_fusion
from bfx import arriba
from bfx import annoFuse
from bfx import rseqc
from bfx import rnaseqc2

from bfx import bash_cmd as bash

log = logging.getLogger(__name__)

class RnaSeqRaw(common.Illumina):
    """
    RNA-Seq Pipeline
    ================

    The standard MUGQIC RNA-Seq pipeline has three protocols (stringtie, variants, cancer), stringtie is the
    default protocol and applicable in most cases.

    All three protocols are based on the use of the [STAR aligner](https://code.google.com/p/rna-star/)
    to align reads to the reference genome. These alignments are used during
    downstream analysis (for example in stringtie protocol, to determine differential expression of genes and transcripts).

    The [StringTie](https://ccb.jhu.edu/software/stringtie/) suite is used for differential transcript expression (DTE) analysis, whereas
    [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and
    [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) are used for the differential gene expression (DGE) analysis.

    The "stringtie" protocol is requireds a design file which will be used to define the comparison groups
    in the differential analyses. The design file format is described
    [here](https://genpipes.readthedocs.io/en/latest/get-started/concepts/design_file.html). In addition,
    [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential
    transcript and gene expression levels and test them for significant differences.

    The variants protocol is used when variant detection, is the main goal of the analysis. GATK best practices workflow
    is used for variant calling. It also contains a step for annotating genes using [gemini](https://gemini.readthedocs.io/en/latest/)

    The cancer protocol contains all the steps in the variants protocol but it is
    specifically designed for cancer data sets due to the
    complexity of cancer samples and additional analyses those projects often entail.
    The goal of the cancer protocol is comparing expression to known benchmarks. In addition to the steps in the variants
    protocol, it contains four specific steps. Three of them (run_star_fusion, run_arriba, run_annofuse)
    are related to detection and annotation of gene fusions. For that, [Star-fusion](https://github.com/STAR-Fusion/STAR-Fusion-Tutorial/wiki),
    [Arriba](https://arriba.readthedocs.io/en/latest/) and [anno-Fuse](https://rdrr.io/github/d3b-center/annoFuse/) are
    used. Furthermore, [RSeQC](http://rseqc.sourceforge.net/) provides RNA-seq quality control metrics to asses the
    quality of the data.

    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=argparse.FileType('r'))
        super(RnaSeqRaw, self).__init__(protocol)

    @property
    def sequence_dictionary(self):
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', param_type='filepath'), variant=False)
        return self._sequence_dictionary

    def build_adapter_file(self, directory, readset):
        adapter_file = os.path.join(directory, "adapter.tsv")
        if readset.run_type == "SINGLE_END":
            if readset.adapter1:
                adapter_job = Job(
                    command="""\
`cat > {adapter_file} << END
>Adapter\n{sequence}\n
END
`""".format(
                        adapter_file=adapter_file,
                        sequence=readset.adapter1
                    )
                )
            else:
                raise Exception("Error: missing adapter1 for SINGLE_END readset \"" + readset.name + "\", or missing adapter_file parameter in config file!")

        elif readset.run_type == "PAIRED_END":
            if readset.adapter1 and readset.adapter2:
                adapter_job = Job(
                    command="""\
`cat > {adapter_file} << END
>Adapter1\n{sequence1}\n
>Adapter2\n{sequence2}\n
END
`""".format(
                        adapter_file=adapter_file,
                        sequence1=readset.adapter1,
                        sequence2=readset.adapter2,
                    )
                )

        return adapter_job


    def skewer_trimming(self):
        """
        [Skewer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182) is used mainly for
        detection and trimming adapter sequences from raw fastq files. Other features of Skewer is listed
        [here](https://github.com/relipmoc/skewer).
        """

        jobs = []

        for readset in self.readsets:
            output_dir = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(output_dir, readset.name)

            adapter_file = config.param('skewer_trimming', 'adapter_file', required=False, param_type='filepath')
            adapter_job = None

            if not adapter_file:
                adapter_file = os.path.join(output_dir, "adapter.tsv")
                adapter_job = adapters.create(
                    readset,
                    adapter_file
                )

            fastq1 = ""
            fastq2 = ""
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dir,
                        "raw_reads",
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                    prefix = os.path.join(
                        self.output_dir,
                        "raw_reads",
                        readset.sample.name,
                        readset.name + "."
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dir,
                        "raw_reads",
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + ".single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        output_dir,
                        remove=True
                    ),
                    adapter_job,
                    skewer.trim(
                        fastq1,
                        fastq2,
                        trim_file_prefix,
                        adapter_file
                        ),
                    bash.ln(
                        trim_file_prefix + "-trimmed-pair1.fastq.gz",
                        trim_file_prefix + ".trim.pair1.fastq.gz",
                        self.output_dir
                        ),
                    bash.ln(
                        trim_file_prefix + "-trimmed-pair2.fastq.gz",
                        trim_file_prefix + ".trim.pair2.fastq.gz",
                        self.output_dir
                        )
                    ],
                    name="skewer_trimming." + readset.name,
                    removable_files=[output_dir],
                    samples=[readset.sample]
                    )
                )

        return jobs

    def star_genome_length(self):
        """
        Calculation for setting genomeSAindexNbases for STAR index
        """
        genome_index = csv.reader(open(config.param('DEFAULT', 'genome_fasta', param_type='filepath') + ".fai", 'r'),
                                  delimiter='\t')
    
        return int(min(14, math.log2(sum([int(chromosome[1]) for chromosome in genome_index])) / 2 - 1))

    def star(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per readset of sequencing
        using the [STAR](https://code.google.com/p/rna-star/) software. It generates a Binary Alignment Map file (.bam).

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        project_index_directory = "reference.Merged"
        project_junction_file = os.path.join("alignment_1stPass", "AllSamples.SJ.out.tab")
        individual_junction_list = []
        genome_length = self.star_genome_length()
        ######
        # pass 1 -alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_1stPass_directory = os.path.join("alignment_1stPass", readset.sample.name, readset.name)
            individual_junction_list.append(os.path.join(alignment_1stPass_directory, "SJ.out.tab"))
    
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                  re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))
    
            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)
    
            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_1stPass_directory,
                genome_index_folder=None,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else ""
            )
            job.name = "star_align.1." + readset.name
            job.samples = [readset.sample]
            jobs.append(job)

        ######
        jobs.append(concat_jobs([
            # pass 1 - contatenate junction
            Job(samples=self.samples),
            star.concatenate_junction(
                input_junction_files_list=individual_junction_list,
                output_junction_file=project_junction_file
            ),
            # pass 1 - genome indexing
            star.index(
                genome_index_folder=project_index_directory,
                junction_file=project_junction_file,
                genome_length=genome_length
            )], name="star_index.AllSamples", samples=self.samples))

        ######
        # Pass 2 - alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_2ndPass_directory = os.path.join("alignment", readset.sample.name, readset.name)
    
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                  re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))
    
            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)
    
            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_2ndPass_directory,
                genome_index_folder=project_index_directory,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else "",
                create_wiggle_track=True,
                search_chimeres=True,
                cuff_follow=True,
                sort_bam=True
            )
            job.samples = [readset.sample]
            job.input_files.append(os.path.join(project_index_directory, "SAindex"))
    
            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            # remove older symlink before otherwise it raise an error if the link already exist (in case of redo)
            if len(readset.sample.readsets) == 1:
                readset_bam = os.path.join(alignment_2ndPass_directory, "Aligned.sortedByCoord.out.bam")
                sample_bam = os.path.join("alignment", readset.sample.name, readset.sample.name + ".sorted.bam")
                job = concat_jobs([
                    job,
                    Job(
                        [readset_bam],
                        [sample_bam],
                        command="ln -s -f " + os.path.relpath(readset_bam,os.path.dirname(sample_bam))
                                + " " + sample_bam,
                        removable_files=[sample_bam]
                    )
                ])
    
            job.name = "star_align.2." + readset.name
            jobs.append(job)

        report_file = os.path.join("report", "RnaSeq.star.md")
        jobs.append(
            Job(
                [os.path.join("alignment", readset.sample.name, readset.name, "Aligned.sortedByCoord.out.bam") for
                 readset in self.readsets],
                [report_file],
                [['star', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable scientific_name="{scientific_name}" \\
  --variable assembly="{assembly}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    scientific_name=config.param('star', 'scientific_name'),
                    assembly=config.param('star', 'assembly'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="star_report",
                samples=self.samples
            )
        )

        return jobs

    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name, "Aligned.sortedByCoord.out.bam")
                          for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                job = gatk4.merge_sam_files(
                    inputs,
                    output
                )
                job.name = "picard_merge_sam_files." + sample.name
                job.samples = [sample]
                jobs.append(job)

        return jobs

    def picard_sort_sam(self):
        """
        The alignment file is reordered (QueryName) using [Picard](http://broadinstitute.github.io/picard/). The QueryName-sorted bam files will be used to determine raw read counts.
        """
    
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)
        
            job = gatk4.sort_sam(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".QueryNameSorted.bam",
                "queryname"
            )
            job.name = "picard_sort_sam." + sample.name
            job.samples = [sample]
            jobs.append(job)
        return jobs

    def picard_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        Will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "sorted.bam"
            output = alignment_file_prefix + "sorted.mdup.bam"

            metrics_file = os.path.join("metrics",
                                        "rna",
                                        sample.name,
                                        sample.name + ".sorted.mdup.metrics")
            jobs.append(concat_jobs([
                bash.mkdir(
                    os.path.join("metrics",
                                 "rna",
                                 sample.name),
                    remove=True
                ),
                gatk4.picard_mark_duplicates(
                    [input],
                    output,
                    metrics_file
                ),
            ],name = "picard_mark_duplicates." + sample.name)
            )

        return jobs

    def metrics_rna_fastqc(self):

        jobs = []
        for sample in self.samples:
            fastqc_directory = os.path.join("metrics", "rna", sample.name, "fastqc")
            input = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam")
            output_dir = os.path.join(fastqc_directory)
            output = os.path.join(fastqc_directory, sample.name + ".sorted.mdup_fastqc.zip")

            adapter_file = config.param('fastqc', 'adapter_file', required=False, param_type='filepath')
            adapter_job = None

            if not adapter_file:
                adapter_job = self.build_adapter_file(output_dir, sample.readsets[0])
                adapter_file = os.path.join(output_dir, "adapter.tsv")

            jobs.append(concat_jobs([
                bash.mkdir(
                    output_dir,
                    remove=True
                ),
                adapter_job,
                fastqc.fastqc(
                    input,
                    None,
                    output_dir,
                    output,
                    adapter_file
                ),
            ], name="fastqc." + sample.name))

        return jobs

    def picard_rna_metrics(self):
        """
        Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions
        metrics are collected using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        reference_file = config.param('picard_rna_metrics', 'genome_fasta', param_type='filepath')
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            [input] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]
            ])
            
            output_directory = os.path.join("metrics", "rna", sample.name)

            job = concat_jobs([
                bash.mkdir(
                    output_directory,
                    remove=True
                ),
                gatk4.collect_multiple_metrics(
                    input,
                    os.path.join(output_directory, sample.name),
                    reference_file
                ),
                gatk4.collect_rna_metrics(
                        input,
                        os.path.join(output_directory, sample.name + ".picard_rna_metrics")
                )
            ], name="picard_rna_metrics." + sample.name)

            jobs.append(job)

        return jobs

    def rseqc(self):
        """
        Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions
        metrics are collected using [Picard](http://broadinstitute.github.io/picard/).
        """
    
        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            [input] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]
            ])
        
            output_directory = os.path.join("metrics", "rna", sample.name)

            mkdir_job = bash.mkdir(
                output_directory,
                remove=True
            )

            jobs.append(
                concat_jobs([
                    mkdir_job,
                    bash.chgdir(
                        output_directory
                    ),
                    rseqc.tin(
                        input,
                        output_directory
                    )
                ],
                    name="rseqc.tin." + sample.name,
                    samples=[sample]
                )
            )
            
        return jobs
    
    def bam_hard_clip(self):
        jobs = []

        for sample in self.samples:
            alignment_input = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam")
            alignment_output = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.hardClip.bam")
            job=pipe_jobs([
                samtools.view(
                    alignment_input,
                    None,
                    "-h"
                ),
                Job(
                    [None],
                    [alignment_output],
                    # awk to transform soft clip into hard clip for tuxedo suite
                    command="""\
                    awk 'BEGIN {{OFS="\\t"}} {{if (substr($1,1,1)=="@") {{print;next}}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {{$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}}; if (C[length(C)]=="S") {{L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }}; gsub(/[0-9]*S/,"",$6); print}}' """.format()
                ),
                samtools.view(
                    "-",
                    alignment_output,
                    "-hbS"
                ),
            ])
            job.name="tuxedo_hard_clip."+ sample.name
            job.samples = [sample]
            jobs.append(job)
        
        return jobs

    def gatk_callable_loci(self):
        """
        Computes the callable region or the genome as a bed track.
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            alignment_directory = os.path.join("alignment", sample.name)
            [input] = self.select_input_files([
                [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                [os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]
            ])

            job = gatk4.callable_loci(
                input,
                alignment_file_prefix + "callable.bed",
                alignment_file_prefix + "callable.summary.txt"
            )
            job.name = "gatk_callable_loci." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def rnaseqc(self):
        """
        Computes a series of quality control metrics using [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc).
        """

        jobs = []
        sample_file = os.path.join("alignment", "rnaseqc.samples.txt")
        sample_rows = [
            [sample.name, os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam"), "RNAseq"] for
            sample in self.samples]
        input_bams = [sample_row[1] for sample_row in sample_rows]
        output_directory = os.path.join("metrics", "rnaseqRep")
        # Use GTF with transcript_id only otherwise RNASeQC fails
        gtf_transcript_id = config.param('rnaseqc', 'gtf_transcript_id', param_type='filepath')

        jobs.append(concat_jobs([
            Job(
                command="mkdir -p " + output_directory,
                removable_files=[output_directory],
                samples=self.samples
            ),
            Job(
                input_bams,
                [sample_file],
                command="""\
echo "Sample\tBamFile\tNote
{sample_rows}" \\
  > {sample_file}""".format(
                    sample_rows="\n".join(["\t".join(sample_row) for sample_row in sample_rows]),
                    sample_file=sample_file)
            ),
            metrics.rnaseqc(
                sample_file,
                output_directory,
                self.run_type == "SINGLE_END",
                gtf_file=gtf_transcript_id,
                reference=config.param(
                    'rnaseqc',
                    'genome_fasta',
                    param_type='filepath'
                ),
                ribosomal_interval_file=config.param(
                    'rnaseqc',
                    'ribosomal_fasta',
                    param_type='filepath')
            ),
            Job([],
                [output_directory + ".zip"],
                command="zip -r {output_directory}.zip {output_directory}".format(
                output_directory=output_directory)
                )
        ], name="rnaseqc")
        )

        trim_metrics_file = os.path.join("metrics", "trimSampleTable.tsv")
        metrics_file = os.path.join("metrics", "rnaseqRep", "metrics.tsv")
        report_metrics_file = os.path.join("report", "trimAlignmentTable.tsv")
        report_file = os.path.join("report", "RnaSeq.rnaseqc.md")
        jobs.append(
            Job(
                [metrics_file],
                [report_file, report_metrics_file],
                [['rnaseqc', 'module_python'], ['rnaseqc', 'module_pandoc']],
                # Ugly awk to merge sample metrics with trim metrics if they exist; knitr may do this better
                command="""\
mkdir -p report && \\
cp {output_directory}.zip report/reportRNAseqQC.zip && \\
python -c 'import csv; csv_in = csv.DictReader(open("{metrics_file}"), delimiter="\t")
print "\t".join(["Sample", "Aligned Reads", "Alternative Alignments", "%", "rRNA Reads", "Coverage", "Exonic Rate", "Genes"])
print "\\n".join(["\t".join([
    line["Sample"],
    line["Mapped"],
    line["Alternative Aligments"],
    str(float(line["Alternative Aligments"]) / float(line["Mapped"]) * 100),
    line["rRNA"],
    line["Mean Per Base Cov."],
    line["Exonic Rate"],
    line["Genes Detected"]
]) for line in csv_in])' \\
  > {report_metrics_file}.tmp && \\
if [[ -f {trim_metrics_file} ]]
then
  awk -F"\t" 'FNR==NR{{raw_reads[$1]=$2; surviving_reads[$1]=$3; surviving_pct[$1]=$4; next}}{{OFS="\t"; if ($2=="Aligned Reads"){{surviving_pct[$1]="%"; aligned_pct="%"; rrna_pct="%"}} else {{aligned_pct=($2 / surviving_reads[$1] * 100); rrna_pct=($5 / surviving_reads[$1] * 100)}}; printf $1"\t"raw_reads[$1]"\t"surviving_reads[$1]"\t"surviving_pct[$1]"\t"$2"\t"aligned_pct"\t"$3"\t"$4"\t"$5"\t"rrna_pct; for (i = 6; i<= NF; i++) {{printf "\t"$i}}; print ""}}' \\
  {trim_metrics_file} \\
  {report_metrics_file}.tmp \\
  > {report_metrics_file}
else
  cp {report_metrics_file}.tmp {report_metrics_file}
fi && \\
rm {report_metrics_file}.tmp && \\
trim_alignment_table_md=`if [[ -f {trim_metrics_file} ]] ; then cut -f1-13 {report_metrics_file} | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%.1f", $4), sprintf("%\\47d", $5), sprintf("%.1f", $6), sprintf("%\\47d", $7), sprintf("%.1f", $8), sprintf("%\\47d", $9), sprintf("%.1f", $10), sprintf("%.2f", $11), sprintf("%.2f", $12), sprintf("%\\47d", $13)}}}}' ; else cat {report_metrics_file} | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%.1f", $4), sprintf("%\\47d", $5), sprintf("%.2f", $6), sprintf("%.2f", $7), $8}}}}' ; fi`
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable trim_alignment_table="$trim_alignment_table_md" \\
  --to markdown \\
  > {report_file}""".format(
                    output_directory=output_directory,
                    report_template_dir=self.report_template_dir,
                    trim_metrics_file=trim_metrics_file,
                    metrics_file=metrics_file,
                    basename_report_file=os.path.basename(report_file),
                    report_metrics_file=report_metrics_file,
                    report_file=report_file
                ),
                report_files=[report_file],
                name="rnaseqc_report",
                samples=self.samples
            )
        )

        return jobs

    def rnaseqc2(self):
            """
            Computes a series of quality control metrics using [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc).
            """
            jobs = []

            # Use GTF with transcript_id only otherwise RNASeQC fails
            for sample in self.samples:
                alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)
                [input] = self.select_input_files([
                    [os.path.join(alignment_file_prefix + ".sorted.mdup.bam")]
                ])

                output_directory = os.path.join("metrics", "rna", sample.name, "rnaseqc2")

                jobs.append(concat_jobs([
                    Job(
                        command="mkdir -p " + output_directory,
                        removable_files=[output_directory],
                        samples=self.samples
                        ),
                    rnaseqc2.run(
                        input,
                        output_directory,
                    ),
            ], name="rnaseqc2." + sample.name))
    
            return jobs

    def split_N_trim(self):
        """
        SplitNtrim. A [GATK](https://software.broadinstitute.org/gatk/) tool called SplitNCigarReads developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions.
        """
    
        jobs = []
    
        nb_jobs = config.param('gatk_split_N_trim', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning(
                "Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
    
        for sample in self.samples:
            alignment_dir = os.path.join("alignment", sample.name)
            split_dir = os.path.join(alignment_dir, "splitNtrim")
            split_file_prefix = os.path.join(split_dir, sample.name + ".")
            alignment_file_prefix = os.path.join(alignment_dir, sample.name + ".")
        
            if nb_jobs == 1:
                job = concat_jobs([
                    bash.mkdir(
                        split_dir,
                        remove=True
                    ),
                    gatk4.split_n_cigar_reads(
                        alignment_file_prefix + "sorted.mdup.bam",
                        split_file_prefix + "sorted.mdup.split.bam"
                    ),
                    Job(
                        [split_file_prefix + "sorted.mdup.split.bam"],
                        [alignment_file_prefix + "sorted.mdup.split.bam"],
                        command="ln -s -f " +
                                os.path.relpath(split_file_prefix + "sorted.mdup.split.bam",
                                                os.path.dirname(
                                                    alignment_file_prefix + "sorted.mdup.split.bam"
                                                )
                                                ) + " " + alignment_file_prefix + "sorted.mdup.split.bam"
                    )
                ], name="gatk_split_N_trim." + sample.name)
                job.samples = [sample]
                jobs.append(job)
        
        
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)
            
                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    job = concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            split_dir,
                            remove=True
                        ),
                        gatk4.split_n_cigar_reads(
                            alignment_file_prefix + "sorted.mdup.bam",
                            split_file_prefix + "sorted.mdup.split." + str(idx) + ".bam",
                            intervals=sequences
                        )
                    ], name="gatk_split_N_trim." + sample.name + "." + str(idx))
                    job.samples = [sample]
                    jobs.append(job)
            
                job = concat_jobs([
                    bash.mkdir(
                        split_dir,
                        remove=True
                    ),
                    gatk4.split_n_cigar_reads(
                        alignment_file_prefix + "sorted.mdup.bam",
                        split_file_prefix + "sorted.mdup.split.others.bam",
                        exclude_intervals=unique_sequences_per_job_others
                    )
                ], name="gatk_split_N_trim." + sample.name + ".others")
                job.samples = [sample]
                jobs.append(job)
    
        return jobs

    def sambamba_merge_splitNtrim_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Sambamba] (http://lomereiter.github.io/sambamba/docs/sambamba-merge.html).
        """
    
        jobs = []
    
        nb_jobs = config.param('gatk_split_N_trim', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning(
                "Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
    
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            split_file_prefix = os.path.join(alignment_directory, "splitNtrim", sample.name + ".")
            output = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.bam")
        
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)
            
                inputs = []
                for idx, sequences in enumerate(unique_sequences_per_job):
                    inputs.append(
                        os.path.join(
                            split_file_prefix + "sorted.mdup.split." + str(idx) + ".bam"
                        )
                    )
                inputs.append(
                    os.path.join(
                        split_file_prefix + "sorted.mdup.split.others.bam"
                    )
                )
            
                job = sambamba.merge(inputs, output)
                job.name = "sambamba_merge_splitNtrim_files." + sample.name
                job.samples = [sample]
                jobs.append(job)
    
        return jobs

    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.
        """
        
        jobs = []
    
        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
    
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            input = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.bam")
            # input = os.path.join(alignment_directory, sample.name + ".keep.merged.dup.sorted.bam")
        
            if nb_jobs == 1:
                realign_prefix = os.path.join(realign_directory, "all")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                sample_output_bam = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")
                job = concat_jobs([
                    bash.mkdir(
                        realign_directory,
                        remove=True
                    ),
                    gatk4.realigner_target_creator(
                        input,
                        realign_intervals,
                        output_dir=self.output_dir,
                    ),
                    gatk4.indel_realigner(
                        input,
                        output=output_bam,
                        target_intervals=realign_intervals,
                        output_dir=self.output_dir,
                    ),
                    # Create sample realign symlink since no merging is required
                    Job(
                        [output_bam],
                        [sample_output_bam],
                        command="ln -s -f " +
                                os.path.relpath(output_bam,
                                                os.path.dirname(
                                                    sample_output_bam
                                                )
                                                ) + " " + sample_output_bam
                    )
                ], name="gatk_indel_realigner." + sample.name
                )
                job.samples = [sample]
                jobs.append(job)
        
            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)
            
                # Create one separate job for each of the first sequences
                for idx, sequence in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(realign_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequence
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    output_bam = realign_prefix + ".bam"
                    job = concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            realign_directory,
                            remove=True
                        ),
                        gatk4.realigner_target_creator(
                            input,
                            realign_intervals,
                            intervals=intervals,
                            output_dir=self.output_dir,
                        ),
                        gatk4.indel_realigner(
                            input,
                            output=output_bam,
                            target_intervals=realign_intervals,
                            intervals=intervals,
                            output_dir=self.output_dir,
                        )
                    ], name="gatk_indel_realigner." + sample.name + "." + str(idx))
                    job.samples = [sample]
                    jobs.append(job)
            
                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                job = concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    bash.mkdir(
                        realign_directory,
                        remove=True
                    ),
                    gatk4.realigner_target_creator(
                        input,
                        realign_intervals,
                        exclude_intervals=unique_sequences_per_job_others,
                        output_dir=self.output_dir,
                    ),
                    gatk4.indel_realigner(
                        input,
                        output=output_bam,
                        target_intervals=realign_intervals,
                        exclude_intervals=unique_sequences_per_job_others,
                        output_dir=self.output_dir,
                    )
                ], name="gatk_indel_realigner." + sample.name + ".others")
                job.samples = [sample]
                jobs.append(job)
    
        return jobs

    def sambamba_merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """
    
        jobs = []
    
        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', param_type='posint')
    
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")
        
            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)
            
                inputs = []
                for idx, sequences in enumerate(unique_sequences_per_job):
                    inputs.append(os.path.join(realign_directory, str(idx) + ".bam"))
                inputs.append(os.path.join(realign_directory, "others.bam"))
            
                job = sambamba.merge(inputs, merged_realigned_bam)
                job.name = "sambamba_merge_realigned." + sample.name
                job.removable_files = [merged_realigned_bam]
                job.samples = [sample]
                jobs.append(job)
    
        return jobs

    def recalibration(self):
        """
        Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration,
        the quality scores in the QUAL field in each read in the output BAM are more accurate in that
        the reported quality score is closer to its actual probability of mismatching the reference genome.
        Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle
        and sequence context, and by doing so, provides not only more accurate quality scores but also
        more widely dispersed ones.
        """
    
        jobs = []
    
        created_interval_lists = []
    
        for sample in self.samples:
            file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.split.realigned.")
            input = file_prefix + "bam"
            print_reads_output = file_prefix + "recal.bam"
            base_recalibrator_output = file_prefix + "recalibration_report.grp"
        
            interval_list = None
        
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)
            
                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(None, coverage_bed, interval_list)
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    job.samples = [sample]
                    jobs.append(job)
                    # created_interval_lists.append(interval_list)
        
            jobs.append(concat_jobs([
                gatk4.base_recalibrator(
                    input,
                    base_recalibrator_output,
                    intervals=interval_list
                ),
            ], name="gatk_base_recalibrator." + sample.name))
        
            jobs.append(concat_jobs([
                gatk4.print_reads(
                    input,
                    print_reads_output,
                    base_recalibrator_output
                ),
            ], name="gatk_print_reads." + sample.name))
    
        return jobs

    def gatk_haplotype_caller(self):
        """
        GATK haplotype caller for snps and small indels.
        """
    
        jobs = []
    
        nb_haplotype_jobs = config.param('gatk_haplotype_caller', 'nb_jobs', param_type='posint')
        if nb_haplotype_jobs > 50:
            log.warning(
                "Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
    
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
            input = os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")
            input = self.select_input_files(
                [[os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.recal.bam")],
                 [os.path.join(alignment_directory, sample.name + ".sorted.mdup.split.realigned.bam")],
                 [os.path.join(alignment_directory, sample.name + ".sorted.bam")]])
        
            if nb_haplotype_jobs == 1:
                job = concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(
                        command="mkdir -p " + haplotype_directory,
                        removable_files=[haplotype_directory]
                    ),
                    gatk4.haplotype_caller(
                        input,
                        os.path.join(haplotype_directory, sample.name + ".hc.vcf.gz")
                    )
                ], name="gatk_haplotype_caller." + sample.name)
                job.samples = [sample]
                jobs.append(job)
        
        
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(
                    self.sequence_dictionary,
                    nb_haplotype_jobs - 1,
                    variant=True
                )
            
                job = []
            
                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    job = concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(
                            command="mkdir -p " + haplotype_directory,
                            removable_files=[haplotype_directory]
                        ),
                        gatk4.haplotype_caller(
                            input, os.path.join(
                                haplotype_directory,
                                sample.name + "." + str(idx) + ".hc.vcf.gz"
                            ),
                            intervals=sequences
                        )
                    ], name="gatk_haplotype_caller." + sample.name + "." + str(idx))
                    job.samples = [sample]
                    jobs.append(job)
            
                # Create one last job to process the last remaining sequences and 'others' sequences
                job = concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(
                        command="mkdir -p " + haplotype_directory,
                        removable_files=[haplotype_directory]
                    ),
                    gatk4.haplotype_caller(
                        input,
                        os.path.join(haplotype_directory,
                                     sample.name + ".others.hc.vcf.gz"
                                     ),
                        exclude_intervals=unique_sequences_per_job_others
                    )
                ], name="gatk_haplotype_caller." + sample.name + ".others")
                job.samples = [sample]
                jobs.append(job)
    
        return jobs

    def merge_hc_vcf(self):
        """
        Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.
        """
    
        jobs = []
        nb_haplotype_jobs = config.param('gatk_haplotype_caller', 'nb_jobs', param_type='posint')
    
        for sample in self.samples:
            haplotype_file_prefix = os.path.join("alignment", sample.name, "rawHaplotypeCaller", sample.name)
            output_haplotype_file_prefix = os.path.join("alignment", sample.name, sample.name)
            if nb_haplotype_jobs == 1:
                vcfs_to_merge = [haplotype_file_prefix + ".hc.vcf.gz"]
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_haplotype_jobs - 1)
            
                vcfs_to_merge = [haplotype_file_prefix + "." + str(idx) + ".hc.vcf.gz" for idx in
                                 range(len(unique_sequences_per_job))]
                vcfs_to_merge.append(haplotype_file_prefix + ".others.hc.vcf.gz")
        
            job = concat_jobs([
                gatk4.cat_variants(
                    vcfs_to_merge,
                    output_haplotype_file_prefix + ".hc.vcf.gz"
                ),
            ], name="gatk_merge_vcfs." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def run_vcfanno(self):
        """
        vcfanno is used to annotate VCF files with preferred INFO fields from anu number of VCF or BED files. For more
        information [visit](https://github.com/brentp/vcfanno)
        """
        
        jobs = []
    
        for sample in self.samples:
            input = os.path.join("alignment", sample.name, sample.name + ".hc.vcf.gz")
            output = os.path.join("alignment", sample.name, sample.name + ".hc.rnaedit.vcf.gz")
            
            job = pipe_jobs([
                vcfanno.run(
                    input,
                    None,
                ),
                htslib.bgzip_tabix(
                    None,
                    output
                )
            ], name="run_vcfanno.rnaedit." + sample.name)
            
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def variant_filtration(self):
        """
        GATK VariantFiltration.
        VariantFiltration is a GATK tool for hard-filtering variant calls based on certain criteria. Records are hard-filtered
        by changing the value in the FILTER field to something other than PASS.
        """
    
        jobs = []
    
        for sample in self.samples:
            input_file_prefix = os.path.join("alignment", sample.name, sample.name)
        
            job = concat_jobs([
                gatk4.variant_filtration(
                    input_file_prefix + ".hc.rnaedit.vcf.gz",
                    input_file_prefix + ".hc.flt.vcf.gz",
                    config.param('gatk_variant_filtration', 'other_options')
                )
            ], name="gatk_variant_filtration." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def decompose_and_normalize(self):
        """
        [vt](https://genome.sph.umich.edu/wiki/Vt#Normalization) is used to normalized and decompose VCF files. For more
        information about normalizing and decomposing visit
        [here](https://research-help.genomicsengland.co.uk/display/GERE/Variant+Normalisation). An indexed file is also
        generated from the output file using [htslib](http://www.htslib.org/download/).
        """
    
        jobs = []
    
        for sample in self.samples:
            input_file_prefix = os.path.join("alignment", sample.name, sample.name)
        
            job = pipe_jobs([
                vt.decompose_and_normalize_mnps(
                    input_file_prefix + ".hc.flt.vcf.gz",
                    None
                ),
                htslib.bgzip_tabix(
                    None,
                    input_file_prefix + ".hc.flt.vt.vcf.gz"
                ),
            ], name="decompose_and_normalize." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def compute_snp_effects(self):
        """
        [SnpEff](https://pcingola.github.io/SnpEff/) is used to variant annotation and effect prediction on genes by
        using an interval forest approach. It annotates and predicts the effects of genetic variants
        (such as amino acid changes).
        """
        jobs = []
    
        for sample in self.samples:
            input_file_prefix = os.path.join("alignment", sample.name, sample.name)
        
            job = concat_jobs([
                snpeff.compute_effects(
                    input_file_prefix + ".hc.flt.vt.vcf.gz",
                    input_file_prefix + ".hc.snpeff.vcf",
                    options=config.param('compute_effects', 'options')
                ),
                htslib.bgzip_tabix(
                    input_file_prefix + ".hc.snpeff.vcf",
                    input_file_prefix + ".hc.snpeff.vcf.gz"
                ),
            ], name="compute_effects." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def gemini_annotations(self):
        """
        [gemini](https://github.com/arq5x/gemini) (GEnome MINIng) is used to integrative exploration of genetic
        variation and genome annotations. For more information
        [visit](https://gemini.readthedocs.io/en/latest/).
        """
    
        jobs = []
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])
    
        for sample in self.samples:
            temp_dir = os.path.join("alignment")
            input_file_prefix = os.path.join(temp_dir, sample.name, sample.name)
        
            job = concat_jobs([
                gemini.gemini_annotations(
                    input_file_prefix + ".hc.snpeff.vcf.gz",
                    input_file_prefix + ".gemini." + gemini_version + ".db",
                    temp_dir
                )
            ], name="gemini_annotations." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def run_star_fusion(self):
        """
        STAR-Fusion is a component of the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). Based on the STAR
        aligner it identifies candidate fusion transcripts supported by Illumina reads.
        https://github.com/STAR-Fusion/STAR-Fusion/wiki
        """
    
        jobs = []
    
        left_fastqs = collections.defaultdict(list)
        right_fastqs = collections.defaultdict(list)
    
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + "-trimmed-")
        
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                  re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
        
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
        
            left_fastqs[readset.sample.name].append(fastq1)
            right_fastqs[readset.sample.name].append(fastq2)
    
        for sample in self.samples:
            output_dir = os.path.join("fusion", sample.name, "star_fusion")
        
            mkdir_job = Job(command="mkdir -p " + output_dir)
        
            job = concat_jobs([
                mkdir_job,
                star_fusion.run(
                    left_fastqs[sample.name],
                    right_fastqs[sample.name],
                    output_dir
                )
            ], name="run_star_fusion." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def run_arriba(self):
        """
        [arriba](https://github.com/suhrig/arriba) is used for the detection of gene fusions from RNA-Seq data.
        arriba is based on the [STAR](https://github.com/alexdobin/STAR) aligner. Apart from gene fusions,
        Arriba can detect other structural rearrangements with potential clinical relevance, including viral integration
        sites, internal tandem duplications, whole exon duplications and intragenic inversions etc...
        """
    
        jobs = []
    
        left_fastqs = collections.defaultdict(list)
        right_fastqs = collections.defaultdict(list)
    
        for readset in self.readsets:
            trim_dir = os.path.abspath("trim")
            trim_file_prefix = os.path.join(trim_dir, readset.sample.name, readset.name + "-trimmed-")
        
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                  re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
        
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
        
            left_fastqs[readset.sample.name].append(fastq1)
            right_fastqs[readset.sample.name].append(fastq2)
    
        for sample in self.samples:
            output_dir = os.path.join("fusion", sample.name, "arriba")
        
            mkdir_job = Job(command="mkdir -p " + output_dir)
        
            chgdir_job = Job(command="cd " + output_dir)
        
            job = concat_jobs([
                mkdir_job,
                chgdir_job,
                arriba.run(
                    left_fastqs[sample.name],
                    right_fastqs[sample.name],
                    output_dir
                )
            ], name="run_arriba." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def run_annofuse(self):
        """
        [annofuse](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03922-7)
        is a R package and it is used to annotate, prioritize, and interactively explore putative oncogenic
        RNA fusions.
        """
        jobs = []

        for sample in self.samples:
            arriba_output = os.path.join("fusion", sample.name, "arriba",
                                         "fusions.tsv")
            star_fusion_output = os.path.join("fusion", sample.name, "star_fusion",
                                              "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
            output = os.path.join("fusion", sample.name, "annoFuse")

            job = concat_jobs([
                bash.mkdir(
                    output,
                    remove=True
                ),
                annoFuse.run(
                    arriba_output,
                    star_fusion_output,
                    os.path.join(output, sample.name)
                )
            ], name="run_annoFuse." + sample.name)
            
            job.samples = [sample]
            jobs.append(job)
            
        return jobs

    def estimate_ribosomal_rna(self):
        """
        Use bwa mem to align reads on the rRNA reference fasta and count the number of read mapped
        The filtered reads are aligned to a reference fasta file of ribosomal sequence. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files: readset Bam files
        """
    
        jobs = []
        for readset in self.readsets:
            readset_bam = os.path.join("alignment", readset.sample.name, readset.name , "Aligned.sortedByCoord.out.bam")
            output_folder = os.path.join("metrics", readset.sample.name, readset.name)
            readset_metrics_bam = os.path.join(output_folder, readset.name + "rRNA.bam")
        
            job = concat_jobs([
                Job(command="mkdir -p " + os.path.dirname(readset_bam) + " " + output_folder),
                pipe_jobs([
                    bvatools.bam2fq(
                        readset_bam
                    ),
                    bwa.mem(
                        "/dev/stdin",
                        None,
                        read_group="'@RG" + \
                                   "\tID:" + readset.name + \
                                   "\tSM:" + readset.sample.name + \
                                   ("\tLB:" + readset.library if readset.library else "") + \
                                   (
                                       "\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                                   ("\tCN:" + config.param('bwa_mem_rRNA', 'sequencing_center') if config.param(
                                       'bwa_mem_rRNA', 'sequencing_center', required=False) else "") + \
                                   "\tPL:Illumina" + \
                                   "'",
                        ref=config.param('bwa_mem_rRNA', 'ribosomal_fasta'),
                        ini_section='bwa_mem_rRNA'
                    ),
                    picard2.sort_sam(
                        "/dev/stdin",
                        readset_metrics_bam,
                        "coordinate",
                        ini_section='picard_sort_sam_rrna'
                    )
                ]),
                tools.py_rrnaBAMcount(
                    bam=readset_metrics_bam,
                    gtf=config.param('bwa_mem_rRNA', 'gtf'),
                    output=os.path.join(output_folder, readset.name + "rRNA.stats.tsv"),
                    typ="transcript")], name="bwa_mem_rRNA." + readset.name)
        
            job.removable_files = [readset_metrics_bam]
            job.samples = [readset.sample]
            jobs.append(job)
            
        return jobs

    def wiggle(self):
        """
        Generate wiggle tracks suitable for multiple browsers.
        """
    
        jobs = []
    
        ##check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "PAIRED_END"
            if readset.run_type == "SINGLE_END":
                library[readset.sample] = "SINGLE_END"
    
        for sample in self.samples:
            bam_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.")
            input_bam = bam_file_prefix + "bam"
            bed_graph_prefix = os.path.join("tracks", sample.name, sample.name)
            big_wig_prefix = os.path.join("tracks", "bigWig", sample.name)
        
            if (config.param('DEFAULT', 'strand_info') != 'fr-unstranded') and library[sample] == "PAIRED_END":
                input_bam_f1 = bam_file_prefix + "tmp1.forward.bam"
                input_bam_f2 = bam_file_prefix + "tmp2.forward.bam"
                input_bam_r1 = bam_file_prefix + "tmp1.reverse.bam"
                input_bam_r2 = bam_file_prefix + "tmp2.reverse.bam"
                output_bam_f = bam_file_prefix + "forward.bam"
                output_bam_r = bam_file_prefix + "reverse.bam"
            
                bam_f_job = concat_jobs([
                    samtools.view(input_bam, input_bam_f1, "-bh -F 256 -f 81"),
                    samtools.view(input_bam, input_bam_f2, "-bh -F 256 -f 161"),
                    gatk4.merge_sam_files([input_bam_f1, input_bam_f2], output_bam_f),
                    Job(command="rm " + input_bam_f1 + " " + input_bam_f2)
                ], name="wiggle." + sample.name + ".forward_strandspec")
                bam_f_job.samples = [sample]
                # Remove temporary-then-deleted files from job output files, otherwise job is never up to date
                bam_f_job.output_files.remove(input_bam_f1)
                bam_f_job.output_files.remove(input_bam_f2)
            
                bam_r_job = concat_jobs([
                    Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " " + os.path.join("tracks",
                                                                                                       "bigWig")),
                    samtools.view(input_bam, input_bam_r1, "-bh -F 256 -f 97"),
                    samtools.view(input_bam, input_bam_r2, "-bh -F 256 -f 145"),
                    gatk4.merge_sam_files([input_bam_r1, input_bam_r2], output_bam_r),
                    Job(command="rm " + input_bam_r1 + " " + input_bam_r2)
                ], name="wiggle." + sample.name + ".reverse_strandspec")
                bam_r_job.samples = [sample]
                # Remove temporary-then-deleted files from job output files, otherwise job is never up to date
                bam_r_job.output_files.remove(input_bam_r1)
                bam_r_job.output_files.remove(input_bam_r2)
            
                jobs.extend([bam_f_job, bam_r_job])
            
                outputs = [
                    [bed_graph_prefix + ".forward.bedGraph", big_wig_prefix + ".forward.bw"],
                    [bed_graph_prefix + ".reverse.bedGraph", big_wig_prefix + ".reverse.bw"],
                ]
            else:
                outputs = [[bed_graph_prefix + ".bedGraph", big_wig_prefix + ".bw"]]
        
            for bed_graph_output, big_wig_output in outputs:
                if "forward" in bed_graph_output:
                    in_bam = bam_file_prefix + "forward.bam"  # same as output_bam_f from previous picard job
                elif "reverse" in bed_graph_output:
                    in_bam = bam_file_prefix + "reverse.bam"  # same as output_bam_r from previous picard job
                else:
                    in_bam = input_bam
                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + os.path.join("tracks", sample.name) + " ", removable_files=["tracks"],
                            samples=[sample]),
                        bedtools.graph(in_bam, bed_graph_output, library[sample])
                    ], name="bed_graph." + re.sub(".bedGraph", "", os.path.basename(bed_graph_output)))
                )
                jobs.append(
                    concat_jobs([
                        Job(command="mkdir -p " + os.path.join("tracks", "bigWig"), samples=[sample]),
                        ucsc.bedGraphToBigWig(bed_graph_output, big_wig_output, False)
                    ], name="wiggle." + re.sub(".bw", "", os.path.basename(big_wig_output)))
                )
    
        return jobs

    def raw_counts(self):
        """
        Count reads in features using [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html).
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)
            input_bam = alignment_file_prefix + ".QueryNameSorted.bam"

            # Count reads
            output_count = os.path.join("raw_counts", sample.name + ".readcounts.csv")
            stranded = "no" if config.param('DEFAULT', 'strand_info') == "fr-unstranded" else "reverse"
            job = concat_jobs([
                Job(command="mkdir -p raw_counts"),
                pipe_jobs([
                        samtools.view(
                                input_bam,
                                options="-F 4"
                        ),
                        htseq.htseq_count(
                        "-",
                        config.param('htseq_count', 'gtf', param_type='filepath'),
                        output_count,
                        config.param('htseq_count', 'options'),
                        stranded
                        )
                ])
            ], name="htseq_count." + sample.name)
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def raw_counts_metrics(self):
        """
        Create rawcount matrix, zip the wiggle tracks and create the saturation plots based on standardized read counts.
        """

        jobs = []

        # Create raw count matrix
        output_directory = "DGE"
        read_count_files = [os.path.join("raw_counts", sample.name + ".readcounts.csv") for sample in self.samples]
        output_matrix = os.path.join(output_directory, "rawCountMatrix.csv")

        job = Job(read_count_files, [output_matrix], [['raw_counts_metrics', 'module_mugqic_tools']],
                  name="metrics.matrix")

        job.command = """\
mkdir -p {output_directory} && \\
gtf2tmpMatrix.awk \\
  {reference_gtf} \\
  {output_directory}/tmpMatrix.txt && \\
HEAD='Gene\\tSymbol' && \\
for read_count_file in \\
  {read_count_files}
do
  sort -k1,1 $read_count_file > {output_directory}/tmpSort.txt && \\
  join -1 1 -2 1 <(sort -k1,1 {output_directory}/tmpMatrix.txt) {output_directory}/tmpSort.txt > {output_directory}/tmpMatrix.2.txt && \\
  mv {output_directory}/tmpMatrix.2.txt {output_directory}/tmpMatrix.txt && \\
  na=$(basename $read_count_file | rev | cut -d. -f3- | rev) && \\
  HEAD="$HEAD\\t$na"
done && \\
echo -e $HEAD | cat - {output_directory}/tmpMatrix.txt | tr ' ' '\\t' > {output_matrix} && \\
rm {output_directory}/tmpSort.txt {output_directory}/tmpMatrix.txt""".format(
            reference_gtf=config.param('raw_counts_metrics', 'gtf', param_type='filepath'),
            output_directory=output_directory,
            read_count_files=" \\\n  ".join(read_count_files),
            output_matrix=output_matrix
        )
        job.samples = self.samples
        jobs.append(job)

        # Create Wiggle tracks archive
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "PAIRED_END"
            if readset.run_type == "SINGLE_END":
                library[readset.sample] = "SINGLE_END"

        wiggle_directory = os.path.join("tracks", "bigWig")
        wiggle_archive = "tracks.zip"
        if config.param('DEFAULT', 'strand_info') != 'fr-unstranded':
            wiggle_files = []
            for sample in self.samples:
                if library[sample] == "PAIRED_END":
                    wiggle_files.extend([os.path.join(wiggle_directory, sample.name) + ".forward.bw",
                                         os.path.join(wiggle_directory, sample.name) + ".reverse.bw"])
        else:
            wiggle_files = [os.path.join(wiggle_directory, sample.name + ".bw") for sample in self.samples]
        jobs.append(Job(wiggle_files, [wiggle_archive], name="metrics.wigzip",
                        command="zip -r " + wiggle_archive + " " + wiggle_directory, samples=self.samples))

        # RPKM and Saturation
        count_file = os.path.join("DGE", "rawCountMatrix.csv")
        gene_size_file = config.param('rpkm_saturation', 'gene_size', param_type='filepath')
        rpkm_directory = "raw_counts"
        saturation_directory = os.path.join("metrics", "saturation")

        job = concat_jobs([
            Job(command="mkdir -p " + saturation_directory),
            metrics.rpkm_saturation(count_file, gene_size_file, rpkm_directory, saturation_directory)
        ], name="rpkm_saturation")
        job.samples = self.samples
        jobs.append(job)

        report_file = os.path.join("report", "RnaSeq.raw_counts_metrics.md")
        jobs.append(
            Job(
                [wiggle_archive, saturation_directory + ".zip", "metrics/rnaseqRep/corrMatrixSpearman.txt"],
                [report_file],
                [['raw_counts_metrics', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp metrics/rnaseqRep/corrMatrixSpearman.txt report/corrMatrixSpearman.tsv && \\
cp {wiggle_archive} report/ && \\
cp {saturation_archive} report/ && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable corr_matrix_spearman_table="`head -16 report/corrMatrixSpearman.tsv | cut -f-16| awk -F"\t" '{{OFS="\t"; if (NR==1) {{$0="Vs"$0; print; gsub(/[^\t]/, "-"); print}} else {{printf $1; for (i=2; i<=NF; i++) {{printf "\t"sprintf("%.2f", $i)}}; print ""}}}}' | sed 's/\t/|/g'`" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    wiggle_archive=wiggle_archive,
                    saturation_archive=saturation_directory + ".zip",
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="raw_count_metrics_report",
                samples=self.samples
            )
        )

        return jobs

    def stringtie(self):
        """
        Assemble transcriptome using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
        Warning: Still in testing.
        """
        jobs = []

        gtf = config.param('stringtie','gtf', param_type='filepath')

        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.hardClip.bam")
            output_directory = os.path.join("stringtie", sample.name)

            job = stringtie.stringtie(
                input_bam,
                output_directory,
                gtf
            )
            job.name = "stringtie." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def stringtie_merge(self): 
        """
        Merge assemblies into a master teranscriptome reference using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
        Warning: still in testing
        """

        output_directory = os.path.join("stringtie", "AllSamples")
        sample_file = os.path.join("stringtie", "stringtie-merge.samples.txt")
        input_gtfs = [os.path.join("stringtie", sample.name, "transcripts.gtf") for sample in self.samples]
        gtf = config.param('stringtie','gtf', param_type='filepath')


        job = concat_jobs([
            Job(
                command="mkdir -p " + output_directory,
                samples=self.samples
            ),
            Job(
                input_gtfs,
                [sample_file],
                command="""\
`cat > {sample_file} << END
{sample_rows}
END
`""".format(
                    sample_rows="\n".join(input_gtfs),
                    sample_file=sample_file)
            ),
            stringtie.stringtie_merge(
                sample_file,
                output_directory,
                gtf)
        ], name="stringtie-merge"
        )

        return [job]

    def stringtie_abund(self):
        """
        Assemble transcriptome and compute RNA-seq expression using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
        Warning: Still in testing.
        """
        jobs = []

        gtf = os.path.join("stringtie", "AllSamples", "merged.gtf")

        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.hardClip.bam")
            output_directory = os.path.join("stringtie", sample.name)

            job = stringtie.stringtie(
                input_bam,
                output_directory,
                gtf,
                abund=True
            )
            job.name = "stringtie_abund." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    def ballgown(self):
        """
        [Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential
        transcript and gene expression levels and test them for significant differences.
        """

        jobs = []

        # Perform ballgown on each design contrast
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts: 
            design_file = os.path.relpath(self.args.design.name, self.output_dir)        
        output_directory = "ballgown" 
        input_abund = [os.path.join("stringtie", sample.name, "abundance.tab") for sample in self.samples]

        ballgown_job = ballgown.ballgown(
            input_abund,
            design_file,
            output_directory
        )
        ballgown_job.name = "ballgown"
        ballgown_job.samples = self.samples
        jobs.append(ballgown_job)

        return jobs

    def gq_seq_utils_exploratory_analysis_rnaseq(self):
        """
        Exploratory analysis using the gqSeqUtils R package.
        """

        jobs = []

        # gqSeqUtils function call
        sample_fpkm_readcounts = [[
            sample.name,
            os.path.join("cufflinks", sample.name, "isoforms.fpkm_tracking"),
            os.path.join("raw_counts", sample.name + ".readcounts.csv")
        ] for sample in self.samples]
        jobs.append(concat_jobs([
            Job(command="mkdir -p exploratory", samples=self.samples),
            gq_seq_utils.exploratory_analysis_rnaseq(
                os.path.join("DGE", "rawCountMatrix.csv"),
                "cuffnorm",
                config.param('gq_seq_utils_exploratory_analysis_rnaseq', 'genes', param_type='filepath'),
                "exploratory"
            )
        ], name="gq_seq_utils_exploratory_analysis_rnaseq"))

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input            = os.path.join("exploratory", "index.tsv"),
                job_name             = "gq_seq_utils_exploratory_analysis_rnaseq_report",
                input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.Rmd"),
                samples              = self.samples,
                render_output_dir    = 'report',
                module_section       = 'report', # TODO: this or exploratory?
                prerun_r             = 'report_dir="report";' # TODO: really necessary or should be hard-coded in exploratory.Rmd?
            )
        )



        report_file = os.path.join("report", "RnaSeq.cuffnorm.md")
        jobs.append(
            Job(
                [os.path.join("cufflinks", "AllSamples","merged.gtf")],
                [report_file],
                command="""\
mkdir -p report && \\
zip -r report/cuffAnalysis.zip cufflinks/ cuffdiff/ cuffnorm/ && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="cuffnorm_report",
                samples=self.samples
            )
        )

        return jobs


    def differential_expression(self):
        """
        Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
        Merge the results of the analysis in a single csv file.
        """

        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)
        output_directory = "DGE"
        count_matrix = os.path.join(output_directory, "rawCountMatrix.csv")

        edger_job = differential_expression.edger(design_file, count_matrix, output_directory)
        edger_job.output_files = [os.path.join(output_directory, contrast.name, "edger_results.csv") for contrast in self.contrasts]
        edger_job.samples = self.samples

        deseq_job = differential_expression.deseq2(design_file, count_matrix, output_directory)
        deseq_job.output_files = [os.path.join(output_directory, contrast.name, "dge_results.csv") for contrast in self.contrasts]
        deseq_job.samples = self.samples

        return [concat_jobs([
            Job(command="mkdir -p " + output_directory),
            edger_job,
            deseq_job
        ], name="differential_expression")]

    def differential_expression_goseq(self):
        """
        Gene Ontology analysis for RNA-Seq using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
        Generates GO annotations for differential gene expression analysis.
        """

        jobs = []

    def ihec_metrics(self):
        """
        Generate IHEC's standard metrics.
        """

        genome = config.param('ihec_metrics', 'assembly')

        return [metrics.ihec_metrics_rnaseq(genome)]


    @property
    def steps(self):
        return [
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.star,
                self.picard_merge_sam_files,
                self.picard_sort_sam,
                self.picard_mark_duplicates,
                self.picard_rna_metrics,
                self.estimate_ribosomal_rna,
                self.bam_hard_clip,
                self.rnaseqc,
                self.wiggle,
                self.raw_counts,
                self.raw_counts_metrics,
                self.stringtie,
                self.stringtie_merge,
                self.stringtie_abund,
                self.ballgown,
                self.differential_expression,
                self.cram_output
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.star,
                self.picard_merge_sam_files,
                self.picard_mark_duplicates,
                self.split_N_trim,
                self.sambamba_merge_splitNtrim_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.recalibration,
                self.gatk_haplotype_caller,
                self.merge_hc_vcf,
                self.run_vcfanno,
                self.variant_filtration,
                self.decompose_and_normalize,
                self.compute_snp_effects,
                self.gemini_annotations,
                self.picard_rna_metrics,
                self.estimate_ribosomal_rna,
                self.rnaseqc2,
                self.gatk_callable_loci,
                self.wiggle,
                self.cram_output
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.star,
                self.picard_merge_sam_files,
                #self.picard_sort_sam,
                self.picard_mark_duplicates,
                self.split_N_trim,
                self.sambamba_merge_splitNtrim_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.recalibration,
                self.gatk_haplotype_caller,
                self.merge_hc_vcf,
                self.run_vcfanno,
                self.variant_filtration,
                self.decompose_and_normalize,
                self.compute_snp_effects,
                self.gemini_annotations,
                self.run_star_fusion,
                self.run_arriba,
                self.run_annofuse,
                self.picard_rna_metrics,
                self.estimate_ribosomal_rna,
                self.rnaseqc2,
                self.rseqc,
                self.gatk_callable_loci,
                self.wiggle,
                self.cram_output
            ]
        ]
class RnaSeq(RnaSeqRaw):
    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-t", "--type", help="RNAseq analysis type", choices=["stringtie","variants","cancer"], default="stringtie")
        super(RnaSeq, self).__init__(protocol)

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        RnaSeq(protocol=['stringtie','variants', 'cancer'])
