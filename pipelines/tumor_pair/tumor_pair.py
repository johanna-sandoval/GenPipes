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
import logging
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs, pipe_jobs
from bfx.sample_tumor_pairs import parse_tumor_pair_file
from bfx.sequence_dictionary import split_by_size, parse_sequence_dictionary_file
import utils.utils

import gzip
from sys import stderr
from pipelines.dnaseq import dnaseq

#utilizes
from bfx import sambamba
from bfx import bcftools
from bfx import tools
from bfx import metric_tools
from bfx import bvatools
from bfx import vt
from bfx import snpeff
from bfx import vawk
from bfx import deliverables
from bfx import bash_cmd as bash

#metrics
from bfx import conpair
from bfx import qualimap
from bfx import adapters
from bfx import fastqc
from bfx import multiqc

#variants
from bfx import htslib
from bfx import samtools
from bfx import varscan
from bfx import gatk
from bfx import gatk4
from bfx import vardict
from bfx import strelka2
from bfx import bcbio_variation_recall
from bfx import cpsr
from bfx import pcgr
from bfx import gemini

#sv
from bfx import delly
from bfx import manta
from bfx import lumpy
from bfx import svtyper
from bfx import wham
from bfx import metasv
from bfx import cnvkit
from bfx import scones
from bfx import sequenza
from bfx import amber
from bfx import cobalt
from bfx import purple
from bfx import svaba
from bfx import annotations

log = logging.getLogger(__name__)

class TumorPair(dnaseq.DnaSeqRaw):
    """
    Tumor Pair Pipeline
    =================

    The Tumor Pair pipeline inherits the initial bam preparation steps of the DNA-Seq pipeline with the exception of the
    indel realignment (IR) step. In the tumor pipeline the IR step utilizes both the normal and tumor bam to further reduce
    false positives (FPs) in and around indels. The tumor pipeline deviates from the DNA-seq pipeline at the variant calling step. 
    At this point, a paired caller is used to call SNVs and Indels from the pairs given as input. Additional, muliple cancer callers 
    are utilized using an ensemble approach and SNVs and Indels seen in at least 2 different callers are retained for further 
    investigation.

    Example command:
    python tumor_pair.py -c a.ini b.base.ini -s x-y,z -r readset.tsv -p pairs.csv
    
    -c ini files: multiple can be specified e.g WGS or exome, or different clusters e.g. base (abacus) or guillimin

    -r readset: derived from GQ lims or made yourself. See : https://bitbucket.org/mugqic/mugqic_pipelines#markdown-header-readset-file

    -p pairs : format - patient_name,normal_sample_name,tumor_sample_name 
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=argparse.FileType('r'))
        self.argparser.add_argument("--profyle", help="adjust deliverables to PROFYLE folder conventions (Default: False)", action="store_true")
        self.argparser.add_argument("-t", "--type", help="Tumor pair analysis type", choices = ["fastpass", "ensemble", "sv"], default="ensemble")
        super(TumorPair, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
            'alignment_directory': 'alignment',
            'metrics_directory': 'metrics',
            'paired_variants_directory': 'pairedVariants',
            'sv_variants_directory': 'SVariants'
        }
        return dirs

    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(
                self.args.pairs.name,
                self.samples,
                self.args.profyle
            )
        return self._tumor_pairs

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(
                config.param('DEFAULT', 'genome_dictionary', param_type='filepath'),
                variant=True
            )
        return self._sequence_dictionary_variant

    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary_variant()):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in self.sequence_dictionary_variant()]
        else:
            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary_variant()])
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary_variant()))))
            windows = []

            for sequence in self.sequence_dictionary_variant():
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in range(1, sequence['length'] + 1, approximate_window_size)]:
                    windows.append(sequence['name'] + ":" + str(start) + "-" + str(end))

        return windows

    def is_gz_file(self, name):
        if not os.path.isfile(name):
            return True
        #if os.stat(name).st_size == 0:
        #    return False

        with gzip.open(name, 'rb') as f:
            try:
                file_content = f.read(1)
                return len(file_content) > 0
            except:
                return False

    def sym_link_fastq_pair(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            inputs["Normal"] = [
                self.select_input_files(
                    [
                        [readset.fastq1],
                        [os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".pair1.fastq.gz")]
                    ]
                ) for readset in tumor_pair.readsets[tumor_pair.normal.name]
            ][0]
            inputs["Normal"].append(
                [
                    self.select_input_files(
                        [
                            [readset.fastq2],
                            [os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".pair2.fastq.gz")]
                        ]
                    ) for readset in tumor_pair.readsets[tumor_pair.normal.name]
                ][0][0]
            )

            inputs["Tumor"] = [
                self.select_input_files(
                    [
                        [readset.fastq1],
                        [os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".pair1.fastq.gz")]
                    ]
                ) for readset in tumor_pair.readsets[tumor_pair.tumor.name]
            ][0]
            inputs["Tumor"].append(
                [
                    self.select_input_files(
                        [
                            [readset.fastq2],
                            [os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".pair2.fastq.gz")]
                        ]
                    ) for readset in tumor_pair.readsets[tumor_pair.tumor.name]
                ][0][0]
            )
            
            for key, input_files in inputs.items():
                for read, input_file in enumerate(input_files):
                    symlink_pair_job = deliverables.sym_link_pair(
                        input_file,
                        tumor_pair,
                        self.output_dir,
                        type="raw_reads",
                        sample=key,
                        profyle=self.args.profyle
                    )
                    dir_name, file_name = os.path.split(symlink_pair_job.output_files[0])
                    # do not compute md5sum in the readset input directory
                    md5sum_job = deliverables.md5sum(
                        symlink_pair_job.output_files[0],
                        file_name + ".md5",
                        dir_name
                    )
                    jobs.append(
                        concat_jobs(
                            [
                                symlink_pair_job,
                                md5sum_job
                            ],
                            name="sym_link_fastq.pairs." + str(read) + "." + tumor_pair.name + "." + key
                        )
                    )

        return jobs

    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.

        Note: modified to use both normal and tumor bams to reduce FPs around indels

        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            pair_directory = os.path.join(self.output_dirs['alignment_directory'], "realign", tumor_pair.name)

            input_normal = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_intervals = os.path.join(pair_directory, "all.intervals")
                bam_postfix = ".realigned.all.bam"
                
                normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned.all.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                normal_output_bam = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                
                tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned.all.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                tumor_output_bam = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    bash.mkdir(
                        pair_directory,
                        remove=True
                    ),
                    bash.chgdir(
                        pair_directory
                    ),
                    gatk.realigner_target_creator(
                        input_normal,
                        realign_intervals,
                        output_dir=self.output_dir,
                        input2=input_tumor
                    ),
                    gatk.indel_realigner(
                        input_normal,
                        input2=input_tumor,
                        output_dir=self.output_dir,
                        output_norm_dep=[normal_bam,normal_index],
                        output_tum_dep=[tumor_bam,tumor_index],
                        target_intervals=realign_intervals,
                        optional=bam_postfix
                    ),
                    # Move sample realign
                    bash.ln(
                        normal_bam,
                        normal_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        normal_index,
                        normal_output_index,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_bam,
                        tumor_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_index,
                        tumor_output_index,
                        self.output_dir
                    ),
                ], name="gatk_indel_realigner." + tumor_pair.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)
                normal_realign_directory = os.path.join(normal_alignment_directory, "realign")
                tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")
                
                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(pair_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequences
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    bam_postfix = ".realigned." + str(idx) + ".bam"
                    normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_index = re.sub("\.bam$", ".bai", normal_bam)
                    tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                    normal_output_bam = os.path.join(normal_realign_directory,
                                                     tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                    tumor_output_bam = os.path.join(tumor_realign_directory,
                                                    tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            pair_directory,
                            remove=True
                        ),
                        bash.mkdir(
                            normal_realign_directory,
                            remove=True
                        ),
                        bash.mkdir(
                            tumor_realign_directory,
                            remove=True
                        ),
                        bash.chgdir(
                            pair_directory
                        ),
                        gatk.realigner_target_creator(
                            input_normal,
                            realign_intervals,
                            output_dir=self.output_dir,
                            input2=input_tumor,
                            intervals=intervals
                        ),
                        gatk.indel_realigner(
                            input_normal,
                            input2=input_tumor,
                            output_dir=self.output_dir,
                            output_norm_dep=[normal_bam,normal_index],
                            output_tum_dep=[tumor_bam,tumor_index],
                            target_intervals=realign_intervals,
                            intervals=intervals,
                            optional=bam_postfix
                        ),
                        bash.ln(
                            normal_bam,
                            normal_output_bam,
                            self.output_dir
                        ),
                        bash.ln(
                            normal_index,
                            normal_output_index,
                            self.output_dir
                        ),
                        bash.ln(
                            tumor_bam,
                            tumor_output_bam,
                            self.output_dir
                        ),
                        bash.ln(
                            tumor_index,
                            tumor_output_index,
                            self.output_dir
                        ),
                    ], name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_intervals = os.path.join(pair_directory, "others.intervals")
                bam_postfix = ".realigned.others.bam"
                normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                normal_output_bam = os.path.join(normal_realign_directory,
                                                 tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_output_bam = os.path.join(tumor_realign_directory,
                                                tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    bash.mkdir(
                        pair_directory,
                        remove=True
                    ),
                    bash.mkdir(
                        normal_realign_directory,
                        remove=True
                    ),
                    bash.mkdir(
                        tumor_realign_directory,
                        remove=True
                    ),
                    bash.chgdir(
                        pair_directory
                    ),
                    gatk.realigner_target_creator(
                        input_normal,
                        realign_intervals,
                        output_dir=self.output_dir,
                        input2=input_tumor,
                        exclude_intervals=unique_sequences_per_job_others
                    ),
                    gatk.indel_realigner(
                        input_normal,
                        input2=input_tumor,
                        output_dir=self.output_dir,
                        output_norm_dep=[normal_bam, normal_index],
                        output_tum_dep=[tumor_bam, tumor_index],
                        target_intervals=realign_intervals,
                        exclude_intervals=unique_sequences_per_job_others,
                        optional=bam_postfix
                    ),
                    bash.ln(
                        normal_bam,
                        normal_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        normal_index,
                        normal_output_index,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_bam,
                        tumor_output_bam,
                        self.output_dir
                    ),
                    bash.ln(
                        tumor_index,
                        tumor_output_index,
                        self.output_dir
                    ),
                ], name="gatk_indel_realigner." + tumor_pair.name + ".others"))

        return jobs

    def sambamba_merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using
        [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
 
            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_jobs - 1)

                normal_inputs = []
                for idx, sequences in enumerate(unique_sequences_per_job):
                    normal_inputs.append(
                        os.path.join(
                            normal_alignment_directory,
                            "realign",
                            tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam"
                        )
                    )
                normal_inputs.append(
                    os.path.join(
                        normal_alignment_directory,
                        "realign",
                        tumor_pair.normal.name + ".sorted.realigned.others.bam"
                    )
                )

                tumor_inputs = []
                for idx, sequences in enumerate(unique_sequences_per_job):
                    tumor_inputs.append(
                        os.path.join(
                            tumor_alignment_directory,
                            "realign",
                            tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam"
                        )
                    )
                tumor_inputs.append(
                    os.path.join(
                        tumor_alignment_directory,
                        "realign",
                        tumor_pair.tumor.name + ".sorted.realigned.others.bam"
                    )
                )

                job = sambamba.merge(
                    normal_inputs,
                    os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")
                )
                job.name = "sambamba_merge_realigned." + tumor_pair.name + "." + tumor_pair.normal.name
                job.samples = [tumor_pair.normal]
                jobs.append(job)

                job = sambamba.merge(
                    tumor_inputs,
                    os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")
                )
                job.name = "sambamba_merge_realigned." + tumor_pair.name + "." + tumor_pair.tumor.name
                job.samples = [tumor_pair.tumor]
                jobs.append(job)

        return jobs

    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using
        [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            [normal_input] = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")],
            ])
            normal_output = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")

            [tumor_input] = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")],
            ])
            tumor_output = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name  + ".sorted.dup.bam")

            job = sambamba.markdup(
                normal_input,
                normal_output,
                config.param('sambamba_mark_duplicates', 'tmp_dir'),
                other_options=config.param('sambamba_mark_duplicates', 'options')
            )
            job.name = "sambamba_mark_duplicates." + tumor_pair.name + "." + tumor_pair.normal.name
            #job.samples = [tumor_pair.normal]
            jobs.append(job)

            job = sambamba.markdup(
                tumor_input,
                tumor_output,
                config.param('sambamba_mark_duplicates', 'tmp_dir'),
                other_options=config.param('sambamba_mark_duplicates', 'options')
            )
            job.name = "sambamba_mark_duplicates." + tumor_pair.name + "." + tumor_pair.tumor.name
            #job.samples = [tumor_pair.tumor]
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

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            normal_prefix = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.")
            tumor_prefix = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.")
            
            normal_input = normal_prefix + "bam"
            tumor_input = tumor_prefix + "bam"
            
            normal_print_reads_output = normal_prefix + "recal.bam"
            tumor_print_reads_output = tumor_prefix + "recal.bam"
            
            normal_base_recalibrator_output = normal_prefix + "recalibration_report.grp"
            tumor_base_recalibrator_output = tumor_prefix + "recalibration_report.grp"
            
            interval_list = None
        
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )
            if coverage_bed:
                interval_list = os.path.join(tumor_alignment_directory, re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed)))
            
                if not os.path.isfile(interval_list):
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(tumor_alignment_directory),
                                tools.bed2interval_list(
                                    coverage_bed,
                                    interval_list
                                )
                            ],
                            name="interval_list." + os.path.basename(coverage_bed)
                        )
                    )
        
            job = gatk4.base_recalibrator(
                normal_input,
                normal_base_recalibrator_output,
                intervals=interval_list
            )
            job.name = "gatk_base_recalibrator." + tumor_pair.name + "." + tumor_pair.normal.name
            jobs.append(job)
        
            job = gatk4.print_reads(
                normal_input,
                normal_print_reads_output,
                normal_base_recalibrator_output
            )
            job.name = "gatk_print_reads." + tumor_pair.name + "." + tumor_pair.normal.name
            jobs.append(job)

            job = gatk4.base_recalibrator(
                tumor_input,
                tumor_base_recalibrator_output,
                intervals=interval_list
            )
            job.name = "gatk_base_recalibrator." + tumor_pair.name + "." + tumor_pair.tumor.name
            jobs.append(job)

            job = gatk4.print_reads(
                tumor_input,
                tumor_print_reads_output,
                tumor_base_recalibrator_output
            )
            job.name = "gatk_print_reads." + tumor_pair.name + "." + tumor_pair.tumor.name
            jobs.append(job)

        return jobs

    def sym_link_final_bam(self):
        """
        Create sym link of final bam for delivery of data to clients
        :return:
        """
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            inputs["Normal"] = [self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")]
            ])][0]

            inputs["Normal"].append(self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bai")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam.bai")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bai")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam.bai")]
            ])[0])

            inputs["Tumor"] = [self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")]
            ])][0]
            
            inputs["Tumor"].append(self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bai")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam.bai")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bai")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam.bai")]
            ])[0])

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            input_file,
                            input_file + ".md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            input_file,
                            tumor_pair,
                            self.output_dir,
                            type="alignment",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            input_file + ".md5",
                            tumor_pair,
                            self.output_dir,
                            type="alignment",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_final_bam.pairs." + str(idx) + "." + tumor_pair.name + "." + key))

        return jobs

    def conpair_concordance_contamination(self):
        """
        Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification
        (= samples coming from the same individual), as well as cross-individual contamination level estimation in
        whole-genome and whole-exome sequencing experiments. Importantly, the method of estimates contamination in
        the tumor samples not affected by copy number changes and is able to detect contamination levels as low as 0.1%.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'])

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")
            pileup_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".gatkPileup")
            pileup_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".gatkPileup")

            concordance_out = os.path.join(metrics_directory, tumor_pair.tumor.name + ".concordance.tsv")
            contamination_out = os.path.join(metrics_directory, tumor_pair.tumor.name + ".contamination.tsv")

            jobs.append(concat_jobs([
                conpair.pileup(
                    input_normal,
                    pileup_normal
                ),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.name + "." + tumor_pair.normal.name))

            jobs.append(concat_jobs([
                conpair.pileup(
                    input_tumor,
                    pileup_tumor
                ),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.name + "." + tumor_pair.tumor.name))

            jobs.append(concat_jobs([
                bash.mkdir(
                    metrics_directory,
                    remove=False
                ),
                conpair.concordance(
                    pileup_normal,
                    pileup_tumor,
                    concordance_out
                ),
                conpair.contamination(
                    pileup_normal,
                    pileup_tumor,
                    contamination_out
                )
            ], name="conpair_concordance_contamination." + tumor_pair.name))

        return jobs

    def rawmpileup_panel(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup_panel', 'nb_jobs', param_type='posint')
            bedfile = config.param('rawmpileup_panel', 'panel')

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")
    
                jobs.append(concat_jobs([
                    bash.mkdir(
                        varscan_directory,
                        remove=True
                    ),
                    samtools.mpileup(
                        [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam"),
                         os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                        input_pair,
                        config.param('rawmpileup_panel', 'mpileup_other_options'),
                        regionFile=bedfile
                    ),
                    ], name="rawmpileup_panel." + tumor_pair.name + ".all")
                )
                
            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")
    
                        jobs.append(concat_jobs([
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            samtools.mpileup(
                                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam"),
                                 os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                                pair_output,
                                config.param('rawmpileup_panel', 'mpileup_other_options'),
                                region=sequence['name'],
                                regionFile=bedfile
                            ),
                            ], name="rawmpileup_panel." + tumor_pair.name + "." + sequence['name'])
                        )
        return jobs

    def paired_varscan2_panel(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup_panel', 'nb_jobs', param_type='posint')

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")

                output = os.path.join(varscan_directory, tumor_pair.name)
                output_snp = os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")

                jobs.append(concat_jobs([
                    bash.mkdir(
                        varscan_directory,
                        remove=True
                    ),
                    varscan.somatic(
                        input_pair,
                        output,
                        config.param('varscan2_somatic_panel', 'other_options'),
                        output_vcf_dep=output_vcf_gz,
                        output_snp_dep=output_snp,
                        output_indel_dep=output_indel
                    ),
                    htslib.bgzip_tabix(
                        output_snp,
                        os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz")
                    ),
                    htslib.bgzip_tabix(
                        output_indel,
                        os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
                    ),
                    pipe_jobs([
                        bcftools.concat(
                            [os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz"),
                             os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")],
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="sed 's/TUMOR/"
                                    + tumor_pair.tumor.name
                                    + "/g' | sed 's/NORMAL/"
                                    + tumor_pair.normal.name + "/g' "
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vcf_gz
                        ),
                    ]),
                ], name="varscan2_somatic_panel." + tumor_pair.name + ".all"))
                
            else:
                
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")
    
                        output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                        output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                        output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                        output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")
    
                        jobs.append(concat_jobs([
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            varscan.somatic(
                                input_pair,
                                output,
                                config.param('varscan2_somatic_panel', 'other_options'),
                                output_vcf_dep=output_vcf_gz,
                                output_snp_dep=output_snp,
                                output_indel_dep=output_indel
                            ),
                            htslib.bgzip_tabix(
                                output_snp,
                                os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz")
                            ),
                            htslib.bgzip_tabix(
                                output_indel,
                                os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")
                            ),
                            pipe_jobs([
                                bcftools.concat(
                                    [os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz"),
                                     os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")],
                                    None
                                ),
                                Job(
                                    [None],
                                    [None],
                                    command="sed 's/TUMOR/"+ tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"+ tumor_pair.normal.name + "/g' "
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_vcf_gz
                                ),
                            ]),
                        ], name="varscan2_somatic_panel." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_varscan2_panel(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup_panel', 'nb_jobs', param_type='posint')

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    pipe_jobs([
                        Job(
                            [os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")],
                            [None],
                            command="zcat " + os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")
                        ),
                        tools.fix_varscan_output(
                            None,
                            None,
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job([None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                            ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                        ),
                    ]),
                    bcftools.view(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                        config.param('merge_varscan2', 'somatic_filter_options')
                    ),
                    htslib.tabix(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                        config.param('merge_varscan2', 'tabix_options', required=False)
                    ),
                    bcftools.view(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz"),
                        config.param('merge_varscan2', 'germline_filter_options')
                    ),
                    htslib.tabix(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz"),
                        config.param('merge_varscan2', 'tabix_options', required=False)
                    ),
                ], name = "merge_varscan2." + tumor_pair.name))

            else:
                all_inputs = [os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")
                              for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary']

                for input_vcf in all_inputs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete panel varscan2 vcf: %s\n" % input_vcf)

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            all_inputs,
                            None),
                        tools.fix_varscan_output(
                            None,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job([None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                            ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                        ),
                    ]),
                    bcftools.view(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                        config.param('merge_varscan2', 'somatic_filter_options')
                    ),
                    htslib.tabix(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                        config.param('merge_varscan2', 'tabix_options', required=False)
                    ),
                    bcftools.view(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz"),
                        config.param('merge_varscan2', 'germline_filter_options')
                    ),
                    htslib.tabix(
                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz"),
                        config.param('merge_varscan2', 'tabix_options', required=False)
                    ),
                ], name="merge_varscan2." + tumor_pair.name))

        return jobs

    def preprocess_vcf_panel(self):
        """
        Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
        vcf FORMAT modification for correct loading into gemini
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")

            prefix = os.path.join(pair_directory, tumor_pair.name)
            output_somatic = prefix + ".varscan2.somatic.vt.vcf.gz"

            output_germline = prefix + ".varscan2.germline.vt.vcf.gz"

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        prefix + ".varscan2.somatic.vcf.gz" ,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        prefix + ".prep.vt.vcf.gz"
                    ),
                ]),
                tools.preprocess_varscan(
                    prefix + ".prep.vt.vcf.gz",
                    output_somatic
                ),
            ], name="preprocess_vcf_panel.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        prefix + ".varscan2.germline.vcf.gz" ,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        prefix + ".germline.prep.vt.vcf.gz"
                    ),
                ]),
                tools.preprocess_varscan(
                    prefix + ".germline.prep.vt.vcf.gz",
                    output_germline
                ),
            ], name="preprocess_vcf_panel.germline." + tumor_pair.name))

        return jobs

    def snp_effect_panel(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            input_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf")
            output_somatic_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf.gz")

            input_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")
            output_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.snpeff.vcf")
            output_germline_gz = os.path.join(pair_directory,
                                              tumor_pair.name + ".varscan2.germline.vt.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(varscan_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(
                    input_somatic,
                    output_somatic,
                    cancer_sample_file=cancer_pair_filename,
                    options=config.param('compute_cancer_effects_somatic', 'options')
                ),
                htslib.bgzip_tabix(
                    output_somatic,
                    output_somatic_gz
                ),
            ], name = "compute_cancer_effects_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(
                    input_germline,
                    output_germline,
                    cancer_sample_file=cancer_pair_filename,
                    options=config.param('compute_cancer_effects_germline', 'options')
                ),
                htslib.bgzip_tabix(
                    output_germline,
                    output_germline_gz
                ),
            ], name = "compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def gemini_annotations_panel(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database
        [Gemini] (http://gemini.readthedocs.org/en/latest/index.html)
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            temp_dir = config.param('DEFAULT', 'tmp_dir')
            gemini_prefix = os.path.join(pair_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(
                    gemini_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",
                    gemini_prefix + ".somatic.gemini.db", temp_dir
                )
            ], name="gemini_annotations.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(
                    gemini_prefix + ".varscan2.germline.vt.snpeff.vcf.gz",
                    gemini_prefix + ".germline.gemini.db",
                    temp_dir
                )
            ], name="gemini_annotations.germline." + tumor_pair.name))

        return jobs

    def sym_link_panel(self):
        """
        Create sym links of panel variants for deliverables to the clients
        """
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] =  [os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel", tumor_pair.name)]

            for key, input_files in inputs.items():
                for idx, sample_prefix in enumerate(input_files):
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(
                            sample_prefix + ".varscan2.vcf.gz",
                            tumor_pair, self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".varscan2.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle),
                        deliverables.sym_link_pair(
                            sample_prefix + ".varscan2.germline.vt.snpeff.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".varscan2.germline.vt.snpeff.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".somatic.gemini.db",
                            tumor_pair,
                            self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".germline.gemini.db",
                            tumor_pair, self.output_dir,
                            type="snv/panel",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_panel." + str(idx) + "." + tumor_pair.name + "." + key))

        return jobs

    def metrics_dna_picard_metrics(self):
        """
        Runs specific QC metrics on DNA data
        Functions: collect_multiple_metrics, CollectOxoGMetrics and collect_sequencing_artifacts_metrics
        [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html)
        """
    
        ffpe = config.param('picard_collect_sequencing_artifacts_metrics', 'FFPE', param_type='boolean')

        ##check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                normal_metrics = os.path.join(tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                normal_metrics = os.path.join(tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            normal_picard_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", normal_metrics, "picard_metrics")
            tumor_picard_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "picard_metrics")

            [normal_input] = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")],

            ])

            [tumor_input] = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")],

            ])
            # log.info(input)
            mkdir_job_normal = bash.mkdir(
                normal_picard_directory,
                remove=True
            )

            jobs.append(
                concat_jobs([
                    mkdir_job_normal,
                    gatk4.collect_multiple_metrics(
                        normal_input,
                        os.path.join(normal_picard_directory, tumor_pair.normal.name + ".all.metrics"),
                        library_type=library[tumor_pair.normal]
                    )
                ],
                    name="picard_collect_multiple_metrics." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal]
                )
            )

            jobs.append(
                concat_jobs([
                    mkdir_job_normal,
                    gatk4.collect_oxog_metrics(
                        normal_input,
                        os.path.join(normal_picard_directory, tumor_pair.normal.name + ".oxog_metrics.txt")
                    )
                ],
                    name="picard_collect_oxog_metrics." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal]
                )
            )
        
            jobs.append(
                concat_jobs([
                    mkdir_job_normal,
                    gatk4.collect_gcbias_metrics(
                        normal_input,
                        os.path.join(normal_picard_directory, tumor_pair.normal.name + ".qcbias_metrics.txt"),
                        os.path.join(normal_picard_directory, tumor_pair.normal.name + ".qcbias_metrics.pdf"),
                        os.path.join(normal_picard_directory, tumor_pair.normal.name + ".qcbias_summary_metrics.txt")
                    )
                ],
                    name="picard_collect_gcbias_metrics." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal]
                )
            )
            # log.info(input)
            mkdir_job_tumor = bash.mkdir(
                tumor_picard_directory,
                remove=True
            )

            jobs.append(
                concat_jobs([
                    mkdir_job_tumor,
                    gatk4.collect_multiple_metrics(
                        tumor_input,
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".all.metrics"),
                        library_type=library[tumor_pair.tumor]
                    )
                ],
                    name="picard_collect_multiple_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs([
                    mkdir_job_tumor,
                    gatk4.collect_oxog_metrics(
                        tumor_input,
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".oxog_metrics.txt")
                    )
                ],
                    name="picard_collect_oxog_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs([
                    mkdir_job_tumor,
                    gatk4.collect_gcbias_metrics(
                        tumor_input,
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".qcbias_metrics.txt"),
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".qcbias_metrics.pdf"),
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".qcbias_summary_metrics.txt")
                    )
                ],
                    name="picard_collect_gcbias_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor]
                )
            )

            if ffpe == True:
                jobs.append(concat_jobs([
                    mkdir_job_normal,
                    gatk4.collect_sequencing_artifacts_metrics(
                        normal_input,
                        os.path.join(normal_picard_directory, tumor_pair.normal.name)
                    )
                ],
                    name="picard_collect_sequencing_artifacts_metrics." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal]
                )
                )
                jobs.append(concat_jobs([
                    mkdir_job_tumor,
                    gatk4.collect_sequencing_artifacts_metrics(
                        tumor_input,
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name)
                    )
                ],
                    name="picard_collect_sequencing_artifacts_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor]
                )
                )

        return jobs

    def metrics_dna_sample_qualimap(self):
        """
        QC alignment metrics generated by
        [Qualimap](http://qualimap.conesalab.org/)
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                normal_metrics = os.path.join(tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                normal_metrics = os.path.join(tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            normal_qualimap_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", normal_metrics,
                                                     "qualimap", tumor_pair.normal.name)
            tumor_qualimap_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name,
                                                     "qualimap", tumor_pair.tumor.name)

            [normal_input] = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
            ])
        
            normal_output = os.path.join(normal_qualimap_directory, "genome_results.txt")

            [tumor_input] = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
            ])

            tumor_output = os.path.join(tumor_qualimap_directory, "genome_results.txt")
            use_bed = config.param('dna_sample_qualimap', 'use_bed', param_type='boolean', required=True)
        
            options = None
            if use_bed:
                bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])
                options = config.param('dna_sample_qualimap', 'qualimap_options') + " --feature-file " + bed
        
            else:
                options = config.param('dna_sample_qualimap', 'qualimap_options')
        
            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        normal_qualimap_directory,
                        remove=False
                    ),
                    qualimap.bamqc(
                        normal_input,
                        normal_qualimap_directory,
                        normal_output,
                        options
                    )
                ],
                    name="dna_sample_qualimap." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal]
                )
            )
            
            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        tumor_qualimap_directory,
                        remove=False
                    ),
                    qualimap.bamqc(
                        tumor_input,
                        tumor_qualimap_directory,
                        tumor_output,
                        options
                    )
                ],
                    name="dna_sample_qualimap." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor]
                )
            )
    
    
        return jobs

    def metrics_dna_fastqc(self):
        """
        QCing metrics generated on the read level using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        """
    
        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                normal_metrics = os.path.join(tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                normal_metrics = os.path.join(tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            normal_fastqc_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", normal_metrics, "fastqc")

            tumor_fastqc_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "fastqc")
  
            [normal_input] = self.select_input_files([
                # [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
            ])
            
            normal_output_dir = normal_fastqc_directory
            normal_file = re.sub(".bam", "", os.path.basename(normal_input))
            normal_output = os.path.join(normal_fastqc_directory, normal_file + "_fastqc.zip")
            
            [tumor_input] = self.select_input_files([
                # [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
            ])
        
            tumor_output_dir = tumor_fastqc_directory
            tumor_file = re.sub(".bam", "", os.path.basename(tumor_input))
            tumor_output = os.path.join(tumor_fastqc_directory, tumor_file + "_fastqc.zip")
        
            adapter_file = config.param('fastqc', 'adapter_file', required=False, param_type='filepath')
            normal_adapter_job = None
            tumor_adapter_job = None
        
            if not adapter_file:
                normal_adapter_job = adapters.create(
                    tumor_pair.normal.readsets[0],
                    os.path.join(normal_output_dir, "adapter.tsv"),
                    fastqc=True
                )
                tumor_adapter_job = adapters.create(
                    tumor_pair.tumor.readsets[0],
                    os.path.join(tumor_output_dir, "adapter.tsv"),
                    fastqc=True
                )
        
            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        normal_output_dir,
                        remove=True
                    ),
                    normal_adapter_job,
                    fastqc.fastqc(
                        normal_input,
                        None,
                        normal_output_dir,
                        normal_output,
                        os.path.join(normal_output_dir, "adapter.tsv")
                    )
                ],
                    name="fastqc." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal]
                )
            )
            
            jobs.append(
                concat_jobs([
                    bash.mkdir(
                        tumor_output_dir,
                        remove=True
                    ),
                    tumor_adapter_job,
                    fastqc.fastqc(
                        tumor_input,
                        None,
                        tumor_output_dir,
                        tumor_output,
                        os.path.join(tumor_output_dir, "adapter.tsv")
                    )
                ],
                    name="fastqc." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor]
                )
            )
    
        return jobs

    def run_pair_multiqc(self):
        """
        Aggregate results from bioinformatics analyses across many samples into a single report
        MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
        perfect for summarising the output from numerous bioinformatics tools
        [MultiQC](https://multiqc.info/)
        """

        jobs = []

        metrics_directory = os.path.join(self.output_dirs['metrics_directory'], "dna")
        input_dep = []
        for tumor_pair in self.tumor_pairs.values():
            input_dep = []
            inputs = []
            if tumor_pair.multiple_normal == 1:
                normal_directory = os.path.join(metrics_directory, tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_directory = os.path.join(metrics_directory, tumor_pair.normal.name)
    
            tumor_directory = os.path.join(metrics_directory, tumor_pair.tumor.name)

            input_normal_oxog = os.path.join(normal_directory, "picard_metrics", tumor_pair.normal.name + ".oxog_metrics.txt")
            input_normal_qcbias = os.path.join(normal_directory, "picard_metrics", tumor_pair.normal.name +".qcbias_metrics.txt")
            input_normal_all_picard = os.path.join(normal_directory, "picard_metrics", tumor_pair.normal.name + ".all.metrics.quality_distribution.pdf")
            input_normal_qualimap = os.path.join(normal_directory, "qualimap", tumor_pair.normal.name, "genome_results.txt")
            
            [input_normal_fastqc] = self.select_input_files([
                [os.path.join(normal_directory, "fastqc", tumor_pair.normal.name + ".sorted.dup_fastqc.zip")],
                [os.path.join(normal_directory, "fastqc", tumor_pair.normal.name + "_fastqc.zip")],
            ])

            input_tumor_oxog = os.path.join(tumor_directory, "picard_metrics", tumor_pair.tumor.name + ".oxog_metrics.txt")
            input_tumor_qcbias = os.path.join(tumor_directory, "picard_metrics", tumor_pair.tumor.name + ".qcbias_metrics.txt")
            input_tumor_all_picard = os.path.join(tumor_directory, "picard_metrics", tumor_pair.tumor.name + ".all.metrics.quality_distribution.pdf")
            input_tumor_qualimap = os.path.join(tumor_directory, "qualimap", tumor_pair.tumor.name, "genome_results.txt")

            [input_tumor_fastqc] = self.select_input_files([
                [os.path.join(tumor_directory, "fastqc", tumor_pair.tumor.name + ".sorted.dup_fastqc.zip")],
                [os.path.join(tumor_directory, "fastqc", tumor_pair.tumor.name + "_fastqc.zip")],
            ])

            input_dep += [
                input_normal_oxog,
                input_normal_qcbias,
                input_normal_all_picard,
                input_normal_qualimap,
                input_normal_fastqc,
                input_tumor_oxog,
                input_tumor_qcbias,
                input_tumor_all_picard,
                input_tumor_qualimap,
                input_tumor_fastqc
            ]

            output = os.path.join(metrics_directory, tumor_pair.name + ".multiqc")

            jobs.append(
                concat_jobs([
                    multiqc.run(
                        input_dep,
                        output
                        )
            ], name="multiqc." + tumor_pair.name))

        return jobs

    def sym_link_report(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] = [os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.name + ".multiqc.html")]

            for key, input_files in inputs.items():
                for idx, report_file in enumerate(input_files):
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(
                            report_file,
                            tumor_pair,
                            self.output_dir,
                            type="metrics",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_fastq.report." + str(idx) + "." + tumor_pair.name + "." + key))

        return jobs

    def rawmpileup(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                bed_file = coverage_bed

            input_normal = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
            ])

            input_tumor = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
            ])

            nb_jobs = config.param('rawmpileup', 'nb_jobs', param_type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
            
            if nb_jobs == 1:
                pair_output = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            samtools.mpileup(
                                [
                                    input_normal[0],
                                    input_tumor[0]
                                ],
                                pair_output,
                                config.param('rawmpileup', 'mpileup_other_options'),
                                regionFile=bed_file
                            )
                        ],
                        name="rawmpileup." + tumor_pair.name
                    )
                )

            else:
                
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        varscan_directory,
                                        remove=True
                                    ),
                                    samtools.mpileup(
                                        [
                                            input_normal[0],
                                            input_tumor[0]
                                        ],
                                        pair_output,
                                        config.param('rawmpileup', 'mpileup_other_options'),
                                        region=sequence['name'],
                                        regionFile=bed_file
                                    )
                                ],
                                name="rawmpileup." + tumor_pair.name + "." + sequence['name']
                            )
                        )

        return jobs

    def paired_varscan2(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data. 
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
        SSC INFO field remove to prevent collison with Samtools output during ensemble                     
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")
            output = os.path.join(varscan_directory, tumor_pair.name)

            nb_jobs = config.param('rawmpileup', 'nb_jobs', param_type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")
    
                output_snp = os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf")
                output_vcf = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")
    
                jobs.append(concat_jobs([
                    bash.mkdir(
                        varscan_directory,
                        remove=True
                    ),
                    varscan.somatic(
                        input_pair,
                        output,
                        config.param('varscan2_somatic', 'other_options'),
                        output_vcf_dep=output_vcf,
                        output_snp_dep=output_snp,
                        output_indel_dep=output_indel
                    ),
                    htslib.bgzip_tabix(
                        output_snp,
                        os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz")
                    ),
                    htslib.bgzip_tabix(
                        output_indel,
                        os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
                    ),
                    pipe_jobs([
                        bcftools.concat(
                            [os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz"),
                             os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")],
                            None
                        ),
                        Job(
                            [None],
                            [output_vcf],
                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                    + tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > "
                                    + output_vcf
                        ),
                    ]),
                    htslib.bgzip_tabix(
                        output_vcf,
                        output_vcf_gz
                    ),
                ], name="varscan2_somatic." + tumor_pair.name))

            else:

                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                        output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                        output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                        output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")
                        output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")

                        jobs.append(concat_jobs([
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            varscan.somatic(
                                input_pair,
                                output,
                                config.param('varscan2_somatic', 'other_options'),
                                output_vcf_dep=output_vcf,
                                output_snp_dep=output_snp,
                                output_indel_dep=output_indel
                            ),
                            htslib.bgzip_tabix(
                                output_snp,
                                os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz")
                            ),
                            htslib.bgzip_tabix(
                                output_indel,
                                os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")
                            ),
                            pipe_jobs([
                                bcftools.concat(
                                    [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz"),
                                     os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")],
                                    None
                                ),
                                Job(
                                    [None],
                                    [output_vcf],
                                    command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                            + tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > "
                                            + output_vcf
                                ),
                            ]),
                            htslib.bgzip_tabix(
                                output_vcf,
                                output_vcf_gz
                            ),
                        ], name="varscan2_somatic." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_varscan2(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup', 'nb_jobs', param_type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            all_inputs = []
            if nb_jobs == 1:
                all_inputs = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")

            else:
                all_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")
                              for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary']

            for input_vcf in all_inputs:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete varscan2 vcf: %s\n" % input_vcf)

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            all_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vt.vcf.gz")

            somtic_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            germline_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.view(
                            all_inputs,
                            None
                        ),
                        tools.fix_varscan_output(
                            None,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        #vt.sort("-", all_output, "-m full"),
                        htslib.bgzip_tabix(
                            None,
                            all_output
                        ),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            all_output,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            all_output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            all_output_vt,
                            None,
                            config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            somtic_output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            all_output_vt,
                            None,
                            config.param('varscan2_readcount_fpfilter', 'germline_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            germline_output_vt
                        ),
                    ]),
            	], name="merge_varscan2." + tumor_pair.name))

            else:
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            all_inputs,
                            None
                        ),
                        tools.fix_varscan_output(
                            None,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                        ),
                        #vt.sort("-", all_output, "-m full"),
                        htslib.bgzip_tabix(
                            None,
                            all_output
                        ),
                ]),
                #htslib.tabix(all_output),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(
                        all_output,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        all_output_vt
                    ),
                ]),
                pipe_jobs([
                    bcftools.view(
                        all_output_vt,
                        None,
                        config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')
                    ),
                    htslib.bgzip_tabix(
                        None,
                        somtic_output_vt
                    ),
                ]),
                pipe_jobs([
                    bcftools.view(
                        all_output_vt,
                        None,
                        config.param('varscan2_readcount_fpfilter', 'germline_filter_options')
                    ),
                    htslib.bgzip_tabix(
                        None,
                        germline_output_vt
                    ),
                ]),
            	], name="merge_varscan2." + tumor_pair.name))
             
        return jobs

    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        """

        jobs = []

        created_interval_lists = []
        
        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")

            input_normal = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
            ])

            input_tumor = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
            ])

            interval_list = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])
            if coverage_bed:
                interval_list = os.path.join(mutect_directory, re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed)))

                if not interval_list in created_interval_lists:
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(mutect_directory),
                                tools.bed2interval_list(
                                    coverage_bed,
                                    interval_list
                                )
                            ],
                            name="interval_list." + os.path.basename(coverage_bed)
                        )
                    )
                    created_interval_lists.append(interval_list)

            if nb_jobs == 1:

                jobs.append(
                    concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                mutect_directory,
                                remove=True
                            ),
                            gatk4.mutect2(
                                input_normal[0],
                                tumor_pair.normal.name,
                                input_tumor[0],
                                tumor_pair.tumor.name,
                                os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"),
                                os.path.join(mutect_directory, tumor_pair.name + ".f1r2.tar.gz"),
                                interval_list=interval_list
                            )
                        ],
                        name="gatk_mutect2." + tumor_pair.name
                    )
                )

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):

                    outprefix = tumor_pair.name + "." + str(idx) + ".mutect2"
                    jobs.append(
                        concat_jobs(
                            [
                                # Create output directory since it is not done by default by GATK tools
                                bash.mkdir(
                                    mutect_directory,
                                    remove=True
                                ),
                                gatk4.mutect2(
                                    input_normal[0],
                                    tumor_pair.normal.name,
                                    input_tumor[0],
                                    tumor_pair.tumor.name,
                                    os.path.join(mutect_directory, outprefix + ".vcf.gz"),
                                    os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".f1r2.tar.gz"),
                                    intervals=sequences,
                                    interval_list=interval_list
                                )
                            ],
                            name="gatk_mutect2." + tumor_pair.name + "." + str(idx)
                        )
                    )

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(
                    concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                mutect_directory,
                                remove=True
                            ),
                            gatk4.mutect2(
                                input_normal[0],
                                tumor_pair.normal.name,
                                input_tumor[0],
                                tumor_pair.tumor.name,
                                os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"),
                                os.path.join(mutect_directory, tumor_pair.name + ".others.f1r2.tar.gz"),
                                exclude_intervals=unique_sequences_per_job_others,
                                interval_list=interval_list
                            )
                        ],
                        name="gatk_mutect2." + tumor_pair.name + ".others"
                    )
                )

        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2
        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
        Generate a somatic vcf containing only PASS variants        
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_flt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.flt.vcf.gz")
            output_vt_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            output_somatic_vt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                if config.param('gatk_mutect2', 'module_gatk').split("/")[2] > "4":
                    jobs.append(concat_jobs([
                        Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                        gatk4.learn_read_orientation_model(
                            [os.path.join(mutect_directory, tumor_pair.name + ".f1r2.tar.gz")],
                            os.path.join(pair_directory, tumor_pair.name + ".f1r2.tar.gz")
                        ),
                        gatk4.filter_mutect_calls(
                            os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"),
                            output_flt,
                            read_orientation=os.path.join(pair_directory, tumor_pair.name + ".f1r2.tar.gz")
                        ),
                        pipe_jobs([
                            vt.decompose_and_normalize_mnps(
                                output_flt,
                                None
                            ),
                            Job(
                                [None],
                                [None],
                                command=" grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -vE 'EBV|hs37d5'"
                                        + " | sed -e 's#/\.##g'"
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_vt_gz
                            ),
                        ]),
                        pipe_jobs([
                            bcftools.view(
                                output_vt_gz,
                                None,
                                config.param('merge_filter_mutect2', 'filter_options')
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic_vt
                            ),
                        ]),
                    ], name="merge_filter_mutect2." + tumor_pair.name))
                
                else:
                    input_vcf = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                    jobs.append(concat_jobs([
                        Job(
                            [input_vcf],
                            [output_gz],
                            command="ln -s -f " + os.path.relpath(input_vcf, os.path.dirname(output_gz)) + " " + output_gz,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        ),
                        #gatk4.filter_mutect_calls(output_gz, output_flt),
                        pipe_jobs([
                            vt.decompose_and_normalize_mnps(
                                output_gz,
                                None
                            ),
                            Job(
                                [None],
                                [None],
                                command="sed 's/TUMOR/" + tumor_pair.tumor.name
                                        + "/g' | sed 's/NORMAL/"
                                        + tumor_pair.normal.name
                                        + "/g' | sed 's/Number=R/Number=./g' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -vE 'EBV|hs37d5'"
                                        + " | sed -e 's#/\.##g'"
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic_vt
                            ),
                        ]),
                    ], name="symlink_mutect_vcf." + tumor_pair.name))

            elif nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(
                    self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputs = []
                for idx, sequences in enumerate(unique_sequences_per_job):
                    inputs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))
                inputs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))

                for input_vcf in inputs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete mutect2 vcf: %s\n" % input_vcf)

                if config.param('gatk_mutect2', 'module_gatk').split("/")[2] > "4":

                    output_stats = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz.stats")
                    stats = []
                    for idx, sequences in enumerate(unique_sequences_per_job):
                        stats.append(
                            os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz.stats"))
                    stats.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz.stats"))

                    output_models = os.path.join(pair_directory, tumor_pair.name + ".read-orientation-model.tar.gz")
                    models = []
                    for idx, sequences in enumerate(unique_sequences_per_job):
                        models.append(
                            os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".f1r2.tar.gz"))
                    models.append(os.path.join(mutect_directory, tumor_pair.name + ".others.f1r2.tar.gz"))

                    jobs.append(concat_jobs([
                        Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                        gatk4.learn_read_orientation_model(
                            models,
                            output_models
                        ),
                        gatk4.cat_variants(
                            inputs,
                            output_gz
                        ),
                        gatk4.merge_stats(
                            stats,
                            output_stats
                        ),
                        gatk4.filter_mutect_calls(
                            output_gz,
                            output_flt,
                            read_orientation=output_models
                        ),
                        pipe_jobs([
                            vt.decompose_and_normalize_mnps(
                                output_flt,
                                None
                            ),
                            Job(
                                [None],
                                [None],
                                command=" grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -vE 'EBV|hs37d5'"
                                        + " | sed -e 's#/\.##g'"
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_vt_gz
                            ),
                        ]),
                        pipe_jobs([
                            bcftools.view(
                                output_vt_gz,
                                None,
                                config.param('merge_filter_mutect2', 'filter_options')
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic_vt
                            ),
                        ]),
                    ], name="merge_filter_mutect2." + tumor_pair.name))

                else:
                    jobs.append(concat_jobs([
                        Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                        pipe_jobs([
                            bcftools.concat(
                                inputs,
                                None,
                                config.param('merge_filter_mutect2', 'bcftools_options')
                            ),
                            Job(
                                [None],
                                [None],
                                command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/"
                                        + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                            ),

                            htslib.bgzip_tabix(
                                None,
                                output_gz
                            ),
                        ]),
                        #gatk4.filter_mutect_calls(output_gz, output_flt),
                        pipe_jobs([
                            vt.decompose_and_normalize_mnps(
                                output_gz,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_vt_gz
                            ),
                        ]),
                        pipe_jobs([
                            bcftools.view(
                                output_vt_gz,
                                None,
                                config.param('merge_filter_mutect2', 'filter_options')
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic_vt
                            ),
                        ]),
                    ], name="merge_filter_mutect2." + tumor_pair.name))

        return jobs

    def strelka2_paired_somatic(self):
        """
        Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs
        This implementation is optimized for somatic calling.
        [Strelka2](https://github.com/Illumina/strelka)
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if (tumor_pair.multiple_normal == 1):
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            somatic_dir = os.path.join(pair_directory, "rawStrelka2_somatic")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            input_normal = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
            ])

            input_tumor = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
            ])

            mantaIndels = None
            if os.path.isfile(os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, "rawManta", "results", "variants", "candidateSmallIndels.vcf.gz")):
                mantaIndels = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, "rawManta", "results", "variants", "candidateSmallIndels.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                local_coverage_bed = os.path.join(pair_directory, os.path.basename(coverage_bed))
                bed_file = local_coverage_bed + ".gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(pair_directory),
                            Job(
                                [coverage_bed],
                                [local_coverage_bed + ".sort"],
                                command="sort -V -k1,1 -k2,2n -k3,3n "
                                        + coverage_bed + " > "
                                        + local_coverage_bed + ".sort ; sleep 15"
                            ),
                            htslib.bgzip(
                                local_coverage_bed + ".sort",
                                bed_file
                            ),
                            htslib.tabix(
                                bed_file,
                                "-p bed"
                            )
                        ],
                        name="bed_index." + tumor_pair.name
                    )
                )

            else:
                bed_file=config.param('strelka2_paired_somatic', 'bed_file')

            output_dep = [
                os.path.join(somatic_dir, "results/variants/somatic.snvs.vcf.gz"),
                os.path.join(somatic_dir, "results/variants/somatic.indels.vcf.gz")
            ]

            jobs.append(
                concat_jobs(
                    [
                        bash.rm(somatic_dir),
                        strelka2.somatic_config(
                            input_normal[0],
                            input_tumor[0],
                            somatic_dir,
                            bed_file,
                            mantaIndels
                        ),
                        strelka2.run(
                            somatic_dir,
                            output_dep=output_dep
                        ),
                    ],
                    name="strelka2_paired_somatic.call."+tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    input_dependency=[input_normal[0], input_tumor[0]],
                    output_dependency=output_dep
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                bcftools.concat(
                                    output_dep,
                                    None
                                ),
                                Job(
                                    [None],
                                    [None],
                                    command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name
                                        + "/g' | sed 's/Number=R/Number=./g' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    output_prefix + ".strelka2.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.vt.vcf.gz"
                                )
                            ]
                        ),
                        tools.fix_genotypes_strelka(
                            output_prefix + ".strelka2.vt.vcf.gz",
                            output_prefix + ".strelka2.somatic.gt.vcf.gz",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name
                        ),
                        bcftools.view(
                            output_prefix + ".strelka2.somatic.gt.vcf.gz",
                            output_prefix + ".strelka2.somatic.vt.vcf.gz",
                            config.param('strelka2_paired_somatic', 'filter_options')
                        )
                    ],
                    name="strelka2_paired_somatic.filter." + tumor_pair.name
                )
            )

        return jobs

    def strelka2_paired_germline(self):
        """
        Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs
        This implementation is optimized for germline calling in cancer pairs.
        [Strelka2](https://github.com/Illumina/strelka)
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if (tumor_pair.multiple_normal == 1):
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            germline_dir = os.path.join(pair_directory, "rawStrelka2_germline")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            input_normal = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
            ])

            input_tumor = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
            ])
        
            input = [input_normal[0], input_tumor[0]]
        
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )
        
            if coverage_bed:
                local_coverage_bed = os.path.join(pair_directory, os.path.basename(coverage_bed))
                bed_file = local_coverage_bed + ".gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(pair_directory),
                            Job(
                                [coverage_bed],
                                [local_coverage_bed + ".sort"],
                                command="sort -V -k1,1 -k2,2n -k3,3n "
                                        + coverage_bed + " > "
                                        + local_coverage_bed + ".sort ; sleep 15"
                            ),
                            htslib.bgzip(
                                local_coverage_bed + ".sort",
                                bed_file
                            ),
                            htslib.tabix(
                                bed_file,
                                "-p bed"
                            )
                        ],
                        name="bed_index." + tumor_pair.name
                    )
                )
            
            else:
                bed_file = config.param('strelka2_paired_germline', 'bed_file')
                
            output_dep = [os.path.join(germline_dir, "results/variants/variants.vcf.gz")]
        
            jobs.append(
                concat_jobs(
                    [
                        bash.rm(germline_dir),
                        strelka2.germline_config(
                            input,
                            germline_dir,
                            bed_file,
                        ),
                        strelka2.run(
                            germline_dir,
                            output_dep=output_dep
                        )
                    ],
                    name="strelka2_paired_germline.call."+tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    input_dependency=input,
                    output_dependency=output_dep
                )
            )
        
            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                Job(
                                    [os.path.join(germline_dir, "results/variants/variants.vcf.gz")],
                                    [None],
                                    command="zcat " + os.path.join(germline_dir, "results/variants/variants.vcf.gz")
                                            + " | sed 's/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                            + " | sed 's/NORMAL/" + tumor_pair.normal.name
                                            + "/g' | sed 's/Number=R/Number=./g' | grep -vE 'GL00|hs37d5' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.germline.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    output_prefix + ".strelka2.germline.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.germline.gt.vcf.gz"
                                )
                            ]
                         ),
                        bcftools.view(
                            output_prefix + ".strelka2.germline.gt.vcf.gz",
                            output_prefix + ".strelka2.germline.vt.vcf.gz",
                            config.param('strelka2_paired_germline', 'filter_options')
                        )
                    ],
                    name="strelka2_paired_germline.filter." + tumor_pair.name
                )
            )
    
        return jobs

    def vardict_paired(self):
        """
        vardict caller for SNVs and Indels.
        Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is
        AUPAC nomenclature
        """

        ##TO DO - the BED system needs to be revisted !! 
        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', param_type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', param_type='filepath')

        interval_list = []

        splitjobs_dir = os.path.join(self.output_dirs['paired_variants_directory'], "splitjobs", "vardict" )
        if use_bed:
            for idx in range(nb_jobs):
                interval_list.append(
                    os.path.join(
                        splitjobs_dir,
                        "exome",
                        "interval_list",
                        str(idx).zfill(4) + "-scattered.interval_list"
                    )
                )

            jobs.append(concat_jobs([
                bash.mkdir(
                    os.path.join(splitjobs_dir,
                                 "exome",
                                 "interval_list"
                                 ),
                    remove=True
                ),
                gatk4.bed2interval_list(
                    genome_dictionary,
                    self.samples[0].readsets[0].beds[0],
                    os.path.join(splitjobs_dir,
                                 "exome",
                                 "interval_list",
                                 config.param('vardict_paired', 'assembly') + ".interval_list"
                                 )
                ),
                gatk4.splitInterval(
                    os.path.join(splitjobs_dir,
                                 "exome",
                                 "interval_list",
                                 config.param('vardict_paired', 'assembly') + ".interval_list"
                                 ),
                    os.path.join(splitjobs_dir, "exome", "interval_list"),
                    nb_jobs,
                    options="--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
                ),
                ], name="vardict_paired.create_splitjobs")
            )
        # else:
        #     for idx in range(nb_jobs):
        #         interval_list.append(os.path.join(splitjobs_dir,
        #                                           "wgs",
        #                                           "interval_list",
        #                                           str(idx).zfill(4) + "-scattered.interval_list"
        #                                           )
        #                              )
        #     jobs.append(concat_jobs([
        #         bash.mkdir(
        #             os.path.join(splitjobs_dir, "wgs", "interval_list"),
        #             remove=True
        #         ),
        #         picard2.scatterIntervalsByNs(
        #             config.param('vardict_paired', 'genome_fasta', type='filepath'),
        #             os.path.join(splitjobs_dir,
        #                          "wgs",
        #                          "interval_list",
        #                          config.param('vardict_paired', 'assembly') + ".interval_list"
        #                          ),
        #             options="OUTPUT_TYPE=ACGT"
        #         ),
        #         gatk4.splitInterval(
        #             os.path.join(splitjobs_dir,
        #                          "wgs",
        #                          "interval_list",
        #                          config.param('vardict_paired', 'assembly') + ".interval_list"
        #                          ),
        #             os.path.join(splitjobs_dir,
        #                          "wgs",
        #                          "interval_list"),
        #             nb_jobs,
        #         ),
        #         ], name="vardict_paired.create_splitjobs")
        #     )
            
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            
            input_normal = self.select_input_files(
                [[os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]])

            input_tumor = self.select_input_files(
                [[os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]])

            if use_bed:
                idx = 0
                for interval in interval_list:
                    bed = re.sub("interval_list$", "bed", interval)
                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx).zfill(4) + ".vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            vardict_directory,
                            remove=True
                        ),
                        gatk4.interval_list2bed(
                            interval,
                            bed
                        ),
                        pipe_jobs([
                        vardict.paired_java(
                            input_normal[0],
                            input_tumor[0],
                            tumor_pair.name,
                            None,
                            bed
                        ),
                        vardict.testsomatic(
                            None,
                            None
                        ),
                        vardict.var2vcf(
                            None,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output
                        ),
                        ]),
                    ],name="vardict_paired." + tumor_pair.name + "." + str(idx).zfill(4))
                    )
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))
            
                jobs.append(concat_jobs([
                    bash.mkdir(
                        vardict_directory,
                        remove=True
                    ),
                    vardict.dict2beds(
                        genome_dictionary,
                        beds
                    ),
                    ], name="vardict.genome.beds." + tumor_pair.name)
                )
                for idx in range(nb_jobs):
                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            vardict_directory,
                            remove=True
                        ),
                        pipe_jobs([
                            vardict.paired_java(
                                input_normal[0],
                                input_tumor[0],
                                tumor_pair.name,
                                None,
                                beds[idx]
                            ),
                            vardict.testsomatic(
                                None,
                                None
                            ),
                            vardict.var2vcf(
                                None,
                                tumor_pair.normal.name,
                                tumor_pair.tumor.name,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output
                            ),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + "." + str(idx))
                    )
                
        return jobs

    def merge_filter_paired_vardict(self):
        """
        The fully merged vcf is filtered using following steps:
        1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
        2. Somatics identified in step 1 must have PASS filter
        """

        jobs = []
        nb_jobs = config.param('vardict_paired', 'nb_jobs', param_type='posint')
        use_bed = config.param('vardict_paired', 'use_bed', param_type='boolean', required=True)

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output_tmp = os.path.join(pair_directory, tumor_pair.name + ".vardict.tmp.vcf.gz")
            output = os.path.join(pair_directory, tumor_pair.name + ".vardict.vcf.gz")
            output_vt = os.path.join(pair_directory, tumor_pair.name + ".vardict.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".vardict.germline.vt.vcf.gz")

            if nb_jobs == 1 and use_bed:
                inputs = os.path.join(vardict_directory, tumor_pair.name + "." + str(0).zfill(4) + ".vardict.vcf.gz")
                jobs.append(concat_jobs([
                    Job(
                        [inputs],
                        [output_tmp],
                        command="ln -s -f " + os.path.relpath(inputs, os.path.dirname(output_tmp)) + " " + output_tmp,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    ),
                    pipe_jobs([
                        Job(
                            [output_tmp],
                            [None],
                            command="zcat " + output_tmp + " | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"),
                        htslib.bgzip_tabix(
                            None,
                            output
                        )
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'somatic_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline_loh
                        ),
                    ]),
                ], name="symlink_vardict_vcf." + tumor_pair.name))
            else:
                inputVCFs = []
                for idx in range(nb_jobs):
                    inputVCFs.append(
                        os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))

                for input_vcf in inputVCFs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete vardict vcf: %s\n" % input_vcf)

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            inputVCFs,
                            None
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                        ),
                        Job(
                            [None],
                            [None],
                            command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output
                        ),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(
                            output,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_vt
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'somatic_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic
                        ),
                    ]),
                    pipe_jobs([
                        bcftools.view(
                            output_vt,
                            None,
                            config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline_loh
                        ),
                    ]),
                ], name="merge_filter_paired_vardict." + tumor_pair.name))

        return jobs

    def ensemble_somatic(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")
            input_strelka2 = os.path.join(input_directory, tumor_pair.name + ".strelka2.somatic.vt.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            inputs_somatic = [input_mutect2, input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_somatic:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete ensemble vcf: %s\n" % input_vcf)

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                bash.mkdir(
                    paired_ensemble_directory,
                    remove=True
                ),
                bcbio_variation_recall.ensemble(
                    inputs_somatic,
                    output_ensemble,
                    config.param('bcbio_ensemble_somatic', 'options')
                ),
            ], name="bcbio_ensemble_somatic." + tumor_pair.name))

        return jobs

    def ensemble_germline_loh(self):
        """
        Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            input_strelka2 = os.path.join(input_directory, tumor_pair.name + ".strelka2.germline.vt.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")
            
            inputs_germline = [input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_germline:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete ensemble vcf: %s\n" % input_vcf)

            output_ensemble = os.path.join(paired_ensemble_directory,
                                           tumor_pair.name + ".ensemble.germline.vt.vcf.gz")

            # if os.path.isdir(os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline.vt-work")):
            #     rm_job = bash.rm(
            #         os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline.vt-work")
            #     )
            #     jobs.append(rm_job)

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                bash.mkdir(
                    paired_ensemble_directory,
                    remove=True
                ),
                bcbio_variation_recall.ensemble(
                    inputs_germline,
                    output_ensemble,
                    config.param('bcbio_ensemble_germline', 'options')
                ),
            ], name="bcbio_ensemble_germline." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: Standard and Somatic annotations
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if (tumor_pair.multiple_normal == 1):
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            annot_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, "rawAnnotation")
            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join( tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
    
                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_somatic_variants,
                        output_somatic_variants,
                        config.param('gatk_variant_annotator_somatic', 'other_options')
                    ),
                ], name="gatk_variant_annotator_somatic." + tumor_pair.name))
                
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_somatic_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot." + str(idx) + ".vcf.gz")

                    jobs.append(concat_jobs([
                        bash.mkdir(
                            annot_directory,
                            remove=True
                        ),
                        gatk.variant_annotator(
                            input_normal,
                            input_tumor,
                            input_somatic_variants,
                            output_somatic_variants,
                            config.param('gatk_variant_annotator_somatic', 'other_options'),
                            intervals=sequences
                        ),
                    ], name="gatk_variant_annotator_somatic." + str(idx) + "." + tumor_pair.name))

                output_somatic_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.others.vcf.gz")

                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_somatic_variants,
                        output_somatic_variants,
                        config.param('gatk_variant_annotator_somatic', 'other_options'),
                        exclude_intervals=unique_sequences_per_job_others
                    ),
                ], name="gatk_variant_annotator_somatic.others." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if (tumor_pair.multiple_normal == 1):
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            annot_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, "rawAnnotation")
            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.vcf.gz")
    
            if nb_jobs == 1:
                output_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
        
                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_germline_variants,
                        output_germline_variants,
                        config.param('gatk_variant_annotator_germline', 'other_options'),
                    ),
                ], name="gatk_variant_annotator_germline." + tumor_pair.name))
    
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_germline_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot." + str(idx) + ".vcf.gz")
            
                    jobs.append(concat_jobs([
                        bash.mkdir(
                            annot_directory,
                            remove=True
                        ),
                        gatk.variant_annotator(
                            input_normal,
                            input_tumor,
                            input_germline_variants,
                            output_germline_variants,
                            config.param('gatk_variant_annotator_germline', 'other_options'),
                            intervals=sequences
                        ),
                    ], name="gatk_variant_annotator_germline." + str(idx) + "." + tumor_pair.name))
        
                output_germline_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot.others.vcf.gz")
        
                jobs.append(concat_jobs([
                    bash.mkdir(
                        annot_directory,
                        remove=True
                    ),
                    gatk.variant_annotator(
                        input_normal,
                        input_tumor,
                        input_germline_variants,
                        output_germline_variants,
                        config.param('gatk_variant_annotator_germline', 'other_options'),
                        exclude_intervals=unique_sequences_per_job_others
                    ),
                ], name="gatk_variant_annotator_germline.others." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_somatic(self):
        """
        Merge annotated somatic vcfs
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_somatic = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                vcfs_to_merge = [os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot." + str(idx) +".vcf.gz")
                                  for idx in range(len(unique_sequences_per_job))]
                
                vcfs_to_merge.append(os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.others.vcf.gz"))
                
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            vcfs_to_merge,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_somatic
                        ),
                    ]),
                ], name="merge_gatk_variant_annotator.somatic." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_germline(self):
        """
        Merge annotated germline and LOH vcfs
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_germline = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
            
            if nb_jobs > 1:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                vcfs_to_merge = [os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation", tumor_pair.name + ".ensemble.germline.vt.annot." + str(idx) + ".vcf.gz")
                                 for idx in range(len(unique_sequences_per_job))]

                vcfs_to_merge.append(os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot.others.vcf.gz"))
        
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(
                            vcfs_to_merge,
                            None
                        ),
                        htslib.bgzip_tabix(
                            None,
                            output_germline
                        ),
                    ]),
                ], name="merge_gatk_variant_annotator.germline." + tumor_pair.name))

        return jobs

    def filter_ensemble_germline(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields
        """
        jobs = []
    
        ensemble_directory = os.path.join(self.output_dir, "pairedVariants", "ensemble")
    
        for tumor_pair in self.tumor_pairs.values():
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz"
            )
            output_2caller = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.2caller.vcf.gz"
            )
            output_filter = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.2caller.flt.vcf.gz"
            )
            
            jobs.append(concat_jobs([
                tools.format2pcgr(
                    input,
                    output_2caller,
                    config.param('filter_ensemble', 'call_filter'),
                    "germline",
                    tumor_pair.tumor.name,
                ),
                pipe_jobs([
                    bcftools.view(
                        output_2caller,
                        None,
                        filter_options=config.param('filter_ensemble', 'filter_options'),
                    ),
                    bcftools.view(
                        None,
                        None,
                        filter_options="-Ov -s ^" + tumor_pair.normal.name
                    ),
                    bcftools.sort(
                        None,
                        output_filter,
                        sort_options="-Oz"
                    ),
                ]),
                htslib.tabix(
                    output_filter,
                    options="-pvcf"
                ),
                ], name="filter_ensemble.germline." + tumor_pair.name))
    
        return jobs

    def report_cpsr(self):
        """
        Creates a cpsr gremline report (https://sigven.github.io/cpsr/)
        input: filtered ensemble gremline vcf
        output: html report and addtionalflat files
        """
        jobs = []
    
        ensemble_directory = os.path.join(self.output_dir, "pairedVariants", "ensemble")
    
        for tumor_pair in self.tumor_pairs.values():
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.2caller.flt.vcf.gz"
            )
            cpsr_directory = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                "cpsr"
            )
        
            jobs.append(concat_jobs([
                bash.mkdir(
                    cpsr_directory,
                ),
                cpsr.report(
                    input,
                    cpsr_directory,
                    tumor_pair.name
                )
            ], name="report_cpsr." + tumor_pair.name))
    
        return jobs

    def filter_ensemble_somatic(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields
        """
        jobs = []
    
        ensemble_directory = os.path.join(self.output_dir, "pairedVariants", "ensemble")
    
        for tumor_pair in self.tumor_pairs.values():
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz"
            )
            output_2caller = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.vcf.gz"
            )
            output_filter = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
            )
        
            jobs.append(concat_jobs([
                tools.format2pcgr(
                    input,
                    output_2caller,
                    config.param('filter_ensemble', 'call_filter'),
                    "somatic",
                    tumor_pair.tumor.name,
                ),
                bcftools.view(
                    output_2caller,
                    output_filter,
                    filter_options=config.param('filter_ensemble', 'filter_options'),
                ),
                htslib.tabix(
                    output_filter,
                    options="-pvcf"
                ),
            ], name="filter_ensemble.somatic." + tumor_pair.name))
    
        return jobs

    def report_pcgr(self):
        """
        Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
        input: filtered ensemble gremline vcf
        output: html report and addtionalflat files
        """
        jobs = []
    
        ensemble_directory = os.path.join(
            self.output_dir,
            "pairedVariants",
            "ensemble"
        )
        assembly = config.param('report_pcgr', 'assembly')

        for tumor_pair in self.tumor_pairs.values():
            cpsr_directory = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                "cpsr"
            )
            input_cpsr = os.path.join(
                cpsr_directory,
                tumor_pair.name + ".cpsr." + assembly + ".json.gz"
            )
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
            )
            input_cna = os.path.join(
                self.output_dir,
                "SVariants",
                tumor_pair.name,
                tumor_pair.name + ".cnvkit.vcf.gz"
            )
            header = os.path.join(
                self.output_dir,
                "SVariants",
                "header"
            )
            output_cna_body = os.path.join(
                self.output_dir,
                "SVariants",
                tumor_pair.name + ".cnvkit.body.tsv"
            )
            output_cna = os.path.join(
                self.output_dir,
                "SVariants",
                tumor_pair.name + ".cnvkit.cna.tsv"
            )
            pcgr_directory = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                "pcgr"
            )
            output = os.path.join(
                pcgr_directory,
                tumor_pair.name + ".pcgr_acmg." + assembly + ".flexdb.html"
            )
     
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            pcgr_directory,
                        ),
                        pcgr.create_header(
                            header,
                        ),
                        bcftools.query(
                            input_cna,
                            output_cna_body,
                            query_options="-f '%CHROM\\t%POS\\t%END\\t%FOLD_CHANGE_LOG\\n'"
                        ),
                        bash.cat(
                            [
                                header,
                                output_cna_body,
                            ],
                            output_cna,
                        ),
                        pcgr.report(
                            input,
                            output_cna,
                            input_cpsr,
                            pcgr_directory,
                            tumor_pair.tumor.name
                        )
                    ],
                    name="report_pcgr." + tumor_pair.name,
                    input_dependency = [header, input, input_cna, input_cpsr, output_cna_body],
                    output_dependency = [header, output_cna_body, output_cna, output]
                )
            )
        
        return jobs

    def compute_cancer_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        if not os.path.exists(ensemble_directory):
            os.makedirs(ensemble_directory)

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            if not os.path.exists(paired_directory):
                os.makedirs(paired_directory)

            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                snpeff.compute_effects(
                    input_somatic,
                    output_somatic,
                    cancer_sample_file=cancer_pair_filename,
                                       options=config.param('compute_cancer_effects_somatic', 'options')
                ),
                htslib.bgzip_tabix(
                    output_somatic,
                    output_somatic + ".gz"
                ),
            ], name="compute_cancer_effects_somatic." + tumor_pair.name))

        return jobs

    def compute_cancer_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)

            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
            output_germline = os.path.join(paired_directory,
                                           tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.vcf")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                snpeff.compute_effects(
                    input_germline,
                    output_germline,
                    options=config.param('compute_cancer_effects_germline', 'options')
                ),
                htslib.bgzip_tabix(
                    output_germline,
                    output_germline + ".gz"
                ),
            ], name="compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def ensemble_somatic_dbnsfp_annotation(self):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        """
    
        jobs = []
    

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")
            output_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.dbnsfp.vcf")
            
            jobs.append(concat_jobs([
                snpeff.snpsift_dbnsfp(
                    input_vcf,
                    output_vcf
                ),
                htslib.bgzip_tabix(
                    output_vcf,
                    output_vcf + ".gz"
                ),
            ], name="dbnsfp_annotation.somatic." + tumor_pair.name))
        # job.samples = self.samples
    
        return jobs

    def ensemble_germline_dbnsfp_annotation(self):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        """
    
        jobs = []
    
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
    
        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.vcf.gz")
            output_vcf = os.path.join(paired_directory,
                                      tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.dbnsfp.vcf")
        
            jobs.append(concat_jobs([
                snpeff.snpsift_dbnsfp(
                    input_vcf,
                    output_vcf
                ),
                htslib.bgzip_tabix(
                    output_vcf,
                    output_vcf + ".gz"
                ),
            ], name="dbnsfp_annotation.germline." + tumor_pair.name))
        # job.samples = self.samples
    
        return jobs

    def sample_gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database :
        [Gemini](http://gemini.readthedocs.org/en/latest/index.html)
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)


            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                gemini.gemini_annotations(
                    gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                    gemini_prefix + ".somatic.gemini." + gemini_version + ".db",
                    self.output_dir
                )
            ], name="gemini_annotations.somatic." + tumor_pair.name))

        return jobs

    def sample_gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database :
        [Gemini](http://gemini.readthedocs.org/en/latest/index.html)
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                bash.mkdir(
                    paired_directory,
                    remove=True
                ),
                gemini.gemini_annotations(
                    gemini_prefix + ".ensemble.germline.vt.annot.snpeff.vcf.gz",
                    gemini_prefix + ".germline.gemini." + gemini_version + ".db",
                    self.output_dir
                )
            ], name="gemini_annotations.germline." + tumor_pair.name))

        return jobs

    def sym_link_ensemble(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] =  [os.path.join(self.output_dirs["paired_variants_directory"], "ensemble", tumor_pair.name, tumor_pair.name)]

            for key, input_files in inputs.items():
                for idx, sample_prefix in enumerate(input_files):
                    jobs.append(concat_jobs([
                        deliverables.md5sum(
                            sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz",
                            sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz.md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz.md5",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.md5sum(
                            sample_prefix + ".ensemble.germline.vt.annot.vcf.gz",
                            sample_prefix + ".ensemble.germline.vt.annot.vcf.gz.md5",
                            self.output_dir
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".ensemble.germline.vt.annot.vcf.gz.md5",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".ensemble.germline.vt.annot.vcf.gz",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                        deliverables.sym_link_pair(
                            sample_prefix + ".ensemble.germline.vt.annot.vcf.gz.tbi",
                            tumor_pair,
                            self.output_dir,
                            type="snv/ensemble",
                            sample=key,
                            profyle=self.args.profyle
                        ),
                    ], name="sym_link_ensemble." + str(idx) + "." + tumor_pair.name + "." + key))

        return jobs

    def combine_tumor_pairs_somatic(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.values()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                Job(
                    [input_merged_vcfs[0]],
                    [output],
                    command="ln -s -f " + os.path.relpath(input_merged_vcfs[0], os.path.dirname(output)) + " " + output
                )
            ], name="gatk_combine_variants.somatic.allPairs"))

        else:

            jobs.append(concat_jobs([
                bash.mkdir(
                    ensemble_directory,
                    remove=True
                ),
                gatk.combine_variants(
                    input_merged_vcfs,
                    output
                )
            ], name="gatk_combine_variants.somatic.allPairs"))

        return jobs

    def combine_tumor_pairs_germline(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name,
                                          tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz") for tumor_pair in
                             self.tumor_pairs.values()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        if len(input_merged_vcfs) == 1:
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        Job(
                            [input_merged_vcfs[0]],
                            [output],
                            command="ln -s -f " + os.path.relpath(input_merged_vcfs[0], os.path.dirname(output)) + " " + output
                        )
                    ],
                    name="gatk_combine_variants.germline.allPairs",
                    samples=sample_list
                )
            )

        else:
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        gatk.combine_variants(
                            input_merged_vcfs,
                            output
                        )
                    ],
                    name="gatk_combine_variants.germline.allPairs",
                    samples=sample_list
                )
            )

        return jobs

    def decompose_and_normalize_mnps_somatic(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    vt.decompose_and_normalize_mnps(
                        input,
                        output
                    )
                ],
                name="decompose_and_normalize_mnps.somatic.allPairs",
                samples=sample_list
            )
        )

        return jobs

    def decompose_and_normalize_mnps_germline(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline.annot.vcf.gz")
        output_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input_vcf, output_vcf)
        job.name = "decompose_and_normalize_mnps.germline.allPairs"

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    vt.decompose_and_normalize_mnps(
                        input_vcf,
                        output_vcf
                    )
                ],
                name="decompose_and_normalize_mnps.somatic.allPairs",
                samples=sample_list
            )
        )

        return jobs

    def all_pairs_compute_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        cancer_pair = open(cancer_pair_filename, 'w')

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    snpeff.compute_effects(
                        input,
                        output,
                        cancer_sample_file=cancer_pair_filename,
                        options=config.param('compute_cancer_effects_somatic', 'options')
                    ),
                    htslib.bgzip_tabix(
                        output,
                        output_gz
                    ),
                ],
                name="compute_effects.somatic.allPairs",
                samples=sample_list
            )
        )

        return jobs

    def all_pairs_compute_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.snpeff.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    snpeff.compute_effects(
                        input,
                        output,
                        options=config.param('compute_cancer_effects_germline', 'options')
                    ),
                    htslib.bgzip_tabix(
                        output,
                        output_gz
                    ),
                ],
                name="compute_effects.germline.allPair",
                samples=sample_list
            )
        )

        return jobs

    def gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    gemini.gemini_annotations(
                        gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                        gemini_prefix + ".somatic.gemini.db",
                        temp_dir
                    )
                ],
                name="gemini_annotations.somatic.allPairs",
                samples=sample_list
            )
        )

        return jobs

    def gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    gemini.gemini_annotations(
                        gemini_prefix + ".ensemble.germline.vt.annot.snpeff.vcf.gz",
                        gemini_prefix + ".germline.gemini.db",
                        temp_dir
                    )
                ],
                name="gemini_annotations.germline.allPairs",
                samples=sample_list
            )
        )

        return jobs

    def sequenza(self):
        """
        Sequenza is a novel set of tools providing a fast python script to genotype cancer samples,
        and an R package to estimate cancer cellularity, ploidy, genome wide copy number profile and infer
        for mutated alleles.

        """
        jobs = []
        nb_jobs = config.param('sequenza', 'nb_jobs', param_type='posint')
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            sequenza_directory = os.path.join(pair_directory, "sequenza")
            rawSequenza_directory = os.path.join(sequenza_directory, "rawSequenza")
            
            inputNormal = self.select_input_files(
                [[os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]])

            inputTumor = self.select_input_files(
                [[os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]])

            rawOutput = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".")
            output = os.path.join(sequenza_directory, tumor_pair.name + ".")
            
            if nb_jobs == 1:
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            sequenza.bam2seqz(
                                inputNormal[0],
                                inputTumor[0],
                                config.param('sequenza', 'gc_file'),
                                rawOutput + "all.seqz.gz",
                                None
                            ),
                            sequenza.bin(
                                rawOutput + "all.seqz.gz",
                                output + "all.binned.seqz.gz",
                            )
                        ],
                        name="sequenza.create_seqz." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )
                
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            sequenza.main(
                                output + "all.binned.seqz.gz",
                                sequenza_directory,
                                tumor_pair.name
                            ),
                            # sequenza.filter(
                            #     os.path.join(sequenza_directory, tumor_pair.name + "_segments.txt"),
                            #     tumor_pair.name, os.path.join(sequenza_directory, tumor_pair.name + ".segments.txt")
                            # ),
                            # sequenza.annotate(
                            #     os.path.join(sequenza_directory, tumor_pair.name + ".segments.txt"),
                            #     os.path.join(sequenza_directory, tumor_pair.name + ".annotated"),
                            #     os.path.join(sequenza_directory, tumor_pair.name + ".tmp")
                            # )
                        ],
                        name="sequenza." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        
                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        rawSequenza_directory,
                                        remove=True
                                    ),
                                    sequenza.bam2seqz(
                                        inputNormal[0],
                                        inputTumor[0],
                                        config.param('sequenza', 'gc_file'),
                                        rawOutput + "seqz." + sequence['name'] + ".gz",
                                        sequence['name']
                                    ),
                                    sequenza.bin(
                                        rawOutput + "seqz." + sequence['name'] + ".gz",
                                        rawOutput + "binned.seqz." + sequence['name'] + ".gz",
                                    ),
                                ],
                                name="sequenza.create_seqz." + sequence['name'] + "." + tumor_pair.name,
                                samples=[tumor_pair.normal, tumor_pair.tumor]
                            )
                        )

                seqz_outputs = [rawOutput + "binned.seqz." + sequence['name'] + ".gz" for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary']

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            Job(
                                seqz_outputs,
                                [output + "binned.merged.seqz.gz"],
                                command = "zcat "
                                        + " \\\n".join(seqz_outputs)
                                        + " \\\n | gawk 'FNR==1 && NR==1{print;}{ if($1!=\"chromosome\" && $1!=\"MT\" && $1!=\"chrMT\" && $1!=\"chrM\") {print $0} }' | "
                                        + " \\\n gzip -cf > "
                                        + output + "binned.merged.seqz.gz"
                                )
                        ],
                        name="sequenza.merge_binned_seqz." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )
    
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                rawSequenza_directory,
                                remove=True
                            ),
                            sequenza.main(
                                output + "binned.merged.seqz.gz",
                                sequenza_directory,
                                tumor_pair.name
                            )
                            #sequenza.filter(
                            #    os.path.join(sequenza_directory, tumor_pair.name + "_segments.txt"),
                            #    tumor_pair.name, os.path.join(sequenza_directory, tumor_pair.name + ".segments.txt")
                            #),
                            #sequenza.annotate(
                            #    os.path.join(sequenza_directory, tumor_pair.name + ".segments.txt"),
                            #    os.path.join(sequenza_directory, tumor_pair.name + ".annotated"),
                            #    os.path.join(sequenza_directory, tumor_pair.name + ".tmp")
                            #)
                        ],
                        name="sequenza." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )

        return jobs

    def sym_link_sequenza(self):
        """
        Sym link of sequenza outputs
        """
        jobs = []

        inputs = dict()

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            inputs["Tumor"] = [
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_chromosome_view.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_genome_view.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CN_bars.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CP_contours.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_ploidy_celularity.tsv")
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv/cnv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link.sequenza." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs

    def purple(self):
        """
        PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

        It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
        somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_dir = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            purple_dir = os.path.join(pair_dir, "purple")
            amber_dir = os.path.join(purple_dir, "rawAmber")
            cobalt_dir = os.path.join(purple_dir, "rawCobalt")
        
            inputNormal = self.select_input_files(
                [[os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]])
        
            inputTumor = self.select_input_files(
                [[os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]])

            somatic_snv = None
            if os.path.join(pair_dir, tumor_pair.name + ".strelka2.somatic.vt.vcf.gz"):
                somatic_snv = os.path.join(pair_dir, tumor_pair.name + ".strelka2.somatic.purple.vcf.gz")
                jobs.append(
                    concat_jobs(
                        [
                        purple.strelka2_convert(
                            os.path.join(pair_dir, tumor_pair.name + ".strelka2.somatic.vt.vcf.gz"),
                            somatic_snv,
                        )
                    ],
                    name="purple.convert_strelka2." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            amber_dir,
                            remove=True
                        ),
                        amber.run(
                            inputNormal[0],
                            inputTumor[0],
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            amber_dir,
                        )
                    ],
                    name="purple.amber." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [                
                        bash.mkdir(
                            cobalt_dir,
                            remove=True
                        ),
                        cobalt.run(
                            inputNormal[0],
                            inputTumor[0],
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            cobalt_dir,
                        ),
                    ],
                    name="purple.cobalt." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        purple.run(
                            amber_dir,
                            cobalt_dir,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            purple_dir,
                            somatic_snv,
                        )  
                    ],
                    name="purple.purity." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )
            
        return jobs

    def delly_call_filter(self):
        """
        Delly2 is an integrated structural variant prediction method that can
        discover, genotype and visualize deletions, tandem duplications, inversions and translocations
        at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends
        and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.
        Structural variants can be visualized using Delly-maze and Delly-suave.
        input: normal and tumor final bams
        Returns:bcf file

        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")

            filename = os.path.join(delly_directory, tumor_pair.name + '.tsv')
            if not os.path.exists(os.path.dirname(filename)):
                os.makedirs(os.path.dirname(filename))
           
            cancer_pair = open(filename, 'w')
            cancer_pair.write(tumor_pair.tumor.name + "\ttumor\n")
            cancer_pair.write(tumor_pair.normal.name + "\tcontrol\n")

            inputNormal = os.path.join(normal_alignment_directory,
                                       tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join(tumor_alignment_directory,
                                      tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            inputs = [inputTumor, inputNormal]

            SV_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            for sv_type in SV_types:
                output_bcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf")
                output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".somatic.flt.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                delly_directory,
                                remove=True
                            ),
                            delly.call(
                                inputs,
                                output_bcf,
                                sv_type
                            ),
                            pipe_jobs([
                                bcftools.view(
                                    output_bcf,
                                    None,
                                    config.param('delly_call_filter_somatic', 'bcftools_options')
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_vcf
                                ),
                            ]),
                        ],
                        name="delly_call_filter." + str(sv_type) + "." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )

        return jobs

    def delly_sv_annotation(self):
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            final_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name, tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")
            output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly.merge.sort.vcf.gz")
            output_flt_vcf = os.path.join(pair_directory, tumor_pair.name + ".delly.merge.sort.flt.vcf.gz")
            
            SV_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            inputBCF = []
            for sv_type in SV_types:
                inputBCF.append(os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf"))

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                bcftools.concat(
                                    inputBCF,
                                    None,
                                    "-O v"
                                ),
                                vt.sort(
                                    "-",
                                    "-",
                                    "-m full"
                                ),
                                htslib.bgzip(
                                    None,
                                    output_vcf
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                bcftools.view(
                                    output_vcf,
                                    None,
                                    "-f PASS"
                                ),
                                htslib.bgzip(
                                    None,
                                    output_flt_vcf
                                )
                            ]
                        ),
                    ],
                    name="sv_annotation.delly.merge_sort_filter." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        vawk.paired_somatic(
                            output_flt_vcf,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            final_directory + ".delly.somatic.vcf"
                        ),
                        htslib.bgzip(
                            final_directory + ".delly.somatic.vcf",
                            final_directory + ".delly.somatic.vcf.gz"
                        ),
                        snpeff.compute_effects(
                            final_directory + ".delly.somatic.vcf",
                            final_directory + ".delly.somatic.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            final_directory + ".delly.somatic.snpeff.vcf",
                            final_directory + ".delly.somatic.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            final_directory + ".delly.somatic.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "DELLY",
                            final_directory + ".delly.somatic.prioritize.tsv"
                        ),
                    ],
                    name="sv_annotation.delly.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        vawk.paired_germline(
                            output_flt_vcf,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            final_directory + ".delly.germline.vcf"
                        ),
                        htslib.bgzip(
                            final_directory + ".delly.germline.vcf",
                            final_directory + ".delly.germline.vcf.gz"
                        ),
                        snpeff.compute_effects(
                            final_directory + ".delly.germline.vcf",
                            final_directory + ".delly.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            final_directory + ".delly.germline.snpeff.vcf",
                            final_directory + ".delly.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            final_directory + ".delly.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "DELLY",
                            final_directory + ".delly.germline.prioritize.tsv"
                        ),
                    ],
                    name="sv_annotation.delly.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )
            
        return jobs
    
    def sym_link_delly(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".delly.somatic.snpeff.annot.vcf",
                pair_directory + ".delly.somatic.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_delly.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".delly.germline.snpeff.annot.vcf",
                pair_directory + ".delly.germline.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_delly.germline." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs
    
        
    def manta_sv_calls(self):
        """
        Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
        analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
        Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
        single efficient workflow.
        Returns:Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
         in VCF 4.1 format.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            manta_directory = os.path.join(pair_directory, "rawManta")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            inputNormal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            manta_somatic_output = os.path.join(manta_directory, "results/variants/somaticSV.vcf.gz")
            manta_germline_output = os.path.join(manta_directory, "results/variants/diploidSV.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                local_coverage_bed = os.path.join(pair_directory, os.path.basename(coverage_bed))
                bed_file = local_coverage_bed + ".gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(manta_directory),
                            Job(
                                [coverage_bed],
                                [local_coverage_bed + ".sort"],
                                command="sort -V -k1,1 -k2,2n -k3,3n " + coverage_bed + " | sed 's#chr##g' > "
                                        + local_coverage_bed + ".sort"
                            ),
                            htslib.bgzip(
                                local_coverage_bed + ".sort",
                                bed_file
                            ),
                            htslib.tabix(
                                bed_file,
                                "-p bed"
                            )
                        ],
                        name="bed_index." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )

            output_dep = [
                manta_somatic_output,
                manta_somatic_output + ".tbi",
                manta_germline_output,
                manta_germline_output + ".tbi"
            ]

            jobs.append(
                concat_jobs(
                    [
                        bash.rm(manta_directory),
                        bash.mkdir(
                            manta_directory,
                            remove=True
                         ),
                         manta.manta_config(
                            inputNormal,
                            inputTumor,
                            manta_directory,
                            bed_file
                        ),
                        manta.manta_run(
                            manta_directory,
                            output_dep=output_dep
                        ),
                        bash.ln(
                            manta_somatic_output,
                            output_prefix + ".manta.somatic.vcf.gz",
                            self.output_dir,
                        ),
                        bash.ln(
                            manta_somatic_output + ".tbi",
                            output_prefix + ".manta.somatic.vcf.gz.tbi",
                            self.output_dir
                        ),
                        bash.ln(
                            manta_germline_output,
                            output_prefix + ".manta.germline.vcf.gz",
                            self.output_dir,
                        ),
                        bash.ln(
                            manta_germline_output + ".tbi",
                            output_prefix + ".manta.germline.vcf.gz.tbi",
                            self.output_dir,
                        ),
                    ],
                    name="manta_sv." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def manta_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(
                    pair_directory + ".manta.somatic.vcf.gz",
                    pair_directory + ".manta.somatic.snpeff.vcf"
                ),
                annotations.structural_variants(
                    pair_directory + ".manta.somatic.snpeff.vcf",
                    pair_directory + ".manta.somatic.snpeff.annot.vcf"
                ),
                vawk.sv(
                    pair_directory + ".manta.somatic.snpeff.annot.vcf",
                    tumor_pair.normal.name,
                    tumor_pair.tumor.name,
                    "MANTA",
                    pair_directory + ".manta.somatic.prioritize.tsv"
                ),
            ], name="sv_annotation.manta_somatic." + tumor_pair.name))

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            pair_directory + ".manta.germline.vcf.gz",
                            pair_directory + ".manta.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            pair_directory + ".manta.germline.snpeff.vcf",
                            pair_directory + ".manta.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            pair_directory + ".manta.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "MANTA",
                            pair_directory + ".manta.germline.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.manta_germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def sym_link_manta(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".manta.somatic.snpeff.annot.vcf",
                pair_directory + ".manta.somatic.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_manta.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".manta.germline.snpeff.annot.vcf",
                pair_directory + ".manta.germline.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_manta.germline." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs

    def lumpy_paired_sv(self):
        """
        A probabilistic framework for structural variant discovery.
        Lumpy traditional with paired ends and split reads on tumor normal pair.
        Returns:bams.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name)
            lumpy_directory = os.path.join(pair_directory, "rawLumpy")
            inputNormal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            discordants_normal = os.path.join(lumpy_directory, normal_alignment_directory + ".discordants.sorted.bam")
            discordants_tumor = os.path.join(lumpy_directory, tumor_alignment_directory + ".discordants.sorted.bam")

            splitters_tumor = os.path.join(lumpy_directory, normal_alignment_directory + ".splitters.sorted.bam")
            splitters_normal = os.path.join(lumpy_directory, tumor_alignment_directory + ".splitters.sorted.bam")

            output_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf")
            gzip_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf.gz")

            genotype_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf")
            genotype_gzip = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            lumpy_directory,
                            remove=True
                        ),
                        pipe_jobs(
                            [
                                samtools.view(
                                    inputNormal,
                                    None,
                                    "-b -F 1294"
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    discordants_normal,
                                    lumpy_directory,
                                    config.param('extract_discordant_reads', 'sambamba_options')
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                samtools.view(
                                    inputTumor,
                                    None,
                                    "-b -F 1294"
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    discordants_tumor,
                                    lumpy_directory,
                                    config.param('extract_discordant_reads', 'sambamba_options')
                                )
                            ]
                        )
                    ],
                    name="extract_discordant_reads." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            lumpy_directory,
                            remove=True
                        ),
                        pipe_jobs(
                            [
                                samtools.view(
                                    inputNormal,
                                    None,
                                    "-h"
                                ),
                                Job(
                                    [None],
                                    [None],
                                    [['lumpy_sv', 'module_lumpy']],
                                    command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"
                                ),
                                samtools.view(
                                    "-",
                                    None,
                                    " -Sb "
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    splitters_normal,
                                    lumpy_directory,
                                    config.param('extract_split_reads', 'sambamba_options')
                                ),
                            ]
                        ),
                        pipe_jobs(
                            [
                                samtools.view(
                                    inputTumor,
                                    None,
                                    "-h"
                                ),
                                Job(
                                    [None],
                                    [None],
                                    [['lumpy_sv', 'module_lumpy']],
                                    command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"
                                ),
                                samtools.view(
                                    "-",
                                    None,
                                    " -Sb "
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    splitters_tumor,
                                    lumpy_directory,
                                    config.param('extract_split_reads', 'options')
                                )
                            ]
                        )
                    ],
                    name="extract_split_reads." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            lumpy_directory,
                            remove=True
                        ),
                        lumpy.lumpyexpress_pair(
                            inputNormal,
                            inputTumor,
                            output_vcf,
                            spl_normal=splitters_normal,
                            spl_tumor=splitters_tumor,
                            dis_normal=discordants_normal,
                            dis_tumor=discordants_tumor
                        ),
                        htslib.bgzip(
                            output_vcf,
                            gzip_vcf
                        ),
                    ],
                    name="lumpy_paired_sv_calls." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                Job(
                                    [gzip_vcf],
                                    [None],
                                    command="zcat " + gzip_vcf + " | grep -v \"^##contig\""
                                ),
                                bcftools.annotate(
                                    None,
                                    None,
                                    config.param('lumpy_paired_sv_calls', 'header_options')
                                ),
                                vt.sort(
                                    "-",
                                    os.path.join(pair_directory, tumor_pair.name + ".lumpy.sorted.vcf"),
                                    "-m full"
                                )
                            ]
                        ),
                        svtyper.genotyper(
                            inputTumor,
                            inputNormal,
                            os.path.join(pair_directory, tumor_pair.name + ".lumpy.sorted.vcf"),
                            genotype_vcf
                        ),
                        htslib.bgzip(
                            genotype_vcf,
                            genotype_gzip
                        )
                    ],
                    name="lumpy_paired_sv_calls.genotype." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def lumpy_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name)
            prefix = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            
            genotype_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf")
            somatic_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.somatic.vcf.gz")
            germline_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.germline.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                vawk.paired_somatic(
                                    genotype_vcf,
                                    tumor_pair.normal.name,
                                    tumor_pair.tumor.name,
                                    None
                                ),
                                htslib.bgzip(
                                    None,
                                    somatic_vcf
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vawk.paired_germline(
                                    genotype_vcf,
                                    tumor_pair.normal.name,
                                    tumor_pair.tumor.name,
                                    None
                                ),
                                htslib.bgzip(
                                    None,
                                    germline_vcf
                                )
                            ]
                        )
                    ],
                    name="sv_annotation.lumpy.genotypes." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            somatic_vcf,
                            prefix + ".lumpy.somatic.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            prefix + ".lumpy.somatic.snpeff.vcf",
                            prefix + ".lumpy.somatic.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            prefix + ".lumpy.somatic.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "LUMPY",
                            prefix + ".lumpy.somatic.prioritize.tsv"
                        ),
                    ],
                    name="sv_annotation.lumpy.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            germline_vcf,
                            prefix + ".lumpy.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            prefix + ".lumpy.germline.snpeff.vcf",
                            prefix + ".lumpy.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            prefix + ".lumpy.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "LUMPY",
                            prefix + ".lumpy.germline.prioritize.tsv"
                        ),
                    ],
                    name="sv_annotation.lumpy.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def sym_link_lumpy(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".lumpy.somatic.snpeff.annot.vcf",
                pair_directory + ".lumpy.somatic.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_lumpy.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".lumpy.germline.snpeff.annot.vcf",
                pair_directory + ".lumpy.germline.prioritize.tsv"
            ]
        
            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_lumpy.germline." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs

    def wham_call_sv(self):
        """
        Wham (Whole-genome Alignment Metrics) to provide a single, integrated framework for both structural variant
        calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts
        to employ SVs in association testing.
        Returns:vcf.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name)
            wham_directory = os.path.join(pair_directory, "rawWham")
            inputNormal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            output_vcf = os.path.join(wham_directory, tumor_pair.name + ".wham.vcf")
            merge_vcf = os.path.join(wham_directory, tumor_pair.name + ".wham.merged.vcf")
            genotyped_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.merged.genotyped.vcf.gz")
            somatic_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.somatic.vcf.gz")
            germline_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.germline.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            wham_directory,
                            remove=True
                        ),
                        wham.call_sv(
                            inputTumor,
                            inputNormal,
                            output_vcf
                        ),
                        pipe_jobs(
                            [
                                wham.merge(
                                    output_vcf,
                                    None
                                ),
                                Job(
                                    [None],
                                    [merge_vcf],
                                    command="sed 's/NONE/" + tumor_pair.tumor.name + "/g' | sed -e 's#\"\"#\"#g' > " + merge_vcf
                                ),
                            ]
                        ),
                    ],
                    name="wham_call_sv.call_merge." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            wham_directory,
                            remove=True
                        ),
                        pipe_jobs(
                            [
                                Job(
                                    [merge_vcf],
                                    [None],
                                    command="cat " + merge_vcf + " | grep -v \"^##contig\""
                                ),
                                bcftools.annotate(
                                    None,
                                    None,
                                    config.param('wham_call_sv', 'header_options')
                                ),
                                vt.sort(
                                    "-",
                                    os.path.join(pair_directory, tumor_pair.name + ".wham.sorted.vcf"),
                                    "-m full"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                svtyper.genotyper(
                                    inputTumor,
                                    inputNormal,
                                    os.path.join(pair_directory, tumor_pair.name + ".wham.sorted.vcf"),
                                    None
                                ),
                                Job(
                                    [None],
                                    [None],
                                    command=" sed -e 's#\"\"#\"#g' "
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    genotyped_vcf
                                )
                            ]
                        )
                    ],
                    name="wham_call_sv.genotype." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def wham_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name)
            
            genotyped_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.merged.genotyped.vcf.gz")

            prefix = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            
            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                vawk.paired_somatic(
                                    genotyped_vcf,
                                    tumor_pair.normal.name,
                                    tumor_pair.tumor.name,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    prefix + ".wham.somatic.vcf.gz"
                                ),
                            ]
                        ),
                        snpeff.compute_effects(
                            prefix + ".wham.somatic.vcf.gz",
                            prefix + ".wham.somatic.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            prefix + ".wham.somatic.snpeff.vcf",
                            prefix + ".wham.somatic.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            prefix + ".wham.somatic.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "WHAM",
                            prefix + ".wham.somatic.prioritize.tsv"
                        ),
                    ], name="sv_annotation.wham.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                vawk.paired_germline(
                                    genotyped_vcf,
                                    tumor_pair.normal.name,
                                    tumor_pair.tumor.name,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    prefix + ".wham.germline.vcf.gz"
                                )
                            ]
                        ),
                        snpeff.compute_effects(
                            prefix + ".wham.germline.vcf.gz",
                            prefix + ".wham.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            prefix + ".wham.germline.snpeff.vcf",
                            prefix + ".wham.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            prefix + ".wham.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "WHAM", prefix + ".wham.germline.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.wham.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def sym_link_wham(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".wham.somatic.snpeff.annot.vcf",
                pair_directory + ".wham.somatic.prioritize.tsv"
            ]
            
            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_wham.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".wham.germline.snpeff.annot.vcf",
                pair_directory + ".wham.germline.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_wham.germline." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs

    def cnvkit_batch(self):
        """
        CNVkit is a Python library and command-line software toolkit to infer and visualize copy number from
        high-throughput DNA sequencing data. It is designed for use with hybrid capture, including both whole-exome and
        custom target panels, and short-read sequencing platforms such as Illumina and Ion Torrent.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            cnvkit_dir = os.path.join(pair_directory, "rawCNVkit")
            inputNormal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            tarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.targetcoverage.cnn")
            antitarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.antitargetcoverage.cnn")
            ref_cnn = os.path.join(cnvkit_dir, tumor_pair.name + ".reference.cnn")
            tumor_cns = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".cns")
            vcf_gz = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            metrics = os.path.join(self.output_dirs['sv_variants_directory'], "cnvkit_reference")
            poolRef = os.path.join(self.output_dir, metrics, "pooledReference.cnn")

            if os.path.isfile(poolRef):
                pool_ref_cnn = poolRef
                ref_cnn = None

            else:
                pool_ref_cnn = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            bed = None

            if coverage_bed:
                bed = coverage_bed

            vardict_vcf = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, tumor_pair.name + ".vardict.germline.vt.vcf.gz")

            input_vcf = None
            normal = None
            tumor = None
            if os.path.isfile(vardict_vcf):
                input_vcf = vardict_vcf
                normal = tumor_pair.normal.name
                tumor = tumor_pair.tumor.name

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.batch(
                            inputTumor,
                            inputNormal,
                            cnvkit_dir,
                            tar_dep=tarcov_cnn,
                            antitar_dep=antitarcov_cnn,
                            target_bed=bed,
                            reference=pool_ref_cnn,
                            output_cnn=ref_cnn
                        )
                    ],
                    name="cnvkit_batch." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.fix(
                            tarcov_cnn,
                            antitarcov_cnn,
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            reference=pool_ref_cnn,
                            ref_cnn=ref_cnn
                        ),
                        cnvkit.segment(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            tumor_cns
                        )
                    ],
                    name="cnvkit_batch.correction." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.call(
                            tumor_cns,
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns")
                        ),
                        pipe_jobs(
                            [
                                cnvkit.export(
                                    os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                                    None,
                                    sample_id=tumor_pair.tumor.name
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    vcf_gz
                                )
                            ]
                        )
                    ],
                    name="cnvkit_batch.call." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.metrics(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                            os.path.join(metrics, tumor_pair.name + ".metrics.tsv")
                        ),
                        cnvkit.scatter(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".scatter.pdf"),
                            input_vcf,
                            normal,
                            tumor
                        ),
                        cnvkit.diagram(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".diagram.pdf")
                        )
                    ],
                    name="cnvkit_batch.metrics." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def cnvkit_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name, tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            pair_directory + ".cnvkit.vcf.gz",
                            pair_directory + ".cnvkit.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            pair_directory + ".cnvkit.snpeff.vcf",
                            pair_directory + ".cnvkit.snpeff.annot.vcf"
                        ),
                    ],
                    name="sv_annotation.cnvkit." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def sym_link_cnvkit(self):
        jobs = []
    
        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dir,"SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [pair_directory + ".cnvkit.snpeff.annot.vcf"]
        
            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_cnvkit.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs
     
    def ensemble_metasv_somatic(self):
        """
        MetaSV: An accurate and integrative structural-variant caller for next generation sequencing
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            ensemble_directory = os.path.join(self.output_dirs['sv_variants_directory'], "ensemble", tumor_pair.name)

            inputTumor = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            isize_file = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "picard_metrics", "picard_metrics.all.metrics.insert_size_metrics")
            gatk_vcf = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vcf.gz")
            lumpy_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.somatic.vcf.gz")
            manta_vcf = os.path.join(pair_directory, tumor_pair.name + ".manta.somatic.vcf.gz")
            wham_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.somatic.vcf.gz")
            delly_vcf= os.path.join(pair_directory, tumor_pair.name + ".delly.somatic.vcf.gz")
            cnvkit_vcf = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            if os.path.isfile(isize_file):
                isize_mean, isize_sd = metric_tools.extract_isize(
                    isize_file
                )

            else:
                isize_mean = 325
                isize_sd = 75
                
            gatk_pass = None
            if os.path.isfile(gatk_vcf):
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                ensemble_directory,
                                remove=True
                            ),
                            bcftools.view(
                                gatk_vcf,
                                gatk_pass,
                                config.param('metasv_ensemble', 'filter_somatic_options')
                            ),
                        ],
                        name="metasv_ensemble.ensemble_pass." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )
            
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        metasv.ensemble(
                            lumpy_vcf,
                            manta_vcf,
                            cnvkit_vcf,
                            wham_vcf,
                            delly_vcf,
                            gatk_pass,
                            inputTumor,
                            tumor_pair.tumor.name,
                            os.path.join(ensemble_directory, "rawMetaSV_somatic"),
                            ensemble_directory,
                            isize_mean=str(isize_mean),
                            isize_sd=str(isize_sd),
                            output_vcf=os.path.join(ensemble_directory, "variants.vcf.gz")
                        ),
                    ],
                    name="metasv_ensemble." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def ensemble_metasv_germline(self):
        """
        MetaSV: An accurate and integrative structural-variant caller for next generation sequencing
        """
        jobs = []
    
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            ensemble_directory = os.path.join(self.output_dirs['sv_variants_directory'], "ensemble", tumor_pair.name)
        
            inputTumor = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            isize_file = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "picard_metrics", "picard_metrics.all.metrics.insert_size_metrics")
            gatk_vcf = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, tumor_pair.name + ".ensemble.germline.vcf.gz")
            gatk_pass = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, tumor_pair.name + ".ensemble.germline.flt.pass.vcf.gz")
            lumpy_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.germline.vcf.gz")
            manta_vcf = os.path.join(pair_directory, tumor_pair.name + ".manta.germline.vcf.gz")
            wham_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.germline.vcf.gz")
            delly_vcf = os.path.join(pair_directory, tumor_pair.name + ".delly.germline.vcf.gz")
            cnvkit_vcf = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            if os.path.isfile(isize_file):
                isize_mean, isize_sd = metric_tools.extract_isize(
                    isize_file
                )
        
            else:
                isize_mean = 325
                isize_sd = 75
        
            gatk_pass = None
            if os.path.isfile(gatk_vcf):
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                ensemble_directory,
                                remove=True
                            ),
                            bcftools.view(
                                gatk_vcf,
                                gatk_pass,
                                config.param('metasv_ensemble', 'filter_germline_options')
                            ),
                        ],
                        name="metasv_ensemble.ensemble_pass." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor]
                    )
                )
        
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        metasv.ensemble(
                            lumpy_vcf,
                            manta_vcf,
                            cnvkit_vcf,
                            wham_vcf,
                            delly_vcf,
                            gatk_pass,
                            inputTumor,
                            tumor_pair.tumor.name,
                            os.path.join(ensemble_directory, "rawMetaSV_germline"),
                            ensemble_directory,
                            isize_mean=str(isize_mean),
                            isize_sd=str(isize_sd),
                            output_vcf=os.path.join(ensemble_directory, "variants.vcf.gz")
                        ),
                    ],
                    name="metasv_ensemble." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )
    
        return jobs

    def metasv_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            ensemble_directory = os.path.join(self.output_dirs['sv_variants_directory'], "ensemble", tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            os.path.join(ensemble_directory, "variants.vcf.gz"),
                            os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.vcf")
                        ),
                        annotations.structural_variants(
                            os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.vcf"),
                            os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.annot.vcf")
                        ),
                    ],
                    name="sv_annotation.metasv_ensemble." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def sym_link_metasv(self):
        jobs = []
    
        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], "ensemble", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [pair_directory + ".metasv.snpeff.annot.vcf"]
        
            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_metasv." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )
                    
        return jobs

    def scones(self):
        """
        This step aims to estimate somatic Copy Number Variation using BVAtools and SCoNEs. BVAtools generate the bined Depth ratio values from the
        tumor and normal BAM files. SCoNEs is tool to deconvolution the logR signal of the tumor-normal coverage into a mixture of baysian sub-signal
        for each copy number state. The result is a set of several deconvolution using  0-7 sub-signal. As each tumor sample is unique the choice of
        the best final model (number of sub-signal) needs to be manually evaluated using the log ratio graphical representation.

        """
        window_size = config.param('scones', 'window', required=True)
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            sv_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            scones_directory = os.path.join(sv_directory, "SCoNEs")
            inputNormal = self.select_input_files(
                [[os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]])[0]
            inputTumor = self.select_input_files(
                [[os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]])[0]

            bined_count_file = os.path.join(scones_directory, tumor_pair.normal.name + ".bin" + window_size + ".tsv")
            bined_count_fix_file = os.path.join(scones_directory, tumor_pair.normal.name + ".bin" + window_size + ".fix.tsv")

            output_scones_basename = os.path.join(scones_directory,
                                                  tumor_pair.normal.name + ".bin" + window_size + "_SCoNEs")
            scones_best_model_basename = output_scones_basename + "_Model_" + config.param('scones', 'best_model',
                                                                                           required=True)
            scones_calls_file = scones_best_model_basename + "_CNVcalls.txt"
            scones_filtered_file = scones_best_model_basename + "_CNVcalls.filtered.tsv"
            scones_annotate_basename = scones_best_model_basename + "_CNVcalls.filtered.anotated"
            scones_annotate_tmp_basename = scones_best_model_basename + "_CNVcalls.filtered.tmp"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            scones_directory,
                            remove=True
                        ),
                        bvatools.bincounter(
                            bam=inputTumor,
                            refbam=inputNormal,
                            out=bined_count_fix_file,
                            window=window_size
                        ),
                        Job(
                            [bined_count_fix_file],
                            [bined_count_file],
                                command="cat <(head -1 " + bined_count_fix_file + ") <(grep -v \"_\" " + bined_count_fix_file
                                    + " | grep -v \"EBV\" ) > " + bined_count_file
                        ),
                    ],
                    name="bvatools_bincounter." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            scones_directory,
                            remove=True
                        ),
                        scones.scones_pair(
                            bined_file=bined_count_file,
                            output_basename=output_scones_basename,
                            window=window_size
                        )
                    ],
                    name="scones_pair." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            scones_directory,
                            remove=True
                        ),
                        scones.scones_filter(
                            scones_calls=scones_calls_file,
                            pair_name=tumor_pair.name,
                            output=scones_filtered_file
                        )
                    ],
                    name="scones_filter." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            scones_directory,
                            remove=True
                        ),
                        scones.scones_annotate(
                            scones_calls_filtered=scones_filtered_file,
                            output_basename=scones_annotate_basename,
                            tmp_basename=scones_annotate_tmp_basename
                        )
                    ],
                    name="scones_annotate." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def svaba_assemble(self):
        """
        SvABA - Structural variation and indel analysis by assembly
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
    
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            svaba_directory = os.path.join(pair_directory, "rawSvaba")

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            somatic_input = tumor_pair.name + ".svaba.somatic.sv.vcf"
            somatic_output = os.path.join(pair_directory, tumor_pair.name + ".svaba.somatic.vcf")

            germline_input = tumor_pair.name + ".svaba.germline.sv.vcf"
            germline_output = os.path.join(pair_directory, tumor_pair.name + ".svaba.germline.vcf")

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.tumor.readsets[0]
            )
            bed = None

            if coverage_bed:
                bed = coverage_bed

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            svaba_directory,
                            remove=True
                        ),
                        Job(
                            command="cd " + svaba_directory
                        ),
                        svaba.run(
                            input_tumor,
                            tumor_pair.name,
                            input_normal,
                            bed
                        ),
                        Job(
                            [somatic_input],
                            [somatic_output],
                            command="sed -e 's#" + input_normal + "#" + tumor_pair.normal.name + "#g' " + somatic_input + " | " + "sed -e 's#" + input_tumor + "#" + tumor_pair.tumor.name + "#g' > " + somatic_output),
                        Job(
                            [germline_input],
                            [germline_output],
                            command="sed -e 's#" + input_normal + "#" + tumor_pair.normal.name + "#g' " + germline_input + " | " + "sed -e 's#" + input_tumor + "#" + tumor_pair.tumor.name + "#g' > " + germline_output)
                    ],
                    name="svaba_run." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def svaba_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        Job(
                            [pair_directory + ".svaba.somatic.vcf"],
                            [pair_directory + ".svaba.somatic.flt.vcf"],
                            command="cat <(grep \"^#\" " + pair_directory
                                    + ".svaba.somatic.vcf) <(grep -v \"^#\" " + pair_directory
                                    + ".svaba.somatic.vcf | cut -f1-9,13-14) > " + pair_directory + ".svaba.somatic.flt.vcf"
                        ),
                        snpeff.compute_effects(
                            pair_directory + ".svaba.somatic.flt.vcf",
                            pair_directory + ".svaba.somatic.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            pair_directory + ".svaba.somatic.snpeff.vcf",
                            pair_directory + ".svaba.somatic.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            pair_directory + ".svaba.somatic.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "SVABA",
                            pair_directory + ".svaba.somatic.prioritize.tsv"
                        ),
                    ],
                    name="sv_annotation.svaba_somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )       

            jobs.append(
                concat_jobs(
                    [
                        Job(
                            [pair_directory + ".svaba.germline.vcf"],
                            [pair_directory + ".svaba.germline.flt.vcf"],
                            command="cat <(grep \"^#\" " + pair_directory + ".svaba.germline.vcf) <(grep -v \"^#\" "
                                    + pair_directory + ".svaba.germline.vcf | cut -f1-9,13-14) > "
                                    + pair_directory + ".svaba.germline.flt.vcf"
                        ),
                        snpeff.compute_effects(
                            pair_directory + ".svaba.germline.flt.vcf",
                            pair_directory + ".svaba.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            pair_directory + ".svaba.germline.snpeff.vcf",
                            pair_directory + ".svaba.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            pair_directory + ".svaba.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "SVABA",
                            pair_directory + ".svaba.germline.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.svaba_germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor]
                )
            )

        return jobs

    def sym_link_svaba(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".svaba.somatic.snpeff.annot.vcf",
                pair_directory + ".svaba.somatic.prioritize.tsv"
            ]
                               
            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_file):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_svaba.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        inputs = dict()
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".svaba.germline.sv.snpeff.annot.vcf",
                pair_directory + ".svaba.germline.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                            ],
                            name="sym_link_svaba.germline." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor]
                        )
                    )

        return jobs

    @property
    def steps(self):
        return [
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.manta_sv_calls,
                self.rawmpileup_panel,
                self.paired_varscan2_panel,
                self.merge_varscan2_panel,
                self.preprocess_vcf_panel,
                self.snp_effect_panel,
                self.gemini_annotations_panel,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_fastqc,
                self.sequenza,
                self.run_pair_multiqc,
                self.sym_link_report,
                self.sym_link_fastq_pair,
                self.sym_link_panel,
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_fastqc,
                self.sequenza,
                self.manta_sv_calls,
                self.strelka2_paired_somatic,
                self.strelka2_paired_germline,
                self.purple,
                self.rawmpileup,
                self.paired_varscan2,
                self.merge_varscan2,
                self.paired_mutect2,
                self.merge_mutect2,
                self.vardict_paired,
                self.merge_filter_paired_vardict,
                self.ensemble_somatic,
                self.gatk_variant_annotator_somatic,
                self.merge_gatk_variant_annotator_somatic,
                self.ensemble_germline_loh,
                self.gatk_variant_annotator_germline,
                self.merge_gatk_variant_annotator_germline,
                self.cnvkit_batch,
                self.filter_ensemble_germline,
                self.filter_ensemble_somatic,
                self.report_cpsr,
                self.report_pcgr,
                self.run_pair_multiqc,
                self.sym_link_fastq_pair,
                self.sym_link_final_bam,
                self.sym_link_report,
                self.sym_link_ensemble,
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.strelka2_paired_somatic,
                self.strelka2_paired_germline,
                self.metrics_dna_picard_metrics,
                self.sequenza,
                self.delly_call_filter,
                self.delly_sv_annotation,
                self.manta_sv_calls,
                self.manta_sv_annotation,
                self.lumpy_paired_sv,
                self.lumpy_sv_annotation,
                self.wham_call_sv,
                self.wham_sv_annotation,
                self.cnvkit_batch,
                self.cnvkit_sv_annotation,
                self.scones,
                self.svaba_assemble,
                self.svaba_sv_annotation,
                self.ensemble_metasv_somatic,
                self.ensemble_metasv_germline,
                self.metasv_sv_annotation,
                self.sym_link_sequenza,
                self.sym_link_metasv,
                self.sym_link_delly,
                self.sym_link_manta,
                self.sym_link_lumpy,
                self.sym_link_wham,
                self.sym_link_cnvkit,
                #self.sym_link_svaba
            ]
        ]


if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        TumorPair(protocol=['fastpass', 'ensemble', 'sv'])
