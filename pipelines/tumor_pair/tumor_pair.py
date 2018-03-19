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
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
import gzip
from sys import stderr
from core.pipeline import *
from bfx.sample_tumor_pairs import *
from bfx.sequence_dictionary import *
from pipelines.dnaseq import dnaseq

#utilizes
from bfx import sambamba
from bfx import vcflib
from bfx import bcftools
from bfx import tools
from bfx import bed_file
from bfx import metric_tools
from bfx import bvatools
from bfx import vt
from bfx import snpeff
from bfx import vawk
from bfx import deliverables

#metrics
from bfx import conpair
from bfx import multiqc

#variants
from bfx import htslib
from bfx import samtools
from bfx import varscan
from bfx import gatk
from bfx import vardict
from bfx import bcbio_variation_recall
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
from bfx import shapeit
from bfx import scnaphase
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
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=file)
        self.argparser.add_argument("--profyle", help="adjust deliverables to PROFYLE folder conventions (Default: False)", action="store_true")
        self.argparser.add_argument("-t", "--type", help="Tumor pair analysis type",choices = ["fastpass", "ensemble", "sv"], default="ensemble")
        super(TumorPair, self).__init__(protocol)

        


    #def __args__(self):
        #self.args.type.help = "Tumor pair analysis type"
        #self.args.type.choices = ["fastpass", "ensemble", "sv"]
        #self.args.type.default="ensemble"
        #return self.args


    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(self.args.pairs.name, self.samples, self.args.profyle)
        return self._tumor_pairs

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(
                config.param('DEFAULT', 'genome_dictionary', type='filepath'), variant=True)
        return self._sequence_dictionary_variant

    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary_variant()):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in
                    self.sequence_dictionary_variant()]
        else:
            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary_variant()])
            approximate_window_size = int(
                math.floor(total_length / (nb_jobs - len(self.sequence_dictionary_variant()))))
            windows = []

            for sequence in self.sequence_dictionary_variant():
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in
                                   range(1, sequence['length'] + 1, approximate_window_size)]:
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

    def build_ped_file(self, directory, tumor_pair):
        ped_file = os.path.join(directory, tumor_pair.name + ".ped")
        ped_job = Job(
            command="""\
`cat > {ped_file} << END
#Family_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\tPhenotype\tEthnicity
1\t{normal}\t-9\t-9\t0\t1\t-9
1\t{tumor}\t-9\t-9\t0\t2\t-9
END`""".format(
            ped_file=ped_file,
            normal=tumor_pair.normal.name,
            tumor=tumor_pair.tumor.name,
            )
        )

        return ped_job

    def sym_link_fastq_pair(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Normal"] = [os.path.join("raw_reads", readset.sample.name, readset.name) for readset in tumor_pair.readsets[tumor_pair.normal.name]]
            inputs["Tumor"] = [os.path.join("raw_reads", readset.sample.name, readset.name) for readset in tumor_pair.readsets[tumor_pair.tumor.name]]

            for key,input in inputs.iteritems():
                for readset in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(readset + ".pair1.fastq.gz", tumor_pair, type="raw_reads", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(readset + ".pair2.fastq.gz", tumor_pair, type="raw_reads", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_fastq.pairs." + readset))

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

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            normal_alignment_directory = os.path.join("alignment", tumor_pair.normal.name)
            normal_realign_directory = os.path.join(normal_alignment_directory, "realign")
            tumor_alignment_directory = os.path.join("alignment", tumor_pair.tumor.name)
            tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_intervals = os.path.join(tumor_realign_directory, "all.intervals")
                bam_postfix = ".all.realigned.bam"
                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.all.realigned.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                normal_output_bam = os.path.join(normal_alignment_directory,
                                                 tumor_pair.normal.name + ".realigned.qsorted.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.all.realigned.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                tumor_output_bam = os.path.join(tumor_alignment_directory,
                                                tumor_pair.tumor.name + ".realigned.qsorted.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory], samples=[tumor_pair.normal]),
                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory], samples=[tumor_pair.tumor]),
                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor),
                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam,
                                         output_tum_dep=tumor_output_bam, target_intervals=realign_intervals,
                                         optional=bam_postfix),
                    # Move sample realign 
                    Job([input_normal], [normal_output_bam],
                        command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                    Job([input_tumor], [tumor_output_bam],
                        command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                ], name="gatk_indel_realigner." + tumor_pair.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary,
                                                                                          nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(tumor_realign_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequences
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    bam_postfix = ".realigned." + str(idx) + ".bam"
                    normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_index = re.sub("\.bam$", ".bai", normal_bam)
                    tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                    normal_output_bam = os.path.join(normal_realign_directory,
                                                     tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                    tumor_output_bam = os.path.join(tumor_realign_directory,
                                                    tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory], samples=[tumor_pair.normal]),
                        Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory], samples=[tumor_pair.tumor]),
                        gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor,
                                                      intervals=intervals),
                        gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam,
                                             output_tum_dep=tumor_output_bam, target_intervals=realign_intervals,
                                             intervals=intervals, optional=bam_postfix),
                        Job([input_normal], [normal_output_bam],
                            command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                        Job([input_tumor], [tumor_output_bam],
                            command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                    ], name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(tumor_realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                bam_postfix = ".realigned.others.bam"
                normal_bam = os.path.join(tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                tumor_bam = os.path.join(tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                normal_output_bam = os.path.join(normal_realign_directory,
                                                 tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_output_bam = os.path.join(tumor_realign_directory,
                                                tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + normal_realign_directory, removable_files=[normal_realign_directory], samples=[tumor_pair.normal]),
                    Job(command="mkdir -p " + tumor_realign_directory, removable_files=[tumor_realign_directory], samples=[tumor_pair.tumor]),
                    gatk.realigner_target_creator(input_normal, realign_intervals, input2=input_tumor,
                                                  exclude_intervals=unique_sequences_per_job_others),
                    gatk.indel_realigner(input_normal, input2=input_tumor, output_norm_dep=normal_output_bam,
                                         output_tum_dep=tumor_output_bam, target_intervals=realign_intervals,
                                         exclude_intervals=unique_sequences_per_job_others, optional=bam_postfix),
                    Job([input_normal], [normal_output_bam],
                        command="mv " + normal_bam + " " + normal_output_bam + " && mv " + normal_index + " " + normal_output_index),
                    Job([input_tumor], [tumor_output_bam],
                        command="mv " + tumor_bam + " " + tumor_output_bam + " && mv " + tumor_index + " " + tumor_output_index)
                ], name="gatk_indel_realigner." + tumor_pair.name + ".others"))

        return jobs

    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "realigned.sorted.bam"
            output = alignment_file_prefix + "sorted.dup.bam"

            job = sambamba.markdup(input, output, os.path.join("alignment", sample.name, sample.name))
            job.name = "sambamba_mark_duplicates." + sample.name
            job.samples=[sample]
            jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.picard_mark_duplicates.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="picard_mark_duplicates_report")
        )

        return jobs

    def sym_link_final_bam(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Normal"] = [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal")]
            inputs["Tumor"] =  [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal")]

            for key, input in inputs.iteritems():
                for sample_bam in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample_bam + ".bam", tumor_pair, type="final_data", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample_bam + ".bai", tumor_pair, type="final_data", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample_bam + ".bam.md5", tumor_pair, type="final_data", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_final_bam.pairs." + tumor_pair.name + "." + key))

        return jobs

    def conpair_concordance_contamination(self):
        """
        
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            metrics_directory = os.path.join("metrics")
            input_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")
            pileup_normal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".gatkPileup")
            pileup_tumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".gatkPileup")

            concordance_out = os.path.join(metrics_directory, tumor_pair.name + ".concordance.tsv")
            contamination_out = os.path.join(metrics_directory, tumor_pair.name + ".contamination.tsv")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal]),
                conpair.pileup(input_normal, pileup_normal),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.normal.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.tumor]),
                conpair.pileup(input_tumor, pileup_tumor),
            ], name="conpair_concordance_contamination.pileup." + tumor_pair.tumor.name))

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + metrics_directory, samples=[tumor_pair.normal, tumor_pair.tumor]),
                conpair.concordance(pileup_normal, pileup_tumor, concordance_out),
                conpair.contamination(pileup_normal, pileup_tumor, contamination_out)
            ], name="conpair_concordance_contamination." + tumor_pair.name))

        return jobs

    def rawmpileup_panel(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            bedfile = config.param('rawmpileup_panel', 'panel')

            for sequence in self.sequence_dictionary_variant():
                pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    samtools.mpileup(
                        [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam"),
                         os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                        pair_output, config.param('rawmpileup_panel', 'mpileup_other_options'), region=sequence['name'],
                        regionFile=bedfile),
                ], name="rawmpileup_panel." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def paired_varscan2_panel(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                output_vcf_gz = os.path.join(varscan_directory,
                                             tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    varscan.somatic(input_pair, output, config.param('varscan2_somatic_panel', 'other_options'),
                                    output_vcf_dep=output_vcf_gz, output_snp_dep=output_snp,
                                    output_indel_dep=output_indel),
                    htslib.bgzip_tabix(output_snp, os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz")),
                    htslib.bgzip_tabix(output_indel, os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")),
                    pipe_jobs([
                        bcftools.concat(
                            [os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz"),
                             os.path.join(varscan_directory,
                                          tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")], None),
                        Job([None], [None],
                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' "),
                        htslib.bgzip_tabix(None, output_vcf_gz),
                    ]),
                ], name="varscan2_somatic_panel." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_varscan2_panel(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")
                          for sequence in self.sequence_dictionary_variant()]

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            somatic_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz")
            germline_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vcf.gz")

            for input_vcf in all_inputs:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete panel varscan2 vcf: %s\n" % input_vcf)

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(all_inputs, None),
                    tools.fix_varscan_output(None, None),
                    Job([None], [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                    Job([None], [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                    htslib.bgzip_tabix(None, all_output),
                ]),
                pipe_jobs([
                    bcftools.view(all_output, None, config.param('merge_varscan2', 'somatic_filter_options')),
                    htslib.bgzip_tabix(None, somatic_output),
                ]),
                pipe_jobs([
                    bcftools.view(all_output, None, config.param('merge_varscan2', 'germline_loh_filter_options')),
                    htslib.bgzip_tabix(None, germline_output),
                ]),
            ], name="merge_varscan2." + tumor_pair.name))

        return jobs

    def preprocess_vcf_panel(self):
        """
        Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
        vcf FORMAT modification for correct loading into gemini
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")

            prefix = os.path.join(pair_directory, tumor_pair.name)
            output_somatic = prefix + ".varscan2.somatic.vt.vcf.gz"

            output_germline = prefix + ".varscan2.germline_loh.vt.vcf.gz"

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(prefix + ".varscan2.somatic.vcf.gz", None),
                    htslib.bgzip_tabix(None, prefix + ".prep.vt.vcf.gz"),
                ]),
                tools.preprocess_varscan(prefix + ".prep.vt.vcf.gz", output_somatic),
            ], name="preprocess_vcf_panel.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(prefix + ".varscan2.germline_loh.vcf.gz", None),
                    htslib.bgzip_tabix(None, prefix + ".germline_loh.prep.vt.vcf.gz"),
                ]),
                tools.preprocess_varscan(prefix + ".germline_loh.prep.vt.vcf.gz", output_germline),
            ], name="preprocess_vcf_panel.germline." + tumor_pair.name))

        return jobs

    def snp_effect_panel(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            input_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf")
            output_somatic_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf.gz")

            input_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")
            output_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.snpeff.vcf")
            output_germline_gz = os.path.join(pair_directory,
                                              tumor_pair.name + ".varscan2.germline_loh.vt.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(varscan_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename,
                                       options=config.param('compute_cancer_effects_somatic', 'options')),
                htslib.bgzip_tabix(output_somatic, output_somatic_gz),
            ], name="compute_cancer_effects_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                snpeff.compute_effects(input_germline, output_germline, cancer_sample_file=cancer_pair_filename,
                                       options=config.param('compute_cancer_effects_germline', 'options')),
                htslib.bgzip_tabix(output_germline, output_germline_gz),
            ], name="compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def gemini_annotations_panel(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if not os.path.exists(varscan_directory):
                os.makedirs(varscan_directory)

            temp_dir = os.path.join(os.getcwd(), pair_directory)
            gemini_prefix = os.path.join(pair_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(gemini_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",
                                          gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)
            ], name="gemini_annotations.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(gemini_prefix + ".varscan2.germline_loh.vt.snpeff.vcf.gz",
                                          gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)
            ], name="gemini_annotations.germline." + tumor_pair.name))

        return jobs

    def set_somatic_and_actionable_mutations_panel(self):
        """

        """

        jobs = []

        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        ped_file = config.param('set_somatic_and_actionable_mutations', 'ped_file', required=False, type='filepath')
        ped_job = None

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join("pairedVariants", tumor_pair.name, "panel")
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)

            if not ped_file:
                ped_job = self.build_ped_file(paired_directory, tumor_pair)
                ped_file = os.path.join(gemini_prefix + ".ped")

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + paired_directory),
                ped_job,
                gemini.set_somatic(ped_file, gemini_prefix + ".somatic.gemini." + gemini_version + ".db",
                                   gemini_prefix + ".varscan2.somatic.gemini.set_somatic.tsv"),
                gemini.actionable_mutations(gemini_prefix + ".somatic.gemini." + gemini_version + ".db",
                                            gemini_prefix + ".varscan2.somatic.gemini.actionable.tsv")
            ], name="set_somatic_and_actionable_mutations." + tumor_pair.name))

        return jobs

    def sym_link_panel(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] =  [os.path.join("pairedVariants", tumor_pair.name, "panel", tumor_pair.name)]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample + ".varscan2.vcf.gz", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.vcf.gz.tbi", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.somatic.vt.snpeff.vcf.gz", tumor_pair, type="panel",sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.somatic.vt.snpeff.vcf.gz.tbi", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.germline_loh.vt.snpeff.vcf.gz", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.germline_loh.vt.snpeff.vcf.gz.tbi", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.somatic.gemini.set_somatic.tsv", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".varscan2.somatic.gemini.actionable.tsv", tumor_pair, type="panel", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_panel." + tumor_pair.name + "." + key))

        return jobs

    def run_pair_multiqc(self):

        jobs = []

        metrics_directory = os.path.join("metrics", "dna")
        inputs = []
        outputs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            input_normal_oxog = os.path.join(metrics_directory, tumor_pair.normal.name, "picard_metrics.oxog_metrics.txt")
            input_normal_qcbias = os.path.join(metrics_directory, tumor_pair.normal.name, "picard_metrics.qcbias_metrics.txt")
            input_normal_all_picard = os.path.join(metrics_directory, tumor_pair.normal.name, "picard_metrics.all.metrics.quality_distribution.pdf")
            input_normal_qualimap = os.path.join(metrics_directory, tumor_pair.normal.name, "qualimap", tumor_pair.normal.name, "genome_results.txt")
            input_normal_fastqc = os.path.join(metrics_directory, tumor_pair.normal.name, "fastqc", tumor_pair.normal.name + ".sorted.dup_fastqc.zip")
            input_normal_flagstat = os.path.join(metrics_directory, tumor_pair.normal.name, "flagstat", tumor_pair.normal.name + ".flagstat")

            input_tumor_oxog = os.path.join(metrics_directory, tumor_pair.tumor.name, "picard_metrics.oxog_metrics.txt")
            input_tumor_qcbias = os.path.join(metrics_directory, tumor_pair.tumor.name, "picard_metrics.qcbias_metrics.txt")
            input_tumor_all_picard = os.path.join(metrics_directory, tumor_pair.tumor.name, "picard_metrics.all.metrics.quality_distribution.pdf")
            input_tumor_qualimap = os.path.join(metrics_directory, tumor_pair.tumor.name, "qualimap", tumor_pair.tumor.name, "genome_results.txt")
            input_tumor_fastqc = os.path.join(metrics_directory, tumor_pair.tumor.name, "fastqc", tumor_pair.tumor.name + ".sorted.dup_fastqc.zip")
            input_tumor_flagstat = os.path.join(metrics_directory, tumor_pair.tumor.name, "flagstat", tumor_pair.tumor.name + ".flagstat")

            input_dep = [input_normal_oxog, input_normal_qcbias, input_normal_all_picard, input_normal_qualimap, input_normal_fastqc, 
                        input_normal_flagstat, input_tumor_oxog, input_tumor_qcbias, input_tumor_all_picard, input_tumor_qualimap, input_tumor_fastqc,
                        input_tumor_flagstat]

            input = [os.path.join(metrics_directory, tumor_pair.normal.name), os.path.join(metrics_directory,tumor_pair.tumor.name)]
            output = os.path.join(metrics_directory, tumor_pair.name + ".multiqc.html")
            inputs.append(input)
            outputs.append(output)

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + os.path.join(metrics_directory, tumor_pair.normal.name) + " " +
                            os.path.join(metrics_directory, tumor_pair.normal.name + "_json"), samples=[tumor_pair.normal, tumor_pair.tumor]),
                Job(command="mkdir -p " + os.path.join(metrics_directory, tumor_pair.tumor.name) + " " +
                            os.path.join(metrics_directory, tumor_pair.tumor.name + "_json"), samples=[tumor_pair.normal, tumor_pair.tumor]),
                multiqc.run(input, output, input_dep=input_dep, sample=metrics_directory),
            ], name="multiqc." + tumor_pair.name))

        #jobs.append(concat_jobs([
        #    multiqc.run(inputs, os.path.join(metrics_directory, "allPairs.multiqc.html") , input_dep=outputs, sample="allPairs"),
        #], name="multiqc.allPairs"))

        return jobs

    def sym_link_report(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] = [os.path.join("metrics", "dna", tumor_pair.name + ".multiqc.html")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="metrics", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_fastq.report." + tumor_pair.name + "." + key))

        return jobs

    def rawmpileup(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    samtools.mpileup([os.path.join("alignment", tumor_pair.normal.name,
                                                   tumor_pair.normal.name + ".sorted.dup.recal.bam"),
                                      os.path.join("alignment", tumor_pair.tumor.name,
                                                   tumor_pair.tumor.name + ".sorted.dup.recal.bam")], pair_output,
                                     config.param('rawmpileup', 'mpileup_other_options'), sequence['name']),
                ], name="rawmpileup." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def rawmpileup_cat(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            mpileup_normal_file_prefix = os.path.join(varscan_directory, tumor_pair.normal.name + ".")
            mpileup_normal_inputs = [mpileup_normal_file_prefix + sequence['name'] + ".mpileup" for sequence in
                                     self.sequence_dictionary_variant()]

            mpileup_tumor_file_prefix = os.path.join(varscan_directory, tumor_pair.tumor.name + ".")
            mpileup_tumor_inputs = [mpileup_tumor_file_prefix + sequence['name'] + ".mpileup" for sequence in
                                    self.sequence_dictionary_variant()]

            normal_output = mpileup_normal_file_prefix + "mpileup"
            tumor_output = mpileup_tumor_file_prefix + "mpileup"

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                Job(mpileup_normal_inputs, [normal_output],
                    command="cat \\\n  " + " \\\n  ".join(mpileup_normal_inputs) + " \\\n  > " + normal_output),
                Job(mpileup_tumor_inputs, [tumor_output],
                    command="cat \\\n  " + " \\\n  ".join(mpileup_tumor_inputs) + " \\\n  > " + tumor_output)
            ], name="rawmpileup_cat." + tumor_pair.name))

        return jobs

    def paired_varscan2(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data. 
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
        SSC INFO field remove to prevent collison with Samtools output during ensemble                     
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            for sequence in self.sequence_dictionary_variant():
                input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory,
                                             tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + varscan_directory, removable_files=[varscan_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),
                    varscan.somatic(input_pair, output, config.param('varscan2_somatic', 'other_options'),
                                    output_vcf_dep=output_vcf, output_snp_dep=output_snp,
                                    output_indel_dep=output_indel),
                    htslib.bgzip_tabix(output_snp, os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz")),
                    htslib.bgzip_tabix(output_indel, os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")),
                    pipe_jobs([
                        bcftools.concat(
                            [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz"),
                             os.path.join(varscan_directory,
                                          tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")], None),
                        Job([None], [output_vcf],
                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' | grep -v \"INFO=<ID=SSC\" | sed -E \"s/SSC=(.*);//g\" > " + output_vcf),
                    ]),
                    htslib.bgzip_tabix(output_vcf, output_vcf_gz),
                ], name="varscan2_somatic." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_varscan2(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            all_inputs = [os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")
                          for sequence in self.sequence_dictionary_variant()]

            for input_vcf in all_inputs:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete varscan2 vcf: %s\n" % input_vcf)

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            all_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vt.vcf.gz")

            somtic_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            germline_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(all_inputs, None),
                    tools.fix_varscan_output(None, None),
                    Job([None], [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                    Job([None], [None],
                        command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                    Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                    htslib.bgzip_tabix(None, all_output),
                ]),
                pipe_jobs([
                    vt.decompose_and_normalize_mnps(all_output, None),
                    htslib.bgzip_tabix(None, all_output_vt),
                ]),
                pipe_jobs([
                    bcftools.view(all_output_vt, None,
                                  config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')),
                    htslib.bgzip_tabix(None, somtic_output_vt),
                ]),
                pipe_jobs([
                    bcftools.view(all_output_vt, None,
                                  config.param('varscan2_readcount_fpfilter', 'germline_loh_filter_options')),
                    htslib.bgzip_tabix(None, germline_output_vt),
                ]),
            ], name="merge_varscan2." + tumor_pair.name))

        return jobs

    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            input_normal = os.path.join("alignment", tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,
                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + mutect_directory, removable_files=[mutect_directory], samples=[tumor_pair.normal, tumor_pair.tumor])
            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(input_normal, input_tumor,
                                 os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"))
                ], name="gatk_mutect2." + tumor_pair.name))

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(
                    self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    outprefix = tumor_pair.name + "." + str(idx) + ".mutect2"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        mkdir_job,
                        gatk.mutect2(input_normal, input_tumor, os.path.join(mutect_directory, outprefix + ".vcf.gz"),
                                     intervals=sequences)
                    ], name="gatk_mutect2." + tumor_pair.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk.mutect2(input_normal, input_tumor,
                                 os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"),
                                 exclude_intervals=unique_sequences_per_job_others)
                ], name="gatk_mutect2." + tumor_pair.name + ".others"))

        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2
        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
        Generate a somatic vcf containing only PASS variants        
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_vt_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            output_somatic_vt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                input_vcf = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                jobs.append(concat_jobs([
                    Job([input_vcf], [output_gz], command="ln -s -f " + input_vcf + " " + output_gz, samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_gz, None),
                        Job([output_vt_gz], [None],
                            command="zcat " + output_vt_gz + " | sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),
                        htslib.bgzip_tabix(None, output_somatic_vt),
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

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        bcftools.concat(inputs, None, config.param('merge_filter_mutect2', 'bcftools_options')),
                        Job([None], [None],
                            command="sed 's/TUMOR/" + tumor_pair.tumor.name + "/g' | sed 's/NORMAL/" + tumor_pair.normal.name + "/g' | sed 's/Number=R/Number=./g' "),
                        htslib.bgzip_tabix(None, output_gz),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output_gz, None),
                        htslib.bgzip_tabix(None, output_vt_gz),
                    ]),
                    pipe_jobs([
                        bcftools.view(output_vt_gz, None, config.param('merge_filter_mutect2', 'filter_options')),
                        htslib.bgzip_tabix(None, output_somatic_vt),
                    ]),
                ], name="merge_filter_mutect2." + tumor_pair.name))

        return jobs

    def samtools_paired(self):
        """
        Samtools caller for SNVs and Indels using verison 0.1.19.
        """

        jobs = []

        nb_jobs = config.param('samtools_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            input_normal = os.path.join("alignment", tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,
                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            paired_sample = [input_normal, input_tumor]

            mkdir_job = Job(command="mkdir -p " + samtools_directory, removable_files=[samtools_directory], samples=[tumor_pair.normal, tumor_pair.tumor])

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    mkdir_job,
                    pipe_jobs([
                        bcftools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options')),
                        bcftools.call("", os.path.join(samtools_directory, tumor_pair.name + ".bcf"),
                                      config.param('samtools_paired', 'bcftools_calls_options')),
                    ]),
                    bcftools.index(os.path.join(samtools_directory, tumor_pair.name + ".bcf")),
                ], name="samtools_paired." + tumor_pair.name))

            else:
                for sequence in self.sequence_dictionary_variant():
                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            bcftools.mpileup(paired_sample, None, config.param('samtools_paired', 'mpileup_other_options'),
                                             sequence['name']),
                            bcftools.call("", os.path.join(samtools_directory, tumor_pair.name + "." + sequence['name'] + ".bcf"),
                                          config.param('samtools_paired', 'bcftools_calls_options')),
                        ]),
                        bcftools.index(os.path.join(samtools_directory, tumor_pair.name + "." + sequence['name'] + ".bcf")),
                    ], name="samtools_paired." + tumor_pair.name + "." + sequence['name']))

        return jobs

    def merge_filter_paired_samtools(self):
        """
        bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.
        The output of bcftools is fed to varfilter, which does an additional filtering of the variants
        and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls
        for all samples in the experiment.
        Additional somatic filters are performed to reduce the number of FPs: 
        1. vcflibs vcfsamplediff tags each variant with <tag>={germline,somatic,loh} to specify the type 
        of variant given the genotype difference between the two samples.
        2. bcftools filter is used to retain only variants with CLR>=15 and have STATUS=somatic from 
        vcfsamplediff
        3. bcftools filter is used to retain only variants that have STATUS=germline or STATUS=loh from
        vcfsamplediff
        """

        jobs = []
        nb_jobs = config.param('merge_filter_paired_samtools', 'approximate_nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            samtools_directory = os.path.join(pair_directory, "rawSamtools")
            output = os.path.join(samtools_directory, tumor_pair.name + ".samtools.bcf")
            output_vcf = os.path.join(pair_directory, tumor_pair.name + ".samtools.vcf.gz")
            output_vcf_vt = os.path.join(pair_directory, tumor_pair.name + ".samtools.vt.vcf.gz")
            output_somatics = os.path.join(pair_directory, tumor_pair.name + ".samtools.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".samtools.germline_loh.vt.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(samtools_directory, tumor_pair.name + ".bcf")
                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    bcftools.concat(inputs, output, config.param('merge_filter_paired_samtools', 'concat_options')),
                    pipe_jobs([
                        bcftools.view(output, None),
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix(None, output_vcf),
                    ]),
                    vt.decompose_and_normalize_mnps(output_vcf, output_vcf_vt),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None,
                                        config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        htslib.bgzip_tabix(None, output_somatics),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None,
                                        config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix(None, output_germline_loh),
                    ]),
                ], name="merge_filter_paired_samtools." + tumor_pair.name))

            else:
                inputs = [os.path.join(samtools_directory, tumor_pair.name + "." + sequence['name'] + ".bcf") for sequence in self.sequence_dictionary_variant()]

                for input_vcf in inputs:
                    if not self.is_gz_file(input_vcf):
                        stderr.write("Incomplete samtools vcf: %s\n" % input_vcf)

                jobs.append(concat_jobs([
                    Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                    bcftools.concat(inputs, output, config.param('merge_filter_paired_samtools', 'concat_options')),
                    pipe_jobs([
                        bcftools.view(output, None),
                        vcflib.vcfsamplediff(tumor_pair.normal.name, tumor_pair.tumor.name, None, None),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix(None, output_vcf),
                    ]),
                    vt.decompose_and_normalize_mnps(output_vcf, output_vcf_vt),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None,
                                        config.param('merge_filter_paired_samtools', 'somatic_filter_options')),
                        htslib.bgzip_tabix(None, output_somatics),
                    ]),
                    pipe_jobs([
                        bcftools.filter(output_vcf_vt, None,
                                        config.param('merge_filter_paired_samtools', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix(None, output_germline_loh),
                    ]),
                ], name="merge_filter_paired_samtools." + tumor_pair.name))

        report_file = os.path.join("report", "DnaSeq.merge_filter_bcf.md")
        jobs.append(
            Job(
                [output_vcf],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_template_dir}/HumanVCFformatDescriptor.tsv \\
  report/ && \\
sed 's/\t/|/g' report/HumanVCFformatDescriptor.tsv | sed '2i-----|-----' >> {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_filter_bcf_report")
        )

        return jobs

    def vardict_paired(self):
        """
        vardict caller for SNVs and Indels.
        Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is AUPAC nomenclature 
        """

        ##TO DO - the BED system needs to be revisted !! 
        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', type='filepath')

        bed_file_list = []
        if use_bed:
            bed = self.samples[0].readsets[0].beds[0]
            bed_intervals, interval_size = bed_file.parse_bed_file(bed)
            last_bed_file = 'vardict.tmp.' + str(nb_jobs - 1) + '.bed'
            if not os.path.exists(last_bed_file):
                bed_file_list = bed_file.split_by_size(bed_intervals, interval_size, nb_jobs, output="./vardict.tmp")
            else:
                for idx in range(nb_jobs):
                    bed_file_list.append(os.path.join("vardict.tmp." + str(idx) + ".bed"))

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            input_normal = os.path.join("alignment", tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,
                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + vardict_directory, removable_files=[vardict_directory], samples=[tumor_pair.normal, tumor_pair.tumor])

            if use_bed:
                idx = 0
                for bf in bed_file_list:
                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            vardict.paired_java(input_normal, input_tumor, tumor_pair.name, None, bf),
                            vardict.testsomatic(None, None),
                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                            htslib.bgzip_tabix(None, output),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + "." + str(idx)))
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))
                if nb_jobs == 1:
                    bedjob = vardict.dict2beds(genome_dictionary, beds)
                    output = os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")
                    jobs.append(concat_jobs([
                        mkdir_job,
                        bedjob,
                        pipe_jobs([
                            vardict.paired_java(input_normal, input_tumor, tumor_pair.name, None, beds.pop()),
                            vardict.testsomatic(None, None),
                            vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                            htslib.bgzip_tabix(None, output),
                        ]),
                    ], name="vardict_paired." + tumor_pair.name + ".0"))
                else:
                    bedjob = vardict.dict2beds(genome_dictionary, beds)
                    jobs.append(concat_jobs([mkdir_job, bedjob], name="vardict.genome.beds." + tumor_pair.name))

                    for idx in range(nb_jobs):
                        output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
                        jobs.append(concat_jobs([
                            mkdir_job,
                            pipe_jobs([
                                vardict.paired_java(input_normal, input_tumor, tumor_pair.name, None, beds[idx]),
                                vardict.testsomatic(None, None),
                                vardict.var2vcf(None, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                                htslib.bgzip_tabix(None, output),
                            ]),
                        ], name="vardict_paired." + tumor_pair.name + "." + str(idx)))
        return jobs

    def merge_filter_paired_vardict(self):
        """
        The fully merged vcf is filtered using following steps:
        1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
        2. Somatics identified in step 1 must have PASS filter
        """

        jobs = []
        nb_jobs = config.param('vardict_paired', 'nb_jobs', type='posint')

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output_tmp = os.path.join(pair_directory, tumor_pair.name + ".vardict.tmp.vcf.gz")
            output = os.path.join(pair_directory, tumor_pair.name + ".vardict.vcf.gz")
            output_vt = os.path.join(pair_directory, tumor_pair.name + ".vardict.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")

            if nb_jobs == 1:
                inputs = os.path.join(vardict_directory, tumor_pair.name + ".0.vardict.vcf.gz")
                jobs.append(concat_jobs([
                    Job([inputs], [output_tmp], command="ln -s -f " + inputs + " " + output_tmp, samples=[tumor_pair.normal, tumor_pair.tumor]),
                    pipe_jobs([
                        Job([output_tmp], [None],
                            command="zcat {output} | awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00'"),
                        htslib.bgzip_tabix(None, output)
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output, None),
                        htslib.bgzip_tabix(None, output_vt),
                    ]),
                    pipe_jobs([
                        bcftools.view(output_vt, None,
                                      config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                        htslib.bgzip_tabix(None, output_somatic),
                    ]),
                    pipe_jobs([
                        bcftools.view(output_vt, None,
                                      config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix(None, output_germline_loh),
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
                    pipe_jobs([
                        Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                        bcftools.concat(inputVCFs, None),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}' | grep -v 'GL00'"),
                        htslib.bgzip_tabix(None, output),
                    ]),
                    pipe_jobs([
                        vt.decompose_and_normalize_mnps(output, None),
                        htslib.bgzip_tabix(None, output_vt),
                    ]),
                    pipe_jobs([
                        bcftools.view(output_vt, None,
                                      config.param('merge_filter_paired_vardict', 'somatic_filter_options')),
                        htslib.bgzip_tabix(None, output_somatic),
                    ]),
                    pipe_jobs([
                        bcftools.view(output_vt, None,
                                      config.param('merge_filter_paired_vardict', 'germline_loh_filter_options')),
                        htslib.bgzip_tabix(None, output_germline_loh),
                    ]),
                ], name="merge_filter_paired_vardict." + tumor_pair.name))

        return jobs

    def ensemble_somatic(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)

            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.somatic.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            inputs_somatic = [input_mutect2, input_vardict, input_varscan2, input_samtools]

            for input_vcf in inputs_somatic:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete ensemble vcf: %s\n" % input_vcf)

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")
            #output_flt = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.flt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[
                os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt-work"), output_ensemble], samples=[tumor_pair.normal, tumor_pair.tumor])

            rm_job = Job(command="rm -Rf " + os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt-work") + " " +
                                 output_ensemble)

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                mkdir_job,
                rm_job,
                bcbio_variation_recall.ensemble(inputs_somatic, output_ensemble,
                                                config.param('bcbio_ensemble_somatic', 'options')),
            ], name="bcbio_ensemble_somatic." + tumor_pair.name))

        return jobs

    def ensemble_germline_loh(self):
        """
        Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls
        Filter ensemble calls to retain only calls overlapping 2 or more callers
        """

        jobs = []
        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join("pairedVariants", tumor_pair.name)

            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")
            input_samtools = os.path.join(input_directory, tumor_pair.name + ".samtools.germline_loh.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline_loh.vt.vcf.gz")
            inputs_germline = [input_vardict, input_varscan2, input_samtools]

            for input_vcf in inputs_germline:
                if not self.is_gz_file(input_vcf):
                    stderr.write("Incomplete ensemble vcf: %s\n" % input_vcf)

            output_ensemble = os.path.join(paired_ensemble_directory,
                                           tumor_pair.name + ".ensemble.germline_loh.vt.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + paired_ensemble_directory, removable_files=[
                os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh.vt-work"),
                output_ensemble], samples=[tumor_pair.normal, tumor_pair.tumor])

            rm_job = Job(command="rm -Rf " + os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline_loh.vt-work") + " " +
                                 output_ensemble)

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                mkdir_job,
                rm_job,
                bcbio_variation_recall.ensemble(inputs_germline, output_ensemble,
                                                config.param('bcbio_ensemble_germline_loh', 'options')),
            ], name="bcbio_ensemble_germline_loh." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: Standard and Somatic annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            input_normal = os.path.join("alignment", tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,
                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name,
                                                  tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            for sequence in self.sequence_dictionary_variant():
                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name,
                                                       tumor_pair.name + ".ensemble.somatic.vt.annot." + sequence[
                                                           'name'] + ".vcf.gz")

                mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_somatic_variants], samples=[tumor_pair.normal, tumor_pair.tumor])

                jobs.append(concat_jobs([
                    mkdir_job,
                    gatk.variant_annotator(input_normal, input_tumor, input_somatic_variants, output_somatic_variants,
                                           intervals=[sequence['name']]),
                ], name="gatk_variant_annotator.somatic." + sequence['name'] + "." + tumor_pair.name))

        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            input_normal = os.path.join("alignment", tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join("alignment", tumor_pair.tumor.name,
                                        tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            input_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name,
                                                       tumor_pair.name + ".ensemble.germline_loh.vt.vcf.gz")

            for sequence in self.sequence_dictionary_variant():
                output_germline_loh_variants = os.path.join(ensemble_directory, tumor_pair.name,
                                                            tumor_pair.name + ".ensemble.germline_loh.vt.annot." +
                                                            sequence['name'] + ".vcf.gz")

                mkdir_job = Job(command="mkdir -p " + ensemble_directory,
                                removable_files=[output_germline_loh_variants], samples=[tumor_pair.normal, tumor_pair.tumor])

                jobs.append(concat_jobs([
                    mkdir_job,
                    gatk.variant_annotator(input_normal, input_tumor, input_germline_loh_variants,
                                           output_germline_loh_variants, intervals=[sequence['name']]),
                ], name="gatk_variant_annotator.germline." + sequence['name'] + "." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_somatic(self):
        """
        Merge annotated somatic vcfs
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            inputs_somatic = [os.path.join(ensemble_directory, tumor_pair.name,
                                           tumor_pair.name + ".ensemble.somatic.vt.annot." + sequence[
                                               'name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]
            ouputs_somatic = os.path.join(ensemble_directory, tumor_pair.name,
                                          tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(inputs_somatic, None),
                    htslib.bgzip_tabix(None, ouputs_somatic),
                ]),
            ], name="merge_gatk_variant_annotator.somatic." + tumor_pair.name))

        return jobs

    def merge_gatk_variant_annotator_germline(self):
        """
        Merge annotated germline and LOH vcfs
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            germline_inputs = [os.path.join(ensemble_directory, tumor_pair.name,
                                            tumor_pair.name + ".ensemble.germline_loh.vt.annot." + sequence[
                                                'name'] + ".vcf.gz") for sequence in self.sequence_dictionary_variant()]
            germline_output = os.path.join(ensemble_directory, tumor_pair.name,
                                           tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz")

            jobs.append(concat_jobs([
                Job(samples=[tumor_pair.normal, tumor_pair.tumor]),
                pipe_jobs([
                    bcftools.concat(germline_inputs, None),
                    htslib.bgzip_tabix(None, germline_output),
                ]),
            ], name="merge_gatk_variant_annotator.germline." + tumor_pair.name))

        return jobs

    def compute_cancer_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        if not os.path.exists(ensemble_directory):
            os.makedirs(ensemble_directory)

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            if not os.path.exists(paired_directory):
                os.makedirs(paired_directory)

            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf")
            output_somatic_gz = os.path.join(paired_directory,
                                             tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_somatic], samples=[tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                snpeff.compute_effects(input_somatic, output_somatic, cancer_sample_file=cancer_pair_filename,
                                       options=config.param('compute_cancer_effects_somatic', 'options')),
                htslib.bgzip_tabix(output_somatic, output_somatic_gz),
            ], name="compute_cancer_effects_somatic." + tumor_pair.name))

        return jobs

    def compute_cancer_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)

            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz")
            output_germline = os.path.join(paired_directory,
                                           tumor_pair.name + ".ensemble.germline_loh.vt.annot.snpeff.vcf")
            output_germline_gz = os.path.join(paired_directory,
                                              tumor_pair.name + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            mkdir_job = Job(command="mkdir -p " + paired_directory, removable_files=[output_germline], samples=[tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                snpeff.compute_effects(input_germline, output_germline,
                                       options=config.param('compute_cancer_effects_germline', 'options')),
                htslib.bgzip_tabix(output_germline, output_germline_gz),
            ], name="compute_cancer_effects_germline." + tumor_pair.name))

        return jobs

    def sample_gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)


            jobs.append(concat_jobs([
                Job(command="mkdir -p " + paired_directory, samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                                          gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)
            ], name="gemini_annotations.somatic." + tumor_pair.name))

        return jobs

    def sample_gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """
        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + paired_directory, samples=[tumor_pair.normal, tumor_pair.tumor]),
                gemini.gemini_annotations(gemini_prefix + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz",
                                          gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)
            ], name="gemini_annotations.germline_loh." + tumor_pair.name))

        return jobs

    def set_somatic_and_actionable_mutations(self):
        """

        """       

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        ped_file = config.param('set_somatic_and_actionable_mutations', 'ped_file', required=False, type='filepath')
        ped_job = None

        for tumor_pair in self.tumor_pairs.itervalues():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)

            if not ped_file:
                ped_job = self.build_ped_file(paired_directory, tumor_pair)
                ped_file = os.path.join(gemini_prefix + ".ped")
               
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + paired_directory, samples=[tumor_pair.normal, tumor_pair.tumor]),
                ped_job,
                gemini.set_somatic(ped_file, gemini_prefix + ".somatic.gemini." + gemini_version + ".db", 
                                   gemini_prefix + ".somatic.gemini.set_somatic.tsv"),
                gemini.actionable_mutations(gemini_prefix + ".somatic.gemini." + gemini_version + ".db",
                                   gemini_prefix + ".somatic.gemini.actionable.tsv")
             ], name="set_somatic_and_actionable_mutations." + tumor_pair.name))

        return jobs

    def sym_link_ensemble(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] =  [os.path.join("pairedVariants", "ensemble", tumor_pair.name, tumor_pair.name)]

            for key,input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample + ".ensemble.somatic.vt.vcf.gz", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.somatic.vt.vcf.gz.tbi", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.somatic.vt.annot.snpeff.vcf.gz", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.somatic.vt.annot.snpeff.vcf.gz.tbi", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".somatic.gemini.set_somatic.tsv", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".somatic.gemini.actionable.tsv", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.germline_loh.vt.vcf.gz", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.germline_loh.vt.vcf.gz.tbi", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                        deliverables.sym_link_pair(sample + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz.tbi", tumor_pair, type="ensemble", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_ensemble." + tumor_pair.name + "." + key))

        return jobs

    def combine_tumor_pairs_somatic(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [
            os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for
            tumor_pair in self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor])

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                mkdir_job,
                Job([input_merged_vcfs[0]], [output],
                    command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)
            ], name="gatk_combine_variants.somatic.allPairs"))

        else:

            jobs.append(concat_jobs([
                mkdir_job,
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="gatk_combine_variants.somatic.allPairs"))

        return jobs

    def combine_tumor_pairs_germline(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name,
                                          tumor_pair.name + ".ensemble.germline_loh.vt.annot.vcf.gz") for tumor_pair in
                             self.tumor_pairs.itervalues()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor])

        if len(input_merged_vcfs) == 1:
            jobs.append(concat_jobs([
                mkdir_job,
                Job([input_merged_vcfs[0]], [output],
                    command="ln -s -f " + os.path.abspath(input_merged_vcfs[0]) + " " + output)
            ], name="gatk_combine_variants.germline_loh.allPairs"))

        else:

            jobs.append(concat_jobs([
                mkdir_job,
                gatk.combine_variants(input_merged_vcfs, output)
            ], name="gatk_combine_variants.germline_loh.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_somatic(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor])

        jobs.append(concat_jobs([
            mkdir_job,
            vt.decompose_and_normalize_mnps(input, output)
        ], name="decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def decompose_and_normalize_mnps_germline(self):
        """
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt)
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.annot.vcf.gz")
        output_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input_vcf, output_vcf)
        job.name = "decompose_and_normalize_mnps.germline.allPairs"

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output_vcf], samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor])

        jobs.append(concat_jobs([
            mkdir_job,
            vt.decompose_and_normalize_mnps(input_vcf, output_vcf)
        ], name="decompose_and_normalize_mnps.somatic.allPairs"))

        return jobs

    def all_pairs_compute_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        cancer_pair = open(cancer_pair_filename, 'w')

        for tumor_pair in self.tumor_pairs.itervalues():
            cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor])

        jobs.append(concat_jobs([
            mkdir_job,
            snpeff.compute_effects(input, output, cancer_sample_file=cancer_pair_filename,
                                   options=config.param('compute_cancer_effects_somatic', 'options')),
            htslib.bgzip_tabix(output, output_gz),
        ], name="compute_effects.somatic.allPairs"))

        return jobs

    def all_pairs_compute_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline_loh.vt.annot.snpeff.vcf.gz")

        mkdir_job = Job(command="mkdir -p " + ensemble_directory, removable_files=[output], samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor])

        jobs.append(concat_jobs([
            mkdir_job,
            snpeff.compute_effects(input, output, options=config.param('compute_cancer_effects_germline', 'options')),
            htslib.bgzip_tabix(output, output_gz),
        ], name="compute_effects.germline.allPair"))

        return jobs

    def gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + ensemble_directory, samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor]),
            gemini.gemini_annotations(gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                                      gemini_prefix + ".somatic.gemini." + gemini_version + ".db", temp_dir)
        ], name="gemini_annotations.somatic.allPairs"))

        return jobs

    def gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html
        """

        jobs = []

        ensemble_directory = os.path.join("pairedVariants", "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + ensemble_directory, samples=[self.tumor_pairs.normal, self.tumor_pairs.tumor]),
            gemini.gemini_annotations(gemini_prefix + ".ensemble.germline_loh.vt.annot.snpeff.vcf.gz",
                                      gemini_prefix + ".germline_loh.gemini." + gemini_version + ".db", temp_dir)
        ], name="gemini_annotations.germline.allPairs"))

        return jobs

    def sequenza(self):
        """
        Sequenza is a novel set of tools providing a fast python script to genotype cancer samples,
        and an R package to estimate cancer cellularity, ploidy, genome wide copy number profile and infer
        for mutated alleles.

        """
        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            sequenza_directory = os.path.join(pair_directory, "sequenza")
            rawSequenza_directory = os.path.join(sequenza_directory, "rawSequenza")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            mkdir_job = Job(command="mkdir -p " + rawSequenza_directory, removable_files=[rawSequenza_directory], samples=[tumor_pair.normal, tumor_pair.tumor])

            for sequence in self.sequence_dictionary_variant():
                normal_mpileup = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.normal.name + "." + sequence['name'] + ".mpileup")
                tumor_mpileup = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup")
                normal_gz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.normal.name + "." + sequence['name'] + ".mpileup.gz")
                tumor_gz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.tumor.name + "." + sequence['name'] + ".mpileup.gz")
                out_seqz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + "." + sequence['name'] + ".seqz.gz")
                binned_seqz = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".binned.seqz." + sequence['name'] + ".gz")

                if os.path.isfile(normal_mpileup) and os.path.isfile(tumor_mpileup):

                    jobs.append(concat_jobs([
                        mkdir_job,
                        Job([normal_mpileup], [normal_gz], command="gzip -cf " + normal_mpileup + " > " + normal_gz),
                        Job([tumor_mpileup], [tumor_gz], command="gzip -cf " + tumor_mpileup + " > " + tumor_gz),
                        pipe_jobs([
                            sequenza.sequenza_seqz(normal_gz, tumor_gz, config.param('sequenza', 'gc_file'), None),
                            Job([None], [out_seqz], command="gzip -cf > " + out_seqz)
                        ]),
                        pipe_jobs([
                            sequenza.sequenza_bin(out_seqz, None),
                            Job([None], [binned_seqz], command="gzip -c > " + binned_seqz),
                        ]),
                    ], name="sequenza.create_seqz." + sequence['name'] + "." + tumor_pair.name))

                else:

                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            samtools.mpileup([inputNormal], None, config.param('sequenza', 'mpileup_options'), sequence['name']),
                            Job([None], [normal_gz], command="gzip -cf > " + normal_gz),
                        ]),
                        pipe_jobs([
                            samtools.mpileup([inputTumor], None, config.param('sequenza', 'mpileup_options'), sequence['name']),
                            Job([None], [tumor_gz], command="gzip -cf > " + tumor_gz),
                        ]),
                    ], name="mpileup_sequenza." + sequence['name'] + "." + tumor_pair.name))

                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            sequenza.sequenza_seqz(normal_gz, tumor_gz, config.param('sequenza', 'gc_file'), None),
                            Job([None], [out_seqz], command="gzip -c > " + out_seqz),
                        ]),
                        pipe_jobs([
                            sequenza.sequenza_bin(out_seqz, None),
                            Job([None], [binned_seqz], command="gzip -c > " + binned_seqz),
                        ]),
                    ], name="sequenza.create_seqz." + sequence['name'] + "." + tumor_pair.name))

            seqz_outputs = [os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".binned.seqz." + sequence['name'] + ".gz") for sequence
                            in self.sequence_dictionary_variant()]
            seqz_input = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".binned.seqz.1.gz")
            header = os.path.join(sequenza_directory, tumor_pair.name + ".header.gz")
            tmp_output = os.path.join(sequenza_directory, tumor_pair.name + ".tmp.seqz.gz")
            merged_seqz = os.path.join(sequenza_directory, tumor_pair.name + ".binned.merged.seqz.gz")

            jobs.append(concat_jobs([
                mkdir_job,
                Job([seqz_input], [header], command="zcat " + seqz_input + " | head -1 | gzip -cf > " + header, removable_files=[header]),
                Job(seqz_outputs, [tmp_output],
                    command="zcat \\\n  " + " \\\n  ".join(seqz_outputs) + " \\\n  | grep -v 'chrom' | grep -v 'MT' | gzip -cf > " + tmp_output,
                    removable_files=[tmp_output]),
                Job([tmp_output], [merged_seqz], command="zcat " + header + " " + tmp_output + " | gzip -cf > " + merged_seqz),
            ], name="sequenza.merge_binned_seqz." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                sequenza.sequenza_main(merged_seqz, sequenza_directory, tumor_pair.name),
            ], name="sequenza." + tumor_pair.name))

        return jobs

    def sym_link_sequenza(self):
        jobs = []

        inputs = dict()

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("pairedVariants", tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory, "sequenza", tumor_pair.name + "_chromosome_view.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_genome_view.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CN_bars.pdf"),
                               os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CP_contours.pdf")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="cnv", sample=key,
                                                   profyle=self.args.profyle),
                    ], name="sym_link_fastq.report." + tumor_pair.name + "." + key))

        return jobs

    def sCNAphase(self):
        """


        """
        jobs = []
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("pairedVariants", tumor_pair.name)
            scnaphase_directory = os.path.join(pair_directory, "sCNAphase")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name,
                                       tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name,
                                      tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            output_normal_vcf = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".vcf.gz")
            output_normal_vcf_vt = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".vt.vcf.gz")
            output_normal_vcf_vt_flt = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".vt.flt.vcf.gz")
            output_tumor_vcf = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".vcf.gz")
            output_tumor_vcf_vt = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".vt.vcf.gz")
            output_tumor_vcf_vt_flt = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".vt.flt.vcf.gz")

            mkdir_job=Job(command="mkdir -p " + scnaphase_directory, removable_files=[scnaphase_directory], samples=[tumor_pair.normal, tumor_pair.tumor]),

            nb_jobs = config.param('samtools_single', 'nb_jobs', type='posint')

            if nb_jobs == 1:
                jobs.append(concat_jobs([
                    mkdir_job,
                    pipe_jobs([
                        samtools.mpileup([inputNormal], None, config.param('samtools_single', 'mpileup_other_options'),
                                         ini_section="samtools_paired"),
                        samtools.bcftools_call_pair("-",
                                                    os.path.join(scnaphase_directory, tumor_pair.normal.name + ".bcf"),
                                                    config.param('samtools_single', 'bcftools_view_options')),
                    ]),
                ], name="samtools_single." + tumor_pair.normal.name))

                jobs.append(concat_jobs([
                    mkdir_job,
                    pipe_jobs([
                        samtools.mpileup([inputTumor], None, config.param('samtools_single', 'mpileup_other_options'),
                                         ini_section="samtools_paired"),
                        samtools.bcftools_call_pair("-",
                                                    os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".bcf"),
                                                    config.param('samtools_single', 'bcftools_view_options')),
                    ]),
                ], name="samtools_single." + tumor_pair.tumor.name))

            else:
                for region in self.generate_approximate_windows(
                        nb_jobs):  # for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            samtools.mpileup([inputNormal], None,
                                             config.param('samtools_single', 'mpileup_other_options'), region,
                                             ini_section="samtools_paired"),
                            samtools.bcftools_call_pair("-", os.path.join(scnaphase_directory,
                                                                          tumor_pair.normal.name + "." + region + ".bcf"),
                                                        config.param('samtools_single', 'bcftools_view_options')),
                        ]),
                    ], name="samtools_single." + tumor_pair.normal.name + "." + region))

                    jobs.append(concat_jobs([
                        mkdir_job,
                        pipe_jobs([
                            samtools.mpileup([inputTumor], None,
                                             config.param('samtools_single', 'mpileup_other_options'), region,
                                             ini_section="samtools_paired"),
                            samtools.bcftools_call_pair("-", os.path.join(scnaphase_directory,
                                                                          tumor_pair.tumor.name + "." + region + ".bcf"),
                                                        config.param('samtools_single', 'bcftools_view_options')),
                        ]),
                    ], name="samtools_single." + tumor_pair.tumor.name + "." + region))

                inputsNormal = [os.path.join(scnaphase_directory, tumor_pair.normal.name + "." + region + ".bcf") for
                                region in
                                self.generate_approximate_windows(nb_jobs)]
                jobs.append(concat_jobs([
                    pipe_jobs([
                        samtools.bcftools_cat_pair(inputsNormal, None),
                        samtools.bcftools_view_pair("-", None),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix(None, output_normal_vcf),
                    ]),
                    vt.decompose_and_normalize_mnps(output_normal_vcf, output_normal_vcf_vt),
                    pipe_jobs([
                        Job([output_normal_vcf_vt], [None],
                            command="zgrep -Pv '\\tN\\t' " + output_normal_vcf_vt + " | grep -v 'INDEL' | grep -v '\.\/' | grep -v '\/\.' "),
                        htslib.bgzip_tabix(None, output_normal_vcf_vt_flt),
                    ]),
                ], name="merge_samtools_single." + tumor_pair.normal.name))

                inputsTumor = [os.path.join(scnaphase_directory, tumor_pair.tumor.name + "." + region + ".bcf") for
                               region in
                               self.generate_approximate_windows(nb_jobs)]
                jobs.append(concat_jobs([
                    pipe_jobs([
                        samtools.bcftools_cat_pair(inputsTumor, None),
                        samtools.bcftools_view_pair("-", None),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"),
                        Job([None], [None],
                            command="awk -F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"),
                        Job([None], [None], command="awk -F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"),
                        htslib.bgzip_tabix(None, output_tumor_vcf),
                    ]),
                    vt.decompose_and_normalize_mnps(output_tumor_vcf, output_tumor_vcf_vt),
                    pipe_jobs([
                        Job([output_tumor_vcf_vt], [None],
                            command="zgrep -Pv '\\tN\\t' " + output_tumor_vcf_vt + " | grep -v 'INDEL' | grep -v '\.\/' | grep -v '\/\.' "),
                        htslib.bgzip_tabix(None, output_tumor_vcf_vt_flt),
                    ]),
                ], name="merge_samtools_single." + tumor_pair.tumor.name))

                shapeit_nprefix = os.path.join(scnaphase_directory, tumor_pair.normal.name + ".")
                shapeit_tprefix = os.path.join(scnaphase_directory, tumor_pair.tumor.name + ".")

                for chr in range(1, 23):
                    jobs.append(concat_jobs([
                        htslib.tabix_split(output_normal_vcf_vt_flt,
                                           os.path.join(shapeit_nprefix + "chr" + str(chr) + ".vcf"), str(chr)),
                        htslib.tabix_split(output_tumor_vcf_vt_flt,
                                           os.path.join(shapeit_tprefix + "chr" + str(chr) + ".vcf"), str(chr)),
                    ], name="tabix_split." + tumor_pair.name + "." + str(chr)))

                for chr in range(1, 23):
                    jobs.append(concat_jobs([
                        shapeit.check(os.path.join(shapeit_nprefix + "chr" + str(chr) + ".vcf"),
                                      os.path.join(shapeit_nprefix + "chr" + str(chr) + ".alignments"), str(chr)),
                    ], name="shapeit.check." + tumor_pair.normal.name + "." + str(chr)))

                    jobs.append(concat_jobs([
                        shapeit.phase(os.path.join(shapeit_nprefix + "chr" + str(chr) + ".vcf"),
                                      os.path.join(
                                          shapeit_nprefix + "chr" + str(chr) + ".alignments.snp.strand.exclude"),
                                      os.path.join(shapeit_nprefix + "chr" + str(chr)),
                                      os.path.join(shapeit_nprefix + "chr" + str(chr) + ".phase"),
                                      str(chr)),
                    ], name="shapeit.phase." + tumor_pair.normal.name + "." + str(chr)))

            cd_job = Job(command="cd " + scnaphase_directory)

            jobs.append(concat_jobs([
                mkdir_job,
                cd_job,
                scnaphase.run(tumor_pair.name, tumor_pair.normal.name, tumor_pair.tumor.name),
            ], name="scnaphase." + tumor_pair.name))

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
        for tumor_pair in self.tumor_pairs.itervalues():

            pair_directory = os.path.join("SVariants", tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")

            os.system("mkdir -p " + delly_directory)

            cancer_pair_filename = os.path.join(delly_directory, tumor_pair.name + '.tsv')
            cancer_pair = open(cancer_pair_filename, 'w')
            cancer_pair.write(tumor_pair.tumor.name + "\ttumor\n")
            cancer_pair.write(tumor_pair.normal.name + "\tcontrol\n")

            inputNormal = os.path.join("alignment", tumor_pair.normal.name,
                                       tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name,
                                      tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            tumorPair = [inputTumor, inputNormal]

            mkdir_job = Job(command="mkdir -p " + delly_directory, removable_files=[delly_directory], samples = [tumor_pair.normal, tumor_pair.tumor])

            SV_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            for sv_type in SV_types:
                output_bcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf")
                output_vcf = os.path.join(delly_directory,
                                          tumor_pair.name + ".delly." + str(sv_type) + ".somatic.flt.vcf.gz")
                output_filter_somatic_bcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(
                    sv_type) + ".somatic.pre.bcf")
                output_filter_germline_bcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(
                    sv_type) + ".germline.pre.bcf")

                jobs.append(concat_jobs([
                    mkdir_job,
                    delly.call(tumorPair, output_bcf, sv_type),
                    pipe_jobs([
                        bcftools.view(output_bcf, None, config.param('delly_call_filter_somatic', 'bcftools_options')),
                        htslib.bgzip_tabix(None, output_vcf),
                    ]),
                    #delly.filter(output_bcf, output_filter_somatic_bcf, sv_type,
                    #             config.param('delly_call_filter_somatic', 'type_options'),
                    #             config.param('delly_call_filter_somatic', sv_type + '_options'),
                    #             sample_file=cancer_pair_filename),
                    #delly.filter(output_bcf, output_filter_germline_bcf, sv_type,
                    #             config.param('delly_call_filter_germline', 'type_options'),
                    #             config.param('delly_call_filter_germline', sv_type + '_options'))
                ], name="delly_call_filter." + str(sv_type) + "." + tumor_pair.name))

        return jobs

    def delly_sv_annotation(self):
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():

            pair_directory = os.path.join("SVariants", tumor_pair.name)
            final_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")
            output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly.merge.sort.vcf.gz")
            output_flt_vcf = os.path.join(pair_directory, tumor_pair.name + ".delly.merge.sort.flt.vcf.gz")
            
            SV_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            inputBCF = []
            for sv_type in SV_types:
                inputBCF.append(os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf"))

            jobs.append(concat_jobs([
                pipe_jobs([
                    bcftools.concat(inputBCF, None, "-O v"),
                    vt.sort("-", "-", "-m full"),
                    htslib.bgzip(None, output_vcf),
                ]),
                pipe_jobs([
                    bcftools.view(output_vcf, None, "-f PASS"),
                    htslib.bgzip(None, output_flt_vcf),
                ]),
            ], name="sv_annotation.delly.merge_sort_filter." + tumor_pair.name))

            jobs.append(concat_jobs([
                vawk.somatic(output_flt_vcf, tumor_pair.normal.name, tumor_pair.tumor.name, final_directory + ".delly.somatic.vcf"),
                snpeff.compute_effects(final_directory + ".delly.somatic.vcf", final_directory + ".delly.somatic.snpeff.vcf"),
                annotations.structural_variants(final_directory + ".delly.somatic.snpeff.vcf",
                                                final_directory + ".delly.somatic.snpeff.annot.vcf"),
                vawk.sv(final_directory + ".delly.somatic.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "DELLY",
                        final_directory + ".delly.somatic.prioritize.tsv"),
            ], name="sv_annotation.delly.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                vawk.germline(output_flt_vcf, tumor_pair.normal.name, tumor_pair.tumor.name, final_directory + ".delly.germline.vcf"),
                snpeff.compute_effects(final_directory + ".delly.germline.vcf", final_directory + ".delly.germline.snpeff.vcf"),
                annotations.structural_variants(final_directory + ".delly.germline.snpeff.vcf", final_directory + ".delly.germline.snpeff.annot.vcf"),
                vawk.sv(final_directory + ".delly.germline.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "DELLY",
                        final_directory + ".delly.germline.prioritize.tsv"),
            ], name="sv_annotation.delly.germline." + tumor_pair.name))
            
        return jobs
        
    def manta_pair_sv_calls(self):
        """
        Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
        analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
        Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
        single efficient workflow.
        Returns:Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
         in VCF 4.1 format.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            manta_directory = os.path.abspath(os.path.join(pair_directory, "rawManta"))
            output_prefix = os.path.abspath(os.path.join(pair_directory, tumor_pair.name))

            mkdir_job = Job(command="mkdir -p " + manta_directory, samples = [tumor_pair.normal, tumor_pair.tumor])

            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            manta_somatic_output = os.path.join(manta_directory, "results/variants/somaticSV.vcf.gz")
            manta_somatic_output_tbi = os.path.join(manta_directory, "results/variants/somaticSV.vcf.gz.tbi")
            manta_germline_output = os.path.join(manta_directory, "results/variants/diploidSV.vcf.gz")
            manta_germline_output_tbi = os.path.join(manta_directory, "results/variants/diploidSV.vcf.gz.tbi")

            output_dep = [manta_somatic_output, manta_somatic_output_tbi, manta_germline_output, manta_germline_output_tbi]

            jobs.append(concat_jobs([
                mkdir_job,
                manta.manta_pair_config(inputNormal, inputTumor, manta_directory),
                manta.manta_run(manta_directory, output_dep=output_dep),
                Job([manta_somatic_output], [output_prefix + ".manta.somatic.vcf.gz"],
                    command="ln -sf " + manta_somatic_output + " " + output_prefix + ".manta.somatic.vcf.gz"),
                Job([manta_somatic_output_tbi], [output_prefix + ".manta.somatic.vcf.gz"],
                    command="ln -sf " + manta_somatic_output_tbi + " " + output_prefix + ".manta.somatic.vcf.gz.tbi"),
                Job([manta_germline_output], [output_prefix + ".manta.germline.vcf.gz"],
                    command="ln -sf " + manta_germline_output + " " + output_prefix + ".manta.germline.vcf.gz"),
                Job([manta_germline_output_tbi], [output_prefix + ".manta.germline.vcf.gz"],
                    command="ln -sf " + manta_germline_output_tbi + " " + output_prefix + ".manta.germline.vcf.gz.tbi"),
            ], name="manta_pair_sv_calls." + tumor_pair.name))

        return jobs

    def manta_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(pair_directory + ".manta.somatic.vcf.gz", pair_directory + ".manta.somatic.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".manta.somatic.snpeff.vcf", pair_directory + ".manta.somatic.snpeff.annot.vcf"),
                vawk.sv(pair_directory + ".manta.somatic.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "MANTA",
                        pair_directory + ".manta.somatic.prioritize.tsv"),
            ], name="sv_annotation.manta_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(pair_directory + ".manta.germline.vcf.gz", pair_directory + ".manta.germline.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".manta.germline.snpeff.vcf", pair_directory + ".manta.germline.snpeff.annot.vcf"),
                vawk.sv(pair_directory + ".manta.germline.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "MANTA",
                        pair_directory + ".manta.germline.prioritize.tsv")
            ], name="sv_annotation.manta_germline." + tumor_pair.name))

        return jobs

    def sym_link_manta(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.abspath(os.path.join("SVariants", tumor_pair.name, tumor_pair.name))
            inputs["Tumor"] = [os.path.join(pair_directory + ".manta.somatic.snpeff.annot.vcf"),
                               pair_directory + ".manta.somatic.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_somatic", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_manta.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] = [os.path.join(pair_directory + ".manta.germline.snpeff.annot.vcf"),
                               pair_directory + ".manta.germline.prioritize.tsv"]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_germline", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_manta.germline." + tumor_pair.name + "." + key))

        return jobs

    def lumpy_paired_sv(self):
        """
        A probabilistic framework for structural variant discovery.
        Lumpy traditional with paired ends and split reads on tumor normal pair.
        Returns:bams.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            lumpy_directory = os.path.join(pair_directory, "rawLumpy")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            discordants_normal = os.path.join(lumpy_directory, tumor_pair.normal.name + ".discordants.sorted.bam")
            discordants_tumor = os.path.join(lumpy_directory, tumor_pair.tumor.name + ".discordants.sorted.bam")

            splitters_tumor = os.path.join(lumpy_directory, tumor_pair.tumor.name + ".splitters.sorted.bam")
            splitters_normal = os.path.join(lumpy_directory, tumor_pair.normal.name + ".splitters.sorted.bam")

            output_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf")
            gzip_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf.gz")

            genotype_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf")

            mkdir_job = Job(command="mkdir -p " + lumpy_directory, samples = [tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                pipe_jobs([
                    samtools.view(inputNormal, None, "-b -F 1294"),
                    sambamba.sort("/dev/stdin", discordants_normal, lumpy_directory, config.param('extract_discordant_reads', 'discordants_sort_option')),
                ]),
                pipe_jobs([
                    samtools.view(inputTumor, None, "-b -F 1294"),
                    sambamba.sort("/dev/stdin", discordants_tumor, lumpy_directory, config.param('extract_discordant_reads', 'discordants_sort_option')),
                ]),
            ], name="extract_discordant_reads." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                pipe_jobs([
                    samtools.view(inputNormal, None, "-h"),
                    Job([None], [None], [['lumpy_sv', 'module_lumpy']], command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"),
                    samtools.view("-", None, " -Sb "),
                    sambamba.sort("/dev/stdin", splitters_normal, lumpy_directory, config.param('extract_split_reads', 'split_sort_option')),
                ]),
                pipe_jobs([
                    samtools.view(inputTumor, None, "-h"),
                    Job([None], [None], [['lumpy_sv', 'module_lumpy']], command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"),
                    samtools.view("-", None, " -Sb "),
                    sambamba.sort("/dev/stdin", splitters_tumor, lumpy_directory, config.param('extract_split_reads', 'split_sort_option')),
                ]),
            ], name="extract_split_reads." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                lumpy.lumpyexpress_pair(inputNormal, inputTumor, output_vcf, spl_normal=splitters_normal, spl_tumor=splitters_tumor,
                                        dis_normal=discordants_normal, dis_tumor=discordants_tumor),
                htslib.bgzip(output_vcf, gzip_vcf),
            ], name="lumpy_paired_sv_calls." + tumor_pair.name))

            jobs.append(concat_jobs([
                pipe_jobs([
                    Job([gzip_vcf],[None], command="zcat " + gzip_vcf + " | grep -v \"^##contig\""),
                    bcftools.annotate(None, None, config.param('lumpy_paired_sv_calls', 'header_options')),
                    vt.sort("-", os.path.join(pair_directory, tumor_pair.name + ".lumpy.sorted.vcf"), "-m full"),
                ]),
                svtyper.genotyper(inputTumor, inputNormal, os.path.join(pair_directory, tumor_pair.name + ".lumpy.sorted.vcf"), genotype_vcf),
            ], name="lumpy_paired_sv_calls.genotype." + tumor_pair.name))

        return jobs

    def lumpy_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            prefix = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            
            genotype_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.genotyped.vcf")
            somatic_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.somatic.vcf.gz")
            germline_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.germline.vcf.gz")

            jobs.append(concat_jobs([
                pipe_jobs([
                    vawk.somatic(genotype_vcf, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                    htslib.bgzip(None, somatic_vcf),
                ]),
                pipe_jobs([
                    vawk.germline(genotype_vcf, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                    htslib.bgzip(None, germline_vcf),
                ]),
            ], name="sv_annotation.lumpy.genotypes." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(somatic_vcf, prefix + ".lumpy.somatic.snpeff.vcf"),
                annotations.structural_variants(prefix + ".lumpy.somatic.snpeff.vcf", prefix + ".lumpy.somatic.snpeff.annot.vcf"),
                vawk.sv(prefix + ".lumpy.somatic.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "LUMPY",
                        prefix + ".lumpy.somatic.prioritize.tsv"),
            ], name="sv_annotation.lumpy.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(germline_vcf, prefix + ".lumpy.germline.snpeff.vcf"),
                annotations.structural_variants(prefix + ".lumpy.germline.snpeff.vcf", prefix + ".lumpy.germline.snpeff.annot.vcf"),
                vawk.sv(prefix + ".lumpy.germline.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "LUMPY",
                        prefix + ".lumpy.germline.prioritize.tsv"),
            ], name="sv_annotation.lumpy.germline." + tumor_pair.name))

        return jobs

    def sym_link_lumpy(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory + ".lumpy.somatic.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".lumpy.somatic.prioritize.tsv")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_somatic", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_lumpy.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory + ".lumpy.germline.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".lumpy.germline.prioritize.tsv")]
        
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_germline", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_lumpy.germline." + tumor_pair.name + "." + key))

        return jobs

    def wham_call_sv(self):
        """
        Wham (Whole-genome Alignment Metrics) to provide a single, integrated framework for both structural variant
        calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts
        to employ SVs in association testing.
        Returns:vcf.

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            wham_directory = os.path.join(pair_directory, "rawWham")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            output_vcf = os.path.join(wham_directory, tumor_pair.name + ".wham.vcf")
            merge_vcf = os.path.join(wham_directory, tumor_pair.name + ".wham.merged.vcf")
            genotyped_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.merged.genotyped.vcf.gz")
            somatic_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.somatic.vcf.gz")
            germline_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.germline.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + wham_directory, samples = [tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                wham.call_sv(inputTumor, inputNormal, output_vcf),
                wham.merge(output_vcf, merge_vcf),
            ], name="wham_call_sv.call_merge." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                vt.sort(merge_vcf, os.path.join(pair_directory, tumor_pair.name + ".wham.sorted.vcf"), "-m full"),
                pipe_jobs([
                    svtyper.genotyper(inputTumor, inputNormal,
                                  os.path.join(pair_directory, tumor_pair.name + ".wham.sorted.vcf"), None),
                    htslib.bgzip_tabix(None, genotyped_vcf),
                ]),
                #pipe_jobs([
                #wham.genotype(merge_vcf, inputTumor, inputNormal, None),
                #htslib.bgzip_tabix(None, genotyped_vcf),
                #]),
                pipe_jobs([
                    vawk.somatic(genotyped_vcf, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                    htslib.bgzip_tabix(None, somatic_vcf),
                ]),
                pipe_jobs([
                    vawk.germline(genotyped_vcf, tumor_pair.normal.name, tumor_pair.tumor.name, None),
                    htslib.bgzip_tabix(None, germline_vcf),
                ]),
            ], name="wham_call_sv.genotype." + tumor_pair.name))

        return jobs

    def wham_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(pair_directory + ".wham.somatic.vcf.gz", pair_directory + ".wham.somatic.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".wham.somatic.snpeff.vcf",
                                                pair_directory + ".wham.somatic.snpeff.annot.vcf"),
                vawk.sv(pair_directory + ".wham.somatic.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "WHAM",
                        pair_directory + ".wham.somatic.prioritize.tsv"),
            ], name="sv_annotation.wham.somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(pair_directory + ".wham.germline.vcf.gz", pair_directory + ".wham.germline.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".wham.germline.snpeff.vcf",
                                                pair_directory + ".wham.germline.snpeff.annot.vcf"),
                vawk.sv(pair_directory + ".wham.germline.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "WHAM",
                        pair_directory + ".wham.germline.prioritize.tsv"),
            ], name="sv_annotation.wham.germline." + tumor_pair.name))

        return jobs

    def sym_link_wham(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory + ".wham.somatic.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".wham.somatic.prioritize.tsv")]
            
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_somatic", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_wham.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            inputs["Tumor"] = [os.path.join(pair_directory + ".wham.germline.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".wham.germline.prioritize.tsv")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_germline", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_wham.germline." + tumor_pair.name + "." + key))

        return jobs

    def cnvkit_batch(self):
        """
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            cnvkit_dir = os.path.join(pair_directory, "rawCNVkit")
            inputNormal = os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            inputTumor = os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            tarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.targetcoverage.cnn")
            antitarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.antitargetcoverage.cnn")
            ref_cnn = os.path.join(cnvkit_dir, tumor_pair.name + ".reference.cnn")
            tumor_cns = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".cns")
            vcf_gz = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            metrics = os.path.join("SVariants", "cnvkit_reference")
            poolRef = os.path.join(metrics, "pooledReference.cnn")

            if os.path.isfile(poolRef):
                pool_ref_cnn = poolRef
                ref_cnn = None

            else:
                pool_ref_cnn = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            bed = []

            if coverage_bed:
                bed = coverage_bed

            else:
                bed = None

            #mutect2_vcf = os.path.join("pairedVariants", tumor_pair.name, tumor_pair.name + ".mutect2.vcf.gz")
            vardict_vcf = os.path.join("pairedVariants", tumor_pair.name,
                                      tumor_pair.name + ".vardict.germline_loh.vt.vcf.gz")

            if os.path.isfile(vardict_vcf):
                input_vcf = vardict_vcf
                normal = tumor_pair.normal.name
                tumor = tumor_pair.tumor.name

            else:
                input_vcf = None
                normal = None
                tumor = None

            mkdir_job = Job(command="mkdir -p " + cnvkit_dir, samples = [tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                cnvkit.batch(inputTumor, inputNormal, cnvkit_dir, tar_dep=tarcov_cnn, antitar_dep=antitarcov_cnn, target_bed=bed, reference=pool_ref_cnn, output_cnn=ref_cnn),
            ], name="cnvkit_batch." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                cnvkit.fix(tarcov_cnn, antitarcov_cnn, os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"), reference=pool_ref_cnn, ref_cnn=ref_cnn),
                cnvkit.segment(os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"), tumor_cns),
            ], name="cnvkit_batch.correction." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                cnvkit.call(tumor_cns, os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns")),
                pipe_jobs([
                    cnvkit.export(os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"), None, sample_id=tumor_pair.tumor.name),
                    htslib.bgzip_tabix(None, vcf_gz),
                ]),
            ], name="cnvkit_batch.call." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                cnvkit.metrics(os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                               os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"), os.path.join(metrics, tumor_pair.name + ".metrics.tsv")),
                cnvkit.scatter(os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                               os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                               os.path.join(cnvkit_dir, tumor_pair.name + ".scatter.pdf"), input_vcf, normal, tumor),
                cnvkit.diagram(os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                               os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                               os.path.join(cnvkit_dir, tumor_pair.name + ".diagram.pdf")),
            ], name="cnvkit_batch.metrics." + tumor_pair.name))

        return jobs

    def cnvkit_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(pair_directory + ".cnvkit.vcf.gz", pair_directory + ".cnvkit.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".cnvkit.snpeff.vcf",
                                                pair_directory + ".cnvkit.snpeff.annot.vcf"),
            ], name="sv_annotation.cnvkit." + tumor_pair.name))

        return jobs

    def sym_link_cnvkit(self):
        jobs = []
    
        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory + ".cnvkit.snpeff.annot.vcf")]
        
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_somatic", sample=key, profyle=self.args.profyle),
                    ], name="sym_link_cnvkit.somatic." + tumor_pair.name + "." + key))
     


    def ensemble_metasv(self):
        """

        """
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            ensemble_directory = os.path.join("SVariants", "ensemble", tumor_pair.name)

            inputTumor = os.path.join("alignment", tumor_pair.tumor.name,
                                      tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            isize_file = os.path.join("metrics", "dna", tumor_pair.tumor.name,
                                      "picard_metrics.all.metrics.insert_size_metrics")
            gatk_vcf = os.path.join("pairedVariants", "ensemble", tumor_pair.name,
                                    tumor_pair.name + ".ensemble.somatic.flt.vcf.gz")
            gatk_pass = os.path.join("pairedVariants", "ensemble", tumor_pair.name,
                                     tumor_pair.name + ".ensemble.somatic.flt.pass.vcf.gz")
            lumpy_vcf = os.path.join(pair_directory, tumor_pair.name + ".lumpy.vcf.gz")
            manta_vcf = os.path.join(pair_directory, tumor_pair.name + ".manta.somatic.vcf.gz")
            wham_vcf = os.path.join(pair_directory, tumor_pair.name + ".wham.merged.genotyped.vcf.gz")
            delly_vcf= os.path.join(pair_directory, tumor_pair.name + ".delly.merge.sort.flt.vcf.gz")
            cnvkit_vcf = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            mkdir_job = Job(command="mkdir -p " + ensemble_directory, samples = [tumor_pair.normal, tumor_pair.tumor])

            if os.path.isfile(isize_file):
                isize_mean, isize_sd = metric_tools.extract_isize(isize_file)

            else:
                raise Exception("Error " + isize_file + " does not exist. Please run metrics step\n")

            if os.path.isfile(wham_vcf):
                input_wham = wham_vcf

            else:
                input_wham = None

            if os.path.isfile(delly_vcf):
                input_delly = delly_vcf

            else:
                input_delly = None

            if os.path.isfile(gatk_vcf):
                jobs.append(concat_jobs([
                    mkdir_job,
                    vcflib.vcffilter(gatk_vcf, gatk_pass, config.param('metasv_ensemble', 'filter_pass_options')),
                ], name="metasv_ensemble.gatk_pass." + tumor_pair.name))

            else:
                gatk_pass = None
            
            
            jobs.append(concat_jobs([
                mkdir_job,
                metasv.ensemble(lumpy_vcf, manta_vcf, cnvkit_vcf, input_wham, input_delly, gatk_pass, inputTumor,
                                tumor_pair.tumor.name,
                                os.path.join(ensemble_directory, "rawMetaSV"), ensemble_directory,
                                isize_mean=str(isize_mean), isize_sd=str(isize_sd),
                                output_vcf=[]),
                # input_wham=input_wham, input_gatk=gatk_pass, output_vcf=[]),
            ], name="metasv_ensemble." + tumor_pair.name))

        return jobs

    def metasv_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            ensemble_directory = os.path.join("SVariants", "ensemble", tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(os.path.join(ensemble_directory, "variants.vcf.gz"),
                                       os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.vcf")),
                annotations.structural_variants(
                    os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.vcf"),
                    os.path.join(ensemble_directory, tumor_pair.name + ".metasv.snpeff.annot.vcf")),
            ], name="sv_annotation.metasv_ensemble." + tumor_pair.name))

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

        for tumor_pair in self.tumor_pairs.itervalues():
            sv_directory = os.path.join("SVariants", tumor_pair.name)
            scones_directory = os.path.join(sv_directory, "SCoNEs")
            inputNormal = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")]])[0]
            inputTumor = self.select_input_files(
                [[os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")],
                 [os.path.join("alignment", tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.bam")]])[0]

            bined_count_file = os.path.join(scones_directory, tumor_pair.normal.name + ".bin" + window_size + ".tsv")

            output_scones_basename = os.path.join(scones_directory,
                                                  tumor_pair.normal.name + ".bin" + window_size + "_SCoNEs")
            scones_best_model_basename = output_scones_basename + "_Model_" + config.param('scones', 'best_model',
                                                                                           required=True)
            scones_calls_file = scones_best_model_basename + "_CNVcalls.txt"
            scones_filtered_file = scones_best_model_basename + "_CNVcalls.filtered.tsv"
            scones_annotate_basename = scones_best_model_basename + "_CNVcalls.filtered.anotated"
            scones_annotate_tmp_basename = scones_best_model_basename + "_CNVcalls.filtered.tmp"

            mkdir_job = Job(command="mkdir -p " + scones_directory, samples = [tumor_pair.normal, tumor_pair.tumor])

            jobs.append(concat_jobs([
                mkdir_job,
                bvatools.bincounter(bam=inputTumor, refbam=inputNormal, out=bined_count_file, window=window_size),
            ], name="bvatools_bincounter." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                scones.scones_pair(bined_file=bined_count_file, output_basename=output_scones_basename,
                                   window=window_size)
            ], name="scones_pair." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                scones.scones_filter(scones_calls=scones_calls_file, pair_name=tumor_pair.name,
                                     output=scones_filtered_file)
            ], name="scones_filter." + tumor_pair.name))

            jobs.append(concat_jobs([
                mkdir_job,
                scones.scones_annotate(scones_calls_filtered=scones_filtered_file,
                                       output_basename=scones_annotate_basename,
                                       tmp_basename=scones_annotate_tmp_basename)
            ], name="scones_annotate." + tumor_pair.name))

        return jobs

    def svaba_assemble(self):
        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name)
            svaba_directory = os.path.join(pair_directory, "rawSvaba")
            abs_alignment = os.path.abspath("alignment")
            input_normal = os.path.join(abs_alignment, tumor_pair.normal.name,
                                        tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join(abs_alignment, tumor_pair.tumor.name,
                                       tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            somatic_input = tumor_pair.name + ".svaba.somatic.sv.vcf"
            somatic_output = os.path.join(os.path.abspath(pair_directory), tumor_pair.name + ".svaba.somatic.vcf")

            germline_input = tumor_pair.name + ".svaba.germline.sv.vcf"
            germline_output = os.path.join(os.path.abspath(pair_directory), tumor_pair.name + ".svaba.germline.vcf")

            mkdir_job = Job(command="mkdir -p " + svaba_directory, samples = [tumor_pair.normal, tumor_pair.tumor])
            cd_job = Job(command="cd " + svaba_directory)

            jobs.append(concat_jobs([
                mkdir_job,
                cd_job,
                svaba.run(input_tumor, tumor_pair.name, input_normal),
                Job([somatic_input], [somatic_output], command="sed -e 's#" + input_normal + "#" + tumor_pair.normal.name + "#g' " + somatic_input + " | "
                                                               "sed -e 's#" + input_tumor + "#" + tumor_pair.tumor.name + "#g' > " + somatic_output),
                Job([germline_input], [germline_output], command="sed -e 's#" + input_normal + "#" + tumor_pair.normal.name + "#g' " + germline_input + " | "
                                                               "sed -e 's#" + input_tumor + "#" + tumor_pair.tumor.name + "#g' > " + germline_output)
            ], name="svaba_run." + tumor_pair.name))

        return jobs

    def svaba_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)

            jobs.append(concat_jobs([
                snpeff.compute_effects(os.path.abspath(pair_directory) + ".svaba.somatic.vcf", pair_directory + ".svaba.somatic.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".svaba.somatic.snpeff.vcf", pair_directory + ".svaba.somatic.snpeff.annot.vcf"),
                vawk.sv(pair_directory + ".svaba.somatic.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "SVABA",
                        pair_directory + ".svaba.somatic.prioritize.tsv"),
            ], name="sv_annotation.svaba_somatic." + tumor_pair.name))

            jobs.append(concat_jobs([
                snpeff.compute_effects(os.path.abspath(pair_directory) + ".svaba.germline.vcf", pair_directory + ".svaba.germline.snpeff.vcf"),
                annotations.structural_variants(pair_directory + ".svaba.germline.snpeff.vcf", pair_directory + ".svaba.germline.snpeff.annot.vcf"),
                vawk.sv(pair_directory + ".svaba.germline.snpeff.annot.vcf", tumor_pair.normal.name, tumor_pair.tumor.name, "SVABA",
                        pair_directory + ".svaba.germline.prioritize.tsv"),
            ], name="sv_annotation.svaba_germline." + tumor_pair.name))

        return jobs

    def sym_link_svaba(self):
        jobs = []

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory + ".svaba.somatic.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".svaba.somatic.prioritize.tsv")]
                               
            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_somatic", sample=key,
                                                   profyle=self.args.profyle),
                    ], name="sym_link_svaba.somatic." + tumor_pair.name + "." + key))

        inputs = dict()
        for tumor_pair in self.tumor_pairs.itervalues():
            pair_directory = os.path.join("SVariants", tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [os.path.join(pair_directory + ".svaba.germline.sv.snpeff.annot.vcf"),
                               os.path.join(pair_directory + ".svaba.germline.prioritize.tsv")]

            for key, input in inputs.iteritems():
                for sample in input:
                    jobs.append(concat_jobs([
                        deliverables.sym_link_pair(sample, tumor_pair, type="sv_germline", sample=key,
                                                   profyle=self.args.profyle),
                    ], name="sym_link_svaba.germline." + tumor_pair.name + "." + key))

        return jobs

    @property
    def steps(self):
        return [
            [
                self.picard_sam_to_fastq,
                self.sym_link_fastq_pair,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.skewer_trimming,
                self.bwa_mem_picard_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.rawmpileup_panel,
                self.paired_varscan2_panel,
                self.merge_varscan2_panel,
                self.preprocess_vcf_panel,
                self.snp_effect_panel,
                self.gemini_annotations_panel,
                self.set_somatic_and_actionable_mutations_panel,
                self.sym_link_panel,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_sambamba_flagstat,
                self.metrics_dna_fastqc,
                self.run_pair_multiqc,
                self.sym_link_report
            ],
            [
                self.picard_sam_to_fastq,
                self.sym_link_fastq_pair,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.skewer_trimming,
                self.bwa_mem_picard_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.sym_link_final_bam,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_sambamba_flagstat,
                self.metrics_dna_fastqc,
                self.run_pair_multiqc,
                self.sym_link_report,
                self.rawmpileup,
                self.paired_varscan2,
                self.merge_varscan2,
                self.paired_mutect2,
                self.merge_mutect2,
                self.samtools_paired,
                self.merge_filter_paired_samtools,
                self.vardict_paired,
                self.merge_filter_paired_vardict,
                self.ensemble_somatic,
                self.gatk_variant_annotator_somatic,
                self.merge_gatk_variant_annotator_somatic,
                self.compute_cancer_effects_somatic,
                self.sample_gemini_annotations_somatic,
                self.set_somatic_and_actionable_mutations,
                self.sym_link_ensemble,
                self.combine_tumor_pairs_somatic,
                self.decompose_and_normalize_mnps_somatic,
                self.all_pairs_compute_effects_somatic,
                self.gemini_annotations_somatic,
                self.ensemble_germline_loh,
                self.gatk_variant_annotator_germline,
                self.merge_gatk_variant_annotator_germline,
                self.compute_cancer_effects_germline,
                self.sample_gemini_annotations_germline,
                self.combine_tumor_pairs_germline,
                self.decompose_and_normalize_mnps_germline,
                self.all_pairs_compute_effects_germline,
                self.gemini_annotations_germline
            ],
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.skewer_trimming,
                self.bwa_mem_picard_sort_sam,
                self.sambamba_merge_sam_files,
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.sequenza,
                self.sCNAphase,
                self.delly_call_filter,
                self.delly_sv_annotation,
                self.manta_pair_sv_calls,
                self.manta_sv_annotation,
                self.lumpy_paired_sv,
                self.lumpy_sv_annotation,
                self.wham_call_sv,
                self.wham_sv_annotation,
                self.cnvkit_batch,
                self.cnvkit_sv_annotation,
                self.scones,
                self.ensemble_metasv,
                self.metasv_sv_annotation,
                self.svaba_assemble,
                self.svaba_sv_annotation,
                self.sym_link_sequenza,
                self.sym_link_manta,
                self.sym_link_lumpy,
                self.sym_link_wham,
                self.sym_link_svaba
            ]
        ]


if __name__ == '__main__':
    TumorPair(protocol=['fastpass','ensemble','sv'])
