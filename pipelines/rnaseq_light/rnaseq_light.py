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
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

# from bfx import bedtools
# from bfx import cufflinks
# from bfx import differential_expression
# from bfx import gq_seq_utils
# from bfx import htseq
# from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import star
# from bfx import bvatools
# from bfx import rmarkdown
from pipelines import common
import utils

from pipelines.rnaseq import rnaseq

log = logging.getLogger(__name__)

class RNAseqLight(rnaseq.RnaSeq):
	def __init__(self):
	# Add pipeline specific arguments
		self.argparser.add_argument("-d", "--design", help="design file", type=file)
		super(RNAseqLight, self).__init__()

	# def picard_sam_to_fastq_star(self):
 #      """
 #        Convert SAM/BAM files from the input readset file into FASTQ format
 #        if FASTQ files are not already specified in the readset file. Do nothing otherwise.
 #        Modified to take Star bam files
 #        """
 #        jobs = []
 #        for readset in self.readsets:
 #            # If readset FASTQ files are available, skip this step
 #            if not readset.fastq1:
 #                if readset.bam:
 #                    if readset.run_type == "PAIRED_END":
 #                        fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
 #                        fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
 #                    elif readset.run_type == "SINGLE_END":
 #                        fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
 #                        fastq2 = None
 #                    else:
 #                        raise Exception("Error: run type \"" + readset.run_type +
 #                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

 #                    job = picard.sam_to_fastq(readset.bam, fastq1, fastq2)
 #                    job.name = "picard_sam_to_fastq." + readset.name
 #                    jobs.append(job)
 #                else:
 #                    raise Exception("Error: BAM file not available for readset \"" + readset.name + "\"!")
 #        return jobs

	# 	INPUT=/home/emercier/projects/Moorehead_RNAseq_PRJBFX_1462/raw_reads/KW390ASoy200/KW390ASoy200.MPS12341978-E08.4045.6.bam FASTQ=/home/emercier/projects/Moorehead_RNAseq_PRJBFX_1462/raw_reads/KW390ASoy200/KW390ASoy200.MPS12341978-E08.4045.6.pair1.fastq.gz



	# 	input_bam=os.path.join("alignment", readset.sample.name ,readset.sample.name + ".sorted.bam")
	# 	output_fasta=[os.path.join("FILEDIRECTORY", readset.sample.name ,readset.sample.name + ".sorted.pair1.fastq.gz"), os.path.join("FILEDIRECTORY", readset.sample.name ,readset.sample.name + ".sorted.pair1.fastq.gz")]
	# 	sam_to_fastq(input_bam, fastq, second_end_fastq=None)



	def kallisto(self):
		"""
			Run Kallisto on fastq files for a fast esimate of abundance.
		"""
		transcriptome_file = config.param('kallisto', 'transcriptome', type="filepath")
		gtf_file = config.param('DEFAULT', 'gtf_transcript_id', type="filepath")

		jobs = []
		for readset in self.readsets:
			trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")

			#PAIRED
			candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
			if readset.fastq1 and readset.fastq2:
				candidate_input_files.append([readset.fastq1, readset.fastq2])
			if readset.bam:
				candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
			[fastq1, fastq2] = self.select_input_files(candidate_input_files)

			job_name = "kallisto." + readset.name
			output_dir=self.output_dir+"/kallisto/" + readset.sample.name
			job = tools.rnaseqLight_kallisto(fastq1, fastq2, transcriptome_file, gtf_file, output_dir, job_name)
			jobs.append(job)

			#SINGLE
		return jobs

	# def gq_seq_utils_exploratory_analysis_rnaseq_light(self):
	#     """
	#     Exploratory analysis using the gqSeqUtils R package adapted for RNAseqLight
	#     """

	#     jobs = []

	#     gqSeqUtils function call
	#     sample_fpkm_readcounts = [[
	#         sample.name,
	#         os.path.join("kallisto", sample.name, "abundance_transcripts.tsv"),
	#         os.path.join("kallisto", sample.name, "abundance_genes.tsv")
	#         # os.path.join("cufflinks", sample.name, "isoforms.fpkm_tracking"),
	#         # os.path.join("raw_counts", sample.name + ".readcounts.csv")
	#     ] for sample in self.samples]
	#     jobs.append(concat_jobs([
	#         Job(command="mkdir -p exploratory"),
	#         gq_seq_utils.exploratory_analysis_rnaseq(
	#             os.path.join("DGE", "rawCountMatrix.csv"),
	#             "cuffnorm",
	#             config.param('gq_seq_utils_exploratory_analysis_rnaseq', 'genes', type='filepath'),
	#             "exploratory"
	#         )
	#     ], name="gq_seq_utils_exploratory_analysis_rnaseq"))

	#     Render Rmarkdown Report
	#     jobs.append(
	#         rmarkdown.render(
	#          job_input            = os.path.join("exploratory", "index.tsv"),
	#          job_name             = "gq_seq_utils_exploratory_analysis_rnaseq_report",
	#          input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.Rmd") ,
	#          render_output_dir    = 'report',
	#          module_section       = 'report', # TODO: this or exploratory?
	#          prerun_r             = 'report_dir="report";' # TODO: really necessary or should be hard-coded in exploratory.Rmd?
	#          )
	#     )

	#     report_file = os.path.join("report", "RnaSeq.cuffnorm.md")
	#     jobs.append(
	#         Job(
	#             [os.path.join("cufflinks", "AllSamples","merged.gtf")],
	#             [report_file],
	#             command="""\
	# 			mkdir -p report && \\
	# 			zip -r report/cuffAnalysis.zip cufflinks/ cuffdiff/ cuffnorm/ && \\
	# 			cp \\
	# 			  {report_template_dir}/{basename_report_file} \\
	# 			  {report_file}""".format(
	#                 report_template_dir=self.report_template_dir,
	#                 basename_report_file=os.path.basename(report_file),
	#                 report_file=report_file
	#             ),
	#             report_files=[report_file],
	#             name="cuffnorm_report")
	#     )

	# 	return jobs


############

	@property
	def steps(self):
		return [
			self.picard_sam_to_fastq,
			self.trimmomatic,
			self.merge_trimmomatic_stats,
			self.kallisto #meme input que star
			# self.gq_seq_utils_exploratory_analysis_rnaseq_light
			]

if __name__ == '__main__':
	RNAseqLight()

