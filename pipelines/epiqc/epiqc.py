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
import sys
import csv
from operator import attrgetter


# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from pipelines import common

from bfx import bigwiginfo
from bfx import chromimpute
from bfx import wigSignalNoise
from bfx import epigeec
from bfx import epiqc_report
from bfx.readset import parse_illumina_readset_file
from bfx.design import parse_chipseq_design_file
import utils.utils
from shutil import copyfile
from bfx import genome

#from pipelines.chipseq import chipseq
from bfx import bash_cmd as bash
log = logging.getLogger(__name__)


class EpiQC(common.Illumina):
    """
        EpiQC Pipeline
        ==============

        EpiQC is a quality control pipeline for signal files (bigwig) generated from ChIP-Seq. The pipeline does a series of calculations on
        these files to determine whether or not they can be considered good. Four metrics are computed from a single bigwig file.
        With BigWigInfo it is determined if there are missing chromosomes, if the chromosome count is lower than 23 it raises a high level alert.
        ChromImpute imputes signal tracks for the first chromosome and using these imputed files EpiQC computes 2 other metrics.
        And finally the pipeline creates a heatmap from the correlation matrix obtained with EpiGeEC.

        You can test this pipeline with ChIP-Seq samples from the IHEC portal :
        https://epigenomesportal.ca/ihec/grid.html?assembly=4&build=2018-10
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)

        super(EpiQC, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {'bigwiginfo_output_directory': 'bigwiginfo',
                'chromimpute_output_directory': 'imputation',
                'bedgraph_converted_directory': 'bedgraph_data',
                'bedgraph_chr_converted_directory': 'bedgraph_chr_data',
                'chromimpute_converted_directory': 'converted',
                'chromimpute_distance_directory': 'distance',
                'chromimpute_traindata_directory': 'traindata',
                'chromimpute_predictor_directory': 'predictor',
                'chromimpute_apply': 'imputed',
                'chromimpute_eval': 'eval',
                'signal_to_noise_output_directory': 'signal_to_noise',
                'epigeec_output_directory': 'epigeec',
                'epigeec_hdf5': 'hdf5',
                'epigeec_filtered': 'filtered',
                'epigeec_output': 'output',
                'report_dir': 'report'
                }

        return dirs

    # @property
    # def contrasts(self):
    #     contrasts = super(EpiQC, self).contrasts
    #
    #     # Parse contrasts to retrieve name and type
    #     for contrast in contrasts:
    #         if re.search("^\w[\w.-]*,[BN]$", contrast.name):
    #             contrast.real_name = contrast.name.split(",")[0]
    #             if contrast.name.split(",")[1] == 'B':
    #                 contrast.type = 'broad'
    #             elif contrast.name.split(",")[1] == 'N':
    #                 contrast.type = 'narrow'
    #         else:
    #             _raise(SanitycheckError(
    #                 "Error: contrast name \"" + contrast.name + "\" is invalid (should be <contrast>,B for broad or <contrast>,N for narrow)!"))
    #
    #     return contrasts

    @property
    def inputinfo_file(self):
        inputinfo_filename = "inputinfofile.ChromImpute.txt"
        return inputinfo_filename

    @property
    def mark_type_conversion(self):
        dirs = {'N': 'narrow',
                'B': 'broad',
                'I': 'Input'
                }
        return dirs

    @property
    def readsets(self):
        flag = False
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_illumina_readset_file(self.args.readsets.name)
                for readset in self.readsets:
                    if not readset.mark_name:
                        _raise(SanitycheckError("Error: missing readset MarkName for " + readset.name))
                        flag = True
                    elif not readset.mark_type:
                        _raise(SanitycheckError("Error: missing readset MarkType for " + readset.name))
                        flag = True
                if flag:
                    exit()
            else:
                self.argparser.error("argument -r/--readsets is required!")

        return self._readsets

    @property
    def contrasts(self):
        flag = False
        if self.args.design:
            self._contrast = parse_chipseq_design_file(self.args.design.name, self.samples)
        else:
            self.argparser.error("argument -d/--design is required!")
        return self._contrast

    def createInputInfoFile(self):
        """
            Creates the input info file for ChromImpute after converting the dataset in a bedgraph format.
        """
        with open(config.param('chromimpute', 'inputinfofile'), "w+") as inputinfofile:
            for contrast in self.contrasts:
                if contrast.treatments:
                    for sample in contrast.treatments:
                        inputinfofile.write("{sample}\t{histone}\t{file}".format(
                            sample=sample.name,
                            histone=contrast.real_name,
                            file=os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                              sample.name + ".bedgraph")))

        # inputinfofile = open(config.param('chromimpute', 'inputinfofile'), "w+")
        # marks = config.param('chromimpute', 'marks') # Looks for the marks in the ini file (marks are seperated with commas)
        # marks = marks.split(",")
        # cpt = 0

        # # Convert bigwig files to bedgraph and create inputinfofile
        # for sample in self.samples:
        #     for readset in sample.readsets:
        #         if readset.bigwig: # Check if the readset has a BIGWIG column != None
        #             bigwig_file = readset.bigwig
        #             if cpt < len(marks):
        #                 mark = marks[cpt]
        #                 cpt += 1
        #         else:                      # If not, we search for the path from a chipseq pipeline
        #             prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
        #             bigwig_file = os.path.join(prefix_path, "tracks", sample.name, "bigWig", sample.name + ".bw") # Create path to bigwig file
        #             mark = config.param('DEFAULT', 'chip_type')

        #         # mark = os.path.basename(bigwig_file+".bedgraph").split(".")[-4] # Use only for IHEC database samples
        #         inputinfofile.write(sample.name+"\t"+mark+"\t"+os.path.basename(bigwig_file+".bedgraph")+"\n")

        # inputinfofile.close()

    def parseInputInfoFile(self, inputinfofile):
        """
            Parses a chromimpute input info file.
            Each element of the array samplesMarksFiles contains the name of the sample, its mark, and the name of the file in that order.
        """
        samplesMarksFiles = []
        for line in inputinfofile:
            sampleMarkFile = [line[0], line[1], line[2]]
            samplesMarksFiles.append(sampleMarkFile)

        return samplesMarksFiles

    def test(self):
        jobs = []
        S1 = {}
        for sample in self.samples:
            for readset in sample.readsets:
                # S1.setdefault(sample.name,[]).append(readset.bigwig)
                S1.setdefault(sample.name, {})[readset.bigwig] = 1
                print(sample.name)
                print(readset.bigwig)

        print ("dict")

        for key, values in S1.iteritems():
            for value in values:
                print(key, '->', value)

        return jobs

    def bigwiginfo(self):
        """
            Runs the tool bigWigInfo on bigwig files.
            It inspects the signal tracks to identify some obvious problems that
            could have an impact on the results quality. bigWigInfo[3] allows to identify things such as
            missing chromosomes and insufficient track coverage, which are usually symptoms of
            improperly generated tracks.
            If an EpiQC readset is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location
             of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

        output_dir = self.output_dirs['bigwiginfo_output_directory']

        log.debug("bigwiginfo_creating_jobs")

        #obtain samples names for only histone marks (not inputs) from design
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig: # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:                      # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(self.args.readsets.name.split("/")[:-1]) # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, "tracks", sample.name, readset.mark_name,"bigWig", sample.name +"."+ readset.mark_name+".bw") # Create path to bigwig file
                    file_ext = bigwig_file.split(".")[-1]
                    #This is not a good way. we need to accept this and write a descriptive warning in job output without raising an error the job
                    #Becuase bigwinfoinfo it self can identify the file type without the extension
                 #   if file_ext not in ("bigWig", "bw"):# != "bigWig" and file_ext != "bw":
                 #       raise Exception("Error : " + os.path.basename(bigwig_file) + " not a bigWig file !")

                    job = concat_jobs([
                            Job(command="mkdir -p " + output_dir),
                            bigwiginfo.bigWigInfo(bigwig_file, self.output_dirs['bigwiginfo_output_directory'])
                        ])
                    job.name = "bigwiginfo." + readset.name
                    job.samples = self.samples
                    jobs.append(job)

        return jobs

    def bigwig_to_bedgraph(self):
        """
            Converts bigwig files in readset to bedgraph files (Used to train ChromImpute)
            If you already have bigwig for samples create a new column and name it as "BIGWIG"
            add the corresponding bigwigs for each samples
            If you have multiple readsets paste the link for same bigwig
            If you do not have bigwigs script searches for bigwigs in chipseq output folder
            Since this pipeline is intended to run after the chipseq pipeline either way you should have bigwig files to start
            #If you have bedgraph files create a new folder
        """
        jobs = []

        # select histone mark sample names removing input chip-seq samples from desgin file
        # print(sample_name_list)
        # for contrast in self.contrasts:
        #     for sample in contrast.treatments:
        #         sample_name_list.append(sample.name)
        # print(sample_name_list)

        # /!\ gz bedgraph to save space

        #         jobs.append(Job(
        #             [],
        #             ['dataset_dir'],
        #             [],
        #             command="""\
        # mkdir -p \\
        #   {dataset_dir} \\
        #   {dataset_error}""".format(
        #     dataset_dir=self.output_dirs['bedgraph_converted_directory'],
        #     dataset_error=self.output_dirs['bedgraph_converted_directory']+"_error"),
        #             name='mkdirs_bedgraph_converted_directory'))

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:  # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(
                            self.args.readsets.name.split("/")[:-1])  # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, "tracks", sample.name, readset.mark_name, "bigWig",
                                                   sample.name + "." + readset.mark_name + ".bw")  # Create path to bigwig file
                    # checking the file extension and writing a warning in job output should be also done here
                    output_bedgraph = os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                                   sample.name + "_" + readset.mark_name + ".bedgraph")
                    output_bedgraph_gz = os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                                      sample.name + "_" + readset.mark_name + ".bedgraph.gz")

                    job = concat_jobs([
                        Job(command="mkdir -p " + self.output_dirs['bedgraph_converted_directory']),
                        bigwiginfo.bigWigToBedGraph(bigwig_file, output_bedgraph),
                        Job(output_files=[output_bedgraph_gz], command="gzip -f " + output_bedgraph)
                    ])
                    job.name = "bigwig_to_bedgraph." + sample.name + "_" + readset.mark_name
                    jobs.append(job)
                # jobs.append(bigwiginfo.bigWigToBedGraph(bigwig_file, self.output_dirs['bedgraph_converted_directory']))

        return jobs

    def chromimpute_preprocess(self):
        """
            Runs the training steps (Convert, ComputeGlobalDist, GenerateTrainData, Train) of the ChromImpute tool on the bigwig files given to the pipeline.
            This step can be skipped if the trained data from the IHEC database is used.
            All the output are stored in the imputation directory.

            Important note :
            The name of the inputinfofile and the dataset, the path to the chromsizes, the resolution and chromosome number has to be specified and the base.ini file under the chromimpute section.

            If the readset file has a BIGWIG column, we use these files for the pipeline
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        # for now there is no way we can start from bedgraph files without creating a folder names "bedgraph_info" and put all the
        # files there
        # there should be a way to define it in the readset file. or use a specific readset file for epiqc
        # for now pipeline does not support that
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)

        bedgraph_converted_files = []

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    bedgraph_converted_files.append(
                        os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                     sample.name + "_" + readset.mark_name + ".bedgraph.gz"))

        job_folder_create = Job(input_files=bedgraph_converted_files,
                                command="""\
mkdir -p \\
{output_dir}/{converteddir} \\
{output_dir}/{compute_global_dist} \\
{output_dir}/{generate_train_data} \\
{output_dir}/{train} \\
{output_dir}/{apply} \\
{output_dir}/{eval}""".format(
                                    output_dir=self.output_dirs['chromimpute_output_directory'],
                                    converteddir=self.output_dirs['chromimpute_converted_directory'],
                                    compute_global_dist=self.output_dirs['chromimpute_distance_directory'],
                                    generate_train_data=self.output_dirs['chromimpute_traindata_directory'],
                                    train=self.output_dirs['chromimpute_predictor_directory'],
                                    apply=self.output_dirs['chromimpute_apply'],
                                    eval=self.output_dirs['chromimpute_eval'])
                                )

        # copy IHEC inputinfor file from the cvmfs and paste it in the chromimpute folder
        job_copy_inputinfo = Job(output_files=[inputinfofile],
                                 command="""\
                        if test -f "{inputinfofile}"; then
                    rm {inputinfofile} 
                    fi &&
        cp {ihec_inputinfofile} {inputinfofile}""".format(
                                     inputinfofile=inputinfofile,
                                     ihec_inputinfofile="/project/6007512/C3G/projects/pubudu/epiQC_main/inputinfofile.txt")
                                 )

        # jobs.append(job_folder_create)
        # job_folder_create.name = "chromimpute_preprocess"
        job = []
        # copy histone mark, sample and file paths in user's samples into inputinfo file (avoid input histone files)
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    input_file = os.path.join(self.output_dirs['bedgraph_converted_directory'],
                                     sample.name + "_" + readset.mark_name + ".bedgraph.gz")
                    job_sample = chromimpute.modify_inputinfofile(input_file, sample.name, readset.mark_name,
                                                              inputinfofile)
                    job_sample.samples = [sample]
                    job.append(job_sample)
                    job_inputinfo = concat_jobs(job)
                # job_inputinfo.name = "chromimpute_preprocess.inputinfo"

        job_create_simlinks = Job(command="""\
    if [ "$(ls -A {output_dir}/{user_converteddir})" ]; then
    rm {output_dir}/{user_converteddir}/*
    fi &&
    ln -s {ihec_converteddir}/* {output_dir}/{user_converteddir}/""".format(
            output_dir=self.output_dirs['chromimpute_output_directory'],
            user_converteddir=self.output_dirs['chromimpute_converted_directory'],
            ihec_converteddir="/project/6007512/C3G/projects/pubudu/epiQC_main/CONVERTEDDIR")
        )

        output_dir = self.output_dirs['chromimpute_output_directory']
        chr_sizes = config.param('DEFAULT', 'chromosome_size')
        chrs = config.param('chromimpute_preprocess', 'chromosomes')
        # get the chromosome from the ini file if one chr specified it gets the chromosome and split the string in the
        # ini file
        # otherwise get all chrs information from dictionary file and creates the chromosome length file
        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', type='filepath'))
            chrs = genome.chr_names_conv(genome_dict)
        else:
            chrs = config.param('chromimpute_preprocess', 'chromosomes').split(",")

        chr_sizes_file = os.path.join(output_dir,
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))

        job_create_chr_sizes = Job(
            output_files=[chr_sizes_file],
            command="""\
        if test -f "{chr_sizes_file}"; then
        rm {chr_sizes_file} &&
        touch {chr_sizes_file}
        else
            touch {chr_sizes_file}
        fi""".format(
            chr_sizes_file=chr_sizes_file)
        )

        ##double check this code, the below jobs seems only working in chr1, not sure though
        job = []
        for chr in chrs:
            job_chr_sizes = chromimpute.generate_chr_sizes(chr_sizes_file, chr_sizes, chr)
            job.append(job_chr_sizes)
            job_chrs_sizes = concat_jobs(job)

            # do not delete below command
            # if someday you want to run convert for each chromosomes parellaly just uncomment below code

        # job_bdg_folder = Job(input_files=bedgraph_converted_files,
        #                      output_files=[os.path.join(self.output_dirs['chromimpute_output_directory'],
        #                                                self.output_dirs['bedgraph_chr_converted_directory'])],
        #                      command="""\
        #         mkdir -p {chromimpute_output_dir}/{input_dir}""".format(
        #                          chromimpute_output_dir=self.output_dirs['chromimpute_output_directory'],
        #                          input_dir=self.output_dirs['bedgraph_chr_converted_directory'])
        #                      )

        jobs.append(concat_jobs([job_folder_create, job_copy_inputinfo, job_inputinfo, job_create_simlinks,
                                 job_create_chr_sizes, job_chrs_sizes], name="chromimpute_preprocess"))
        # todo #add samples info to the above job,

        # do not delete below code
        # if someday you want to run convert for each chromosomes parellaly just uncomment below code

        # if user request not to run all the chromosomes, bedgraph files for specific chromosomes should be created
        # if(chrs!="All"):
        #     chrs = config.param('chromimpute_preprocess', 'chromosomes').split(",")
        #     job =[]
        #     outputdir = os.path.join(self.output_dirs['chromimpute_output_directory'],
        #                                self.output_dirs['bedgraph_chr_converted_directory'])
        #     for chr in chrs:
        #         for contrast in self.contrasts:
        #             for sample in contrast.treatments:
        #                 input_file = os.path.join(self.output_dirs['bedgraph_converted_directory'],
        #                                           sample.name + ".bedgraph.gz")
        #
        #                 output_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
        #                                            self.output_dirs['bedgraph_chr_converted_directory'],
        #                                           chr +"_"+ sample.name + ".bedgraph.gz")
        #
        #
        #
        #                 job_bdg_chr = chromimpute.convert_chr_bedgraph(input_file, output_file, chr, outputdir)
        #                 job_bdg_chr.samples = [sample]
        #                 job_bdg_chr.name = "bigwig_to_bedgraph.select_chr_" + sample.name + "_"+ chr
        #                 jobs.append(job_bdg_chr)
        # job_inputinfo.name = "chromimpute_preprocess.inputinfo"
        # jobs.extend([job_folder_create, job_inputinfo])
        # jobs.extend([job_folder_create, job_inputinfo])
        return jobs

    # def create_bedgraph_for_chr(self, chr):

    def chromimpute(self):
        self.bigwig_to_bedgraph()
        self.chromimpute_preprocess()

    def chromimpute_convert(self):
        """
            Creates a job for chromimpute Convert for each unique mark/ sample combination in the dataset
            If you got index out of bound exception, check whether your reference genome version of the bedgraph file is similar to chr_sizes_file
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)

        # check ini file and user request specific chromosoems change the input folder accrdingly
        chrs = config.param('chromimpute_preprocess', 'chromosomes')

        # if someday you want to run convert for each chromosomes parellaly just uncomment below code
        # if chrs=="All":
        #     input_dir = self.output_dirs['bedgraph_converted_directory']
        #     genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', type='filepath'))
        #     all_chrs = genome.chr_names_conv(genome_dict)
        # else:
        #     input_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['bedgraph_chr_converted_directory'])
        #     all_chrs = config.param('chromimpute_preprocess', 'chromosomes').split(",")

        if chrs == "All":
            genome_dict = os.path.expandvars(config.param('DEFAULT', 'genome_dictionary', type='filepath'))
            all_chrs = genome.chr_names_conv(genome_dict)
        else:
            all_chrs = config.param('chromimpute_preprocess', 'chromosomes').split(",")

        input_dir = self.output_dirs['bedgraph_converted_directory']
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_converted_directory'])

        # for contrast in self.contrasts:
        #     if contrast.treatments:
        #         output_files = []
        #         for sample in contrast.treatments:
        #             with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:
        #                 for line in chrominfofile:
        #                     output_files.append(os.path.join(output_dir, "%s_%s.bedgraph.gz.wig.gz" % (line.split("\t")[0], sample.name)))
        #                 job = chromimpute.convert(input_dir, output_dir, output_files, inputinfofile, contrast.real_name, sample.name)
        #                 jobs.append(job)
        chr_sizes_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    output_files = []
                    for chr in all_chrs:
                        output_files.append(os.path.join(output_dir, chr + "_" + sample.name + "_" + readset.mark_name +".bedgraph.gz.wig.gz"))

                    input_file = os.path.join(input_dir,
                                     sample.name + "_" + readset.mark_name + ".bedgraph.gz")

                    job = chromimpute.convert(input_dir, input_file, output_dir, output_files, inputinfofile,
                                              readset.mark_name, sample.name, chr_sizes_file)
                    job.samples = [sample]
                    jobs.append(job)
                # treatment_files = [os.path.join(self.output_dirs['bedgraph_converted_directory'], sample.name + ".bedgraph.gz") for sample in contrast.treatments]

        return jobs

    def chromimpute_compute_global_dist(self):
        """
            Creates a job for chromimpute ComputeGlobalDist for each unique mark in the inputinfo file
        """
        jobs = []
        chr_sizes_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))


        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_distance_directory'])

        histone_marks = []
        # create a job for each histone maek in inputinfor file together with user's histone marks
        # get unique histone marks from the inputinfo file
        with open(inputinfofile, "r") as histone_marks_total:
            for line in histone_marks_total:
                histone_marks.append(line.strip().split("\t")[1])
        histone_marks = list(set(histone_marks))



        # this check all the converted file for converted step as input file
        # so it takes time
        # if you want to do this uncomment this and comment the code below this
        # for histone in histone_marks:
        #     with open(inputinfofile, "r") as inputinfo:
        #         for inputinfoline in inputinfo:
        #             #if histone mark in inputinfo file present in design file it only includes as an input fileof the convert global distance step
        #             if inputinfoline.split("\t")[1]==histone:
        #
        #                 with open(chr_sizes_file, "r") as chrominfofile:
        #                     for line in chrominfofile:
        #
        #                         input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'],
        #                                     self.output_dirs['chromimpute_converted_directory'], "%s_%s.wig.gz" %
        #                                             (line.strip().split("\t")[0], inputinfoline.strip().split("\t")[2])))
        #
        #                 output_files.append(
        #                     os.path.join(output_dir, "%s_%s.txt" % (inputinfoline.split("\t")[0], histone)))
        #
        #     job = chromimpute.compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile,
        #                                           histone, chr_sizes_file)

        # this only get first line in the chromosome file(i.e chr1) as inputs to the converted step. if u want to get
        # all the chrms uncomment above comment and comment below code
        for histone in histone_marks:
            input_files = []
            output_files = []
            with open(chr_sizes_file, "r") as chrominfofile:
                line = chrominfofile.readline()

                with open(inputinfofile, "r") as inputinfo:
                    for inputinfoline in inputinfo:

                        # if histone mark in inputinfo file present in design file it only includes as an input fileof the convert global distance step
                        if inputinfoline.split("\t")[1] == histone:

                            input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'],
                                                            self.output_dirs['chromimpute_converted_directory'],
                                                            "%s_%s.wig.gz" %
                                                            (line.strip().split("\t")[0],
                                                             inputinfoline.strip().split("\t")[2])))

                            output_files.append(
                                os.path.join(output_dir, "%s_%s.txt" % (inputinfoline.split("\t")[0], histone)))

            job = chromimpute.compute_global_dist(input_files, output_dir, output_files, converteddir, inputinfofile,
                                                  histone, chr_sizes_file)
            jobs.append(job)

        return jobs

    def chromimpute_generate_train_data(self):
        """
            Creates a job for chromimpute GenerateTrainData for each unique mark in the dataset
        """
        jobs = []

        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], self.inputinfo_file)
        temp2_inputinfofile_path = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                                "temp2_" + self.inputinfo_file)
        temp_inputinfofile_path = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                               "temp_" + self.inputinfo_file)
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                   self.output_dirs['chromimpute_distance_directory'])
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_traindata_directory'])
        chr_sizes_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))
        #get all the unique histone marks in user's readset file
        histone_marks =[]
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    histone_marks.append(readset.mark_name)

        histone_marks = set(histone_marks)

        train_data_path = config.param('chromimpute_generate_train_data', 'pre_trained_data_path')

        # create temp inputinfo file with interneal indexes used in traindata step
        # this file will be used to get the integer index for output file

        # job_temp_inputinfo = chromimpute.temp_inputinfo( inputinfofile, temp_inputinfofile)
        # jobs.append(job_temp_inputinfo)

        train_user_data = config.param('DEFAULT', 'train_only_user_data')

        if train_user_data == "F":
            if os.path.exists(temp2_inputinfofile_path):
                os.remove(temp2_inputinfofile_path)

            if os.path.exists(temp_inputinfofile_path):
                os.remove(temp_inputinfofile_path)

            temp_inputinfofile = open(temp2_inputinfofile_path, "a")
            lines = open(inputinfofile, 'r').readlines()

            for line in sorted(lines, key=lambda line: line.split()[0]):
                temp_inputinfofile.write(line)
            temp_inputinfofile.close()
            temp_inputinfofile = open(temp_inputinfofile_path, "a")
            with open(temp2_inputinfofile_path, "r") as inputinfo:
                i = 1
                index_i = 0
                for inputinfoline in inputinfo:
                    if (i == 1):
                        sample = inputinfoline.split("\t")[0]
                        index_i = i - 1
                        temp_inputinfofile.write("%s\t%i%s" % (inputinfoline.strip(), index_i, os.linesep))
                        i += 1
                    else:
                        if sample == inputinfoline.split("\t")[0]:
                            temp_inputinfofile.write("%s\t%i%s" % (inputinfoline.strip(), index_i, os.linesep))
                            i += 1
                        else:
                            sample = inputinfoline.split("\t")[0]
                            index_i += 1
                            temp_inputinfofile.write("%s\t%i%s" % (inputinfoline.strip(), index_i, os.linesep))
                            i += 1
            os.remove(temp2_inputinfofile_path)
            temp_inputinfofile.close()
            if train_data_path == "":

                #        jobs.append(chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir,
                #        distancedir, inputinfofile, contrast.real_name))

                for histone in histone_marks:
                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            input_files = []
                            output_files = []
                            with open(temp_inputinfofile_path, "r") as inputinfo:
                                for inputinfoline in inputinfo:
                                    # if histone mark in inputinfo file present in readset file it only includes as an input fileof the convert global distance step
                                    if inputinfoline.split("\t")[1] == histone:
                                        input_files.append(os.path.join(converteddir,
                                                                        "%s_%s.wig.gz" %
                                                                        (chr_name,
                                                                         inputinfoline.strip().split("\t")[2])))

                                        input_files.append(os.path.join(distancedir, "%s_%s.txt" % (
                                            inputinfoline.split("\t")[0], histone)))

                                        output_files.append(
                                            os.path.join(output_dir, "%s_traindata_%s_%i_0.txt.gz" % (
                                                chr_name, histone, int(inputinfoline.strip().split("\t")[3]))))
                                        output_files.append(
                                            os.path.join(output_dir, "attributes_%s_%i_0.txt.gz" % (
                                                histone, int(inputinfoline.strip().split("\t")[3]))))

                    # input_files.append(temp_inputinfofile)
                    job = chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir,
                                                          distancedir,
                                                          inputinfofile, histone, chr_sizes_file, chr_name)
                    jobs.append(job)

            else:
                job_create_simlinks = Job(command="""\
                    if [ "$(ls -A {output_dir}/{usertraindatadir})" ]; then
                    rm {output_dir}/{usertraindatadir}/*
                    fi &&
                    ln -s {ihec_traindatadir}/* {output_dir}/{usertraindatadir}/""".format(
                    output_dir=self.output_dirs['chromimpute_output_directory'],
                    usertraindatadir=self.output_dirs['chromimpute_traindata_directory'],
                    ihec_traindatadir=train_data_path)
                )
                jobs.append(
                    concat_jobs([job_create_simlinks], name="chromimpute_generate_train_data.use_pre_trained_data"))
        else:
            log.info("Currently, only training user data is not supported. Skipping ....")
            # todo
            # mula idanma hadan enna one meka
            # user data witharak train karanna

        # chroms = []
        # with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:#with open(config.param('chromimpute', 'chrominfofile')) as chrominfofile:             os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")]
        #     for line in chrominfofile:
        #         chroms.append(line.split("\t")[0])

        # for contrast in self.contrasts:
        #     if contrast.treatments:
        #         input_files = []
        #         output_files = []
        #         indice = 0
        #         for sample in contrast.treatments:
        #             with open(os.path.join(os.environ[config.param('chromimpute', 'chromosome_size').split("/")[0].replace("$", "")], config.param('chromimpute', 'chromosome_size').replace(config.param('chromimpute', 'chromosome_size').split("/")[0]+"/", "")), "r") as chrominfofile:
        #                 for line in chrominfofile:
        #                     input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'], "%s_%s.bedgraph.gz.wig.gz" % (line.split("\t")[0], sample.name)))
        #             input_files.append(os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_distance_directory'], "%s_%s.txt" % (sample.name, contrast.real_name)))
        #             output_files.append(os.path.join(output_dir, "attributes_%s_%i_0.txt.gz" % (contrast.real_name, indice)))
        #             output_files.append(os.path.join(output_dir, "traindata_%s_%i_0.txt.gz" % (contrast.real_name, indice)))
        #             indice += 1
        #         # for chrom in chroms:
        #         jobs.append(chromimpute.generate_train_data(input_files, output_dir, output_files, converteddir, distancedir, inputinfofile, contrast.real_name))

        return jobs

    def chromimpute_train(self):
        """
            Creates a job for chromimpute Train for every sample mark combination given in the dataset
        """
        jobs = []

        temp_inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                          "temp_" + self.inputinfo_file)
        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                     self.inputinfo_file)
        traindatadir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_traindata_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_predictor_directory'])

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    input_files = []
                    output_files = []
                    with open(temp_inputinfofile, "r") as inputinfo:
                        for inputinfoline in inputinfo:
                            inputinfo_histone = inputinfoline.strip().split("\t")[1]
                            inputinfo_sample = inputinfoline.strip().split("\t")[0]
                            if readset.mark_name == inputinfo_histone:
                                if sample.name != inputinfo_sample:
                                    output_files.append(os.path.join(output_dir, "useattributes_%s_%s_%i_0.txt.gz" % (
                                        sample.name, readset.mark_name, int(inputinfoline.strip().split("\t")[3]))))

                                input_files.append(
                                    os.path.join(traindatadir, "attributes_%s_%i_0.txt.gz" % (
                                        readset.mark_name, int(inputinfoline.strip().split("\t")[3]))))

                            # adding chr name will make things complex. so haven't added it. may be later
                            # input_files.append(
                            #   os.path.join(output_dir, "%s_traindata_%s_%i_0.txt.gz" % (
                            #      chr_name, histone, int(inputinfoline.strip().split("\t")[3]))))

                    jobs.append(
                    chromimpute.train(input_files, output_dir, output_files, traindatadir, inputinfofile, sample.name,
                                      readset.mark_name))

        return jobs

    def chromimpute_train_step(self):
        """
            Runs the training steps (Convert, ComputeGlobalDist, GenerateTrainData, Train) of the ChromImpute tool on the bigwig files given to the pipeline.
            This step can be skipped if the trained data from the IHEC database is used.
            All the output are stored in the imputation directory.

            Important note :
            The resolution and chromosome name to be run has to be specified and the base.ini file under the chromimpute section.

            If the readset file has a BIGWIG column, we use these files for the pipeline
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        # TODO idea: Delete predictordir to save space after run

        jobs = []

        # inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'], "inputinfofile.ChromImpute.txt")
        # path_dataset = os.path.join(os.getcwd(), config.param('chromimpute', 'dataset'))

        # self.createInputInfoFile()
        # read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter='\t')
        # samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        log.debug("chromimpute_convert")
        jobs.extend(self.chromimpute_convert())
        log.debug("chromimpute_compute_global_dist")
        jobs.extend(self.chromimpute_compute_global_dist())
        log.debug("chromimpute_generate_train_data")
        jobs.extend(self.chromimpute_generate_train_data())
        log.debug("chromimpute_train")
        jobs.extend(self.chromimpute_train())

        return jobs

    def chromimpute_apply(self):
        """
            Creates a job for chromimpute Apply for every sample and mark given in the dataset
        """
        jobs = []

        temp_inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                          "temp_" + self.inputinfo_file)
        inputinfofile = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                     self.inputinfo_file)
        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        distancedir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                   self.output_dirs['chromimpute_distance_directory'])
        predictordir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_predictor_directory'])
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_apply'])
        chr_sizes_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    input_files = []
                    output_files = []
                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            with open(temp_inputinfofile, "r") as inputinfo:
                                for inputinfoline in inputinfo:
                                    inputinfo_histone = inputinfoline.strip().split("\t")[1]
                                    inputinfo_sample = inputinfoline.strip().split("\t")[0]
                                    if readset.mark_name == inputinfo_histone:
                                        if sample.name != inputinfo_sample:
                                            input_files.append(
                                                os.path.join(predictordir, "useattributes_%s_%s_%i_0.txt.gz" % (
                                                    sample.name, readset.mark_name,
                                                    int(inputinfoline.strip().split("\t")[3]))))

                                # output_files.append(os.path.join(output_dir, "classifier_%s_%s_%i_0.txt.gz" % (
                                #   sample.name, contrast.real_name, int(inputinfoline.strip().split("\t")[3]))))
                            output_files.append(
                                os.path.join(output_dir, "%s_impute_%s_%s.wig.gz" % (chr_name, sample.name,
                                                                                     readset.mark_name)))

                            # adding chr name will make things complex. so haven't added it. may be later
                            # input_files.append(
                            #   os.path.join(output_dir, "%s_traindata_%s_%i_0.txt.gz" % (
                            #      chr_name, histone, int(inputinfoline.strip().split("\t")[3]))))

                            jobs.append(
                                chromimpute.apply(input_files, output_dir, converteddir, distancedir, predictordir,
                                                  inputinfofile, sample.name, readset.mark_name, chr_name,
                                                  chr_sizes_file))

        return jobs

    # def chromimpute_eval(self, samplesMarksFile):
    #     """
    #         Creates a job for chromimpute Eval for every sample mark given in the dataset
    #     """
    #     jobs = []
    #
    #     input_base = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_apply']) # input_base is the name of the folder containing imputed files
    #     converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_converted_directory'])
    #     percent1 = config.param('chromimpute', 'percent1')
    #     percent2 = config.param('chromimpute', 'percent2')
    #
    #
    #     for contrast in self.contrasts:
    #         if contrast.treatments:
    #             input_dir = input_base+"_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]
    #             output_path = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs['chromimpute_eval'], "eval_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".txt")
    #             jobs.append(chromimpute.eval(input_dir, percent1, percent2, converteddir, sampleMarkFile[2]+".wig.gz", input_base, "impute_"+sampleMarkFile[0]+"_"+sampleMarkFile[1]+".wig.gz", output_path))
    #
    #     return jobs

    def chromimpute_eval(self):
        """
            Creates a job for chromimpute Eval for every sample mark given in the dataset
        """
        jobs = []

        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])
        imputeddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                  self.output_dirs['chromimpute_apply'])
        chr_sizes_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))

        percent1 = config.param('chromimpute_eval', 'percent1')
        percent2 = config.param('chromimpute_eval', 'percent2')
        output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'], self.output_dirs[
            'chromimpute_eval'])

        for sample in self.samples:
            for readset in sample.readsets:
                input_files = []
                if readset.mark_type != "I":
                    with open(chr_sizes_file, "r") as chrominfofile:
                        line = chrominfofile.readline()
                        chr_name = line.strip().split("\t")[0]

                        input_files.append(os.path.join(imputeddir, "%s_impute_%s_%s.wig.gz" % (
                        chr_name, sample.name, readset.mark_name)))

                        input_files.append(os.path.join(converteddir, "%s_%s_%s.bedgraph.gz.wig.gz" % (
                        chr_name, sample.name, readset.mark_name)))

                        imputed_file = "impute_%s_%s.wig.gz" % (
                        sample.name, readset.mark_name)

                        converted_file = "%s_%s.bedgraph.gz.wig.gz" % (
                        sample.name, readset.mark_name)

                    output_file = os.path.join(output_dir, "eval_%s_%s.tsv" % (sample.name,
                                                                               readset.mark_name))

                    jobs.append(chromimpute.eval(input_files, imputed_file, converted_file, output_file, converteddir,
                                                  imputeddir, percent1, percent2, chr_sizes_file, sample.name,
                                                  readset.mark_name))

        return jobs

    def signal_to_noise(self):
        """
            Uses the converted files obtained from chromimpute and calculates the ratio between the top 10% & top 5% of the signals over the sum the the signals
            Outputs the results in a tsv file
        """
        jobs = []
        output_dir = self.output_dirs['signal_to_noise_output_directory']
        chr_sizes_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                      "_".join((config.param('DEFAULT', 'scientific_name'), config.param('DEFAULT',
                                                                                                         'assembly'),
                                                "chrom_sizes.txt")))
   #     jobs.append()

        converteddir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                    self.output_dirs['chromimpute_converted_directory'])

        for sample in self.samples:
            for readset in sample.readsets:
                input_files = []
                if readset.mark_type != "I":
                    with open(chr_sizes_file, "r") as chrominfofile:
                        for line in chrominfofile:
                            chr_name = line.strip().split("\t")[0]
                            # If not, we search for the path from a chipseq pipeline
                            converted_bedgraph_file = os.path.join(converteddir,
                                                                   "%s_%s_%s.bedgraph.gz.wig.gz" % (chr_name, sample.name, readset.mark_name))

                            output_file = os.path.join(output_dir,
                                                       "%s_%s_%s.tsv" % (chr_name, sample.name, readset.mark_name))
                            signal_noise_folder = Job(
                                command="mkdir -p {output_dir}".format(
                                    output_dir=output_dir
                                ))
                            signal_noise_job = Job(
                                [converted_bedgraph_file],
                                [output_file],
                                [['signal_noise', 'module_python']],
                                command="""\
                            python /home/pubudu/projects/rrg-bourqueg-ad/pubudu/epiqc/epiqc_2021/genpipes/bfx/signal_noise.py \\
                              -i {input_file} \\
                              -p1 {percent1} \\
                              -p2 {percent2} \\
                              -o {output_dir}""".format(
                                    input_file=converted_bedgraph_file,
                                    percent1=config.param('signal_noise', 'percent1'),
                                    percent2=config.param('signal_noise', 'percent2'),
                                    output_dir=output_file
                                ))

                            jobs.append(concat_jobs([signal_noise_folder, signal_noise_job], name='signal_noise.' + "%s_%s_%s.tsv" % (chr_name, sample.name, readset.mark_name)))

        return jobs

   # job_chr_sizes = chromimpute.generate_chr_sizes(chr_sizes_file, chr_sizes, chr)
   # job.append(job_chr_sizes)
    #job_chrs_sizes = concat_jobs(job)
    def epigeec_tohdf5(self):
        jobs = []

        input_dir = 'epigeec'
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        job_epigeec_hdf5 =[]
        hdf5_folder = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=output_dir
            ))
        for sample in self.samples:
            for readset in sample.readsets:
                input_files = []
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:  # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(self.args.readsets.name.split("/")[:-1])  # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, "tracks", sample.name, readset.mark_name, "bigWig",
                                                   sample.name + "." + readset.mark_name + ".bw")  # Create path to bigwig file
                    file_ext = bigwig_file.split(".")[-1]


                    job_epigeec_hdf5.append(epigeec.tohdf5(output_dir, bigwig_file))
                    jobs_epigeec_hdf5 = concat_jobs(job_epigeec_hdf5)
        jobs.append(concat_jobs([hdf5_folder, jobs_epigeec_hdf5],name = "epigeec_tohdf5"))

        return jobs

    def epigeec_filter(self, hdf5_files):
        jobs = []
        filter_job = []
        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_hdf5'])
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_filtered'])
        filter_folder_job = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=output_dir
            ))



        for hdf5_file in hdf5_files:
            path_to_hdf5_file = os.path.join(input_dir, hdf5_file)
            filter_job.append(epigeec.filter(output_dir, path_to_hdf5_file))
            filter_jobs = concat_jobs(filter_job)

        jobs.append(concat_jobs([filter_folder_job, filter_jobs],name="epigeec_filter"))

        return jobs

    def epigeec_correlate(self, skip_filter_step, hdf5_files):
        jobs=[]
        job=[]
        hdf5_file_list_name = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 "hdf5_list.txt")

        input_files =[]
        output_dir = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 self.output_dirs['epigeec_output'])
        output_folder_job = Job(
            command="mkdir -p {output_dir}".format(
                output_dir=output_dir
            ))

        job_create_hdf5_list_file = Job(
            output_files=[hdf5_file_list_name],
            command="""\
        if test -f "{hdf5_file_list_name}"; then
        rm {hdf5_file_list_name} &&
        touch {hdf5_file_list_name}
        else
            touch {hdf5_file_list_name}
        fi""".format(
            hdf5_file_list_name=hdf5_file_list_name)
        )

        for hdf5_file in hdf5_files:
            if skip_filter_step:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 self.output_dirs['epigeec_hdf5'], hdf5_file)
                input_files.append(path_to_hdf5_file)
            else:
                path_to_hdf5_file = os.path.join(self.output_dirs['epigeec_output_directory'],
                                                 self.output_dirs['epigeec_filtered'], hdf5_file)
                input_files.append(path_to_hdf5_file)

            job_hdf5_list_file = epigeec.generate_hdf5_list(hdf5_file_list_name, path_to_hdf5_file)
            job.append(job_hdf5_list_file)
            jobs_hdf5_list_file = concat_jobs(job)


        job_matrix = epigeec.correlate(input_files, output_dir, hdf5_file_list_name)
        jobs.append(concat_jobs([output_folder_job, job_create_hdf5_list_file, jobs_hdf5_list_file, job_matrix], name="epigeec_matrix"))

        return jobs

    def epigeec(self):
        """
            Runs the epigeec tool on the bigwig files given to the pipeline
            Bigwig files are first converted to the hdf5 format, then the correlation matrix is computed

            If a readset file is given to the pipeline, we search for the BIGWIG column to get the files
            If epiqc is ran after a chipseq pipeline, the path to the bigwig files is reconstructed through the location of the chipseq readset file (HAS TO BE IN THE SAME FOLDER AS THE OUTPUT OF THE CHIPSEQ PIPELINE)
        """
        jobs = []

        skip_filter_step = config.param('epigeec', 'select') == '' and config.param('epigeec', 'exclude') == ''

        hdf5_files = []

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.mark_type != "I":
                    if readset.bigwig:  # Check if the readset has a BIGWIG column != None
                        bigwig_file = readset.bigwig
                    else:  # If not, we search for the path from a chipseq pipeline
                        prefix_path = "/".join(self.args.readsets.name.split("/")[:-1])  # Find path to chipseq folder
                        bigwig_file = os.path.join(prefix_path, "tracks", sample.name, readset.mark_name, "bigWig",
                                                   sample.name + "." + readset.mark_name + ".bw")  # Create path to bigwig file

                    hdf5_files.append(os.path.basename(bigwig_file) + ".hdf5")

        jobs.extend(self.epigeec_tohdf5())

        if not skip_filter_step:  # We skip this step if there are no filter files

            jobs.extend(self.epigeec_filter(hdf5_files))


        jobs.extend(self.epigeec_correlate(skip_filter_step, hdf5_files))

        return jobs

    def epiqc_report(self):
        """
            Creates a report file for each bigwig file.

            Alert levels :
                High Level Alert:
                    Chromosome count is under 23
                Medium Level Alert:
                    Whole genome bases covered under 25,000,000
                    GeEC average correlation score under 50% within the same consortium tracks
                    Signal in top 10% bins below 30%
                    ChromImpute OBSERVED_1.0_IMPUTE_5.0 below 30%
                Low Level Alert:
                    Whole genome bases covered under 75,000,000
                    Signal in top 5% bins below 20%
                    ChromImpute BOTH_1.0 below 20%
        """
        jobs = []

        jobs.append(Job(
            [],
            ['epiqc_report'],
            [],
            command=
            "mkdir \
              {output_dir}".format(
                output_dir=self.output_dirs['report_dir']),
            name="mkdir_report_epiqc"))

        chromimpute_output_dir = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                              self.output_dirs['chromimpute_eval'])

        read_inputinfofile = csv.reader(open(config.param('chromimpute', 'inputinfofile'), 'rb'), delimiter='\t')
        samplesMarksFiles = self.parseInputInfoFile(read_inputinfofile)

        marks = config.param('chromimpute', 'marks')
        marks = marks.split(",")
        cpt = 0

        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bigwig is not None:
                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'],
                                                   "bigwiginfo_" + os.path.basename(readset.bigwig) + ".txt")
                    signal_noise_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'],
                                                     os.path.basename(readset.bigwig) + ".bedgraph.wig.gz.tsv")
                    mark = marks[cpt]  # Corresponds to the marks specified in the epiqc.base.ini
                    cpt += 1
                else:
                    bigwiginfo_file = os.path.join(self.output_dirs['bigwiginfo_output_directory'],
                                                   "bigwiginfo_" + sample.name + ".bw.txt")
                    signal_noise_file = os.path.join(self.output_dirs['signal_to_noise_output_directory'],
                                                     sample.name + ".bw.bedgraph.wig.gz.tsv")
                    mark = config.param('DEFAULT', 'chip_type')

                report_file = os.path.join(self.output_dirs['report_dir'],
                                           "report_" + sample.name + "_" + mark + ".txt")

                jobs.append(Job(
                    ['epiqc_report', bigwiginfo_file],
                    [],
                    [],
                    name="report_bigwiginfo_" + os.path.basename(bigwiginfo_file),
                    command=
                    "python ../genpipes/bfx/epiqc_report.py \
                      -b {bigwiginfo_file} \
                      -cb {chromCount} \
                      -bc1 {low_alert_bases_covered} \
                      -bc2 {medium_alert_bases_covered} \
                      -o {output_file}".format(
                        bigwiginfo_file=bigwiginfo_file,
                        chromCount=config.param('epiqc_report', 'chromcount_threshold'),
                        low_alert_bases_covered=config.param('epiqc_report', 'low_alert_bases_covered'),
                        medium_alert_bases_covered=config.param('epiqc_report', 'medium_alert_bases_covered'),
                        output_file=report_file)))

                eval_file = "eval_" + sample.name + "_" + mark + ".txt"
                eval_file = os.path.join(self.output_dirs['chromimpute_output_directory'],
                                         self.output_dirs['chromimpute_eval'], eval_file)

        #                 jobs.append(Job(
        #                     ['epiqc_report', eval_file],
        #                     [],
        #                     [],
        #                     name = "report_eval_"+os.path.basename(eval_file),
        #                     command =
        # "python ../genpipes/bfx/epiqc_report.py \
        #   -c {eval_file} \
        #   -p1 {percent1} \
        #   -p2 {percent2} \
        #   -ct1 {chromimpute_threshold_M} \
        #   -ct2 {chromimpute_threshold_L} \
        #   -o {output_file}".format(
        #                     eval_file = eval_file,
        #                     percent1 = config.param('chromimpute', 'percent1'),
        #                     percent2 = config.param('chromimpute', 'percent2'),
        #                     chromimpute_threshold_M=config.param('epiqc_report', 'chromimpute_threshold_M'),
        #                     chromimpute_threshold_L=config.param('epiqc_report', 'chromimpute_threshold_L'),
        #                     output_file = report_file)))

        #             jobs.append(Job(
        #                     ['epiqc_report', signal_noise_file],
        #                     [],
        #                     [],
        #                     name = "report_signal_noise_" + os.path.basename(signal_noise_file),
        #                     command =
        # "python ../genpipes/bfx/epiqc_report.py \
        #   -s {signal_noise_file} \
        #   -s1 {percent1} \
        #   -s2 {percent2} \
        #   -st1 {signal_noise_threshold_M} \
        #   -st2 {signal_noise_threshold_L} \
        #   -o {output_file}".format(
        #                     signal_noise_file = signal_noise_file,
        #                     percent1 = config.param('signal_noise', 'percent1'),
        #                     percent2 = config.param('signal_noise', 'percent2'),
        #                     signal_noise_threshold_M = config.param('epiqc_report', 'signal_noise_threshold_M'),
        #                     signal_noise_threshold_L = config.param('epiqc_report', 'signal_noise_threshold_L'),
        #                     output_file = report_file)))

        input_dir = os.path.join(self.output_dirs['epigeec_output_directory'], self.output_dirs['epigeec_output'])
        jobs.append(Job(
            ['epiqc_report', input_dir],
            [],
            [],
            name="report_heatmap",
            command=
            "python ../genpipes/bfx/epiqc_report.py \
              -e {correlation_matrix} \
              -o {output_dir}".format(
                correlation_matrix=input_dir + "/correlation_matrix.tsv",
                output_dir=self.output_dirs['report_dir'])))

        return jobs

    @property
    def steps(self):
        # TODO : - Create steps table
        return [
           # self.test,
            self.bigwiginfo,
            self.bigwig_to_bedgraph,
            self.chromimpute_preprocess,
            self.chromimpute_convert,
            self.chromimpute_compute_global_dist,
            self.chromimpute_generate_train_data,
            self.chromimpute_train,
            self.chromimpute_apply,
            self.chromimpute_eval,
           # self.chromimpute_train_step,
           # self.chromimpute_compute_metrics,
             self.signal_to_noise,
            #self.epigeec_tohdf5,
         self.epigeec,
            # self.epiqc_report
        ]


if __name__ == "__main__":
    EpiQC()
