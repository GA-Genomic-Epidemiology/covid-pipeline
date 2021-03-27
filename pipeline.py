#!/usr/bin/env python3

__author__ = "Lavanya Rishishwar, Anna Gaines, Emily Norris, Andrew Conley"
__copyright__ = "Copyright 2021, Lavanya Rishishwar, Anna Gaines, Emily Norris, Andrew Conley"
__credits__ = ["Lavanya Rishishwar", "Emily Norris", "Anna Gaines", "Andrew Conley"]
__license__ = "GPL"
__version__ = "0.2"
__maintainer__ = "Lavanya Rishishwar, Anna Gaines, Emily Norris, Andrew Conley"
__email__ = "lavanyarishishwar@gmail.com"
__status__ = "Development"

import datetime
import glob
import logging
import os.path
import random
import re
import shlex
import string
import shutil
import subprocess
import sys
import time
from argparse import ArgumentParser, HelpFormatter

import pathos.multiprocessing as mp

PROGRAM_NAME = "pipeline.py"
VERSION = __version__
PROGRAM_DESCRIPTION = "SARS-CoV-2 sequence read analysis pipeline"
TIMESTAMP = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


class Support:
    @staticmethod
    def error_out(message: str = None, messages: list = None):
        if message is not None:
            sys.exit(Colors.FAIL + f"Error: {message}" + Colors.ENDC)
        elif messages is not None:
            sys.exit(Colors.FAIL + f"Encountered the following errors, can not continue:\n" + "\n".join(
                messages) + Colors.ENDC)
        else:
            sys.exit("The program encountered an error and has to exit.")

    @staticmethod
    def validate_file(the_file: str):
        return os.path.isfile(the_file)

    @staticmethod
    def validate_file_size(the_file: str):
        return os.stat(the_file).st_size > 0

    @staticmethod
    def validate_file_and_size_or_error(the_file: str, error_prefix: str = 'The file',
                                        presence_suffix: str = 'doesn\'t exist',
                                        size_suffix: str = 'is size 0'):
        if not Support.validate_file(the_file=the_file):
            print(' '.join([error_prefix, the_file, presence_suffix]), 0, Colors.FAIL)
            Support.error_out()

        if not Support.validate_file_size(the_file=the_file):
            print(' '.join([error_prefix, the_file, size_suffix]), 0, Colors.FAIL)
            Support.error_out()

    @staticmethod
    def validate_dir(the_dir: str):
        return os.path.isdir(the_dir)

    # TODO: This is a hastily written method, needs error fixing
    @staticmethod
    def validate_dir_or_error(the_dir: str, error_prefix: str = "The dir", presence_suffix: str = "doesn't exist"):
        if not Support.validate_dir(the_dir=the_dir):
            print(' '.join([error_prefix, the_dir, presence_suffix]), 0, Colors.FAIL)
            Support.error_out()

    # TODO: Implement checks for dependency progams
    @staticmethod
    def check_dependencies(program_list: list = None):
        pass

    @staticmethod
    def run_command(command_str: str = None, command_list: list = None, shell=False, split=True):
        if command_str is None and command_list is None:
            raise ValueError("Support.run_command() was called without any command to execute.")
        try:
            if command_str is not None:
                logging.debug(f"Attempting to run: {command_str}")
                if split:
                    command_str = shlex.split(command_str)
                proc = subprocess.Popen(command_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        shell=shell, universal_newlines=True)
                output, stderr = proc.communicate()
                logging.debug(f"Standard Error stream:\n{stderr}")
            else:
                logging.debug(f"Attempting to run: " + " ".join([str(x) for x in command_list]))
                proc = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        shell=shell, universal_newlines=True)
                output, stderr = proc.communicate()
                logging.debug(f"Standard Error stream:\n{stderr}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Encountered an error executing the command: ")
            if command_str is not None:
                logging.error(command_str)
            else:
                logging.error(command_list)
            logging.error(f"Error details:")
            logging.error(f"Exit code={e.returncode}")
            logging.error(f"Error message={e.output}")
            sys.exit(1)
        logging.debug("Command executed without raising any exceptions")
        return output

    @staticmethod
    def validate_filename(filename: str):
        if re.match(r"^[a-zA-Z0-9_.-]+$", filename):
            return True
        else:
            return False

    @staticmethod
    def validate_output_prefix(out_prefix: str):
        parent, prefix = os.path.split(out_prefix)
        if parent != "":
            if not Support.validate_dir(parent):
                Support.safe_dir_create(parent)
        return Support.validate_filename(prefix)

    @staticmethod
    def safe_dir_create(this_dir: str):
        try:
            # TODO: Migrate to newer Pathlib based dir creation
            os.makedirs(this_dir, exist_ok=True)
        except IOError:
            print(f"I don't seem to have access to create directory '{this_dir}'. Are the permissions correct?")
            sys.exit(1)

    @staticmethod
    def safe_dir_rm(this_dir: str):
        if not os.path.isdir(this_dir):
            logging.debug(f"{this_dir} doesn't exist")
            return
        try:
            shutil.rmtree(this_dir, ignore_errors=True)
        except IOError:
            print(f"I don't seem to have access to remove the directory '{this_dir}'. Are the permissions correct?")
            sys.exit(1)

    @staticmethod
    def init_logger(log_file, verbosity):
        """Configures the logging for printing
        Returns
        -------
        None
            Logger behavior is set based on the Inputs variable
        """
        try:
            logging.basicConfig(filename=log_file, filemode="w", level=logging.DEBUG,
                                format=f"[%(asctime)s] %(message)s",
                                datefmt="%m-%d-%Y %I:%M:%S %p")
        except FileNotFoundError:
            print(f"The supplied location for the log file '{log_file}'" +
                  f"doesn't exist. Please check if the location exists.")
            sys.exit(1)
        except IOError:
            print(f"I don't seem to have access to make the log file." +
                  f"Are the permissions correct or is there a directory with the same name?")
            sys.exit(1)

        if verbosity:
            console = logging.StreamHandler()
            console.setLevel(logging.INFO)
            formatter = logging.Formatter(fmt=f"[%(asctime)s] %(message)s", datefmt="%m-%d-%Y %I:%M:%S %p")
            console.setFormatter(formatter)
            logging.getLogger().addHandler(console)


class Analysis:
    # TODO: package these files
    reference_file = "NC_045512.2.fasta"
    reference_mmi = reference_file.replace(".fasta", ".mmi")
    reference_seqid = "NC_045512.2"
    adapter_file = "adapters.fasta"
    swift_regions = "25-29852"
    ftypes = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    program_list = ["fastqc", "fastp", "trimmomatic", "bwa2", "bwa", "sambamba", "bcftools", "tabix"]
    out_subdir = {'qa': '1-QualityAssessment', 'qc': '2-QualityControl', 'map': '3-Mapping',
                  'cov': '4-CoverageCalculation', 'varcall': '4-VariantCalling', 'varfilt': '5-VariantFiltering',
                  'cons': '6-ConsensusCalling', 'pangolin': '7-Pangolin', 'nextclade': '8-NextClade'}

    def __init__(self, opts):
        # General attributes
        self.sub_command = opts.sub_command
        self.directory = opts.directory
        self.metafile = opts.metafile
        self.seqprot = opts.seqprot
        self.out_prefix = opts.out_prefix
        self.depth_filter = opts.depth_filter
        self.min_read_len = opts.min_read_len
        self.min_read_qual = opts.min_read_qual
        self.min_base_qual = opts.min_base_qual
        self.min_genome_cov = opts.min_genome_cov
        self.trimming_program = opts.trimming_program
        self.mapping_program = opts.mapping_program

        self.threads = opts.threads
        self.log_file = opts.log_file
        self.verbosity = opts.verbosity

        # Any errors we encounter
        self.errors = []
        # The actual queue for the analysis
        self.analysis = []
        # Samples collected
        self.samples = {}
        # Output report
        self.report = {}

        # Verbosity levels and colors
        self.error_color = Colors.FAIL
        self.warning_color = Colors.WARNING
        self.main_process_color = Colors.OKGREEN
        self.sub_process_color = Colors.OKBLUE

        self.time = time.time()
        # Initialize the logger
        Support.init_logger(log_file=self.log_file, verbosity=self.verbosity)

    def __str__(self):
        long_string = f"""
        Class constants: 
                None

        Instance variables:
                * General program parameters
                log_file                      =  {self.log_file}
                threads                       =  {self.threads}
                verbosity                     =  {self.verbosity}
                * Subcommand to be executed
                sub_command                   =  {self.sub_command}
                * Input options
                directory                     =  {self.directory}
                metafile                      =  {self.metafile}
                seqprot                       =  {self.seqprot}
                * Output options
                out_prefix                    =  {self.out_prefix}
                * Fine tuning options
                depth_filter                  =  {self.depth_filter}
                min_read_len                  =  {self.min_read_len}
                min_read_qual                 =  {self.min_read_qual}
                min_base_qual                 =  {self.min_base_qual}
                min_genome_cov                =  {self.min_genome_cov}
                trimming_program              =  {self.trimming_program}
                mapping_program               =  {self.mapping_program}
        """
        return long_string

    def validate_options(self):
        if self.sub_command == "all":
            self.validate_all_arguments()
            self.analysis += ['get_files', 'perform_qa', 'perform_qc', 'map_reads',
                              'calculate_coverage_metrics', 'call_variants', 'filter_variants', 'call_consensus',
                              'run_pangolin', 'run_nextclade', 'print_report']

    # TODO: Implement validations for file presence and stuff
    def validate_all_arguments(self):
        if self.directory is None:
            self.errors += [self.sub_command + ' requires --directory']
            return

        if self.mapping_program == "bwa":
            bwa_ext = ["amb", "ann", "bwt", "fai", "pac", "sa"]
            for ext in bwa_ext:
                if not os.path.isfile(f"{Analysis.reference_file}.{ext}"):
                    self.analysis += ['preprocess_ref_files']
                    break
        else:
            if not os.path.isfile(Analysis.reference_mmi) or not os.path.isfile(f"{Analysis.reference_file}.fai"):
                self.analysis += ['preprocess_ref_files']

        if os.path.isdir(self.out_prefix):
            logging.warning(self.warning_color +
                            f"Output directory '{self.out_prefix}' exists, will be overwritten" + Colors.ENDC)
            Support.safe_dir_rm(self.out_prefix)

        Support.safe_dir_create(self.out_prefix)

    def summarize_run(self):
        logging.info(self.main_process_color + str(self) + Colors.ENDC)
        logging.info(self.main_process_color + f"Analysis to perform: " + " -> ".join(
            [x.replace("_", " ") for x in self.analysis]) + Colors.ENDC)

    def go(self):
        self.summarize_run()
        time_dict = {}
        total_time = 0
        while True:
            step = self.analysis[0]
            self.analysis = self.analysis[1:]
            function = getattr(self, step)

            start_time = time.time()
            function()
            time_dict[step] = time.time() - start_time
            total_time += time_dict[step]
            if len(self.analysis) == 0:
                logging.info(f"Full process took {round(time.time() - self.time, 2)} seconds")
                break
        with open(os.path.join(self.out_prefix, "runtime.log"), "w") as pf:
            pf.write("Step\tTime(sec)\tPercentTotal\n")
            for step in time_dict:
                pf.write(f"{step}\t{round(time_dict[step], 2)}\t{round(100 * time_dict[step] / total_time, 2)}\n")

    # TODO: Make fasta index
    def preprocess_ref_files(self):
        if self.mapping_program == "bwa":
            Support.run_command(f"bwa index {Analysis.reference_file}")
        else:
            Support.run_command(f"minimap2 -d {Analysis.reference_mmi} {Analysis.reference_file}")
        Support.run_command(f"samtools index {Analysis.reference_file}")

    # TODO: Implement safe file pulls
    # TODO: Implement lane merging
    def get_files(self):
        files = []
        for ext in Analysis.ftypes:
            files = files + glob.glob(os.path.join(self.directory, f"*{ext}"))
        for file in files:
            basename = os.path.basename(file)
            sample_id, snum, lane, readn, suffix = basename.split(".")[0].split("_")

            # If sample is undertermined, skip
            if sample_id == "Undetermined":
                continue

            # Sample id will be sample_id_snum
            sample_id = f"{sample_id}_{snum}"
            if sample_id not in self.samples:
                self.samples[sample_id] = {"r1": "", "r2": ""}

            if readn == "R1":
                if self.samples[sample_id]["r1"] != "":
                    logging.error(f"""Found another R1 for {sample_id}, possibly from a different lane.
                                    This feature isn't implemented yet.""")
                else:
                    self.samples[sample_id]["r1"] = file
            elif readn == "R2":
                if self.samples[sample_id]["r2"] != "":
                    logging.error(f"""Found another R2 for {sample_id}, possibly from a different lane.
                                    This feature isn't implemented yet.""")
                else:
                    self.samples[sample_id]["r2"] = file
            else:
                logging.warning(f"Discarding read: {file}")

        for sample_id in self.samples:
            if self.samples[sample_id]["r1"] == "" or self.samples[sample_id]["r2"] == "":
                logging.debug(f"{sample_id} only has one file from the pair")
                self.errors.append(f"{sample_id} only has one file from the pair")
            else:
                self.report[sample_id] = {}

        if len(self.errors) > 0:
            Support.error_out(messages=self.errors)

        logging.info(f"Read {len(self.samples)} samples.")
        logging.info("Samples found: " + ", ".join(list(self.samples.keys())))

    # TODO: polish code
    def perform_qa(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['qa'])
        logging.info(f"Starting quality assessment (FastQC), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        cmd_queue = []
        for sample_id in self.samples:
            logging.info(f"Queuing {sample_id}")
            r1 = self.samples[sample_id]["r1"]
            r2 = self.samples[sample_id]["r2"]
            cmd = f"fastqc -o {outdir} {r1} {r2}"
            # cmd = f"fastqc -o {outdir} -t {self.threads} {r1} {r2}"
            # Support.run_command(command_str=cmd)
            cmd_queue.append(cmd)
        pool = mp.Pool(processes=self.threads)
        pool.map(lambda x: Support.run_command(command_str=x), cmd_queue)
        logging.info(f"Quality assessed for all samples")

        # TODO: Implement multiqc
        # logging.info(f"Compiling QA reports")
        # cmd = f"multiqc {outdir}"
        # Support.run_command(cmd)

    # TODO: polish code
    def perform_qc(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['qc'])
        logging.info(f"Starting quality control ({self.trimming_program}), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            logging.info(f"Trimming {sample_id}")
            self.samples[sample_id]["qc_r1"] = f"{outdir}/{sample_id}_R1.trimmed.fastq.gz"
            self.samples[sample_id]["qc_r2"] = f"{outdir}/{sample_id}_R2.trimmed.fastq.gz"
            r1 = self.samples[sample_id]["r1"]
            r2 = self.samples[sample_id]["r2"]
            trim_r1 = self.samples[sample_id]["qc_r1"]
            trim_r2 = self.samples[sample_id]["qc_r2"]

            if self.trimming_program == "trimmomatic":
                cmd = f"trimmomatic PE -threads {self.threads} -trimlog {outdir}/{sample_id}.trimlog "
                cmd += f"{r1} {r2} {trim_r1} /dev/null {trim_r2} /dev/null "
                cmd += f"ILLUMINACLIP:{Analysis.adapter_file}:2:30:10:1:true"
                cmd += f"SLIDINGWINDOW:10:{self.min_base_qual} MINLEN:{self.min_read_len} AVGQUAL:{self.min_read_qual}"
            else:
                cmd = f"fastp -i {r1} -o {trim_r1} -I {r2} -O {trim_r2} --adapter_fasta {Analysis.adapter_file} -5 -W 10"
                cmd += f" -M {self.min_read_qual} -e {self.min_read_qual} -l {self.min_read_len} -w {self.threads}"
            Support.run_command(command_str=cmd)

    # TODO: polish code
    def map_reads(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['varcall'])
        logging.info(f"Starting read-to-genome mapping ({self.mapping_program}), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            logging.info(f"Mapping {sample_id}")
            self.samples[sample_id]["bam"] = f"{outdir}/{sample_id}.bam"
            r1 = self.samples[sample_id]["qc_r1"]
            r2 = self.samples[sample_id]["qc_r2"]
            sam = f"{outdir}/{sample_id}.sam"
            ubam = f"{outdir}/{sample_id}.unsorted.bam"
            bam = self.samples[sample_id]["bam"]

            if self.mapping_program == "bwa":
                cmd = f"bwa mem {Analysis.reference_file} {r1} {r2} -M -t {self.threads} -o {sam}"
            else:
                cmd = f"minimap2 -ax sr -o {sam} -t {self.threads} {Analysis.reference_mmi} {r1} {r2}"
            Support.run_command(command_str=cmd)

            # TODO: Calculate mapped reads %

            logging.info(f"Converting {sample_id} SAM to BAM")
            Support.run_command(command_str=f"""sambamba view -F 'not (unmapped or mate_is_unmapped)'
             -f bam -S {sam} -o {ubam} -t {self.threads}""")
            os.remove(sam)

            logging.info(f"Sorting {sample_id} BAM")
            Support.run_command(command_str=f"sambamba sort -o {bam} -t {self.threads} {ubam}")
            os.remove(ubam)

    # TODO: polish code
    def calculate_coverage_metrics(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['cov'])
        logging.info(f"Starting coverage calculation (HTSlib), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            logging.info(f"Calculating depth for {sample_id}")
            bam = self.samples[sample_id]["bam"]
            depth_file = os.path.join(outdir, f"{sample_id}.depth")

            cmd = f"samtools depth -d 0 -a -r '{Analysis.reference_seqid}:{Analysis.swift_regions}' {bam} > {depth_file}"

            Support.run_command(command_str=cmd, shell=True, split=False)

            total_bases = 0
            depth_coverage = 0
            genome_coverage = 0

            with open(depth_file, "r") as f:
                for line in f:
                    depth = int(line.rstrip().split("\t")[2])
                    if depth >= self.depth_filter:
                        genome_coverage += 1

                    depth_coverage += depth
                    total_bases += 1

            self.report[sample_id]["Depth Coverage"] = round(depth_coverage / total_bases, 2)
            self.report[sample_id]["Genome Coverage"] = round(genome_coverage * 100 / total_bases, 2)

            if self.report[sample_id]["Depth Coverage"] < self.depth_filter or \
                    self.report[sample_id]["Depth Coverage"] < self.min_genome_cov:
                self.report[sample_id]["Coverage Pass"] = "Fail"
            else:
                self.report[sample_id]["Coverage Pass"] = "Pass"

    # TODO: polish code
    def call_variants(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['varcall'])
        logging.info(f"Starting variant calling (BWA/HTSlib), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            if self.report[sample_id]["Coverage Pass"] == "Fail":
                continue
            self.samples[sample_id]["vcf"] = f"{outdir}/{sample_id}.vcf.gz"
            bam = self.samples[sample_id]["bam"]
            vcf = self.samples[sample_id]["vcf"]

            logging.info(f"Calling variants {sample_id}")
            mpileup_cmd = f"bcftools mpileup -Ou -Q 0 -B -a FORMAT/AD,FORMAT/DP -L 10000 -f {Analysis.reference_file} {bam}"
            call_cmd = f"bcftools call -mv -Oz -o {vcf} --ploidy 1"
            Support.run_command(command_str=f"{mpileup_cmd} | {call_cmd}", shell=True, split=False)

            logging.info(f"Indexing {sample_id}")
            Support.run_command(command_str=f"tabix -p vcf {vcf}")

    # TODO: polish code
    # TODO: Implement push to a big database
    def filter_variants(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['varfilt'])
        logging.info(f"Filtering variants (HTSlib), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            if self.report[sample_id]["Coverage Pass"] == "Fail":
                continue
            logging.info(f"Filtering variants for {sample_id}")
            self.samples[sample_id]["filtvcf"] = f"{outdir}/filt.{sample_id}.vcf.gz"

            in_vcf = self.samples[sample_id]["vcf"]
            out_vcf = self.samples[sample_id]["filtvcf"]

            cmd = f"bcftools filter -i'FORMAT/DP>{self.depth_filter} && QUAL>30' "
            cmd += f"-r'{Analysis.reference_seqid}:{Analysis.swift_regions}' {in_vcf} -O z -o {out_vcf}"
            Support.run_command(command_str=cmd)

            cmd = f"bcftools index {out_vcf}"
            Support.run_command(command_str=cmd)

    # TODO: polish code
    def call_consensus(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['cons'])
        logging.info(f"Consensus calling (HTSlib), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            if self.report[sample_id]["Coverage Pass"] == "Fail":
                continue
            logging.info(f"Calling consensus sequence for {sample_id}")
            self.samples[sample_id]["cons"] = f"{outdir}/{sample_id}.fasta"

            if "filtvcf" in self.samples[sample_id]:
                in_vcf = self.samples[sample_id]["filtvcf"]
            else:
                in_vcf = self.samples[sample_id]["vcf"]
            out_fa = self.samples[sample_id]["cons"]
            cmd = f"bcftools consensus -f {Analysis.reference_file} -o {out_fa} {in_vcf}"
            Support.run_command(command_str=cmd)

    # TODO: polish code
    def run_pangolin(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['pangolin'])
        outfile = "pangolin_consensus.csv"
        logging.info(f"Identifying lineages (Pangolin), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            if self.report[sample_id]["Coverage Pass"] == "Fail":
                self.report[sample_id]["Lineage"] = "."
                self.report[sample_id]["Lineage_Probability"] = "."
                self.report[sample_id]["Pangolin_QC"] = "."
                self.report[sample_id]["Pangolin_Version"] = "."
                continue

            logging.info(f"Identifying lineage for {sample_id}")
            self.samples[sample_id]["pangolin-dir"] = f"{outdir}/pangolin-{sample_id}"

            in_fa = self.samples[sample_id]["cons"]
            save_dir = self.samples[sample_id]["pangolin-dir"]

            cmd = f"pangolin {in_fa} -o {save_dir} --outfile {outfile}"
            Support.run_command(command_str=cmd)

            with open(f"{save_dir}/{outfile}", "r") as f:
                info = f.readlines()[1].split(",")
                self.report[sample_id]["Lineage"] = info[1]
                self.report[sample_id]["Lineage_Probability"] = info[2]
                self.report[sample_id]["Pangolin_QC"] = info[4]
                self.report[sample_id]["Pangolin_Version"] = info[3]
                logging.info(f"{sample_id} has lineage " + self.report[sample_id]["Lineage"])

    # TODO: polish code
    def run_nextclade(self):
        outdir = os.path.join(self.out_prefix, Analysis.out_subdir['nextclade'])
        logging.info(f"Identifying clades (Nextclade), results will be stored here: {outdir}")
        Support.safe_dir_create(outdir)
        for sample_id in self.samples:
            if self.report[sample_id]["Coverage Pass"] == "Fail":
                self.report[sample_id]["Clade"] = "."
                continue

            logging.info(f"Identifying clade for {sample_id}")
            self.samples[sample_id]["nextclade"] = f"{outdir}/nextclade-{sample_id}.tsv"

            in_fa = self.samples[sample_id]["cons"]
            out_tsv = self.samples[sample_id]["nextclade"]

            cmd = f"nextclade --input-fasta {in_fa} --output-tsv {out_tsv}"
            Support.run_command(command_str=cmd)

            with open(out_tsv, "r") as f:
                clade = f.readlines()[1].split("\t")[1]
                self.report[sample_id]["Clade"] = clade
                logging.info(f"{sample_id} has clade {clade}")

    # TODO: Pretty print
    def print_report(self):
        outfile = os.path.join(self.out_prefix, "report.csv")
        samples = list(self.samples.keys())
        headers = list(self.report[samples[0]].keys())
        with open(outfile, "w") as f:
            f.write("Sample," + ",".join(headers) + "\n")
            for sample_id in samples:
                f.write(f"{sample_id}," + ",".join([str(self.report[sample_id][x]) for x in headers]) + "\n")

        logging.info(f"Final report file available at: {outfile}")


if __name__ == '__main__':
    parser = ArgumentParser(prog=PROGRAM_NAME,
                            add_help=False,
                            description=f'''
                            {PROGRAM_NAME} - {PROGRAM_DESCRIPTION}
                            Version: {VERSION}
                            ''',
                            formatter_class=lambda prog: HelpFormatter(prog, width=120, max_help_position=120))

    parser.add_argument('--help', '-h', '--h', action='store_true', default=False)
    parser.add_argument('--version', action='version', version=f"{PROGRAM_NAME} {VERSION}")

    subparsers = parser.add_subparsers(title=f"{PROGRAM_NAME} commands")

    # Make the "all" sub-command parser
    all_parser = subparsers.add_parser('all', help='Perform complete analysis workflow',
                                       formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                  max_help_position=120))
    all_parser.set_defaults(sub_command='all')

    # Add all input arguments
    all_input_group = all_parser.add_argument_group('Required Input options')
    all_input_group.add_argument('--directory', '-d', '--d', required=False, default=None, metavar='<DIR>',
                                 help='Directory containing .fastq(.gz) files', dest="directory")
    all_input_group.add_argument('--metafile', '-m', '--m', required=False, default=None, metavar='<FILE>',
                                 help='Metafile containing sample to data mapping', dest="metafile")
    all_input_group.add_argument('--seq-prot', '-s', '--s', required=False, default="Swift",
                                 metavar='<PROT>', dest="seqprot",
                                 help='Which sequencing protocols was used? Options: Swift. (Default: %(default)s)')

    # Add all output arguments
    all_output_group = all_parser.add_argument_group('Output options')
    all_output_group.add_argument('--out-prefix', '-o', '--o', required=False, default=f"output-{TIMESTAMP}",
                                  metavar='<OUTPREFIX>', type=str, dest="out_prefix",
                                  help="""Output directory where all the files will be created. 
                                  Output report will be created as <OUTPUTPREFIX>.report""")

    # Add all finer parameters arguments
    all_filter_group = all_parser.add_argument_group('Optional arguments for finer level parameter control')
    all_filter_group.add_argument('--depth-filter', required=False, default=30, metavar='<DP>', type=int,
                                  help='Minimum read depth for variants (Default: %(default)s)',
                                  dest="depth_filter")
    all_filter_group.add_argument('--min-read-len', required=False, default=101, metavar='<READLEN>', type=int,
                                  help='Minimum sequencing read length (Default: %(default)s)',
                                  dest="min_read_len")
    all_filter_group.add_argument('--min-read-quality', required=False, default=18, metavar='<READQUAL>', type=int,
                                  help='Minimum sequencing read quality for read filtering (Default: %(default)s)',
                                  dest="min_read_qual")
    all_filter_group.add_argument('--min-base-quality', required=False, default=18, metavar='<BASEQUAL>', type=int,
                                  help='Minimum sequencing base quality for read trimming (Default: %(default)s)',
                                  dest="min_base_qual")
    all_filter_group.add_argument('--min-genome-cov', required=False, default=90, metavar='<GENCOV>', type=float,
                                  help='Minimum %% of genome covered (Default: %(default)s)',
                                  dest="min_genome_cov")
    all_filter_group.add_argument('--trimming-program', required=False, default="fastp", metavar='<TRIMPROG>', type=str,
                                  help='Which trimming program to use? Options: fastp,trimmomatic (Default: %(default)s)',
                                  dest="trimming_program")
    all_filter_group.add_argument('--mapping-program', required=False, default="minimap2", metavar='<MAPPROG>',
                                  type=str,
                                  help='Which mapper to use? Options: minimap2,bwa (Default: %(default)s)',
                                  dest="mapping_program")

    # Add all performance arguments
    all_runtime_group = all_parser.add_argument_group('Runtime options')
    all_runtime_group.add_argument('--log-file', "-l", "--l", required=False, default="run.log",
                                   metavar='run.log', type=str, dest="log_file",
                                   help='File containing log information (default: %(default)s)')
    all_runtime_group.add_argument('--threads', "-t", required=False, default=4, metavar='N', type=int,
                                   help='Number of concurrent threads (default: %(default)s)', dest="threads")
    all_runtime_group.add_argument('--verbose', "-v", required=False, action="store_true",
                                   help='Print program progress information on screen', dest="verbosity")

    options, unknown_arguments = parser.parse_known_args()

    if len(unknown_arguments) > 0:
        print(Colors.HEADER)
        print("User specificed unknown arguments: " + " ".join([str(x) for x in unknown_arguments]))
        print("Please see the correct usage below:")
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.help or "sub_command" not in options:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if options.sub_command == 'all' and options.help:
        print(Colors.HEADER)
        all_parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    # Start the analysis
    analysis = Analysis(options)
    analysis.validate_options()

    if len(analysis.errors) > 0:
        print(Colors.HEADER)
        if options.sub_command == 'all':
            all_parser.print_help()
        print(Colors.ENDC)
        print(Colors.FAIL + '\n\nErrors:')
        print("\n".join(analysis.errors))
        print(Colors.ENDC)
        sys.exit()

        # If we're still good, start the actual analysis
    analysis.go()
