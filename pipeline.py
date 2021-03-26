#!/usr/bin/env python3

__author__ = "Lavanya Rishishwar, Anna Gaines, Emily Norris, Andrew Conley"
__copyright__ = "Copyright 2021, Lavanya Rishishwar, Anna Gaines, Emily Norris, Andrew Conley"
__credits__ = ["Lavanya Rishishwar"]
__license__ = "GPL"
__version__ = "0.1"
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
import subprocess
import sys
from argparse import ArgumentParser, HelpFormatter

# import pathos.multiprocessing as mp

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
    def check_dependencies():
        pass

    @staticmethod
    def run_command(command_str: str = None, command_list: list = None, shell=False):
        if command_str is None and command_list is None:
            raise ValueError("Support.run_command() was called without any command to execute.")
        try:
            if command_str is not None:
                logging.info(f"Attempting to run: {command_str}")
                proc = subprocess.Popen(shlex.split(command_str), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        shell=shell, universal_newlines=True)
                output, stderr = proc.communicate()
                logging.debug(f"Error stream: {stderr}")
            else:
                logging.info(f"Attempting to run: " + " ".join([str(x) for x in command_list]))
                proc = subprocess.Popen(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        shell=shell, universal_newlines=True)
                output, stderr = proc.communicate()
                logging.debug(f"Error stream: {stderr}")
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
        logging.info("Command executed without raising any exceptions")
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
            os.makedirs(this_dir)
        except IOError:
            print(f"I don't seem to have access to output prefix directory. Are the permissions correct?")
            sys.exit(1)

    @staticmethod
    def safe_dir_rm(this_dir: str):
        try:
            os.rmdir(this_dir)
        except IOError:
            print(f"I don't seem to have access to output prefix directory. Are the permissions correct?")
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
    ftypes = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]

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
        self.threads = opts.threads
        self.log_file = opts.log_file
        self.verbosity = opts.verbosity

        # Samples collected
        self.samples = {}
        # Any errors we encounter
        self.errors = []
        # The actual queue for the analysis
        self.analysis = []

        # Verbosity levels and colors
        self.error_color = Colors.FAIL
        self.warning_color = Colors.WARNING
        self.main_process_color = Colors.OKGREEN
        self.sub_process_color = Colors.OKBLUE

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
        """
        return long_string

    def validate_options(self):
        if self.sub_command == "all":
            self.validate_all_arguments()
            self.analysis += ['get_files']

    # TODO: Implement validations for file presence and stuff
    def validate_all_arguments(self):
        if self.directory is None:
            self.errors += [self.sub_command + ' requires --directory']
            return

    def summarize_run(self):
        logging.info(self.main_process_color + str(self) + Colors.ENDC)
        logging.info(self.main_process_color + f"Analysis to perform: " + ",".join(self.analysis) + Colors.ENDC)

    def go(self):
        self.summarize_run()
        while True:
            step = self.analysis[0]
            self.analysis = self.analysis[1:]
            function = getattr(self, step)
            function()
            if len(self.analysis) == 0:
                break

    # TODO: Make fasta index
    def preprocess_ref_files(self):
        pass

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

        if len(self.errors) > 0:
            Support.error_out(messages=self.errors)

        logging.info(f"Read {len(self.samples)} samples.")
        logging.info(self.samples)

    # TODO: polish code
    def perform_qa(self):
        # for sample in self.samples:
        #     Support.run_command(f"fastqc {} {}")
        pass

    # TODO: polish code
    def perform_qc(self):
        pass

    # TODO: polish code
    def call_variants(self):
        pass

    # TODO: polish code
    # TODO: Implement push to a big database
    def filter_variants(self):
        pass

    # TODO: polish code
    def call_consensus(self):
        pass

    # TODO: polish code
    def run_pangolin(self):
        pass

    # TODO: polish code
    def run_nextclade(self):
        pass


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
    all_output_group.add_argument('--out-prefix', required=False, default=f"output-{TIMESTAMP}", metavar='<OUTPREFIX>',
                                  type=str, dest="out_prefix",
                                  help="""Output directory where all the files will be created. 
                                  Output report will be created as <OUTPUTPREFIX>.report""")

    # Add all finer parameters arguments
    all_filter_group = all_parser.add_argument_group('Optional arguments for finer level parameter control')
    all_filter_group.add_argument('--depth-filter', required=False, default=30, metavar='<DP>', type=int,
                                  help='Minimum read depth for variants (Default: %(default)s)',
                                  dest="depth_filter")
    all_filter_group.add_argument('--min-read-len', required=False, default=120, metavar='<READLEN>', type=int,
                                  help='Minimum sequencing read length (Default: %(default)s)',
                                  dest="min_read_len")
    all_filter_group.add_argument('--min-read-quality', required=False, default=18, metavar='<READQUAL>', type=int,
                                  help='Minimum sequencing read quality for read filtering (Default: %(default)s)',
                                  dest="min_read_qual")
    all_filter_group.add_argument('--min-base-quality', required=False, default=18, metavar='<BASEQUAL>', type=int,
                                  help='Minimum sequencing base quality for read trimming (Default: %(default)s)',
                                  dest="min_base_qual")

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
