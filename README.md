# covid-pipeline
Pipeline for analyzing SARS-CoV-2 sequence reads  

Current routine only handles Swift Biosciences data is handled.  

Codebase stage: development  
Developers and maintainers: Lavanya Rishishwar, Andrew Conley, Emily Norris, Anna Gaines  
Testers: Lavanya Rishishwar  


## Installation
We are working on improving the user experience with this as-of-unnamed software.  This section will be updated in near future with easier installation procedure for the core pipeline wand it's dependencies.

### Dependencies
This pipeline requires the following OS/programs/modules to run:
* Linux/Mac OS - tested on RedHat OS, expected to run on Mac environments, will possibly work on WSL though Mac and WSL haven't been tested so far
* fastqc
* multiqc
* sambamba
* bcftools
* tabix
* fastp or trimmomatic
* minimap2 or bwa
* primerclip
* pangolin 
* 
* Python version 3.7 or higher
* Python libraries = Pathos



#### Conda based installation
The following commands will install the dependencies and python modules, and clone the git repo.  I will be adding an environment file in near future.
```
# Create a conda environment
git clone https://github.com/cov-lineages/pangolin
cd pangolin
conda env update --name covidpipeline -f environment.yml

# Activate the newly created conda environment
conda activate covidpipeline

# Finish installing pangolin
python setup.py install
cd ..
rm -rf pangolin

# Install core packages
conda install -c conda-forge pathos nodejs
conda install -c bioconda fastqc sambamba bcftools tabix fastp trimmomatic minimap2 bwa primerclip

# these software are separated out otherwise conda dependency conflicts...
pip install multiqc
npm install --global @neherlab/nextclade

# Download this pipeline
git clone https://github.com/GA-Genomic-Epidemiology/covid-pipeline

# Give it executable permissions
cd covid-pipeline
chmod +x pipeline.py

# Run the pipeline
./pipeline.py

# Deactivate environment, when not using pastrami
conda deactivate covidpipeline
```

## Quickstart guide
This section will be populated in near future with small example datasets and commands for analyzing them.

## Basic usage
### General program structure overview
As of now, the pipeline only has a single subcommand implemented - all.  In future, additional subcommands will be added to provide users finer level control of the output.

```
usage: pipeline.py [--help] [--version] {all} ...

pipeline.py - SARS-CoV-2 sequence read analysis pipeline Version: 0.2

optional arguments:
  --help, -h, --h
  --version        show program's version number and exit

pipeline.py commands:
  {all}
    all            Perform complete analysis workflow
```

The all subcommand has a lot more options implemented within it:

```
usage: pipeline.py all [-h] [--directory <DIR>] [--metafile <FILE>] [--seq-prot <PROT>] [--out-prefix <OUTPREFIX>]
                       [--depth-filter <DP>] [--min-read-len <READLEN>] [--min-read-quality <READQUAL>]
                       [--min-base-quality <BASEQUAL>] [--min-genome-cov <GENCOV>] [--trimming-program <TRIMPROG>]
                       [--mapping-program <MAPPROG>] [--log-file run.log] [--threads N] [--verbose]

optional arguments:
  -h, --help                                                 show this help message and exit

Required Input options:
  --directory <DIR>, -d <DIR>, --d <DIR>                     Directory containing .fastq(.gz) files
  --metafile <FILE>, -m <FILE>, --m <FILE>                   Metafile containing sample to data mapping
  --seq-prot <PROT>, -s <PROT>, --s <PROT>                   Which sequencing protocols was used? Options: Swift.
                                                             (Default: Swift)

Output options:
  --out-prefix <OUTPREFIX>, -o <OUTPREFIX>, --o <OUTPREFIX>  Output directory where all the files will be created.
                                                             Output report will be created as <OUTPUTPREFIX>.report

Optional arguments for finer level parameter control:
  --depth-filter <DP>                                        Minimum read depth for variants (Default: 30)
  --min-read-len <READLEN>                                   Minimum sequencing read length (Default: 101)
  --min-read-quality <READQUAL>                              Minimum sequencing read quality for read filtering
                                                             (Default: 18)
  --min-base-quality <BASEQUAL>                              Minimum sequencing base quality for read trimming (Default:
                                                             18)
  --min-genome-cov <GENCOV>                                  Minimum % of genome covered (Default: 90)
  --trimming-program <TRIMPROG>                              Which trimming program to use? Options: fastp,trimmomatic
                                                             (Default: fastp)
  --mapping-program <MAPPROG>                                Which mapper to use? Options: minimap2,bwa (Default:
                                                             minimap2)

Runtime options:
  --log-file run.log, -l run.log, --l run.log                File containing log information (default: run.log)
  --threads N, -t N                                          Number of concurrent threads (default: 4)
  --verbose, -v                                              Print program progress information on screen

```

### A typical run
The ***all*** subcommand is the only available subcommand as of this date.  If you have your input FASTQ(.gz) files located in a directory called *input_fastq*, you can run the pipeline as:

```
./pipeline.py all -d input_fastq -v
# -v flag prints program activity on the screen so you can assess what step the program is at currently.
```
**NOTE**: By default, 4 threads are utilized.  Change this setting by the -t option to adapt to your computational infrastructure.
