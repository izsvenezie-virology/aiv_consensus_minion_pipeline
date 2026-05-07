# AIV consensus MinION pipeline

This is a pipeline developed to create the consensus sequence of an Avian Influenza Virus (AIV) sample sequenced on MinION platform. It relies on a reference-based approach by automatically detecting a suitable reference, aligning data againts it with [minimap2](https://pubmed.ncbi.nlm.nih.gov/29750242/) and creating consensus sequence with a script developed at IZSVe.

## Installation

This pipeline is available for Linux and requires blast, chopper, perl, minimap2 and samtools as dependencies. The easiest way to install it is to clone current github repository with the following command
```
git clone https://github.com/izsvenezie-virology/aiv_consensus_minion_pipeline
```
and to create a mamba environment using the yml file provided.
```
mamba env create -p aiv_consensus_minion_pipeline_env -f env.yml
```
You can check if installation was successful in your system by activating the environment with the following command
```
mamba activate aiv_consensus_minion_pipeline_env
```
and then executing the pipeline on example data provided
```
bash minion_aiv_consensus.sh examples/sample.fastq.gz output_dir_example 1
```

## Usage

### Input
The pipeline tool requires a single FASTQ file as input, containing all the reads basecalled by MinION for a sample. In case reads belonging to a sample are divided into multiple files, you can concatenate them together using zcat utility, e.g. with the following command
```
zcat FOLDER_RUN/fastq_pass/barcodeXXX/*.gz | gzip >input_data.fastq.gz
```
where "XXX" in "barcodeXXX" is the barcode used for the wanted sample.

### Basic usage
With the mamba environment active, you can use the pipeline on a single sample to produce its consensus by running:
```
bash minion_aiv_consensus.sh input_data.fastq.gz OUTPUT.FOLDER 1
```
where
* "input_data.fastq.gz" is the fastq file containing sequenced data for your sample
* "OUTPUT.FOLDER" is the directory where computation will be done and consensus sequences will be saved in "consensus.fa" file
* "1" is the number of cpus used by the pipeline, the more the faster the analysis will be done. It is highly recommended to use the maximum number of CPUs available by your system to speed up the analysis.

## Cite the pipeline
Please cite the GitHub repository:

[https://github.com/izsvenezie-virology/aiv_consensus_minion_pipeline](https://github.com/izsvenezie-virology/aiv_consensus_minion_pipeline)

## License
GenIn is licensed under the GNU Affero v3 license (see [LICENSE](LICENSE)).

## Fundings

This work was supported by ??????.

>Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HEDEA). 
>Neither the European Union nor the granting authority can be held responsible for them
