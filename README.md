# Challenge 1 Evaluation

The evaluation of submitted transcript models will be done using [SQANTI3](https://github.com/ConesaLab/SQANTI3) descriptors and number of LRGASP-agreed evaluation metrics. The general procedure to assess these evaluation metrics on your data is to run **sqanti3_lrgasp.py**, an adapted version of the original code that will generate automatically an HTML report with the results of the evaluation. Please, download/clone the entire [sqanti3_evaluation](https://github.com/LRGASP/lrgasp-challenge-1-evaluation.git) directory, including the **utilities** subfolder.

## Setting up the environment

In order to install all the dependencies needed by **sqanti3_lrgasp.py**, please use the [YML](https://github.com/LRGASP/lrgasp-challenge-1-evaluation/blob/main/sqanti3_lrgasp.yml) file to build a conda environment. 

```
conda env create -f sqanti3_lrgasp.yml
source activate sqanti3_lrgasp
```

SQANTI3 also takes advantage of some scripts of [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/wiki#install). Please install it with the *sqanti3_lrgasp* environment activated and add `cDNA_Cupcake`and `cDNA_Cupcake/sequence` to your `PYTHONPATH`.

```
(sqanti3_lrgasp)$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
(sqanti3_lrgasp)$ cd cDNA_Cupcake
(sqanti3_lrgasp)$ python setup.py build
(sqanti3_lrgasp)$ python setup.py install

(sqanti3_lrgasp)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/sequence/
(sqanti3_lrgasp)$ export PYTHONPATH=$PYTHONPATH:<path_to>/cDNA_Cupcake/

```

Remember to activate the *sqanti3_lrgasp* environment and setting up the `PYTHONPATH`every time you are running the evaluation.

## Run SQANTI3

When running [SQANTI3](https://github.com/ConesaLab/SQANTI3), your transcript-models will be compared against the GENCODE annotation using the last version available of the reference genome. We highly recommend to use the mouse and human genome annotation files available at the [LRGASP submission github](https://github.com/LRGASP/lrgasp-submissions/blob/master/docs/reference-genomes.md), as they already include the spike-ins information (SIRVs and ERCC).

LRGASP will be using CAGE peak data, polyA motif list and Illumina junction coverage to evaluate your transcript models using SQANTI3. We therefore recommend you run **sqanti3_lrgasp.py** enabling these analyses. To do so:

-   **CAGE peak data**:  Download BED files of CAGE peak data for [human](https://github.com/LRGASP/lrgasp-challenge-1-evaluation/blob/main/utilities/refTSS.human.bed) and [mouse](https://github.com/LRGASP/lrgasp-challenge-1-evaluation/blob/main/utilities/refTSS.mouse.bed) and provide them to **sqanti3_lrgasp.py** using the `--cage_peak` option
-   **polyA motif list**: This is a TXT file with the most common polyA motifs for human and mouse that can be downloaded from [here](https://github.com/LRGASP/lrgasp-challenge-1-evaluation/blob/main/utilities/polyA_list.txt). Include this file when running **sqanti3_qc.py** using the `--polyA_motif_list` option.
-   **SJ coverage**:  As SJ information is dependent on the sample being analyzed, it is necessary to run previously STAR to map the Illumina reads against the genome and identify possible splice-junctions using the `--twopassMode`. Then, the resulting _*SJ.oyut.tab_ file can be input to **sqanti3_lrgasp.py** with the parameter `-c`. This is an example of how we normally run STAR for this SJ detection:

1. Create a genome index with STAR without providing the reference annotation. We don't want to bias the SJ-detection towards the annotated splice sites.

```
STAR --runThreadN <num_threads> --runMode genomeGenerate --genomeDir <star_index> --genomeFastaFiles <reference_genome_FASTA> --outTmpDir <index_dir_tmp> 
```

2. Map short-reads using `--twopassMode`

```
STAR --runThreadN <num_threads> --genomeDir <star_index> --readFilesIn <read1> <read2> --outFileNamePrefix <output_prefix> --alignSJoverhangMin 8  --alignSJDBoverhangMin 1 --outFilterType BySJout --outSAMunmapped Within --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.04 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --sjdbScore 1 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --twopassMode Basic
```

It is also neccessary to provide a metadata file in JSON format. I should be like the experiment JSON file that is required to complete a submission. [Here](https://lrgasp.github.io/lrgasp-submissions/docs/metadata.html#experimentjson) you can find which information is expected to be provided through the metadata file.

### Example

This is an example of how to run the **sqanti3_lrgasp.py** script:

```
python sqanti3_lrgasp.py human_submitted.gtf lrgasp_gencode_v38.gtf lrgasp_grch38_sirvs.fasta \
	--gtf --json human_example.json  --cage_peak refTSS.human.bed \
	--polyA_motif_list polyA_list.txt -c my_test.SJ.out.tab \
	-d /my_output/directory -o human_submission_test
```

This will generate in `/my_output/directory` the two main SQANTI3 outputs (`*_classification.txt`and `*_junctions.txt`) and a HTML file that will be called, in this case, `human_submission_test_Evaluation_report.html`.


