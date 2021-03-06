# Drop-seq Tools Pipeline

A pipeline for processing raw single-cell RNA-seq generated by the Drop-seq protocol.

## Quickstart

To run this pipeline, you'll need [Nextflow](https://www.nextflow.io/) installed on your machine. Install it with the following command:
```
curl -s https://get.nextflow.io | bash
```
With Nextflow installed and accessible in your `PATH`, the pipeline can be run with the following command:
```
nextflow run gates-mri-bioinformatics/drop-seq-tools-pipeline -user <Github username>
```
You'll be prompted for your Github password. After providing it, Nextflow will automatically pull scripts from this repository
and run them in the environment specified by the configuration file `nextflow.config` (see more details on this file below). There's no need to clone this repository.

## Parameters

Pipeline parameters can be set in a `nextflow.config` file in the directory from which `nextflow` is launched:

```
params.output_directory = "/path/to/my/results"
params.fastq_directory = "/path/to/my/reads"
params.reference_fasta = "/path/to/my/sequence"
params.reference_gtf = "/path/to/my/annotations"
params.star_index = "/path/to/my/index
```

Alternatively, the required parameters can also be provided via the command line:
```
nextflow run gates-mri-bioinformatics/drop-seq-tools-pipeline -user liameabbott \
--output-directory /path/to/my/results \
--fastq-directory /path/to/my/reads \
--reference-fasta /path/to/my/sequence \
--reference-gtf /path/to/my/annotations
--star-index /path/to/my/index
```

Note the inclusion of the `params.` prefix in the config file version, and the use of underscores (`_`) in the config file version versus hyphens (`-`) in the command-line version.

#### Required:
```
--output-directory:
    Path to a directory to write pipeline output files.
--fastq-directory:
    Path to a directory with paired FASTQ files containing raw reads.
--reference-fasta:
    Path to a FASTA file containing a reference genome sequence.
--reference-gtf:
    Path to a GTF file containing annotations for the reference genome sequence.
--star-index:
    Path to a `STAR` index created from the `--reference-fasta` genome sequence.
```

#### Optional:
```
--mt-sequence (default "null"):
    The string used to specify mitochondrial contig in the reference files.
--cell-barcode-base-range (default "1-12"):
    The range of bases in the barcode read used to define the cell barcode.
--umi-barcode-base-range  (default "13-20"):
    The range of bases in the barcode read used to define the UMI barcode.
--barcoded-read (default "1"):
    The number of the barcode read in the read pair (either "1" or "2").
--barcode-min-base-quality (default "10"):
    The minimum acceptable base quality in the barcodes.
--barcode-n-bases-below-min-quality (default "1"):
    The number of bases in the barcode allowed to be below the quality threshold before the read is discarded.
--trim-starting-sequence (default "AAGCAGTGGTATCAACGCAGAGTGAATGGG"):
    The adapter sequence to look for at the beginning of the read.
--trim-starting-n-bases (default "5"):
    The number of bases at the beginning of the read that must match "--trim-starting-sequence" before trimming occurs.
--trim-starting-n-mismatches (default "0"):
    The number of mismatches allowed in the starting sequence before determining that the read should be trimmed.
--trim-polyA-adapter (default "~XM~XCGTACTCTGCGTTGATACCACTGCTT"):
    The adapter sequence. "^XM" references the value of the "XM" tag applied in the barcode tagging stage. "~XM" references the reverse complement of the tag.
--star-n-cores (default "1"):
    The number of CPU cores to use in `STAR` alignment.
--star-limit-out-sj-collapsed (default "5000000"):
    STAR parameter - maximum number of collapsed splice junctions.
--repair-min-umis-per-cell (default "20"):
    The minimum number of UMIs to consider a cell barcode for collapsing during barcode repair steps.
```


### Execution

With [Docker](https://www.docker.com/) installed on your machine, you can execute the pipeline in a containerized environment with all required dependencies installed:
```
nextflow run gates-mri-bioinformatics/drop-seq-tools-pipeline -user <Github username> -profile docker
```

