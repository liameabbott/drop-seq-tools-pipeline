
nextflow.enable.dsl=2

log.info """\

====================================
|| drop-seq-tools-pipeline - v0.1 ||
====================================
| execution                        
------------------------------------
| script            : ${workflow.scriptFile}
| repository        : ${workflow.repository}
| commit            : ${workflow.commitId}  
| project_directory : ${workflow.projectDir}
| launch_directory  : ${workflow.launchDir}
| work_directory    : ${workflow.workDir}
| profile           : ${workflow.profile} 
| docker_image      : ${workflow.container}
====================================
| required parameters              
------------------------------------
| output_directory : ${params.output_directory}
| fastq_directory  : ${params.fastq_directory}
| reference_fasta  : ${params.reference_fasta}
| reference_gtf    : ${params.reference_gtf}
| star_index       : ${params.star_index}
====================================
| optional parameters              
------------------------------------
| mt_sequence                       : ${params.mt_sequence}
| cell_barcode_base_range           : ${params.cell_barcode_base_range}
| umi_barcode_base_range            : ${params.umi_barcode_base_range}
| barcoded_read                     : ${params.barcoded_read}
| barcode_min_base_quality          : ${params.barcode_min_base_quality}
| barcode_n_bases_below_min_quality : ${params.barcode_n_bases_below_min_quality}
| trim_starting_sequence            : ${params.trim_starting_sequence}
| trim_starting_n_bases             : ${params.trim_starting_n_bases}
| trim_starting_n_mismatches        : ${params.trim_starting_n_mismatches}
| trim_polyA_adapter                : ${params.trim_polyA_adapter}
| star_n_cores                      : ${params.star_n_cores}
| star_limit_out_sj_collapsed       : ${params.star_limit_out_sj_collapsed}
| repair_min_umis_per_cell          : ${params.repair_min_umis_per_cell}
------------------------------------

"""

// import metadata preparation modules
include {
	normalize_fasta;
	create_sequence_dictionary;
	convert_to_refFlat;
	reduce_gtf;
	create_intervals;
	convert_intervals_to_bed;
	convert_intervals_to_bed_mt;
} from './prepare-metadata.nf'

// import read processing modules
include {
	convert_reads;
	prepare_reads;
	align_annotate_reads;
	repair_reads;
} from './process-reads.nf'

process get_picard_jar {
	output:
	path("picard.jar")

	"""
	wget -O picard.jar \
		https://github.com/broadinstitute/picard/releases/download/${params.picard_version}/picard.jar
	"""
}

// run the whole pipeline
workflow {

	// PREPARE METADATA

	// fetch Picard Tools jarfile
	picard_jar = get_picard_jar()

	// load reference sequence
	reference_fasta = Channel
		.fromPath("${params.reference_fasta}")
		.ifEmpty{ error "No reference sequence found at ${params.reference_fasta}." }

	// load reference annotation
	reference_gtf = Channel
		.fromPath("${params.reference_gtf}")
		.ifEmpty{ error "No reference annotation found at ${params.reference_gtf}." }

	// ensure all reference FASTA lines are of equal length
	normalized_fasta = normalize_fasta(picard_jar, reference_fasta)

	// create a dictionary for the reference sequence
	reference_dict = create_sequence_dictionary(picard_jar, normalized_fasta)

	// create Drop-seq tools refFlat annotation file
	reference_refFlat = convert_to_refFlat(reference_gtf, reference_dict)

	// create reduced GTF annotation file
	reduced_gtf = reduce_gtf(reference_gtf, reference_dict)

	// create intervals files from reduced GTF
	intervals = create_intervals(reference_dict, reduced_gtf)

	// convert intervals files to BED format
	beds = convert_intervals_to_bed(picard_jar, intervals)

	// add mitochondrial bed file, if intervals file exists
	if ( file("reference.mt.intervals").exists() )
		beds = beds.concat(convert_intervals_to_bed_mt(picard_jar, intervals))


	// PROCESS READS

	// load paired-end FASTQ read files
	raw_read_pairs = Channel
		.fromFilePairs("${params.fastq_directory}/*_{1,2}*.fastq{.gz,}")
		.ifEmpty{ error "No read pairs found in ${params.fastq_directory}." }
		.map{ sample_name, reads -> tuple(sample_name, reads[0], reads[1]) }

	// convert reads from FASTQ to BAM
	converted_bam = convert_reads(picard_jar, raw_read_pairs)

	// tag reads with barcodes; filter and trim reads
	prepared_bam = prepare_reads(converted_bam)

	// load STAR index
	star_index = Channel
		.fromPath("${params.star_index}", checkIfExists: true)

	// align and annotate reads
	aligned_bam = align_annotate_reads(
		picard_jar,
		prepared_bam,
		star_index,
		normalized_fasta,
		reference_dict,
		reference_refFlat)

	// repair read barcodes
	repaired_bam = repair_reads(aligned_bam)


	// COMPUTE QC METRICS

}

