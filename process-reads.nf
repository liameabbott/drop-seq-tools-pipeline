
fastq_directory = params.fastq_directory
output_directory = params.output_directory

process convert_reads {
	publishDir "${output_directory}/${sample_name}/bams", \
		pattern: "*.bam", mode: 'copy', overwrite: true

	input:
	path(picard_jar)
	tuple val(sample_name), path(fastq1), path(fastq2)

	output:
	tuple val(sample_name), path("${sample_name}.converted.bam")

	"""
	java -jar ${picard_jar} FastqToSam \
	FASTQ=${fastq1} \
	FASTQ2=${fastq2} \
	OUTPUT=${sample_name}.converted.bam \
	SAMPLE_NAME=${sample_name} \
	SORT_ORDER=queryname \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp
	"""
}

process prepare_reads {
	publishDir "${output_directory}/${sample_name}/bams", \
		pattern: "${sample_name}.prepared.bam", mode: "copy", overwrite: true
	publishDir "${output_directory}/${sample_name}/summaries", \
		pattern: "${sample_name}.*.summary.txt", mode: "copy", overwrite: true

	input:
	tuple val(sample_name), path(converted_bam)

	output:
	tuple val(sample_name), path("${sample_name}.prepared.bam")

	"""
	TagBamWithReadSequenceExtended \
	INPUT=${converted_bam} \
	OUTPUT=/dev/stdout \
	SUMMARY=${sample_name}.cb-tag.summary.txt \
	BASE_RANGE=${params.cell_barcode_base_range} \
	BARCODED_READ=${params.barcoded_read} \
	DISCARD_READ=false \
	BASE_QUALITY=${params.barcode_min_base_quality} \
	NUM_BASES_BELOW_QUALITY=${params.barcode_n_bases_below_min_quality} \
	TAG_NAME=XC \
	TAG_QUALITY=XQ \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp | \

	TagBamWithReadSequenceExtended \
	INPUT=/dev/stdin \
	OUTPUT=/dev/stdout \
	SUMMARY=${sample_name}.umi-tag.summary.txt \
	BASE_RANGE=${params.umi_barcode_base_range} \
	BARCODED_READ=${params.barcoded_read} \
	DISCARD_READ=true \
	BASE_QUALITY=${params.barcode_min_base_quality} \
	NUM_BASES_BELOW_QUALITY=${params.barcode_n_bases_below_min_quality} \
	TAG_NAME=XM \
	TAG_QUALITY=XQ \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp | \

	FilterBam \
	INPUT=/dev/stdin \
	OUTPUT=/dev/stdout \
	SUMMARY=filter.summary.txt \
	TAG_REJECT=XQ \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp | \

	TrimStartingSequence \
	INPUT=/dev/stdin \
	OUTPUT=/dev/stdout \
	OUTPUT_SUMMARY=${sample_name}.trim-starting.summary.txt \
	SEQUENCE=${params.trim_starting_sequence} \
	MISMATCHES=${params.trim_starting_n_mismatches} \
	NUM_BASES=${params.trim_starting_n_bases} \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp | \

	PolyATrimmer \
	INPUT=/dev/stdin \
	OUTPUT=${sample_name}.prepared.bam \
	OUTPUT_SUMMARY=${sample_name}.trim-polyA.summary.txt \
	USE_NEW_TRIMMER=true \
	ADAPTER=${params.trim_polyA_adapter} \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process align_annotate_reads {
	publishDir "${output_directory}/${sample_name}/bams", \
		pattern: "${sample_name}.annotated.bam", mode: "copy", overwrite: true

	input:
	path(picard_jar)
	tuple val(sample_name), path(prepared_bam)
	path(star_index)
	path(normalized_fasta)
	path(reference_dict)
	path(reference_refFlat)

	output:
	tuple val(sample_name), path("${sample_name}.annotated.bam")

	"""
	java -jar ${picard_jar} SamToFastq \
	INPUT=${prepared_bam} \
	FASTQ=/dev/stdout \
	QUIET=true \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp | \

	STAR \
	--genomeDir ${star_index} \
	--runThreadN ${params.star_n_cores} \
	--readFilesIn /dev/stdin \
	--outFileNamePrefix ${sample_name}-star-alignment/ \
	--limitOutSJcollapsed ${params.star_limit_out_sj_collapsed}

	java -jar ${picard_jar} SortSam \
	INPUT=${sample_name}-star-alignment/Aligned.out.sam \
	OUTPUT=${sample_name}.sorted.bam \
	SORT_ORDER=queryname \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	java -jar ${picard_jar} MergeBamAlignment \
	REFERENCE_SEQUENCE=${normalized_fasta} \
	UNMAPPED_BAM=${prepared_bam} \
	ALIGNED_BAM=${sample_name}.sorted.bam \
	INCLUDE_SECONDARY_ALIGNMENTS=false \
	CLIP_ADAPTERS=false \
	OUTPUT=${sample_name}.merged.bam \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	TagReadWithGeneFunction \
	INPUT=${sample_name}.merged.bam \
	OUTPUT=${sample_name}.annotated.bam \
	ANNOTATIONS_FILE=${reference_refFlat} \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process repair_reads {
	publishDir "${output_directory}/${sample_name}/bams", \
		pattern: "${sample_name}.repaired.bam", mode: "copy", overwrite: true
	publishDir "${output_directory}/${sample_name}/summaries", \
		pattern: "${sample_name}.*.{report,summary,stats}.txt", mode: "copy", overwrite: true

	input:
	tuple val(sample_name), path(annotated_bam)

	"""
	DetectBeadSubstitutionErrors \
	INPUT=${annotated_bam} \
	OUTPUT=${sample_name}.substitution.bam \
	OUTPUT_REPORT=${sample_name}.substitution-repair.report.txt \
	OUTPUT_SUMMARY=${sample_name}.substitution-repair.summary.txt \
	MIN_UMIS_PER_CELL=${params.repair_min_umis_per_cell} \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	DetectBeadSynthesisErrors \
	INPUT=${sample_name}.substitution.bam \
	OUTPUT=${sample_name}.repaired.bam \
	REPORT=${sample_name}.synthesis-repair.report.txt \
	OUTPUT_STATS=${sample_name}.synthesis-repair.stats.txt \
	SUMMARY=${sample_name}.synthesis-repair.summary.txt \
	MIN_UMIS_PER_CELL=${params.repair_min_umis_per_cell} \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}