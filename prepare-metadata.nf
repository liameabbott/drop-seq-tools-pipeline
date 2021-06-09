
metadata_dir = "${projectDir}/test/data/metadata"

process normalize_fasta {

	input:
	path(picard_jar)
	path(reference_fasta)

	output:
	path("reference.fasta")

	"""
	java -jar ${picard_jar} NormalizeFasta \
	INPUT=${reference_fasta} \
	OUTPUT=reference.fasta \
	LINE_LENGTH=100 \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process create_sequence_dictionary {
	
	input:
	path(picard_jar)
	path(normalized_fasta)

	output:
	path("reference.dict")

	"""
	java -jar ${picard_jar} CreateSequenceDictionary \
	REFERENCE=${normalized_fasta} \
	OUTPUT=reference.dict \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process convert_to_refFlat {
	
	input:
	path(reference_gtf)
	path(reference_dict)

	output:
	path("reference.refFlat")

	"""
	ConvertToRefFlat \
	ANNOTATIONS_FILE=${reference_gtf} \
	SEQUENCE_DICTIONARY=${reference_dict} \
	OUTPUT=reference.refFlat \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process reduce_gtf {
	
	input:
	path(reference_gtf)
	path(reference_dict)

	output:
	path("reference.reduced.gtf")

	"""
	ReduceGtf \
	GTF=${reference_gtf} \
	SEQUENCE_DICTIONARY=${reference_dict} \
	OUTPUT=reference.reduced.gtf \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process create_intervals {
	
	input:
	path(reference_dict)
	path(reduced_gtf)

	output:
	path("reference.*.intervals")

	"""
	CreateIntervalsFiles \
	SEQUENCE_DICTIONARY=${reference_dict} \
	REDUCED_GTF=${reduced_gtf} \
	PREFIX=reference \
	MT_SEQUENCE=${params.mt_sequence} \
	OUTPUT=./ \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process convert_intervals_to_bed {
	
	input:
	path(picard_jar)
	path(intervals)

	output:
	path("reference.*.bed")

	"""
	java \
	-jar ${picard_jar} IntervalListToBed \
	INPUT=reference.rRNA.intervals \
	OUTPUT=reference.rRNA.bed \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	java \
	-jar ${picard_jar} IntervalListToBed \
	INPUT=reference.consensus_introns.intervals \
	OUTPUT=reference.consensus_introns.bed \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	java \
	-jar ${picard_jar} IntervalListToBed \
	INPUT=reference.exons.intervals \
	OUTPUT=reference.exons.bed \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	java \
	-jar ${picard_jar} IntervalListToBed \
	INPUT=reference.intergenic.intervals \
	OUTPUT=reference.intergenic.bed \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 

	java \
	-jar ${picard_jar} IntervalListToBed \
	INPUT=reference.genes.intervals \
	OUTPUT=reference.genes.bed \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}

process convert_intervals_to_bed_mt {

	input:
	path(picard_jar)
	path(intervals)

	output:
	path("reference.mt.bed")

	"""
	java \
	-jar ${picard_jar} IntervalListToBed \
	INPUT=reference.mt.intervals \
	OUTPUT=reference.mt.bed \
	USE_JDK_DEFLATER=true \
	USE_JDK_INFLATER=true \
	TMP_DIR=./.tmp 
	"""
}
