#!/usr/bin/env nextflow

params.cpus = 8
params.bowtie_index = "$baseDir/data/FN433596"

params.read = "$baseDir/data"
params.results = "$baseDir/results"
params.out = "$params.results"

outDir = file(params.out)

readChannel = Channel.fromFilePairs("${params.read}/*{1,2}.fastq.gz")

process mapping {

	conda "bioconda::bowtie2"

	publishDir "$outDir", mode: "copy"

	input:
	set pair_id, file(reads) from readChannel

	output:
	set pair_id, file("*.sam") into mappingChannel

	script:
	"""
	bowtie2 -q -1 ${reads[0]} -2 ${reads[1]} -x ${params.bowtie_index} -S ${pair_id}.sam -p ${params.cpus}
	"""
}


process samtools_view {

	conda "bioconda::samtools"

	publishDir "$outDir", mode: "copy"
	
	input:
	set pair_id, file(reads) from mappingChannel

	output:
	set pair_id, file("*.bam") into bamChannel

	script:
	"""
	samtools view -S -@ ${params.cpus} -b -o output.bam ${reads[0]}
	"""
}


process samtools_sort {

	conda "bioconda::samtools"

	publishDir "$outDir", mode: "copy"
	
	input:
	set pair_id, file(reads) from bamChannel

	output:
	set pair_id, file("*.bam") into sortbamChannel

	script:
	"""
	samtools sort -@ ${params.cpus} -o sorted_output.bam ${reads[0]}
	"""
}


process bedtools {

	conda "bioconda::bedtools"

	publishDir "$outDir", mode: "copy"
	
	input:
	set pair_id, file(reads) from sortbamChannel

	output:
	set pair_id, file("*.gcbout") into bedtoolsChannel

	script:
	"""
	bedtools genomecov -ibam ${reads[0]} -d > output.gcbout
	"""
}


process coverageStats {

	conda "python=3.6 numpy"

	publishDir "$outDir", mode: "copy"
	
	input:
	set pair_id, file(reads) from bedtoolsChannel

	output:
    stdout gc_result

	script:
	"""
	bed2coverage ${reads[0]}
	"""
}

gc_result.subscribe { println it }