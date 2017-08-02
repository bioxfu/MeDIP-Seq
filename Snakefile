configfile: "config.yaml"

rule all:
	input:
		expand('fastqc/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
		expand('bam/{sample}.bam', sample=config['samples']),
		expand('bam/{sample}.bam.bai', sample=config['samples']),
		expand('bam/{sample}.bamqc', sample=config['samples']),
		expand('track/{sample}.tdf', sample=config['samples']),
		'RData/MeDIP_edgeR_DMR_anno.RData',
		'table/MeDIP_edgeR_DMR_anno.xlsx',
		'doc/report.html'

rule fastqc:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/{sample}_R1_fastqc.html',
		'fastqc/{sample}_R2_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc {input}'

rule trimmomatic:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic PE -threads 3 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule bowtie2:
	input:
		r1 = 'clean/{sample}_R1_paired.fastq.gz',
		r2 = 'clean/{sample}_R2_paired.fastq.gz'
	output:
		bam = 'bam/{sample}.bam'
	params:
		prefix = 'bam/{sample}',
		cpu = config['cpu'],
		index = config['index']
	shell:
		"bowtie2 -p {params.cpu} -x {params.index} -1 {input.r1} -2 {input.r2}|samtools view -Sh -q 30 -F 4 -|grep -v 'XS:'|samtools view -Shub|samtools sort - -T {params.prefix} -o {output.bam}"

rule bam_idx:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bai = 'bam/{sample}.bam.bai'
	shell:
		'samtools index {input.bam} {output.bai}'

rule bam_qc:
	input:
		bam = 'bam/{sample}.bam'
	output:
		bamqc = 'bam/{sample}.bamqc'
	params:
		cpu = config['cpu']
	shell:
		"qualimap bamqc --java-mem-size=10G -nt {params.cpu} -bam {input.bam} -outdir {output.bamqc}"

rule bam2count:
	input:
		bam = 'bam/{sample}.bam'
	output:
		cnt = 'bam/{sample}.cnt'
	shell:
		"samtools view -c -F 4 {input.bam} > {output.cnt}"
		
rule bam2bedgraph:
	input:
		bam = 'bam/{sample}.bam',
		cnt = 'bam/{sample}.cnt'
	output:
		bg = 'track/{sample}.bedgraph'
	shell:
		"X=`awk '{{print 1/$1*1000000}}' {input.cnt}`; "
		"bedtools genomecov -ibam {input.bam} -bga -scale $X -split > {output.bg}"

rule bedgraph2tdf:
	input:
		bg = 'track/{sample}.bedgraph'
	output:
		tdf = 'track/{sample}.tdf'
	params:
		IGV = config['IGV']
	shell:
		"igvtools toTDF {input.bg} {output.tdf} {params.IGV}"

rule MeDIP_DMR:
	input:
		['bam/{sample}.bam'.format(sample=x) for x in config['samples']]
	output:
		'RData/MeDIP_edgeR_DMR.RData'
	shell:
		'Rscript script/MeDIP_DMR.R {output}'

rule MeDIP_DMR_anno:
	input:
		'RData/MeDIP_edgeR_DMR.RData'
	output:
		'RData/MeDIP_edgeR_DMR_anno.RData',
		'table/MeDIP_edgeR_DMR_anno.xlsx'
	params:
		txdb = config['txdb'],
		gene2go = config['gene2go']
	shell:
		'Rscript script/MeDIP_DMR_anno.R {input} {params} {output}'

rule report:
	input:
		rmd = 'doc/report.Rmd'
	output:
		html = 'doc/report.html'
	params:
		sample = config['samples'][0]
	shell:
		"snakemake --config samples='{params.sample}' --dag | dot -Tsvg > doc/workflow.svg; "
		"Rscript script/make_report.R {input}"

