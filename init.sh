source activate gmatic
conda env export > doc/environment.yml

if [ ! -d fastqc ]; then
	mkdir fastqc clean bam track table figure RData
fi
