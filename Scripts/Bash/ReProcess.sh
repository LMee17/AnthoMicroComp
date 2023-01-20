#usr/bin/bash

#USAGE: bash Reprocess.sh <directory of fastq files to map> <sample/genome Key> <no threads>

#Script to rerun the mapping part of the pipeline, with relaxed mapping parameters

#NOTES
#Input is considered to be a directory of fastq files with its batch ID as directory name
#Sample/Genome key is two columned input = col 1: Sample SRR IDs col2: closest genome to sample is to be mapped against
#No threads = number of threads as input to STAR 

#make a list of the input samples
ls $1*f*q | awk -F "/" '{print$NF}' | awk -F "_" '{print$1}'  | awk -F "." '{print$1}' | uniq > samples.tmp

#extract batch name from the input directory
batch=$(echo $1 | sed 's/.$//')
echo $batch

#make output directories ready to populate
mkdir $1"$batch"_Reads/
mkdir $1Logs/

#run through all samples to allow user to check for issues
while read s; do
	echo $s
done < samples.tmp

#get input from the user to confirm that formatting so far is correct
echo "Are these the correct samples and are we ready to continue?"
read reply
if [[ $reply == y* ]]; then
	echo "Great. Onwards"
else
	echo "Please check the read directory and try again."
	exit
fi

#for each sample,..
while read s; do
	echo $s
	#gen = genome to be mapped against, second column of the sample/genome key
	gen=$(grep $s $2 | awk '{print$2}')
	echo $gen
	#list number of fastq files (ie paired or single)
	readz=$(ls $1$s* | wc -l)
	echo $readz
	#if the sample is paired reads...
	if [[ readz -gt 1 ]]; then
		#make sure paired reads are in the right order
		seqkit pair -1 $1$s*1.fastq -2 $1$s*2.fastq -u
		#run mapping with relaxed parameters (
		STAR --runThreadN $3 --genomeDir ../$gen --outReadsUnmapped Fastx --outFileNamePrefix $s. --readFilesIn $1$s*1.paired.fastq $1$s*2.paired.fastq --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
		#pair these
		seqkit pair -1 $s*.out.mate1 -2 $s*.out.mate2 -u
		mv $s*.out.paired.mate1 $s.unmapped_R1.fq
		mv $s*.out.paired.mate2 $s.unmapped_R2.fq
		rm $s.Unmapped.out.mate*
	else
		#map single read directly to genome and keep unmapped reads
		STAR --runThreadN $3 --genomeDir ../$gen --outReadsUnmapped Fastx --outFileNamePrefix $s. --readFilesIn $1$s*fastq --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
		mv $s*.out.mate1 $s.unmapped.fq
	fi
	#clean up
	mv *unmapped*.fq $1"$batch"_Reads/
	rm -r $s*._STARtmp/
	rm *sam
	mv $s*Log* $1Logs/
	rm $1*fastq
	rm $s.SJ.out.tab $1
done < samples.tmp

rm *tmp
