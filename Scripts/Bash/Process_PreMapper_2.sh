#usr/bin/bash

#USAGE: bash Process_Prep <sample and genome key> <SRA directory> <directory of genomes> <Batch Name> <no Threads>

#Script that takes a list of samples and the respective closest genome to said sample's species as input
#Then this script will download the sample from the Sequence Read Archive, unpack it, and map it against
#the reference genome indicated (assumed to be found in the directory of genomes given as argument $3 above).
#Unmapped reads are retained for the next step in the pipeline.
#Samples that have <50% reads mapped to the given genome for reasons of reads being "too short"
#will be indicated so that they can be mapped again with more relaxed parameters (ReProcess.sh)

#NOTES
#Sample and genome key must have two columns: col 1 = SRA IDs, col 2 = closest host genome.
#All genomes to be used in the map run should be kept in the directory of genomes (path provided as argument $3)
#Unmapped reads will be kept in an output directory named after $4

#take note of how many samples are to be iterated through
tot=$(wc -l $1 | awk '{print$1}')
#set counter
count=1

#make list of samples to be iterated through
awk '{print$1}' $1 > samps.tmp

#make output directory space
mkdir $4/
mkdir $4/"$4"_Reads/
mkdir $4/IndSampleLogs/
mkdir $4/BatchLogs/

#begin to write the script's report
printf "%s\t%s\t%s\n" Sample Total_Reads Unmapped_Reads >> $4.Map.report

#whilst iterating through the samples....
while read s; do
	echo "$count/$tot: Processing $s ... "
	#download SRA from NCBI
	echo "Fetching $s ..."
	prefetch $s
	#unpack to get to the fastq files
        echo "Unpacking $s..."
	fasterq-dump $2sra/$s.sra >> $3.log
	#remove .sra
	rm $2sra/$s.sra
	#mapping
	#retrieve genome information
	gen=$(grep $s $1 | awk '{print$2}')
	echo "Mapping against $gen genome ..."
	#check if there is already a index folder
	if [[ -d $gen/ ]]; then
	echo "$gen has already been indexed. Beginning mapping ..."
		#Check for single or double reads
		readz=$(ls ./$s* | wc -l)
		if [[ readz -gt 1 ]]; then
			#record the time for the script report
			tik=$(date '+%c')
			echo "$tik: Mapping $s against $gen (paired reads)...." >> $4.log
			#map fastq file against the indicated genome
			STAR --runThreadN $5 --genomeDir $gen/ --outReadsUnmapped Fastx --outFileNamePrefix $s. --readFilesIn $s*1.f*q $s*2.f*q
			#record the time of mapping completion
			tok=$(date '+%c')
			echo "$tok: Mapping complete..." >> $4.log
			echo "Mapping complete..."
			#ensure that paired reads are in correct orders
			seqkit pair -1 $s*out.mate1 -2 $s*out.mate2 -u
			mv $s*.out.paired.mate1 $s.unmapped_R1.fq
			mv $s*.out.paired.mate2 $s.unmapped_R2.fq
			#record the number of unmapped reads
			unmapR=$(grep -c "@" $s.unmapped*fq | awk -F ":" '{print$2}' | uniq)
			#check to see if an unpaired file already exists for this sample
			if [[ -f "$s*unpaired*" ]]; then
				#if it does, check to see if there is already a directory to save output
				if [[ -d "$4/Reads/Unpaired" ]]; then
					#remove the out.mate files
					rm $s.out.mate*
				else
					#if no directory exists, make one
					mkdir $4/"$4"_Reads/Unpaired/
				fi
			#clean up
			mv $s*unpaired* $4/Reads/Unpaired
                        mv $s*.out.mate* $4/Reads/Unpaired
			fi
		else
			#single reads
			#record time at beginning of mapping
			echo "$tik: Mapping $s against $gen (single read)...." >> $4.log
			STAR --runThreadN $5 --genomeDir $gen/ --outReadsUnmapped Fastx --outFileNamePrefix $s. --readFilesIn $s*f*q
			#record time at end of mapping
			tok=$(date '+%c')
                        echo "$tok: Mapping complete..." >> $4.log
			echo "Mapping complete..."
			mv $s*.out.mate1 $s.unmapped.fq
			unmapR=$(grep -c "@" $s.unmapped*fq)
		fi
	else
		#genome requires indexing before it can be mapped against
		echo "$gen has not been indexed. Beginning indexing ..."
		#record time at beginning of indexing
		tik=$(date '+%c')
		echo "$tik:Required $gen genome not indexed. Beginning index..." >> $4.log
		#make a directory ready for the index to be compiled
		mkdir $gen/
		#index the genome
		STAR --runMode genomeGenerate --runThreadN $5 --genomeDir $gen/ --genomeFastaFiles Genomes/$gen* --limitGenomeGenerateRAM 36807840469 --genomeSAindexNbases 13
		#record time as genome finishes being indexed
		tik=$(date '+%c')
		echo "$tik: Indexing complete" >> $4.log
		#Check for single or double reads
                readz=$(ls ./$s* | wc -l)
		if [[ readz -gt 1 ]]; then
			#record time before mapping is commenced
                        tik=$(date '+%c')
                        echo "$tik: Mapping $s against $gen (paired reads)...." >> $4.log
			#map fastq file against genome of choice
			STAR --runThreadN $5 --genomeDir $gen/ --outReadsUnmapped Fastx --outFileNamePrefix $s. --readFilesIn $s*1.f*q $s*2.f*q
                        #record time as mapping completes
			tok=$(date '+%c')
                        echo "$tok: Mapping complete..." >> $4.log
                        echo "Mapping complete..."
			#ensure that paired reads are in correct orders
                        seqkit pair -1 $s*out.mate1 -2 $s*out.mate2 -u
                        mv $s*.out.paired.mate1 $s.unmapped_R1.fq
                        mv $s*.out.paired.mate2 $s.unmapped_R2.fq
			unmapR=$(grep -c "@" $s.unmapped*fq | awk -F ":" '{print$2}' | uniq)
			#check to see if directories need to be created or if they already exist
                        if [[ -f "$s*unpaired*" ]]; then
                                if [[ -d "$4/Reads/Unpaired" ]]; then
                                        rm $s.out.mate*
                                else
                                        mkdir $4/Reads/Unpaired/
                                fi
                        mv $s*unpaired* $4/Reads/Unpaired
			mv $s*.out.mate* $4/Reads/Unpaired
                        fi

                else
                        echo "$tik: Mapping $s against $gen (single read)...." >> $4.log
			STAR --runThreadN $5 --genomeDir $gen/ --outReadsUnmapped Fastx --outFileNamePrefix $s. --readFilesIn $s*f*q
                        tok=$(date '+%c')
                        echo "$tok: Mapping complete..." >> $4.log
                        echo "Mapping complete..."
			mv $s*.out.mate1 $s.unmapped.fq
			unmapR=$(grep -c "@" $s.unmapped*fq)
                fi
	fi
	#record the total number of input reads for report
	totR=$(grep "Number of input reads" $s.Log.final.out | awk -F '\t' '{print$2}')
	printf "%s\t%d\t%d\n" $s $totR $unmapR >> $4.Map.report
	#check for too many unmapped reads: possibly something wrong with mapping step
	qccheck=$(grep "% of reads unmapped: too short" $s.Log.final.out | awk -F '\t' '{print$2}' | awk -F "." '{print$1}')
	if [[ $qccheck -gt 50 ]]; then
		#check to see if an 'Issues' directory already exists: if it does, do nothing, otherwise create one
		if [[ -d "$4/Issues/" ]]; then
			:
		else
			mkdir $4/Issues/
		fi
		#record where there has been a problem
		echo "There appears to have been an issue with mapping: $qccheck % of reads failed to map"
		printf "%s\t%s\n" $s "Mapping issue" >> $4.problem.log
		rm -r $s._STARtmp
		rm $s*sam
		rm $s*unmapped*
		rm $s*Unmapped*
		mv $s* $4/Issues/
	else
		echo "Mapping was ran successfully"
		mv $s.unmapped* $4/"$4"_Reads/
		rm $s*Unmapped*mate*
		rm $s*f*q
		rm -r $s._STARtmp/
		rm $s*sam
		mv $s* $4/IndSampleLogs/
	fi
	echo "$s complete"
	let count=count+1
done < samps.tmp

#Clean up
mv $4.log $4/BatchLogs/
mv $4.Map.report $4/BatchLogs/
mv $1 $4/BatchLogs/
if [[ -f "$4.problem.log" ]]; then
	mv $4.problem.log $4/BatchLogs/
fi
rm samps.tmp

echo "Batch run complete"
