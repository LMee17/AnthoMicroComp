#user/bin/bash

#USAGE bash CZidUpload.sh <directory of reads> <Name of Project on Czid> <Czid metadata file>

#Script to upload a batch of unmapped reads to CZID.org
#Assumed that directory of reads will be in side a larger batch directory and that
#this script is ran from above that directory (ie read path = Batch/Batch_Reads/"
#My directories were names BatchID_reads and that is also assumed down the script
#If this is not the case remove the extra awk steps when setting batch variable.

#logging in to CZ ID must be completed before the script is commenced using the below command
# > czid login
# or
# > czid login --persistent

#check that user is logged into czid
echo "Have you logged into CZID?"
read reply
#assumes that any affirmative replies will begin with "y" (ie yes, yeah, y)
if [[ $reply == y* ]]; then
        echo "Great. Onwards"
else
#exit script if this is not the case
        echo "This won't work until you're logged in."
        echo "Please login using czid login"
        exit
fi

#extract batch name
batch=$(echo $1 | awk -F "/" '{print$2}' | awk -F "_" '{print$1}')
echo $batch

#compile list of samples from the reads in the batch directory
ls $1*fq | awk -F "/" '{print$NF}' | awk -F "." '{print$1}' | uniq > samples.txt

#iterate through each individual sample's unmapped read fastq file
while read s; do
	#store the sample ID
	sam=$(echo $r | awk -F "/" '{print$NF}' | awk -F "." '{print$1}')
	echo "Uploading $sam ..."
	#check for paired reads / single read
	count=$(ls $1$s*fq | wc -l)
	if [ $count -eq 2 ]; then
		#upload to CZID project
		#the output of this process is saved into a log file to look for issues should uploading fail
		echo "Uploading $s (paired reads) to CZID..."
		czid short-read-mngs upload-sample -p "$2" -s "$s" --metadata-csv $3 $1$s*R1.f*q $1$s*R2.f*q
		#the most common issue is czid.org logging out the user whilst uploading. Here we check for this happening
                czcheck=$(tail -n 1 $batch.log | grep -c "not authenticated with czid")
                #if this has occurred ...
                if [[ $czcheck -eq 1 ]]; then
                	echo "CZID has logged out the user. Halting batch script."
                        #log the sample that was half - processed in the Failed output file with the apt reason
                        printf "%s\t%s\n" $s "CZID log-in timed out" >> "$batch"_Failed.txt
                        #log this failed file, and any that were set to follow it, into another batch file
                        #then exit the script
          	        echo $s > $batch.UnProcessed.txt
                        sed "0,/$s/d" samples.txt >> $batch.UnProcessed.txt
                        exit
              	else
			#assume upload was successful
                        echo "$s successfully uploaded to CZID"
                        echo "Removing sample to conserve SSD..."
                        #remove the fastq files from the working directory
                        #please note .sra files will persist in the SRA directory and should be removed after this stage of the process is complete
                       	rm $s*
		fi
	else
		#uploading single read fastq files
        	echo "Uploading $s (single read) to CZID..."
           	czid short-read-mngs upload-sample  -p "$2" -s "$s" --metadata-csv $3 $1$s*.f*q >> $batch.log
                czcheck=$(tail -n 1 $batch.log | grep -c "not authenticated with czid")
                if [[ $czcheck -eq 1 ]]; then
                	echo "CZID has logged out the user. Halting batch script."
                        printf "%s\t%s\n" $s "CZID log-in timed out" >> "$batch"_failed.txt
                        sed "0,/$s/d" samples.txt >> $batch.Unprocessed.txt
                        exit
              	else
                      	echo "$s successfully uploaded to CZID"
                        echo "Removing sample to conserve SSD ..."
                        rm $s*
               	fi
	fi
done < samples.txt

#clean up
rm samples.txt

#script run complete
