untrimmed=(
20151030_seq/KA5146901
);

write_stats_file="20151021_seq/stats_simple.csv"
ext=".fastq.gz"
counter=1

printf "\n****************************************************************\n"
printf "************ PCR Bias Test for NextFLEX 4N Libraries ***********\n"
printf "****************************************************************\n\n"

for fq in ${untrimmed[*]}
do
	#printf "  ($counter of ${#untrimmed[*]}) Trimming reads (50%%) . . . . .            \r"
	cutadapt -a TGGAATTCTCGGGTGCCAAG -o $fq.trim1.fq --too-short-output $fq.short --discard-untrimmed --minimum-length 24 $fq$ext &> $fq.cutadapt.out	
	#printf "  ($counter of ${#untrimmed[*]}) Performing PCR bias test (100%%) . . . . . \r"
	printf "  ($counter of ${#untrimmed[*]}) Performing PCR bias test (100%%) . . . . . \n"
	python pcr_bias_test.py $fq.trim1.fq
	printf "\n"
	#printf "  FastQ file $counter of ${#untrimmed[*]} processed successfully!              \n"
	counter=$(($counter+1))
done
printf "\n"