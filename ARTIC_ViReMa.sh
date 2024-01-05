#!/bin/bash
usage="ARTIC-Recombination Batch Script
        Written by Andrew Routh 2024

USAGE: ./batchscript [OPTIONS] FASTQfile

e.g. ./TCS_batch /path/to/data/mydata_R1.fastq

Required Arguments:
    File	Enter full path of R1 file

Optional Arguments:
    -h show this help text

    -p Perform custom stages; select combination of P, M, R, C, D. No whitespace, e.g. 'PM'.
            (default = PM: Preprocess and map)
        B Perform bowtie2/pilon reconstruction
        M Merge R1 and R2 data
        V Do ViReMa mapping
	N Normalize Recombination Rates

    -t set threads (default: 1)

    -g provide base genome
    "

GENOME='NC_045512.2.fasta'
STAGING='PMV'
THREADS=1
while getopts 'hp:t:g:' option; do
  case "$option" in
    h ) echo "$usage"
       exit
       ;;
    p ) STAGING=$OPTARG
       ;;
    t ) THREADS=$OPTARG
       ;;
    g ) GENOME=$OPTARG
       ;;

   \? ) printf "unrecognised option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND -1))

##REQUIRED INPUT
FileR1=$1
FileR2=$2
DirRoot=${FileR1%%_R1_00*}
Root=${DirRoot##*/}
WKDIR=$0
ScriptPath=${WKDIR%/*}'/Scripts/'
echo $File
echo $DirRoot
echo $Root


##Initial Mapping
if [[ "$STAGING" == *"B"* ]]; then
	bowtie2 -p 16 -x NC_045512.2.fasta -1 $1 -2 $2 | samtools view -buSh - | samtools sort -@ 16 - -o $Root'_bwt2.bam'
	samtools index $Root'_bwt2.bam' 
	pilon --fix bases --genome $GENOME --flank 1 --mindepth 10 --frags $Root'_bwt2.bam' --vcf --changes --output $Root --outdir $Root'_pilon'
fi

#ViReMa Merge Data
if [[ "$STAGING" == *"M"* ]]; then
	gunzip -dc $1 > $Root'_merge.fastq'
	gunzip -dc $2 >> $Root'_merge.fastq'
	sed -i 's/\ 1\:N\:/_1\:N\:/g' $Root'_merge.fastq'
	sed -i 's/\ 2\:N\:/_2\:N\:/g' $Root'_merge.fastq'
	gzip $Root'_merge.fastq'
fi

##Run ViReMa
if [[ "$STAGING" == *"V"* ]]; then
	python3 $ScriptPath'ViReMa_0.29/ViReMa.py' $Root'_pilon/'$Root'.fasta' $Root'_merge.fastq.gz' $Root'_ViReMa.sam' --Output_Dir $Root'_ViReMa' -BED12 --N 2 --X 3 --MicroInDel_Length 5 --Defuzz 0 --p $THREADS --Output_Tag $Root --MaxIters 25 -Overwrite --Chunk 5000000
fi
	
if [[ "$STAGING" == *"N"* ]]; then
	python3 $ScriptPath'Transpose_to_WA1-Coords.py' $Root'_pilon/'$Root'.changes' $Root'_ViReMa/BED_Files/'$Root'_Virus_Recombination_Results.bed' $Root'_ViReMa/BED_Files/'$Root'_Virus_Recombination_Results_WA1coords.bed'
	python3 $ScriptPath'Plot_CS_Freq.py' $Root'_ViReMa/'$Root'_ViReMa' $Root'_ViReMa/BED_Files/'$Root'_Virus_Recombination_Results.bed' $Root'_pilon/'$Root'.fasta' --MicroInDel_Length 25 -CoVData -Ends --MinCov 100 --MinCount 3
fi



