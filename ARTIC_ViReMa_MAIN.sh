#!/bin/bash
usage="ARTIC-Recombination Batch Script
        Written by Andrew Routh 2024

USAGE: ./batchscript [OPTIONS] FASTQfile

e.g. ./ARTIC_ViReMa_Main.sh /path/to/data/mydata_R1.fastq

Required Arguments:
File	Enter full path of R1 file. Expected to end in '_1.prep.fastq.gz'

Optional Arguments:
    -h show this help text

    -p Perform custom stages; select combination of P, R, C, D. No whitespace, e.g. 'PM'.
            (default = PM: Preprocess and map)
        B Perform bowtie2/pilon reconstruction
        V Do ViReMa mapping
	C Repeat ViReMa Compilation Steps - if SAM file already present (must ungzipped).
	N Normalize Recombination Rates

    -t set threads (default: 1)

    -g provide base genome
    "

GENOME='NC_045512.2.fasta'
STAGING='PV'
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
DirRoot=${FileR1%%_1.prep.fastq.gz}
Root=${DirRoot##*/}
WKDIR=$0
ScriptPath=${WKDIR%/*}'/Scripts/'
echo $File
echo $DirRoot
echo $Root

##Initial Mapping
if [[ "$STAGING" == *"B"* ]]; then
	bowtie2 -p 16 -x NC_045512.2.fasta -1 FileR1 -2 FileR2 | samtools view -F 4 -buSh - | samtools sort -@ 16 - -o $Root'_bwt2.bam'
	samtools index $Root'_bwt2.bam' 
	pilon --fix bases --genome $GENOME --flank 5 --mindepth 25 --frags $Root'_bwt2.bam' --vcf --changes --output $Root --outdir $Root'_pilon'
	GENOME=$Root'_pilon/'$Root'.fasta'
	grep PASS $Root'_pilon/'$Root'.vcf' | awk '{OFS=""}{if($5 ~ /[A,T,G,C]/)print $4, $2, $5}' > $Root'.changes.txt'
fi

##Run ViReMa
if [[ "$STAGING" == *"V"* ]]; then
	python3 $ScriptPath'ViReMa_0.32/ViReMa.py' $GENOME FileR1 $Root'_1_ViReMa.sam' --Output_Dir $Root'_1_ViReMa' -BAM -BED12 --N 1 --X 3 --MicroInDel_Length 5 --Defuzz 0 --p $THREADS --Output_Tag $Root'_1' --MaxIters 50 --Chunk 5000000 -Stranded -FuzzEntry --Coverage_Offset 10
	python3 $ScriptPath'ViReMa_0.32/ViReMa.py' $GENOME FileR2 $Root'_2_ViReMa.sam' --Output_Dir $Root'_2_ViReMa' -BAM -BED12 --N 1 --X 3 --MicroInDel_Length 5 --Defuzz 0 --p $THREADS --Output_Tag $Root'_2' --MaxIters 50 --Chunk 5000000 -Stranded -FuzzEntry --Coverage_Offset 10
fi

##Make Normalized and coordinated BED files
if [[ "$STAGING" == *"N"* ]]; then
	#python3 $ScriptPath'Combine_unstranded_annotations.py' $Root'_ViReMa/BED_Files/'$Root'_Virus_Recombination_Results.bed' $Root'_ViReMa/BED_Files/'$Root'_Virus_Recombination_Results_noDir.bed' -BED12 -Stranded
	
 	python3 Combine_unstranded_annotations.py $Root'_1_ViReMa/BED_Files/'$Root'_1_Virus_Recombination_Results.bed' temp1.bed -CountCombine
	python3 Combine_unstranded_annotations.py $Root'_2_ViReMa/BED_Files/'$Root'_2_Virus_Recombination_Results.bed' temp2.bed -CountCombine
	cat temp1.bed temp2.bed > temp3.bed
	python3 Combine_unstranded_annotations.py temp3.bed $Root'_ViReMa_comb_Recombination_Results.bed' -CountCombine
	python3 $ScriptPath'Transpose_to_WA1-Coords.py' $Root'.changes.txt' $Root'_ViReMa_comb_Recombination_Results_noDir.bed' $Root'_ViReMa_comb_Recombination_Results_noDir_WA1coords.bed'
	python3 Normalize_BED_toCount.py $Root --MicroInDel_Length 25
 
	#python3 $ScriptPath'Plot_CS_Freq.py' $Root'_ViReMa/'$Root'_ViReMa' $Root'_ViReMa/BED_Files/'$Root'_Virus_Recombination_Results_noDir_WA1coords.bed' $Root'_pilon/'$Root'.fasta' --MicroInDel_Length 25 -CoVData -Ends --MinCov 100 --MinCount 3
fi
