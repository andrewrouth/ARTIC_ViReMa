ARTIC-Recombination Batch Script
Written by Andrew Routh 2024

USAGE: ./batchscript [OPTIONS] FASTQfile1 FASTQfile2

e.g. ./ARTIC_ViReMa.sh /path/to/data/mydata_R1_001.fastq.gz /path/to/data/mydata_R2_001.fastq.gz

Required Arguments:
    FASTQfile1	Enter full path of R1 file
    FASTQfile2	Enter full path of R2 file

Optional Arguments:
    -h show this help text

    -p Perform custom stages; select combination of P, M, R, C, D. No whitespace, e.g. 'PM'.
            (default = PM: Preprocess and map)
        P Only perform data preprocessing
        M Only default data mapping
        V Do ViReMa mapping

    -t set threads (default: 1)

    -g provide base genome