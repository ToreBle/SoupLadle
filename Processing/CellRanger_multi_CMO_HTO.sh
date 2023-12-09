#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=48G
#SBATCH --job-name=cellRanger
#SBATCH --time=48:00:00
#SBATCH --output=output.%J.out


#Cellranger multi for CellPlex demultiplexing
cellranger multi --id KH108_CMO --csv ./Samples/CMO_CellPLex_KH111.csv

#Cellranger multi for Hashtag demultiplexing
cellranger multi --id KH108_HTO --csv ./Samples/HTO_Hashtag_KH112.csv
