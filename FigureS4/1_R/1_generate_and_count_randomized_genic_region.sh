#!/bin/bash

# main script
window=1000

for SLURM_ARRAY_TASK_ID in $(seq 0 9999)
do
    Rscript shuffleRNA.r ${SLURM_ARRAY_TASK_ID}
    awk -v OFS="\t" '$1=$1' temp_${SLURM_ARRAY_TASK_ID}.bed > temp_${SLURM_ARRAY_TASK_ID}_tsv.bed
    sort-bed temp_${SLURM_ARRAY_TASK_ID}_tsv.bed > temp_${SLURM_ARRAY_TASK_ID}_tsv_sorted.bed

    bedtools closest -a hg38.rmsk.bed -b temp_${SLURM_ARRAY_TASK_ID}_tsv_sorted.bed -d -io -id -D ref > temp_${SLURM_ARRAY_TASK_ID}_upstream.tsv
    bedtools closest -a hg38.rmsk.bed -b temp_${SLURM_ARRAY_TASK_ID}_tsv_sorted.bed -d -io -iu -D ref > temp_${SLURM_ARRAY_TASK_ID}_downstream.tsv

    echo "Generated a Alu-circRNA list..!"

    python3 count.py -i ${SLURM_ARRAY_TASK_ID} -w ${window} > result_${SLURM_ARRAY_TASK_ID}.csv
    rm temp_${SLURM_ARRAY_TASK_ID}.bed 
    rm temp_${SLURM_ARRAY_TASK_ID}_tsv.bed 
    rm temp_${SLURM_ARRAY_TASK_ID}_tsv_sorted.bed
done
# Finish
echo "job is finished"
echo "$(date)"
