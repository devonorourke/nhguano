#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/paper3/qiime/individ_libs/lib12
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="lib12_dada2"
#SBATCH --output="p12.dada2.log"

module purge
module load anaconda/colsa
conda activate qiime2-2019.4

################################################################################
## Denoise Cutadapt-trimmed reads with DADA2
################################################################################

LIB=p12                                                                             ## edit here

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$LIB".demux.qza \
  --p-trim-left-f 5 \
  --p-trim-left-r 5 \
  --p-trunc-len-f 175 \
  --p-trunc-len-r 175 \
  --p-n-threads 24 \
  --o-table "$LIB".raw_table.qza \
  --o-representative-sequences "$LIB".raw_repSeqs.qza \
  --o-denoising-stats "$LIB".denoisingStats.qza

qiime metadata tabulate \
  --m-input-file "$LIB".denoisingStats.qza \
  --o-visualization "$LIB".denoisingStats.qzv
