#!/bin/bash

#SBATCH --cpus-per-task=24
#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/paper3/fastq/
#SBATCH --partition=macmanes,shared
#SBATCH --job-name="lib12"
#SBATCH --output="p12.cutadapt_and_import.log"

module purge
module load anaconda/colsa
conda activate arrR

## need to modify output names for the `L00*` value for UNH-seq'd to match 001/002; all NAU are 001
################################################################################
## Take raw reads and trim with Cutadapt using an anchored adapters
################################################################################
LIB=p12                                                                             ## edit here
LIBALT=lib12                                                                        ## edit here

RAWDIR=/mnt/lustre/macmaneslab/devon/guano/Data/individLibs/"$LIB"/fastq
OUTDIR=/mnt/lustre/macmaneslab/devon/guano/paper3/fastq/"$LIBALT"

for SAMPLE in $(ls $RAWDIR | cut -f 1 -d '_' | sort -u); do
  cutadapt --cores=24 \
  -a "GGTCAACAAATCATAAAGATATTGG;optional...GGATTTGGAAATTGATTAGTWCC" \
  -A "GGWACTAATCAATTTCCAAATCC;optional...CCAATATCTTTATGATTTGTTGACC" \
  --minimum-length 160 --maximum-length 200 \
  -o "$OUTDIR"/"$SAMPLE"_1.fastq.gz -p "$OUTDIR"/"$SAMPLE"_2.fastq.gz \
  "$RAWDIR"/"$SAMPLE"_L002_R1_001.fastq.gz "$RAWDIR"/"$SAMPLE"_L002_R2_001.fastq.gz;
done

################################################################################
## Create manifest file and import into QIIME
################################################################################
conda deactivate
conda activate qiime2-2019.4

LIB=p12                                                                             ## edit here
LIBALT=lib12                                                                        ## edit here

OUTDIR=/mnt/lustre/macmaneslab/devon/guano/paper3/fastq/"$LIBALT"
cd $OUTDIR

## create manifest file
pwd > "$LIB"_pwd.tmptxt
find . -name "*.gz" | sort -u | cut -d '/' -f 2 > "$LIB"_filenames.tmptxt
cut -f 1 -d "_" "$LIB"_filenames.tmptxt > "$LIB"_col1.tmptxt
wc -l "$LIB"_col1.tmptxt | cut -f 1 -d ' ' > "$LIB"_lines.tmptxt
seq $(echo $(cat "$LIB"_lines.tmptxt)) | xargs -Iz echo $(cat "$LIB"_pwd.tmptxt) > "$LIB"_dirpath.tmptxt
paste "$LIB"_dirpath.tmptxt "$LIB"_filenames.tmptxt -d "/" > "$LIB"_col2.tmptxt
for i in $(cat "$LIB"_filenames.tmptxt); do
if [[ $i == *_1.fastq.gz ]];
then
  echo "forward"
else
  echo "reverse"
fi;
done > "$LIB"_col3.tmptxt
paste "$LIB"_col1.tmptxt "$LIB"_col2.tmptxt "$LIB"_col3.tmptxt -d "," > "$LIB"_manifest.tmptxt
echo 'sample-id,absolute-filepath,direction' | cat - "$LIB"_manifest.tmptxt > ../$(echo "$LIB").manifest.file
rm *.tmptxt

## import into qiime
QIIMEDIR=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/individ_libs/"$LIBALT"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../"$LIB".manifest.file \
  --output-path "$QIIMEDIR"/"$LIB".demux.qza \
  --input-format PairedEndFastqManifestPhred33

################################################################################
## Generate .qzv file to visualize per-base quality scores among trimmed dataset
## These data are useful for setting library-specific denoising parameters
################################################################################

qiime demux summarize \
  --i-data "$QIIMEDIR"/"$LIB".demux.qza \
  --p-n 50000 \
  --o-visualization "$QIIMEDIR"/"$LIB".demux_sumry.qzv
