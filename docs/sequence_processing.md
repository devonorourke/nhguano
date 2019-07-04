# Overview
We utilized the [QIIME 2](https://qiime2.org/) suite of tools to perform much of the sequencing processing tasks in this project. A new virtual environment for QIIME 2 v-2019.4 installation. Most of the steps in this sequence processing workflow were executed within the QIIME environment, though the initial primer trimming was performed with an updated standalone version of Cutadapt.

# Primer trimming with Cutadapt
Unjoined paired end sequences were trimmed with Cutadapt v-2.3:
> `$RAWDIR` points to one of the (1 of 9 libraries) directories of sequence data
> `$OUTDIR` points to an output directory where trimmed reads were stored

```
for SAMPLE in $(ls $RAWDIR | cut -f 1 -d '_' | sort -u); do
  cutadapt --cores=24 \
  -a "GGTCAACAAATCATAAAGATATTGG;optional...GGATTTGGAAATTGATTAGTWCC" \
  -A "GGWACTAATCAATTTCCAAATCC;optional...CCAATATCTTTATGATTTGTTGACC" \
  --minimum-length 160 --maximum-length 200 \
  -o "$OUTDIR"/"$SAMPLE"_1.fastq.gz -p "$OUTDIR"/"$SAMPLE"_2.fastq.gz \
  "$RAWDIR"/"$SAMPLE"_L002_R1_001.fastq.gz "$RAWDIR"/"$SAMPLE"_L002_R2_001.fastq.gz;
done
```

These raw data were then imported as a QIIME-formatted artifact. This required creating a `manifest file`:
> `$OUTDIR` refers to the path set from previous command (where Cutadapt output fq files were output)
> `$LIB` refers to the individual library name a sample was sequenced within (e.g. "Lib12")
> `$LIBALT` refers to the individual library directory a sample was sequenced within (e.g. /path/to/Lib12)

```
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
```

This manifest file was then used in the QIIME import function to create a `.qza` artifact containing all the unjoined paired end data:
```
QIIMEDIR=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/individ_libs/"$LIBALT"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../"$LIB".manifest.file \
  --output-path "$QIIMEDIR"/"$LIB".demux.qza \
  --input-format PairedEndFastqManifestPhred33
```

We then visualized the per-base sequence quality using a subset of the merged data:
```
qiime demux summarize \
  --i-data "$QIIMEDIR"/"$LIB".demux.qza \
  --p-n 50000 \
  --o-visualization "$QIIMEDIR"/"$LIB".demux_sumry.qzv
```

These visualizations can be loaded in the [QIIME viewer website](https://view.qiime2.org/).
# Denoising with DADA2

Drop out the shell script. paste here.
Mention that .qzv files indicated that some NTCs remained; some even had moderate read depths (x samples with > 1000 reads). Link to qza files in Repo.

# QIIME 2 sequence processing
The steps in sequence processing began with merging the individually-denoised DADA2 ASV tables and representative sequence artifact files,


1. Merging all DADA2 libraries
2. Running decontam R script. Point to decontam_workflow.md file (link).
3. Filtering out non NH-bat samples, NTCs, mock samples.
```
## filtering table
qiime feature-table filter-samples \
  --i-table all.raw_table.qza --o-filtered-table study.raw_table.qza \
  --m-metadata-file "$META" --p-where "SampleID='subject-1'"

## filtering repSeqs
qiime feature-table filter-seqs \
  --i-data all.raw_repSeqs.qza --i-table study.raw_table.qza --o-filtered-data study.raw_repSeqs.qza
```

4. Figure out sampling depth with alpha rarefaction viz. Using default 10 iterations, but setting 1000 reads as minimum.

```
qiime diversity alpha-rarefaction \
  --p-metrics observed_otus --p-min-depth 1000 --p-max-depth 10000 --p-iterations 10 \
  --i-table study.raw_table.qza --o-visualization study.raw.alphaRareViz.qzv
```
