# Overview
We utilized the [QIIME 2](https://qiime2.org/) suite of tools to perform much of the sequencing processing tasks in this project. A new virtual environment for QIIME 2 v-2019.4 installation. Most of the steps in this sequence processing workflow were executed within the QIIME environment, though the initial primer trimming was performed with an updated standalone version of Cutadapt.

# Primer trimming with Cutadapt
Unjoined paired end sequences were trimmed with Cutadapt v-2.3. We used an anchored approach, though our read data are such that the 5' adapter is not expected to be present in either the forward or reverse read.
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
> `$QIIMEDIR` refers to the path to the output path for the resulting QIIME artifact outputs

```
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

These visualizations can be loaded in the [QIIME viewer website](https://view.qiime2.org/). A directory containing the individual library `.qzv` files is [available here](https://github.com/devonorourke/nhguano/tree/master/data/qiime_qzv/demux_sumry). These per-base characteristics informed the next stage of the sequence processing: setting trimming parameters for DADA2 denoising.

# Denoising with DADA2
The per-base visualizations suggested that the entirety of the expected read length was of high quality. Because our expected amplicon length was just 181 bp, we set the 3' truncation length for a given read to 181 bases. Files containing the denoised representative sequences (`"$LIB".raw_repSeqs.qza`), a table of the read abundances per sample per ASV (i.e. an OTU table), and a summary file containing the per-sample read abundances (`"$LIB".denoisingStats.qza`) were generated for each library. The summary `.qza` file was used as input to create the similarly named `.qzv` file used as input for visualization in the [qiime viewer](view.qiime2.org):
```
## generate repseqs
qiime dada2 denoise-paired \
  --p-n-threads 24 --p-trunc-len-f 181 --p-trunc-len-r 181 \
  --i-demultiplexed-seqs "$LIB".demux.qza --o-denoising-stats "$LIB".denoisingStats.qza \
  --o-table "$LIB".raw_table.qza --o-representative-sequences "$LIB".raw_repSeqs.qza \

## generate summary visualization
qiime metadata tabulate \
--m-input-file "$LIB".denoisingStats.qza --o-visualization "$LIB".denoisingStats.qzv  
```

- DADA2-summary visualization files are [available here](https://github.com/devonorourke/nhguano/tree/master/data/qiime_qzv/dada_sumry)  
- DADA2-summary stat files are [available here](https://github.com/devonorourke/nhguano/tree/master/data/qiime_qza/dada2_denoisingStats)  

## Combining DADA2 datasets
Because each library was separately processed in DADA2 we combined all ASV table and representative sequence `.qza` file into a single pair of artifacts:
> `$PFX` refers to the parent directory path to the individual libraries with DADA2-processed .qza tables and representative sequence files

```
# tables
qiime feature-table merge --i-tables "$PFX"/lib12/p12.raw_table.qza \
  --i-tables "$PFX"/lib31/p31.raw_table.qza --i-tables "$PFX"/lib32/p32.raw_table.qza \
  --i-tables "$PFX"/lib41/p41.raw_table.qza --i-tables "$PFX"/lib42/p42.raw_table.qza \
  --i-tables "$PFX"/lib51/p51.raw_table.qza --i-tables "$PFX"/lib52/p52.raw_table.qza \
  --i-tables "$PFX"/lib71/p71.raw_table.qza --i-tables "$PFX"/lib72/p72.raw_table.qza \
  --o-merged-table tmp.raw_table.qza

# sequences
qiime feature-table merge-seqs --i-data "$PFX"/lib12/p12.raw_repSeqs.qza \
  --i-data "$PFX"/lib31/p31.raw_repSeqs.qza --i-data "$PFX"/lib32/p32.raw_repSeqs.qza \
  --i-data "$PFX"/lib41/p41.raw_repSeqs.qza --i-data "$PFX"/lib42/p42.raw_repSeqs.qza \
  --i-data "$PFX"/lib51/p51.raw_repSeqs.qza --i-data "$PFX"/lib52/p52.raw_repSeqs.qza \
  --i-data "$PFX"/lib71/p71.raw_repSeqs.qza --i-data "$PFX"/lib72/p72.raw_repSeqs.qza \
  --o-merged-data tmp.raw_repSeqs.qza
```

The `.qza` files for ASV data:
- ASV table: [tmp.raw_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) 
- ASV sequence: [tmp.raw_repSeqs.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza)

In addition, the fasta file is exported in a traditional text format:
- ASV sequence (text): [allSamps_ASVseqs.fasta.gz](https://github.com/devonorourke/nhguano/blob/master/data/text_tables/asv_data/allSamps_ASVseqs.fasta.gz)

## Next steps in analysis
1. The `tmp.raw_table.qza` file served as input into the contamination overview outlined in the [decontam workflow document](https://github.com/devonorourke/nhguano/blob/master/docs/decontam_workflow.md). That analysis included evaluating sequence variants for potential wet-bench cross contamination and sequencing platform contamination. 
2. Because we lacked evidence of pervasive contamination, we proceeded to clustering the exact sequence variants (ASVs) into clusters of representative variants sharing (at least) 98.5% similarity. These clusters were then classified using both alignment and kmer-based approaches, allowing the identification of bat-host species identification, and arthropod prey items. These OTUs (and their taxonomic assignments) were used in the subsequent diversity analyses presented in the manuscript. These analyses are described in the [diversity analyses](https://github.com/devonorourke/nhguano/blob/master/docs/diversity_analyses.md) document.
