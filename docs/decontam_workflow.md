# Contamination overview
Raw sequence data was processed as described in the [sequence_processing](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md) document. Demultiplexed paired-end read data from multiple libraries were trimmed with Cutadapt, denoised with DADA2, and merged into a single QIIME-formatted ASV table and fasta file. These data&mdash;the QIIME-formatted [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) and the [fasta](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/repSeqs/all.raw_repSeqs.qza) files merged from all libraries&mdash;served as input for these analyses.

We generated sequencing data from 9 Illumina libraries. These data consisted of true samples, positive biological mock samples, and negative controls (DNA extraction blanks that were amplified in conjunction with true samples). Because negative controls (herein termed NTCs) generated some sequencing data we explored whether there was evidence for pervasive reagent contamination, specific sequencing batch contamination, or particular amplicon sequence variants (ASVs) that were likely contaminants. Notably, contamination can originate in two distinct ways:
1. **Platform-based contamination** Also known as "index-bleed", or "cross-talk", these data are a result of the sequencing of the libraries whereby the indexes assigned to one sample are mistakenly assigned to a different sample. As a result, the ASVs present in one sample are derived from a different sample. Because we included a biological mock community in every library, we were able to empirically evaluate the likelihood of such cross talk (insofar for that the one mock sample is representative of the entire community of samples).
2. **Wet-bench contamination** can arise from the various molecular workflows conducted by the researcher in preparing the libraries: DNA extraction and PCR amplification are the principal sources of potential contamination in our experiments. DNA was extracted from samples using a 96-well plate format that has multiple opportunities for contamination: the initial bead-beating step requires the user to place individual guano samples in a 96-well plate with liquid reagents and shake the plate up to 30 Hz for 20 minutes. It's possible that the silicone mat used to secure and separate individual samples within the wells leaks liquid among the wells. It's also possible that some of this liquid and debris can be mixed between wells during the step following the shaking where the user removes the silicone mat. In addition, multiple rounds of pipetting various supernatants into and out of the wells poses opportunities for multichannel tips to drag across non-target wells. We later amplified our DNA via PCR, which required additional pipetting and added another layer of reagent (principally primers and PCR master mix) contamination. Finally, reagent contamination from the kits directly (also known as "kitomes") are possible, though less likely for our system given that these are arthropod amplicons and most kitomes are suspected to have microbial contaminants.

Our contamination evaluation explores the distributions of read abundances per sample and per ASV (i.e. how many sequence counts were present per sample? per ASV among all samples?), the distribution of unique ASVs per sample, and the prevalence of individual ASVs among all samples. We applied the R library [decontam](https://github.com/benjjneb/decontam) to evaluate the likelihood of ASV contamination, as describe below:

# Negative control samples generally produce few ASVs or sequence counts per sample

While a complete [script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R) defines all programs and code required for the following analysis, we wanted to provide a more extended visualization and commentary here.

Data was processed in RStudio v-1.2.1335 with R v-3.5.3. We imported the [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) with the [qiime2R](https://github.com/jbisanz/qiime2R) package, and created a [phyloseq](https://joey711.github.io/phyloseq/) object that contained sequence counts per ASV per sample, as well as [metadata](https://github.com/devonorourke/nhguano/blob/master/data/metadata/allbat_meta.csv) information such as the sequencing library (`SeqBatch`) a sample was derived, as well as the DNA plate the guano sample was extracted (`DNAplate`). The `SampleType` was defined as either **sample** (guano sample), **ncontrol** (NTCs), or **mock** (biological mock sample).

While there were 10,797 unique sequence variants identified across all nine libraries, the majority of these were observed in just a single sample:
  - **7,587** ASVs were singletons - a sequence variant detected in just one sample. The vast majority of these singletons generated very few reads: the mean (62) or median(13) sequence counts among these singleton ASVs indicate that these variants were rare largely because they failed to generate a significant number of amplicons among any sample.  

There were many other ASVs that were also quite rare:
  - **1,189** ASVs were observed in at least 5 samples
  - **128** ASVs were observed in at least 50 samples

In fact, we observed just **39** ASVs in at least 100 samples. This is important in our understanding of contamination in general: we included 205 negative control samples among the 9 sequenced libraries, and all 205 generated at least some amount of sequence data. But this does not in and of itself indicate that these samples were contaminated with guano! We would expect some degree of cross-talk from the sequencing platform (for extended discussion of the various contamination sources and more links on the subject, check out [this blog post](http://fiererlab.org/2018/08/15/garbage-in-garbage-out-wrestling-with-contamination-in-microbial-sequencing-projects/)). What we want to understand isn't whether or not an ASV is observed in a **ncontrol** sample. What we're interested is whether the prevalence of a particular ASV is more likely to be detected in a control than a true sample. Furthermore, we'd like to better understand the distribution of read abundances of those ASVs in **ncontrol** samples: very low abundances are potentially a result of cross-talk. They may also be a result of very low-level contamination, but given that our diversity estimates include measures that incorporate read abundances, these rare and low abundant sequences are likely not have a major impact in our downstream analyses.

Therefore, we want to determine if:
1. The abundant ASVs are more likely to be observed in the **ncontrol** than guano **sample** `SampleType` data
2. The **ncontrol** samples generated relatively more or less sequence counts than true samples
3. The relative proportion of sequence counts in **mock** samples are of the _expected_ sequence types versus any unexpected ASVs

We first partitioned the samples into their respective `SampleType` groups (**sample**, **ncontrol**, and **mock**) to determine how the per-sample ASV prevalence (number of ASVs detected per sample) and per-sample read abundances were distributed. As seen in the figure below, the observed richness of ASVs (y axis) and the number of raw sequence counts (x axis) is typically far greater among true guano **samples** than the majority of **ncontrol** samples. In addition, the biological **mock** samples produce relatively high read abundances (reflecting the fact that purified samples as opposed to guano samples often amplify better) and have most sequence richness compared to many other guano samples. This is expected because mock samples were supposed to only have ~ 20 unique sequence variants.

![imagehere: contam_eval_allSampls_Counts-ASVs_scatterplot](https://github.com/devonorourke/nhguano/blob/master/figures/contam_eval_allSampls_Counts-ASVs_scatterplot.png)

A few negative control samples are concerning in that they produce large numbers of total reads per sample, however we see that the vast majority of negative control samples are producing relatively few total reads, and are generally containing fewer than 10 unique ASVs. In other words, the negative controls represent the left tail of the read depth and ASV prevalence distributions. In fact, if we filter _out_ all samples that have fewer than 500 total sequences per sample, we find that the majority of negative control samples are discarded: just **24** of the original **206** control samples generated at least 500 reads. The following plot below labels a few of the remaining NTC samples:
> all points shown are from same dataset; we simply dropped the samples with less than 500 total reads

![imagehere: contam_eval_allSampls_Counts-ASVs_scatterplot_filt500ReadMin](https://github.com/devonorourke/nhguano/blob/master/figures/contam_eval_allSampls_Counts-ASVs_scatterplot_filt500ReadMin.png)

The labels for each negative control reflect the DNA plate number and DNA well position. For example, `negoro14D04` indicates that this particular NTC was extracted from plate **14** in well **D04**. It's possible that these few outlier NTCs were the result of trace amounts of bat guano being loaded into their respective wells: a piece of guano is the size of a large rice grain, and it's quite plausible among the >2000 guano samples that a few pieces broke into a well and were not detected by the researcher. Interestingly, _plate 14_ has two of the top 5 most abundant negative control samples (and there were 38 unique DNA plates with control samples), and perhaps this is an example in which the bead beating portion of the extraction process contributed to some contamination. One strong indication of this would be whether the community composition of the observed sequence variants are more similar to each other on this plate than other NTC samples. We'll explicitly test this in this workflow later.

Overall, this broad perspective illustrates that while we have many negative control samples producing some amount of sequence data, the vast majority of these **ncontrol** samples produce very few reads or ASVs. In fact, among the 9 libraries of data, the 206 NTC samples constitute just 0.43% of the overall sequence data. Thus, if contamination is persistent, it is at very low abundances.  

# Relatively few ASVs are identified by Decontam as likely contaminants

One approach to classifying a particular ASV as a suspected contaminant would be to highlight those sequence variants that are more prevalent in control samples than in true samples. The [decontam paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0605-2) suggests that the reason is because the lack of competing DNA in control samples would allow for particular ASVs to be repeatedly amplified. In other words, the lack of variation in control samples would tend to produce the same ASVs among all control samples, while the highly variable sequence variants among true samples would tend to compete for sequence space, and thus a particular ASV would not be identified in proportionally as many true samples as in a control.   

The decontam function creates a score by generating a chi-square statistic from a 2x2 presence-absence table for each feature. The score statistic reflects the one-sided tail of the chi-square distribution, where smaller decontam scores indicate the contaminant model (an ASV is likely a contaminant) is more likely. We followed the authors suggestion of testing multiple thresholds with which an ASV would be flagged as a contaminant; a threshold of 0.5 would indicate ASVs that are more prevalent in a higher proportion of negative controls than the guano samples. As depicted in the histogram below, we found that the distribution of these decontam scores reflected that few sequence variants were likely contaminants:

![imagehere:decontam_prevHist_basic](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_prevHist_basic.png)

There are indeed _some_ ASVs in the left tail of that distribution, indicating that the program has targeted a few sequence variants as likely contaminants. We can see how many are flagged as contaminants at a given threshold in the following table:

![imagehere:decontam_prevTable_basic](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_prevTable_basic.png)

A few notes about this table:
1. Those marked as `FALSE` may be a situation where an ASV is present in both negative and true samples, but is better reflected by the non-contaminat model. It's also possible, however, that the ASV isn't detected in a negative control sample at all. These certainly wouldn't be contaminant ASVs then! In fact, there are **2,803** ASVs that are private to the true samples (that is, they are not detected in any negative control sample).
2. Those marked as `TRUE` may be a situation where an ASVs is present in both negative and true samples, but is better reflected by the contaminant model. It's also possible, however, that an ASV isn't detected in a true sample at all. These certainly wouldn't be contaminants _either_ because they aren't in our samples of concern! There are, in fact, 12 ASVs marked as `TRUE` by the decontam program that fit this description, but the majority that are flagged are indeed present in at least some number of true samples.

Notably, the Decontam package also contains a parameter to control for batch effects. In our case we wanted to test for both DNA plate ("DNAplate") and sequencing ("SeqBatch") batch effects. We found that there were similar distributions of the left tail of Decontam scores for each batch type, but unique profiles to the right of the distribution:

![imagehere:decontam_prevHist_batch](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_prevHist_batch.png)

The overall number of ASVs identified as suspected contaminants were fairly similar among the `DNAplate` and `SeqBatch` batch types:

![imagehere:decontam_prevTable_batch](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_prevTable_batch.png)

Nevertheless, the distributional differences for each batch type had me wondering whether or not the _same_ ASVs were flagged as contaminants or not. We wouldn't expect that to necessarily be true: as noted by the Decontam authors, this program is not intended to model cross-talk, which is most likely the driving force behind any `SeqBatch` effect. Instead, the `DNAplate` model is better suited to model the ASV prevalence within a set of samples that shared a common DNA extraction condition. However, the statistical power to detect such differences is diminished when grouping by this batch type: there were between 1-10 (mean 5.42, median 6) NTC samples per DNA plate analyzed.

The following table demonstrates how many ASVs are shared or are distinct to a given set. With three different _'batch'_ types to consider in the Decontam model ("basic" == no batch, "DNAplate" for the DNA extraction plate type, and "SeqBatch" for the Illumina sequencing run) there seven distinct sets (listed in the table as a "Group"). We find that most ASV identified as a contaminant are generally shared among the three different batch types at a given threshold:

![imagehere:TableOf_ASVcounts_byThreshold_byGroup.png](https://github.com/devonorourke/nhguano/blob/master/figures/TableOf_ASVcounts_byThreshold_byGroup.png)

Because there were distinct ASVs being flagged depending on the batch type, I wanted to explore _which_ ASVs those where, and how frequently they were detected among all the samples. I chose a threshold of 0.2 for the ASVs identified as contaminants in either the `SeqBatch` or `DNAplate` batch method. The following plot illustrates that the particular batch type we're employing makes a big difference if we were going to blindly remove a particular ASV from the dataset that was flagged by Decontam:

![decontam_Reads-ASVs_ContamOrNot](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_Reads-ASVs_ContamOrNot.png)

A few observations from this figure:
1. The most prevalent ASVs in our dataset are routinely identified as contaminants by the DNAplate method, but are typically not by the sequencing batch method. This makes me suspicious that there is not a sufficient statistical power to support the DNAplate-based batch filter. For example, it could be that a partiuclar ASV is marked as a contaminant in 3 of 40 DNA plates, but not so in the other 37. Nevertheless, we're seeing some ASVs in hundreds of samples, and these are likely _not_ contaminants: they're the kinds of sequence variants we expect in the bat diet (each ASV was classified and assigned taxonomic identity, explained in the last section of this document).
2. The mock ASVs sequenced were often flagged as likely contaminants themselves. Given that these samples were not part of the DNA extraction process, the more abundant samples are only a result of sequencing cross-talk (and we know Decontam isn't built for that). Given that the mock samples often generated among the greatest per-sample read abundances per Illumina run, it's not surprising that we might be seeing the mock ASVs showing up in negative control samples at low abundances. This is explored in more detail in the next section.
3. There are negative control samples that are _not_ identified as contaminants. These are clearly not of concern. It's unclear to me why an ASV would be unique to a control sample and not present in any mock sample or guano sample, but these are not of a concern with our data and will be filtered out.

Because the Decontam package works by a presence-absence detection method, I wondered whether these ASVs being marked as contaminants were potentially overly sensitive, and whether accounting for differences in read abundances may be worth investigating. A negative control sample may have low read abundances whereas a true guano sample may have higher read abundances, for instance, yet these low reads are still "detected" in many control samples. As mentioned in the beginning of this document, we expected contaminants to occur more often in negative control samples than true samples, but we're not accounting for any sort of read abundances. If it turns out that most negative controls have very low sequence abundances of these ASVs that are highly prevalent, then perhaps these aren't very concerning. Likewise, it may be that when we visualize the per-sample ASV prevalence (like in the plot above) and find that a few ASVs have large read abundances, it could be that just one or two negative control samples generated that majority of those sequences. Thus, I wanted to examine the read counts on a per-sample basis, and selected a few ASVs that were highly prevalent, and marked either as contaminants or not among the two batch methods used in Decontam. As you can see in the plot below, the ASVs that are highly prevalent in our dataset are often considered contaminants (i.e. `TRUE` in the plot below) when using the `DNAplate` batch method in Decontam, but are `FALSE` when using the `SeqBatch` method:

![imagehere:decontam_Prevalence_ncontrolANDsample](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_Prevalence_ncontrolANDsample.png)

I've highlighted 9 different ASVs above:
- The contamination status is different for 3 ASVs depending on the batch type (ASVs 2, 5, and 7)
- Three ASVs are both `FALSE` (not identified as contaminant): ASVs 15, 18, and 20
- Three ASVs are both `TRUE` (identified as contaminants in both batch types): ASVs 8, 14, and 16

It's interesting to see that there is a generally linear trend among the presence of an ASV in both control and true samples. This indicates that we are more likely to detect an ASV in a control sample when it is found in many true samples. On one hand, that would be expected if it was a pervasive contaminant, but interestingly, these ASVs in most control samples are rarely in multiple taxa. Among the hundreds of ASVs detected in these control samples, just 13 occur in at least 10 samples (out of a possible 206 control samples). This suggests that while ASVs may occur in multiple samples, it's probably not because of some pervasive contaminant, rather, it's just because those commonly observed ASVs are already in the 96-well plate and have a greater likelihood of being picked up in a negative control well during DNA extraction.

Using those same 9 ASVs, let's look at the distributions of the per-sample read abundances each of those ASVs in the among negative control samples relative to true samples:

![imagehere:decontam_selectASVabundance-perSample_contamComparison](https://github.com/devonorourke/nhguano/blob/master/figures/decontam_selectASVabundance-perSample_contamComparison.png)

1. There doesn't appear to be a large difference in read abundances among the ASVs that are suspected as contaminants using the `DNAplate` batch method in Decontam, but not identified as contaminants using the `SeqBatch` class as the batch factor (top row of "A" in panel, purple square, ASVs 2,5,7). The sequence variants that were `TRUE` for both like ASVs 8 or 14 (middle row, blue square) appear to have a few more detections than the ASVs that were considered not to be indicators of contamination for either batch group (bottom row, red square), but compared with the sheer number of true samples ("B" panel in figure) these differences are questionable. We find that ASV-16 have generally higher read abundances, but these were not classified with any confidence to a particular arthropod, and likely represent an undetected chimeric sequence (note that the DADA2 pipeline we used to denoise data does contain a denovo chimera detection step).
2. Among control samples, we rarely detect many samples with high read abundances, but these same ASVs frequently are detected in higher amounts in true samples. This indicates to me that the presence-absence framework used in our analysis may not be wrong in identifying ASVs as potential contaminants, but the these ASVs aren't likely to cause any problems by leaving them into our subsequent diversity estimates that incorporate abundance information.
3. As mentioned above, we tend to see ASVs in our negative controls with the highest prevalence when they are also among the most prevalent ASVs in true samples. What's interesting about the above plot is that it suggests that for those ASVs which were represented in many controls samples like ASV2 (57 samples), ASV5 (33 samples, and ASV7 (30 samples), we rarely see many samples with high read abundances, but there are a few specific samples that generate the majority of these reads. If there were consistent and persistent contamination, I would have suspected that ASVs would have generated more similar patterns of sequence counts. However, elevated ASVs randomly increased in single negative control samples is exactly what you'd expect in situations where an empty well has no alternative DNA tempaltes to compete with during DNA amplification. An indeed, in an example like sample `negoro14G07` which was one of the outliers noted in the earlier plot describing read abundances and ASV prevalence per sample, over 75% of all reads are attributed to `ASV2`. In fact, the majority of reads are assigned to a taxonomic group of the same Genus of beetle (_Phyllophaga_) suggesting that a small fragment of that insect either was dislodged from a guano pellet during the loading of the 96-well plate, or, that a single DNA template was amplified and generated the observed sequence variants. We observe the same effect with the negative control that generated the most sequence data - `negoro37F08` - whereby a single ASV (`ASV83`) generated over 78% of all sequences observed. Moreover, the species assigned to that ASV, _Phyllophaga crenulata_, was assigned to 12 different ASVs.

Overall, the analyses thus far have suggested that:
1. ASVs are most prevalent in negative control samples when they are also highly prevalent among true samples.
2. Negative controls rarely produce substantial read abundances collectively (per sample) or individually per sequence variant.
3. There are rare cases in which a particular ASV generates a substantial number of reads or is highly prevalent, but even among those instances the number of samples that have substantial read counts is very rare.
4. There is little evidence that there are specific sequencing run (SeqBatch) or DNA extraction (DNAplate) batch effects. Most ASVs identified by Decontam as contaminants are shared among both of these batch methods, though the DNAplate batch type generated the most number of suspected ASVs.


We wanted to next explore whether the composition of the sequence variants was associated with these batch effects. This required rarefying the data to conduct a series of distance estimates using nonphylogenetic and phylogenetic metrics.

# Contamination events are likley localized to individual plates, random across a given plate, and not pervasive within a plate
## Generating the required data: rarefying data, calculating distances, and PCoA
We'll rarefy our data prior to calculating any community composition distances. QIIME 2 has an [alpha rarefaction function](https://docs.qiime2.org/2019.4/plugins/available/diversity/alpha-rarefaction/) whereby we can visualize how the observed richness of a sample (or other alpha diversity metrics) changes with sampling depth. Because this analysis is focused on looking at how negative control samples may or may not associate with a particular group (i.e. DNA extraction plate, or Sequencing batch) we care more about preserving as many NTC samples as possible, so we're going to focus our alpha diversity analysis only on those samples. Thus, we're first going to filter out all other non-NTC samples from the `tmp.raw_table.qza` initially imported in the [decontam R script](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R), then run the function on those NTC samples for a range of sampling depths. We'll then choose a single depth, rarefy _all_ samples at that depth, then use that rarefied table to conduct the distance estimates. All of this will happen using QIIME functions.

We know that most NTCs have low abundances, so we'll investigate the observed richness of ASVs in a narrow range: 100 to 1000 sequence counts:

> `$TABLE` refers to the [tmp.raw_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) table
> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)

```
## Filtering data to retain only NTC samples
qiime feature-table filter-samples \
  --i-table $TABLE \
  --m-metadata-file $META \
  --p-where "SampleType='ncontrol'" \
  --o-filtered-table tmp.neg_table.qza

## Summary of remaining samples:
qiime feature-table summarize --i-table tmp.neg_table.qza --o-visualization tmp.neg_sumrytable.qzv

## Alpha diversity visualization
qiime diversity alpha-rarefaction \
--i-table tmp.neg_table.qza --o-visualization tmp.neg_alphaviz.qza\
--p-metrics observed_otus --p-max-depth 1000 --p-min-depth 100
```

The [tmp.neg_sumrytable.qzv](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/tmp.neg_sumrytable.qzv) file can be loaded into the [QIIME viewer online](view.qiime2.org) and there is an interactive feature that illustrates how the number of features and samples are lost as a result of the sampling depth. The [tmp.neg_alphaviz.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/tmp.neg_alphaviz.qza) file contains a different summary of the number of observed ASVs at each of the specified sampling depths.  

The `tmp.neg_sumrytable.qzv` file (visualized in the table below) suggests that a sampling depth between 400-600 sequences will retain a balance between preserving the most number of sequence variants and the greatest number of samples. However this is only relative to the negative control samples; clearly, we'd want to have a higher number of sequences for the true samples.

| Sampling depth (bp) | % ASVs remaining | # Samples remaining |
| -------------- | ------ | --------- |
| 100 | 7% | 78 |
| 200 | 8.3% | 47 |
| 300 | 9.5% | 36 |
| 400 | 11.7% | 33 |
| 500 | 10.6% | 24 |
| 600 | 11.3% | 21 |
| 700 | 11.3% | 18 |
| 800 | 9.1% | 13 |
| 900 | 8.8% | 11 |
| 1000 | 8% | 9 |

We can see that the overall observed richness doesn't change much across that range of 400-600 bases in the `tmp.neg_alphaviz.qza` file. Here's a screenshot of what that file looks like in the viewer:

![imagehere:contam_alpharareViz](https://github.com/devonorourke/nhguano/blob/master/figures/contam_alpharareViz.png)

Given these data, we selected a sampling depth of 500 reads in attempts to preserve diversity without losing too many NTC samples. We next rarefied the _entire_ dataset at that depth in QIIME. We also created another summary visualization of the remaining table to clarify which samples were retained in the analysis.

```
qiime feature-table rarefy \
  --i-table $TABLE \
  --p-sampling-depth 500 \
  --o-rarefied-table contam_rfyd_table.qza

qiime feature-table summarize --i-table contam_rfyd_table.qza --o-visualization contam_rfyd_sumrytable.qzv
```

As we can see in the [contam_rfyd_sumrytable.qzv](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/contam_rfyd_sumrytable.qzv) file, there are 1,515 samples remaining at this sampling depth, 24 of which are NTCs. The resulting [contam_rfyd_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/contam_rfyd_table.qza) ASV table contains the rarefied data, but we're going to filter out any samples with just a single ASV before we can use the data to calculate distances and ordinate the data:

```
qiime feature-table filter-samples --p-min-features 2 \
  --i-table contam_rfyd_table.qza --o-filtered-table contam_rfyd_filtd_table.qza

qiime feature-table summarize --i-table contam_rfyd_filtd_table.qza --o-visualization contam_rfyd__filtd_sumrytable.qzv
```

As we can see in the [contam_rfyd_filtd_sumrytable.qzv](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qzv/contam_analysis/contam_rfyd_filtd_sumrytable.qzv) visualization file, the filtered [contam_rfyd_filtd_table.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/contam_rfyd_filtd_table.qza) ASVtable now has just **1,498** samples instead of the original **1,515** samples. Notably, all 24 NTCs remain.  

Next, we're going to be using Unifrac to estimate community compositional differences between samples, so we need to create a tree using the ASV sequences in the dataset. We'll again rely on a QIIME function to [align ASVs to build a tree](https://docs.qiime2.org/2019.4/plugins/available/phylogeny/align-to-tree-mafft-fasttree/). This is done with the original [tmp.raw_repSeqs.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza) QIIME-formatted fasta file equivalent.

> `$READS` refers to the `tmp.raw_repSeqs.qza` file linked above

```
## create tree with FastTree
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences $READS \
--p-n-threads 24 \
--i-sequences tmp.raw_repSeqs.qza \
--p-n-threads 24 \
--o-alignment raw.ASV_alignment.qza --o-masked-alignment raw.ASV_alignment_masked.qza \
--o-tree raw.ASVtree_unrooted.qza --o-rooted-tree raw.ASVtree_rooted.qza
```

The rooted [raw.ASVtree_rooted.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/trees) file is used for the Unifrac distance estimates next.

QIIME has a pair of beta diversity functions for [non-phylogenetic](https://docs.qiime2.org/2019.4/plugins/available/diversity/beta/) and [phylogenetic](https://docs.qiime2.org/2019.4/plugins/available/diversity/beta-phylogenetic/) metrics. We're using four metrics here in this comparison:
- Dice-Sorensen (a.k.a. _Observed otus_): unweighted abundance, unweighted phylogenetic
- Bray-Curtis: weighted abundance, unweighted phylogenetic
- unweighted Unifrac: unweighted abundance, weighted phylogenetic
- weighted Unifrac: weighted abundance, weighted phylogenetic  

Each distance estimate is calculated as follows:
> `$TABLE` refers to the `contam_rfyd_filtd_table.qza` file
> `$TREE` refers to the `raw.ASVtree_rooted.qza` file

```
## non-phylogenetic
qiime diversity beta --i-table $TABLE --p-metric dice --o-distance-matrix contam_ds_dist.qza
qiime diversity beta --i-table $TABLE --p-metric braycurtis --o-distance-matrix contam_bc_dist.qza

## phylogenetic
qiime diversity beta-phylogenetic --i-table $TABLE --i-phylogeny "$TREE" --p-metric unweighted_unifrac --o-distance-matrix contam_uu_dist.qza
qiime diversity beta-phylogenetic --i-table $TABLE --i-phylogeny "$TREE" --p-metric weighted_unifrac --o-distance-matrix contam_wu_dist.qza
```

We then use each of these distance estimates in a Principal Components Analysis:
> `$DISTDIR` refers to the [directory with each distance artifact](https://github.com/devonorourke/nhguano/data/qiime_qza/distmat/contam_evals) - the `*dist.qza` files from the previous output  
> `$PCOADIR` refers to the [output director where each PCoA artifact is found](https://github.com/devonorourke/nhguano/data/qiime_qza/pcoa/contam_evals) - the `*pcoa.qza` files

```
qiime diversity pcoa --i-distance-matrix "$DISTDIR"/contam_ds_dist.qza --o-pcoa ds_pcoa.qza
qiime diversity pcoa --i-distance-matrix "$DISTDIR"/contam_bc_dist.qza --o-pcoa bc_pcoa.qza
qiime diversity pcoa --i-distance-matrix "$DISTDIR"/contam_uu_dist.qza --o-pcoa uu_pcoa.qza
qiime diversity pcoa --i-distance-matrix "$DISTDIR"/contam_wu_dist.qza --o-pcoa wu_pcoa.qza
```

These PCoA `.qza` files were imported into the [R script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R) to produce the subsequent ordinations presented in the following section. The `*dist.qza` objects were used to conduct a permutational analysis of variance using the [QIIME adonis plugin](https://docs.qiime2.org/2019.4/plugins/available/diversity/adonis/), which is a wrapper for the [Vegan Adonis R script](http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html) that allows us to assess the proportion of variation associated with a series of main effect variables. For our analysis we examined the model of (y ~ SampleType * SeqBatch * DNAplate).  

> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)
> `$ADONISDIR` refers to the directory with the [output of the Adonis directory]((https://github.com/devonorourke/nhguano/data/qiime_qzv/contam_analysis/adonis)) containing the `.qzv` visualization files

```
qiime diversity adonis \
  --i-distance-matrix contam_ds_dist.qza --o-visualization "$ADONISDIR"/contam_ds_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
qiime diversity adonis \
  --i-distance-matrix contam_bc_dist.qza --o-visualization "$ADONISDIR"/contam_bc_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
qiime diversity adonis \
  --i-distance-matrix contam_uu_dist.qza --o-visualization "$ADONISDIR"/contam_uu_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
qiime diversity adonis \
  --i-distance-matrix contam_wu_dist.qza --o-visualization "$ADONISDIR"/contam_wu_adonis.qzv \
  --m-metadata-file $META --p-formula "SampleType*SeqBatch*DNAplate"
```

## Diversity analyses
PERMANOVA reports for each of the four distance metrics are available as `.qzv` files in the [Adonis directory of the Repo]((https://github.com/devonorourke/nhguano/data/qiime_qzv/contam_analysis/adonis)) and can be viewed in the [QIIME viewer online](view.qiime2.org). We find that while there are significant main effects for all three groups (**SampleType**, **SeqBatch**, and **DNAplate**), the strength of these relationships are extremely low with a single exception: the main effect of DNAplate. This finding suggests that when contamination does persist, it is a limited to the samples within the single plate itself, thus project-wide pervasive contamination is not of a great concern.  

The per-sample dissimilarities were subsequently ordinated and visualized in the following plots. Both present the first two principle component axes for each of the four distance measures.

The first plot examines the relationship of community composition between negative control samples (purple color) and positive control samples (grey), where each label of text on the plot represents the Sequencing batch the sample was from:

![imagehere:contam_pcoa_4metric_bySeq](https://github.com/devonorourke/nhguano/blob/master/figures/contam_pcoa_4metric_bySeq.png)

The trends in these relationshps are strongest among the unweighted metrics, particularly the Dice-Sorensen metric. Because sequencing batches generally had samples that were extracted in similar locations and dates we would expect there to be some minor relationship to sequencing batches even among negative controls. However when abundances are taken into account, we find less of a relationship within sequencing batches. Further, when phylogenetic information is added with abundance information we find that there are negative control samples from each sequencing run in the same two dimensional space as another sequencing run. With so little overall variation explained by the unweighted metrics, we find no reason to be concerned about contamination due to the sequencing batch.  

The next plot highlights how negative control samples relate to their respective positive controls with respect to the DNA plate a sample was extracted from (again, purple indicates a negative control sample, and a gray sample is a guano sample):  

![imagehere:contam_pcoa_4metric_byDNA](https://github.com/devonorourke/nhguano/blob/master/figures/contam_pcoa_4metric_byDNA.png)

We find that there is a greater similarity among samples with common DNAplate numbers than with sequencing batch numbers. This is precicely what we would expect if there was a minor amount of variation occurring during DNA extraction - negative control samples that looked more like _other_ plates would be an indication of reagent contamination occurring across multiple extraction experiments and would be even more concerning. Thus we'd expect both guano and negative control samples to cluster together given that the guano samples were generally extracted in batches related to the site and week they were obtained. Nevertheless, while some samples are more similar to each other within the same DNAplate, the similarity of multiple negative control samples within a single DNAplate often vary in the same ordination space. If the entirety of an extraction was contaminated, we would expect all the negative control samples to cluster together more tightly to each other than they do to the other true samples, yet we don't have any strong evidence for that.

To show one example of this, we selected a single DNAplate (**DNAplate 33**) and illustrate that the specific well location of a negative control sample often does not match up in composition to the expected neighboring wells:

![imagehere:contam_pcoa_4metric_byDNAwel33only](https://github.com/devonorourke/nhguano/blob/master/figures/contam_pcoa_4metric_byDNAwell33only.png)

These data for just a single DNAplate indicate two further pieces of evidence suggesting that we do not need to be concerned about particular ASVs contaminating entire DNAplates:
1. We do not see that the NTC samples are clustering together, thus there are no specific ASVs that appear to be dominating the negative control composition  
2. The NTCs do not necessarily contain similar compositions to all the neighboring wells where they were extracted. If there was local contamination, we would expect that a given NTC sample would be similar to the surrounding wells. For example, the NTC sample in well location **B02** would be surrounded by plate well locations: A01, B01, C01, A02, C02, A03, B03, and C03. We find some evidence for that (ex. C02) but other instances where it's not apparent at all (ex. C01, A03, B03).

Collectively these analyses point to random and infrequent instances where some minor amount of contamination may have occurred during DNA extraction, we expect that these instances are rare and will likely not contribute to any batch effects when analyzing the various sets of samples that were extracted across different plates. Furthermore, we find little evidence for any sequencing group bias. Overall we see no reason to remove any ASVs from these samples, and will simply remove the negative control samples from subsequent analysis.


# Mock community samples suggest there is some degree of cross-talk, but it is extremely low
We included a positive control in each sequencing run. These biological mock samples consisted of about 20 unique sequence variants spanning 10 distinct arthropod Orders as described in Michelle Juisnio's [paper](https://doi.org/10.1111/1755-0998.12951). We were interested in using these mock samples to evaluate two different contamination features:
1. Reagent contamination. We do not expect the same contaminants identified in the Decontam workflow (see above) to be present in the mock samples because these mock samples were not extracted in conjunction with negative controls or true samples (the mock samples consist of equimolar pools of plasmids containing the individual sequence variants).  In addition, mock samples were amplified using similar primer constructs as the guano samples, except these reactions took place using batches of reagents that were separate from those used in the true samples.  
2. Cross talk. We _would_ expect some unexpected sequence variants to be present in the mock samples because of the sequencing process itself. These would typically be low abundance reads in the mock sample, though the particular unexpected ASVs in these mocks would most likely be derived from the most highly abundant ASVs in other samples.  

We applied the following code to perform the necessary tasks to complete this evaluation:
1. Subset the original ASV table (consisting of DADA2-filtered reads from all samples) to create a  mock-only ASV table  
2. Use QIIME's [quality control]() plugin to assign ASVs in the mock ASV table as either "expected" or "unexpected". We aligned all ASV sequences observed in each mock samle to a [fasta file containing the known mock sequences](), and considered any sequence within 99% identity as "expected", while any other sequence as "unexpected". We then created a list of the expected and unexpected ASVid's, and added applied that information to the mock ASV table in the [R script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R) to create the resulting visualizations. In addition, the samples that were verified as expected mock samples were ultimately removed from analysis (as mock sample ASVs can also result in cross-talk _into_ true samples).  
> `$ALLTABLE` refers to the initial [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) imported in this workflow containing all DADA2-processed samples
> `$ALLSEQS` refers to the initial [ASV fasta](https://github.com/devonorourke/nhguano/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza) file associated with the ASV table
> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)
> `$MOCKSEQ` refers to the QIIME-formatted mock sequences file [partialCOImock.seqs.qza](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/mock/partialCOImock.seqs.qza)

```
## filtering data table to retain only mock samples and sequences
qiime feature-table filter-samples \
  --i-table "$TABLE" \
  --m-metadata-file "$META" --p-where "SampleType='mock'" \
  --o-filtered-table mock.table.qza

## filter the original fasta file to retain only the ASVs remaining in that mock table
qiime feature-table filter-seqs \
  --i-data tmp.raw_repSeqs.qza \
  --i-table mock.table.qza \
  --o-filtered-data mock.seqs.qza

## run VSEARCH to align all ASVs in mock samples to known mock reference sequences
MOCKSEQ=/mnt/lustre/macmaneslab/devon/guano/mockFastas/partialCOImock.seqs.qza
ALLSEQS=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/tmp.raw_repSeqs.qza
qiime quality-control exclude-seqs --p-method vsearch \
  --i-query-sequences $ALLSEQS \
  --i-reference-sequences $MOCKSEQ \
  --p-perc-identity 0.98 \
  --o-sequence-hits mock.expectSeqs.qza \
  --o-sequence-misses mock.notexpectSeqs.qza
```

The ASVs identified in the [mock.expectSeqs.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/mock/mock.expectSeqs.qza) file served as input for the final section of the [R script for this workflow](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R). We assigned any ASVid identified by the above alignment as a MockASV in the following plot.

![imagehere:contam_mockIndexBleed](https://github.com/devonorourke/nhguano/blob/master/figures/contam_mockIndexBleed.png)

Interestingly, there were **347 ASVs** that aligned within 98% identity and 97% coverage** among the 10,797 ASVs in the entire dataset. However, as the plot above indicates, these ASVs were uniformly extremely low abundant sequences relative to the expected mock sequences: the highest _unexpected_ ASV contained just 14 reads, with a mean of just 8.3 reads per unexpected ASV. This is exactly what we would expect: a very small proportion of cross talk in terms of read abundance, but given the large number of unique sequence variants among the entirety of the guano samples in the dataset, it's unsurprising that a few hundred of these are misassigned to a mock sample. On a per-sample basis, there are 321 NH guano samples that contain at least one of these mock-associated ASVs, but there are always very few sequence counts attributed to any one of those 25 mock-associated ASVs. While there was a maximum of 249 sequences detected for a single mock-associated ASV in a single sample, the average number with which any one mock-associated ASV is detected in a sample is just 2 sequences. Among all 25 of these ASVs, we detect a mean of 27 sequences and a median of 15 sequences. These data suggest that while cross-talk is quite pervasive from our mock samples into the broader population of sequence samples, but this is expected given that these mock samples were typically more deeply sequenced per run than other guano samples. Moreover, the overal abundance of these mock ASVs are extremely low, suggesting that abundance-based diversity metrics will not be influenced by this degree of cross-talk.

In addition, we can clearly see that just a few of these ASVs are repeatedly detected within the mock samples themselves, and these are very likely our expected mock sequences. Just 24 ASVs are present in at least 5 of 10 mock samples, and there are just 25 distinct ASVs with at least 200 reads per sample. We used the earlier classification results to examine these 25 ASVs and found that every one matched one of the expected mock taxonomies. Thus, these sequence variants likely represent those directly from the mock sample and should be removed from the final dataset. As a second comparison we also queried these 25 ASV sequences with the nr Database in NCBI using their online BLAST application; as expected, nearly all sequences were matches for the expected taxonomies, and the two instances where they were not observed to be matches were when no suitable match was found (thus they were not the _incorrect_ assignment). The list of these [25 ASV identifiers is found here](https://github.com/devonorourke/nhguano/data/fasta/prevalentMockASVs.txt).

Collectively we find that  there is clear evidence for per-library cross-talk. though the abundances with which those sequence-based contamination events occur is very low. This suggests that we are best suited to utilize diversity metrics that incorporate abundance information in our analyses, otherwise we will likely be inflating our alpha diversity estimates, and likely inflating the similarity between samples that are otherwise less related. We therefore proceed by removing all mock and negative control samples, as well as removing the 25 ASVs identified as likely mock sequences in out dataset.

# Data filtering

## Removing control samples and ASVs
We first remove all control samples from the dataset. The [25 ASVs associated with mock samples](https://github.com/devonorourke/nhguano/data/fasta/prevalentMockASVs.txt) are also discarded, as well as any ASVs that were associated exclusively with negative controls.  
> `$ALLTABLE` refers to the initial [ASV table](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/ASVtable/tmp.raw_table.qza) imported in this workflow containing all DADA2-processed samples
> `$ALLSEQS` refers to the initial [ASV fasta](https://github.com/devonorourke/nhguano/data/qiime_qza/repSeqs/tmp.raw_repSeqs.qza) file associated with the ASV table
> `$META` refers to the QIIME-formatted metadata file [qiime_allbat_meta.tsv](https://github.com/devonorourke/nhguano/blob/master/data/metadata/qiime_allbat_meta.tsv)
> `$MOCKASV` refers to the [25 mock ASV sequences](https://github.com/devonorourke/nhguano/data/fasta/prevalentMockASVs.txt)

```
## remove all negative and positive control samples from ASV table
qiime feature-table filter-samples \
  --i-table "$ALLTABLE" --o-filtered-table tmpfilt1_table.qza \
  --m-metadata-file "$META" \
  --p-where "StudyID='oro15' OR StudyID='oro16'"  

## Remove the ASVs associated with the known mock sequences
qiime feature-table filter-features \
  --i-table tmpfilt1_table.qza --o-filtered-table tmpfilt2_table.qza \
  --m-metadata-file "$MOCKASV" --p-exclude-ids

## Drop any samples that no long have any ASVs
qiime feature-table filter-features \
  --i-table tmpfilt2_table.qza \
  --p-min-samples 1 \
  --o-filtered-table tmpfilt3_table.qza

## Drop any ASVs that no long have any samples  
qiime feature-table filter-samples \
  --i-table tmpfilt3_table.qza \
  --p-min-features 1 \
  --o-filtered-table sampleOnly_table.qza

rm tmpfilt*_table.qza
```

## Removing bat (host) sequences
We next take the [sampleOnly_table.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_table.qza) table and identify all bat host sequences in the samples. We're going to use two separate databases to query this dataset:

1. The first database consists of a selection of host references designed specifically for this project. These sequences contain references for bat species known to inhabit New England and New York. We also included other reference sequences assigned to organisms that our lab had been conducting DNA extractions with at the same time this project was being conducted. Full details describing the database design are available - see [hostCOI_database_design](https://github.com/devonorourke/nhguano/blob/master/docs/hostCOI_database_design.md).
2. A second database consisting of millions of COI sequences from arthropod and non-arthropod animals, as well as non-animal COI from subjects like fungi and microeukaryotes. The construction of this database is described in the [broadCOI_database_design.md](https://github.com/devonorourke/nhguano/blob/master/docs/broadCOI_database_design.md) file.

All files associated with these databases are hosted in an [Open Science Framework project](https://osf.io/qju3w/).

We used the QIIME 2 VSEARCH plugin to align representative sequences to the host and broad COI databases (separately) using similar parameters. We also trained a Naive Bayes classifier in QIIME 2 using the same broad COI database, then assigned taxonomy to the same ASV dataset a second time to identify if any ASVs were missed by either of the VSEARCH classification methods. We executed the following code to train the Naive Bayes classifier and classify the sequences using the Naive Bayes machine learning classifier as well as the VSEARCH alignment method:
> `$HOSTDBSEQ` refers to the QIIME-formatted sequence file for the host COI database ([host_seqs.qza](https://osf.io/p7sze/))  
> `$HOSTDBTAX` refers to the QIIME-formatted taxonomy file for the host COI database ([host_taxonomy.qza](https://osf.io/aq6ex/))  
> `$BIGDBSEQ` refers to the QIIME-formatted sequence file for the broad COI database ([bigCOI.derep.seqs.qza](https://osf.io/2sh7n/))  
> `$BIGDBTAX` refers to the QIIME-formatted taxonomy file for the broad COI database ([bigCOI.derep.tax.qza](https://osf.io/aqvcj/))  
> `$READS` refers to the original [fasta](https://github.com/devonorourke/nhguano/blob/master/data/qiime_qza/repSeqs/all.raw_repSeqs.qza) `all.raw_repSeqs.qza` artifact  
> `$BIGCLFYR` refers to the QIIME-formatted Naive Bayes classifier file [nbClassifer_hostDBonly.qza](https://osf.io/vj6xn/)

```
## classify ASVs with host database
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$READS" --o-classification tmp.raw_hostDB_VStax.qza \
  --i-reference-reads "$HOSTDBSEQ" --i-reference-taxonomy "$HOSTDBTAX" \
  --p-maxaccepts 1000 --p-perc-identity 0.97 --p-query-cov 0.89 --p-strand both --p-threads 24

## classify ASVs with broad COI database
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$READS" --o-classification tmp.raw_bigDB_VStax.qza
  --i-reference-reads "$BIGDBSEQ" --i-reference-taxonomy "$BIGDBTAX" \
  --p-maxaccepts 1000 --p-perc-identity 0.97 --p-query-cov 0.89 --p-strand both --p-threads 24 \

## train the Naive Bayes classifier using the broad COI database
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$BIGDBSEQ" --i-reference-taxonomy "$BIGDBTAX" \
  --o-classifier nbClassifer_hostDBonly.qza

## classify ASVs with Naive Bayes classifier     
qiime feature-classifier classify-sklearn \
  --i-reads "$READS" --i-classifier "$BIGCLFYR" \
  --p-n-jobs 1 --p-reads-per-batch 2000 \
  --o-classification tmp.raw_bigDB_NBtax.qza
```

Each classification file output was exported into a `.tsv` format:
```
qiime tools export --input-path {some.qza} --output-path tmp && mv ./tmp/* . && mv taxonomy.tsv {some.tsv}
```

The resulting [tmp.raw_bigDB_NBtax.tsv](https://github.com/devonorourke/nhguano/blob/master/data/tax/tmp.raw_bigDB_NBtax.tsv), [tmp.raw_bigDB_VStax.tsv](https://github.com/devonorourke/nhguano/blob/master/data/tax/tmp.raw_bigDB_VStax.tsv), [tmp.raw_hostDB_VStax.tsv](https://github.com/devonorourke/nhguano/blob/master/data/tax/tmp.raw_hostDB_VStax.tsv) files were then available for analysis in the [decontam R script](https://github.com/devonorourke/nhguano/blob/master/scripts/r_scripts/decontam_efforts.R). The equivalent `.qza` files are available in [this directory](https://github.com/devonorourke/nhguano/data/qiime_qza/taxonomy).

We found substantial evidence that indicated that the only New Hampshire bat sampled in this study is _Myotis lucifugus_. Both the Naive Bayes and VSEARCH classifiers identified a common set of 27 ASVs assigned to several bat species. However, these other bats were present in just 1 or 2 samples generating just 22 reads total among the entire 9 libraries. The little brown bat, on the other hand, was detected in 595 samples, and generated over 1.6 million sequences. No other expected bat species was detected in our study, confirming that these guano samples are likely exclusively from little brown bats. Interestingly, the host DB method assigned just 3 unique ASVs totaling just 10 reads, and perhaps is an indication that the _lucifugus_ reference in that database is not similar enough to the NH species types we have in the broader COI reference. Nevertheless, these analyses confirm that our New Hampshire guano is highly likely to have originated from a single species: _Myotis lucifugus_.

All bat-associated from either database were removed from the `nocontrol_nomockASV_table.qza`:
> The `$BATASV` file refers to the [bat-associated list of ASVs](https://github.com/devonorourke/nhguano/data/host/batASVs.txt) generated at the end of the R script
> `$TABLE` refers to [sampleOnly_table.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_table.qza)
> `$NHMETA` refers to the [qiime_NHbat_meta.tsv](https://github.com/devonorourke/nhguano/data/metadata/qiime_NHbat_meta.tsv) metadata file

```
qiime feature-table filter-features \
  --i-table $TABLE --o-filtered-table tmpfilt4_table.qza \
  --m-metadata-file "$BATASV" --p-exclude-ids

qiime feature-table filter-samples \
  --i-table tmpfilt4_table.qza --o-filtered-table sampleOnly_nobatASV_table.qza \
  --m-metadata-file "$NHMETA"
```

Our final step in a filtering analysis requires making decisions on what ASVs should be retained based upon taxonomic information.

## Filtering out non-arthropod sequences

- Filtered just the ASVs that were in Naive Bayes with at least Family-name, but missing from VSEARCH in R script.
- Subset the original full fasta file:
```
seqkit grep raw_repSeqs.fasta --pattern-file vsearch_missingFaminfo_asvs.txt -w 0 > vsearch_missingFaminfo.fasta
```
- ran standalone vsearch to generate the %id values
```
READS=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/vsearch_missingFaminfo.fasta
REFS=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/bigDB/bigCOI.derep.fasta.gz
vsearch --usearch_global $READS --db $REFS \
--id 0.8 --query_cov 0.89 --strand both --maxaccepts 100 --threads 24 --blast6out vsearch_missingFam_vsearchOut.tsv

cat vsearch_missingFam_vsearchOut.tsv | cut -f 1,2,3,4,7,8 | gzip > vsearch_missingFam_vsearchOut.tsv.gz
```

Can sort the output to keep:
1. Only values with %qcov > 0.89
2. Only the top hit (regardless if there is a tie)

AACATTATATTTTATTTTTGGAATTTGAGCAGGTATAGTAGGAACTTCTTTAAGATTATTAATTCGAGCAGAATTAGGAAATCCTGGATCTTTAATTGGTGATGATCAAATTTATAATACTATTGTAACAGCCCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGG


TATTCTTTATTTTTTATTTGCCATCTGAGCAGGAATAATTGGATCATCCATAAGTATAATTATTCGACTAGAATTAGGATCATGTAATTCTTTAATTAATAATGATATAATTTATAATATTCTAGTAACAAGACACGGTTTTATTATAATTTTTTTTATAATTATACCTATTATAATCGGG
