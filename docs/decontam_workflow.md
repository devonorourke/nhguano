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
1. The most prevalent ASVs in our dataset are routinely identified as contaminants by the DNAplate method, but are typically not by the sequencing batch method. This makes me suspicious that there is not a sufficient statistical power to support the DNAplate-based batch filter. For example, it could be that a partiuclar ASV is marked as a contaminant in 3 of 40 DNA plates, but not so in the other 37. Nevertheless, we're seeing some ASVs in hundreds of samples, and these are likely _not_ contaminants: they're the kinds of sequence variants we expect in the bat diet (each ASV was classified and assigned taxonomic identity, explained in the next section).
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

# Diversity analyses using all data

We'll rarefy our data prior to calculating any community composition distances. QIIME 2 has an [alpha rarefaction function](https://docs.qiime2.org/2019.4/plugins/available/diversity/alpha-rarefaction/) whereby we can visualize how the observed richness of a sample (or other alpha diversith metrics) changes with sampling depth. Because this analysis is focused on looking at how negative control samples may or may not associate with a particular group (i.e. DNA extraction plate, or Sequencing batch) we care more about preserving as many NTC samples as possible, so we're going to focus our alpha diversity analysis only on those samples. Thus, we're first going to filter out all other non-NTC samples from the `tmp.raw_table.qza` initially imported in the decontam R script, then run the function on those NTC samples for a range of sampling depths. We'll then choose a single depth, rarefy _all_ samples at that depth, then use that new table in the R script to carry out the analysis.

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

The `tmp.neg_sumrytable.qzv` file (visualized in the table below) suggests that a sampling depth between 400-600 sequences will retain a balance betwen preserving the most number of sequence variants and the greatest number of samples. However this is only relative to the negative control samples; clearly, we'd want to have a higher number of sequences for the true samples.

| Sampling_depth | % ASVs remaining in dataset | # Samples |
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

What's

# QIIME 2 data filtering
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

## Next steps in analysis
The filtered dataset was used as input into the [diversity analyses](https://github.com/devonorourke/nhguano/blob/master/docs/decontam_workflow.md) pipeline, which included normalizing the data by rarefying, and estimating measures of alpha and beta diversity.
