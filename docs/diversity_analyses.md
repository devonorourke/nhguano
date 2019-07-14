
# Overview
We initially processed all samples following steps described in the [sequence_processing](https://github.com/devonorourke/nhguano/blob/master/docs/sequence_processing.md) document, and performed a detailed analysis of potential contamination among negative control and positive control samples as outlined in the [decontam_workflow](https://github.com/devonorourke/nhguano/blob/master/docs/decontam_workflow.md) file. The final ASV table had been filtered to remove both bat (host) DNA and non-arthropod sequence variants.

# Determining sampling depth for normalization
Because the previous filtering regime did not remove samples with low read abundances (any samples with >= 1 read were retained) we're going to first determine how the sampling depth will effect the number of samples retained, the overall richness of the dataset, and the per-sample richness.
> `$TABLE` refers to the final taxonomy-filtered [ASV table](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_arthOnly_table.qza) produced at the conclusion of the [decontam_workflow](https://github.com/devonorourke/nhguano/blob/master/docs/decontam_workflow.md) document
> `$NHMETA` refers to the [qiime_NHbat_meta.tsv](https://github.com/devonorourke/nhguano/data/metadata/qiime_NHbat_meta.tsv) metadata file

```
## generating summary table
qiime feature-table summarize \
--i-table $TABLE --o-visualization sampleOnly.summaryTable.qzv --m-sample-metadata-file $NHMETA

## boostrap estimates of richness per sample across a range of sampling depths
qiime diversity alpha-rarefaction \
  --i-table $TABLE --o-visualization sampleOnly.alphaRareViz.qzv --m-metadata-file $NHMETA \
  --p-metrics observed_otus --p-min-depth 500 --p-max-depth 5000 --p-iterations 15
```

We used the [sampleOnly.summaryTable.qzv](https://github.com/devonorourke/nhguano/data/qiime_qzv/table_sumry/sampleOnly.summaryTable.qzv) summary visualization to determine how the sampling depth would influence the number of samples that would be retained for different variables like sampling sites and weeks. This table demonstrates the tradeoff we expected: because we pooled hundreds of samples on a single MiSeq run the average read depth per sample is low (per sample sequencing counts median = 400), but there are hundreds of samples with thousands of reads. Because we are not particularly concerned with rare variants or low abundance taxa and instead are more interested in retaining as many samples per group as possible, we sought to find a sampling depth that preserved a balance between retaining as many samples as possible at a point where observed ASV richness began to plateau on the alpha rarefaction curve. The [sampleOnly.alphaRareViz.qzv](https://github.com/devonorourke/nhguano/data/qiime_qzv/alpha_viz/sampleOnly.alphaRareViz.qzv) visualization suggests that a relatively low sampling depth between 1000-2000 reads will retain most of the observed variation among all but a few sites (one site, EPS (Epsom, NH) had far more diversity than other samples); importantly this range of sampling depth will preserve many of our samples. Samples were then rarefied to our selected sampling depth of 1000 reads.

```
qiime feature-table rarefy \
  --i-table $TABLE \
  --p-sampling-depth 1000 \
  --o-rarefied-table sampleOnly_rfyd_table.qza
```

The rarefied table ([sampleOnly_rfyd_table.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_rfyd_table.qza)) is used next to calculate estimates of diversity.

# Alpha diversity estimates

We explored three estimates of within sample ASV diversity: richness (observed ASVs), Shannon's entropy, and Faith's phylogenetic diversity.
`$TREE` refers to the rooted tree file ([raw.ASVtree_rooted.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/trees)) created previously in the `decontam workflow` document
`$TABLE` refers to the rarefied ASV table ([sampleOnly_rfyd_table.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_rfyd_table.qza))

```
qiime diversity alpha-phylogenetic --i-table "$TABLE" --i-phylogeny "$TREE" --p-metric faith_pd --o-alpha-diversity alpha.vals_fp.qza
qiime diversity alpha --i-table "$TABLE" --p-metric observed_otus --o-alpha-diversity alpha.vals_ob.qza
qiime diversity alpha --i-table "$TABLE" --p-metric shannon --o-alpha-diversity alpha.vals_sh.qza
```

All `alpha.vals*qza` artifacts are available at [this directory](https://github.com/devonorourke/nhguano/data/qiime_qza/alpha). We partitioned our analyses to focus on particular 2016 sites that contained the greatest proportion of samples across that years sampling dates (April to October). These `.qza` files input to the [NH_diversity.R](https://github.com/devonorourke/nhguano/scripts/r_scripts/alphadiv.R) script produce the tables and figures presented in this manuscript.


## Beta diversity estimates

We compared community composition among two sets of samples: select 2016 samples and particular sites repeatedly sampled in 2015 and 2016. We assessed differences in ASV composition using four distance metrics:
- Dice-Sorensen (unweighted abundance, unweighted phylogenetic)  
- Bray-Curtis (weighted abundance, unweighted phylogenetic)  
- Unweighted Unifrac (unweighted abundance, weighted phylogenetic)  
- Weighted Unifrac (weighted abundance, weighted phylogenetic)  

The following code was executed to generate the distance matrices. We then applied a Principal Correspondence Analysis for each distance matrix:  
```
## distance estimates
qiime diversity beta-phylogenetic --i-table "$TABLE" --i-phylogeny "$TREE" --p-metric unweighted_unifrac --o-distance-matrix uu_dist.qza
qiime diversity beta-phylogenetic --i-table "$TABLE" --i-phylogeny "$TREE" --p-metric weighted_unifrac --o-distance-matrix wu_dist.qza
qiime diversity beta --i-table "$TABLE" --p-metric dice --o-distance-matrix ds_dist.qza
qiime diversity beta --i-table "$TABLE" --p-metric braycurtis --o-distance-matrix bc_dist.qza

## pcoa
qiime diversity pcoa --i-distance-matrix uu_dist.qza --o-pcoa uu_pcoa.qza
qiime diversity pcoa --i-distance-matrix wu_dist.qza --o-pcoa wu_pcoa.qza
qiime diversity pcoa --i-distance-matrix ds_dist.qza --o-pcoa ds_pcoa.qza
qiime diversity pcoa --i-distance-matrix bc_dist.qza --o-pcoa bc_pcoa.qza
```

The `$TABLE` used in the distance estimate was first filtered according to each sub study. In the case of the 2016 dataset, we filtered the rarefied [sampleOnly_rfyd_table.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_rfyd_table.qza) table to contain only the samples selected in that particular investigation.
> `$STUDY1META` refers to the select list of samples filtered in the [NH_diversity.R](https://github.com/devonorourke/nhguano/scripts/r_scripts/alphadiv.R) script that pertain to nine sites from 2016: the [alpha_study1names.txt](https://github.com/devonorourke/nhguano/data/metadata/alpha_study1names.txt) file
> `$RARETABLE` refers to the original rarefied table: [sampleOnly_rfyd_table.qza](https://github.com/devonorourke/nhguano/data/qiime_qza/ASVtable/sampleOnly_rfyd_table.qza)

```
qiime feature-table filter-samples \
--i-table "$RARETABLE" --m-metadata-file "$STUDY1META" --o-filtered-table sampleOnly_select2016_rfyd_table.qza
```

Thus, the above beta distance and PCoA calculations are representing a generic function: the `$TABLE` would be altered depending upon which study was being analyzed. All distance metrics are available at [this directory](https://github.com/devonorourke/nhguano/data/qiime_qza/distmat/select2016), while all PCoA artifacts are [available here](https://github.com/devonorourke/nhguano/data/qiime_qza/pcoa/select2016).  





Last section on machine learning??

# consider using abundances and not normalizing
for alpha diversity, could look at richness here: https://forum.qiime2.org/t/q2-breakaway-community-tutorial/5756

For beta diversity, Nick mentioned using a DEseq wald test. See: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variations-to-the-standard-workflow
There's also a Corncob tutorial here: https://github.com/bryandmartin/corncob/blob/master/vignettes/corncob-intro.Rmd
Corncob paper here: https://arxiv.org/pdf/1902.02776.pdf
