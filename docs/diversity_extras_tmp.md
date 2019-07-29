
## for DEICODE, need to use non-rarefied table, but want to focus on specific samples:
NONRARETABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_arthOnly_table.qza
FILTLIST=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/alpha_study1names.txt
qiime feature-table filter-samples \
--i-table $NONRARETABLE \
--m-metadata-file $FILTLIST \
--o-filtered-table sampleOnly_select2016_nonrfyd_table.qza

## install deicode with conda
## see tutorial here: https://forum.qiime2.org/t/robust-aitchison-pca-beta-diversity-with-deicode/8333
conda install -c conda-forge deicode

## run deicode to generate pcoa biplot
DEICODETABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_select2016_nonrfyd_table.qza
qiime deicode rpca \
    --i-table $DEICODETABLE \
    --p-min-feature-count 2 \
    --p-min-sample-count 1000 \
    --o-biplot deicode_biplot.qza \
    --o-distance-matrix deicode_distmat.qza

DEICODETABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_select2016_nonrfyd_table.qza
qiime deicode rpca \
    --i-table $DEICODETABLE \
    --p-min-feature-count 10 \
    --p-min-sample-count 1000 \
    --o-biplot deicode_biplot_min10feat.qza \
    --o-distance-matrix deicode_distmat_min10feat.qza


# consider using abundances and not normalizing
for alpha diversity, could look at richness here: https://forum.qiime2.org/t/q2-breakaway-community-tutorial/5756

For beta diversity, Nick mentioned using a DEseq wald test. See: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#variations-to-the-standard-workflow

There's also a Corncob tutorial here: https://github.com/bryandmartin/corncob/blob/master/vignettes/corncob-intro.Rmd

Corncob paper here: https://arxiv.org/pdf/1902.02776.pdf
