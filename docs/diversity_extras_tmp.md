NONRARETABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_arthOnly_table.qza

TABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_rfyd_table.qza
FILTTABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_select2016_rfyd_table.qza
FILTLIST=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/alpha_study1names.txt
TREE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/trees/raw.ASVtree_rooted.qza
BETADIR=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/betas/select2016
qiime feature-table filter-samples \
--i-table $TABLE \
--m-metadata-file $FILTLIST \
--o-filtered-table sampleOnly_select2016_rfyd_table.qza
  ## sampleOnly_select2016_rfyd_table.qza has 529 samples; matches number of samples listed in FILTLIST...

qiime diversity beta-phylogenetic --i-table "$FILTTABLE" --i-phylogeny "$TREE" --p-metric unweighted_unifrac --o-distance-matrix s16_dist_uu.qza
qiime diversity beta-phylogenetic --i-table "$FILTTABLE" --i-phylogeny "$TREE" --p-metric weighted_unifrac --o-distance-matrix s16_dist_wu.qza
qiime diversity beta --i-table "$FILTTABLE" --p-metric dice --o-distance-matrix s16_dist_ds.qza
qiime diversity beta --i-table "$FILTTABLE" --p-metric braycurtis --o-distance-matrix s16_dist_bc.qza

qiime diversity pcoa --i-distance-matrix "$BETADIR"/s16_dist_uu.qza --o-pcoa s16_pcoa_uu.qza
qiime diversity pcoa --i-distance-matrix "$BETADIR"/s16_dist_wu.qza --o-pcoa s16_pcoa_wu.qza
qiime diversity pcoa --i-distance-matrix "$BETADIR"/s16_dist_ds.qza --o-pcoa s16_pcoa_ds.qza
qiime diversity pcoa --i-distance-matrix "$BETADIR"/s16_dist_bc.qza --o-pcoa s16_pcoa_bc.qza

## for biplot, need a relative frequency table
qiime feature-table relative-frequency --i-table $FILTTABLE --o-relative-frequency-table sampleOnly_select2016_rfyd_relfreq_table.qza
RELFREQTABLE=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/reads/sampleOnly_select2016_rfyd_relfreq_table.qza
PCOADIR=/mnt/lustre/macmaneslab/devon/guano/paper3/qiime/select_libs/pcoa
qiime diversity pcoa-biplot --i-pcoa s16_pcoa_wu.qza --o-biplot s16_pcoabiplot_wu.qza --i-features $RELFREQTABLE

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



Phyllophaga  drakii, Phyllophaga  forsteri, Phyllophaga  tristis


Phyllophaga  anxia
Phyllophaga  drakii
Phyllophaga  forsteri
Phyllophaga  hirticula
Phyllophaga  marginalis
Phyllophaga  tristis
Phyllophaga  fervida
Phyllophaga  crenulata

Phyllophaga  fraterna
Phyllophaga  fusca
Phyllophaga  gracilis
Phyllophaga  longispina
Phyllophaga  occidentalis

Phyllophaga crassissima
Phyllophaga fervida
Phyllophaga foxii
......Phyllophaga hirsuta
Phyllophaga ilicis
