## Script requires a qiime2 installed
## I run this with a virtual Conda environment named 'rescript_2020.6'

## module load anaconda3
## conda activate rescript_2020.6

## import tree file for all datasets
## convert to QIIME-formatted .qza file
wget https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/trees/finalfiltd_rootedtree.nwk
qiime tools import --input-path finalfiltd_rootedtree.nwk --output-path finalfiltd_rooted_tree.qza --type 'Phylogeny[Rooted]'


ORIDIR=$(pwd)

for PREFIX in $(printf "win456\nwin567\n"); do
  ## create tmp directory to store all files
  mkdir -p "$PREFIX"_tmpdir
  cd "$PREFIX"_tmpdir
  ln -s ../finalfiltd_rooted_tree.qza .

  ## import metadata
  wget $(echo https://raw.githubusercontent.com/devonorourke/nhguano/master/data/metadata/"$PREFIX"meta_forQIIME.txt)

  ## import OTU table
    wget https://raw.githubusercontent.com/devonorourke/nhguano/master/data/text_tables/otu_tables/"$PREFIX"_OTUtable.tsv

  ## convert OTU table to biom format, then import for QIIME
    biom convert -i "$PREFIX"_OTUtable.tsv -o "$PREFIX"_OTUtable.biom --to-hdf5
    qiime tools import --input-path "$PREFIX"_OTUtable.biom --output-path "$PREFIX"_OTUtable.qza --type FeatureTable[Frequency]
  ## calculate distances using four different metrics
    qiime diversity beta --i-table "$PREFIX"_OTUtable.qza --p-metric 'dice' --o-distance-matrix "$PREFIX"_dist_dice.qza
    qiime diversity beta --i-table "$PREFIX"_OTUtable.qza --p-metric 'braycurtis' --o-distance-matrix "$PREFIX"_dist_bray.qza
    qiime diversity beta-phylogenetic --i-table "$PREFIX"_OTUtable.qza --i-phylogeny finalfiltd_rooted_tree.qza --p-metric 'unweighted_unifrac' --o-distance-matrix "$PREFIX"_dist_uuni.qza
    qiime diversity beta-phylogenetic --i-table "$PREFIX"_OTUtable.qza --i-phylogeny finalfiltd_rooted_tree.qza --p-metric 'weighted_unifrac' --o-distance-matrix "$PREFIX"_dist_wuni.qza

## run beta group significance for each distance; once for 'Window' and once for 'Site' variables
  ## export each beta group significance visualization, collect raw text data, and add
  #### for site factorgroup
    for METRIC in $(printf "dice\nbray\nuuni\nwuni\n"); do
      qiime diversity beta-group-significance \
        --i-distance-matrix "$PREFIX"_dist_"$METRIC".qza \
        --m-metadata-file "$PREFIX"meta_forQIIME.txt \
        --m-metadata-column 'Site' --o-visualization "$PREFIX"_betasig_site_"$METRIC".qzv;
      qiime tools export \
        --input-path "$PREFIX"_betasig_site_"$METRIC".qzv \
        --output-path "$PREFIX"_betasig_site_"$METRIC"_output
      cut -f 2-6 ./"$PREFIX"_betasig_site_"$METRIC"_output/raw_data.tsv | \
      awk -v var=$METRIC 'NR>1 {print $0"\t"var}' | \
      awk -v var=$PREFIX '{print $0"\t"var"\t""site"}' > "$PREFIX"_betasig_site_"$METRIC"_data.tsv
    #### for window factorgroup
      qiime diversity beta-group-significance \
        --i-distance-matrix "$PREFIX"_dist_"$METRIC".qza \
        --m-metadata-file "$PREFIX"meta_forQIIME.txt \
        --m-metadata-column 'Window' --o-visualization "$PREFIX"_betasig_window_"$METRIC".qzv;
      qiime tools export \
        --input-path "$PREFIX"_betasig_window_"$METRIC".qzv \
        --output-path "$PREFIX"_betasig_window_"$METRIC"_output
      cut -f 2-6 ./"$PREFIX"_betasig_window_"$METRIC"_output/raw_data.tsv | \
      awk -v var=$METRIC 'NR>1 {print $0"\t"var}' | \
      awk -v var=$PREFIX '{print $0"\t"var"\t""window"}' > "$PREFIX"_betasig_window_"$METRIC"_data.tsv
      
    done

    ## combine all outputs into single file for particular group (win456 or win567):
  printf 'SubjectID1\tSubjectID2\tGroup1\tGroup2\tDistance\tMetric\tWindowGroup\tFactorGroup\n' > "$PREFIX"_betasig_allMetrics_data.tsv

  cat "$PREFIX"_betasig_site_dice_data.tsv "$PREFIX"_betasig_site_bray_data.tsv \
  "$PREFIX"_betasig_site_uuni_data.tsv "$PREFIX"_betasig_site_wuni_data.tsv \
  "$PREFIX"_betasig_window_dice_data.tsv "$PREFIX"_betasig_window_bray_data.tsv \
  "$PREFIX"_betasig_window_uuni_data.tsv "$PREFIX"_betasig_window_wuni_data.tsv \
  >> "$PREFIX"_betasig_allMetrics_data.tsv

  cd $ORIDIR

done
