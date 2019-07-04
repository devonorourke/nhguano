# Overview

# Determining sampling depth for normalization
4. Figure out sampling depth with alpha rarefaction viz. Using default 10 iterations, but setting 1000 reads as minimum.

```
qiime diversity alpha-rarefaction \
  --p-metrics observed_otus --p-min-depth 1000 --p-max-depth 10000 --p-iterations 10 \
  --i-table study.raw_table.qza --o-visualization study.raw.alphaRareViz.qzv
```
