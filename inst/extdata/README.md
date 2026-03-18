This directory is intentionally lightweight in BrainEnrich 2.x.

BrainEnrich keeps a small offline fixture bundle here for tests, examples, and
vignette builds:

- `geneExp/desikan_r0.4.csv.bz2`
- `geneExp/desikan_r0.6.csv.bz2`
- `geneSets/SynGO.rds`
- `geneSets/GO_MF.rds`

The full curated expression matrices, gene-set collections, and manifest live
in the companion package `BrainEnrichData`. BrainEnrich resolves resources from
an installed `BrainEnrichData`, a sibling development checkout, these packaged
fixtures, or a cache download from the matching BrainEnrichData release.
