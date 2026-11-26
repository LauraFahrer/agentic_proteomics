# Evaluation

## Using Identical or Similar Tools

## Using Self-Written Code

## Conclusion

1. Subtype classification:

| Analysis | Result 0006 | Result 0007 | Verdict |
| --- | --- | --- | ---|
| NMF Clustering | - (partially) reproduced<br>- AGS cellline in wrong cluster<br>- PP NMF ARI: 0.73, the other 1.0<br>- NMF parameter for most stable decomposition NOT reproduced | | |



2. Functional analysis:

| Analysis | Result 0006 | Result 0007 | Verdict |
| --- | --- | --- | ---|
| GO | - distribution on cellular components mostly reproduced |
| GSEA | - same processes identified<br>- significance depend on variant of permutation method (undefined in paper)<br>- weaker significance for PP YP | | |
| WGCNA | - parameters not reported in paper<br>- strong deviation of identified protein modules<br>- strongest modules associated with subtypes | | |


3. Kinase-activity analsis:
| PTM-SEA | - CDK4- and CDK6-signatures in EMT significantly higher<br>- confirmes paper findings | | |
| Kinome clustering/differential kinases | - confirms EMT/Metabloite classification | | |
| Outlier kinases | - mostly reproduced<br>- e.g. ERBB2 in N87 und FGFR2 in KatoIII | | |