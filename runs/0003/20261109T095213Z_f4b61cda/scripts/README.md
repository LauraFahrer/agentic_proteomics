# Analysis scripts — PXD048538 run 0003 (signaling / PTM analysis)

NOTE: the quantms release e60b9e9aa4c8 has NO PTM layer. The phosphoproteome (PP),
tyrosine (YP) and whole-proteome (WP) matrices are from the publication Supplementary
Table S3 (MOL2-18-1437-s005.xls), fetched by 00_download.sh.

    export PXD048538_DATA=/path/to/data
    bash   00_download.sh "$PXD048538_DATA"   # data, Table S3 supplement, UniProt FASTA
    pip install phosx                          # Kinase Library PSSMs for kinase inference
    python 01_phospho_build.py                 # PP/YP line-mean + protein-normalised PP (fig 01)
    python 02_diff_phospho.py                  # differential phospho + seqrnk files (fig 02)
    # kinase activity: for each contrast and layer run phosX on its seqrnk, e.g.
    #   phosx -s -n 200 -c 4 -o results/tables/kinase/EMT_vs_Metabolism_ST_kinases.tsv \
    #         results/seqrnk/EMT_vs_Metabolism_ST.seqrnk
    #   phosx -t -n 200 -c 4 -o results/tables/kinase/EMT_vs_Metabolism_Y_kinases.tsv \
    #         results/seqrnk/EMT_vs_Metabolism_Y.seqrnk   (and the EBV_vs_nonEBV pair)
    python 03_kinase_activity.py               # parse phosX, kinase figure (fig 03)
    python 04_hypotheses.py                    # signaling-module figure + summary (fig 04)

- phospho_common.py : PP/YP/WP loaders, 10-mer window builder, Excel gene->date fix.
- common.py, stats.py : shared annotation/palette and moderated-t / GSEA helpers.

phosX is the slow step (~15 min per Ser/Thr contrast). Ser/Thr: Johnson 2023 PSSMs;
Tyr: Yaron-Barir 2024 PSSMs. Kinases are ranked by the continuous KS statistic because
the -log10(p) Activity Score saturates at the permutation floor.
