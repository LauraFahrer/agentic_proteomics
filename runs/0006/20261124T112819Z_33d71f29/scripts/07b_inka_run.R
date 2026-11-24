#!/usr/bin/env Rscript
# Step 07b - run the original InKA code (digitalbiologylab/inka-intensity @782c1ab,
# inka_v1.6.0.R / InKAscore 1.5.0) on the YP inputs built by 07a. Unmodified script; it is
# sourced inside the work directory with its bundled resources (PSP Jul-2016, NetworKIN 2014,
# HGNC 2016, phomics activation-loop table) symlinked next to it.
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "common.R"))
W <- file.path(DC, "inka", "work")
res_dir <- file.path(SRC, "inka-intensity", "inka_20220330")
for (f in c("UniProt_data_08Jun2016.Rdata", "PSP_human_03Jul2016.Rdata", "PSP_KSR_human_03Jul2016.Rdata",
            "Kinases_HGNC2017.01.20_v2.Rdata", "complete-HGNC-21Apr2016.csv", "phomics_global_output.txt",
            "TOP_NWK_nwk_proteome_unprot_ref_2014.Rdata", "inka_v1.6.0.R")) {
  if (!file.exists(file.path(W, f))) file.symlink(file.path(res_dir, f), file.path(W, f))
}
setwd(W)
TXT_DIR <- "./"
commandArgs <- function(trailingOnly = FALSE) if (trailingOnly) c("--no.network") else base::commandArgs(FALSE)
source("inka_v1.6.0.R")
