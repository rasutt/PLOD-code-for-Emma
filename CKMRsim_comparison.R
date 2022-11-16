# Code showing half-sibling pair vs unrelated pair PLODs from CKMRsim package
# equal to mine up to rounding on bi-allelic srwqc dataset.

# Load CKMRsim package
library(CKMRsim)

# Can recompute, or skip and load results at bottom of script

# Set dataset name
dataset_name <- "srwqc"

# Load saved genind object
load(paste0("datasets/", dataset_name, "_genind.rdata"))

# Store the genind object in known variable
dataset.gi <- srwqc.gi
rm(srwqc.gi)

# Read and summarise data for finding PLODs
source("code/prepare_data.R")

# Find half-sibling versus unrelated pair PLODs and save them in the analysis
# folder. ~14secs
# source("code/find_plods.R")
load(paste0("analyses/", dataset_name, "/plods_df_and_hsp_up_plods.rdata"))

# Make locus names
loci_ales = names(ale_freqs)
loci = substr(loci_ales, 1, nchar(loci_ales) - 3)

# Combine data required for CKMRsim allele frequencies input dataset
ckmr_sim_AFs = data.frame(
  Chrom = "Unk", 
  Pos = rep(1:(length(ale_freqs) / 2), each = 2), 
  Locus = paste0("L", gsub(":", ".", loci, fixed = T)), 
  Allele = as.character(0:1),
  Freq = ale_freqs, 
  AlleIdx = NA, 
  LocIdx = NA, 
  row.names = 1:length(ale_freqs)
)
ckmr_sim_AFs = reindex_markers(ckmr_sim_AFs)

# Make CKMR object to interface with the package functions. ~6.7mins (4.3mins
# when more memory free I think) vs < 1sec for mine, but incorporates true and
# assumed genotype error models for simulation and likelihoods.  Setting epsilon
# = 0 for assumed error model to compare with my code.
s = Sys.time()
ckmr_obj = create_ckmr(
  D = ckmr_sim_AFs, kappa_matrix = kappas[c("HS", "U"), ],
  ge_mod_assumed_pars_list = list(epsilon = 0)
)
print(Sys.time() - s)

# Get original genotypes with numbers of alleles of each type at each locus
data_gts = tab(dataset.gi)

# Make matrix for new form with columns for type of each allele at each locus
ale_cpy_gts = array(dim = dim(data_gts))
colnames(ale_cpy_gts) = paste0("L", loci, c("_1", "_2"))

# Fill new matrix depending on original
ale_0_inds = seq(1, ncol(data_gts), 2)
ale_1_inds = seq(2, ncol(data_gts), 2)
ale_cpy_gts[, ale_0_inds][data_gts[, ale_0_inds] > 0] = "0"
ale_cpy_gts[, ale_0_inds][data_gts[, ale_0_inds] == 0] = "1"
ale_cpy_gts[, ale_1_inds][data_gts[, ale_1_inds] > 0] = "1"
ale_cpy_gts[, ale_1_inds][data_gts[, ale_1_inds] == 0] = "0"

# Combine sample names as first column
ale_cpy_gts_nm = data.frame(cbind(rownames(data_gts), ale_cpy_gts))
colnames(ale_cpy_gts_nm)[1] = "Indiv"

# Make long form genotype dataset with columns for sample name, locus, gene
# copy, and allele
library(tidyr)
long_gts = pivot_longer(
  data = ale_cpy_gts_nm,
  cols = 2:ncol(ale_cpy_gts_nm),
  names_to = c("Locus", "gene_copy"),
  names_sep = "_",
  values_to = "Allele"
)

# Compute PLODs. Does some data transformation, then lookup and
# summation over pairs and loci using C/C++. 9.3secs vs 14secs for mine.
s = Sys.time()
hsp_up_plods_ckmrsim = 
  pairwise_kin_logl_ratios(long_gts, long_gts, ckmr_obj, "HS", "U")
print(Sys.time() - s)

# Combine results to compare
comp_df1 = merge(
  hsp_up_plods_ckmrsim, plods_df, 
  by.x = c("D2_indiv", "D1_indiv"), by.y = c("samp_name_1", "samp_name_2")
)[-5:-6]
comp_df2 = merge(
  hsp_up_plods_ckmrsim, plods_df, 
  by.x = c("D1_indiv", "D2_indiv"), by.y = c("samp_name_1", "samp_name_2")
)[-5:-6]
comp_df = rbind(comp_df1, comp_df2)

# Normalise log-likelihood ratios from CKMRsim with number of shared loci and
# compare with my plods
comp_df$CKMRsimPLOD = comp_df$logl_ratio / comp_df$num_loc

# Save results
save(ale_freqs, ckmr_obj, long_gts, comp_df, file = "CKMRsim_comparison.Rdata")
# load(file = "CKMRsim_comparison.Rdata")

# Compare plods from ckmrsim and my code, same up to rounding error
head(comp_df)
nrow(comp_df) == nrow(plods_df)
all(abs(comp_df$CKMRsimPLOD - comp_df$PLOD) < 1e-10)


