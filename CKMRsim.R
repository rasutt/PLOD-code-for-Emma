# Code to try out CKMRsim package

# Load CKMRsim package
library(CKMRsim)

# Try our bi-allelic dataset srwqc for now.  I've found the allele frequencies
# using my own code.

# Make locus names
str(ale_freqs)
any(is.na(ale_freqs))
loci_ales = names(ale_freqs)
loci = substr(loci_ales, 1, nchar(loci_ales) - 3)
tail(loci)

# Combine data required for allele frequencies input dataset
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
head(ckmr_sim_AFs)
ckmr_sim_AFs = reindex_markers(ckmr_sim_AFs)
head(ckmr_sim_AFs)

# Make a ckmr object to interface with the package functions.  4.3mins for the
# srwqc dataset, which doesn't include finding PLODs.  Mine takes 14secs. Theirs
# does multi-core on non-windows, and incorporates true and assumed genotype
# error models, for simulation and likelihood calculation.  Using epsilon = 0
# for the assumed error model resolves the difference between our results.
s = Sys.time()
ckmr_obj = create_ckmr(
  D = ckmr_sim_AFs, kappa_matrix = kappas[c("HS", "U"), ],
  ge_mod_assumed_pars_list = list(epsilon = 0)
)
print(Sys.time() - s)
summary(ckmr_obj)

# Now I need my genotypes in the dumb long format required for this package

# Get original genotypes with numbers of alleles of each type at each locus
data_gts = tab(dataset.gi)
str(data_gts)
data_gts[1:5, 1:5]

# Make matrix for new form with columns for type of each allele at each locus
ale_cpy_gts = array(dim = dim(data_gts))
colnames(ale_cpy_gts) = paste0("L", loci, c("_1", "_2"))
ale_cpy_gts[1:5, 1:5]

# Fill new matrix depending on original
ale_0_inds = seq(1, ncol(data_gts), 2)
ale_1_inds = seq(2, ncol(data_gts), 2)
ale_cpy_gts[, ale_0_inds][data_gts[, ale_0_inds] > 0] = "0"
ale_cpy_gts[, ale_0_inds][data_gts[, ale_0_inds] == 0] = "1"
ale_cpy_gts[, ale_1_inds][data_gts[, ale_1_inds] > 0] = "1"
ale_cpy_gts[, ale_1_inds][data_gts[, ale_1_inds] == 0] = "0"
data_gts[1:50, 1:4]
ale_cpy_gts[1:50, 1:4]

# Combine sample names as first column
ale_cpy_gts_nm = data.frame(cbind(rownames(data_gts), ale_cpy_gts))
colnames(ale_cpy_gts_nm)[1] = "Indiv"
ale_cpy_gts_nm[1:4, 1:7]

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
head(long_gts)
sum(is.na(long_gts$Allele))

# Save data to speed up later attempts
# save(ale_freqs, ckmr_obj, long_gts, file = "CKMRsim_data.Rdata")
# load(file = "CKMRsim_data.Rdata")

# Compute PLODs, this actually just does some data transformation, and then
# lookup and summation (using C/C++) of plods over pairs and loci
hsp_up_plods_ckmrsim = 
  pairwise_kin_logl_ratios(long_gts, long_gts, ckmr_obj, "HS", "U")
head(hsp_up_plods_ckmrsim)

# Check same number of plods
nrow(plods_df)
nrow(hsp_up_plods_ckmrsim)

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
nrow(comp_df)

# Normalise log-likelihood ratios from CKMRsim with number of shared loci and
# compare with my plods
comp_df$logl_ratio = comp_df$logl_ratio / comp_df$num_loc
head(comp_df)
all(abs(comp_df$logl_ratio - comp_df$PLOD) < 1e-10)

# The results are the same up to rounding error when the parameter epsilon in
# the assumed genotype error model is set to zero.

