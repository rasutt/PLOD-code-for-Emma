# Read and summarise data for finding PLODs

# Find numbers of ales per locus in genind genotypes
ns_ales <- nAll(dataset.gi)

# Find number of loci
n_loci <- length(ns_ales)

# Get sample names
samp_name <- rownames(tab(dataset.gi))

# Find number of samples
n_samps <- length(samp_name)

# Find all pairs of sample indices, and number of pairs of samples
samp_prs_inds <- combn(n_samps, 2)
samp_inds_1 <- samp_prs_inds[1, ]
samp_inds_2 <- samp_prs_inds[2, ]
n_pairs <- choose(n_samps, 2)

# Find allele frequencies over all samples
ale_freqs <- colMeans(tab(dataset.gi), na.rm = T) / 2

# Find numbers of loci with each number of alleles. Table converts things to
# factors first and works on levels, so won't have zero counts
ns_ales_tab_full <- table(ns_ales)

# Stop if only one allele at every locus
if (all(names(ns_ales_tab_full) == "1")) stop("All loci monomorphic.")

# Remove monomorphic loci (just one allele), as they do not add any information
ns_ales_tab <- ns_ales_tab_full[names(ns_ales_tab_full) != "1"]

# Find numbers of alleles at any locus
ns_ales_any = as.numeric(names(ns_ales_tab))

# Find number of numbers of alleles at any locus
n_ns_ales <- length(ns_ales_any)

# Make lists for PLOD arrays and keys indices
plod_lst <- keys_inds_lst <- vector("list", n_ns_ales)

# Set scalars for expected values to zero
ev_sp <- ev_pop <- ev_up <- 0

# Loop over numbers of alleles at any locus
for (i in seq_along(ns_ales_any)) {
  n_ales = ns_ales_any[i]
  
  # Find number of possible genotypes
  n_gts <- n_ales * (n_ales + 1) / 2
  
  # Find possible genotypes 
  gts <- combn(n_ales, 2)
  gts <- cbind(gts, rbind(1:n_ales, 1:n_ales))
  ales_1 <- gts[1, ]
  ales_2 <- gts[2, ]
  
  # Find and store keys indices for looking up possible genotypes
  keys_inds <- 1:n_gts
  names(keys_inds) <- colSums(2^gts)
  keys_inds_lst[[i]] <- keys_inds
  
  # Find loci indices, and number
  loci_sz_bool <- ns_ales == n_ales
  loci_sz_inds <- as.character(locFac(dataset.gi)) %in% 
    names(ns_ales[loci_sz_bool])
  n_loci_sz <- sum(loci_sz_bool)
  loci_sz_inds_rep <- rep(1:n_loci_sz, each = n_gts)
  
  # Find allele frequencies, remember more than two alleles now
  ale_frqs_sz <- matrix(
    ale_freqs[loci_sz_inds], nrow = n_loci_sz, ncol = n_ales, byrow = T
  )
  
  # Find possible genotype probabilities
  gt_prbs_mat <- matrix(
    ale_frqs_sz[cbind(loci_sz_inds_rep, ales_1)] *
      ale_frqs_sz[cbind(loci_sz_inds_rep, ales_2)] *
      (1 + (ales_1 != ales_2)), 
    nrow = n_gts, ncol = n_loci_sz
  )
  
  # Create array for half-sibling versus unrelated pair plods and add possible
  # second genotype probabilities
  gt_2_prbs_mat <- array(
    rep(gt_prbs_mat, each = n_gts), dim = c(n_gts, n_gts, n_loci_sz)
  )
  hsp_up_plods_ary <- gt_2_prbs_mat
  
  # Find possible genopairs
  gts_1 <- gts[, rep(1:n_gts, n_gts)]
  gts_2 <- gts[, rep(1:n_gts, each = n_gts)]
  
  # Find allele equalities among possible genopairs
  cis_eqs <- gts_1 == gts_2
  cis_eqs_1_mat <- matrix(cis_eqs[1, ], nrow = n_gts)
  cis_eqs_2_mat <- matrix(cis_eqs[2, ], nrow = n_gts)
  
  trans_eqs_gts_2_htro <- gts_1 == gts_2[c(2, 1), ] & 
    rep(gts_2[1, ] != gts_2[2, ], each = 2)
  trans_eqs_gts_2_htro_1_mat <- matrix(trans_eqs_gts_2_htro[1, ], nrow = n_gts)
  trans_eqs_gts_2_htro_2_mat <- matrix(trans_eqs_gts_2_htro[2, ], nrow = n_gts)
  
  # Add conditional probabilities of possible second genotypes, given
  # parent-offspring with first, to PLODs, when alleles shared
  hsp_up_plods_ary[rep(cis_eqs_1_mat, n_loci_sz)] <- 
    hsp_up_plods_ary[rep(cis_eqs_1_mat, n_loci_sz)] + 0.5 *
    ale_frqs_sz[cbind(rep(1:n_loci_sz, each = sum(cis_eqs_1_mat)), 
                      gts_2[2, cis_eqs[1, ]])]
  hsp_up_plods_ary[rep(cis_eqs_2_mat, n_loci_sz)] <- 
    hsp_up_plods_ary[rep(cis_eqs_2_mat, n_loci_sz)] + 0.5 *
    ale_frqs_sz[cbind(rep(1:n_loci_sz, each = sum(cis_eqs_2_mat)), 
                      gts_2[1, cis_eqs[2, ]])]
  
  # Don't add twice when second genotype homozygous
  hsp_up_plods_ary[rep(trans_eqs_gts_2_htro_1_mat, n_loci_sz)] <- 
    hsp_up_plods_ary[rep(trans_eqs_gts_2_htro_1_mat, n_loci_sz)] + 0.5 *
    ale_frqs_sz[cbind(rep(1:n_loci_sz, each = sum(trans_eqs_gts_2_htro_1_mat)), 
                      gts_2[1, trans_eqs_gts_2_htro[1, ]])]
  hsp_up_plods_ary[rep(trans_eqs_gts_2_htro_2_mat, n_loci_sz)] <- 
    hsp_up_plods_ary[rep(trans_eqs_gts_2_htro_2_mat, n_loci_sz)] + 0.5 *
    ale_frqs_sz[cbind(rep(1:n_loci_sz, each = sum(trans_eqs_gts_2_htro_2_mat)), 
                      gts_2[2, trans_eqs_gts_2_htro[2, ]])]
  
  # Find conditional second genotype probabilities given parent-offspring with
  # first
  cnd_gt_2_prbs_pop_ary <- hsp_up_plods_ary - gt_2_prbs_mat
  
  # Take logs of HSP vs UP pseudo likelihood ratios but skip dividing by 2 for
  # now
  hsp_up_plods_ary <- log(hsp_up_plods_ary) - 
    rep(log(gt_prbs_mat), each = n_gts)

  # Add PLOD array to list
  plod_lst[[i]] <- hsp_up_plods_ary
  
  # Make array of first genotype probabilities
  gt_1_prbs_mat <- array(gt_prbs_mat[, rep(1:n_loci_sz, each = n_gts)], 
                         dim = c(n_gts, n_gts, n_loci_sz))
  
  # Find and add expected values
  ev_sp <- ev_sp + sum(hsp_up_plods_ary[cbind(
    rep(1:n_gts, n_loci_sz), rep(1:n_gts, n_loci_sz), 
    rep(1:n_loci_sz, each = n_gts))] * gt_prbs_mat)
  ev_pop <- ev_pop + sum(gt_1_prbs_mat * cnd_gt_2_prbs_pop_ary * 
                           hsp_up_plods_ary)
  ev_up <- ev_up + sum(gt_1_prbs_mat * gt_2_prbs_mat * hsp_up_plods_ary)
}

# Divide expected values by number of loci, divide by two, and combine in list
ev_sp <- ev_sp / n_loci - log(2)
ev_pop <- ev_pop / n_loci - log(2)
ev_up <- ev_up / n_loci - log(2)
ev_hsgpop <- (ev_pop + ev_up) / 2
ev_tcggpop <- (ev_pop + 3 * ev_up) / 4
ev_fcgtcgggpop <- (ev_pop + 7 * ev_up) / 8
evs <- c(ev_up, ev_fcgtcgggpop, ev_tcggpop, ev_hsgpop, ev_pop, ev_sp)
names(evs) <- c("Unrelated", "First cousin", "Avuncular", "Half-sibling",
                    "Parent-offspring", "Self")

# Display progress
cat("Found", n_samps, "samples, with", n_loci, "loci \n")
cat("Numbers of loci with each number of alleles: \n")
print(ns_ales_tab_full)
