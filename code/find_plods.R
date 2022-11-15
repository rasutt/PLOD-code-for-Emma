# Find and save half-sibling versus unrelated pair PLODs

# Set size of batches of loci to keep memory usage < 1GB
# n_pairs * btch_sz * 8 < 1Gb = 1e9 ~ 2^30 => btch_sz ~< 2^27 / n_pairs
btch_sz <- 2^(26 - ceiling(log(n_pairs, 2)))

# Find numbers of batches
ns_btchs <- ceiling(ns_ales_tab / btch_sz)

# Create directory for analysis
dir.create(path = paste0("analyses/", dataset_name))

# Make vectors for PLODs and numbers of shared markers
hsp_up_plods <- ns_shrd_mrkrs <- numeric(n_pairs)

# Set batch counter to zero
btch_cnt <- 0

# Set timer
s_time <- proc.time()[3]

# Display progress
cat("Finding", n_pairs, "PLODs, in", sum(ns_btchs), "batches, over", n_loci,
    "loci \n")
cat("Batch: ")

# Loop over numbers of alleles at any locus
for (i in seq_along(ns_ales_any)) {
  n_ales = ns_ales_any[i]
  
  # Find loci indices, and number
  loci_sz_inds <- as.character(locFac(dataset.gi)) %in% 
    names(ns_ales[ns_ales == n_ales])
  n_loci_sz <- ns_ales_tab[as.character(n_ales)]
  
  # Set timer for this number of alleles
  s_time_sz <- proc.time()[3]
  
  # Get observed genotypes and transform into keys indices for looking up
  # probabilities
  gt_keys_mat <- apply(array(tab(dataset.gi)[, loci_sz_inds], 
                             dim = c(n_samps, n_ales, n_loci_sz)) * 
                         2^rep(1:n_ales, each = n_samps), c(1, 3), sum)
  gt_keys_inds_mat <- matrix(keys_inds_lst[[i]][as.character(gt_keys_mat)], 
                             nrow = n_samps)
  
  # Loop over batches of loci
  for(btch_ind in 1:ns_btchs[i]) {
    # Increment batch counter
    btch_cnt <- btch_cnt + 1
    
    # Display progress
    cat(btch_cnt, "")
    
    # Find locus indices
    loci_inds <- 
      ((btch_ind - 1) * btch_sz + 1):min(btch_ind * btch_sz, n_loci_sz)
    
    # Lookup plods
    hsp_up_plods_obs_mat <- matrix(plod_lst[[i]][
      cbind(as.vector(gt_keys_inds_mat[samp_inds_1, loci_inds]), 
            as.vector(gt_keys_inds_mat[samp_inds_2, loci_inds]), 
            rep(loci_inds, each = n_pairs))], nrow = n_pairs)
    
    # Find and add numbers of shared markers
    ns_shrd_mrkrs <- ns_shrd_mrkrs + 
      (length(loci_inds) - rowSums(is.na(hsp_up_plods_obs_mat)))
    
    # Find and add plods over shared markers
    hsp_up_plods <- hsp_up_plods + rowSums(hsp_up_plods_obs_mat, na.rm = T)
  }
}

# Divide by numbers of shared markers. Include division by two skipped for hsp
# probs earlier. P(gts|HSP) = (P(gts|UP) + P(gts|POP)) / 2 at each locus, so
# subtract log(2) * n_loci and divide by n_loci.
hsp_up_plods <- hsp_up_plods / ns_shrd_mrkrs - log(2)

# Combine PLODs with sample names and populations
plods_df <- data.frame(
  samp_name_1 = samp_name[samp_inds_1], 
  samp_name_2 = samp_name[samp_inds_2],
  ppln_1 = ppln[samp_inds_1], 
  ppln_2 = ppln[samp_inds_2],
  PLOD = hsp_up_plods,
  ns_shrd_mrkrs = ns_shrd_mrkrs
)

# Add sample names to plods
names(hsp_up_plods) <- 
  paste0(samp_name[samp_inds_1], ".", samp_name[samp_inds_2])

# Save PLODs
save(evs, plods_df, hsp_up_plods, 
     file = paste0("analyses/", dataset_name, 
                   "/plods_df_and_hsp_up_plods.rdata"))

# Display progress
cat("Done \n")
cat("Time taken:", proc.time()[3] - s_time, "seconds \n")
