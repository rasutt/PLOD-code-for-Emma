# Copy of pairwise_kin_logl_ratios function from CKMRsim, with my own formatting
# and comments.

# Function to find kinpair log-likelihood ratios from genotype datasets, and
# object created from allele frequencies.  Optionally only keep results for
# keep_top pairs with top results.
function (
    D1, D2, CK, numer, denom = NULL, keep_top = NULL, 
    num_cores = parallel::detectCores()
) {
  # Make 'integer genotype matrices' from long-form genotype datasets and allele
  # frequencies
  M1 <- create_integer_genotype_matrix(D1, CK$orig_data)
  M2 <- create_integer_genotype_matrix(D2, CK$orig_data)
  
  # Get possible probabilities, and shared loci, for numerator kinship in ratio
  numer_flat <- flatten_ckmr(CK, numer)
  logl_flat <- numer_flat
  
  # Get possible probabilities, and shared loci, for denominator kinship, and
  # compute log-ratios
  if (!is.null(denom)) {
    denom_flat <- flatten_ckmr(CK, denom)
    logl_flat$probs <- log(numer_flat$probs/denom_flat$probs)
    rm(denom_flat)
  }
  else {
    logl_flat$probs <- log(numer_flat$probs)
  }
  rm(numer_flat)
  
  # Set number of cores for parallel computation in windows systems
  if (.Platform$OS.type == "windows") {
    num_cores <- 1
  }
  
  # Loop over genotypes in second dataset, finding actual plods with all
  # genotypes in first for each one.  Seems to call C/C++ code to lookup and sum
  # log-likelihood ratios at each locus.
  idx <- 1:nrow(M2)
  names(idx) <- idx
  comps <- parallel::mclapply(idx, mc.cores = num_cores, FUN = function(i) {
    tmp <- comp_ind_pairwise(
      S = M1, T = M2, t = i, values = logl_flat$probs, 
      nGenos = logl_flat$nGenos, Starts = logl_flat$base0_locus_starts
    )
    if (!is.null(keep_top)) {
      ret <- tmp[rev(top_index(tmp$value, min(keep_top, nrow(tmp)))), ]
    }
    else {
      ret <- tmp
    }
    ret
  }) %>% dplyr::bind_rows(.id = "D2_indiv") %>% tibble::as_tibble() %>% 
    dplyr::rename(D1_indiv = ind) %>% 
    dplyr::mutate(D2_indiv = as.integer(D2_indiv))
  
  # If two copies of same genotype dataset delete self-comparisons and
  # duplicated comparisons (could also check this at the start to avoid extra
  # computation for larger datasets)
  if (identical(D1, D2)) {
    message("D1 and D2 are identical: dropping self comparisons and keeping only
            first instance of each pair")
    comps <- comps %>% dplyr::filter(D2_indiv < D1_indiv)
  }
  
  # Format and return results
  comps %>% 
    dplyr::mutate(
      D2_indiv = rownames(M2)[D2_indiv], 
      D1_indiv = rownames(M1)[D1_indiv]
    ) %>% dplyr::rename(logl_ratio = value)
}
