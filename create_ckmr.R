# Copy of create_ckmr function from CKMRsim, with my own formatting
# and comments.

# Function to create ckmr_class object, including possible genotype-pair
# probabilities at each locus given range of kinships, from dataset of allele
# frequencies, and kappa matrix for required kinships.  Also includes objects
# for simulation of new genotype-pairs given certain relationships, based on
# specified genotype error models.
function (
    D, kappa_matrix = kappas[c("MZ", "PO", "FS", "HS", "U"), ], 
    ge_mod_assumed = ge_model_TGIE, ge_mod_true = ge_model_TGIE, 
    ge_mod_assumed_pars_list = list(epsilon = 0.1), 
    ge_mod_true_pars_list = list(epsilon = 0.1)
) {
  # Check for missing allele frequencies and abort if any found
  if (any(is.na(D$Allele))) 
    stop("Some entries in Allele column in D are NA.  This cannot be so.  
         Please fix.  Aborting...")
  
  # Make list of lists, one list for each locus, of matrices representing
  # possible genotype-pair probabilities for each kinship requested
  mhlist <- long_markers_to_X_l_list(D = D, kappa_matrix = kappa_matrix)
  
  # Add matrices for each locus giving probability of observing each allele,
  # given each true allele, given an assumed genotype error model, for
  # calculation of plods, and a true one for simulation of genotype-pairs.
  mhlist2 <- insert_C_l_matrices(
    XL = mhlist, ge_mod_assumed = ge_mod_assumed, 
    ge_mod_true = ge_mod_true, 
    ge_mod_assumed_pars_list = ge_mod_assumed_pars_list, 
    ge_mod_true_pars_list = ge_mod_true_pars_list
  )
  
  # Add matrices representing probabilities of genotype-pairs for each
  # relationship, given assumed and true genotype error models.  Results of
  # multiplying matrices found in previous two steps.
  mhlist3 <- insert_Y_l_matrices(mhlist2)
  ret <- list(orig_data = D, loci = mhlist3)
  ckmr_class(ret)
}
