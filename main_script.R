# By Robin Aldridge-Sutton, for Emma Carroll, 3/11/22

# This is code to read in a dataset of genotyped samples, calculate the
# half-sibling versus unrelated pair PLODs for all pairs, and plot them to see
# how well they identify parent-offspring and half-sibling pairs.

# The dataset should be saved in a genind object called dataset_name.gi, in a
# file called dataset_name_genind.rdata.  

# Put the file in the datasets folder, set the dataset_name in the first line of
# code below, and the genind object variable name in the third and fourth lines
# of code. Set the working directory to the location of main_script.R, and run
# the code.

# A folder called dataset_name will be created in the analyses folder.  Use a
# different dataset name each time if you don't want to overwrite the old
# results!  The PLODs will be saved there as a named vector, and in a dataframe
# with their sample names and populations.  Their expected values at common
# kinships will also be saved there.

# Plots of the PLODs and their expected values at common kinships, both
# altogether, and separated by population, will also be saved there.

# There is already an analysis folder for each of the datasets you have given
# me.  The ones called srwhiMAF and srw.mhap have analyses based on all the
# loci, just loci with major allele frequencies below 0.8, and just loci with
# more than 4 alleles per locus, each of which numbered about 5000.  The one
# called srwqc is just based on all the loci.

# Set the dataset name, you need a different one for each analysis if you want
# to keep the old ones!
dataset_name <- "Change this to your dataset_name"

# Load saved genind object
load(paste0("datasets/", dataset_name, "_genind.rdata"))

# Store the genind object in known variable
dataset.gi <- Change_this_to_your_dataset_name.gi
rm(Change_this_to_your_dataset_name.gi)

# Read and summarise data for finding PLODs
source("code/prepare_data.R")

# Find half-sibling versus unrelated pair PLODs and save them in the analysis
# folder
source("code/find_plods.R")

# View the expected values and first few PLODs
head(hsp_up_plods)
head(plods_df)
evs

# Create plots of PLODs and save them in the analysis folder
source("code/plot_plods.R")
