# Save plots of PLODs to pdfs

# Open pdf file to autosave plot
pdf(paste0("analyses/", dataset_name, "/hsp_up_plods.pdf"))

# Set two plots per page
par(mfrow = c(2, 1))

# Plot plods
hist(
  hsp_up_plods, 
  main = "HSP vs UP PLODs for all samples",
  sub = paste(dataset_name, "->", n_loci, "loci"),
  xlab = "PLOD",
  breaks = 200,
)

# Plot expected values
abline(v = ev_up, col = 2)

# Add legend
legend(
  "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 0.75,
  legend = c(
    "Expected value given kinship", "Unrelated", "First cousin",
    "Avuncular", "Half-sibling", "Parent-offspring", "Self"
  )
)

# Plot uncommon values
hist(
  hsp_up_plods, 
  main = "Uncommon values suggest likely close-kin",
  xlab = "PLOD",
  breaks = 200, 
  ylim = c(0, 100)
)

# Plot expected values
abline(v = evs, col = 2:7)

# Save pdf of plot
dev.off()

# Loop over populations
for (pln.ind in 1:n_pops){
  # Find population name
  pln.nm <- levels(ppln)[pln.ind]
  
  # Open pdf file to autosave plot
  pdf(paste0("analyses/", dataset_name, "/hsp_up_plods_", pln.nm, ".pdf"))
  
  # Set two plots per page
  par(mfrow = c(2, 1))
  
  # Plot plods
  hist(
    hsp_up_plods[ppln[samp_inds_1] == pln.nm & ppln[samp_inds_2] == pln.nm], 
    main = paste("HSP vs UP PLODs within", pln.nm),
    sub = paste(dataset_name, "->", n_loci, "loci"),
    xlab = "PLOD",
    breaks = 200,
  )
  
  # Plot expected values
  abline(v = ev_up, col = 2)
  
  # Add legend
  legend(
    "topright", col = 1:7, lty = c(0, rep(1, 6)), cex = 0.75,
    legend = c(
      "Expected value given kinship", "Unrelated", "First cousin",
      "Avuncular", "Half-sibling", "Parent-offspring", "Self"
    )
  )
  
  # Plot uncommon values
  hist(
    hsp_up_plods[ppln[samp_inds_1] == pln.nm & ppln[samp_inds_2] == pln.nm], 
    main = "Uncommon values suggest likely close-kin",
    xlab = "PLOD",
    breaks = 200, 
    ylim = c(0, 50)
  )
  
  # Plot expected values
  abline(v = evs, col = 2:7)
  
  # Save pdf of plot
  dev.off()
}
