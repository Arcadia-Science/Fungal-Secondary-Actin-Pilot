library(tidyverse)
library(ggtree)
library(corHMM)

setwd('~/Documents/ArcadiaScience/github/Fungal-Secondary-Actin-Pilot/auxin_response_promoter/')

# Read in the data and tree
dat <- read.table('Organized_trait_data_andFSA_94species_Austin.csv', sep = ",", header = T, row.names = 1)
tree <- ladderize(read.tree('tree_trait_mapping.nwk'))

# Recode dat$auxinResponsivePromoter_count as a binary presence/absence, and recode the FSA variable in an equivalent manner
dat$FSA <- as.factor(gsub("Yes", "Present", dat$FSA) |> gsub(pattern = "No", replacement = "Absent"))
dat$auxinResponsivePromoter <- dat$auxinResponsivePromoter_count
dat$auxinResponsivePromoter[which(dat$auxinResponsivePromoter > 0)] <- "Present"
dat$auxinResponsivePromoter[which(dat$auxinResponsivePromoter == 0)] <- "Absent"
dat$auxinResponsivePromoter <- as.factor(dat$auxinResponsivePromoter)

# Now, reduce this down to just the species for which we have data for auxinResponsivePromoter
dat <- dat[!is.na(dat$auxinResponsivePromoter),]

# Prune the tree down to these species:
tree <- keep.tip(tree, dat$tip.label)

# Sort data to match tip labels
dat <- dat[match(tree$tip.label, dat$tip.label),]

### Markov models
# Now, let's explore the use of models of discrete trait evolution - Markov models. 
# These models (like the ones used for nucleotide or amino acid substitution) 
# estimate the rates at which some discrete character transitions from one state 
# to another. We can use these models to infer ancestral character states, and 
# to assess whether discrete characters have correlated evolutionary histories. 

# Now, let's try corHMMs test of correlated characters, to see
# whether this approach can detect an association between the
# two presence/absence states.
fsa_aux_corr_test <-
  corHMM::fitCorrelationTest(tree,
                             data = dat[c('tip.label', 'FSA', 'auxinResponsivePromoter')])

# Let's compare model fits:
corhmm_res <- data.frame(
  model = names(fsa_aux_corr_test),
  n_rate_cats = NA,
  n_params = NA,
  neg_log_lik = NA,
  AIC = NA,
  AICc = NA
)
for(model in 1:length(ls(fsa_aux_corr_test))){
  mod_name <- names(fsa_aux_corr_test[model])
  mod_cats <- fsa_aux_corr_test[[model]]$rate.cat
  mod_nparam <- length(na.omit(c(fsa_aux_corr_test[[model]]$solution)))
  mod_lik <- round(fsa_aux_corr_test[[model]]$loglik, 2)
  mod_aic <- round(fsa_aux_corr_test[[model]]$AIC, 2)
  mod_aicc <- round(fsa_aux_corr_test[[model]]$AICc, 2)
  
  corhmm_res[model,] <-
    c(mod_name, mod_cats, mod_nparam,
      mod_lik, mod_aic, mod_aicc)
}

# Lets see!
corhmm_res

# So, this would suggest that when looking across all species, a model in which 
# the presence/absence of FSA is independent of whether a species has a responsive 
# Auxin promoter. 

# Let's visualize the best-fit model:
best_fit <- which(corhmm_res$AICc == min(corhmm_res$AICc))

#Legend
#                1                 2                 3                 4
#  "absent_absent"  "absent_present"  "present_absent" "present_present"

# Lets rename the states something more useful
states <- c("FSA: 0, Auxin: 0", "FSA: 0, Auxin: 1", 
            "FSA: 1, Auxin: 0", "FSA: 1, Auxin: 1")
rate_states <- c("0_0", "0_1", "1_0", "1_1")
colnames(fsa_aux_corr_test[[best_fit]]$states) <- states
dimnames(fsa_aux_corr_test[[best_fit]]$solution) <- 
  list(rate_states, rate_states)

# First the estimated rate matrix:
corHMM::plotMKmodel(fsa_aux_corr_test[[1]])

# Now, look at inferred ancestral state reconstructions
corHMM::plotRECON(fsa_aux_corr_test[[best_fit]]$phy,
                  fsa_aux_corr_test[[best_fit]]$states,
                  show.tip.label = T, method = "marginal",
                  piecolors = c("#292928", "#596F74",
                                "#F7B846", "#F1E8DA"))

# Determine the color of the symbols based on the values in fsa_binary
fsa_color <- ifelse(dat$FSA == levels(dat$FSA)[1], 'white', 'black')
aux_color <- ifelse(dat$auxinResponsivePromoter == levels(dat$auxinResponsivePromoter)[1], 'white', 'black')

# Add observed states to the tips
tiplabels(pch = 21, bg = fsa_color,
          cex = 1, offset = 130)
tiplabels(pch = 21, bg = aux_color,
          cex = 1, offset = 140)

# and do this again, saving this out to a PDF
pdf('fsa_aux_binary_ancstate_recon.pdf',
    height = 9, width = 10)
corHMM::plotRECON(fsa_aux_corr_test[[best_fit]]$phy,
                  fsa_aux_corr_test[[best_fit]]$states,
                  piecolors = c("#292928", "#596F74",
                                "#F7B846", "#F1E8DA"),
                  show.tip.label = T)
fsa_color <- ifelse(dat$FSA == levels(dat$FSA)[1], 'white', 'black')
aux_color <- ifelse(dat$auxinResponsivePromoter == levels(dat$auxinResponsivePromoter)[1], 'white', 'black')
tiplabels(pch = 21, bg = fsa_color,
          cex = 1, offset = 150)
tiplabels(pch = 21, bg = aux_color,
          cex = 1, offset = 165)
dev.off()
# Save the fitted models and summary table out to file
saveRDS(fsa_aux_corr_test, file = 'fsa_aux_binary_corrtest_fitted_models.Rds')
write.table(corhmm_res, file = 'fsa_aux_binary_model_fits.tsv', sep = "\t", col.names = T, row.names = F, quote = F)












# And one that allows for hidden/latent states (i.e. unmeasured,
# covarying biological features that impact the rates of
# character state transitions).
fsa_hla_hmm2 <-
  corHMM(phy = ladderize(tree),
         data = pdat[c('species', 'fsa_binary', 'hla_binary')],
         rate.cat = 2, n.cores = 7, model = "ARD")

# Check the model fits:
fsa_hla_mm1$AICc
fsa_hla_hmm2$AICc

# Clearly the model with hidden states fits much better - this
# is not uncommon, as otherwise these markov models attribute
# all variation in transition rates to the sampled character.
# This assumption, while understandable, is not biologically
# realistic.

# Now, let's visualize this model, knowing that the character
# are as follows (ordered FSA_HLA):
# 1 = "absent_absent"
# 2 = "absent_present"
# 3 = "present_absent"
# 4 = "present_present"
corHMM::plotMKmodel(fsa_hla_hmm2)

# Great - now lets visualize the ancestral state reconstructions!
# Note: this is a huge tree, and so we'll be better off visualizing in a large pdf.
corHMM::plotRECON(fsa_hla_hmm2$phy,
                  fsa_hla_hmm2$states,
                  show.tip.label = F)

# Add observed states to the tips
# Determine the color of the symbols based on the values in fsa_binary
fsa_color <- ifelse(fsa_binary == levels(fsa_binary)[1], 'lightgrey', 'black')
hla_color <- ifelse(hla_binary == levels(hla_binary)[1], 'red', 'black')

# Add observed states to the tips
tiplabels(pch = 16, col = fsa_color,
          cex = 0.25, offset = 10)
tiplabels(pch = 16, col = hla_color,
          cex = 0.25, offset = 15)



#### Presence/Absence
# First fit models for a simplified character - the presence or absence of FSAs or HLAs.

# Fit a markov model with two states, binary presence absence for
# both FSA and auxinResponsivePromoter In this model, transitions can occur 
# between each of the four states (Absent_Absent, Absent_Present, Present_Present, 
# Present_Absent), excluding "double" state transitions (i.e. Absent_Absent -> 
# Present_Present).

# fsa_aux_mm1 <-
#   corHMM(phy = ladderize(tree),
#          data = dat[c('tip.label', 'FSA', 'auxinResponsivePromoter')],
#          rate.cat = 1, n.cores = 7, model = "ARD")

