library(phylolm)
library(phyr)
library(tidyverse)
library(ape)
library(ggtree)
library(phytools)

setwd('~/Documents/ArcadiaScience/Projects/manon_fungi_actin/')

# Read in the fsa/hla data, again removing the mollusc
dat <- read_csv('data/hla_fsa_dataset_for_tree_Austin.csv')
dat <- dat[-which(dat$Species == "Siphonaria sp. JEL0065"),]

# reformat, and also calculate the % Actins that are FSA
pdat <- data.frame(
  species = gsub(" ", "_", dat$Species),
  fsa_count = dat$fsa_count,
  hla_count = dat$hla_count,
  fsa_prop = dat$fsa_count / (dat$fsa_count + dat$hla_count),
  phylum = as.factor(dat$Phylum),
  phylum_simple = dat$Phylum
)
# Replace everything other than Ascomycota and Basidiomycota with "other", 
# since the former two dominate the dataset
pdat$phylum_simple[-which(pdat$phylum %in% c('Ascomycota', 'Basidiomycota'))] <- 'Other'
pdat$phylum_simple <- factor(pdat$phylum_simple, levels = c('Other', 'Ascomycota', 'Basidiomycota'))
pdat$fsa_prop[which(!is.finite(pdat$fsa_prop))] <- 0
rownames(pdat) <- 
  gsub(" ", "_", dat$Species)

# Read in the cladogram with species we want to retain
tree <- read.tree('data/tree_for_hla_and_fsa.nwk')
tree <- drop.tip(tree, "Siphonaria_sp._JEL0065")

# Make sure the tree is ultrametric and fully bifurcating, with no 0-length branches
tree$edge.length <- tree$edge.length+0.00001
tree <- multi2di(tree)
tree <- force.ultrametric(tree, method = 'extend')

# Read in the time calibrated tree of mushroom forming fungi:
mush_tt <- read.tree('~/Downloads/doi_10.5061_dryad.gc2k9r9__v1/Genome_chrono/r8s_default_calibration_scheme.tre')
mush_tt <- force.ultrametric(mush_tt, method = 'extend')
# rename tips in the mushroom time tree to match the timetree tree more closely:
corr_tiplabs <- sub("^(.*?_.*?)_.*", "\\1", tree$tip.label)
keep <- which(mush_tt$tip.label %in% sub("^(.*?_.*?)_.*", "\\1", tree$tip.label))
newnames <- tree$tip.label[which(sub("^(.*?_.*?)_.*", "\\1", tree$tip.label) %in% mush_tt$tip.label)]
newnames <- newnames[match(mush_tt$tip.label[keep], sub("^(.*?_.*?)_.*", "\\1", newnames))]
mush_tt$tip.label[keep] <- replace

# Now, see what congruifying the timetree tree looks like:
cong_tree <- 
  geiger::congruify.phylo(reference = mush_tt, target = tree, 
                          scale = 'treePL', ncores = 10)

# Great. This might be a bit bettwe behaved than the one from timetree. 
# Use this one instead. 
tree <- cong_tree$phy

# Begin by simply visualizing the trait on the tree with some inferred ancestral character states
fsa <- as.numeric(pdat$fsa_count)
hla <- as.numeric(pdat$hla_count)
prop <- as.numeric(pdat$fsa_prop)
names(fsa) <- rownames(pdat)
names(hla) <- rownames(pdat)
names(prop) <- rownames(pdat)
fsa_fit <- phytools::fastAnc(tree, log10(fsa+1))
hla_fit <- phytools::fastAnc(tree, log10(hla+1))
prop_fit <- phytools::fastAnc(tree, prop)

# Make a dataframe with trait values at the tips
td <- 
  data.frame(
    node = nodeid(tree, names(fsa)),
    fsa_count = log10(fsa+1),
    hla_count = log10(hla+1),
    prop_fsa = prop)

# Make a dataframe with estimated trait values at the nodes
nd <- data.frame(node = names(fsa_fit), 
                 fsa_count = c(fsa_fit),
                 hla_count = c(hla_fit),
                 prop_fsa = c(prop_fit))

# Combine these with the tree data for plotting with ggtree
d <- rbind(td, nd)
d$node <- as.numeric(d$node)

actin_tree <- full_join(tree, d, by = 'node')

fsa_tree_p <-
  ggtree(actin_tree, 
         aes(color = fsa_count), 
         ladderize = TRUE, continuous = "color", size = 0.5) +
  scale_color_viridis_c(option = "inferno", 
                        name = "# FSA Copies", 
                        breaks = c(0, 0.30103, 0.4771213,
                                   0.60206, 0.69897),
                        labels = c(0, 1, 2, 3, 4)) + 
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        legend.key.width = unit(0.9,"cm"), 
        legend.title = element_text(size = 14)) +
  coord_cartesian(clip="off")

hla_tree_p <-
  ggtree(actin_tree, 
         aes(color = hla_count), 
         ladderize = TRUE, continuous = "color", size = 0.5) +
  scale_color_viridis_c(option = "inferno", 
                        name = "# HLA Copies", 
                        breaks = c(0, 0.30103, 0.4771213,
                                   0.60206, 0.69897),
                        labels = c(0, 1, 2, 3, 4)) + 
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        legend.key.width = unit(0.9,"cm"), 
        legend.title = element_text(size = 14)) +
  coord_cartesian(clip="off")
prop_fsa_tree_p <-
  ggtree(actin_tree, 
         aes(color = prop_fsa), 
         ladderize = TRUE, continuous = "color", size = 0.5) +
  scale_color_viridis_c(option = "inferno", 
                        name = "Prop. Actins FSA") + 
  theme(legend.position = 'top',
        legend.text = element_text(size = 12),
        legend.key.width = unit(0.9,"cm"), 
        legend.title = element_text(size = 14)) +
  coord_cartesian(clip="off")

actin_plt <- cowplot::plot_grid(fsa_tree_p, hla_tree_p, prop_fsa_tree_p, nrow = 1)

ggsave(actin_plt, file = 'fsa_hla_actin_ancestral_reconstructions.pdf',
       height = 8, width = 12)

########################
# Great, now formally assess the extent to which FSA counts are 
# predicted by the counts of HLA

# First, a simple glm, not accounting for phylogeny
fit_gaus <- glm(fsa_count ~ hla_count, data = pdat)

# Then, one that explicitly accounts for phylogeny:
fit_pgaus <- pglmm_compare(fsa_count ~ hla_count, data = pdat, phy = tree)

# Note, we can fit models using either a poisson or binomial distribution. 
# Poisson is useful for count data, but assumes the mean and variance of the 
# distribution is the same, something that is not true of our data. 
# Binomial distributions are useful for binary success/failure type data. 
fit_ppois <- pglmm_compare(fsa_count ~ hla_count, data = pdat, phy = tree, family = "poisson")

# As you can see, the poisson fits the data more poorly than does the gaussian.
# Consequently, we just proceed with the gaussian
summary(fit_pgaus)
summary(fit_ppois)

# Now plot this, adding both regression 
# slopes to show the effect of accounting for phylogeny and include 
# model coefficients and p-values

# Define the intercept, beta (slope), and p-value text for each model
non_phylo_text <- paste("Non-phylogenetic Model\n",
                        "Intercept: ", round(fit_gaus$coefficients[1], 4), "\n",
                        "Beta: ", round(fit_gaus$coefficients[2], 4), "\n",
                        "p-value: ", round(summary(fit_gaus)$coefficients[2,4], 5))

phylo_text <- paste("Phylogenetic Model\n",
                    "Intercept: ", round(fit_pgaus$B[1], 4), "\n",
                    "Beta: ", round(fit_pgaus$B[2], 4), "\n",
                    "p-value: ", round(fit_pgaus$B.pvalue[2], 5))

fsa_hsa_corrplot <- 
  ggplot(data = pdat, aes(y = fsa_count, x = hla_count)) + 
  geom_point(alpha = 0.5, 
             position = position_jitter(width = 0.1, 
                                        height = 0.1)) + 
  geom_abline(aes(intercept = fit_gaus$coefficients[1], slope = fit_gaus$coefficients[2], color='Non-phylogenetic'), size = 1) +
  geom_abline(aes(intercept = fit_pgaus$B[1], slope = fit_pgaus$B[2], color='Phylogenetic'), size = 1) +
  scale_color_manual(name='Regression Model',
                     values=c('Non-phylogenetic'='black', 
                              'Phylogenetic'='red')) +
  theme_bw(base_size = 14) + 
  xlab('# HLA Copies') +
  ylab('# FSA Copies') +
  annotate("text", x = min(pdat$hla_count), y = max(pdat$fsa_count), hjust = 0, vjust = 1,
           label = non_phylo_text, color = "black", size = 3) +
  annotate("text", x = min(pdat$hla_count), y = max(pdat$fsa_count) - 0.75, hjust = 0, vjust = 1,
           label = phylo_text, color = "red", size = 3) +
  theme(legend.position = 'top')

ggsave(fsa_hsa_corrplot, file = 'fsa_hsa_correlation_plt.pdf', height = 7, width = 7)

# Now, what if we account for an effect of phylum on these predictions (i.e. 
# do species within different phyla experience a different relationship 
# between FSA/HLA counts?)
fit_pgaus_phylum <- pglmm_compare(fsa_count ~ hla_count * phylum_simple, data = pdat, phy = tree)
summary(fit_pgaus_phylum)

# Again, check to see whether we observe the same pattern with a non-phylogenetic
# glm. 
# First, a simple glm, not accounting for phylogeny
fit_gaus_phylum <- glm(fsa_count ~ hla_count * phylum_simple, data = pdat)
summary(fit_gaus_phylum)

# So, it does seem to be that phylum plays a role here. Specifically, we see that
# the relationship is still significant for "other". Although the relationship
# does not differ between Basidiomycota and Other (still positive), Ascomycota
# actually has a significantly different relationship that either Basidiomycota
# and everything else - the relationship is now negative! Furthermore, this is 
# distinct from what we observe when looking at a non-phylogenetic model!

# Together, this suggests that in Ascomycota, there is a tradeoff between the 
# counts of HLAs and FSAs

# Again, plot this. 
other_int <- fit_gaus_phylum$coefficients[1]
other_slope <- fit_gaus_phylum$coefficients[2]
other_pint <- fit_pgaus_phylum$B[1]
other_pslope <- fit_pgaus_phylum$B[2]

fsa_hsa_phylum_corrplot <- 
  ggplot(data = pdat, aes(y = fsa_count, x = hla_count, color = phylum_simple)) + 
  geom_point(alpha = 0.5, 
             position = position_jitter(width = 0.1, 
                                        height = 0.1)) + 
  geom_abline(aes(intercept = other_int, slope = other_slope, color='Other', linetype = "Non-phylogenetic"), size = 1) +
  geom_abline(aes(intercept = other_int+fit_gaus_phylum$coefficients[3], slope = other_slope+fit_gaus_phylum$coefficients[5], color='Ascomycota', linetype = "Non-phylogenetic"), size = 1) +
  geom_abline(aes(intercept = other_int+fit_gaus_phylum$coefficients[4], slope = other_slope+fit_gaus_phylum$coefficients[5], color='Basidiomycota', linetype = "Non-phylogenetic"), size = 1) +
  geom_abline(aes(intercept = other_pint[1], slope = other_pslope, color='Other', linetype = "Phylogenetic"), size = 1) +
  geom_abline(aes(intercept = other_pint+fit_pgaus_phylum$B[3], slope = other_pslope+fit_pgaus_phylum$B[5], color='Ascomycota', linetype = "Phylogenetic"), size = 1) +
  geom_abline(aes(intercept = other_pint+fit_pgaus_phylum$B[4], slope = other_pslope+fit_pgaus_phylum$B[6], color='Basidiomycota', linetype = "Phylogenetic"), size = 1) +
  scale_color_manual(name='Phylum',
                     values=c('Other' = 'black', 
                              'Ascomycota' = 'red',
                              "Basidiomycota" = "blue")) +
  scale_linetype_manual(name='Model type',
                        values=c('Non-phylogenetic' = 2, 
                                 'Phylogenetic' = 1)) +
  theme_bw(base_size = 14) + 
  xlab('# HLA Copies') +
  ylab('# FSA Copies') +
  theme(legend.position = 'top', 
        legend.box = "vertical", 
        legend.direction = "horizontal",
        legend.spacing.y = unit(0.01, "cm"))

ggsave(fsa_hsa_phylum_corrplot, file = 'fsa_hsa_phylum_correlation_plt.pdf', height = 7, width = 7)

##############################################################
# Now let's just model this as a binary presence/absence
##############################################################
# Prepare the data
pdat$fsa_binary <- as.factor(pdat$fsa_count > 0)
pdat$hla_binary <- as.factor(pdat$hla_count > 0)
pdat$fsa_binary <- 
  gsub("TRUE", "present", pdat$fsa_binary) |>
  gsub(pattern = "FALSE", replacement = "absent") |> 
  as.factor()
pdat$hla_binary <- 
  gsub("TRUE", "present", pdat$hla_binary) |>
  gsub(pattern = "FALSE", replacement = "absent") |> 
  as.factor()
fsa_binary <- pdat$fsa_binary
hla_binary <- pdat$hla_binary
names(fsa_binary) <- rownames(pdat)
names(hla_binary) <- rownames(pdat)

# Get the common ancestor node for the two large phyla - Basidiomycota and 
# Ascomycota - everything else will be treated in a third group:
phyla <- unique(pdat$phylum_simple)
phyla_mrca <- c()
for(p in 1:length(phyla)){
  phyla_mrca[p] <- getMRCA(tree, rownames(pdat)[which(pdat$phylum_simple == phyla[p])])
  names(phyla_mrca)[p] <- as.character(phyla[p])
}

phyla_tree <- 
  paintSubTree(ladderize(tree), node = 1, 
               state = 'Other', 
               anc.state = 'Other', 
               stem = T)
phyla_tree <- 
  paintSubTree(phyla_tree, node = phyla_mrca[['Ascomycota']], 
               state = 'Ascomycota', 
               anc.state = 'Ascomycota')
phyla_tree <- 
  paintSubTree(phyla_tree, node = phyla_mrca[['Basidiomycota']], 
               state = 'Basidiomycota', 
               anc.state = 'Basidiomycota')

# Now plot to see how this looks
plotSimmap(phyla_tree, fsize = 0.001)


library(corHMM)
# Fit a markov model with two states, binary presence absence for
# both FSA and HLA. In this model, transitions can occur between 
# each of the four states (0_0, 0_1, 1_1, 1_0), excluding
# "double" state transitions (i.e. 0_0 -> 1_1). 
fsa_hla_mm1 <- 
  corHMM(phy = ladderize(tree), 
         data = pdat[c('species', 'fsa_binary', 'hla_binary')],
         rate.cat = 1, n.cores = 7, model = "ARD")

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
corHMM::plotRECON(fsa_hla_hmm2$phy, 
                  fsa_hla_hmm2$states, 
                  show.tip.label = F)

# Now, let's try corHMMs test of correlated characters, to see
# whether this approach can detect an association between the 
# two presence/absence states. 
fsa_hla_corr_test <- 
  corHMM::fitCorrelationTest(ladderize(tree),
                             data = pdat[c('species', 'fsa_binary', 'hla_binary')],
                             simplified_models = TRUE)























tmp <- tree
tmp$edge.length <- tmp$edge.length + 0.00001
fsa_bin_mk <- 
  geiger::fitDiscrete(phy = phyla_tree, dat = fsa_binary, model = "ARD", transform = "lambda", ncores = 10)
fsa_bin_qmat <- geiger::as.Qmatrix.gfit(fsa_bin_mk)

fsa_bin_mk <- 
  phytools::fitMk.parallel(phyla_tree, fsa_binary, model = "ARD", 
                  ncores = 10)
# Great, now fit a model wherein the rates of character state transition 
# (HLA/FSA presence/absence differs by phyla)
fsa_bin_phylum_mk <- 
  phytools::fitmultiMk(phyla_tree, fsa_binary, model = "ARD", 
                       ncores = 10)


fsa_bin_simhist <- 
  phytools::sim.history(tree, Q = matrix(fsa_bin_qmat, nrow = 2), 
                        nsim = 100)
pd <- summary(fsa_bin_simhist, plot=FALSE)

nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)


# function to compute the states
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}
XX<-sapply(fsa_bin_simhist,foo)
pies<-t(apply(XX,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=c("0","1"),Nsim=100))
# done computing the states

# code to plot the tree
mtrees<-rescaleSimmap(fsa_bin_simhist,1)
cols<-c("blue","red"); names(cols)<-c(0,1)
par(mar=rep(0.1,4))
plot.phylo(mtrees[[1]],plot=FALSE,no.margin=T)
plotSimmap(mtrees[[1]],cols,pts=FALSE,lwd=3,ftype="off", add=TRUE)
text(1,1:length(mtrees[[1]]$tip),mtrees[[1]]$tip.label, pos=4,font=3)
nodelabels(pie=pies,cex=0.6,piecol=cols)

#Let's extract some basic information from our maps.
#Get numbers of changes for each trait
# function counts transitions from a mapped history
# written by Liam J. Revell 2013
countSimmap<-function(tree,states=NULL,message=TRUE){
  n<-sum(sapply(tree$maps,length))-nrow(tree$edge)
  if(is.null(states)) states<-colnames(tree$mapped.edge)
  m<-length(states)	
  TT<-matrix(NA,m,m,dimnames=list(states,states))
  foo<-function(map,a,b){
    if(length(map)==1) zz<-0
    else {
      zz<-0; i<-2
      while(i<=length(map)){
        if(names(map)[i]==b&&names(map)[i-1]==a) zz<-zz+1
        i<-i+1
      }
    } 
    return(zz)
  }
  for(i in 1:m) for(j in 1:m)
    if(i==j) TT[i,j]<-0
  else TT[i,j]<-sum(sapply(tree$maps,foo,a=states[i],b=states[j]))
  return(list(N=n,Tr=TT,message=c(
    "N is the total number of character changes on the tree",
    "Tr gives the number of transitions from row state->column state")))
}

countSimmap(fsa_bin_simhist)
phenogram(GAanolisMicroSimmap[[1]], micro, colors=micro_colors)
map.overlap(GAanolisMicroSimmap[[1]],GAanolisMicroSimmap[[2]])
map.overlap(GAanolisMacroSimmap[[10]],GAanolisMacroSimmap[[21]])

obj <- describe.simmap(fsa_bin_simhist, plot = FALSE)
plot(obj, type="fan", fsize = 0.001, )



# fit the MK model for each trait individually
fsa_bin_mk_phytools <- 
  phytools::fitMk.parallel(tree, fsa_binary, model = "ARD", 
                           ncores = 10)
hla_bin_mk <- 
  phytools::fitMk.parallel(tree, hla_binary, model = "ARD", 
                           ncores = 10)


# And use sfreemap to analytically generate simmaps (predictions of ancestral 
# discrete character states under fitted evolutionary model)
library(sfreemap)
fsa_bin_sfreemap <- 
  sfreemap(tree, tip_states = fsa_binary, 
           model = "ARD")
hla_bin_sfreemap <- 
  sfreemap(tree, tip_states = hla_binary, 
           model = "ARD")
fsa_bin_map <- 
  map_posterior_distribution(fsa_bin_sfreemap[[1]], 
                             fsa_bin_sfreemap, 
                             parallel = TRUE)

plot_distribution_tree(fsa_bin_sfreemap)


fsa_hla_pagel_depx <- 
  phytools::fitPagel(tree, x = fsa_binary, y = hla_binary, 
                     method = 'fitDiscrete', model = "ARD", depvar = "x")
fsa_hla_pagel_depxy <- 
  phytools::fitPagel(tree, x = fsa_binary, y = hla_binary, 
                     method = 'fitDiscrete', model = "ARD", depvar = "xy")
fsa_hla_pagel_depxy <- 
  phytools::fitPagel(tree, x = fsa_binary, y = hla_binary, 
                     method = 'fitDiscrete', model = "ARD", depvar = "y")
# And visualize: simulate the discrete character on the tree 
# using the fitted model many times. Do this in parallel,
# as the simmap function scales poorly to large trees
library(snow)
cl <- makeSOCKcluster(rep("localhost", 12))
smap_fsa_bin <- 
  clusterApply(cl, x = replicate(10, fsa_binary, simplify = FALSE), 
               fun = make.simmap, tree = tree, 
               Q = as.Qmatrix(fsa_bin_mk), pi = c(1,0), 
               nsim = 10)

smap_fsa_bin <- do.call(c, smap_fsa_bin)
if(!("multiSimmap" %in% class(smap_habitat.mc))) 
  class(smap_habitat.mc) <- 
  c("multiSimmap", class(smap_habitat.mc))
stopCluster(cl)
# 
# fsa_simmap <- make.simmap(tree, fsa_binary, model = "ARD", 
#                           nsim = 100)
# hla_simmap <- make.simmap(tree, hla_binary, model = "ARD", 
#                           nsim = 100)

plot_simmap <- function(simmap = NULL, states = NULL, state_name = NULL, contmap = NULL, colors = NULL){
  tree_p <- 
    ggtree(simmap[[1]], branch.length = 'none')
  
  anc_states <- as.data.frame(summary(simmap)$ace)
  anc_states$node <- rownames(anc_states)
  node_indices <- match(tree_p$data$label[1:length(simmap[[1]]$tip.label)], anc_states$node)
  spp_nodes <- na.omit(match(anc_states$node, tree_p$data$label[1:length(simmap[[1]]$tip.label)]))
  anc_states$node[node_indices] <- spp_nodes
  nstates <- length(unique(states))
  if(is.null(colors)){
    if(nstates == 2){
      colors <- RColorBrewer::brewer.pal(3, "Set2")
      colors <- colors[c(1,3)]
    }else{
      colors <- RColorBrewer::brewer.pal(nstates, "Set2")
    }
  }
  node_pies <- nodepie(anc_states, cols = 1:nstates,
                       outline.color = 'black', 
                       outline.size = 0.1)
  
  for(i in 1:length(node_pies)){
    node_pies[[i]] <- 
      node_pies[[i]] + 
      scale_fill_manual(values = colors)
  }
  
  taxa_tibble <- tibble(
    tip = simmap[[1]]$tip.label,
    state = states)
  
  if(is.null(contmap)){
    tree_p <- tree_p %<+% taxa_tibble + 
      geom_tippoint(aes(fill = state), pch = 21) +
      scale_fill_manual(name = state_name, values = colors) + 
      guides(fill = guide_legend(override.aes = list(size = 5), nrow = 2)) + 
      geom_inset(node_pies, vjust = 0.1, hjust = 0.02, 
                 width = .055, height = .055) +
      theme(legend.position = 'top',
            plot.margin = margin(5,140,10,5),
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14)) +
      geom_tiplab(offset = 0.2, size = 3, color = "black") + 
      coord_cartesian(clip="off")
  }else{
    tree_p <- contmap %<+% taxa_tibble + 
      geom_tippoint(aes(fill = state), pch = 21) +
      scale_fill_manual(name = state_name, values = colors) + 
      guides(fill = guide_legend(override.aes = list(size = 5), nrow = 2)) + 
      geom_inset(node_pies, vjust = 0.1, hjust = 0.02, 
                 width = .055, height = .055) +
      theme(legend.position = 'top',
            plot.margin = margin(5,140,10,5),
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14))
  }
  
  return(tree_p)
}

# Now, we can use the observed and simulated discrete 
# characters to ask whether FSA presence/absence is 
# predicted by HLA presence/absense. 
# Specifically, we will use Pagels's test of 




fsa_cols <- rev(RColorBrewer::brewer.pal(4, "Set2"))
hla_cols <- RColorBrewer::brewer.pal(4, "Paired")
trophic_tree_p <- plot_simmap(trophic_simmap, trophic, "Trophic Mode", contmap = dt_nov_plt, colors = trophic_cols)
environ_tree_p <- plot_simmap(environ_simmap, environ, "Environment", contmap = dt_nov_plt, colors = environ_cols)















# Create a Phylogenetic Variance-Covariance Matrix:
A_inv <- solve(vcv(tree))

# Convert to numeric
pdat$fsa_count <- as.integer(pdat$fsa_count)
pdat$hla_count <- as.integer(pdat$hla_count)

# Set the prior
prior_p <- list(
  R = list(V = 2,  # Single variance for the Poisson residuals
           nu = 0.002),
  G = list(list(V = 50,  # Single variance for the random effect
                nu = 0.002, 
                alpha.mu = 0, alpha.V = 1000, 
                Ginv = A_inv))
)
model_g <- MCMCglmm(
  fsa_count ~ hla_count,  # Remove the trait specification
  random = ~species,
  data = pdat,
  prior = prior_p,
  family = "gaussian",
  burnin = 5000,
  nitt = 55000,
  thin = 100,
  verbose = TRUE
)
model_p <- MCMCglmm(
  fsa_count ~ hla_count,  # Remove the trait specification
  random = ~species,
  data = pdat,
  prior = prior_p,
  family = "poisson",
  burnin = 5000,
  nitt = 55000,
  thin = 100,
  verbose = TRUE
)

model_p <- MCMCglmm(
  fsa_count ~ hla_count,  # Remove the trait specification
  random = ~species,
  data = pdat,
  prior = prior_p,
  family = "poisson",
  burnin = 50000,
  nitt = 550000,
  thin = 1000,
  verbose = TRUE
)

# Now plot! Add both regression slopes to show the effect of 
# accounting for phylogeny
ggplot(data = pdat, aes(y = fsa_count, x = hla_count)) + 
  geom_point(alpha = 0.5, 
             position = position_jitter(width = 0.15, 
                                        height = 0.15)) + 
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2], color = "black", size = 1) +
  geom_abline(intercept = summary(model_p)$solutions[1], slope = summary(model_p)$solutions[2], color = "red", size = 1) +
  theme_bw()


prior_hup <- list(
  R = list(V = diag(2),  # Increasing the variances in the diagonal
           nu = 0.002),
  G = list(list(V = 50,  # Increasing the variance for the random effect
                nu = 0.001, 
                alpha.mu = 0, alpha.V = 1000, 
                Ginv = A_inv)))

model_hup <- MCMCglmm(
  fsa_count ~ trait - 1 + 
    at.level(trait, 1):hla_count,
  rcov = ~idh(trait):units, 
  random = ~species,
  data = pdat,
  prior = prior_hup,
  family = "hupoisson", 
  burnin = 50000,
  nitt = 550000, 
  thin = 1000,
  verbose = TRUE
)

pdat <- data.frame(
  sp = gsub(" ", "_", dat$Species),
  species = gsub(" ", "_", dat$Species),
  fsa_count = dat$fsa_count,
  hla_count = dat$hla_count,
  phylum = dat$Phylum
)
rownames(pdat) <- 
  gsub(" ", "_", dat$Species)

# Identify terminal branches and add the tiny length
terminal_nodes <- which(tree$edge[, 2] <= Ntip(tree))
tree$edge.length[terminal_nodes] <- tree$edge.length[terminal_nodes] + 0.00001

# zero inflated model with phylogenetic effect using MCMCglmm
library(MCMCglmm)

# Specify the prior
prior_zip_nophy <- list(
  R = list(V = diag(2),  # Increasing the variances in the diagonal
           nu = 0.002))
prior_p <- list(
  R = list(V = 2,  # Single variance for the Poisson residuals
           nu = 0.002),
  G = list(list(V = 50,  # Single variance for the random effect
                nu = 0.002, 
                alpha.mu = 0, alpha.V = 1000, 
                Ginv = A_inv))
)

# Fit the models (zero-inflated poisson, hurdle poisson, and regular poisson):
model_zip <- MCMCglmm(
  fsa_count ~ trait - 1 + 
    at.level(trait, 1):hla_count,
  rcov = ~idh(trait):units, 
  random = ~species,
  data = pdat,
  prior = prior_zip,
  family = "zipoisson", 
  burnin = 10000,
  nitt = 50000,
  verbose = TRUE
)

model_hup <- MCMCglmm(
  fsa_count ~ trait - 1 + 
    at.level(trait, 1):hla_count,
  rcov = ~idh(trait):units, 
  random = ~species,
  data = pdat,
  prior = prior_zip,
  family = "hupoisson", 
  burnin = 10000,
  nitt = 50000,
  verbose = TRUE
)

model_hup_nophy <- MCMCglmm(
  fsa_count ~ trait - 1 + 
    at.level(trait, 1):hla_count,
  rcov = ~idh(trait):units, 
  data = pdat,
  prior = prior_zip_nophy,
  family = "hupoisson", 
  burnin = 10000,
  nitt = 50000,
  verbose = TRUE
)

model_p <- MCMCglmm(
  fsa_count ~ hla_count,  # Remove the trait specification
  random = ~species,
  data = pdat,
  prior = prior_p,
  family = "poisson",
  burnin = 10000,
  nitt = 50000,
  verbose = TRUE
)

DICs <- c(model_p$DIC, model_zip$DIC, model_hup$DIC)
model_fits <- 
  data.frame(model = c("poisson", 'zero-inflated poisson', 'hurdle poisson'),
             DIC = DICs, deltaDIC = DICs - min(DICs))
zip_preds <- predict(model_zip)
hup_preds <- predict(model_hup)
hup_phy_preds <- predict(model_hup_phy)

compdat <- caper::comparative.data(tree, pdat, names.col = 'species')
m <- caper::pgls(fsa_count ~ hla_count, compdat)

# and format things for visualizing the credible interval:
hup_preds <- predict.MCMCglmm(model_hup, newdata = pdat, type = "response", level = 0.95, verbose = T)
pdat$predicted <- hup_preds$fit
pdat$lwr <- hup_preds$lower
pdat$upr <- hup_preds$upper

ggplot(data = pdat, aes(y = fsa_count, x = hla_count)) + 
  geom_point(alpha = 0.5, 
             position = position_jitter(width = 0.15, 
                                        height = 0.15)) + 
  stat_smooth(method = 'lm') +
  geom_abline(intercept = -0.29, slope = 0.3915, color = "red", size = 1) +
  geom_abline(intercept = -0.29, slope = 0.3915, color = "red", size = 1) +
  theme_bw()

# Now conduct ancestral state reconstruction of the two traits and visualize:
fsa <- as.numeric(pdat$fsa_count)
hla <- as.numeric(pdat$hla_count)
names(fsa) <- rownames(pdat)
names(hla) <- rownames(pdat)
fsa_fit <- phytools::fastAnc(tree, fsa)
hla_fit <- phytools::fastAnc(tree, hla)

fsa_sym <- fitMk.parallel(tree, fsa, ncores = 14)
fsa_simmap <- make.simmap(tree, fsa, model = "SYM", nsim = 100)


library(corHMM)
dat_hmm <- pdat[,c('species', 'fsa_count', 'hla_count')]
dat_hmm$fsa_count <- 
  factor(dat_hmm$fsa_count, 
         levels = unique(sort(dat_hmm$fsa_count)), 
         ordered = T)
dat_hmm$hla_count <- 
  factor(dat_hmm$hla_count, 
         levels = unique(sort(dat_hmm$hla_count)), 
         ordered = T)
# Get the full matrix of potential state transitions
full_stat_mat <- getStateMat4Dat(dat_hmm)
mk_1rate <- corHMM(tree, dat_hmm, rate.cat = 1, model = 'SYM', n.cores = 14)

# Now drop extreme state transitions. In one, we will only allow for single
get_values_to_remove <- function(state_matrix, max_state_changes) {
  # Extract the state values using the provided legend
  states <- t(sapply(strsplit(state_matrix$legend, "_"), as.numeric))
  
  # Compute the state changes for each transition
  state_changes <- matrix(0, nrow(state_matrix$rate.mat), ncol(state_matrix$rate.mat))
  for (i in 1:nrow(state_matrix$rate.mat)) {
    for (j in 1:ncol(state_matrix$rate.mat)) {
      state_changes[i, j] <- sum(abs(states[i,] - states[j,]))
    }
  }
  
  # Override max_state_changes restriction for the nearest states
  for (i in 1:nrow(state_matrix$rate.mat)) {
    # Indices of the nearest states
    nearest_states <- order(state_changes[i,])[1]
    # If the nearest state is more than max_state_changes away
    if (state_changes[i, nearest_states] > max_state_changes) {
      state_changes[i, nearest_states] <- max_state_changes
    }
  }
  
  # Find the indices where the transitions exceed the max state changes
  indices_to_remove <- which(state_changes > max_state_changes, arr.ind = TRUE)
  
  # Extract the actual values at those cells
  values_to_remove <- state_matrix$rate.mat[indices_to_remove]
  
  # Filter out zeros
  values_to_remove <- values_to_remove[values_to_remove != 0]
  
  # Return the values
  return(values_to_remove)
}

state_1_mat <- dropStateMatPars(full_stat_mat$rate.mat, get_values_to_remove(state_matrix = full_stat_mat, 1))
state_2_mat <- dropStateMatPars(full_stat_mat$rate.mat, get_values_to_remove(state_matrix = full_stat_mat, 2))
state_3_mat <- dropStateMatPars(full_stat_mat$rate.mat, get_values_to_remove(state_matrix = full_stat_mat, 3))

# Fit the model with symmetrical rates for the single state transition matrix, 
# and two rate categories modeled
ss_hmm <- 
  corHMM(tree, dat_hmm, rate.mat = state_1_mat, rate.cat = 1, model = "SYM", 
         node.states = "marginal", fixed.nodes=FALSE, p=NULL, root.p="yang", 
         ip=NULL, nstarts=0, n.cores=14, get.tip.states = FALSE, 
         lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
         upper.bound = 100, opts=NULL)




model_hup_phy <- MCMCglmm(
  fsa_count ~ trait - 1 + 
    at.level(trait, 1):hla_count +
    at.level(trait, 1):phylum,
  rcov = ~idh(trait):units, 
  random = ~species,
  data = pdat,
  prior = prior_zip,
  family = "hupoisson", 
  burnin = 10000,
  nitt = 50000,
  verbose = TRUE
)


ggplot(data = dat, aes(y = fsa_count, x = hla_count)) + 
  geom_point(alpha = 0.5, 
             position = position_jitter(width = 0.15, 
                                        height = 0.15)) + 
  geom_abline(intercept = 1.9990835, slope = -0.3992514, color = "red", size = 1) +
  stat_smooth(method = 'lm') +
  theme_bw()

fsa <- pdat$fsa_count
hla <- pdat$hla_count
names(fsa) <- rownames(pdat)
names(hla) <- rownames(pdat)

fsa_contmap <- contMap(tree, fsa, fsize=c(0,0.8), lwd = 0.75)
hla_contmap <- contMap(tree, hla, fsize=c(0,0.8), lwd = 0.75)

cpdat <- caper::comparative.data(tree, pdat, names.col = sp)
fit_pgls <- caper::pgls(fsa_count ~ hla_count, data = cpdat)








# Read in information about the Kohler time calibrated fungal tree
kohler_tax <- read.table('kohler_et_al_strain_abbrevs.txt', sep = "\t", header = T)
kohler_tree <- read.tree('kohler_et_al_fungal_chronogram.newick')

kohler_tax <- kohler_tax[match(kohler_tree$tip.label, kohler_tax$abbreviation),]

# Keep a single representative per-family
reps <- c('Aspnid1', 'Agabi_varb', 'Phybl2', 'Batde5')
rep_phyl <- c('Ascomycota', 'Basidiomycota', 'Mucoromycotina', 'Chytridiomycetes')
kohler_tree <- keep.tip(kohler_tree, reps)

# Rename the kohler tree with family name
kohler_tree$tip.label <- rep_phyl
kohler_tree <- force.ultrametric(kohler_tree, method = 'extend')
# Now congruify
cong_kohler <- 
  geiger::congruify.phylo(reference = kohler_tree, scale = 'treePL',
                          target = tree, taxonomy = pdat[,c('species', 'phylum')], 
                          ncores = 10)

#mush_tt <- keep.tip(mush_tt, mush_tt$tip.label[which(mush_tt$tip.label %in% sub("^(.*?_.*?)_.*", "\\1", tree$tip.label))])
mush_tt$tip.label[which(mush_tt$tip.label %in% sub("^(.*?_.*?)_.*", "\\1", tree$tip.label))] <- 
  
  
  
  drop <- !which(mush_tt$tip.label %in% sub("^(.*?_.*?)_.*", "\\1", tree2$tip.label))
mush_tt
tree2$tip.label <- sub("^(.*?_.*?)_.*", "\\1", tree2$tip.label)

# Make sure the tree is ultrametric and fully bifurcating, with no 0-length branches
tree <- force.ultrametric(tree, method = 'extend')
tree$edge.length <- tree$edge.length + 0.001
tree <- multi2di(tree)

