############################################################
#workspace prep
############################################################

#expand memory if necessary
memory.limit(10000)

#set directory
setwd(" ")

#load packages
library(brms)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(mice)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

#MCMC settings
n_iter <- 3000
n_warm <- 1000
n_chains <- 4

#load datasets
dataR = read.csv("dataR.csv")
dataG = read.csv("dataG.csv")
dataGR = read.csv("dataGR.csv")

#add redundant column for phylogenetic effects
dataR$phylo = dataR$Genus_species
dataG$phylo = dataG$Genus_species
dataGR$phylo = dataGR$Genus_species

#load phylogenies
#multiple to capture uncertainty
library(ape)
trees = read.nexus("vert phylo.nex")

#create consensus phylogeny
c_tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")

#directory for results (saving brms and stan models 
#will take GBs of memory)
setwd(" ")

############################################################
#Save image of phylogeny
############################################################

#without tips
png("phyloivso_radial.png", width = 5, height = 5, units = "in", res = 600)
plot.phylo(c_tree, type = "radial", show.tip.label = FALSE, edge.width = 2)
dev.off()

#with contour
library(phytools)
temp = aggregate(IVSO_prop ~ Genus_species, data =  dataGR, FUN = mean)
x = temp[,2]
names(x) = temp[,1]

png("phyloivso_radial_cont.png", width = 10, height = 10, units = "in", res = 600)
obj = contMap(phytools::force.ultrametric(c_tree, method="extend"),
              x, plot=FALSE, type ="fan", lwd=0.5)
n = length(obj$cols)
colfunc = colorRampPalette(c("purple", "red"))
obj$cols[1:n]<-colfunc(n)
plot(obj, type ="fan", ftype = "off", lwd = 3, spread.labels = TRUE)        
dev.off()

############################################################
#histograms of raw data
############################################################
#create copy of dataframe for manipulation
df1 = dataGR

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]


#representation of superfamilies
library(ggplot2)

p.supf =
ggplot(df1, aes(x = superfamily))+geom_bar(aes(y = (..count..)/sum(..count..)))+
  ylab("% of field studies\n")+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 35, vjust = 0.6),
          axis.title.x = element_blank(),
          #axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))

ggsave("p_supf.png", p.supf)  

#organize predictors
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)

#get total counts of main (most frequent) SO
table(df1$Main1)
table(df1$Main2) #varies because of unresolved main SO

dfp = data.frame(cat = factor(c("solitary","MF","MFF","FMM","FFMM"), levels = 
                                c("solitary","MF","MFF","FMM","FFMM")),
                 valh = c(16, 96, 140, 9, 269), 
                 vall = c(15, 88, 119, 7, 240)) 

colors = brewer.pal(n = 8, "Dark2")

hist1 =
ggplot(dfp, aes(x = cat, y = valh, alpha = 0.8)) + 
  geom_col(color = "black", fill = "grey")+
  xlab("\nPrimary social organization")+
  ylab("Number of populations\n")+
  #scale_fill_manual(values = colors[1:5])+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( size = "none", fill = "none", alpha = "none")+
  geom_col(inherit.aes=FALSE, data = dfp, aes(x = cat, y = vall, fill = cat, alpha = 0.8))

ggsave("figure1b_hist1.png", hist1, width = 6, height = 7)


dfp2 = data.frame(value = df1$IVSO_prop)
hist2 = 
  ggplot(dfp2, aes(x = value)) + geom_histogram(binwidth = 0.05, fill = colors[6], color = "black")+
  xlab("\nIntra-population variation in social organization (IVSO)")+
  ylab("Number of populations\n")+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( size = "none", fill = "none", alpha = "none")

library(cowplot)

#combine
f1b = plot_grid(hist1, hist2, align = "v", ncol = 1)
ggsave("figure1b_hist.png", f1b, width = 5, height = 8)


############################################################
#random sample size to capture uncertainty in phylogeny
############################################################

#determine n of random samples from the phylogenetic
#tree to reduce error in estimated probabilities of SO 
#and IVSO; the comparison is made using a naive model 
#for computational efficiency

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 + (1|gr(phylo, cov = A)) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
            1 + (1|gr(phylo, cov = A)) )

#create copy of dataframe for manipulation
df1 = dataGR

#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#Values are 0 (no units observed), not missing or NA
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0

#create response variables for brms
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1)) #observation-level random effect

#estimate random effects model with 1 random phylogeny
m.re1 = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), 
            family = c(multinomial, binomial),
            data = df1, data2 = list(A = A),
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", 
                      paste0(text="prior(\"exponential(2)\", 
                      class = \"sd\", resp=\"SOcounts\", dpar=\"mu", 
                      levels(as.factor(df1$Main1))[-5], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0, 
            control=list(adapt_delta=0.90, max_treedepth=12))
  #save
  saveRDS(m.re1, "m_re1.RDS")
  m.re1 = readRDS("m_re1.RDS")

#estimate random effects model with 5 random phylogenies
  #create data lists
  datal = rep(list(df1),5)
  phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
  phylol = lapply(phylol, ape::vcv, cor = TRUE)
  phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
                list(A = phylol[[3]]), list(A = phylol[[4]]),
                list(A = phylol[[5]]))
  
  m.re5 = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE), 
            family = c(multinomial, binomial),
            data = datal, data2 = phylol,
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", 
                      paste0(text="prior(\"exponential(2)\", 
                      class = \"sd\", resp=\"SOcounts\", dpar=\"mu", 
                      levels(as.factor(df1$Main1))[-5], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0, 
            control=list(adapt_delta=0.90, max_treedepth=12))
  #save
  saveRDS(m.re5, "m_re5.RDS")
  m.re5 = readRDS("m_re5.RDS")

#estimate random effects model with 10 random phylogenies
  #create data lists
  datal = rep(list(df1),10)
  phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),10)
  phylol = lapply(phylol, ape::vcv, cor = TRUE)
  phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
                list(A = phylol[[3]]), list(A = phylol[[4]]),
                list(A = phylol[[5]]), list(A = phylol[[6]]),
                list(A = phylol[[7]]), list(A = phylol[[8]]),
                list(A = phylol[[9]]), list(A = phylol[[10]]))
  
  m.re10 = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE), 
                       family = c(multinomial, binomial),
                       data = datal, data2 = phylol,
                       prior = c(prior("normal(0,1)", class = "Intercept"),
                                 prior("exponential(2)", class = "sd", resp = "IVSOint"),
                                 eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                                 class = \"sd\", resp=\"SOcounts\", dpar=\"mu", 
                                 levels(as.factor(df1$Main1))[-5], "\")",collapse = ","), ")")))), 
                       warmup = n_warm, iter=n_iter, chains = n_chains, init = 0, 
                       control=list(adapt_delta=0.90, max_treedepth=12))
  #save
  saveRDS(m.re10, "m_re10.RDS")
  m.re10 = readRDS("m_re10.RDS")
  
#compare and plot phylogenetic predictions from naive model
  #get predictions
  p.re1 = fitted(m.re1, newdata = data.frame(SO_tot = 1), re_formula = NA, robust = TRUE)
  p.re5 = fitted(m.re5, newdata = data.frame(SO_tot = 1), re_formula = NA, robust = TRUE)
  p.re10 = fitted(m.re10, newdata = data.frame(SO_tot = 1), re_formula = NA, robust = TRUE)

  #estimates across N
  colnames(p.re1)[1:2] = c("median_1","mad_1") #median absolute deviation (robust SD)
  colnames(p.re5)[1:2] = c("median_5","mad_5")
  colnames(p.re10)[1:2] = c("median_10","mad_10")
  Ndf = rbind(p.re1[,1:2,1:6], p.re5[,1:2,1:6], p.re10[,1:2,1:6])
  write.csv(Ndf, "Ndf.csv")
  
  #difference in estimators between N1, N5, and N10 random trees
  Ndf_med = Ndf[c(1,3,5),]
  Ndf_mad = Ndf[c(2,4,6),]
   #N1 -
   round(Ndf_med[1,] - Ndf_med[2,], 3) #N5
   round(Ndf_med[1,] - Ndf_med[3,], 3) #N10
   #N5 -
   round(Ndf_med[2,] - Ndf_med[3,], 3) #N10
   #Absolute delta values remain < 0.01
   
   

##############################################################
#Phylogenetic predictions for superfamilies (no ecological effects)
##############################################################

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 + effort + (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
            1 + effort + (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

#create copy of dataframe for manipulation
df1 = dataGR

#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#prep response and predictor variables
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))
df1$effort =
  as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))

#create data lists
datal = rep(list(df1),5)
phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
phylol = lapply(phylol, ape::vcv, cor = TRUE)
phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
              list(A = phylol[[3]]), list(A = phylol[[4]]),
              list(A = phylol[[5]]))

#estimate model
m.superf = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE), 
          family = c(multinomial, binomial),
          data = datal, data2 = phylol,
          prior = c(prior("normal(0,1)", class = "Intercept"),
                    prior("normal(0,1)", class = "b"),
                    prior("exponential(2)", class = "sd", resp = "IVSOint"),
                    eval(parse(text=paste0("c(", 
                    paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", 
                    levels(as.factor(df1$Main1))[-5], "\")",collapse = ","), ")")))),
          warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
          control=list(adapt_delta=0.99, max_treedepth=10))

#save
saveRDS(m.superf, "m_superf.RDS")
m.superf = readRDS("m_superf.RDS")

#summary of superfamily predictions
names = unique(cbind(df1$Genus_species, df1$superfamily))

pred = fitted(m.superf, newdata = data.frame(phylo = names[,1], superfamily = names[,2], SO_tot = 1),
              re_formula = ~ (1|superfamily) + (1|gr(phylo, cov = A)), robust = FALSE)

pred = fitted(m.superf, newdata = data.frame(phylo = names[,1], Genus_species = names[,1], superfamily = names[,2], SO_tot = 1),
              re_formula = ~ (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species), summary = FALSE)

w.sf = reshape2::melt(pred)
w.sf$sf = rep(names[,2], each = 8000)

library(tidybayes)
superf.p =
ggplot(w.sf, aes(x = value, y = sf, color = sf))+
  facet_wrap(.~ Var3)+
  stat_pointinterval(.width = c(0.9), size = 2.5)+
  stat_pointinterval(w.sf, mapping = aes(x = value, y = sf, group = Var2, color = sf),
                     .width=0, position = position_jitter(), alpha = 0.2, size = 1)+
  scale_color_manual(values = c("#0A2F51","#0E4D64","#137177","#188977","#1D9A6C","#39A96B"))+
  ylab(" ")+
  xlab("Probability of SO for an average social unit")+
  theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( color = "none")

#save
save_plot("superfamily_pred.png", superf.p, base_height = 5.5, base_width = 7.5)

############################################################
#Ancestral state reconstruction (MULTINOMIAL)
############################################################

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 + Activity_pattern + Locom + logmean_bodysize + effort + 
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 1 +
              Activity_pattern + Locom +
              logmean_bodysize + effort + (1|superfamily) +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

#data prep (1 of 10 possible datasets)
df1 = df_10[[1]]

#prepare vcv matrix
A = vcv(c_tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#organize predictors
df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                              ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral", 
                                     df1$Activity_pattern))
df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
df1$Locomotion = as.factor(df1$Locomotion)
df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
df1$logmean_bodysize = 
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                 1, mean, na.rm = TRUE) )))

df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))

#estimate model
m.asr = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
          data = df1, data2 = list(A = A),
          prior = c(prior("normal(0,1)", class = "Intercept"),
                    prior("normal(0,1)", class = "b"), 
                    prior("exponential(2)", class = "sd", resp = "IVSOint"),
                    eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", 
                    levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
          warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
          control=list(adapt_delta=0.99, max_treedepth=12))

#save
saveRDS(m.asr, "m_asr.RDS")
m.asr = readRDS("m_asr.RDS")

#summary of ASR predictions
pred = fitted(m.asr, newdata = data.frame(effort = 0, SO_tot = 1, Nbr_social_units = 1, Activity_pattern = "Nocturnal",
                                        Locom = "AR", logmean_bodysize = -1.5), re_formula = NA, robust = FALSE)
N10_ASRpred = data.frame(pred = names(pred[,"Estimate",]), median = pred[,"Estimate",], sd = pred[,"Est.Error",])
N10_ASRpred[,2:3] = round(N10_ASRpred[,2:3], 2)
N10_ASRpred

pred = data.frame( fitted(m.asr, newdata = data.frame(effort = 0, SO_tot = 1, Nbr_social_units = 1, Activity_pattern = "Nocturnal",
                                                      Locom = "AR", logmean_bodysize = -1.5),
                          re_formula = NA, robust = FALSE, summary=FALSE))

colnames(pred) = c("Solitary", "MF", "MFF", "FMM", "FFMM", "IVSOint")

#probs
pred.df = data.frame(
  median = apply(pred, 2, function(x) median(x)),
  sd = apply(pred, 2, function(x) mad(x)),
  lci = apply(pred, 2, function(x) quantile(x, c(0.05,0.95)))[1,], 
  uci = apply(pred, 2, function(x) quantile(x, c(0.05,0.95)))[2,])
round(pred.df,2)

write.csv(round(pred.df,2),"asr_pred.csv")

#MF versus other primary SOs
diff = data.frame(pred$MF - pred)
diff.df = data.frame(
  median = apply(diff, 2, function(x) median(x)),
  sd = apply(diff, 2, function(x) mad(x)),
  lci = apply(diff, 2, function(x) quantile(x, c(0.05,0.95)))[1,], 
  uci = apply(diff, 2, function(x) quantile(x, c(0.05,0.95)))[2,],
  pp = apply(diff, 2, function(x) sum(x>0)/length(x)))
round(diff.df,2)

write.csv(round(diff.df,2),"asr_diff.csv")

IVSO = 1 - pred$MF
hist(IVSO)
median(IVSO); quantile(IVSO, c(0.25, 0.95))


############################################################
#plot results
pred = data.frame( fitted(m.asr, newdata = 
                          data.frame(effort = 0, SO_tot = 1, 
                                     Activity_pattern = "Nocturnal",
                                     Locom = "AR", logmean_bodysize = -1.5), 
                          re_formula = NA, robust = TRUE,summary = FALSE))

colnames(pred) = c("Solitary", "MF", "MFF", "FMM", "FFMM", "Overall IVSO")

#wide to long
library(tidyr)
lnd1 = tidyr::gather(pred)
lnd1$key = factor(lnd1$key, 
                  levels = c("Solitary","MF", "MFF", "FMM", "FFMM", "Overall IVSO"))

library(ggplot2)
library(tidybayes)

colors = brewer.pal(n = 8, "Dark2")

pred.med = data.frame(med = apply(pred, 2, median), key = factor(colnames(pred), 
                  levels = c("Solitary","MF", "MFF", "FMM", "FFMM", "Overall IVSO")))

asr.primary = 
  ggplot(data = lnd1, aes(x = value, color = key, fill = key, group = key))+
  geom_density(alpha = 0.25, aes(x = value, y = ..scaled..), size = 1)+
  geom_vline(data = pred.med, aes(xintercept = med, group = key), linetype = "dashed")+
  facet_wrap(. ~ key, nrow = 1)+
  scale_x_continuous(limits=c(0,1), expand = c(0.02,0), labels = c("0", "0.25", "0.5", "0.75", "1") )+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(values = colors[1:6]) +
  scale_fill_manual(values = colors[1:6]) +
  xlab("\nProbability for a social unit within an ancestral population")+
  ylab("Relative posterior density\n")+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( size = "none", color = guide_legend(nrow=1), fill = 
            guide_legend(nrow = 1, override.aes=list(color=NA)))


ggsave("asr_p.png", asr.primary, width = 10, height = 3)

#asr.p = plot_grid(asr.primary, asr.ivso, align = "hv", nrow = 1)
#asr.p2 = ggdraw(add_sub(asr.p, "\nProbability in an ancestral population", fontface = "bold", size = 12, 
#                        vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5))
#ggsave("asr_p2.pdf", asr.p2, width = 9, height = 4)
#ggsave("asr_p2.png", asr.p2, width = 9, height = 4)
#ggsave("asr_p1.png", asr.primary, width = 4, height = 2.5)

#saveRDS(asr.p2, "asr_p2.RDS")
#asr.p2 = readRDS("asr_p2.RDS")


############################################################
#Compare ecological inference between datasets (G and R)
############################################################

...

############################################################
#Collective ecological effects (multivariate)
############################################################
#the model is first written in brms
#and then edited in Stan to allow for directly estimating
#variance explained in solitary (the base category)
#using a non-standard but mathematically equivalent
#parameterization with k rather than k-1 linear predictors.

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~ 
            mo(Habitat_heterogenity) + Habitat_cat + 
            foraging_style + Locom + Activity_pattern +
            fruitprop + folivprop + seedprop + animalprop +
            logmean_bodysize + effort + (1|superfamily) + 
            (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

IVSO_m = bf(IVSOint | trials(Nbr_social_units) ~ 1 +
             mo(Habitat_heterogenity) + Habitat_cat + 
             foraging_style + Locom + Activity_pattern + 
             fruitprop + folivprop + seedprop + animalprop +
             logmean_bodysize + effort + (1|superfamily) + 
             (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

#data prep
df1 = df_10[[1]]
#prepare vcv matrix
A = vcv(c_tree, corr = TRUE)
#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]
#organize predictors
df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                              ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral", 
                                     df1$Activity_pattern))
df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
df1$Locomotion = as.factor(df1$Locomotion)
df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
df1$Locom = as.factor(df1$Locom)
df1$logmean_bodysize = 
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                 1, mean, na.rm = TRUE) )))

df1$foraging_style = ifelse(df1$foraging_style == "one", "solitary",
                            ifelse(df1$foraging_style == "unknown", NA,
                                   ifelse(df1$foraging_style == " ", NA,
                                          df1$foraging_style )))
df1$foraging_style = relevel(as.factor(df1$foraging_style), ref = "solitary")
df1$fruitprop = df1$mean_Fruits_all/100
df1$folivprop = df1$mean_Leaves_all/100
df1$flowerprop = df1$mean_Flowers_all/100
df1$seedprop = df1$mean_Seeds_all/100
df1$animalprop = df1$mean_Animal_all/100

df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0

df1[df1$Habitat_cat=="","Habitat_cat"] = NA
df1$Habitat_cat = as.factor(df1$Habitat_cat)

#generated dataset with imputation of missing values
imp = mice(df1[,c("Genus_species","phylo","superfamily","effort",
                  "Solitary","MF","MFF","FMM","FFMM",
                  "Nbr_social_units", "IVSOint","Main1",
                  "Habitat_heterogenity","Habitat_cat",
                  "foraging_style","logmean_bodysize",
                  "Locom","Activity_pattern","fruitprop","folivprop",
                  "flowerprop","seedprop","animalprop")], m = 1)
df1.mi = complete(imp)

df1.mi$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1.mi$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1.mi$MainSO = df1$MainSO
df1.mi$obs = seq(1:nrow(df1))

#estimate model
m.full = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
            data = df1.mi, data2 = list(A = A),
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("normal(0,1)", class = "b"), 
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu",
                                                          levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
            control=list(adapt_delta=0.90, max_treedepth=10))

#save
saveRDS(m.full, "m_full.RDS")
m.full = readRDS("m_full.RDS")
summary(m.full, robust = TRUE)

#without imputation
#estimate model
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$MainSO = df1$MainSO
df1$obs = seq(1:nrow(df1))

m.full_nomi = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
            data = df1, data2 = list(A = A),
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("normal(0,1)", class = "b"), 
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu",
                   levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
            control=list(adapt_delta=0.90, max_treedepth=10))

#save
saveRDS(m.full_nomi, "m_full_nomi.RDS")
m.full_nomi = readRDS("m_full_nomi.RDS")
summary(m.full_nomi)

#reparameterize in Stan

#get model code
txt = make_stancode(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
              data = df1.mi, data2 = list(A = A),
              prior = c(prior("normal(0,1)", class = "Intercept"),
                        prior("normal(0,1)", class = "b"), 
                        prior("exponential(2)", class = "sd", resp = "IVSOint"),
                        eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1.mi$MainSO))[-1], "\")",collapse = ","), ")")))))
write(txt,"brms_stan1.stan")

#load modified code
mod_r = rstan::stan_model("brms_reparam1.stan")

#get data file
data = make_standata(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
              data = df1.mi, data2 = list(A = A),
              prior = c(prior("normal(0,1)", class = "Intercept"),
                        prior("normal(0,1)", class = "b"), 
                        prior("exponential(2)", class = "sd", resp = "IVSOint"),
                        eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))))

#add extra variables for solitary variable
data$N_21 = data$N_1
data$N_22 = data$N_2
data$N_23 = data$N_3
data$N_24 = data$N_4
data$M_21 = data$M_1
data$M_22 = data$M_2
data$M_23 = data$M_3
data$M_24 = data$M_4
data$J_21_SOcounts = data$J_1_SOcounts
data$J_22_SOcounts = data$J_2_SOcounts
data$J_23_SOcounts = data$J_3_SOcounts
data$J_24_SOcounts = data$J_4_SOcounts
data$Z_21_muSolitary_SOcounts_1 = data$Z_1_muMF_SOcounts_1
data$Z_22_muSolitary_SOcounts_1 = data$Z_2_muMF_SOcounts_1
data$Z_23_muSolitary_SOcounts_1 = data$Z_3_muMF_SOcounts_1
data$Z_24_muSolitary_SOcounts_1 = data$Z_4_muMF_SOcounts_1
data$K_muSolitary_SOcounts = data$K_muMF_SOcounts
data$X_muSolitary_SOcounts = data$X_muMF_SOcounts
data$Ksp_muSolitary_SOcounts = data$Ksp_muMF_SOcounts
data$Imo_muSolitary_SOcounts = data$Imo_muMF_SOcounts
data$Jmo_muSolitary_SOcounts = data$Jmo_muMF_SOcounts
data$Xmo_muSolitary_SOcounts_1 = data$Xmo_muMF_SOcounts_1
data$con_simo_muSolitary_SOcounts_1 = data$con_simo_muMF_SOcounts_1
data$Lcov_23 = data$Lcov_3

#run model

#MCMC settings
n_iter <- 3000
n_warm <- 1000
n_chains <- 4

#estimate
full.mod <- rstan::sampling(mod_r, data=data, init=0, iter = n_iter, warmup = n_warm, seed = 9,
                            chains = n_chains, cores = n_chains, control = list(adapt_delta=0.90, max_treedepth=10))
saveRDS(full.mod, "full_mod.RDS")
full.mod = readRDS("full_mod.RDS")


launch_shinystan(full.mod)

############################################################
#Plot variance explained
############################################################
post = (rstan::extract(full.mod))

#variance explained by fixed effects
fe_v =  data.frame(Solitary = apply(post$muSolitary_SOcounts, 1, var) - apply(post$muSolitary_SOcounts2, 1, var),
                  MF = apply(post$muMF_SOcounts, 1, var) - apply(post$muMF_SOcounts2, 1, var),
                  MFF = apply(post$muMFF_SOcounts, 1, var) - apply(post$muMFF_SOcounts2, 1, var),
                  FMM = apply(post$muFMM_SOcounts, 1, var) - apply(post$muFMM_SOcounts2, 1, var),
                  FFMM = apply(post$muFFMM_SOcounts, 1, var) - apply(post$muFFMM_SOcounts2, 1, var),
                  IVSO = apply(post$mu_IVSOint, 1, var) - apply(post$mu_IVSOint2, 1, var))
#change negative values
fe_v[fe_v<0] = 0


species_v =  data.frame(Solitary = post$sd_21^2,
                        MF = post$sd_1^2,
                      MFF = post$sd_5^2,
                      FMM = post$sd_9^2,
                      FFMM = post$sd_13^2,
                      IVSO = post$sd_17^2)

phylo_v =  data.frame(Solitary = post$sd_23^2,
                        MF = post$sd_3^2,
                        MFF = post$sd_7^2,
                        FMM = post$sd_11^2,
                        FFMM = post$sd_15^2,
                        IVSO = post$sd_19^2)

superfamily_v =  data.frame(Solitary = post$sd_24^2,
                      MF = post$sd_4^2,
                      MFF = post$sd_8^2,
                      FMM = post$sd_12^2,
                      FFMM = post$sd_16^2,
                      IVSO = post$sd_20^2)

pop_v =  data.frame(Solitary = post$sd_22^2,
                            MF = post$sd_2^2,
                            MFF = post$sd_6^2,
                            FMM = post$sd_10^2,
                            FFMM = post$sd_14^2,
                            IVSO = post$sd_18^2)

#latent scale variance
logit_v = (pi/sqrt(3))^2

#calculate
phylo_r2 = phylo_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)

#summarize
apply(phylo_r2, 2, mean); apply(phylo_r2, 2, sd)

fe_r2 = fe_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(fe_r2, 2, mean); apply(fe_r2, 2, sd)

species_r2 = species_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(species_r2, 2, mean); apply(species_r2, 2, sd)

superf_r2 = superfamily_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(superf_r2, 2, mean); apply(superf_r2, 2, sd)

pop_r2 = logit_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(pop_r2, 2, mean); apply(pop_r2, 2, sd)

###########################################################################
#plot

#wide to long
library(tidyr)
ldf_phylo = tidyr::gather(phylo_r2)
ldf_fe = tidyr::gather(fe_r2)
ldf_species = tidyr::gather(species_r2)
ldf_superf = tidyr::gather(superf_r2)
ldf_pop = tidyr::gather(pop_r2)
ldf_phylo$type = rep("Phylogeny",nrow(ldf_phylo))
ldf_fe$type = rep("Ecological predictors",nrow(ldf_fe))
ldf_species$type = rep("Species residual",nrow(ldf_species))
ldf_superf$type = rep("Superfamily residual",nrow(ldf_superf))
ldf_pop$type = rep("Population residual",nrow(ldf_pop))

ldf = rbind(ldf_phylo, ldf_fe, ldf_species, ldf_superf, ldf_pop)
ldf$key[ldf$key=="IVSO"] = "Overall IVSO"
ldf$key = factor(ldf$key, levels = c("Solitary", "MF","MFF" , "FMM", "FFMM", "Overall IVSO"))
ldf$type = factor(ldf$type, levels = c("Phylogeny", "Ecological predictors",  "Population residual",
                                       "Species residual", "Superfamily residual"))
library(ggplot2)
library(tidybayes)

cc <- scales::seq_gradient_pal("lightblue", "darkblue")(seq(0,1,length.out=5))

R2.plot = 
  ggplot(ldf, aes(x = value, y = type, color = type, fill = type))+ 
  coord_cartesian(xlim=c(0,1))+
  stat_pointinterval(.width = c(0.90), size = 10, position = position_dodge(width = 0.5)) +
  facet_wrap(. ~ key, nrow=1)+
  scale_color_manual(values = cc)+
  scale_fill_manual(values = cc)+
  scale_y_discrete(limits=rev(levels(ldf$type)))+
  scale_x_continuous(expand=c(0,0),breaks=c(0,0.5,1), labels =c("0","0.5","1"))+
  xlab(" \nProportion of variance explained")+
  theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=10, face = "bold"),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.spacing = unit(1.3, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top",
        plot.margin = unit(c(0.1,0.4,0.7,0.05), "cm"))

#save
save_plot("variance_explained.png", R2.plot, base_height = 3.5, base_width = 12)
saveRDS(R2.plot, "r2_plot.RDS")
R2.plot = readRDS("r2_plot.RDS")


############################################################
#Specific ecological effects (univariate)
############################################################
#each section fits a model for a single predictor
#and then creates a plot for the total effect

df1 = df_10[[1]]
#prepare vcv matrix
A = vcv(c_tree, corr = TRUE)
#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]
#organize predictors
df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                              ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral", 
                                     df1$Activity_pattern))
df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
df1$Locomotion = as.factor(df1$Locomotion)
df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
df1$Locom = as.factor(df1$Locom)
df1$logmean_bodysize = 
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                 1, mean, na.rm = TRUE) )))

df1$foraging_style = ifelse(df1$foraging_style == "one", "solitary",
                            ifelse(df1$foraging_style == "unknown", NA,
                                   ifelse(df1$foraging_style == " ", NA,
                                          df1$foraging_style )))
df1$foraging_style = relevel(as.factor(df1$foraging_style), ref = "solitary")
df1$fruitprop = df1$mean_Fruits_all/100
df1$folivprop = df1$mean_Leaves_all/100
df1$flowerprop = df1$mean_Flowers_all/100
df1$seedprop = df1$mean_Seeds_all/100
df1$animalprop = df1$mean_Animal_all/100

df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0

df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))

############################################################
#activity
#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 + Activity_pattern + effort +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
            1 + Activity_pattern + effort + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.act = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
             data = df1, data2 = list(A = A),
             prior = c(prior("normal(0,1)", class = "Intercept"),
                       prior("normal(0,1)", class = "b"), 
                       prior("exponential(2)", class = "sd", resp = "IVSOint"),
                       eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
             warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
             control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.act, "m_act.RDS")
m.act = readRDS("m_act.RDS")
summary(m.act, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.act, resp = "SOcounts", effect = "Activity_pattern",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_Activity_pattern:cats__`$data

c2  = plot(conditional_effects(m.act, resp = "IVSOint", effect = "Activity_pattern",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_Activity_pattern$data

df2$cats__ = "Overall IVSO"
df2$effect1__ = df2$Activity_pattern
df2$effect2__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

df3$response = factor(df3$effect1__, levels = c("Solitary", "MF", "MFF", "FMM", "FFMM","Overall IVSO"))

act.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df3$lower__, ymax = df3$upper__, width = 0,
                position=position_dodge(width=0.5))+
  xlab("\nActivity pattern")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

#IVSO effect
#c2  = plot(conditional_effects(m.act, resp = "IVSOint", effect = "Activity_pattern",
#                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
#df = c2$IVSOint.IVSOint_Activity_pattern$data
#
#df$response = "Overall IVSO"
#df$response = factor(df$response, levels = c("Overall IVSO", "Solitary", "MF", "MFF", "FMM", "FFMM"))
#
#act.p2 = 
#  ggplot(df, aes(x = effect1__, y = estimate__, group = effect1__, color = response))+
#  scale_y_continuous(lim=c(0,1))+
#  geom_point(position=position_dodge(width=0.5), size = 2)+
#  geom_errorbar(ymin = df$lower__, ymax = df$upper__, width = 0,
#                position=position_dodge(width=0.5))+
#  scale_color_manual(values = colors[6], drop = TRUE)+
#  xlab("\nActivity pattern")+
#  ylab("Probability\n")+
#  theme(legend.position = "top",
#        legend.title = element_blank(),
#        legend.key = element_rect(fill = "white"),
#        legend.text = element_text(size = 8),
#        axis.title.y = element_blank(),
#        axis.text.y = element_text(size = 8),
#        axis.text.x = element_text(size = 8),
#        axis.title.x = element_blank(),
#        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
#        panel.background= element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.spacing = unit(2, "lines"))
#combine
#act.p3 = plot_grid(act.p1, act.p2, align = "hv", nrow = 1)
#act.p3 = ggdraw(add_sub(act.p3, "\nActivity pattern", fontface = "bold", size = 10, 
#                        vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5))
#ggsave("act_p3.png", act.p3, width = 9, height = 3)

#calculate comparisons of interest
newdf = data.frame(Activity_pattern = unique(df1$Activity_pattern), SO_tot = 1, effort = 0)
pred = fitted(N10_cmact, newdata = newdf, re_formula = NA, summary = FALSE, scale = "response")
comp1 = apply(pred, 3, FUN = function(x) x[,1] - x[,2] ) #nocturnal - diurnal
comp2 = apply(pred, 3, FUN = function(x) x[,3] - x[,2] ) #cathemeral - diurnal
comp3 = apply(pred, 3, FUN = function(x) x[,1] - x[,3] ) #nocturnal - diurnal

act_comp = 
  data.frame(
    outcome = names(pred[1,1,]),
    diff_noc_diu = paste0(round(apply(comp1,2,median),2), "(",
                          apply(apply(comp1, 2, quantile, c(0.05, 0.95)),2, 
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_cat_diu = paste0(round(apply(comp2,2,median),2), "(",
                          apply(apply(comp2, 2, quantile, c(0.05, 0.95)),2, 
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_noc_cat = paste0(round(apply(comp3,2,median),2), "(",
                          apply(apply(comp3, 2, quantile, c(0.05, 0.95)),2, 
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")") )

write.csv(act_comp,"act_comp.csv")

############################################################
#locomotion
SO_m = bf(SO_counts | trials(SO_tot) ~ 
              1 + Locom + effort +
             (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
              1 + Locom + effort +
             (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.loc = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
            data = df1, data2 = list(A = A),
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("normal(0,1)", class = "b"), 
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
            control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.loc, "m_loc.RDS")
m.loc = readRDS("m_loc.RDS")
summary(m.loc, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.loc, resp = "SOcounts", effect = "Locom",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_Locom:cats__`$data

c2  = plot(conditional_effects(m.loc, resp = "IVSOint", effect = "Locom",
                               conditions = data.frame(Nbr_social_units=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_Locom$data

colors = brewer.pal(n = 8, "Dark2")

df2$cats__ = "Overall IVSO"
df2$effect1__ = df2$Locom
df2$effect2__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

df3$response = factor(df3$effect1__, levels = c("Solitary", "MF", "MFF", "FMM", "FFMM","Overall IVSO"))


loc.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_x_discrete(labels = c("Arboreal", "Both", "Terrestrial"))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df3$lower__, ymax = df3$upper__, width = 0,
                position=position_dodge(width=0.5))+
  xlab("Locomotion")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

#IVSO effect
#c2  = plot(conditional_effects(m.loc, resp = "IVSOint", effect = "Locom",
#                               conditions = data.frame(Nbr_social_units=1), robust = TRUE, prob = 0.9), plot = FALSE)
#df = c2$IVSOint.IVSOint_Locom$data
#
#df$response = "Overall IVSO"
#df$response = factor(df$response, levels = c("Overall IVSO", "Solitary", "MF", "MFF", "FMM", "FFMM"))
#
#loc.p2 = 
#  ggplot(df, aes(x = effect1__, y = estimate__, group = effect1__, color = response))+
#  scale_y_continuous(lim=c(0,1))+
#  scale_x_discrete(labels = c("Arboreal", "Both", "Terrestrial"))+
#  geom_point(position=position_dodge(width=0.5), size = 2)+
#  geom_errorbar(ymin = df$lower__, ymax = df$upper__, width = 0,
#                position=position_dodge(width=0.5))+
#  scale_color_manual(values = colors[6], drop = TRUE)+
#  xlab("Locomotion")+
#  theme(legend.position = "top",
#        legend.title = element_blank(),
#        legend.key = element_rect(fill = "white"),
#        legend.text = element_text(size = 8),
#        axis.title.y = element_blank(),
#        axis.text.y = element_text(size = 8),
#        axis.text.x = element_text(size = 8),
#        axis.title.x = element_blank(),
#        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
#        panel.background= element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.spacing = unit(2, "lines"))

#combine
#loc.p3 = plot_grid(loc.p1, loc.p2, align = "hv", nrow = 1)
#loc.p3 = ggdraw(add_sub(loc.p3, "\nLocomotion", fontface = "bold", size = 10, 
#                        vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5))
#ggsave("loc_p3.png", loc.p3, width = 9, height = 3)

#calculate comparisons of interest
newdf = data.frame(Locom_s = unique(df1$Locom_s), Nbr_social_units = 1)
pred = fitted(N10_cmloc, newdata = newdf, re_formula = NA, summary = FALSE, scale = "response")
comp1 = apply(pred, 3, FUN = function(x) x[,1] - x[,2] ) #AR - BOTH
comp2 = apply(pred, 3, FUN = function(x) x[,3] - x[,4] ) #BOTH - T
comp3 = apply(pred, 3, FUN = function(x) x[,1] - x[,4] ) #AR - T

loc_comp = 
  data.frame(
    outcome = names(pred[1,1,]),
    diff_AR_BOTH = paste0(round(apply(comp1,2,median),2), "(",
                          apply(apply(comp1, 2, quantile, c(0.05, 0.95)),2, 
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_BOTH_T = paste0(round(apply(comp2,2,median),2), "(",
                         apply(apply(comp2, 2, quantile, c(0.05, 0.95)),2, 
                               FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                         ")"),
    diff_AR_T = paste0(round(apply(comp3,2,median),2), "(",
                       apply(apply(comp3, 2, quantile, c(0.05, 0.95)),2, 
                             FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                       ")") )

write.csv(loc_comp,"loc_comp.csv")

############################################################
#body size
SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 + logmean_bodysize + effort +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
              1 + logmean_bodysize + effort +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.bs = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
            data = df1, data2 = list(A = A),
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("normal(0,1)", class = "b"), 
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
            control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.bs, "m_bs.RDS")
m.bs = readRDS("m_bs.RDS")
summary(m.bs, robust = TRUE)

############################################################
#MainSO
c1  = plot(conditional_effects(m.bs, resp = "SOcounts", effect = "logmean_bodysize",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_logmean_bodysize:cats__`$data

c2  = plot(conditional_effects(m.bs, resp = "IVSOint", effect = "logmean_bodysize",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_logmean_bodysize$data

#df$response = "x"
df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

bs.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.1)+
  xlab("\n Log mean body size")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

ggsave("bs_p1.png", bs.p1, width = 10, height = 3)


############################################################
#combine plots for Figure 3
library(cowplot)

#combine
f2a = R2.plot
f2b = plot_grid(act.p1,loc.p1, align = "h", nrow = 1)
f2b = plot_grid(f2b, bs.p1, nrow = 2)
f2c = asr.primary
figure2 = plot_grid(f2a, f2b, f2c, ncol = 1, rel_heights = c(0.25,0.5,0.25))
ggsave("figure2_primateivso.png", figure2, width = 10.5, height = 11.5)


############################################################
#habitat effects
df1.1 = df1[df1$Habitat_cat!="",]

SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 +  mo(Habitat_heterogenity) + Habitat_cat + effort +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
              1 + mo(Habitat_heterogenity) + Habitat_cat + effort +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.hab = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
           data = df1.1, data2 = list(A = A),
           prior = c(prior("normal(0,1)", class = "Intercept"),
                     prior("normal(0,1)", class = "b"), 
                     prior("exponential(2)", class = "sd", resp = "IVSOint"),
                     eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
           warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
           control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.hab, "m_hab.RDS")
m.hab = readRDS("m_hab.RDS")
summary(m.hab, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.hab, resp = "SOcounts", effect = "Habitat_cat",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_Habitat_cat:cats__`$data

colors = brewer.pal(n = 8, "Dark2")

hab.p1 = 
  ggplot(df, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df$lower__, ymax = df$upper__, width = 0,
                position=position_dodge(width=0.5))+
  xlab("\n Habitat type")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:5])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

#IVSO effect
c2  = plot(conditional_effects(m.hab, resp = "IVSOint", effect = "Habitat_cat",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df = c2$IVSOint.IVSOint_Habitat_cat$data

df$response = "Overall IVSO"
df$response = factor(df$response, levels = c("Overall IVSO", "Solitary", "MF", "MFF", "FMM", "FFMM"))

hab.p2 = 
  ggplot(df, aes(x = effect1__, y = estimate__, group = effect1__, color = response))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df$lower__, ymax = df$upper__, width = 0,
                position=position_dodge(width=0.5))+
  scale_color_manual(values = colors[6], drop = TRUE)+
  xlab("\nHabitat type")+
  ylab("Probability\n")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))
#combine
act.p3 = plot_grid(hab.p1, hab.p2, align = "hv", nrow = 1)
act.p3 = ggdraw(add_sub(act.p3, "\nHabitat type", fontface = "bold", size = 10, 
                        vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5))
ggsave("hab_p3.png", act.p3, width = 9, height = 3)

c1  = plot(conditional_effects(m.hab, resp = "SOcounts", effect = "Habitat_heterogenity",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_Habitat_heterogenity:cats__`$data

c2  = plot(conditional_effects(m.hab, resp = "IVSOint", effect = "Habitat_heterogenity",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_Habitat_heterogenity$data

df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

het.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.1)+
  xlab("\n Habitat heterogeneity")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

ggsave("het_p1.png", het.p1, width = 10, height = 3)

############################################################
#fruit diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 +  fruitprop + effort + 
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
            1 +  fruitprop + effort + 
            (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.fruit = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
            data = df1, data2 = list(A = A),
            prior = c(prior("normal(0,1)", class = "Intercept"),
                      prior("normal(0,1)", class = "b"), 
                      prior("exponential(2)", class = "sd", resp = "IVSOint"),
                      eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
            warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
            control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.fruit, "m_fruit.RDS")
m.fruit = readRDS("m_fruit.RDS")
summary(m.fruit, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.fruit, resp = "SOcounts", effect = "fruitprop",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_fruitprop:cats__`$data

c2  = plot(conditional_effects(m.fruit, resp = "IVSOint", effect = "fruitprop",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_fruitprop$data

df$response = "x"
df2$response = "x"
df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

fruit.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.1)+
  xlab("\n Proportion of fruit (diet)")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

ggsave("fruit_p1.png", fruit.p1, width = 10, height = 3)


#calculate comparisons of interest
newdf = data.frame(fruitprop = c(0,1), SO_tot = 1, effort = 0)
pred = fitted(m.fruit, newdata = newdf, re_formula = NA, summary = FALSE, scale = "response")
comp1 = apply(pred, 3, FUN = function(x) x[,2] - x[,1] ) #1 - 0

median(comp1[,"IVSOint"]); sd(comp1[,"IVSOint"]); quantile(comp1[,"IVSOint"], c(0.05, 0.95))

bs_comp = 
  data.frame(
    outcome = names(pred[1,1,]),
    diff_low_avg = paste0(round(apply(comp1,2,median),2), "(",
                          apply(apply(comp1, 2, quantile, c(0.05, 0.95)),2, 
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_avg_high = paste0(round(apply(comp2,2,median),2), "(",
                           apply(apply(comp2, 2, quantile, c(0.05, 0.95)),2, 
                                 FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                           ")"),
    diff_low_high = paste0(round(apply(comp3,2,median),2), "(",
                           apply(apply(comp3, 2, quantile, c(0.05, 0.95)),2, 
                                 FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                           ")") )

write.csv(bs_comp,"bs_comp.csv")



############################################################
#folivory diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 +  folivprop + effort + 
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
              1 +  folivprop + effort + 
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.foliv = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
              data = df1, data2 = list(A = A),
              prior = c(prior("normal(0,1)", class = "Intercept"),
                        prior("normal(0,1)", class = "b"), 
                        prior("exponential(2)", class = "sd", resp = "IVSOint"),
                        eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
              warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
              control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.foliv, "m_foliv.RDS")
m.foliv = readRDS("m_foliv.RDS")
summary(m.foliv, robust = TRUE)
############################################################
c1  = plot(conditional_effects(m.foliv, resp = "SOcounts", effect = "folivprop",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_folivprop:cats__`$data

c2  = plot(conditional_effects(m.foliv, resp = "IVSOint", effect = "folivprop",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_folivprop$data

df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

folv.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.1)+
  xlab("\n Proportion of foliage (diet)")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

ggsave("folv_p1.png", folv.p1, width = 10, height = 3)

############################################################
#seed diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~ 
              1 +  seedprop + effort + 
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
              1 +  seedprop + effort + 
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.seed = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
              data = df1, data2 = list(A = A),
              prior = c(prior("normal(0,1)", class = "Intercept"),
                        prior("normal(0,1)", class = "b"), 
                        prior("exponential(2)", class = "sd", resp = "IVSOint"),
                        eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
              warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
              control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.seed, "m_seed.RDS")
m.seed = readRDS("m_seed.RDS")
summary(m.seed, robust = TRUE)
############################################################
c1  = plot(conditional_effects(m.seed, resp = "SOcounts", effect = "seedprop",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_seedprop:cats__`$data

c2  = plot(conditional_effects(m.seed, resp = "IVSOint", effect = "seedprop",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_seedprop$data

df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

seed.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.1)+
  xlab("\n Proportion of seeds (diet)")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

ggsave("seed_p1.png", seed.p1, width = 10, height = 3)

############################################################
#animal diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~ 
            1 +  animalprop + effort + 
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
              1 +  animalprop + effort + 
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.animal = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
             data = df1, data2 = list(A = A),
             prior = c(prior("normal(0,1)", class = "Intercept"),
                       prior("normal(0,1)", class = "b"), 
                       prior("exponential(2)", class = "sd", resp = "IVSOint"),
                       eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
                    class = \"sd\", resp=\"SOcounts\", dpar=\"mu", levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
             warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
             control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.animal, "m_animal.RDS")
m.animal = readRDS("m_animal.RDS")
summary(m.animal, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.animal, resp = "SOcounts", effect = "animalprop",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_animalprop:cats__`$data

c2  = plot(conditional_effects(m.animal, resp = "IVSOint", effect = "animalprop",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_animalprop$data

df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")

animal.p1 = 
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.1)+
  xlab("\n Proportion of animal (diet)")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

ggsave("anim_p1.png", animal.p1, width = 10, height = 3)

############################################################
#foraging style
SO_m = bf(SO_counts | trials(SO_tot) ~ 
    1 + foraging_style + effort +
    (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~ 
    1 + foraging_style + effort +
    (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.forg = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
              data = df1, data2 = list(A = A),
              prior = c(prior("normal(0,1)", class = "Intercept"),
              prior("normal(0,1)", class = "b"), 
              prior("exponential(2)", class = "sd", resp = "IVSOint"),
              eval(parse(text=paste0("c(", paste0(text="prior(\"exponential(2)\", 
              class = \"sd\", resp=\"SOcounts\", dpar=\"mu", 
              levels(droplevels(df1$MainSO))[-1], "\")",collapse = ","), ")")))),
              warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0, 
              control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.forg, "m_forg.RDS")
m.forg = readRDS("m_forg.RDS")
summary(m.forg, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.forg, resp = "SOcounts", effect = "foraging_style",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_foraging_style:cats__`$data

colors = brewer.pal(n = 8, "Dark2")

forg.p1 = 
  ggplot(df, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df$lower__, ymax = df$upper__, width = 0,
                position=position_dodge(width=0.5))+
  xlab("\nActivity pattern")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:5])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

#IVSO effect
c2  = plot(conditional_effects(m.forg, resp = "IVSOint", effect = "foraging_style",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df = c2$IVSOint.IVSOint_foraging_style$data

df$response = "Overall IVSO"
df$response = factor(df$response, levels = c("Overall IVSO", "Solitary", "MF", "MFF", "FMM", "FFMM"))

forg.p2 = 
  ggplot(df, aes(x = effect1__, y = estimate__, group = effect1__, color = response))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df$lower__, ymax = df$upper__, width = 0,
                position=position_dodge(width=0.5))+
  scale_color_manual(values = colors[6], drop = TRUE)+
  #scale_color_discrete(drop = FALSE)+
  xlab("\nForaging Style")+
  ylab("Probability\n")+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))
#combine
forg.p3 = plot_grid(forg.p1, forg.p2, align = "hv", nrow = 1)
forg.p3 = ggdraw(add_sub(forg.p3, "\nForaging strategy", fontface = "bold", size = 10, 
                        vpadding=grid::unit(0,"lines"),y=5, x=0.55, vjust=4.5))
ggsave("forg_p3.png", forg.p3, width = 9, height = 3)




#please email jordanscottmartin@gmail.com if you have any
#questions or comments

