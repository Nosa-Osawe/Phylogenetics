library(pheatmap)
library(V.PhyloMaker)
library(ape)
library(tidyverse)
library(readxl)
library(vegan)
library(phytools)
library(emmeans)
library(nlme)
library(lme4)
library(lmerTest)


plant <- read_excel("C:\\Users\\DELL\\Documents\\Git in R\\Phylogenetics\\Data\\Herbicide_Phylogenetic.xlsx", 
                    sheet = "Sheet1")
view(plant)
sum(is.na(plant$Q11))

plant$Q11 <- as.numeric(plant$Q11)

plant.list <- plant %>% 
  filter(Week==0, Belt=="B1") %>% 
  select(Scientific_name, Genus,Family ) %>% 
  rename(species = Scientific_name,
         genus = Genus, 
         family = Family) %>% 
  as.data.frame()


plant_tree <-  phylo.maker(plant.list, output.tree = TRUE, output.sp.list = TRUE)

plot(plant_tree$scenario.3, show.tip.label = TRUE)

# Compute the phylogenetic distance matrix
distance <- cophenetic(plant_tree$scenario.3)

# heatmap of phylogenetic distances
pheatmap(distance, 
         color = colorRampPalette(c("white", "red"))(50),
         main = "Phylogenetic Distance Matrix")

# Pre-spray ground cover, collapsed (averaged) by row [Q]--making 17 surveys

pre_spay <- plant %>% 
  filter(Spray=="Pre-spray") %>% 
  pivot_longer(
    cols = -c(1:9),
    names_to = "Q",
    values_to = "Rows"
  ) %>% 
  as.data.frame() %>%
  group_by(Q, Scientific_name) %>% 
  summarise(mean_ground = mean(Rows)) %>% 
  as.data.frame()

species <-  pre_spay %>% 
  filter(Q=="Q1") %>% 
  select(Scientific_name) %>% 
  as.data.frame()
species<-species$Scientific_name


# Q!
q1_cover <- pre_spay %>% 
  filter(Q=="Q1") %>% 
  select(mean_ground) %>% 
  as.data.frame()

q1_cover<- as.vector(t(q1_cover))
names(q1_cover) <- species 

# Blomberg’s K for phylogenetic signal
q1.signal <- phylosig(plant_tree$scenario.3, q1_cover, method = "K")
print(q1.signal)

# Q2
q2_cover <- pre_spay %>% 
  filter(Q=="Q2") %>% 
  select(mean_ground) %>% 
  as.data.frame()

q2_cover<- as.vector(t(q2_cover))
names(q2_cover) <- species 

# Blomberg’s K for phylogenetic signal
q2.signal <- phylosig(plant_tree$scenario.3, q2_cover, method = "K")
print(q2.signal)

######################################################################################

# Based on the methodology, we have to do away with some rows.

G.cover <- plant %>% 
  pivot_longer(
    cols = -c(1:9),
    names_to = "Q",
    values_to = "Rows"
  ) %>% 
  filter(!Q %in% c("Q3", "Q6", "Q9", "Q12", "Q15")) %>%
  select(Week, Scientific_name, Q, Rows) %>% 
  mutate(Week = factor(Week)) %>% 
  mutate(Rep = case_when(
    Q %in% c("Q1", "Q2")  ~ "R1",
    Q %in% c("Q4", "Q5")  ~ "R2",
    Q %in% c("Q7", "Q8")  ~ "R3",
    Q %in% c("Q10", "Q11")  ~ "R4",
    Q %in% c("Q13", "Q14") ~ "R5",
    Q %in% c("Q16", "Q17") ~ "R6"
  )) %>% 
  group_by(Week, Rep, Scientific_name) %>% 
  summarise(mean_ground = mean(Rows)) %>% 
  as.data.frame()

# Just to get the species name-- nothing more
species <-  G.cover %>% 
  filter(Rep=="R1", Week=="0") %>% 
  select(Scientific_name) %>% 
  as.data.frame()
species<-species$Scientific_name


# The first replicate: R1
R1_cover <- G.cover %>% 
  filter(Rep=="R1", Week=="0") %>% 
  select(mean_ground) %>% 
  as.data.frame()

R1_cover<- as.vector(t(R1_cover))
names(R1_cover) <- species 

# Blomberg’s K for phylogenetic signal
R1.signal <- phylosig(plant_tree$scenario.3, R1_cover, 
                      method = "K", nsim = 999)
print(R1.signal)


# Create a function to make things work faster and avoid too much repetition:

Phylo_sig <- function(Rep, Week) {
  R_cover <- G.cover %>% 
    filter(Rep == !!Rep, Week == !!Week) %>% # where Rep matches the input [Rep]
    select(mean_ground) %>% 
    as.data.frame()
  
  R_cover <- as.vector(t(R_cover))
  names(R_cover) <- species 
  
  R_signal <- phylosig(plant_tree$scenario.3, R_cover, 
                       method = "K", nsim = 999)

  print(R_signal)}


# Pre-spray: Week = 0

Phylo_sig(Rep = "R1", Week = "0")
Phylo_sig(Rep = "R2", Week = "0")
Phylo_sig(Rep = "R3", Week = "0")
Phylo_sig(Rep = "R4", Week = "0")
Phylo_sig(Rep = "R5", Week = "0")
Phylo_sig(Rep = "R6", Week = "0")


# Week = 3
Phylo_sig(Rep = "R1", Week = "3")
Phylo_sig(Rep = "R2", Week = "3")
Phylo_sig(Rep = "R3", Week = "3")
Phylo_sig(Rep = "R4", Week = "3")
Phylo_sig(Rep = "R5", Week = "3")
Phylo_sig(Rep = "R6", Week = "3")


# Week = 7
Phylo_sig(Rep = "R1", Week = "7")
Phylo_sig(Rep = "R2", Week = "7")
Phylo_sig(Rep = "R3", Week = "7")
Phylo_sig(Rep = "R4", Week = "7")
Phylo_sig(Rep = "R5", Week = "7")
Phylo_sig(Rep = "R6", Week = "7")


# Week = 9
Phylo_sig(Rep = "R1", Week = "9")
Phylo_sig(Rep = "R2", Week = "9")
Phylo_sig(Rep = "R3", Week = "9")
Phylo_sig(Rep = "R4", Week = "9")
Phylo_sig(Rep = "R5", Week = "9")
Phylo_sig(Rep = "R6", Week = "9")


# Week = 11
Phylo_sig(Rep = "R1", Week = "11")
Phylo_sig(Rep = "R2", Week = "11")
Phylo_sig(Rep = "R3", Week = "11")
Phylo_sig(Rep = "R4", Week = "11")
Phylo_sig(Rep = "R5", Week = "11")
Phylo_sig(Rep = "R6", Week = "11")

# So, I have to manually fetch all # Blomberg’s K  value,
# because I cant find a R script to build it up.
# I'd bring back the result from Excel file!! Thanks for your understanding!





signals <- read_excel("C:\\Users\\DELL\\Documents\\Git in R\\Phylogenetics\\Data\\Herbicide_Phylogenetic.xlsx", 
                    sheet = "signals")
view(signals)

colour_choice <- c("orange", "red", "black", "brown", "darkgreen","purple")

signals %>% 
 # mutate(Week= factor(Week)) %>% 
  pivot_longer(cols = -1,
               names_to = "Replicates",
               values_to = "Signals"
               ) %>% 
  mutate(Week= as.numeric(Week),
         Replicates= factor(Replicates)) %>% 
  ggplot(aes(y= Signals, x= Week, colour = Replicates))+
  geom_point(alpha = 0.7)+
  geom_line(linewidth = 0.7)+
  scale_fill_manual(values = colour_choice)+
  scale_colour_manual(values = colour_choice)+
  scale_x_continuous(breaks = c(0,3,7,9,11)) +
  theme_classic()+ 
  labs( y = "Phylogenetic signal K")

signal_L <- signals %>% 
  pivot_longer(cols = -1,
               names_to = "Replicates",
               values_to = "Signals"
  ) %>% 
  mutate(Week= as.numeric(Week),
         Replicates= factor(Replicates)) %>% 
  as.data.frame()

lm_model <- lm(Signals ~ Week, data = signal_L)
summary(lm_model)



lmer_model <- lmer(Signals ~ Week + (1 | Replicates), data = signal_L)
summary(lmer_model)


summary(lmer_model)

anova(lmer_model,lm_model)


lme_model <-lme(
  Signals ~ Week,                      
  random = ~ 1 | Replicates,           
  correlation = corAR1(form = ~ Week | Replicates),  
  data = signal_L
)
summary(lme_model) 


cover_L <- G.cover %>% 
  pivot_wider(
    names_from = "Scientific_name",
    values_from = mean_ground
  ) %>% 
  group_by(Week, Rep) %>% 
  summarise(total= sum(across(where(is.numeric)))) %>%  
  as.data.frame() 

signal.cover <- cbind(signal_L, cover_L$total)
signal.cover <- signal.cover %>% 
  rename(Cover ="cover_L$total")

# Relationship between plant cover and phylogenetic consevatin signal

lm_mod <- lm(Signals ~ Cover, data = signal.cover)
summary(lm_mod)

poly_mod <- lm(Signals ~ poly(Cover, 3), data = signal.cover)
summary(poly_mod)


 
lm_mixed_nlme <- lme(Signals ~ Cover, random = ~ 1 | Replicates, 
                       data = signal.cover)
summary(lm_mixed_nlme)

poly_mixed_nlme <- lme(Signals ~ poly(Cover, 2), random = ~ 1 | Replicates, 
                       data = signal.cover)
summary(poly_mixed_nlme)


AIC(lm_mixed_nlme,poly_mixed_nlme )

signal.cover %>% 
 ggplot(aes(y= Signals, x= Cover, fill = Replicates, colour = Replicates))+
  geom_point(aes(y= Signals, x= Cover), size= 2, alpha = 0.7)+
  scale_fill_manual(values = colour_choice)+
  scale_colour_manual(values = colour_choice)+
  geom_line(aes(y = fitted(poly_mixed_nlme)), 
            linewidth = 1, alpha= 0.5)+
  scale_x_continuous(breaks = c(0,25,50,75,100)) +
  labs(y = "Phylogenetic signal",
       x = "Plant cover (%)")+
  theme_light()

### Test for difference in Phylogenetic conservation between pre-spay and post-spray

week0 <- signal.cover %>% 
  filter(Week=="0")

Week3 <- signal.cover %>% 
  filter(Week=="3")

t.test(week0$Signals, Week3$Signals)
t.test(week0$Cover, Week3$Cover)


#############################################################################

# Load libraries
library(ape)        # For phylogenetic tree manipulation
library(phytools)   # For phylogenetic signal (Blomberg's K)
library(nlme)       # For PGLS regression

install.packages("caper")
library(caper)      # For comparative phylogenetic analyses
