
# Based on the methodology, we have to do away with some rows.
groundcover<- plant %>% 
  pivot_longer(
    cols = -c(1:9),
    names_to = "Q",
    values_to = "Rows"
  ) %>% 
  filter(!Q %in% c("Q3", "Q4", "Q7", "Q8", "Q11","Q12", "Q15","Q16","Q17")) %>%
  select(Week, Scientific_name, Q, Rows, Belt) %>% 
  mutate(Week = factor(Week)) %>% 
  mutate(Rep = case_when(
    Q %in% c("Q1", "Q2")  ~ "R1",
    Q %in% c("Q5", "Q6")  ~ "R2",
    Q %in% c("Q9", "Q10")  ~ "R3",
    Q %in% c("Q13", "Q14")  ~ "R4"
  )) %>% 
   mutate(Plot = str_c(Belt, Rep)) %>% 
   group_by(Week, Plot, Scientific_name) %>%
   summarise(cover = mean(Rows)) %>% 
   as.data.frame()

view(groundcover)  
  
# to get the species name:

species <-  groundcover %>% 
  filter(Plot=="B1R1", Week=="0") %>% 
  select(Scientific_name) %>% 
  as.data.frame()
species<-species$Scientific_name



# Create a function to make things work faster and avoid too much repetition:

Phylo_sig <- function(Plot, Week) {
  R_cover <- G.cover %>% 
    filter(Plot == !!Plot, Week == !!Week) %>% 
    select(mean_ground) %>% 
    as.data.frame()
  
  R_cover <- as.vector(t(R_cover))
  names(R_cover) <- species 
  
  R_signal <- phylosig(plant_tree$scenario.3, R_cover, 
                       method = "K", nsim = 999)
  
  print(R_signal)}



unique(groundcover$Plot)
# Pre-spray: Week = 0

Phylo_sig(Plot = "B1R1", Week = "0")
Phylo_sig(Plot = "B1R2", Week = "0")
Phylo_sig(Plot = "B1R3", Week = "0")
Phylo_sig(Plot = "B1R4", Week = "0")

Phylo_sig(Plot = "B2R1", Week = "0")
Phylo_sig(Plot = "B2R2", Week = "0")
Phylo_sig(Plot = "B2R3", Week = "0")
Phylo_sig(Plot = "B2R4", Week = "0")

 

