#READING IN THE DATA---------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/TIMAVID/The-utility-of-evolutionary-allometries-in-estimating-body-sizes-of-phrynosomatid-lizards/master/Phrynosomatid_Measurements.csv")
Lizard_Measurements <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of measured specimens 
head(Lizard_Measurements)

f3 <- curl("https://raw.githubusercontent.com/TIMAVID/Allometric-patterns-in-phrynosomatid-lizards-and-the-implications-for-reconstruc-body-size-for-fossi/master/Lizard%20Repeated%20Measurements.csv")
Lizard_Repeated_Measurements <- read.csv(f3, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of repeatedly measured specimens  
head(Lizard_Repeated_Measurements)

#FILTERING THE DATA TO CREATE DATASETS ------------
library(dplyr)
library(forcats)
Phrynosomatids <- filter(Lizard_Measurements,grepl('Sceloporus|Holbrookia|Cophosaurus|Urosaurus|Uta|Uma|Callisaurus|Phrynosoma|Petrosaurus',Specimen))
Phrynosomatids <- filter(Phrynosomatids,!grepl('Fossil',genus))
set.seed(123)
split<- sample(c(rep(0, 0.8 * nrow(Phrynosomatids)), rep(1, 0.2 * nrow(Phrynosomatids))))
table(split)
Phrynosomatid_train <- Phrynosomatids[split == 0, ] 
Phrynosomatid_test <- Phrynosomatids[split== 1, ]    

Fossils <- filter(Lizard_Measurements,grepl('Fossil',genus))
non_S.occident_phrynos <- filter(Phrynosomatids,!grepl('Sceloporus occidentalis',Specimen))
set.seed(123)
split2<- sample(c(rep(0, 0.6 * nrow(non_S.occident_phrynos)), rep(1, 0.4 * nrow(non_S.occident_phrynos))))
table(split2)
non_S.occident_phrynos_train <- non_S.occident_phrynos[split2== 1, ]    

Sceloporus_occidental <- filter(Phrynosomatids, grepl('Sceloporus occidentalis',Specimen))
non_S.occidentalis.Scelop <- filter(non_S.occident_phrynos, grepl('Sceloporus',Specimen))
non_Scelop_phrynos <- filter(Phrynosomatids,!grepl('Sceloporus',Specimen))


Sceloporus <- filter(Phrynosomatids,grepl('Sceloporus',Specimen))

Sceloporine <- filter(Phrynosomatids,grepl('Sceloporus|Uta|Urosaurus|Petrosaurus',Specimen))
set.seed(123)
split1<- sample(c(rep(0, 0.8 * nrow(Sceloporine)), rep(1, 0.2 * nrow(Sceloporine))))
table(split1)
Sceloporine_train <- Sceloporine[split1 == 0, ] 
Sceloporine_test <- Sceloporine[split1== 1, ] 

Phrynosomatine <-filter(Phrynosomatids,grepl('Holbrookia|Cophosaurus|Uma|Callisaurus|Phrynosoma',Specimen))
set.seed(123)
split2<- sample(c(rep(0, 0.8 * nrow(Phrynosomatine)), rep(1, 0.2 * nrow(Phrynosomatine))))
table(split2)
Phrynosomatine_train <- Phrynosomatine[split2 == 0, ] 
Phrynosomatine_test <- Phrynosomatine[split2== 1, ]    

Cophosaurus_texanus <- filter(Phrynosomatids, grepl('Cophosaurus texanus',Specimen))

AvgSceloporus <-droplevels(Sceloporus)
AvgSceloporus$species<-AvgSceloporus$species %>% fct_collapse('Sceloporus occidentalis' = c("Sceloporus occidentalis","Sceloporus occidentalis biseriatus"))
AvgSceloporus$species<-AvgSceloporus$species %>% fct_collapse('Sceloporus graciosus' = c("Sceloporus graciosus","Sceloporus graciosus vandenburgianus"))
levels(AvgSceloporus$species)
AvgSceloporus <- AvgSceloporus %>% 
  group_by(species) %>%
  summarise(across(Maxilla_LDR:SVL, mean, na.rm = TRUE))

Sceloporus_occidental$species<-Sceloporus_occidental$species %>% fct_collapse('Sceloporus occidentalis' = c("Sceloporus occidentalis","Sceloporus occidentalis biseriatus"))
Sceloporus$species<-Sceloporus$species %>% fct_collapse('Sceloporus occidentalis' = c("Sceloporus occidentalis","Sceloporus occidentalis biseriatus"))
Sceloporus$species<-Sceloporus$species %>% fct_collapse('other Sceloporus' = c("Sceloporus graciosus","Sceloporus graciosus vandenburgianus","Sceloporus clarkii","Sceloporus graciosus","Sceloporus grammicus",   
                                                                               "Sceloporus jarrovii"  ,   "Sceloporus licki"    ,    "Sceloporus magister" ,   
                                                                               "Sceloporus olivaceous"  , "Sceloporus orcutti" ,    
                                                                                  "Sceloporus poinsettii"  , "Sceloporus undulatus"  , 
                                                                               "Sceloporus undulatus tristichus" ,  "Sceloporus virgatus" ))
non_S.occidentalis.Scelop$species<-non_S.occidentalis.Scelop$species %>% fct_collapse('other Sceloporus' = c("Sceloporus graciosus","Sceloporus graciosus vandenburgianus","Sceloporus clarkii","Sceloporus graciosus","Sceloporus grammicus",   
                                                                               "Sceloporus jarrovii"  ,   "Sceloporus licki"    ,    "Sceloporus magister" ,   
                                                                               "Sceloporus olivaceous"  , "Sceloporus orcutti" ,    
                                                                                  "Sceloporus poinsettii"  , "Sceloporus undulatus"  , 
                                                                               "Sceloporus undulatus tristichus" ,  "Sceloporus virgatus" ))


# NEW Bones function TO CREATE LINEAR MODELS FOR PREDICTING SVL------------
varlist <- names(Lizard_Measurements)[2:16] #all the different measurements

bones.lm<-function(variables,data){
  require(ggplot2)
  require(gridExtra)
  figs<-lapply(variables, function(x) {
    ggplot(data = data, 
           aes(log(get(x)),log(SVL),)) + geom_point()+ ggtitle(x)+ theme_classic()+ ylab("log(Measurement)")})
  do.call(grid.arrange, c(figs, ncol=3, top = deparse(substitute(data))))
  
  models <- lapply(variables, function(x) { #function to perform linear regression on all measurements for each dataset
    lm(substitute(log(SVL) ~ log(i), list(i = as.name(x))), data = data)
  })
  names(models) <- variables
  (sum <- (lapply(models, summary)))
  b<-(lapply(sum, function (x)x$coefficients[1]))
  m<-(lapply(sum, function (x)x$coefficients[2]))
  R<-(lapply(sum, function (x)x$adj.r.squared))
  P<-(lapply(sum, function (x)x$coefficients[2,4])) #may need to p.adjust
  MSE <- (lapply(sum, function (x)(mean(x$residuals^2))))
  Con <- (lapply(models, confint))
  C<-(lapply(Con, function (x)x[2,]))
  out<-list(m,b,R,P,MSE,models,C)
  names(out)<-c("slope","Y-intercept","adj.R-squared","P-value","MSE","models","ConfidenceInt")
  out <- do.call(rbind, out)
  out <- t(out)
  out <- as.data.frame(out)
  return(out)
}


# NEW estimating SVL function-------------
estimate_SVL <- function(lm, data #need to check specific columns [,2:16] are ones to estimate from
) #function to estimate SVL and calculate difference between actual and estimated SVL
{
  t<- t(lm)
  t <- as.data.frame(t)
  slopes <-as.vector(t[1,])
  intercepts <- as.vector(t[2,])
  estimates <- mapply(function(m,b,data){{(exp(as.numeric(m)*
                                                 log(data) + as.numeric(b)))}}, 
                      data= data[,2:16], m=slopes, b=intercepts)
  diff <- estimates - data$SVL
  percentdiff <- (abs(estimates - data$SVL)/data$SVL)*100
  data<-list(estimates,data$SVL,diff,percentdiff, data$Specimen, data$species, data$genus)
  names(data)<-c("Estimated_SVL","SVL","Difference","Percentage difference", "Specimen", "Species", "Genus")
  return((data))
}

# FUNCTIONS TO CREATE LINEAR MODELS FOR EXAMINING ALLOMETRY AT DIFFERENT TAXONOMIC LEVELS------------
plot.alom.species<-function(variables,data){
  require(ggplot2)
  require(gridExtra)
  figs<-lapply(variables, function(x) {
    ggplot(data = data,
           aes(log(SVL), log(get(x)))) + geom_point(aes(color=species))+
      scale_colour_manual(values=cbPalette,breaks=c("Sceloporus occidentalis", "other Sceloporus"))+ggtitle(x)+
      theme_classic(base_size = 8)+ ylab("log(Measurement)")+
      geom_smooth(method='lm')+theme(legend.position="none")+ coord_fixed(ratio = 1)
    # +geom_abline(intercept=seq(-100, 100, 1),
    #               slope=1,
    #               colour="red",linetype='dashed')
  })
  do.call(grid.arrange, c(figs, ncol=5, top = deparse(substitute(data))))
  
  models <- lapply(variables, function(x) { #function to perform linear regression on all measurements for each dataset
    lm(substitute(log(i) ~ log(SVL), list(i = as.name(x))), data = data)
  })
  names(models) <- variables
  (sum <- (lapply(models, summary)))
  b<-(lapply(sum, function (x)x$coefficients[1]))
  m<-(lapply(sum, function (x)x$coefficients[2]))
  R<-(lapply(sum, function (x)x$adj.r.squared))
  P<-(lapply(sum, function (x)x$coefficients[2,4])) #may need to p.adjust
  Con <- (lapply(models, confint))
  C<-(lapply(Con, function (x)x[2,]))
  out<-list(m,b,R,P,models,C)
  names(out)<-c("slope","Y-intercept","adj.R-squared","P-value","models","ConfidenceInt")
  out <- do.call(rbind, out)
  out <- t(out)
  out <- as.data.frame(out)
  return(out)
}
cbPalette <- c( "#6D0C96","#EFA134")

plot.alom.genus<-function(variables,data){
  require(ggplot2)
  require(gridExtra)
  figs<-lapply(variables, function(x) {
    ggplot(data = data, 
           aes(log(SVL), log(get(x)))) + geom_point(aes(color=genus))+
      scale_colour_manual(values=cbPalette2, breaks=c("Callisaurus", "Cophosaurus", "Petrosaurus", "Phrynosoma","Sceloporus","Uma","Urosaurus","Uta", "Holbrookia"))+
      ggtitle(x)+
      theme_classic(base_size = 8)+ ylab("log(Measurement)")+
      geom_smooth(method='lm')+theme(legend.position="none")+ coord_fixed(ratio = 1)
    # +geom_abline(intercept=seq(-100, 100, 1),
    #               slope=1,
    #               colour="red",linetype='dashed')
    
  })
  do.call(grid.arrange, c(figs, ncol=5, top = deparse(substitute(data))))
  
  models <- lapply(variables, function(x) { #function to perform linear regression on all measurements for each dataset
    lm(substitute(log(i) ~ log(SVL), list(i = as.name(x))), data = data)
  })
  names(models) <- variables
  (sum <- (lapply(models, summary)))
  b<-(lapply(sum, function (x)x$coefficients[1]))
  m<-(lapply(sum, function (x)x$coefficients[2]))
  R<-(lapply(sum, function (x)x$adj.r.squared))
  P<-(lapply(sum, function (x)x$coefficients[2,4])) #may need to p.adjust
  Con <- (lapply(models, confint))
  C<-(lapply(Con, function (x)x[2,]))
  out<-list(m,b,R,P,models,C)
  names(out)<-c("slope","Y-intercept","adj.R-squared","P-value","models","ConfidenceInt")
  out <- do.call(rbind, out)
  out <- t(out)
  out <- as.data.frame(out)
  return(out)
}
cbPalette2 <- c("#D4185D", "#1174CA", "#FFF500", "#5A3E06","#EFA134", "#43B546", "#EF44D5", "#3CE4C8", "#ADADAD")


## ALL LINEAR REGRESSION MODELS USING S. OCCIDENTALIS DATASET--------------------
Sceloporus_occidental_lm <-bones.lm(varlist,Sceloporus_occidental) #liner models of all measurements

Sceloporus_occidental_alom <-plot.alom.species(varlist,Sceloporus_occidental) #liner models of all measurements

#lm(Sceloporus occidentalis ) -predict> non-S. occidentalis Sceloporus

non.occi.sceloporus_estimates <- estimate_SVL(Sceloporus_occidental_lm, non_S.occidentalis.Scelop)

#lm(Sceloporus occidentalis) -predict> non-S. occidentalis phrynosomatids

non_S.occident_phrynos_estimates <- estimate_SVL(Sceloporus_occidental_lm, non_S.occident_phrynos)


## ALL LINEAR REGRESSION MODELS USING SCELOPORUS DATASET--------------------
Sceloporus_lm <-bones.lm(varlist,Sceloporus) #liner models of all measurements

Sceloporus_alom <-plot.alom.species(varlist,Sceloporus) #liner models of all measurements

#lm(Sceloporus) -predict> non-Sceloporus phrynosomatids

non_Scelop_phrynos_estimates <- estimate_SVL(Sceloporus_lm, non_Scelop_phrynos)


## ALL LINEAR REGRESSION MODELS USING non-S.OCCIDENTALIS Sceloporus DATASET--------------------
non_S.occidentalis.Scelop_lm <-bones.lm(varlist,non_S.occidentalis.Scelop) #liner models of all measurements

non_S.occidentalis.Scelop_alom <-plot.alom.species(varlist,non_S.occidentalis.Scelop) #liner models of all measurements

#lm(non-S. occidentalis Sceloporus) -predict> Sceloporus occidentalis

Sceloporus_occidental_estimates <- estimate_SVL(non_S.occidentalis.Scelop_lm, Sceloporus_occidental)

## ALL LINEAR REGRESSION MODELS USING Sceloporine DATASET--------------------
Sceloporine_lm <-bones.lm(varlist,Sceloporine_train) #liner models of all measurements

Sceloporine_alom <-plot.alom.genus(varlist,Sceloporine) #liner models of all measurements

#lm(Sceloporines) -predict> subset Sceloporines

Sceloporine_estimates <- estimate_SVL(Sceloporine_lm, Sceloporine_test)

#lm(Sceloporines) -predict> subset Phrynosomatids
Phrynosomatid_test2 <- rbind(Sceloporine_test, Phrynosomatine_test)
Phrynosomatid_estimates <- estimate_SVL(Sceloporine_lm, Phrynosomatid_test2)


## ALL LINEAR REGRESSION MODELS USING Phrynosomatinae DATASET--------------------
Phrynosomatine_lm <-bones.lm(varlist,Phrynosomatine_train) #liner models of all measurements

Phrynosomatine_alom <-plot.alom.genus(varlist,Phrynosomatine) #liner models of all measurements

#lm(Phrynosomatines) -predict> subset Phrynosomatines

Phrynosomatine_estimates <- estimate_SVL(Phrynosomatine_lm, Phrynosomatine_test)

#lm(Phrynosomatines) -predict> subset Phrynosomatids

Phrynosomatid_test2 <- rbind(Sceloporine_test, Phrynosomatine_test)
Phrynosomatid_estimates2 <- estimate_SVL(Phrynosomatine_lm, Phrynosomatid_test2)


## ALL LINEAR REGRESSION MODELS USING Cophosaurus_texanus DATASET--------------------
Cophosaurus_texanus_alom <-plot.alom.species(varlist,Cophosaurus_texanus) #liner models of all measurements

## ALL LINEAR REGRESSION MODELS USING non-Sceloporus phrynosomatids DATASET--------------------
non_Scelop_phrynos_lm <-bones.lm(varlist,non_Scelop_phrynos) #liner models of all measurements

non_Scelop_phrynos_alom <-plot.alom.genus(varlist,non_Scelop_phrynos) #liner models of all measurements

#lm(non-Sceloporus phrynosomatids) -predict> Sceloporus

Sceloporus_estimates <- estimate_SVL(non_Scelop_phrynos_lm, Sceloporus)


## ALL LINEAR REGRESSION MODELS USING non-S.occidentalis phrynosomatids DATASET--------------------
non_S.occident_phrynos_lm <-bones.lm(varlist,non_S.occident_phrynos) #liner models of all measurements

non_S.occident_phrynos_lm_sub <-bones.lm(varlist,non_S.occident_phrynos_train) #liner models of all measurements

non_S.occident_phrynos_alom <-plot.alom.genus(varlist,non_S.occident_phrynos) #liner models of all measurements

#lm(non-S. occidentalis phrynosomatids) -predict> Sceloporus occidentalis

Sceloporus_occidental_estimates2 <- estimate_SVL(non_S.occident_phrynos_lm, Sceloporus_occidental)

Sceloporus_occidental_estimates_sub <- estimate_SVL(non_S.occident_phrynos_lm_sub, Sceloporus_occidental)


## ALL LINEAR REGRESSION MODELS USING AVERAGE SCELOPORUS DATASET--------------------
AvgSceloporus_lm <-bones.lm(varlist,AvgSceloporus) #liner models of all measurements


## ALL LINEAR REGRESSION MODELS USING ALL PHRYNOSOMATID DATASET--------------------
Phrynosomatids_lm <-bones.lm(varlist,Phrynosomatids) #liner models of all measurements

Phrynosomatids_alom <-plot.alom.genus(varlist,Phrynosomatids) #liner models of all measurements

#lm(Phrynosomatids) -predict> subset Phrynosomatids

Phrynosomatid_estimates3 <- estimate_SVL(Phrynosomatids_lm, Phrynosomatid_test)

#lm(Phrynosomatids) -predict> fossils------------

Fossil_estimates <- estimate_SVL(Phrynosomatids_lm, Fossils)
library(Hmisc)
Fossil_estimates<- llist(Fossil_estimates[c(1,5)])

Fossil_estimates <- unlist(Fossil_estimates,recursive=FALSE)
Fossil_estimates <- lapply(Fossil_estimates, data.frame, stringsAsFactors = FALSE)
Fossil_estimates<-bind_cols(Fossil_estimates, .id = "column_label")

library(dplyr)
library(janitor)

Fossil_estimates <- Fossil_estimates %>%
  mutate_all(funs(na_if(., ""))) %>%
  remove_empty("cols")

# write.table(Fossil_estimates, file = "Fossil_estimates", sep = ",", quote = FALSE, row.names = T)


#All the percent differences from estimates-----------------------------
All_estimates<- llist(non.occi.sceloporus_estimates[c(4)],non_S.occident_phrynos_estimates[c(4)], non_Scelop_phrynos_estimates[c(4)], 
                      Sceloporus_occidental_estimates[c(4)], Sceloporus_occidental_estimates2[c(4)],Sceloporus_occidental_estimates_sub[c(4)], Sceloporine_estimates[c(4)], Phrynosomatine_estimates[c(4)],
                      Sceloporus_estimates[c(4)], Phrynosomatid_estimates[c(4)] ,Phrynosomatid_estimates2[c(4)], Phrynosomatid_estimates3[c(4)])

#FUNCTION: to plot all the percent differences-----------------------------
plot.percent.diff<-function(data){ 
  yo<-names(data)
  
  new.pp <- unlist(data,recursive=FALSE)
  dfs <- lapply(new.pp, data.frame, stringsAsFactors = FALSE)
  dfs <- lapply(dfs, function(x) {tibble::rownames_to_column((data.frame(t(x))),"Measurement")})
  dfs<- lapply(dfs, function(x) {tidyr::gather(x,obs, value, 2:length(x))})
  
  ggBox <- function(x, name) {
    ggplot(data = x, aes(x=Measurement, y=value)) +
      geom_boxplot() + ggtitle(name)+
      theme_classic()+ ylab("Percent diff") +
      scale_x_discrete(guide = guide_axis(n.dodge=3))+
      geom_hline(yintercept=10)+
      geom_hline(yintercept=15)
    #+ylim(0, 50)
  }
  figs<-Map(f = ggBox, dfs, yo)
  do.call(grid.arrange, c(figs, ncol=3, top = deparse(substitute(data))))
}


plot.percent.diff(All_estimates)


#Change All_estimates into a list of dataframes 
diff.dfs <- unlist(All_estimates,recursive=FALSE)
diff.dfs <- lapply(diff.dfs, data.frame, stringsAsFactors = FALSE)
diff.dfs <- lapply(diff.dfs, function(x) {tibble::rownames_to_column((data.frame(t(x))),"Measurement")})
diff.dfs<- lapply(diff.dfs, function(x) {tidyr::gather(x,obs, value, 2:length(x))})

ChangeType <- function(DF){ #change Measurements to a factor
  DF[,1] <- as.factor(DF[,1])
  DF #return the data.frame 
}
diff.dfs <- lapply(diff.dfs, ChangeType) # store the returned value to df.list, thus updating your existing data.frame

#FUNCTION: Relative standard deviation-------------------------------
RSD <- function(x, na.rm = TRUE){ 
  100*sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
}

#FUNCTION: summaries the percent differences and return mean, median and sd-----------------------
Diff.summ <- function(diff.data){ 
  require(dplyr)
  diff.data %>% group_by(Measurement) %>% summarise_at(c("value"),funs(median, mean, sd), na.rm = TRUE)
}

diff.summary <- lapply(diff.dfs, Diff.summ)

diff.summary <- do.call(rbind, diff.summary)

# write.table(diff.summary, file = "diff summary", sep = ",", quote = FALSE, row.names = T)


#FUNCTION: to plot all the percent difference means-------------------------------
# plot.percent.diff.mean<-function(data){ 
#   yo<-names(data)
#   
#   ggBox <- function(x, name) {
#     ggplot(data = x, aes(x=Measurement, y=mean)) +
#       geom_point() + ggtitle(name)+
#       theme_classic()+ ylab("Mean % diff") +
#       scale_x_discrete(guide = guide_axis(n.dodge=3))+
#       geom_hline(yintercept=10)+
#       geom_hline(yintercept=15) +
#       geom_errorbar(width=.1, aes(ymin=mean, ymax=mean+sd)) 
#     #+ylim(0, 50)
#   }
#   figs<-Map(f = ggBox, data, yo)
#   do.call(grid.arrange, c(figs, ncol=2, top = deparse(substitute(data))))
# }
# 
# plot.percent.diff.mean(diff.summary)

#FUNCTION: to plot all the differences with color-----------------------------
plot.diff.color<-function(data){ 
  yo<-names(data)
  require(purrr)
  require(ggplot2)
  require(data.table)
  require(gridExtra)
  dfs <- map(data, as.data.table)
  dfs<- lapply(dfs, function(x) 
  {tidyr::gather(x,Measurement, value, 
                 Difference.Maxilla_LDR:Difference.Ilium_crest_GL)})
  
  ggBox <- function(x, name) {
    ggplot(data = x, aes(x=Measurement, y=value)) +
      geom_point(aes(color=Genus),position=position_jitter(width=.1,height=0)) + 
      scale_colour_manual(values=cbPalette2, breaks=c("Callisaurus", "Cophosaurus", "Petrosaurus", "Phrynosoma","Sceloporus","Uma","Urosaurus","Uta", "Holbrookia"))+ ggtitle(name)+
      ylab("Diff(Est.-Act.)") +
      scale_x_discrete(guide = guide_axis(angle = 90)) +
      theme_classic()
  }
  figs<-Map(f = ggBox, dfs, yo)
  do.call(grid.arrange, c(figs, ncol=2, top = deparse(substitute(data))))
}

Some.diff <- llist(non_Scelop_phrynos_estimates[c(3,6,7)], non_S.occident_phrynos_estimates[c(3,6,7)])

plot.diff.color(Some.diff)

#FUNCTION:summaries the percent differences and return mean,sd, RSD------------------
Diff.all.summ <- function(diff.data){ 
  require(dplyr)
  diff.data %>% summarise_at(c("value"),funs(mean, sd, RSD), na.rm = TRUE)
}

diff.all.summary <- lapply(diff.dfs, Diff.all.summ)

#out <- do.call(rbind, diff.all.summary)
#write.table(out, file = "%diff summary", sep = ",", quote = FALSE, row.names = T)





# Adjusting lm p values using fdr---------------------------------
lm.results <- llist(Sceloporus_occidental_lm[c(1:4)],non_S.occidentalis.Scelop_lm[c(1:4)],
                    AvgSceloporus_lm[c(1:4)], Sceloporus_lm[c(1:4)],
                    Sceloporine_lm[c(1:4)], Phrynosomatine_lm[c(1:4)], non_S.occident_phrynos_lm[c(1:4)], 
                    non_Scelop_phrynos_lm[c(1:4)], Phrynosomatids_lm[c(1:4)])

q.values <-lapply(lm.results, function(x){
  q <- p.adjust(unlist(x[4]),method = "fdr")
  return(q)
})

q.values <- as.data.frame(q.values)

lm.results <- unlist(lm.results,recursive=FALSE)
lm.results <- lapply(lm.results, data.frame, stringsAsFactors = FALSE)
lm.results <- lapply(lm.results, function(x) {tibble::rownames_to_column((data.frame(t(x))),"Measurement")})
lm.results <- as.data.frame(lm.results)

# write.table(lm.results, file = "lm.results", sep = ",", quote = FALSE, row.names = F)
# write.table(q.values, file = "q.values", sep = ",", quote = FALSE, row.names = T)


# Adjusting alom p values using fdr---------------------------------
alom.results <- llist(Sceloporus_occidental_alom[c(1:4,6)],non_S.occidentalis.Scelop_alom[c(1:4,6)], Cophosaurus_texanus_alom[c(1:4,6)],
                      #Sceloporus_alom[c(1:4,6)],
                      Sceloporine_alom[c(1:4,6)], Phrynosomatine_alom[c(1:4,6)], 
                      non_S.occident_phrynos_alom[c(1:4,6)]
                      #,non_Scelop_phrynos_alom[c(1:4,6)], Phrynosomatids_alom[c(1:4,6)]
                      )

alom_q.values <-lapply(alom.results, function(x){
  q <- p.adjust(unlist(x[4]),method = "fdr")
  return(q)
})

alom_q.values <- as.data.frame(alom_q.values)

alom.results <- unlist(alom.results,recursive=FALSE)
alom.results <- lapply(alom.results, data.frame, stringsAsFactors = FALSE)
alom.results <- lapply(alom.results, function(x) {tibble::rownames_to_column((data.frame(t(x))),"Measurement")})
alom.results <- as.data.frame(alom.results)

# write.table(alom.results, file = "alom.results", sep = ",", quote = FALSE, row.names = F)
# write.table(alom_q.values, file = "alom_q.values", sep = ",", quote = FALSE, row.names = T)

# Predicting SVL of NEW Sceloporus specimens------------------------------------------
f4 <- curl("https://raw.githubusercontent.com/TIMAVID/Allometric-patterns-in-phrynosomatid-lizards-and-the-implications-for-reconstruc-body-size-for-fossi/master/New%20Sceloporus%20to%20test%20-%20Sceloporus%20to%20test.csv")
New_Sceloporus_to_test_Sceloporus_to_test <- read.csv(f4, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of repeatedly measured specimens  
head(New_Sceloporus_to_test_Sceloporus_to_test)


New_Sceloporus_to_test_Sceloporus_to_test <- New_Sceloporus_to_test_Sceloporus_to_test %>%
  dplyr::rename(Specimen = Measurement,
                'Occipital_complex_WC' = 'Occipital.complex_WC',
                'Ilium_crest_GL' = 'Ilium.crest_GL')
New_Sceloporus_to_test_Sceloporus_to_test$species <-gsub("_.*","", New_Sceloporus_to_test_Sceloporus_to_test$Specimen) #makes species column
New_Sceloporus_to_test_Sceloporus_to_test$genus <-gsub(" .*","", New_Sceloporus_to_test_Sceloporus_to_test$Specimen) #makes species column

NEWscelop.estimates <- estimate_SVL(Sceloporus_lm, New_Sceloporus_to_test_Sceloporus_to_test)

NEWscelopALL.estimates <- estimate_SVL(Phrynosomatids_lm, New_Sceloporus_to_test_Sceloporus_to_test)

NEWscelopAvg.estimates <- estimate_SVL(AvgSceloporus_lm, New_Sceloporus_to_test_Sceloporus_to_test)


NEW.Sceloporus_occidental <- filter(New_Sceloporus_to_test_Sceloporus_to_test, grepl('Sceloporus occidentalis',Specimen))
NEWscelop.occidental.estimates <- estimate_SVL(Sceloporus_occidental_lm, NEW.Sceloporus_occidental)


New.estimates <- llist(NEWscelop.occidental.estimates[c(4)], NEWscelop.estimates[c(4)], NEWscelopAvg.estimates[c(4)], NEWscelopALL.estimates[c(4)])

#Change All NEW Sceloporus specimens estimates into a list of dataframes 
New.estimates <- unlist(New.estimates,recursive=FALSE)
New.estimates <- lapply(New.estimates, data.frame, stringsAsFactors = FALSE)
New.estimates <- lapply(New.estimates, function(x) {tibble::rownames_to_column((data.frame(t(x))),"Measurement")})
New.estimates<- lapply(New.estimates, function(x) {tidyr::gather(x,obs, value, 2:length(x))})


New.estimates <- lapply(New.estimates, ChangeType) # store the returned value to df.list, thus updating your existing data.frame


NEWdiff.all.summary <- lapply(New.estimates, Diff.summ)

out <- do.call(rbind, NEWdiff.all.summary)
# write.table(out, file = "%NEWdiff summary", sep = ",", quote = FALSE, row.names = T)


# Measurement error--------------------------
library(readr)
library(dplyr)
library(broom)
Lizard_Repeated_Measurements <- read_csv("Lizard Repeated Measurements.csv")

Lizard_Repeated_Measurements$Specimen <- as.factor(Lizard_Repeated_Measurements$Specimen)


# dt_result = Lizard_Repeated_Measurements %>% group_by(Specimen) %>% do(tidy(t.test(Parietal_GW~User, data=.)))
# dt_result



RSD <- function(x){ #RELATIVE STANDARD DEVIATION
  100*sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
}

Measure_error_USER <- function(data){ #summaries the percent differences and return mean and sd
  require(dplyr)
  dat<- data %>% group_by(Specimen,User) %>%
    dplyr::summarise(across(Maxilla_LDR:SVL, list(mean = mean, sd = sd, RSD = RSD)))
}

Measure_error_TOTAL <- function(data){ #summaries the percent differences and return mean and sd
  require(dplyr)
  dat<- data %>% group_by(Specimen) %>%
    dplyr::summarise(across(Maxilla_LDR:SVL, list(mean = mean, sd = sd, RSD = RSD)))
}

Repeat.meaures_USER<- Measure_error_USER(Lizard_Repeated_Measurements)
Repeat.meaures_USER <-data.frame((Repeat.meaures_USER))

# write.table(Repeat.meaures_USER, file = "Repeat.meaures_USER", sep = ",", quote = FALSE, row.names = T)

library(tidyr)
library(tidyverse)
library(broom)
library(fs)
library(lubridate)

Repeat.meaures_TOTAL<- Measure_error_TOTAL(Lizard_Repeated_Measurements)
data_long <- gather(Repeat.meaures_USER, measurement, value, Maxilla_LDR_mean:Ilium_crest_GL_RSD, factor_key=TRUE)
data_long <- select(data_long, -SVL_sd, -SVL_RSD)
data_long <- filter(data_long,grepl('RSD',measurement)) %>% filter(!is.na(value))
data_long$measurement <- droplevels(data_long$measurement)

data_nest = data_long %>% group_by(measurement,User) %>% nest()

cor_fun <- function(df) cor.test(df$value, df$SVL_mean, method = "pearson") %>% tidy()

data_nest <- mutate(data_nest, model = map(data, cor_fun))
head(data_nest)
corr_pr <- select(data_nest, -data) %>% unnest()

corr_pr <- mutate(corr_pr, sig = ifelse(p.value <0.05, "Sig.", "Non Sig."))



library(psych)
gp.cor <- data_long %>%
  split(.$measurement) %>%  
  map(~corr.test(x = .x %>% select(SVL_mean, value),
                 use = "everything",
                 method = "pearson",
                 adjust = "none",
                 alpha = 0.05,
                 ci = TRUE, minlength = 5)
  )
map(gp.cor, ~.x$r)


# Repeat.meaures_TOTAL <-data.frame(t(Repeat.meaures_TOTAL))

# write.table(Repeat.meaures_TOTAL, file = "Repeat.meaures_TOTAL", sep = ",", quote = FALSE, row.names = T)



