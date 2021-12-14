### Script name: Check_ped_Cal_Amat.R
### Purpose: Manipulate and check pedigree file and calculate the triplet format of pedigree relationship matrix (A-matrix)  
### Author: Sajjad Toghiani
### Released Date: 13 Dec 2020
### Last Modified: 13 Dec 2020

rm(list = ls())

## specify the package names and loading them using packman package 
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(nadiv,dplyr,Matrix,data.table,gdata,optiSel,gtools,ggplot2,readr)

### Call ExamplePed data set built-in optiSel package
### This data set gives a small subset of the pedigree of Hinterwald cattle suitable for demonstration purposes.
### The data frames contains columns for "Indiv" (individual ID),"Sire","Dam","Sex","Breed","Born"(year of birth),"BV"(simulated breeding value)
data(ExamplePed)
head(ExamplePed)

## replace 0 as missing with NA
ExamplePed <- ExamplePed %>% dplyr::na_if(0)

## Count number of NA as missing for each column in pedigree
apply(ExamplePed, 2, function(x) length(which(is.na(x))))

##find duplicate animals in pedigree
ExamplePed %>%  group_by(Indiv) %>% dplyr::filter(n() > 1)

##summarize duplicate animals
ExamplePed %>% group_by(Indiv) %>% dplyr::filter(n()>1) %>% summarize(n=n())

## Remove duplicated rows based on animal ID in pedigree
ped <- ExamplePed %>% dplyr::distinct(Indiv, .keep_all = TRUE)


## Prepares a pedigree by sorting and adding 'founders'
full.ped<-optiSel::prePed(ped)
dim(full.ped)
head(full.ped)

row.names(full.ped)<-NULL   ##reset row index of the pedigree file

## determine the proportion of known ancestors of each specified individual in each ancestral generation
compl_indiv <- optiSel::completeness(full.ped,keep=full.ped$Indiv, by="Indiv")
head(compl_indiv)

##the mean completeness of the pedigrees of specified individuals within sexes
compl_sex <- optiSel::completeness(full.ped,keep=full.ped$Indiv, by="Sex")
head(compl_sex)

## plot the mean completeness of pedigree over generations within sexes
ggplot(compl_sex, aes(x=Generation, y=Completeness, col=Sex)) + geom_line()

########################################################################################
####### The completeness of the pedigree of each individual can be summarized    #######
########################################################################################

## The summarized of complete pedigree contains
## 1) equiGen (Number of equivalent complete generations): 
## the sum over all known ancestors of the terms computed as the sum of (1/2)^n, 
## where n is the number of generations separating the individual to each known ancestor
## 2) fullGen (Number of fully traced generations):
## the farthest generation in which all the ancestors are known
## 3) maxGen (Number of maximum generations traced): 
## the number of generations separating the individual from its farthest ancestor
## 4) PCI (Index of pedigree completeness, which is the harmonic mean of the pedigree completenesses of the parents (MacCluer et al, 1983)
## 5) Inbreeding (Inbreeding coefficient)

compl_pedigree <- summary(full.ped, keep.only=full.ped$Indiv)
head(compl_pedigree)
summary(compl_pedigree)

## Calculating pedigree-based relationship matrix (A matrix) using nadiv package 
Amat<-nadiv::makeA(full.ped[,c(1:3)])

## Total elements of a A-Matrix
total<-dim(Amat)[1]*dim(Amat)[2]

## Number of Non-Zero elements of a Sparse A-Matrix
nonzero<-nnzero(Amat, na.counted = TRUE)

## Number of Zero elements of a Sparse A-Matrix
zero<-total-nonzero

## Build (row, col) pairs from non-zero entries of class "dgCMatrix" (Diagonal+Upper Diagonal)
## https://www.gormanalysis.com/blog/sparse-matrix-construction-and-use-in-r/
df <- as.data.frame(summary(Amat))
df$Row_id <- rownames(Amat)[df$i]
df$Col_id <- colnames(Amat)[df$j]

## manipulate dataframe (df)
df <- df %>% select(Row_id,Col_id,x) %>%
  rename(Amat=x) %>%
  mutate_if(is.numeric, round, digits = 3)

head(df)
tail(df)

## Build triplet dataframe from "dgCMatrix" class (Upper Diagonal matrix)
getUpper <- as.data.frame(subset(summary(Amat), j > i))
getUpper$Row_id <- rownames(Amat)[getUpper$i]
getUpper$Col_id <- colnames(Amat)[getUpper$j]

getUpper <- getUpper %>% select(Row_id,Col_id,x) %>%
  rename(Amat=x) %>%
  mutate_if(is.numeric, round, digits = 3)

head(getUpper)
tail(getUpper)

## Build triplet dataframe from "dgCMatrix" class (Diagonal matrix)
getDiag <- as.data.frame(subset(summary(Amat), j == i))
getDiag$Row_id <- rownames(Amat)[getDiag$i]
getDiag$Col_id <- colnames(Amat)[getDiag$j]

getDiag <- getDiag %>% select(Row_id,Col_id,x) %>%
  rename(Amat=x) %>%
  mutate_if(is.numeric, round, digits = 3)

head(getDiag)
tail(getDiag)

##create triplet dataframe for Lower Diagonal matrix from reverse the i and j column of Upper Diagonal
getLower <- getUpper %>% select(Col_id,Row_id,Amat) %>%
  rename(Row_id=Col_id,Col_id=Row_id)

head(getLower)
tail(getLower)

## rbind all three set of triplet dataframes from the A-matrix
df2<-data.table::rbindlist(list(getDiag,getUpper,getLower))  %>% dplyr::arrange(Row_id,Col_id)


## maximum length of string/value in all columns
max_len_df<- apply(df2, 2, function(x) max(nchar(x)))

# Using the write.fwf from the gdata package if you want to create a Fixed Width File.
## The col names don't seem to be saved to the correct position,
## so to work around this you can format them to the width you want the columns to be using formatC
colnames(df2) <- formatC(colnames(df2), format = "d",width = 16, flag = "0")
gdata::write.fwf(df2,"Amat_long_format", width = rep(16,3),rownames = F,quote = F,justify="right")


