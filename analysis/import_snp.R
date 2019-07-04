# fagus -------------------------------------------------------------------
import.snp.fagus <- function(threshold_loci, threshold_ind, maf){

library(tidyverse)
library(poppr)
library(pegas)
  
snp <- read.csv("../data/Genotyping-2137.011-02 Grid_reformated.csv", 
                header = T, 
                na.strings = c("?", "Uncallable", "Bad", "Missing")
                # ,stringsAsFactors = T
                , check.names = F # default is TRUE, and changes the dashes to dots - which created problems
)
# snp



 ## Transform sample names 
   #### This will allow hierarchical analysis when applicable 
   #### Country / Species / Pop / Plot 

snp <- snp %>% 
  mutate(Genotype = str_replace_all(Genotype, "NR", "NR_")) %>% 
  mutate(Genotype = str_replace_all(Genotype, "^AFS", "GR_FSY_A_")) %>% #GR_Adult
  mutate(Genotype = str_replace_all(Genotype, "^RFS", "GR_FSY_NR_")) %>%  #GR_Regen
  mutate(Genotype = str_replace_all(Genotype, "_A_", "_A_1_")) %>% 
  
  mutate(Genotype = str_replace_all(Genotype, "GR_FSY_NR", "GR_FSY_NR_1"))



 
 ### Create a df following the guidelines of the loci format 

snp_loci_format <- as.data.frame(sapply(snp, gsub, pattern=':', replacement='/') )
# snp_loci_format <- as.data.frame(sapply(snp_loci_format, gsub, pattern=".", replacement="_"))
snp_loci_format <- snp_loci_format[,-1]
rownames(snp_loci_format) <- snp[,1]

# snp_loci_format



## Add a second column for population

pop1 <- replicate(147, "DE_FSY_A")
pop2 <- replicate(72, "DE_FSY_NR")
pop3 <- replicate(147, "GR_FSY_A")
pop4 <- replicate(72, "GR_FSY_NR")
pop5 <- replicate(147, "SI_FSY_A")
pop6 <- replicate(72, "SI_FSY_NR")


pop <- c(pop1, pop2, pop3, pop4, pop5, pop6)

# pop

snp_loci_format <- add_column(snp_loci_format, pop, .before = "4_272")

# snp_loci_format


 ## Create genind object 

library(pegas)

data <- as.loci(snp_loci_format, 
                col.pop = 1
                ,allele.sep = "/")
# data

obj_origin <- loci2genind(data,)
# obj_origin


 ### stratify data set 

strata_df <- as.data.frame(snp_loci_format$pop[-1]) # [-1] because 1 ind was removed from DE_A
colnames(strata_df) <- "strata"
strata_df <- separate(strata_df, col = strata, sep="_", 
                      into = c("Country", 
                               "Species", 
                               "Pop"))
strata(obj_origin) <- strata_df

setPop(obj_origin) <- ~Country/Pop
obj_origin


 ### Remove uninformative (monomorphic) loci 

maf <- 0 

obj <- informloci(obj_origin, MAF = maf)


 ### Filter out missing data 

threshold_loci <- 0.1
threshold_ind <- 0.1

obj <- missingno(obj, type = "loci", cutoff = threshold_loci)

obj <- missingno(obj, type = "genotypes", cutoff = threshold_ind)

test <- seppop(obj) %>% 
  lapply(propTyped, by="loc") %>% 
  lapply(sort)

# capture.output(
# info_table(obj, type = "missing", plot = TRUE, plotlab = F),
# file='NUL')




maf <- 5 / nInd(obj)

obj <- informloci(obj, MAF = maf)


return(obj)
}
