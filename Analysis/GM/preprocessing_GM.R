## start with the supplementary
## helper function to create the input file for GM analysis

library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)

Supplementary_Data_2_08272026 <- read_excel("Supplementary Data 2_08272026.xlsx", 
                                            skip = 1)
# Data will have the following columns in long format
#      CHROM, from #CHROM
#      POS, # from Annotation_hg38
#      REF, # from REF
#      ALT, # from ALT
#      EventCount, # from Rep1_Day14	Rep1_Day5	Rep1_Lib	Rep2_Day14	Rep2_Day5	Rep2_Lib	Rep3_Day14	Rep3_Day5	Rep3_Lib convert to long format
#      EventType, # from variant class, rename Nonsense to Synonymous and Slient to StopGain
#      SpliceMax, # maximal value of columns 
#      Rep, # 1to3
#      Time, # 0=Lib, 1=Day5, 2=Day14
#      Exon # from SGE
#      use4norm # assgin 1 to all



test <- Supplementary_Data_2_08272026 %>% dplyr::select(`#CHROM`, POS_hg38,REF, ALT,  AA,
                                                        Rep1_Day14,	Rep1_Day5,	Rep1_Lib,	Rep2_Day14,	Rep2_Day5,	Rep2_Lib,	Rep3_Day14,	Rep3_Day5,	Rep3_Lib,
                                                        `variant class`, 
                                                        SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL,
                                                        SGE) %>% 
  rename(Exon = SGE, POS = POS_hg38, CHROM = `#CHROM`) %>% mutate(EventType0 = `variant class`) %>% mutate(AAPOS = as.numeric(str_extract(AA, "\\d+"))) %>%
  mutate(EventType = ifelse(EventType0 == "Nonsense", "Synonymous",
                                           ifelse(EventType0 == "Silent", "StopGain", EventType0))) %>%
  mutate(SpliceMax = pmax( SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL, na.rm= T)) %>%
  mutate(use4norm =1) 

## reformat to long format
test0 <- test %>% 
  pivot_longer(
    cols = starts_with("Rep"),
    names_to = c("Rep", "TimeLabel"),
    names_pattern = "Rep(\\d+)_(.*)",
    values_to = "EventCount"
  ) %>%
  mutate(
    Rep = as.integer(Rep),
    Time = case_when(
      TimeLabel == "Lib" ~ 0,
      TimeLabel == "Day5" ~ 1,
      TimeLabel == "Day14" ~ 2
    )) %>% 
  select(CHROM, POS, REF, ALT, AAPOS, Exon, EventType, SpliceMax, Rep, Time, EventCount,use4norm)

test <- test0
table(test$Time, test$Rep, test$Exon, useNA = 'always')


## remove E10
# Remove variants from PAM site for Exon 3, varying by replicate.
## relabel clivar variants

cv <- read_csv("/fslustre/qhs/ext_na_jie_mayo_edu/MASSAGE_dev/data/RAD51D/internal/clinVarExclusions.csv")

out <- test %>% mutate(uPOS = paste(POS, REF, ALT, sep = "_")) %>% filter(Exon != "E10")

out_filtered_pam <- out %>%
  filter(
    # Keep rows UNLESS they match one of the PAM site conditions in Exon 3
    !(Exon == "E3" & Rep == 1 & AAPOS == 81) &
      !(Exon == "E3" & Rep == 2 & AAPOS == 49) &
      !(Exon == "E3" & Rep == 3 & AAPOS == 65)
  )

# Recreate EventType 
out_final <- out_filtered_pam %>%
  mutate(
    EventType = case_when(
      # Condition 1: cv$variant & EventType=="StopGain" -> "StopGain_cv_exclude"
      uPOS %in% cv$variant & EventType == "StopGain" ~ "StopGain_exclude",
      # Condition 2: cv$variant & EventType=="Synonymous" -> "Synonymous_cv_exclude"
      uPOS %in% cv$variant & EventType == "Synonymous" ~ "Synonymous_exclude",
      # Default: if no conditions match, keep the original EventType
      TRUE ~ EventType
    )
  )

# Display the head of the processed data frame to review changes
head(out_final)
table(out_final$Exon, out_final$Rep, out_final$Time)
table(out_final$Exon, out_final$EventType)

write.table(out_final, "./RAD51D_input_noE10_github_gm.tsv", 
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)