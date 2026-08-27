## to prepare for input data from provided supplementary
## aggregated_RAD51D_E1_10_SNV_9_12_2025.tsv

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
