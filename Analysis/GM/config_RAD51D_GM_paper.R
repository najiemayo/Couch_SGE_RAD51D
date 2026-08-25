## this contains the options that users need to provide

## data input file in tsv file format that contains the following columns
## CHROM      POS REF               ALT   EventType Exon Rep Time EventCount use4norm
DATAFILE = "./RAD51D_input_noE10.tsv"
OUTDIR <-  "./result/"

RATIO_COMP = "T2_T1" ##"T2_T1" ## T2_T0
ANALYSIS_MODEL = "MM" ## MM varcall...

## MM OPTIONS
LOWESS_ADJ = FALSE ## TRUE #TRUE FALSE

LOWESS_SUBSET = "all" ## "all" "ModFindlay2"
NORM_BASE = 'all' ## 'all' or "SNV"
ModFindlay_CF =0.8 ## ModFindlay2 cutoff ie
LOWESS_SPAN = 0.15 ## lowess span
POS_GROUP = "StopGain" ## positive group in EventType 
NEG_GROUP = "Synonymous" ## negative group in EventType
TYPE_COLUMN = "EventType" ## Type column name default EventType

MM_LAMBDA = "0.2" ## "data", or "0.1", or "0.2"
MM_OPTION = "normal" ## normal or robust
MM_REPORT = TRUE ##TRUE FALSE
