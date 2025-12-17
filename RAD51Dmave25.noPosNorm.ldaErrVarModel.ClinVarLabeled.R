#########################################
##  2025 Analysis of RAD51D MAVE Data  ##
##  Last Modified:    07/22/25 by ESI  ##
#########################################


library(knitr)
knit("RAD51Dmave25.noPosNorm.ldaErrVarModel.ClinVarLabeled.Rtex")

uvDB1<-read.delim("aggregated_RAD51D_E1_10_no3_SNV.tsv")
uvDB2<-read.delim("aggregated_RAD51D_E1_10_no3_V3.tsv")

uvDB<-read.delim("aggregated_RAD51D_E1_9_V3.tsv")
uvDB<-read.delim("aggregated_E4_8_SNV.tsv")

uvDB<-read.delim("aggregated_RAD51D_E1_10_SNV_9_12_2025.tsv")

35119543_T_G (E1 aa24); 35106459_G_A (E6 aa168); 35103282_G_A (E8, aa237)
no.label<-c("35119543_T_G","35106459_G_A","35103282_G_A")
##(2): here is the list of variants of E3 that are missing from the VarCall results,

lost<-c("35118520_G_C","35118520_G_A","35118520_G_T",
        "35118568_G_A","35118568_G_C","35118568_G_T",
        "35118616_G_T","35118616_G_C","35118616_G_A")
## Verify that the lost variants are now included (they are, among 160 new).
oldtabl<-read.csv("oldMAVEpostProbs.csv")
table(rownames(oldtabl) %in% rownames(tabl))
plot(oldtabl$logBF,tabl[rownames(oldtabl),"logBF"])
dim(oldtabl)
dim(tabl)


colnames(uvDB)[c(1:7,9:21,26,29,30,31,32)]  ##annot
colnames(uvDB)[-c(1:7,9:21,26,29,30,31,32)] ## measurement
## Annotation Structure. One Row Per Unique Variant
annot<-uvDB[,c(1:7,9:21,26,29,30,31,32)]
dim(annot)
annot<-unique(annot)
dim(annot)
length(unique(annot$uPOS))
rownames(annot)<-annot$uPOS
## Measurement Structure. One Row Per Measurement.
x<-uvDB[,-c(1:7,9:21,26,29,30,32,33)]
## Day by Rep Sets:
uv14.1<-x[x$Time=="D14" & x$Rep=="1",]
uv14.2<-x[x$Time=="D14" & x$Rep=="2",]
uv14.3<-x[x$Time=="D14" & x$Rep=="3",]
uv5.1<-x[x$Time=="D5" & x$Rep=="1",]
uv5.2<-x[x$Time=="D5" & x$Rep=="2",]
uv5.3<-x[x$Time=="D5" & x$Rep=="3",]
uv0.1<-x[x$Time=="Lib" & x$Rep=="1",]
uv0.2<-x[x$Time=="Lib" & x$Rep=="2",]
uv0.3<-x[x$Time=="Lib" & x$Rep=="3",]
if (nrow(uv14.1)==length(unique(uv14.1$uPOS))) rownames(uv14.1)<-uv14.1$uPOS
if (nrow(uv14.2)==length(unique(uv14.2$uPOS))) rownames(uv14.2)<-uv14.2$uPOS
if (nrow(uv14.3)==length(unique(uv14.3$uPOS))) rownames(uv14.3)<-uv14.3$uPOS
if (nrow(uv5.1)==length(unique(uv5.1$uPOS))) rownames(uv5.1)<-uv5.1$uPOS
if (nrow(uv5.2)==length(unique(uv5.2$uPOS))) rownames(uv5.2)<-uv5.2$uPOS
if (nrow(uv5.3)==length(unique(uv5.3$uPOS))) rownames(uv5.3)<-uv5.3$uPOS
if (nrow(uv0.1)==length(unique(uv0.1$uPOS))) rownames(uv0.1)<-uv0.1$uPOS
if (nrow(uv0.2)==length(unique(uv0.2$uPOS))) rownames(uv0.2)<-uv0.2$uPOS
if (nrow(uv0.3)==length(unique(uv0.3$uPOS))) rownames(uv0.3)<-uv0.3$uPOS
## Start with annotation structure (one row per observed variant)
uvDB<-annot
## Add in Day- and Replicate-Specfic Read Counts:
uvDB$R1_D14<-rep(NA,nrow(uvDB))
uvDB$R2_D14<-rep(NA,nrow(uvDB))
uvDB$R3_D14<-rep(NA,nrow(uvDB))
uvDB$R1_D5<-rep(NA,nrow(uvDB))
uvDB$R2_D5<-rep(NA,nrow(uvDB))
uvDB$R3_D5<-rep(NA,nrow(uvDB))
uvDB$R1_lib<-rep(NA,nrow(uvDB))
uvDB$R2_lib<-rep(NA,nrow(uvDB))
uvDB$R3_lib<-rep(NA,nrow(uvDB))
uvDB$libRepA<-rep(NA,nrow(uvDB))
uvDB$libRepB<-rep(NA,nrow(uvDB))
##
uvDB[rownames(uv14.1),"R1_D14"]<-uv14.1$EventCount
uvDB[rownames(uv14.2),"R2_D14"]<-uv14.2$EventCount
uvDB[rownames(uv14.3),"R3_D14"]<-uv14.3$EventCount
##
uvDB[rownames(uv5.1),"R1_D5"]<-uv5.1$EventCount
uvDB[rownames(uv5.2),"R2_D5"]<-uv5.2$EventCount
uvDB[rownames(uv5.3),"R3_D5"]<-uv5.3$EventCount
##
uvDB[rownames(uv0.1),"R1_lib"]<-uv0.1$EventCount
uvDB[rownames(uv0.2),"R2_lib"]<-uv0.2$EventCount
uvDB[rownames(uv0.3),"R3_lib"]<-uv0.3$EventCount
##
uvDB[rownames(uv0.1),"libRepA"]<-uv0.1$Rep
uvDB[rownames(uv0.2),"libRepB"]<-uv0.2$Rep
##
lost<-c("35118520_G_C","35118520_G_A","35118520_G_T",
        "35118568_G_A","35118568_G_C","35118568_G_T",
        "35118616_G_T","35118616_G_C","35118616_G_A")
table(uvDB$uPOS %in% lost); table(unique(uvDB$uPOS) %in% lost)
uvDB$Exon<-factor(uvDB$Exon)
cnames<-colnames(uvDB)
uvDB<-uvDB[,!(cnames %in% c("sample_id","Rep","Time"))]





















uvDB[(!is.na(uvDB$AApos))&(uvDB$Exon=="E3")&(uvDB$Rep==1)&(uvDB$AApos==81),]
table(apply(is.na(uvDB),1,sum))
table((!dropR1)&(!dropR2)&(!dropR3))

tbl<-table(uvDB$variant[uvDB$Exon=="E3"],uvDB$replicate[uvDB$Exon=="E3"])
table(uvDB$Exon,uvDB$replicate) ## after drop
   R1  R2  R3
  E1  301 296 298
  E10 306 306 306
  E2  241 239 240
  E3  395 395 395
  E4  303 303 303
  E5  461 462 460
  E6  345 345 345
  E7  330 330 330
  E8  270 270 270
  E9  507 521 519

uvDB$uPOS<-paste0(uvDB$POS,"_",uvDB$REF,"_",uvDB$ALT)
dim(uvDB)
length(unique(uvDB$uPOS))
table(uvDB$Exon,uvDB$Rep,useNA="always")
table(uvDB$Time,uvDB$Rep,useNA="always")

  pdf("ResidualQQ.pdf",width=8,height=8)
  qqplot2(x=expected.stresids,tdf=101,
          titl="Normal QQ Plot of Expected Standardized Residuals")
  dev.off()

rg<-rgamma(n=100000,25,1)
rg<-(1/abs(rg))
rt<-rnorm(n=100000,0,rg)
hist(rt,nclass=500)
qqplot2(x=rt[1:1000],titl="",tdf=24)
24*24*var(rt)  ## -> sd=1/24
24/23 ## t scale=1 variance (=df/(df-2))
## sp prior on sigma^2_eta looks like 1/2 t(0,1/24)

uvDB<-uvDB[uvDB$EventType %in% c("Missense","Synonymous","StopGain"),]

tabl.lda<-utils::read.csv("../rad51d.mac.mave25/MAVEpostProbs.csv",row.names="variant")
table(rownames(tabl)==rownames(tabl.lda))
plot(tabl$lPostOdds,tabl.lda$lPostOdds)


pdf("ModelComparison.pdf",height=10,width=8.5)
plot(tabl$lPostOdds,tabl.lda$lPostOdds,las=1,
     xlab="Standard Model",ylab="Error Variance Model",
     main="logPostOdds Standard Model VS Error Variance Model")
abline(a=0,b=1,lwd=2,col=2)
abline(h=0,lty=3)
abline(v=0,lty=3)

tabl$eta.ae<-((tabl$eta.ul - tabl$eta.ll)/2)
tabl.lda$eta.ae<-((tabl.lda$eta.ul - tabl.lda$eta.ll)/2)
plot(tabl$eta.ae,tabl.lda$eta.ae,las=1,
     xlab="Standard Model",ylab="Error Variance Model",
     main="Variant Effects Error Allowance\n Standard Model VS Error Variance Model")
abline(a=0,b=1,lwd=2,col=2)
abline(h=0,lty=3)
abline(v=0,lty=3)

plot(tabl$eta,tabl.lda$eta,las=1,
     xlab="Standard Model",ylab="Error Variance Model",
     main="Variant Effects (Eta's)\n Standard Model VS Error Variance Model")
abline(a=0,b=1,lwd=2,col=2)
abline(h=0,lty=3)
abline(v=0,lty=3)
dev.off()







