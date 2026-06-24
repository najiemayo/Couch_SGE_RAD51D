## this is for the manuscript figure
## 
library(readxl)
library(ggplot2)
library(ggh4x)
library(forcats)
library(readr)
library(tidyr)
library(dplyr)


resdir <- 'respath'
allfs <- read_excel("For figure 2_092325.xlsx", sheet = "SNV")
allfs <- allfs %>% mutate(uPOS = paste(POS_hg38, REF, ALT, sep = "_")) %>% mutate(EventType := !!sym('variant class'))%>%
  mutate(GRCh38Location = POS_hg38)
if(length(which(is.na(allfs$uPOS))) > 0) allfs <- allfs[-which(is.na(allfs$uPOS)), ] ## remove empty rows
allfs$Eta <- as.numeric(allfs$"Model based functional score")
allfs$EventType[which(allfs$EventType %in% c("missense", "Start codon missense", "Nonsense_lost"))] <- "Missense"
allfs$EventType[which(allfs$EventType %in% c("5' UTR", "3' UTR"))] <- "UTR"

##Please merge the "Start codon missense" and "Nonsense_lost" into the "missense" category.
#Please merge the "5' UTR" and "3' UTR" into a new category ‘UTR”.


allfs <- allfs %>% filter(EventType %in% c("Canonical splice", "Intronic", "Missense", "Silent", "Nonsense", "UTR"))

allfs$EventType <- factor(allfs$EventType, levels = c("Canonical splice", "Intronic", "Missense", "Silent", "Nonsense", "UTR"))
allfs$AApos <- as.numeric(str_extract(allfs$AA, "\\d+"))
allfs$Exon <- as.factor(allfs$Exon)

allfs$functional.prediction <- factor(allfs$`Functional category`, levels = c("P strong", "P moderate", "P supporting", "VUS", "B supporting", "B moderate", "B strong"))

## Figure 2A

p <- dt2plot %>%
  ggplot( aes(x=Eta, fill=EventType)) +  
  geom_histogram( binwidth = 0.025, position="stack") +
  scale_fill_manual(values= mycol5[(dt2plot$EventType)]) +
  ggtitle ( "") + theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        legend.title=element_blank(), 
                        legend.position = "inside", legend.position.inside = c(.3, .8)) + scale_y_continuous(expand = c(0, 0)) +
  xlab("Model based functional score") + ylab("Number of variants")
p


ggsave(paste0(resdir, "/Figure2A_eta_hist.tif"),plot = p, width = 5, height = 5, units = "in", device = "tiff", dpi = 300, limitsize = FALSE)

## F2B
### scatter plot
library(ggbreak)
posb <- allfs %>% select(Exon, GRCh38Location)%>% group_by(Exon) %>% slice_max(order_by = GRCh38Location, n =1) %>% unique() %>% arrange(GRCh38Location)
posb2 <- allfs %>% select(Exon, GRCh38Location)%>% group_by(Exon) %>% slice_min(order_by = GRCh38Location, n =1) %>% unique() %>% arrange(GRCh38Location)
library(vctrs)
gap <- 20
mybreak2 <- vec_interleave(posb$GRCh38Location[1:9], posb2$GRCh38Location[2:10]) 
mybreak3 <- c(posb2$GRCh38Location[1]+gap,
              vec_interleave(posb$GRCh38Location[1:9]-gap, posb2$GRCh38Location[2:10]+gap) ,
              posb$GRCh38Location[10]-gap)
mybreak4 <- c(posb2$GRCh38Location[1],
              vec_interleave(posb$GRCh38Location[1:9], posb2$GRCh38Location[2:10]) ,
              posb$GRCh38Location[10])

ylabsize <- 10; xlabsize <- 10; legdsize <- 12
dt2plot <- allfs %>% arrange(GRCh38Location) # GRCh38Location

p <- ggplot(dt2plot, aes(x = GRCh38Location, y = Eta, colour = EventType)) +
  geom_point() + 
  scale_colour_manual(values = mycol5) +
  xlab("GRCh38 Genomic Coordinate") + ylab("Functional score")+
  theme_classic()+
  
  scale_x_break(breaks = mybreak2, ticklabels = NULL, expand = c(0.01, 0.01)) +
  scale_x_reverse() + 
  theme( axis.ticks.x=element_blank(),
         axis.text.x= element_blank(),
         axis.text.y=element_text(size = ylabsize*1.4),
         axis.title.y = element_text(size = ylabsize*1.5),
         axis.title.x = element_text(size = ylabsize*2))+
  annotate("text", x = mybreak3, y = 0, label = mybreak4, angle = 90, hjust = 0, size = ylabsize *0.5)+
  annotate(geom = "text", x = posb$GRCh38Location, y = 1.5, hjust =0, label = paste0(posb$Exon, " "), size = ylabsize*0.5)+
  coord_cartesian(ylim = c(-0.3, 1.8), clip = "off")+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.5, fill = NA))+
  theme(legend.position = "right", legend.title = element_blank(), legend.text=element_text(size = legdsize))

p

ggsave(paste0(resdir, "/Figure2B_eta_byexon_POS.tif"),plot = p, width = 20, height = 5, units = "in", device = "tiff", dpi = 300, limitsize = FALSE)



### F2C
## barplot
ylabsize <- 10; xlabsize <- 10; legdsize <- 12
fs <- allfs

dt2plot <- fs %>% group_by(EventType, functional.prediction) %>% count()
df <- dt2plot %>% group_by(EventType) %>% mutate(percent = prop.table(n))
df$total <- ave(df$n, df$EventType, FUN = sum)
df$cat <- paste0(df$EventType, ", n=", df$total)
unique(df$cat)

table(allfs$EventType)
df$cat <- factor(df$cat, levels = c("Silent, n=728", "Intronic, n=427", "Missense, n=2015", "Canonical splice, n=108",  "Nonsense, n=93", "UTR, n=60"))
# plot horizontal bar chart with percentage labels
mycol7 <- c("#800000","#FFA500", "#FFD700", "gray",  "#00b4c5" ,"#0073e6", "#0059a3") ## "#CC5500",  "#800020",  
names(mycol7) <- levels(df$functional.prediction)


p <- ggplot(df, aes(x = cat, y = percent*100, fill = functional.prediction)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values =  mycol7[(df$functional.prediction)]) +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("") + ylab("Percentage %")+
  theme_classic()+ #theme(text = element_text(family = "Calibri"))+
  theme(axis.text=element_text(size=ylabsize*1.8),
        axis.title=element_text(size=ylabsize*1.8))+
  theme(legend.position = "right", legend.title = element_blank(), legend.text=element_text(size = legdsize))
p


ggsave(paste0(resdir, "/Figure2C_barplot_percentage_all.tif"),plot = p, width = 9, height = 5, units = "in", device = "tiff", dpi = 300, limitsize = FALSE)


### F2D
## barplot
allfs$Exon <- factor(allfs$Exon, levels = c("E1",  "E2", "E3", "E4",  "E5",  "E6",  "E7",  "E8",  "E9",  "E10"))
dt2plot <- allfs  %>% group_by(Exon, EventType, functional.prediction) %>% count() 
df <- dt2plot %>% group_by(Exon, EventType) %>% mutate(percent = prop.table(n))
df$cat <- factor(df$EventType, levels = c("Silent", "Intronic", "Missense", "Canonical splice",  "Nonsense", "UTR"))

p <- ggplot(df, aes(x = Exon, y = percent*100, fill = functional.prediction)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = mycol7[(df$functional.prediction)]) +
  facet_wrap(.~fct_rev(cat), ncol = 1)+
  xlab("Exon") + ylab("Percentage %")+
  theme(strip.text.x = element_text(size = 65))+
  theme_classic()+ #theme(text = element_text(family = "Calibri"))+
  theme(axis.text=element_text(size=ylabsize),
        axis.title=element_text(size=ylabsize))+
  theme(legend.position = "right", legend.title = element_blank(), legend.text=element_text(size = legdsize)) 
p


ggsave(paste0(resdir, "/Figure2E_barplot_percentage_byexon.tif"),plot = p, width = 10, height = 5, units = "in", device = "tiff", dpi = 300, limitsize = FALSE)


## Figure4
## Figure 4A
# 
## heatmap
# # combining by AApos, heatmap all exons
dt2plot0 <- read_excel("RAD51D MAVE heatmap source data_092225updated.xlsx")


dt2plot0 <- dt2plot0 %>% mutate(AltAA = substr(AA, nchar(AA), nchar(AA))) %>% mutate(AApos = as.numeric(gsub("\\D", "", AA)))
unique(dt2plot0$AltAA)
dt2plot0$AltAA[which(dt2plot0$AltAA == "X")] <- "stop"

# dt2plot0$AltAA <- factor(dt2plot0$AltAA, levels = c("Int", "A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N", "Q", "C", "G", "P", "R", "H",  "K","D",  "E", "stop"))
dt2plot0$AltAA <- factor(dt2plot0$AltAA, levels = c( "A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N", "Q", "C", "G", "P", "R", "H",  "K","D",  "E", "stop"))
dt2plot0$functional.prediction <- factor(dt2plot0$`Functional category`, levels = c("P strong", "P moderate","P supporting", "B strong","B moderate", "B supporting", "VUS"))

# pick the stronger level

dt2plot <- dt2plot0 %>% select(AltAA, AApos, functional.prediction) %>% group_by(AApos, AltAA) %>% arrange((functional.prediction)) %>% slice(1) %>% ungroup() %>% as.data.frame()


## fill in the gaps among AA
gaps <- setdiff(1:max(dt2plot$AApos), unique(dt2plot$AApos))
add2 <- data.frame(AltAA = rep(levels(dt2plot$AltAA), length(gaps)), functional.prediction = "gap", AApos = rep(gaps, each = nlevels(dt2plot$AltAA)))

dt2plot <- rbind(dt2plot, add2)

ylabs2 <- rev(c("A", "V", "I", "L", "M", "F", "Y", "W", "S", "T", "N", "Q", "C", "G", "P", "R", "H",  "K","D",  "E", "stop"))
dt2plot$AltAA <- droplevels(dt2plot$AltAA)
dt2plot$AApos <- factor(dt2plot$AApos)

dt2plot$functional.prediction <- droplevels(dt2plot$functional.prediction)
dt2plot$functional.prediction <- factor(dt2plot$functional.prediction, levels = c("P strong", "P moderate", "P supporting", "VUS", "B supporting", "B moderate", "B strong", "gap"))


AAposs <- unique(dt2plot$AApos)
AAposs <- AAposs[order(AAposs)]
mybreaks <- AAposs[c(seq(50, length(AAposs), 50))]
mybreaks <- c(1, mybreaks,  AAposs[length(AAposs)])

##mycol8 <- c("#800000","#FFA500", "#FFD700", "gray",  "#00b4c5" ,"#0073e6", "#0059a3", "white") ## "#CC5500",  "#800020",  

mycol8 <- c("#B22222","#FFA500", "#FFD700", "gray",  "#00FFFF","#66ddff", "#00b4c5", "white") ## 
names(mycol8) <- levels(dt2plot$functional.prediction)
desired_levels <- levels(dt2plot$functional.prediction)[-length(levels(dt2plot$functional.prediction))]

dt2plot <- dt2plot %>% arrange(AApos)
p1 <- ggplot(dt2plot, aes(x = AApos, y = (AltAA), fill = functional.prediction)) +
  #Specify colours
  scale_fill_manual(values = mycol8[(dt2plot$functional.prediction)], breaks = desired_levels) +
  scale_x_discrete(position = "bottom", name = "Amino Acid Position", expand = c(0,0), breaks = mybreaks, ) +
  geom_tile() +
  theme_classic()+ theme(text = element_text(family = "Calibri")) +
  scale_y_discrete(limits=rev)+ ylab("Amino Acid Substitution")+
  theme(axis.text=element_text(size=ylabsize*1.5),
        axis.title.x = element_text(size = 1.5 * 11),  # X label 1.5× default
        axis.title.y = element_text(size = 1.5 * 11))+ # Y label 1.5× default)+
  theme(legend.position = "right", legend.title=element_blank(), legend.text=element_text(size = legdsize*1.5)) +
  theme(axis.line.y = element_blank(), axis.text.y = element_text(hjust = 0.5, margin = margin(r = 2)), 
        axis.text.x = element_text(vjust = 1)) # remove y axis line, align letters in center
