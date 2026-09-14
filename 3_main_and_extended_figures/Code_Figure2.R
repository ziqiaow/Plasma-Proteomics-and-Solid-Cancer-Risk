# Required packages
library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggbreak)

folder <- paste0("./collaboration/ARIC_Cancer/all/original/results/")
in.dir = "./collaboration/ARIC_Cancer/all/original/"


#Model 2
# Get list of all .RData files
rdata_files <- list.files(folder, pattern = "\\.RData$", full.names = TRUE)
rdata_files = rdata_files[c(1,3:8)]
# Extract cancer name (after second-to-last underscore, remove trailing _anno)
cancer_names <-  c("Bladder","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate")
#sub("^.*?_.*?_.*?_.*?_(.*?)_.*$", "\\1", rdata_files)


load(rdata_files[1])


res_final <- res_final %>% mutate(
  chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")

) %>%
  # Drop rows where CHR contains "chr" or is not numeric
  filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
  mutate(
    # Convert remaining CHR to numeric
    chr = as.numeric(chr),
    transcription_start_site  = as.numeric(transcription_start_site)
  )

res_final$cancer = cancer_names[1]
tmp = res_final

for(i in 2:length(rdata_files)){
  
  load(rdata_files[i])
  
  
  res_final <- res_final %>% mutate(
    chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")
    
  ) %>%
    # Drop rows where CHR contains "chr" or is not numeric
    filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
    mutate(
      # Convert remaining CHR to numeric
      chr = as.numeric(chr),
      transcription_start_site  = as.numeric(transcription_start_site)
    )
  
  res_final$cancer = cancer_names[i]
  tmp = rbind(tmp,res_final)
  
  
  
}
res_final = tmp
save(res_final,file=paste0(in.dir,"summary_stats_combined_manhattan.RData"))


  

#Model 1
folder <- paste0("./collaboration/ARIC_Cancer/all/original/results/")
in.dir = "./collaboration/ARIC_Cancer/all/original/"

# Get list of all .RData files
rdata_files <- list.files(folder, pattern = "_smk_anno\\.RData$", full.names = TRUE)
#rdata_files = rdata_files[c(10:18)]
# Extract cancer name (after second-to-last underscore, remove trailing _anno)
cancer_names <-  c("Bladder","Colon","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate","Rectal")
#sub("^.*?_.*?_.*?_.*?_(.*?)_.*$", "\\1", rdata_files)


load(rdata_files[1])


res_final <- res_final %>% mutate(
  chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")
  
) %>%
  # Drop rows where CHR contains "chr" or is not numeric
  filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
  mutate(
    # Convert remaining CHR to numeric
    chr = as.numeric(chr),
    transcription_start_site  = as.numeric(transcription_start_site)
  )

res_final$cancer = cancer_names[1]
tmp = res_final

for(i in 2:length(rdata_files)){
  
  load(rdata_files[i])
  
  
  res_final <- res_final %>% mutate(
    chr = str_replace_all(chromosome_name, regex("^chr", ignore_case = TRUE), "")
    
  ) %>%
    # Drop rows where CHR contains "chr" or is not numeric
    filter(!chromosome_name %in% c("X", "Y", "X,Y")) %>%
    mutate(
      # Convert remaining CHR to numeric
      chr = as.numeric(chr),
      transcription_start_site  = as.numeric(transcription_start_site)
    )
  
  res_final$cancer = cancer_names[i]
  tmp = rbind(tmp,res_final)
  
  
  
}
res_final = tmp
save(res_final,file=paste0(in.dir,"summary_stats_combined_manhattan_minimal_smk.RData"))



#make the plot
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_manhattan.RData")
dat <- res_final %>%
  mutate(CHR = chr, BP = transcription_start_site , P = as.numeric(res_final$`Pr(>|z|)`)) %>%
  arrange(CHR, BP)


library(readr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(latex2exp)

thres <- 0.05

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 7),
  text = element_text(size = 6)
  # axis.title.x = element_text(size = 10),
  # axis.text.x = element_text(size = 8),
  # axis.title.y = element_text(size = 10),
  # axis.text.y = element_text(size = 8),
  # legend.title = element_text(size = 10)
  # legend.text = element_text(size = 8)
)

dat$CHR <- factor(dat$CHR, c(1:22))

# We will use ggrepel for the annotation
library(ggrepel)

# Prepare the dataset
dat <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)



nCHR <- 22

# Prepare X axis
axis.set <- dat %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

sig <- dat$fdr < 0.05
dat$`Entrez Gene Name`[which(sig == T)]
# [1] "DRAXIN"  "FCMR"    "PCYOX1"  "ARL5A"   "POMGNT2" "ALCAM"   "ACP3"    "NLGN1"   "IGFBP7"  "KLKB1"   "TGFBI"  
# [12] "HAVCR1"  "APOM"    "IGFBP3"  "MSR1"    "EBAG9"   "PSAT1"   "AIF1L"   "MMRN2"   "IGF2"    "ARFIP2"  "IL18BP" 
# [23] "MMP7"    "MMP7"    "A2M"     "FAM234B" "LRP1"    "IGF1"    "NRXN3"   "NRXN3"   "IGFALS"  "ADGRG1"  "KRT16"  
# [34] "ADAM11"  "CHST9"   "PLAUR"   "KLK3"    "RSPO4"   "WFDC2"   "WFDC8"  

label <- dat$`Entrez Gene Name`[which(sig == T)] #MMP7 appears to be significant in both kidney and liver cancer
#NRXN3 is significant in liver, but measured by 2 aptamers

dat$ID <- paste0(dat$`Entrez Gene Name`,"\n(",dat$cancer,")")

myColors <- c("#eeee00","#cc9955","#006600","#ff00bb","#0000ff","#995522","#33cccc")#gtex.colors$V2
names(myColors) <- c("Bladder","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate")


dat$cancer[!sig & (dat$CHR %in% ((1:11)*2))] <- "black"
dat$cancer[!sig & (dat$CHR %in% ((1:11)*2-1))] <- "grey"
myColors <- c(myColors,"#252525","#969696")
names(myColors)[length(myColors)-1] <- "black"
names(myColors)[length(myColors)] <- "grey"

dat$point_alpha <- 0
dat$point_alpha[sig] <- 1

labels_df = dat[which(sig == T),]
labels_df <- data.frame(label=labels_df$ID,
                        logP=-log10(labels_df$P),
                        BPcum=labels_df$BPcum,
                        CHR=labels_df$CHR,
                        cancer = labels_df$cancer)
labels_df <- labels_df[order(labels_df$BPcum),]
labels_df.original = labels_df

dat <- rbind(dat[!sig,],dat[sig,])

set.seed(09202025)
tmp1 <- which( (dat$cancer %in% c("grey","black")) & (dat$P > thres) )
tmp1 <- sample(tmp1, 4955)
tmp2 <- which(dat$P <= thres)
tmp3 <- which(!(dat$cancer %in% c("grey","black")))
dat_plot <- dat[c(tmp1, tmp2, tmp3),]

M_eff = 2865.267
p.original <- 0.05/M_eff

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 7),
  text = element_text(size = 6)
  # axis.title.x = element_text(size = 10),
  # axis.text.x = element_text(size = 8),
  # axis.title.y = element_text(size = 10),
  # axis.text.y = element_text(size = 8),
  # legend.title = element_text(size = 10)
  # legend.text = element_text(size = 8)
)


p1 <- ggplot(dat_plot, aes(x = BPcum, y = -log10(P), 
                           color = as.factor(cancer), size = -log10(P))) +
  geom_point(aes(alpha = point_alpha), size=0.8) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat$BPcum),max(dat$BPcum))) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, 60 )) +
  scale_y_continuous(breaks =c(0,5,10,20,40,60),limits = c(0,70))  +#trans = "sqrt",
  scale_color_manual(name = "gtex.colors", values = myColors)+
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.original),
             linetype='dashed', col="black", size=0.3) +
  # scale_y_break(c(50, 55))  +# creates a gap in the y-axis
  guides(color = F, alpha = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 6, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 6, vjust = 0.5),
    axis.title = element_text(size=7),
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7)
  ) + 
  ggrepel::geom_label_repel(data = labels_df,#[-c(19:20,22:24,36:39),],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 1.5, segment.size = 0.2,
                            point.padding = 0.3,  segment.linetype = 1,  direction = "both",      # allow sideways placement
                            nudge_y = 1,             # slight vertical offset
                            # segment.curvature = -0.2,
                            # segment.ncp = 3,
                            # segment.angle = 20,
                            # segment.curvature = -1e-20,
                            #arrow = arrow(length = unit(0.015, "npc")),
                            min.segment.length = 0, force = 1,
                            box.padding = 0.5,max.overlaps = Inf)+ 
  My_Theme



## cancer color legends

tmp <- ggplot(dat_plot[!(dat_plot$cancer %in% c("black","grey")), ], aes(x = BPcum, y = -log10(P), 
                                                                         color = as.factor(cancer))) +
  geom_point() + 
  scale_color_manual(name = "Cancer Type", values = myColors)+
  theme_minimal() +
  theme(
    legend.key.size = unit(5, "mm"),
    legend.text=element_text(size=10),
    legend.title = element_text(size=10),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+ 
  My_Theme+
  guides(color=guide_legend(ncol = 1))

p2 <- as_ggplot(get_legend(tmp))


#Model 1 plot
load("./collaboration/ARIC_Cancer/all/original/summary_stats_combined_manhattan_minimal_smk.RData")
dat <- res_final %>%
  filter(!res_final$cancer %in% c("Colon", "Rectal")) %>%
  mutate(CHR = chr, BP = transcription_start_site , P = as.numeric(`Pr(>|z|)`)) %>%
  arrange(CHR, BP)


library(readr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(latex2exp)

thres <- 0.05

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 7),
  text = element_text(size = 6)
  # axis.title.x = element_text(size = 10),
  # axis.text.x = element_text(size = 8),
  # axis.title.y = element_text(size = 10),
  # axis.text.y = element_text(size = 8),
  # legend.title = element_text(size = 10)
  # legend.text = element_text(size = 8)
)

dat$CHR <- factor(dat$CHR, c(1:22))

# We will use ggrepel for the annotation
library(ggrepel)

# Prepare the dataset
dat <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)



nCHR <- 22

# Prepare X axis
axis.set <- dat %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

sig <- dat$fdr < 0.05
dat$`Entrez Gene Name`[which(sig == T)]
# [1] "DRAXIN"  "FCMR"    "PCYOX1"  "ARL5A"   "POMGNT2" "ALCAM"   "ACP3"    "NLGN1"   "IGFBP7"  "KLKB1"   "TGFBI"  
# [12] "HAVCR1"  "APOM"    "IGFBP3"  "MSR1"    "EBAG9"   "PSAT1"   "AIF1L"   "MMRN2"   "IGF2"    "ARFIP2"  "IL18BP" 
# [23] "MMP7"    "MMP7"    "A2M"     "FAM234B" "LRP1"    "IGF1"    "NRXN3"   "NRXN3"   "IGFALS"  "ADGRG1"  "KRT16"  
# [34] "ADAM11"  "CHST9"   "PLAUR"   "KLK3"    "RSPO4"   "WFDC2"   "WFDC8"  

label <- dat$`Entrez Gene Name`[which(sig == T)] #MMP7 appears to be significant in both kidney and liver cancer
#NRXN3 is significant in liver, but measured by 2 aptamers

dat$ID <- paste0(dat$`Entrez Gene Name`,"\n(",dat$cancer,")")

myColors <- c("#eeee00","#cc9955","#006600","#ff00bb","#0000ff","#995522","#33cccc")#gtex.colors$V2
names(myColors) <- c("Bladder","Colorectal","Kidney","Liver","Lung","Pancreatic","Prostate")


dat$cancer[!sig & (dat$CHR %in% ((1:11)*2))] <- "black"
dat$cancer[!sig & (dat$CHR %in% ((1:11)*2-1))] <- "grey"
myColors <- c(myColors,"#252525","#969696")
names(myColors)[length(myColors)-1] <- "black"
names(myColors)[length(myColors)] <- "grey"

dat$point_alpha <- 0
dat$point_alpha[sig] <- 1

labels_df = dat[which(sig == T),]
labels_df <- data.frame(label = labels_df$ID,
                        logP = -log10(labels_df$P),
                        BPcum = labels_df$BPcum,
                        CHR = labels_df$CHR)
labels_df <- labels_df[order(labels_df$BPcum),]
#label those in the original analysis
id = match(labels_df.original$label,labels_df$label)
id = id[complete.cases(id)]
#and label those with -log10(P) > 20
M_eff = 2865.267
p.original <- 0.05/M_eff
#id2 = which(labels_df$logP > -log10(p.original))
id2 = which(labels_df$logP > 5)
id = unique(c(id,id2))
labels_df = labels_df[id2,] #only plot those top significant ones
which(duplicated(labels_df$label))
#because the aptamer measured ICAM5 for lung cancer twice, we only label one
labels_df = labels_df[-26,]

dat <- rbind(dat[!sig,],dat[sig,])

set.seed(09202025)
tmp1 <- which( (dat$cancer %in% c("grey","black")) & (dat$P > thres) )
tmp1 <- sample(tmp1, 4955)
tmp2 <- which(dat$P <= thres)
tmp3 <- which(!(dat$cancer %in% c("grey","black")))
dat_plot <- dat[c(tmp1, tmp2, tmp3),]


My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 7),
  text = element_text(size = 6)
  # axis.title.x = element_text(size = 10),
  # axis.text.x = element_text(size = 8),
  # axis.title.y = element_text(size = 10),
  # axis.text.y = element_text(size = 8),
  # legend.title = element_text(size = 10)
  # legend.text = element_text(size = 8)
)


p3 <- ggplot(dat_plot, aes(x = BPcum, y = -log10(P), 
                           color = as.factor(cancer), size = -log10(P))) +
  geom_point(aes(alpha = point_alpha), size=0.8) +
  scale_alpha_continuous(range = c(0.3, 1)) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat$BPcum),max(dat$BPcum))) +
  #scale_y_continuous(expand = c(0,0), limits = c(0, 60 )) +
  scale_y_reverse(expand=c(0,0),breaks =c(0,5,10,20,40,60),limits = c(70,0))+
  #scale_y_continuous(breaks =c(5,10,20,40,60,80,100))  +#trans = "sqrt",
  scale_color_manual(name = "gtex.colors", values = myColors)+
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.original),
             linetype='dashed', col="black", size=0.3) +
  # scale_y_break(c(50, 55))  +# creates a gap in the y-axis
  guides(color = F, alpha = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 6, vjust = 0.5,color = "transparent"),
    axis.text.y = element_text(angle = 0, size = 6, vjust = 0.5),
    axis.title = element_text(size=7),
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7)
  ) + 
  ggrepel::geom_label_repel(data = labels_df,#[-c(19:20,22:24,36:39),],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 1.5, segment.size = 0.2,
                            point.padding = 0.3,  segment.linetype = 1,  direction = "both",      # allow sideways placement
                           # nudge_y = 1,             # slight vertical offset
                            # segment.curvature = -0.2,
                            # segment.ncp = 3,
                            # segment.angle = 20,
                            # segment.curvature = -1e-20,
                            #arrow = arrow(length = unit(0.015, "npc")),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.6,max.overlaps = Inf)+
  My_Theme


p <- ggarrange(ggarrange(p1, p3,
                         nrow = 2, labels = c("a", "b"),
                         heights = c(0.5,0.5)),
               p2,
               ncol = 2,
               labels = c(NA, NA),
               widths = c(0.85,0.15)
)

ggsave(filename=paste0("Figure2.pdf"),
       plot=p, device="pdf",
       path=in.dir,
       width=220, height=150, units="mm", dpi=320)

ggsave(filename=paste0("Figure2.png"),
       plot=p, device="png",
       path=in.dir,
       width=220, height=150, units="mm", dpi=320)

