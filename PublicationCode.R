# Packages #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages(c("tidyverse", "ggtext", "ggpubr",
                   "ggbreak", "ggrepel", "patchwork",
                   "openxlsx"))
library(tidyverse)
library(limma)
library(edgeR) # we used v4.3.2 
library(viridisLite)
library(viridis)
library(ggbreak)
library(ggtext)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(openxlsx)

bw <- theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text = element_text(color = "black"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            legend.text = element_text(color = "black"),
            legend.title = element_text(color = "black"),
            plot.title = element_text(color = "black"),
            plot.subtitle = element_text(color = "black")) 

# Working Directory #####
setwd("C:/Users/birch/Desktop/Grad School Research/Polarella")

# Necessary Files for Full Pipeline #####
# For working files with processed data go to line XXX

# Read count table preprocessed
# Some of the contig IDs are duplicates -> add a column for a working ID
geneIDkNum <- read.table("featureCount_ALL_geneid.txt", sep = "\t", header = TRUE, row.names = "Geneid")
geneIDkNum$Col <- c(1:nrow(geneIDkNum))
forRejoin <- geneIDkNum %>%
  select(c("Chr", "Col"))
geneIDkNumEDGER <- geneIDkNum %>%
  select(-"Chr") %>%
  relocate(last_col(), .before = everything())

# Metadata
metadata13 <- read.table("metadata_etnp18.txt", sep = "\t", header = TRUE, row.names = "Sample")
surface <- c('Cast22_60m','Cast59_30m','Cast50_40m','Cast20_54m',
             'Cast34_60m','Cast52_65m','Cast98_25m','Cast85_20m',
             'Cast93_10m','Cast23_40m','Cast55_120m','Cast6_60m',
             'Cast63_10m','Cast73_16m')
oxycline <- c('Cast39_85m','Cast41_90m','Cast97_40m','Cast80_50m',
              'Cast48_90m','Cast31_90m','Cast47_113m','Cast37_92m',
              "Cast59_900m")
anoxic <- c('Cast35_120m','Cast59_250m','Cast27_275m','Cast18_110m',
            'Cast19_90m', 'Cast72_45m','Cast76_36m','Cast101_80m','Cast82_48m')
station1 <- metadata13 %>%
  filter(Station == "1") %>%
  rownames()
station2 <- metadata13 %>%
  filter(Station == "2") %>%
  rownames()
station3 <- metadata13 %>%
  filter(Station == "3") %>%
  rownames()

# Taxonomy BASTA 
# We are focused on Polarella -> create a clean Polarella file
basta <- read.delim("final.contigs-gmst_50hits_50identity.basta", header = FALSE, sep = "\t") 
bastaClean <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "genus")) %>%
  filter(genus == "Polarella")

# Panther/InterproScan (Annotations)
# From prior analysis, we know the genes we want to focus on
pantherAnno <- read.delim("final.contigs-gmst_PANTHER_simple.txt", header = FALSE, sep = "\t")
pantherAnnoSimp <- pantherAnno %>% 
  select (-c(V2))
colnames(pantherAnnoSimp) <- c("Chr",  "ID", "SignatureDescription", 
                               "CONTevalue","AccessIPR", "Descriptive")
genesOfInterest <- c("PTHR42704", "PTHR21649", "PTHR11002", "PTHR28286", 
                     "PTHR42715", "PTHR11177", "PTHR10353", "PTHR23516", 
                     "PTHR12121",  "PTHR33175", "PTHR11630",
                     "PTHR48085", "PTHR23420", "PTHR23515", "PTHR11814", 
                     "PTHR31632", "PTHR22950", "PTHR10464", "PTHR11101", "PTHR19372", 
                     "PTHR43310", "PTHR30520", "PTHR11730", "PTHR19375", "PTHR45527", 
                     "PTHR11544", "PTHR43208",  "PTHR11709", "PTHR11693",
                     "PTHR12400")
commonName <- as.data.frame(c("RuBisCO", "Chl A/B Binding Protein", 
                              "Carbonic Anhydrase", "Rhodopsin", "Beta-Glucosidase (GH3)",                                       
                              "Chitinase (GH18)", "Glycosyl Hydrolase (GH1)", "Molybdate Transporter",
                              "Carbon Catabolite Repressor P4",                
                              "DNA-Binding Protein HU", "DNA Replicating - MCM",
                              "Cd/Zn-Transporting ATPase",  "Adenosylhomocysteinase", 
                              "High-Affinity NO3- Transporter", "Sulfate Transporter",                                  
                              "Iron Transporter FTH1", "Amino Acid Transporter", "Urea Transporter",
                              "Phosphate Transporter", "Sulfite Reductase *", "Sulfate Transporter YBAR",                     
                              "Formate/Nitrite Transporter", "Ammonium Transporter", "Heat Shock Protein 70kDA",
                              "Nonribosomal Peptide Synthase", "Cold Shock Domain",                   
                              "ABC Transporter", 
                              "Auxiliary Activities 1", "ATP Synthase Gamma Subunit", "Polyphosphate Kinase*"))
colnames(commonName) <- "Name"
commonName$ID <- genesOfInterest
selectedGenes <- pantherAnnoSimp %>%
  inner_join(bastaClean, by = "Chr") %>%
  filter(ID %in% genesOfInterest) %>%
  left_join(commonName)

# Eggnog (Annotations)
# Exploration with this final was done in prior analysis, it goes unused with this file as
# annotation IDs from PATNTHER is what was used for joins later in the code (for no particular reason)
eggnog <- read.delim("final.contigs-gmst.emapper.annotations", header = FALSE, sep = "\t", skip = 5)
eggnogSimp <- eggnog %>%
  select(-c(V2, V5, V6, V7, V9, V19, V18, V19, V20))
colnames(eggnogSimp) <- c("Chr", "Eval", "OrthoScore", "ProteinName", "GO",
                          "EC", "KO", "KO_Path", "KO_Mod", "KO_React", 
                          "KO_rclass", "KO_TC", "EggDescript")

# Cast Profiles from Stations
cast18 <- read.table("NutrientENTP2018/ctd_18.txt", skip =1, sep = "")
cast59 <- read.table("NutrientENTP2018/ctd_59.txt", skip =1, sep = "")
cast93 <- read.table("NutrientENTP2018/ctd_93.txt", skip =1, sep = "")

# Nutrient
No3No2 <- read.xlsx("NO3NO2.xlsx")
Nh4 <- read.xlsx("NH4.xlsx")

# dbcan results (extracellular CAZYs)
dbcan <- read.delim("dbcan_res.txt")
dbcan$Chr <- row.names(dbcan)
filterMe <- c("k137_2263365", "k137_1511638", "k137_5267244", "k137_3849297")
dbcan2 <- dbcan %>%
  filter(Chr %in% filterMe)

# Getting RPKM from EdgeR #####
# PolyA --> This is for DE as polyA is generally better for eukaryotes 
# Factor sample depths
depth <- factor(metadata13$Category_new)
# Construct edgeR object
edOb <- DGEList(geneIDkNumEDGER[,-c(1,2,3,4,5)], group=depth, genes=geneIDkNumEDGER[,c(1,5), drop = FALSE])
# Remove genes with low counts
keep <- rowSums(cpm(edOb) > 0.1) >= 3
table(keep) # how many are TRUE
# FALSE    TRUE 
# 1220150 4777327 
edOb2 <- edOb[keep, , keep.lib.sizes=FALSE]
# Normalization for composition bias
edObNorm <- calcNormFactors(edOb2, method = "TMM")
# Design matrix
design <- model.matrix(~0+depth)
colnames(design) <- levels(depth)
design
# Estimate Dispersion
edObDis <- estimateDisp(edObNorm, design, robust=TRUE)
# Full RPKM Table
rpkmData <- as.data.frame(rpkm(edObDis))
rpkmData$Col <- edObDis$genes$Col
# Save the RPKM Data
write.csv(rpkmData, file = "rpkmData.csv",  row.names = FALSE)
# Estimate Quasi-Likelihood (QL) dispersions
edObFit <- glmQLFit(edObDis, design, robust=TRUE)
head(edObFit$coefficients)

# Differential expression analysis # A vs S #####
diffXaAS <- makeContrasts(AvsS = Anoxic - Surface,
                          levels = design)
# QLFTest = Quasi-Likelihood F Test
res1AS <- glmQLFTest(edObFit, contrast=diffXaAS) 
topTags(res1AS) # view the top DE genes
res1.tableAS <-  topTags(res1AS, n = "Inf")$table
res1.tableAS$Chr <- res1AS$genes$Chr
resTableAS <- as.data.frame(res1.tableAS)
resTableAS <- resTableAS[ ,c(1,3,7)]
colnames(resTableAS) <- c("Col", "LogFC_AvsS", "FDR_AvsS")
# Count how many genes are differentially regulated
diffRegAS <- decideTestsDGE(res1AS)
summary(diffRegAS)
# 1*Anoxic -1*Surface
# Down               2300325
# NotSig             1006309
# Up                 1470693

# Differential expression analysis # O vs S #####
diffXaOS <- makeContrasts(OvsS = Oxycline - Surface,
                          levels = design)
# QLFTest = Quasi-Likelihood F Test
res1OS <- glmQLFTest(edObFit, contrast=diffXaOS) 
topTags(res1OS) # view the top DE genes
res1.tableOS <-  topTags(res1OS, n = "Inf")$table
res1.tableOS$Chr <- res1OS$genes$Chr
resTableOS <- as.data.frame(res1.tableOS)
resTableOS <- resTableOS[ ,c(1,3,7)]
colnames(resTableOS) <- c("Col", "LogFC_OvsS", "FDR_OvsS")
# Count how many genes are differentially regulated
diffRegOS <- decideTestsDGE(res1OS)
summary(diffRegOS)
# 1*Oxycline -1*Surface
# Down                  733048
# NotSig               3376377
# Up                    667902

# Differential expression analysis # O vs A #####
diffXaOA <- makeContrasts(OvsA = Oxycline - Anoxic,
                          levels = design)
# QLFTest = Quasi-Likelihood F Test
res1OA <- glmQLFTest(edObFit, contrast=diffXaOA) 
topTags(res1OA) # view the top DE genes
res1.tableOA <-  topTags(res1OA, n = "Inf")$table
res1.tableOA$Chr <- res1OA$genes$Chr
resTableOA <- as.data.frame(res1.tableOA)
resTableOA <- resTableOA[ ,c(1,3,7)]
colnames(resTableOA) <- c("Col", "LogFC_OvsA", "FDR_OvsA")
# Count how many genes are differentially regulated
diffRegOA <- decideTestsDGE(res1OA)
summary(diffRegOA)
# -1*Anoxic 1*Oxycline # pay attention to the negative here!
# Down                1060115
# NotSig              2143832
# Up                  1573380

# Now we have logFC and FDR for all genes in the dataset 
# Make a big table that contains contig ID
EdgeR_Res <- resTableAS %>%
  left_join(resTableOA) %>%
  left_join(resTableOS) %>%
  left_join(forRejoin)
# Save this result so we can analyze other genes if need be
write.csv(EdgeR_Res, "EdgeR_Res.csv")

# Post-EdgeR Pipeline Read In #####
# Still need to read in files above for the joining of data
rpkmData <- read.csv("rpkmData.csv")
EdgeR_Res <- read.csv("EdgeR_Res.csv")

# Finding amount of DE Polarella genes #####
PolarellaDEG <- EdgeR_Res %>%
  inner_join(bastaClean, by = "Chr") %>%
  select(-X) %>%
  inner_join(rpkmData, by = "Col") %>%
  left_join(pantherAnnoSimp)
PolarellaDEGAS <- PolarellaDEG %>%
  filter(FDR_AvsS < 0.05) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(surface, anoxic)))))
nrow(PolarellaDEGAS)
nrow(PolarellaDEGAS %>%
       filter(LogFC_AvsS > 0))
PolarellaDEGOS <- PolarellaDEG %>%
  filter(FDR_OvsS < 0.05) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(surface, oxycline)))))
nrow(PolarellaDEGOS)
nrow(PolarellaDEGOS %>%
       filter(LogFC_OvsS > 0))
PolarellaDEGOA <- PolarellaDEG %>%
  filter(FDR_OvsA < 0.05) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(oxycline, anoxic)))))
nrow(PolarellaDEGOA)
nrow(PolarellaDEGOA %>%
       filter(LogFC_OvsA > 0))

# Heatmap, LogFC Plot, and LogFC Table Maniplutations #####
# The "col" ID column will be what we use to match to original contigs
# and only keep the columns from the cast
# From our dbcan result, we need to preserve only SOME of the AA1 genes (the ones that are extracellular)
# The ones we are keeping are extracellular. All other extracellular genes are already pulled from prior work!
PolarellaDEG_2 <- PolarellaDEG %>%
  filter(!(ID == 'PTHR11709' & !(Chr %in% dbcan2$Chr))) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(oxycline, surface, anoxic))))) %>%
  arrange(desc(sumRPKM)) %>%
  group_by(ID) %>%
  top_n(1, (sumRPKM)) %>% 
  filter(sumRPKM != 0) %>%
  arrange(desc(sumRPKM)) %>%
  ungroup 
PolarellaDEGAS_2 <- PolarellaDEGAS %>%
  filter(!(ID == 'PTHR11709' & !(Chr %in% dbcan2$Chr))) %>%
  inner_join(selectedGenes) %>%
  arrange(desc(sumRPKM)) %>%
  group_by(ID) %>%
  top_n(1, (sumRPKM)) %>% 
  filter(sumRPKM != 0) %>%
  arrange(desc(sumRPKM)) %>%
  ungroup 
PolarellaDEGOS_2 <- PolarellaDEGOS %>%
  filter(!(ID == 'PTHR11709' & !(Chr %in% dbcan2$Chr))) %>%
  inner_join(selectedGenes) %>%
  arrange(desc(sumRPKM)) %>%
  group_by(ID) %>%
  top_n(1, (sumRPKM)) %>% 
  filter(sumRPKM != 0) %>%
  arrange(desc(sumRPKM)) %>%
  ungroup 
PolarellaDEGOA_2 <- PolarellaDEGOA %>%
  filter(!(ID == 'PTHR11709' & !(Chr %in% dbcan2$Chr))) %>%
  inner_join(selectedGenes) %>%
  arrange(desc(sumRPKM)) %>%
  group_by(ID) %>%
  top_n(1, (sumRPKM)) %>% 
  filter(sumRPKM != 0) %>%
  arrange(desc(sumRPKM)) %>%
  ungroup 

# Keeping contigs shared between respective condition comparisons
Contig_ASvOA <- intersect(PolarellaDEGAS_2$Chr, PolarellaDEGOA_2$Chr)
Contig_ASvOS <- intersect(PolarellaDEGAS_2$Chr, PolarellaDEGOS_2$Chr) 
Contig_OAvOS <- intersect(PolarellaDEGOA_2$Chr, PolarellaDEGOS_2$Chr) 

# Combine Same Contig Ones
AS_selected_OA <- PolarellaDEGAS_2 %>%
  filter(Chr %in% Contig_ASvOA) %>%
  select(-all_of(oxycline))
OA_reduced_AS <- PolarellaDEGOA_2 %>%
  filter(Chr %in% Contig_ASvOA) %>%
  select(all_of(oxycline))
RPKM_sameContig_ASvOA <- bind_cols(AS_selected_OA, OA_reduced_AS)
RPKM_sameContig_ASvOA$Comparison <- "AvS | OvA"
AS_selected_OS <- PolarellaDEGAS_2%>%
  filter(Chr %in% Contig_ASvOS) %>%
  select(-all_of(oxycline))
OS_reduced_AS <- PolarellaDEGOS_2%>%
  filter(Chr %in% Contig_ASvOS) %>%
  select(all_of(oxycline))
RPKM_sameContig_ASvOS <- bind_cols(AS_selected_OS, OS_reduced_AS)
RPKM_sameContig_ASvOS$Comparison <- "AvS | OvS"
RPKM_sameContig <- rbind(RPKM_sameContig_ASvOA, RPKM_sameContig_ASvOS)

# Contigs not shared by respective conditions
RPKM_noSharedContig_AS <- PolarellaDEGAS_2 %>%
  filter(!Chr %in% RPKM_sameContig$Chr)
RPKM_noSharedContig_AS$Comparison <- "AvS"
RPKM_noSharedContig_OS <- PolarellaDEGOS_2 %>%
  filter(!Chr %in% RPKM_sameContig$Chr)
RPKM_noSharedContig_OS$Comparison <- "OvS"
RPKM_noSharedContig_OA <- PolarellaDEGOA_2 %>%
  filter(!Chr %in% RPKM_sameContig$Chr)
RPKM_noSharedContig_OA$Comparison <- "OvA"
RPKM_diffContig <- rbind(RPKM_noSharedContig_AS, RPKM_noSharedContig_OA) %>%
  rbind(RPKM_noSharedContig_OS)

# Combined RPKM for heatmap
RPKM_combined <- rbind(RPKM_diffContig, RPKM_sameContig)

# Adding gene group
transporters <- c("PTHR23515", "PTHR11730", "PTHR43029", "PTHR11814", "PTHR43310",
                  "PTHR31632", "PTHR30520", "PTHR43809", "PTHR11101", "PTHR19372",
                  "PTHR10464", "PTHR22950",  "PTHR23516", "PTHR43208", "PTHR48085", "PTHR12400")
stressresponse <- c("PTHR23410", "PTHR23253", "PTHR10891", "PTHR11630", "PTHR33175",
                    "PTHR19211", "PTHR23420", "PTHR11127", "PTHR31062","PTHR10031",
                    "PTHR12399", "PTHR11620", "PTHR12850", "PTHR43652", "PTHR12121")
phototrophy <- c("PTHR21649", "PTHR42704", "PTHR11002", "PTHR28286", "PTHR11693")
heterotrophy <- c("PTHR38674", "PTHR11177", "PTHR10353", 
                  "PTHR31263", "PTHR12121",  "PTHR42715",
                  "PTHR34939", "PTHR11002","PTHR10587",  "PTHR43400", "PTHR11709")
other <- c()
RPKM_combined_GeneGroup <- RPKM_combined %>%
  mutate(type = case_when(
    ID %in% transporters ~ "Transporters",
    ID %in% stressresponse ~ "Stress Response",
    ID %in% phototrophy ~ "Phototrophy",
    ID %in% heterotrophy ~ "PH"
  )) %>%
  mutate(type = replace_na(type, "Stress Response")) %>%
  unite(New_Name, Name, Comparison, sep = " ") %>%
  mutate(LogFC_AvsS = if_else(grepl("AvS", New_Name), LogFC_AvsS, NA)) %>%
  mutate(LogFC_OvsS = if_else(grepl("OvS", New_Name), LogFC_OvsS, NA)) %>%
  mutate(LogFC_OvsA = if_else(grepl("OvA", New_Name), LogFC_OvsA, NA)) %>%
  filter(!(type == "Phototrophy" & grepl("OvS", New_Name))) %>%
  filter(!(type == "PH" & grepl("OvS", New_Name))) %>%
  filter(!(type == "Stress Response" & grepl("OvS", New_Name)))

# Order of Genes
stress_order <- c("2933779", "1605549", "5912680", "751180",
                  "3026034", "2305313", "1690207")
trans_order <- c("1026009", "778892",
                 "1396758", "2794147", "3616774", "41103",
                 "4392548", "2822257", "4005969", "2830222",
                 "3519852", "4065917", "4525068", "5822789",
                 "2266959", 	
                 "3965802", "5374355", "2916846", "912361", 
                 "1656590", "253618", "1388862", "2595805", "907432", "5368877")
het_order <- c("3976743", "1141532", "406252", "3898829",
               "3185097", "2569928", "191792")
photo_order <- c("634406", "3024000", "5609137", "2660400",
                 "5988908", "4228591", "3154977", "5662314")
order_heatmap <- as.data.frame(c(stress_order, trans_order, het_order, photo_order))

# Heatmap pivot
# READ ME: there is an amendment to this at the very end with StX_PAR data
# this information is from "O2, Fluorescence, and PAR Depth Manipluations and Plots"
# section below
Heatmap <- RPKM_combined_GeneGroup %>%
  pivot_longer(cols=c(row.names(metadata13)), names_to = 'Sample', values_to = 'RPKM') %>%
  mutate(logRPKM = log10(RPKM + 0.1)) %>%
  mutate(Group = case_when(
    Sample %in% surface ~ "Surface",
    Sample %in% oxycline ~ "Oxycline",
    Sample %in% anoxic ~ "Anoxic")) %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St_1_",
    Sample %in% station2 ~ "St_2_",
    Sample %in% station3 ~ "St_3_"
  )) %>%
  unite(Samp, Station, Sample, sep = "") %>%
  arrange(match(Col, order_heatmap$`c(stress_order, trans_order, het_order, photo_order)`))

# Function for bold in figure
bold_specific_strings <- function(text) {
  strings_to_bold <- c("All", "AvS", "OvS", "OvA")  # Replace with your specific strings
  
  for (str in strings_to_bold) {
    text <- gsub(str, paste0("<b>", str, "</b>"), text, fixed = TRUE)
  }
  
  return(text)
}

extract_depth <- function(samp) {
  as.numeric(sub(".*_(\\d+)m$", "\\1", samp))
}
St1_vector <- unique(Heatmap$Samp[grepl("St_1", Heatmap$Samp) & extract_depth(Heatmap$Samp) <= St1_PAR])
St2_vector <- unique(Heatmap$Samp[grepl("St_2", Heatmap$Samp) & extract_depth(Heatmap$Samp) <= St2_PAR])
St3_vector <- unique(Heatmap$Samp[grepl("St_3", Heatmap$Samp) & extract_depth(Heatmap$Samp) <= St3_PAR])
 
bold_sample_strings <- function(text) {
  strings_to_bold <- c(St1_vector, St2_vector, St3_vector)  # Replace with your specific strings
  
  for (str in strings_to_bold) {
    text <- gsub(str, paste0("<b>", str, "</b>"), text, fixed = TRUE)
  }
  
  return(text)
}

# Heatmap #####
heatmap <- ggplot(Heatmap, aes((Samp), fct_inorder(New_Name))) +
  geom_tile(width = 0.9, height = 0.9, aes(fill = logRPKM),
            color = "white", linewidth = 1.5, linetype = 1)+
  scale_fill_viridis(option = "viridis", direction = -1, limits = c(-1, 3))  +
  facet_grid(factor(type, levels = c("Stress Response", "Transporters", "PH", "Phototrophy")) ~ Group,
             scales = "free", space = "free") +
  bw + theme(axis.title.y = element_blank(),
             axis.text.x = element_markdown(angle = 45, hjust = 1, size = 18),
             axis.text.y = element_markdown(size = 18, vjust = 0.5, hjust=1, angle = 25),
             axis.title.x = element_text(size = 22),
             strip.text = element_text(size = 22),
             strip.background = element_blank(),
             legend.position = "top",
             legend.direction = "horizontal",
             legend.key.size = unit(2, 'cm'),
             legend.text = element_text(size = 19),
             legend.title = element_text(size = 19)) +
  scale_y_discrete(labels = function(x) bold_specific_strings(x)) +
  scale_x_discrete(labels = function(x) bold_sample_strings(x)) +
  labs(x = "Sample")
ggsave("heatmap5.png", heatmap, height = 24, width = 18.5, units = "in")
ggsave("heatmapFINAL.tiff", heatmap, height = 24, width = 18.5, units = "in", dpi = 300)

# Log FC Plot #####
logFCAS <- ggplot(RPKM_combined_GeneGroup, aes((LogFC_AvsS), fct_inorder(New_Name), fill = LogFC_AvsS)) +
  geom_bar(stat = "identity") + scale_fill_gradientn(colors = c("gold3", "gold2", "gold1",
                                                                "greenyellow", "forestgreen", "aquamarine3",
                                                                "dodgerblue" , "blue", "navyblue"),
                                                     limits = c(-15, 15),
                                                     breaks = c(-15, -10, -5 , 0, 5, 10, 15)) +
  facet_grid(type ~ . , scales = "free", space = "free",
             labeller = labeller(label_wrap_gen(width = 10), label_parsed)) +
  xlim(-14, 14) +
  geom_vline(xintercept=0, color="black", linewidth=0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 55, hjust = 1, color = "black", face = "bold"),
                     strip.background = element_blank(),
                     axis.text.y = element_text(size = 12, vjust = 0.5, hjust=1, angle = 35),
                     strip.text.x = element_text(size = 12),
                     strip.text.y = element_text(size = 12),
                     strip.text.y.right = element_blank(),
                     legend.position = "none",
                     plot.margin = margin(b = 30, unit = "pt")
  )  +  
  scale_y_discrete(labels = function(x) bold_specific_strings(x)) +
  theme(axis.text.y = element_markdown(), element_text(size = 12, vjust = 0.5, hjust=1, angle = 35)) + ggtitle("Anoxic vs Surface")
logFCOS <- ggplot(RPKM_combined_GeneGroup, aes((LogFC_OvsS), fct_inorder(New_Name), fill = LogFC_OvsS)) +
  geom_bar(stat = "identity") + scale_fill_gradientn(colors = c("gold3", "gold2", "gold1",
                                                                "greenyellow", "forestgreen", "aquamarine3",
                                                                "dodgerblue" , "blue", "navyblue"),
                                                     limits = c(-15, 15),
                                                     breaks = c(-15, -10, -5 , 0, 5, 10, 15)) +
  facet_grid(type ~ . , scales = "free", space = "free",
             labeller = labeller(label_wrap_gen(width = 10), label_parsed)) +
  xlim(-14, 14) +
  geom_vline(xintercept=0, color="black", linewidth=0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 55, hjust = 1, color = "black", face = "bold"),
                     strip.background = element_blank(),
                     axis.text.y = element_blank(),
                     strip.text.x = element_text(size = 12),
                     strip.text.y = element_blank(),
                     plot.margin = margin(b = 30, unit = "pt"),
                     legend.position = "none"
  )  + ggtitle("Oxycline vs Surface") 
logFCOA <- ggplot(RPKM_combined_GeneGroup, aes((LogFC_OvsA), fct_inorder(New_Name), fill = LogFC_OvsA)) +
  geom_bar(stat = "identity") + scale_fill_gradientn(colors = c("gold3", "gold2", "gold1",
                                                                "greenyellow", "forestgreen", "aquamarine3",
                                                                "dodgerblue" , "blue", "navyblue"),
                                                     limits = c(-15, 15),
                                                     breaks = c(-15, -10, -5 , 0, 5, 10, 15)) +
  facet_grid(type ~ . , scales = "free", space = "free",
             labeller = labeller(label_wrap_gen(width = 10), label_parsed)) +
  xlim(-14, 14) +
  geom_vline(xintercept=0, color="black", linewidth=0.5) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 55, hjust = 1, color = "black", face = "bold"),
                     strip.background = element_blank(),
                     strip.text.x = element_text(size = 12),
                     strip.text.y = element_text(size = 12, color = "black", face = "bold"),
                     axis.text.y.left = element_blank(),
                     plot.margin = margin(b = 30, unit = "pt")
  ) + 
  scale_y_discrete(labels = function(x) bold_specific_strings(x)) +
  theme(axis.text.y = element_markdown(), element_text(size = 12, vjust = 0.5, hjust=1, angle = 35)) + labs(fill = "log2FC") +  ggtitle("Oxycline vs Anoxic")

logFCplot <- logFCAS + logFCOS + logFCOA 
ggsave("logFCplot.png", logFCplot, height = 26.5, width = 15, units = "in")

# Log FC Table #####
RPKM_combined_GeneGroup2 <- RPKM_combined %>%
  mutate(type = case_when(
    ID %in% transporters ~ "Transporters",
    ID %in% stressresponse ~ "Stress Response",
    ID %in% phototrophy ~ "Phototrophy",
    ID %in% heterotrophy ~ "Potential Heterotrophy"
  )) %>%
  mutate(type = replace_na(type, "Stress Response")) %>%
  unite(New_Name, Name, Comparison, sep = " ") 
logFC_AS_Table <- RPKM_combined_GeneGroup2  %>%
  filter(str_detect(New_Name, "AvS")) %>%
  mutate(New_Name = str_remove_all(New_Name, "AvS|OvS|OvA| \\| ")) %>%
  select(Chr, LogFC_AvsS, New_Name, FDR_AvsS) %>%
  mutate(Chr_AvsS = Chr) %>%
  select(-Chr) %>%
  na.omit()
logFC_OS_Table <- RPKM_combined_GeneGroup2  %>%
  filter(str_detect(New_Name, "OvS")) %>%
  mutate(New_Name = str_remove_all(New_Name, "AvS|OvS|OvA| \\| ")) %>%
  select(Chr, LogFC_OvsS, New_Name, FDR_OvsS) %>%
  mutate(Chr_OvsS = Chr) %>%
  select(-Chr) %>%
  na.omit()
logFC_OA_Table <- RPKM_combined_GeneGroup2  %>%
  filter(str_detect(New_Name, "OvA")) %>%
  mutate(New_Name = str_remove_all(New_Name, "AvS|OvS|OvA| \\| ")) %>%
  select(Chr, LogFC_OvsA, New_Name, FDR_OvsA) %>%
  mutate(Chr_OvsA = Chr) %>%
  select(-Chr) %>%
  na.omit()
logFC_Table <- full_join(logFC_AS_Table, logFC_OS_Table) %>%
  full_join(logFC_OA_Table) %>%
  mutate(Gene = New_Name) %>%
  select(-New_Name) 
write.csv(logFC_Table, "logFC_contig_table.csv")

# Top Polarella RPKM Overall, Sup Table #####
RPKM <- rpkmData %>%
  inner_join(forRejoin) %>%
  inner_join(bastaClean) %>%
  inner_join(pantherAnnoSimp) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(surface, anoxic, oxycline))))) %>%
  arrange(desc(sumRPKM)) %>%
  group_by(ID) %>%
  top_n(1, (sumRPKM)) %>% 
  filter(sumRPKM != 0) %>%
  arrange(desc(sumRPKM)) %>%
  ungroup %>%
  head(10) %>%
  select(sumRPKM, Chr, SignatureDescription, Descriptive)
Name <- c("DNA Binding Protein HU", "Protease S8 Tripeptidyl Peptidase",
          "60S Ribosomal Protein L28", "Ribonucleotide Reductase Small Subunit",
          "Adenosylhomocysteinase", "Archaeal/bacterial/fungal rhodopsins", 
          "Cold Shock Domain Protein", "14-3-3 Protein", "Sodium/Chloride Dependent Transporter",
          "RuBiSCo") 
RPKM$Name <- Name
RPKM_SupTable <- RPKM %>% select(sumRPKM, Name, Chr)
write.csv(RPKM_SupTable, "TopRPKM_SupTable.csv")

# Polarella glacilias specific genes #####
# File for Polarella glacialis to 95% confidence at the species level
pg95 <- read.delim("polarella_95.txt", header = FALSE)
pg95_DEG <- PolarellaDEG %>%
  filter(Chr %in% pg95$V1)
# Genes DEG Functionality
knownFunction_P_Genus <- PolarellaDEG %>%
  left_join(pantherAnnoSimp) %>%
  filter(SignatureDescription != "unknown") %>%
  filter(SignatureDescription != "UNCHARACTERIZED") %>%
  filter(SignatureDescription != "-")
# Genus, Unique Function
uniqueFunction_P_Genus <- count(knownFunction_P_Genus, SignatureDescription)  
# Species DEG Functionality
knownFunction_P_Sp <- pg95_DEG %>%
  left_join(pantherAnnoSimp) %>%
  filter(SignatureDescription != "unknown") %>%
  filter(SignatureDescription != "UNCHARACTERIZED") %>%
  filter(SignatureDescription != "-")
# Species, Unique Function
uniqueFunction_P_Sp <- count(knownFunction_P_Sp, SignatureDescription)  

# Polarella all RPKM Depth Profile #####
# Polarella RPKM in each cast 
# Need to use edgeR again as this RPKM data will use the ribo0 data (prok + euk, to match DNA)
ribo0 <- read.table("ribo0_co_spades_featureCount_ALL.txt", sep = "\t", header = TRUE, row.names = "Geneid") %>%
  rename_with(
    ~ str_extract(., "Cast\\d+_\\d+m"),
    .cols = 6:37)
ribo0$Col <- c(1:nrow(ribo0))
forRejoin_ribo0 <- ribo0 %>%
  select(c("Chr", "Col")) 
ribo0EDGER <- ribo0 %>%
  select(-"Chr") %>%
  relocate(last_col(), .before = everything())
bastaRibo <- read.delim("transcripts_50hits_p1_i50_p51.basta", header = FALSE, sep = "\t") 
bastaRiboClean <- bastaRibo %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "genus")) %>%
  filter(genus == "Polarella")  %>%
  mutate(Chr = str_extract(Chr, "_g\\d+_"))
# Factor sample depths
depth <- factor(metadata13$Category_new)
# Construct edgeR object
edOb <- DGEList(ribo0EDGER[,-c(1,2,3,4,5)], group=depth, genes=ribo0EDGER[,c(1,5), drop = FALSE])
# Remove genes with low counts
keep <- rowSums(cpm(edOb) > 0.1) >= 3
table(keep) 
edOb2 <- edOb[keep, , keep.lib.sizes=FALSE]
# Normalization for composition bias
edObNorm <- calcNormFactors(edOb2, method = "TMM")
# Design matrix
design <- model.matrix(~0+depth)
colnames(design) <- levels(depth)
design
# Estimate Dispersion
edObDis <- estimateDisp(edObNorm, design, robust=TRUE)
# Full RPKM Table
rpkmData_ribo0 <- as.data.frame(rpkm(edObDis))
rpkmData_ribo0$Col <- edObDis$genes$Col
rownames(rpkmData_ribo0) <- sub("_.*", "", rownames(rpkmData_ribo0))
rpkmData_ribo0_profile <- rpkmData_ribo0 %>%
  left_join(forRejoin_ribo0) %>%
  mutate(Chr = str_extract(Chr, "_g\\d+_"))
# Save the RPKM Data
write.csv(rpkmData_ribo0_profile, file = "rpkmData_ribo0.csv",  row.names = FALSE)
# Read Back for Analysis
rpkmData_ribo0_profile <- read.csv("rpkmData_ribo0.csv")
RPKM_Depth_Polarella <- rpkmData_ribo0_profile %>%
  mutate(Chr = str_extract(Chr, "_g\\d+_")) %>%
  filter(Chr %in% bastaRiboClean$Chr) %>%
  select(surface, anoxic, oxycline) %>%
  summarise(across(everything(), sum, .names = "{.col}")) %>%
  pivot_longer(cols=c(row.names(metadata13)), names_to = 'Sample', values_to = 'RPKM_Polarella') %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St 1",
    Sample %in% station2 ~ "St 2",
    Sample %in% station3 ~ "St 3"
  )) %>%
  mutate(Depth = as.numeric(str_extract(Sample, "(?<=_)\\d+(?=m)"))) %>%
  arrange(Depth)
RPKM_Depth_All <- rpkmData_ribo0_profile %>%
  select(surface, anoxic, oxycline) %>%
  summarise(across(everything(), sum, .names = "{.col}")) %>%
  pivot_longer(cols=c(row.names(metadata13)), names_to = 'Sample', values_to = 'RPKM_All') %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St 1",
    Sample %in% station2 ~ "St 2",
    Sample %in% station3 ~ "St 3"
  )) %>%
  mutate(Depth = as.numeric(str_extract(Sample, "(?<=_)\\d+(?=m)"))) %>%
  arrange(Depth)
RPKM_Depth <- left_join(RPKM_Depth_Polarella, RPKM_Depth_All) %>%
  mutate(
    RPKM_Polarella = as.numeric(unlist(RPKM_Polarella)),
    RPKM_All = as.numeric(unlist(RPKM_All))) %>%
  mutate(RPKM = (RPKM_Polarella / RPKM_All) * 100) 
polar_kaiju <- read.delim(file = "polarella_kaiju.txt", header = FALSE)
polar_kaiju_simple <- polar_kaiju[,c(1:2)] %>%
  mutate(kaiju_percent = V2)
colnames(polar_kaiju_simple) <- c("Sample", "kaiju_fraction", "kaiju_percent")
polar_percent <- left_join(polar_kaiju_simple, RPKM_Depth, by = "Sample") %>%
  arrange(Depth)
rpkm_percent <- ggplot(polar_percent) +
  geom_path(aes(x = RPKM, y = Depth, color = as.character(Station))) + 
  geom_point(aes(x = RPKM, y = Depth, color = as.character(Station))) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, by=50)) + 
  scale_x_continuous(limits = c(0, 1),position = "top") +
  labs(x = expression(atop("% of Total RNA Reads", italic("Polarella"))),
       y = NULL,
       color = "Station") + annotate("text", x = 0, y = 0, label = "(E)") + bw +
  theme(text = element_text(size = 15))
kaiju_percent <- ggplot(polar_percent) +
  geom_path(aes(x = kaiju_percent, y = Depth, color = as.character(Station))) + 
  geom_point(aes(x = kaiju_percent, y = Depth, color = as.character(Station))) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
    scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, by=50)) + 
    scale_x_continuous(limits = c(0, 0.04),position = "top") +
    labs(x = expression(atop("% of Metagenomic Reads", italic("Polarella"))),
         y = "Depth (m)",
         color = "Station") + annotate("text", x = 0, y = 0, label = "(D)") + bw +
  theme(text = element_text(size = 15))
# Depth of only P. glacialis 
pg95ribo <- read.delim("95ribo.txt", header = FALSE) %>%
  mutate(Chr = str_extract(V1, "_g\\d+_")) 
RPKM_Depth_Pg <- rpkmData_ribo0_profile %>%
  mutate(Chr = str_extract(Chr, "_g\\d+_")) %>%
  filter(Chr %in% pg95ribo$Chr) %>%
  select(surface, anoxic, oxycline) %>%
  summarise(across(everything(), sum, .names = "{.col}")) %>%
  pivot_longer(cols=c(row.names(metadata13)), names_to = 'Sample', values_to = 'RPKM_Polarella') %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St 1",
    Sample %in% station2 ~ "St 2",
    Sample %in% station3 ~ "St 3"
  )) %>%
  mutate(Depth = as.numeric(str_extract(Sample, "(?<=_)\\d+(?=m)"))) %>%
  arrange(Depth)
RPKM_Depth_All_PG <- rpkmData_ribo0_profile %>%
  select(surface, anoxic, oxycline) %>%
  summarise(across(everything(), sum, .names = "{.col}")) %>%
  pivot_longer(cols=c(row.names(metadata13)), names_to = 'Sample', values_to = 'RPKM_All') %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St 1",
    Sample %in% station2 ~ "St 2",
    Sample %in% station3 ~ "St 3"
  )) %>%
  mutate(Depth = as.numeric(str_extract(Sample, "(?<=_)\\d+(?=m)"))) %>%
  arrange(Depth)
RPKM_Depth_PG <- left_join(RPKM_Depth_Pg, RPKM_Depth_All_PG) %>%
  mutate(
    RPKM_Polarella = as.numeric(unlist(RPKM_Polarella)),
    RPKM_All = as.numeric(unlist(RPKM_All))) %>%
  mutate(RPKM_Pg = (RPKM_Polarella / RPKM_All) * 100) 
rpkm_percent_pg <- ggplot(RPKM_Depth_PG) +
  geom_path(aes(x = RPKM_Pg, y = Depth, color = as.character(Station))) + 
  geom_point(aes(x = RPKM_Pg, y = Depth, color = as.character(Station))) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, by=50)) + 
  scale_x_continuous(limits = c(0, 0.4), position = "top") +
  labs(x = expression(atop("% of Total RNA Reads", italic("Polarella glacialis"))),
       y = NULL,
       color = "Station") + annotate("text", x = 0, y = 0, label = "(F)") + bw +
  theme(text = element_text(size = 15))

# O2, Fluorescence, and PAR Depth Manipluations and Plots #####
colnames_fordata <- c("Time_seconds", "Lat", "Long", "Pressure_db", "Depth_m", 
                      "Density_kgm3", "Salinity_psu", "Temp_C", "C_Star_m-1", "Fluorescence",
                      "PAR_quantacm2s", "Oxygen_umolL", "V4", "V5", "Flags")
colnames(cast18) <- colnames_fordata
colnames(cast59) <- colnames_fordata
colnames(cast93) <- colnames_fordata
cast18$Station <- "St 1"
cast59$Station <- "St 2"
cast93$Station <- "St 3"
cast18_1 <- cast18 %>%
  select(c("Depth_m", "Station", "Oxygen_umolL", "Fluorescence" ,
           "PAR_quantacm2s")) %>%
  filter(Depth_m < 300)
cast59_1 <- cast59 %>%
  select(c("Depth_m", "Station", "Oxygen_umolL", "Fluorescence",
           "PAR_quantacm2s")) %>%
  filter(Depth_m < 300)
cast93_1 <- cast93 %>%
  select(c("Depth_m", "Station", "Oxygen_umolL", "Fluorescence",
           "PAR_quantacm2s")) %>%
  filter(Depth_m < 300)
oxygenFluorPar <- rbind(cast59_1, cast93_1)
oxygenFluorPar2 <- rbind(oxygenFluorPar, cast18_1)

# Depth that anoxia is reached
St1_anoxia <- as.numeric(oxygenFluorPar2 %>%
                          filter(Station == "St 1") %>%
                          filter(Oxygen_umolL < 1) %>%
                          arrange(Depth_m) %>%
                          slice(1) %>%
                          pull(Depth_m))
St2_anoxia <- as.numeric(oxygenFluorPar2 %>%
                           filter(Station == "St 2") %>%
                           filter(Oxygen_umolL < 1) %>%
                           arrange(Depth_m) %>%
                           slice(1) %>%
                           pull(Depth_m))
St3_anoxia <- as.numeric(oxygenFluorPar2 %>%
                           filter(Station == "St 3") %>%
                           filter(Oxygen_umolL < 1) %>%
                           arrange(Depth_m) %>%
                           slice(1) %>%
                           pull(Depth_m))
# Depth of 0.1% PAR
St1_PAR <- as.numeric(oxygenFluorPar2 %>%
                           filter(Station == "St 1") %>%
                           filter(PAR_quantacm2s < (max(PAR_quantacm2s)/1000)) %>%
                           arrange(Depth_m) %>%
                           slice(1) %>%
                           pull(Depth_m))
St2_PAR <- as.numeric(oxygenFluorPar2 %>%
                           filter(Station == "St 2") %>%
                           filter(PAR_quantacm2s < (max(PAR_quantacm2s)/1000)) %>%
                           arrange(Depth_m) %>%
                           slice(1) %>%
                           pull(Depth_m))
St3_PAR <- as.numeric(oxygenFluorPar2 %>%
                           filter(Station == "St 3") %>%
                           filter(PAR_quantacm2s < (max(PAR_quantacm2s)/1000)) %>%
                           arrange(Depth_m) %>%
                           slice(1) %>%
                           pull(Depth_m))

# Plots
O2 <- ggplot(oxygenFluorPar2) +
  geom_path(aes(Oxygen_umolL, Depth_m, color = as.character(Station))) + 
  geom_point(aes(Oxygen_umolL, Depth_m, color = as.character(Station))) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(-2, 250), position = "top") +
  geom_hline(yintercept=St3_anoxia, linetype="dashed", 
             color = "blue", size=1) +
  geom_hline(yintercept=St2_anoxia, linetype="dashed", 
             color = "forestgreen", size=1) +
  geom_hline(yintercept=St1_anoxia, linetype="dashed", 
             color = "goldenrod1", size=1) +
  labs(title = NULL, x = expression("Oxygen (\U003BCmol/L)"), y = "Depth (m)") +
  annotate("text", x = 0, y = 0, label = "(A)") +
  bw +
  theme(text = element_text(size = 15))
PAR <- ggplot(oxygenFluorPar2) +
  geom_path(aes(PAR_quantacm2s, Depth_m, color = as.character(Station))) + 
  geom_point(aes(PAR_quantacm2s, Depth_m, color = as.character(Station))) +
  scale_color_manual(name = "Station", values = c("goldenrod1", "forestgreen", "blue")) +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(-2, 250), position = "top") +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=1) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=1) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=1) +
  labs(title = NULL, x = expression("PAR (W m"^-2*")"), y = NULL) +
  annotate("text", x = 0, y = 0, label = "(B)") +
  bw +
  theme(legend.position = "top",
        legend.title = element_text(size = 15),  # Increase title size
        legend.text = element_text(size = 15),
        text = element_text(size = 15))
Flur <- ggplot(oxygenFluorPar2) +
  geom_path(aes(Fluorescence, Depth_m, color = as.character(Station))) + 
  geom_point(aes(Fluorescence, Depth_m, color = as.character(Station))) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(0, 5), position = "top") +
  labs(title = NULL, x = expression("Fluorescence (mg m"^-3*")"), y = NULL, colour = "Station") +
  annotate("text", x = 0, y = 0, label = "(C)") +
  bw +
  theme(text = element_text(size = 15))
O2FlurPAR <- ((O2 | PAR | Flur) / 
              (kaiju_percent | rpkm_percent | rpkm_percent_pg))
ggsave("O2FlurPAR.png", O2FlurPAR, height = 9, width = 9, units = "in")
ggsave("O2FlurPAR.tiff", O2FlurPAR, height = 9, width = 9, units = "in", dpi = 300)

# Finding Percent Polarella from Genus Data in RPKM Table #####
RPKMbasta <- rpkmData %>%
  left_join(forRejoin)
bClean <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "genus"))
RPKMbasta_2 <- as.data.frame(RPKMbasta %>%
                               left_join(bClean)) 
RPKMbasta_2[RPKMbasta_2 == ""] <- NA
RPKMbasta_2[RPKMbasta_2 == "unknown"] <- NA
genus <- RPKMbasta_2 %>%
  count(genus) 
genus$percent <- genus$n / nrow(RPKMbasta_2) * 100
genus <- arrange(genus, percent)
genus2 <- genus[1743:1751, ]
genus3 <- c(genus = "all others", freq = (sum(genus$n) - sum(genus2$n)), percent = (100 - sum(genus2$percent)))
genus4 <- rbind(genus2, genus3)

# Depth Profiles with Trophic-Mode Genes #####
# Preparring Table
# Read in DeepLOC and DBCAN table with DE extracellular genes
deepLocDBCAN <- read.csv("DeepLocDBCANextra.csv")
RPKM_TrophicMode_Depth <- PolarellaDEG %>%
  filter(ID %in% phototrophy | Chr %in% deepLocDBCAN$Chr) %>%
  group_by(ID) %>%
  left_join(selectedGenes) %>%
  summarise(across(all_of(row.names(metadata13)), sum),
            across(!all_of(row.names(metadata13)), first)) %>%
  ungroup() %>%
  pivot_longer(cols=c(row.names(metadata13)), names_to = 'Sample', values_to = 'RPKM_All') %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St 1",
    Sample %in% station2 ~ "St 2",
    Sample %in% station3 ~ "St 3"
  )) %>%
  mutate(Depth = as.numeric(str_extract(Sample, "(?<=_)\\d+(?=m)"))) %>%
  arrange(Depth) 
subsetRPKM_trophic <- RPKM_TrophicMode_Depth %>%
  select(Name) %>%
  distinct() %>% arrange(Name)

# Photosynthetic Genes
CarbonicAn <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[4,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,80)) + 
  annotate("text", x = 75, y = 0, label = "(A)", size = 7) +
  labs(title = "Carbonic \nAndyrase", x = "RPKM", y = "Depth (m)") + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
ChlAB <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[6,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,350)) + 
  annotate("text", x = 330, y = 0, label = "(B)", size = 7) +
  labs(title = "Chl A/B Binding \nProtein", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
RubISCo <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[9,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,400)) + 
  annotate("text", x = 375, y = 0, label = "(C)", size = 7) +
  labs(title = "RuBiSCo", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
Rhodopsin <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[8,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,1750)) + 
  annotate("text", x = 1680, y = 0, label = "(D)", size = 7) +
  labs(title = "Rhodopsin", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
ATPsyn <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[2,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,40)) + 
  annotate("text", x = 36, y = 0, label = "(E)", size = 7) +
  labs(title = "ATP Synthase \nGamma Subunit", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))

# Heterotrophic Genes
ABCTrans <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[1,])) %>%
  ggplot() +
    geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
    geom_point(aes(RPKM_All, Depth, color = Station)) +
    scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
    scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
    geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
    geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
    geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
    geom_hline(yintercept=St3_PAR, linetype="dotted", 
            color = "blue", size=0.5) +
    geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
    geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
    scale_x_continuous(lim=c(0,5)) + 
    annotate("text", x = 4.5, y = 0, label = "(F)", size = 7) +
    labs(title = "ABC Transporter", x = "RPKM", y = "Depth (m)") + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
GH3 <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[3,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  annotate("text", x = 8, y = 0, label = "(G)", size = 7) +
  labs(title = "Beta-Glucosidase \n(GH3)", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
GH18 <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[5,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,25)) + 
  annotate("text", x = 22.5, y = 0, label = "(H)", size = 7) +
  labs(title = "Chitinase \n(GH18)", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
GH1 <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[7,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), linewidth = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) + 
  geom_hline(yintercept=St3_anoxia, linetype="twodash", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_anoxia, linetype="twodash", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_anoxia, linetype="twodash", 
             color = "goldenrod1", size=0.5) +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=0.5) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=0.5) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=0.5) +
  scale_x_continuous(lim=c(0,10)) + 
  annotate("text", x = 9, y = 0, label = "(I)", size = 7) +
  labs(title = "Glycosyl \nHydrolase (GH1)", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
    axis.title = element_text(size = 23),  
    plot.title = element_text(size = 25, face = "bold"))
# For legend
# Create a dummy dataset
dummy_data <- data.frame(
  x = c(1, 2, 3),
  y = c(1, 2, 3),
  group = factor(c("St. 1", "St. 2", "St. 3"))  # Example groups for the legend
)
# Create a plot with only the legend
p <- ggplot(dummy_data, aes(x = x, y = y, color = group)) +
  geom_point(size = 5, alpha = 0) +  
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue")) +  
  labs(color = "Station") +      
  theme_void() +                       
  theme(legend.position = c(0.5, 0.5), 
        legend.justification = c("center", "center"),
        legend.key.size = unit(1.5, 'cm'),  
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 30)) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) 

RPKM_TrophicMode_Depth_Plot <- (CarbonicAn | ChlAB  | RubISCo | Rhodopsin | ATPsyn) /
  (ABCTrans | GH3 | GH18 | GH1 | p)
ggsave("TMDP.png", RPKM_TrophicMode_Depth_Plot, height = 12, width = 18, units = "in")
ggsave("TMDP.tiff", RPKM_TrophicMode_Depth_Plot, height = 12, width = 18, units = "in", dpi = 300)

# Nutrient Depth Profiles #####

NH4 <- Nh4 %>%
  pivot_longer("NH4_nM", names_to = "Const", values_to = "Conc") %>%
  drop_na()
NOX <- No3No2 %>%
  pivot_longer(c("NO3_uM", "NO2_uM"), names_to = "Const", values_to = "Conc") %>%
  drop_na()
NOX$Station <- gsub("PS","", NOX$Station)
NH4$Station <- as.integer(NH4$Station)
NOX$Station <- as.integer(NOX$Station)
nitrogenDepth <- rbind(NH4, NOX)
ST1nitrogenDepth <- nitrogenDepth %>%
  filter(Station == 1) %>%
  filter(Depth_m < 151) 
ST2nitrogenDepth <- nitrogenDepth %>%
  filter(Station == 2) %>%
  filter(Depth_m < 901) %>%
  filter(Depth_m != 300) %>%
  filter(Depth_m != 260)
ST3nitrogenDepth <- nitrogenDepth %>%
  filter(Station == 3) %>%
  filter(Depth_m < 101) 
nitroAgain <- rbind(ST1nitrogenDepth, ST2nitrogenDepth)
nitrogenDepth2 <- rbind(nitroAgain, ST3nitrogenDepth)
nitrogenDepth3 <- nitrogenDepth2 %>%
  pivot_wider(names_from = Const, values_from = Conc) %>%
  arrange(Depth_m)
NH4 <- nitrogenDepth3 %>%
  drop_na(NH4_nM) %>%
  ggplot() +
  geom_point(aes(NH4_nM, Depth_m, color = as.character(Station))) +
  geom_path(aes(NH4_nM, Depth_m, color = as.character(Station))) + 
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue")) +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(0, 700), position = "top") +
  labs(title = NULL, x = expression("Ammonium (nM)"), y = "Depth (m)",colour = "Station") +
  annotate("text", x = 675, y = 0, label = "(A)") + guides(linetype = guide_legend(order = 2),
                                                           color = "none") +
  bw + theme(text = element_text(size = 20))
NO2 <- nitrogenDepth3 %>%
  drop_na(NO2_uM) %>%
  ggplot() +
  geom_point(aes(NO2_uM, Depth_m, color = as.character(Station))) +
  geom_path(aes(NO2_uM, Depth_m, color = as.character(Station))) + 
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue")) +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(0, 7), position = "top") +
  labs(title = NULL, x = expression("Nitrite (\U003BCM)"), y = NULL,colour = "Station") +
  annotate("text", x = 6.75, y = 0, label = "(B)") + guides(linetype = guide_legend(order = 2),
                                                           color = "none") +
  bw + theme(text = element_text(size = 20))
NO3 <- nitrogenDepth3 %>%
  drop_na(NO3_uM) %>%
  ggplot() +
  geom_point(aes(NO3_uM, Depth_m, color = as.character(Station))) +
  geom_path(aes(NO3_uM, Depth_m, color = as.character(Station))) + 
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue")) +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(0, 90), position = "top") +
  labs(title = NULL, x = expression("Nitrate (\U003BCM)"), y = NULL,colour = "Station") +
  annotate("text", x = 87.5, y = 0, label = "(C)") +
  bw + theme(text = element_text(size = 20))
nutDepthProfiles <- (NH4 | NO2 | NO3)
ggsave("nutDepthProfiles.png", nutDepthProfiles, height = 4.5, width = 9, units = "in")
ggsave("nutDepthProfiles.tiff", nutDepthProfiles, height = 4.5, width = 9, units = "in", dpi = 300)
