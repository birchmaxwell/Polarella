# Pub Code 2025.03.31
# For Analysis ensure packages are installed and files for analysis are downloaded
# Read in lines 5-194 for all files
# DEG analysis itself is 195-317
# Post DEG - analysis files for plots and data manipulation is 318- end of doc 

# Packages #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages(c("tidyverse", "ggtext", "ggpubr",
                   "ggbreak", "ggrepel", "patchwork",
                   "openxlsx", "svglite", "car"))
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
library(ggbreak)
library(svglite)
library(car)

# plot theme
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

# Read count table preprocessed PolyA
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
# Taxonomy BASTA PolyA
# We are focused on Polarella -> create a clean Polarella file
basta <- read.delim("final.contigs-gmst_50hits_50identity.basta", header = FALSE, sep = "\t") 
# Read in additional taxonomy from MarFERReT
MarFERReT_polyA_genus <- read.csv("polyA_Mar_Genus.csv") %>%
  filter(genus == "Polarella")
colnames(MarFERReT_polyA_genus) <- c("Chr", "genus")
basta2 <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "genus")) 
bastaClean <- basta2 %>%
  filter(if_any(everything(), ~. == "Polarella"))
bastaClean <- rbind(bastaClean, MarFERReT_polyA_genus)
# Polarella at 95%
pg95 <- read.delim("polarella_95.txt", header = FALSE)
# For bar plots additional MarFERReT Data
MarFERReT_polyA_genus2 <- read.csv("polyA_Mar_Genus.csv")
MarFERReT_polyA_family2 <- read.csv("polyA_Mar_Family.csv")
MarFERReT_polyA_order2 <- read.csv("polyA_Mar_Order.csv")
MarFERReT_polyA_class2 <- read.csv("polyA_Mar_Class.csv")
MarFERReT_polyA_phylum2 <- read.csv("polyA_Mar_Phylum.csv")
colnames(MarFERReT_polyA_genus2) <- c("Chr", "genus")
colnames(MarFERReT_polyA_family2) <- c("Chr", "family")
colnames(MarFERReT_polyA_order2) <- c("Chr", "order")
colnames(MarFERReT_polyA_class2) <- c("Chr", "class")
colnames(MarFERReT_polyA_phylum2) <- c("Chr", "phylum")

# For ribo0 plots
ribo0_Mar <- read.csv("ribo0_Mar_Genus.csv")
colnames(ribo0_Mar) <- c("Chr", "genus")
ribo0_Mar <- ribo0_Mar %>%
  mutate(Chr = as.character(Chr))
ribo0 <- read.table("ribo0_co_spades_featureCount_ALL.txt", sep = "\t", header = TRUE, row.names = "Geneid") %>%
  rename_with(
    ~ str_extract(., "Cast\\d+_\\d+m"),
    .cols = 6:37)
bastaRibo <- read.delim("transcripts_50hits_p1_i50_p51.basta", header = FALSE, sep = "\t")
bastaRiboClean <- bastaRibo %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "genus")) %>%
  filter(genus == "Polarella")  %>%
  mutate(Chr = str_extract(Chr, "_g\\d+_")) %>%
  left_join(ribo0_Mar, by = "Chr", suffix = c("_bClean", "_MarFERReT")) %>%
  mutate(
    genus = ifelse(is.na(genus_bClean), genus_MarFERReT, genus_bClean)
  ) %>%
  select(Chr, genus)
pg95ribo <- read.delim("95ribo.txt", header = FALSE) %>%
  mutate(Chr = str_extract(V1, "_g\\d+_")) 
# For DNA plots
polar_kaiju <- read.delim(file = "polarella_kaiju.txt", header = FALSE)
polar_kaiju_simple <- polar_kaiju[,c(1:2)] %>%
  mutate(kaiju_percent = V2)
colnames(polar_kaiju_simple) <- c("Sample", "kaiju_fraction", "kaiju_percentC")
kaijuUnclass <- read.csv("kaijuUnPolar.csv") %>% select(-X)
# Panther/InterproScan (Annotations)
# From prior prelim analysis, we know the functional annotations we want to focus on
# prelim analysis contained extensive annotations
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
                     "PTHR12400", "PTHR14218", "PTHR43607", "PTHR43389")
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
                              "Auxiliary Activities 1", "ATP Synthase Gamma Subunit", "Polyphosphate Kinase*",
                              "Protease S8 tripeptidyl peptidase 1", "V-Type proton ATPase, subunit A", "V-Type proton ATPase, subunit B"))
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
# Read in DeepLOC and DBCAN table with DE extracellular genes
deepLocDBCAN <- read.csv("DeepLocDBCANextra.csv")

##### PolyA RPKM Generate and DE Analysis #####
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
# write.csv(rpkmData, file = "rpkmData.csv",  row.names = FALSE)
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
# write.csv(EdgeR_Res, "EdgeR_Res.csv")
##### Ribo RPKM #####
ribo0$Col <- c(1:nrow(ribo0))
forRejoin_ribo0 <- ribo0 %>%
  select(c("Chr", "Col")) 
ribo0EDGER <- ribo0 %>%
  select(-"Chr") %>%
  relocate(last_col(), .before = everything())
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
# write.csv(rpkmData_ribo0_profile, file = "rpkmData_ribo0.csv",  row.names = FALSE)

##### READING in RPKM Files and DEG RES FILES #####
rpkmData_ribo0_profile <- read.csv("rpkmData_ribo0.csv")
rpkmData <- read.csv("rpkmData.csv")
EdgeR_Res <- read.csv("EdgeR_Res.csv")

##### O2, Fluor, PAR ####
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
  geom_point(aes(Oxygen_umolL, Depth_m, color = as.character(Station), shape = as.character(Station)), size = 2) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 250, y = 0, label = "(A)", 
           hjust = 1, vjust = 1, size = 4) +
  bw +
  theme(text = element_text(size = 15))
PAR <- ggplot(oxygenFluorPar2 %>% filter(PAR_quantacm2s < 300)) +
  geom_path(aes(PAR_quantacm2s, Depth_m, color = as.character(Station))) + 
  geom_point(aes(PAR_quantacm2s, Depth_m, color = as.character(Station), shape = as.character(Station)), size = 2) +
  scale_shape_manual(name = "Station", values = c(15, 16,17)) +
  scale_color_manual(name = "Station", values = c("goldenrod1", "forestgreen", "blue")) +
  scale_y_reverse(lim=c(325,0), breaks=seq(0, 325, 50)) +
  scale_x_continuous(limits = c(-2, 325), position = "top") +
  geom_hline(yintercept=St3_PAR, linetype="dotted", 
             color = "blue", size=1) +
  geom_hline(yintercept=St2_PAR, linetype="dotted", 
             color = "forestgreen", size=1) +
  geom_hline(yintercept=St1_PAR, linetype="dotted", 
             color = "goldenrod1", size=1) +
  labs(title = NULL, x = expression("PAR (W m"^-2*")"), y = NULL) +
  annotate("text", x = 325, y = 0, label = "(B)", 
           hjust = 1, vjust = 1, size = 4) +
  bw +
  theme(legend.position = "top",
        legend.title = element_text(size = 15),  # Increase title size
        legend.text = element_text(size = 15),
        text = element_text(size = 15))
Flur <- ggplot(oxygenFluorPar2) +
  geom_path(aes(Fluorescence, Depth_m, color = as.character(Station))) + 
  geom_point(aes(Fluorescence, Depth_m, color = as.character(Station), shape = as.character(Station)), size = 2) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, 50)) +
  scale_x_continuous(limits = c(0, 5), position = "top") +
  labs(title = NULL, x = expression("Fluorescence (mg m"^-3*")"), y = NULL, colour = "Station") +
  annotate("text", x = 5, y = 0, label = "(C)", 
           hjust = 1, vjust = 1, size = 4) +
  bw +
  theme(text = element_text(size = 15))

##### Finding amount of DE Polarella genes #####
PolarellaDEG <- EdgeR_Res %>%
  inner_join(bastaClean, by = "Chr") %>%
  select(-X) %>%
  inner_join(rpkmData, by = "Col") %>%
  left_join(pantherAnnoSimp)
write.csv(PolarellaDEG, "PolarellaDEG_SupplementalFile.csv")
PolarellaDEGAS <- PolarellaDEG %>%
  filter(FDR_AvsS < 0.05) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(surface, anoxic)))))
nrow(PolarellaDEGAS)
nrow(PolarellaDEGAS %>%
       filter(LogFC_AvsS > 0))
nrow(PolarellaDEGAS %>%
       filter(LogFC_AvsS < 0))
PolarellaDEGOS <- PolarellaDEG %>%
  filter(FDR_OvsS < 0.05) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(surface, oxycline)))))
nrow(PolarellaDEGOS)
nrow(PolarellaDEGOS %>%
       filter(LogFC_OvsS > 0))
nrow(PolarellaDEGOS %>%
       filter(LogFC_OvsS < 0))
PolarellaDEGOA <- PolarellaDEG %>%
  filter(FDR_OvsA < 0.05) %>%
  mutate(sumRPKM = rowSums(select(., all_of(c(oxycline, anoxic)))))
nrow(PolarellaDEGOA)
nrow(PolarellaDEGOA %>%
       filter(LogFC_OvsA > 0))
nrow(PolarellaDEGOA %>%
       filter(LogFC_OvsA < 0))

##### Heatmap, LogFC Plot, and LogFC Table Maniplutations #####
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
                  "PTHR34939", "PTHR11002","PTHR10587",  "PTHR43400", "PTHR11709",
                  "PTHR14218", "PTHR43607", "PTHR43389")
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
het_order <- c("471220", "3976743", "1141532", "406252", "3898829",
               "3185097", "2569928", "191792", "4822737", "65353")
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
ggsave("heatmapFINAL_Review.tiff", heatmap, height = 24, width = 18.5, units = "in", dpi = 300)
ggsave("heatmapFINAL_Review.png", heatmap, height = 24, width = 18.5, units = "in", dpi = 300)
ggsave("heatmapFINAL_Review.pdf", heatmap, height = 24, width = 18.5, units = "in", dpi = 300)

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

##### Top Polarella RPKM Overall, Sup Table #####
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

##### Polarella glacilias specific genes and unique function genes #####
# File for Polarella glacialis to 95% confidence at the species level
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

##### Polarella all Depth Profile, Ribo0, DNA (kaiju) #####
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
polar_kaiju_simple2 <- polar_kaiju_simple %>%
  left_join(kaijuUnclass) %>%
  mutate(kaiju_percent = kaiju_percentC + percentWhole) %>%
  select(c(Sample, kaiju_percent))
polar_percent <- left_join(polar_kaiju_simple2, RPKM_Depth, by = "Sample") %>%
  arrange(Depth) 
RPKM_Depth_PolarellaG95 <- rpkmData_ribo0_profile %>%
  mutate(Chr = str_extract(Chr, "_g\\d+_")) %>%
  filter(Chr %in% bastaRiboClean$Chr) %>%
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
RPKM_Depth_PG <- left_join(RPKM_Depth_PolarellaG95, RPKM_Depth_All) %>%
  mutate(RPKM_PolarellaG = as.numeric(unlist(RPKM_Polarella)),
         RPKM_All = as.numeric(unlist(RPKM_All))) %>%
  mutate(RPKM_Pg = (RPKM_PolarellaG / RPKM_All) * 100) 
rpkm_percent_pg <- ggplot(RPKM_Depth_PG) +
  geom_path(aes(x = RPKM_Pg, y = Depth, color = as.character(Station))) + 
  geom_point(aes(x = RPKM_Pg, y = Depth, color = as.character(Station), shape = as.character(Station)), size = 2) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, by=50)) + 
  scale_x_continuous(limits = c(0, 0.04), position = "top") +
  labs(x = expression(atop(" ", italic("Polarella glacialis"))),
       y = NULL,
       color = "Station") + annotate("text", x = 0.04, y = 0, label = "(F)", 
                                     hjust = 1, vjust = 1, size = 4) + bw +
  theme(text = element_text(size = 15))
rpkm_percent <- ggplot(polar_percent) +
  geom_path(aes(x = RPKM, y = Depth, color = as.character(Station))) + 
  geom_point(aes(x = RPKM, y = Depth, color = as.character(Station), shape = as.character(Station)), size = 2) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, by=50)) + 
  scale_x_continuous(limits = c(0, 1),position = "top") +
  labs(x = expression(atop("                                                 % of rRNA-depleted RNA Reads", italic("Polarella"))),
       y = NULL,
       color = "Station") + annotate("text", x = 1, y = 0, label = "(E)", 
                                     hjust = 1, vjust = 1, size = 4) + bw +
  theme(text = element_text(size = 15))
kaiju_percent <- ggplot(polar_percent) +
  geom_path(aes(x = kaiju_percent, y = Depth, color = as.character(Station))) + 
  geom_point(aes(x = kaiju_percent, y = Depth, color = as.character(Station), shape = as.character(Station)), size = 2) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none") +
  scale_y_reverse(lim=c(300,0), breaks=seq(0, 300, by=50)) + 
  scale_x_continuous(limits = c(0, 0.04),position = "top") +
  labs(x = expression(atop("% of Metagenomic Reads", italic("Polarella"))),
       y = "Depth (m)",
       color = "Station") + annotate("text", x = 0.04, y = 0, label = "(D)", 
                                     hjust = 1, vjust = 1, size = 4) + bw +
  theme(text = element_text(size = 15))
percentOfRibo0 <- data.frame(Sample = RPKM_Depth_PG$Sample,
                             PG95 = RPKM_Depth_PG$RPKM_PolarellaG,
                             PG = RPKM_Depth_PG$RPKM_Polarella) %>%
  mutate(percent = PG95/PG*100)
column_mean <- mean(percentOfRibo0$percent, na.rm = TRUE)
column_sd <- sd(percentOfRibo0$percent, na.rm = TRUE)

# Plots with O2, PAR, Flur, DNA, mRNA #####
O2FlurPAR <- ((O2 | PAR | Flur) / 
                (kaiju_percent | rpkm_percent | rpkm_percent_pg))
ggsave("O2FlurPAR.tiff", O2FlurPAR, height = 9, width = 9, units = "in", dpi = 300)
ggsave("O2FlurPAR_FINAL_review.png", O2FlurPAR, height = 9, width = 9, units = "in", dpi = 300)
ggsave("O2FlurPAR_FINAL_review.pdf", O2FlurPAR, height = 9, width = 9, units = "in", dpi = 300)

##### Top Taxonomy at different taxa levels #####
# Genus #####
RPKMbasta <- rpkmData %>%
  left_join(forRejoin)
bClean <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "genus"))
bClean2 <- bClean %>%
  left_join(MarFERReT_polyA_genus2, by = "Chr", suffix = c("_bClean", "_MarFERReT")) %>%
  mutate(
    genus = ifelse(is.na(genus_bClean), genus_MarFERReT, genus_bClean)
  ) %>%
  select(Chr, genus)
RPKMbasta_2 <- as.data.frame(RPKMbasta %>%
                               left_join(bClean2)) %>%
  select(-c(Col, Chr))
RPKMbasta_2[RPKMbasta_2 == ""] <- NA
RPKMbasta_2[RPKMbasta_2 == "unknown"] <- NA
RPKMbasta_2[is.na(RPKMbasta_2)] <- "Unclassified"
genus <- RPKMbasta_2 %>%
  count(genus) 
genus$percent <- genus$n / nrow(RPKMbasta_2) * 100
genus <- arrange(genus, percent)
genus2 <- genus[1891:1899, ]
genus3 <- c(genus = "all others (n = 1890)", freq = (sum(genus$n) - sum(genus2$n)), percent = (100 - sum(genus2$percent)))
genus4 <- rbind(genus2, genus3) 
write.csv(genus4, file = "genusAll.csv",  row.names = FALSE)
valid_genus <- c("Karlodinium", "Emiliania", "Perkinsus", "Phytophthora", 
                 "Aureococcus", "Pelagomonas", "Polarella", "Symbiodinium", "Unclassified")
RPKMbasta_2_grouped <- RPKMbasta_2 %>%
  mutate(Genus_group = ifelse(genus %in% valid_genus, genus, "All Other (n = 1890)   ")) %>%
  group_by(Genus_group) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup() 
RPKMbasta_2_percentage <- RPKMbasta_2_grouped %>%
  select(Genus_group, starts_with("Cast")) %>%  # Select only genus_group and Cast columns
  mutate(across(starts_with("Cast"), 
                ~ . / sum(.) * 100,  # Calculate percentage for each column
                .names = "percent_{.col}")) %>%  # Naming the new columns as 'percent_Cast...'
  select(-starts_with("Cast"))
colnames(RPKMbasta_2_percentage) <- sub("_", "", colnames(RPKMbasta_2_percentage), fixed = TRUE)
colnames(RPKMbasta_2_percentage) <- gsub("percent", "", colnames(RPKMbasta_2_percentage))
colnames(RPKMbasta_2_percentage) <- gsub("group", "", colnames(RPKMbasta_2_percentage))
custom_order <- c("Symbiodinium", "Polarella", "Pelagomonas", "Aureococcus", "Phytophthora",
                  "Perkinsus", "Emiliania", "Karlodinium", "All Other (n = 1890)   ", "Unclassified")
RPKM_by_station_normalized <- RPKMbasta_2_percentage %>%
  pivot_longer(cols = starts_with("Cast"),  
               names_to = "Sample",  
               values_to = "Relative RPKM within Cast")  %>%
  mutate(Category_new = case_when(
    Sample %in% surface ~ "Surface",      # If Cast is in the surface vector
    Sample %in% oxycline ~ "Oxycline",    # If Cast is in the oxycline vector
    Sample %in% anoxic ~ "Anoxic"))  %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St_1_",
    Sample %in% station2 ~ "St_2_",
    Sample %in% station3 ~ "St_3_"
  )) %>%
  unite(Samp, Station, Sample, sep = "")
RPKM_by_station_normalized$Genus <- factor(
  RPKM_by_station_normalized$Genus,
  levels = custom_order)
# For genus only counts of Polarella per cast
RPKMbasta_2_number <- RPKMbasta_2 %>%
  filter(genus == "Polarella")
convert_cast_columns <- function(df) {
  cast_cols <- grep("^Cast", names(df), value = TRUE)
  result <- lapply(df[cast_cols], as.vector)
  return(result)
}
castNumber <- convert_cast_columns(RPKMbasta_2_number)
count_non_zero <- function(list_of_vectors) {
  lapply(list_of_vectors, function(x) sum(x != 0, na.rm = TRUE))
}
non_zero_counts <- count_non_zero(castNumber)
non_zero_df <- data.frame(
  cast = names(non_zero_counts),
  count = unlist(non_zero_counts))
mean_Polarella <- RPKMbasta_2_percentage %>% 
  filter(Genus == "Polarella") %>% select(-Genus) %>%
  rowwise() %>%
  mutate(RowMean = mean(c_across(everything()), na.rm = TRUE)) %>%
  mutate(RowSD = sd(c_across(everything()), na.rm = TRUE)) %>%
  pivot_longer(cols = everything()) %>%
  mutate(PercentEuk = value,
         cast = name) %>%
  select(-c(name, value))

# Family #####
bCleanFAMILY <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "family"))
bClean2FAMILY <- bCleanFAMILY %>%
  left_join(MarFERReT_polyA_family2, by = "Chr", suffix = c("_bClean", "_MarFERReT")) %>%
  mutate(
    family = ifelse(is.na(family_bClean), family_MarFERReT, family_bClean)
  ) %>%
  select(Chr, family)
RPKMbasta_2FAMILY <- as.data.frame(RPKMbasta %>%
                                     left_join(bClean2FAMILY)) %>%
  select(-c(Col, Chr))
RPKMbasta_2FAMILY[RPKMbasta_2FAMILY == ""] <- NA
RPKMbasta_2FAMILY[RPKMbasta_2FAMILY == "unknown"] <- NA
RPKMbasta_2FAMILY[is.na(RPKMbasta_2FAMILY)] <- "Unclassified"
family <- RPKMbasta_2FAMILY %>%
  count(family) 
family$percent <- family$n / nrow(RPKMbasta_2FAMILY) * 100
family <- arrange(family, percent)
family2 <- family[1185:1193, ]
family3 <- c(family = "all others (n = 1184)", freq = (sum(family$n) - sum(family2$n)), percent = (100 - sum(family2$percent)))
family4 <- rbind(family2, family3)
valid_family <- c("Bathycoccaceae", "Trypanosomatidae", "Perkinsidae", "Kareniaceae", 
                  "Peronosporaceae", "Pelagomonadaceae", "Suessiaceae", "Symbiodiniaceae", "Unclassified")
RPKMbasta_2FAMILY_grouped <- RPKMbasta_2FAMILY %>%
  mutate(family_group = ifelse(family %in% valid_family, family, "All Other (n = 1184)      ")) %>%
  group_by(family_group) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup() 
RPKMbasta_2FAMILY_percentage <- RPKMbasta_2FAMILY_grouped %>%
  select(family_group, starts_with("Cast")) %>%  # Select only family_group and Cast columns
  mutate(across(starts_with("Cast"), 
                ~ . / sum(.) * 100,  # Calculate percentage for each column
                .names = "percent_{.col}")) %>%  # Naming the new columns as 'percent_Cast...'
  select(-starts_with("Cast"))
colnames(RPKMbasta_2FAMILY_percentage) <- sub("_", "", colnames(RPKMbasta_2FAMILY_percentage), fixed = TRUE)
colnames(RPKMbasta_2FAMILY_percentage) <- gsub("percent", "", colnames(RPKMbasta_2FAMILY_percentage))
colnames(RPKMbasta_2FAMILY_percentage) <- gsub("group", "", colnames(RPKMbasta_2FAMILY_percentage))
custom_orderFAMILY <- c("Symbiodiniaceae", "Suessiaceae", "Pelagomonadaceae", "Peronosporaceae", "Kareniaceae",
                        "Perkinsidae", "Trypanosomatidae", "Bathycoccaceae", "All Other (n = 1184)      ", "Unclassified")
RPKM_by_station_normalizedFAMILY <- RPKMbasta_2FAMILY_percentage %>%
  pivot_longer(cols = starts_with("Cast"),  
               names_to = "Sample",  
               values_to = "Relative RPKM within Cast")  %>%
  mutate(Category_new = case_when(
    Sample %in% surface ~ "Surface",      # If Cast is in the surface vector
    Sample %in% oxycline ~ "Oxycline",    # If Cast is in the oxycline vector
    Sample %in% anoxic ~ "Anoxic"))  %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St_1_",
    Sample %in% station2 ~ "St_2_",
    Sample %in% station3 ~ "St_3_"
  )) %>%
  unite(Samp, Station, Sample, sep = "")
RPKM_by_station_normalizedFAMILY$family <- factor(
  RPKM_by_station_normalizedFAMILY$family,
  levels = custom_orderFAMILY)

# Order ####
bCleanORDER <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "order"))
bClean2ORDER <- bCleanORDER %>%
  left_join(MarFERReT_polyA_order2, by = "Chr", suffix = c("_bClean", "_MarFERReT")) %>%
  mutate(
    Order = ifelse(is.na(order_bClean), order_MarFERReT, order_bClean)
  ) %>%
  select(Chr, Order)
RPKMbasta_2ORDER <- as.data.frame(RPKMbasta %>%
                                    left_join(bClean2ORDER)) %>%
  select(-c(Col, Chr))
RPKMbasta_2ORDER[RPKMbasta_2ORDER == ""] <- NA
RPKMbasta_2ORDER[RPKMbasta_2ORDER == "unknown"] <- NA
RPKMbasta_2ORDER[is.na(RPKMbasta_2ORDER)] <- "Unclassified"
ORDER <- RPKMbasta_2ORDER %>%
  count(Order) 
ORDER$percent <- ORDER$n / nrow(RPKMbasta_2ORDER) * 100
ORDER <- arrange(ORDER, percent)
ORDER2 <- ORDER[605:613, ]
ORDER3 <- c(Order = "all others (n = 605)", freq = (sum(ORDER$n) - sum(ORDER2$n)), percent = (100 - sum(ORDER2$percent)))
ORDER4 <- rbind(ORDER2, ORDER3)

valid_ORDER <- c("Suessiales", "Pelagomonadales", "Peronosporales", "Gymnodiniales", 
                 "Mamiellales", "Perkinsida", "Isochrysidales", "Prymnesiales", "Unclassified")
RPKMbasta_2ORDER_grouped <- RPKMbasta_2ORDER %>%
  mutate(ORDER_group = ifelse(Order %in% valid_ORDER, Order, "All Other (n = 605)      ")) %>%
  group_by(ORDER_group) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup() 
RPKMbasta_2ORDER_percentage <- RPKMbasta_2ORDER_grouped %>%
  select(ORDER_group, starts_with("Cast")) %>%  # Select only ORDER_group and Cast columns
  mutate(across(starts_with("Cast"), 
                ~ . / sum(.) * 100,  # Calculate percentage for each column
                .names = "percent_{.col}")) %>%  # Naming the new columns as 'percent_Cast...'
  select(-starts_with("Cast"))
colnames(RPKMbasta_2ORDER_percentage) <- sub("_", "", colnames(RPKMbasta_2ORDER_percentage), fixed = TRUE)
colnames(RPKMbasta_2ORDER_percentage) <- gsub("percent", "", colnames(RPKMbasta_2ORDER_percentage))
colnames(RPKMbasta_2ORDER_percentage) <- gsub("group", "", colnames(RPKMbasta_2ORDER_percentage))
custom_orderORDER <- c("Suessiales", "Pelagomonadales", "Peronosporales", "Gymnodiniales", 
                       "Mamiellales", "Perkinsida", "Isochrysidales", "Prymnesiales", "All Other (n = 605)      ", "Unclassified")

RPKM_by_station_normalizedORDER <- RPKMbasta_2ORDER_percentage %>%
  pivot_longer(cols = starts_with("Cast"),  
               names_to = "Sample",  
               values_to = "Relative RPKM within Cast")  %>%
  mutate(Category_new = case_when(
    Sample %in% surface ~ "Surface",      # If Cast is in the surface vector
    Sample %in% oxycline ~ "Oxycline",    # If Cast is in the oxycline vector
    Sample %in% anoxic ~ "Anoxic"))  %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St_1_",
    Sample %in% station2 ~ "St_2_",
    Sample %in% station3 ~ "St_3_"
  )) %>%
  unite(Samp, Station, Sample, sep = "")
RPKM_by_station_normalizedORDER$ORDER <- factor(
  RPKM_by_station_normalizedORDER$ORDER,
  levels = custom_orderORDER)

# Class #####
bCleanCLASS <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "class"))
bClean2CLASS <- bCleanCLASS %>%
  left_join(MarFERReT_polyA_class2, by = "Chr", suffix = c("_bClean", "_MarFERReT")) %>%
  mutate(
    Class = ifelse(is.na(class_bClean), class_MarFERReT, class_bClean)
  ) %>%
  select(Chr, Class)
RPKMbasta_2CLASS <- as.data.frame(RPKMbasta %>%
                                    left_join(bClean2CLASS)) %>%
  select(-c(Col, Chr))
RPKMbasta_2CLASS[RPKMbasta_2CLASS == ""] <- NA
RPKMbasta_2CLASS[RPKMbasta_2CLASS == "unknown"] <- NA
RPKMbasta_2CLASS[is.na(RPKMbasta_2CLASS)] <- "Unclassified"
CLASS <- RPKMbasta_2CLASS %>%
  count(Class) 
CLASS$percent <- CLASS$n / nrow(RPKMbasta_2CLASS) * 100
CLASS <- arrange(CLASS, percent)
CLASS2 <- CLASS[232:240, ]
CLASS3 <- c(Class = "all others (n = 231)", freq = (sum(CLASS$n) - sum(CLASS2$n)), percent = (100 - sum(CLASS2$percent)))
CLASS4 <- rbind(CLASS2, CLASS3)
valid_CLASS <- c("Gammaproteobacteria", "Mamiellophyceae", "Actinopteri", "Alphaproteobacteria", 
                 "Insecta", "Magnoliopsida", "Pelagophyceae", "Dinophyceae", "Unclassified")
RPKMbasta_2CLASS_grouped <- RPKMbasta_2CLASS %>%
  mutate(CLASS_group = ifelse(Class %in% valid_CLASS, Class, "All Other (n = 231)")) %>%
  group_by(CLASS_group) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup() 
RPKMbasta_2CLASS_percentage <- RPKMbasta_2CLASS_grouped %>%
  select(CLASS_group, starts_with("Cast")) %>%  # Select only CLASS_group and Cast columns
  mutate(across(starts_with("Cast"), 
                ~ . / sum(.) * 100,  # Calculate percentage for each column
                .names = "percent_{.col}")) %>%  # Naming the new columns as 'percent_Cast...'
  select(-starts_with("Cast"))
colnames(RPKMbasta_2CLASS_percentage) <- sub("_", "", colnames(RPKMbasta_2CLASS_percentage), fixed = TRUE)
colnames(RPKMbasta_2CLASS_percentage) <- gsub("percent", "", colnames(RPKMbasta_2CLASS_percentage))
colnames(RPKMbasta_2CLASS_percentage) <- gsub("group", "", colnames(RPKMbasta_2CLASS_percentage))
custom_orderCLASS <- c("Dinophyceae", "Pelagophyceae", "Magnoliopsida", "Insecta", "Alphaproteobacteria",
                       "Actinopteri", "Mamiellophyceae", "Gammaproteobacteria", "All Other (n = 231)", "Unclassified")
RPKM_by_station_normalizedCLASS <- RPKMbasta_2CLASS_percentage %>%
  pivot_longer(cols = starts_with("Cast"),  
               names_to = "Sample",  
               values_to = "Relative RPKM within Cast")  %>%
  mutate(Category_new = case_when(
    Sample %in% surface ~ "Surface",      # If Cast is in the surface vector
    Sample %in% oxycline ~ "Oxycline",    # If Cast is in the oxycline vector
    Sample %in% anoxic ~ "Anoxic"))  %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St_1_",
    Sample %in% station2 ~ "St_2_",
    Sample %in% station3 ~ "St_3_"
  )) %>%
  unite(Samp, Station, Sample, sep = "")
RPKM_by_station_normalizedCLASS$CLASS <- factor(
  RPKM_by_station_normalizedCLASS$CLASS,
  levels = custom_orderCLASS)
mean_DINO <- RPKMbasta_2CLASS_grouped %>% 
  filter(CLASS_group == "Dinophyceae") %>%
  select(-CLASS_group)
mean_POLAR <- RPKMbasta_2_grouped %>% 
  filter(Genus_group == "Polarella") %>%
  select(-Genus_group)
meanMean <- rbind(mean_DINO, mean_POLAR) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(Percent = V2/V1 * 100) %>%
  rownames_to_column()
colnames(meanMean) <- c("cast", "V1", "V2", "PercentDino")
non_zero_df2 <- non_zero_df %>% 
  left_join(meanMean) %>%
  select(-c(V1, V2)) %>%
  left_join(mean_Polarella)
column_means <- sapply(non_zero_df2[ ,c(2:4)], mean, na.rm = TRUE)
column_sds <- sapply(non_zero_df2[ ,c(2:4)], sd, na.rm = TRUE)
result <- data.frame(
  Column = names(non_zero_df2[ ,c(2:4)]),
  Mean = column_means,
  SD = column_sds) %>% t() %>% as.data.frame() %>%
  rownames_to_column() %>%
  filter(rowname != "Column") %>%
  mutate(cast = rowname) %>%
  select(-rowname)
non_zero_df3 <- rbind(non_zero_df2, result)
surface <- c('Cast22_60m','Cast59_30m','Cast50_40m','Cast20_54m',
             'Cast34_60m','Cast52_65m','Cast98_25m','Cast85_20m',
             'Cast93_10m','Cast23_40m','Cast55_120m','Cast6_60m',
             'Cast63_10m','Cast73_16m')
oxycline <- c('Cast39_85m','Cast41_90m','Cast97_40m','Cast80_50m',
              'Cast48_90m','Cast31_90m','Cast47_113m','Cast37_92m',
              "Cast59_900m")
anoxic <- c('Cast35_120m','Cast59_250m','Cast27_275m','Cast18_110m',
            'Cast19_90m', 'Cast72_45m','Cast76_36m','Cast101_80m','Cast82_48m')
non_zero_df4 <- non_zero_df3 %>%
  mutate("Oxygen.Regime" = case_when(
    cast %in% surface ~ "Surface",
    cast %in% oxycline ~ "Oxycline",
    cast %in% anoxic ~ "Anoxic")) %>%
  mutate(PercentEuk = as.numeric(PercentEuk))
# For stats to see if oxygen regimes are statistically different
anoxic_data <- non_zero_df4[non_zero_df4$Oxygen.Regime == "Anoxic", ]
shapiro.test(anoxic_data$PercentEuk) # p > 0.05
oxycline_data <- non_zero_df4[non_zero_df4$Oxygen.Regime == "Oxycline", ]
shapiro.test(oxycline_data$PercentEuk) # p > 0.05
surface_data <- non_zero_df4[non_zero_df4$Oxygen.Regime == "Surface", ]
shapiro.test(surface_data$PercentEuk) # p > 0.05
# All data is normally distributed 
leveneTest(PercentEuk ~ Oxygen.Regime, data = non_zero_df4)
# Variances are equal
# Assumptions for parametric anova satisfied 
anova_result <- aov(PercentEuk ~ Oxygen.Regime, data = non_zero_df4)
summary(anova_result)
TukeyHSD(anova_result)
# p > 0.05 in all comparisons. Accept null hypothesis that relative expression is 
# equal in all oxygen regimes 
# add this data to the stats section of the manuscript
write.csv(non_zero_df3, "non_zeroCountsCast.csv", row.names = FALSE)

# Phylum #####
# For 16 Classes, Phylum was "unknown". We fixed this manually via NCBI:
phylum_lookup <- data.frame(
  class = c("Pelagophyceae", "Dinophyceae", "Bigyra", "Eustigmatophyceae", 
            "Cryptophyceae", "Xanthophyceae", "Synurophyceae", "Phaeophyceae", 
            "Aphelidea", "Choanoflagellata", "Ichthyosporea", "Polycystinea", 
            "Raphidophyceae", "Filasterea", "Naldaviricetes", "Acantharea"),  
  phylum = c("Stramenopiles", "Alveolata", "Stramenopiles", "Stramenopiles",
             "unknown", "Stramenopiles", "Stramenopiles", "Stramenopiles",
             "Aphelidea", "unknown", "unknown", "Rhizaria",
             "Stramenopiles", "unknown", "unknown", "Rhizaria"))
bCleanPHYLUM <- basta %>%
  separate(V2, c("kingdom", "phylum", 
                 "class", "order", "family", 
                 "genus", "species"), ";",
           extra = "drop", fill = "right") %>%
  rename_at("V1",~"Chr") %>%
  select(c("Chr", "phylum", "class")) %>%
  left_join(phylum_lookup, by = "class") %>%  # Join by class
  mutate(phylum = ifelse(phylum.x == "unknown", phylum.y, phylum.x)) %>%  # Update phylum
  select(-phylum.y, -phylum.x) %>%
  select(-class)
bClean2PHYLUM <- bCleanPHYLUM %>%
  left_join(MarFERReT_polyA_phylum2, by = "Chr", suffix = c("_bClean", "_MarFERReT")) %>%
  mutate(
    Phylum = ifelse(is.na(phylum_bClean), phylum_MarFERReT, phylum_bClean)
  ) %>%
  select(Chr, Phylum)
RPKMbasta_2PHYLUM <- as.data.frame(RPKMbasta %>%
                                    left_join(bClean2PHYLUM)) %>%
  select(-c(Col, Chr))
RPKMbasta_2PHYLUM[RPKMbasta_2PHYLUM == ""] <- NA
RPKMbasta_2PHYLUM[RPKMbasta_2PHYLUM == "unknown"] <- NA
RPKMbasta_2PHYLUM[is.na(RPKMbasta_2PHYLUM)] <- "Unclassified"
PHYLUM <- RPKMbasta_2PHYLUM %>%
  count(Phylum) 
PHYLUM$percent <- PHYLUM$n / nrow(RPKMbasta_2PHYLUM) * 100
PHYLUM <- arrange(PHYLUM, percent)
PHYLUM2 <- PHYLUM[127:135, ]
PHYLUM3 <- c(Phylum = "all others (n = 126)", freq = (sum(PHYLUM$n) - sum(PHYLUM2$n)), percent = (100 - sum(PHYLUM2$percent)))
PHYLUM4 <- rbind(PHYLUM2, PHYLUM3)
valid_PHYLUM <- c("Haptophyta", "Proteobacteria", "Stramenopiles", "Arthropoda", 
                  "Chordata", "Streptophyta", "Oomycota", " Alveolata", "Unclassified")
RPKMbasta_2PHYLUM_grouped <- RPKMbasta_2PHYLUM %>%
  mutate(PHYLUM_group = ifelse(Phylum %in% valid_PHYLUM, Phylum, "All Other (n = 126)")) %>%
  group_by(PHYLUM_group) %>%
  summarize(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup() 
RPKMbasta_2PHYLUM_percentage <- RPKMbasta_2PHYLUM_grouped %>%
  select(PHYLUM_group, starts_with("Cast")) %>%  # Select only PHYLUM_group and Cast columns
  mutate(across(starts_with("Cast"), 
                ~ . / sum(.) * 100,  # Calculate percentage for each column
                .names = "percent_{.col}")) %>%  # Naming the new columns as 'percent_Cast...'
  select(-starts_with("Cast"))
colnames(RPKMbasta_2PHYLUM_percentage) <- sub("_", "", colnames(RPKMbasta_2PHYLUM_percentage), fixed = TRUE)
colnames(RPKMbasta_2PHYLUM_percentage) <- gsub("percent", "", colnames(RPKMbasta_2PHYLUM_percentage))
colnames(RPKMbasta_2PHYLUM_percentage) <- gsub("group", "", colnames(RPKMbasta_2PHYLUM_percentage))
custom_orderPHYLUM <- c("Haptophyta", "Proteobacteria", "Stramenopiles", "Arthropoda", 
                        "Chordata", "Streptophyta", "Oomycota", " Alveolata", "All Other (n = 126)", "Unclassified")
RPKM_by_station_normalizedPHYLUM <- RPKMbasta_2PHYLUM_percentage %>%
  pivot_longer(cols = starts_with("Cast"),  
               names_to = "Sample",  
               values_to = "Relative RPKM within Cast")  %>%
  mutate(Category_new = case_when(
    Sample %in% surface ~ "Surface",      # If Cast is in the surface vector
    Sample %in% oxycline ~ "Oxycline",    # If Cast is in the oxycline vector
    Sample %in% anoxic ~ "Anoxic"))  %>%
  mutate(Station = case_when(
    Sample %in% station1 ~ "St_1_",
    Sample %in% station2 ~ "St_2_",
    Sample %in% station3 ~ "St_3_"
  )) %>%
  unite(Samp, Station, Sample, sep = "")
RPKM_by_station_normalizedPHYLUM$PHYLUM <- factor(
  RPKM_by_station_normalizedPHYLUM$PHYLUM,
  levels = custom_orderPHYLUM)

# Stacked Bar Plots #####
RPKMnormal_plot <- ggplot(RPKM_by_station_normalized, aes(x = Samp, y = `Relative RPKM within Cast`, fill = `Genus`)) + 
  scale_fill_manual(values = c("goldenrod1", "lavender", "forestgreen", "dodgerblue",
                               "goldenrod4", "darkgreen", "darkblue",
                               "yellow", "green", "lightblue")) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +  # Corrected theme element
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  facet_grid( ~ Category_new,
              scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(80, 100, by = 5)) + 
  coord_cartesian(ylim = c(80, 100)) +
  bw + theme(axis.text.x = element_markdown(angle = 55, hjust = 1, size = 12),
             axis.title.x = element_markdown(size = 14),
             axis.text.y = element_markdown(size = 14),
             legend.text = element_text(size = 14),
             legend.title = element_text(size = 14),
             strip.background = element_blank(),
             strip.text = element_blank(),
             plot.margin = margin(t = 10, r = 5, b = 5, l = 20)) +
  scale_x_discrete(labels = function(x) bold_sample_strings(x)) +
  labs(title =  NULL,
       x = "Sample       ",
       y = "   ",
       fill = "Taxonomy - Genus Level")
RPKMnormal_plotFAMILY <- ggplot(RPKM_by_station_normalizedFAMILY, aes(x = Samp, y = `Relative RPKM within Cast`, fill = `family`)) + 
  scale_fill_manual(values = c("goldenrod1", "lavender", "forestgreen", "dodgerblue",
                               "goldenrod4", "darkgreen", "darkblue",
                               "yellow",  "green", "lightblue")) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +  # Corrected theme element
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  facet_grid( ~ Category_new,
              scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(80, 100, by = 5)) + 
  coord_cartesian(ylim = c(80, 100)) +
  bw + theme(axis.text.x = element_markdown(angle = 55, hjust = 1, size = 12),
             axis.title.x = element_markdown(size = 14),
             axis.text.y = element_markdown(size = 14),
             legend.text = element_text(size = 14),
             legend.title = element_text(size = 14),
             strip.background = element_blank(),
             strip.text = element_blank(),
             plot.margin = margin(t = 10, r = 5, b = 5, l = 20)) +
  scale_x_discrete(labels = function(x) bold_sample_strings(x)) +
  labs(title = NULL,
       x = "Sample       ",
       y = "   ",
       fill = "Taxonomy - Family Level")
RPKMnormal_plotORDER <- ggplot(RPKM_by_station_normalizedORDER, aes(x = Samp, y = `Relative RPKM within Cast`, fill = `ORDER`)) + 
  scale_fill_manual(values = c("lavender", "goldenrod1", "forestgreen", "dodgerblue",
                               "goldenrod4", "darkgreen", "darkblue",
                               "yellow",  "green", "lightblue")) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +  # Corrected theme element
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  facet_grid( ~ Category_new,
              scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(60, 100, by = 5)) + 
  coord_cartesian(ylim = c(60, 100)) +
  bw + theme(axis.text.x = element_markdown(angle = 55, hjust = 1, size = 12),
             axis.title.x = element_markdown(size = 14),
             axis.text.y = element_markdown(size = 14),
             legend.text = element_text(size = 14),
             legend.title = element_text(size = 14),
             strip.background = element_blank(),
             strip.text = element_blank(),
             plot.margin = margin(t = 10, r = 5, b = 5, l = 20)) +
  scale_x_discrete(labels = function(x) bold_sample_strings(x)) +
  labs(title = NULL,
       x = "Sample       ",
       y = "   ",
       fill = "Taxonomy - Order Level")
RPKMnormal_plotCLASS <- ggplot(RPKM_by_station_normalizedCLASS, aes(x = Samp, y = `Relative RPKM within Cast`, fill = `CLASS`)) + 
  scale_fill_manual(values = c("lavender", "goldenrod1", "forestgreen", "dodgerblue",
                               "goldenrod4", "darkgreen", "darkblue",
                               "yellow",  "green", "lightblue")) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +  # Corrected theme element
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  facet_grid( ~ Category_new,
              scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(60, 100, by = 5)) + 
  coord_cartesian(ylim = c(60, 100)) +
  bw + theme(axis.text.x = element_markdown(angle = 55, hjust = 1, size = 12),
             axis.title.x = element_markdown(size = 14),
             axis.text.y = element_markdown(size = 14),
             legend.text = element_text(size = 14),
             legend.title = element_text(size = 14),
             strip.background = element_blank(),
             strip.text = element_text(size = 14),
             plot.margin = margin(t = 10, r = 5, b = 5, l = 20)) +
  scale_x_discrete(labels = function(x) bold_sample_strings(x)) +
  labs(title = NULL,
       x = "Sample       ",
       y = "   ",
       fill = "Taxonomy - Class Level")
RPKMnormal_plotPHYLUM <- ggplot(RPKM_by_station_normalizedPHYLUM, aes(x = Samp, y = `Relative RPKM within Cast`, fill = `PHYLUM`)) + 
  scale_fill_manual(values = c("goldenrod1", "forestgreen", "dodgerblue",
                               "goldenrod4", "darkgreen", "darkblue",
                               "yellow",  "lavender", "green", "lightblue")) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +  # Corrected theme element
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  facet_grid( ~ Category_new,
              scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(60, 100, by = 5)) + 
  coord_cartesian(ylim = c(60, 100)) +
  bw + theme(axis.text.x = element_markdown(angle = 55, hjust = 1, size = 12),
             axis.title.x = element_markdown(size = 14),
             axis.text.y = element_markdown(size = 14),
             legend.text = element_text(size = 14),
             legend.title = element_text(size = 14),
             strip.text = element_text(size = 14),
             strip.background = element_blank(),
             plot.margin = margin(t = 10, r = 5, b = 5, l = 20)) +
  scale_x_discrete(labels = function(x) bold_sample_strings(x)) +
  labs(title = NULL,
       x = "Sample       ",
       y = "   ",
       fill = "Taxonomy - Phylum Level")

# Depth Profiles with Trophic-Mode Genes #####
RPKM_TrophicMode_Depth <- PolarellaDEG %>%
  filter(ID %in% c(phototrophy, "PTHR14218") | Chr %in% deepLocDBCAN$Chr) %>%
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
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 80, y = 0, label = "(A)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Carbonic \nAndyrase", x = "RPKM", y = "Depth (m)") + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
ChlAB <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[6,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 350, y = 0, label = "(B)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Chl A/B Binding \nProtein", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
RubISCo <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[10,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue"), guide = "none")  +
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
  annotate("text", x = 400, y = 0, label = "(C)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "RuBiSCo", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
Rhodopsin <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[9,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 1750, y = 0, label = "(D)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Rhodopsin", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
ATPsyn <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[2,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 40, y = 0, label = "(E)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "ATP Synthase \nGamma Subunit", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
# Heterotrophic Genes
ABCTrans <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[1,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 5, y = 0, label = "(F)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "ABC Transporter", x = "RPKM", y = "Depth (m)") + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
GH3 <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[3,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 10, y = 0, label = "(G)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Beta-Glucosidase \n(GH3)", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
GH18 <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[5,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17)) +
  scale_color_manual(values = c("goldenrod1", "forestgreen", "blue")) +
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
  annotate("text", x = 25, y = 0, label = "(H)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Chitinase \n(GH18)", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom", 
        legend.direction = "horizontal",
        legend.key.size = unit(5, "cm"),
        legend.text = element_text(size = 22),  
        legend.title = element_text(size = 22))
GH1 <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[7,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  annotate("text", x = 10, y = 0, label = "(I)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Glycosyl \nHydrolase (GH1)", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
peptidase <- RPKM_TrophicMode_Depth %>%
  filter(Name == as.character(subsetRPKM_trophic[8,])) %>%
  ggplot() +
  geom_path(aes(RPKM_All, Depth, color = Station), size = 1.1) + 
  geom_point(aes(RPKM_All, Depth, color = Station, shape = Station), size = 4) +
  scale_shape_manual(name = "Station", values = c(15, 16,17), guide = "none") +
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
  scale_x_continuous(lim=c(0,650)) + 
  annotate("text", x = 650, y = 0, label = "(J)", 
           hjust = 1, vjust = 1, size = 6) +
  labs(title = "Tripeptidyl \npeptidase I", x = "RPKM", y = NULL) + bw +
  theme(axis.text = element_text(size = 23),
        axis.text.x = element_markdown(angle = 45, hjust = 1, size = 23),
        axis.title = element_text(size = 23),  
        plot.title = element_text(size = 25, face = "bold"))
RPKM_TrophicMode_Depth_Plot <- (CarbonicAn | ChlAB  | RubISCo | Rhodopsin | ATPsyn) /
  (ABCTrans | GH3 | GH18 | GH1 | peptidase)
ggsave("TMDP_REVISE.tiff", RPKM_TrophicMode_Depth_Plot, height = 12, width = 18, units = "in", dpi = 300)
ggsave("TMDP_FINAL_review.png", RPKM_TrophicMode_Depth_Plot, height = 12, width = 18, units = "in", dpi = 300)
ggsave("TMDP_FINAL_review.pdf", RPKM_TrophicMode_Depth_Plot, height = 12, width = 18, units = "in", dpi = 300)

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
ggsave("nutDepthProfiles.tiff", nutDepthProfiles, height = 4.5, width = 9, units = "in", dpi = 300)
ggsave("nutDepthProfiles_FINAL_review.png", nutDepthProfiles, height = 4.5, width = 9, units = "in", dpi = 300)
ggsave("nutDepthProfiles_FINAL_review.pdf", nutDepthProfiles, height = 4.5, width = 9, units = "in", dpi = 300)
