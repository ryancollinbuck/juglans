########################################################  ANALYSES  ###############################################
qrsh -l h_rt=24:00:00,h_vmem=5G -pe shared 4
module load gcc/10.2.0
module load gdal/3.1.3
module load geos/3.11.1
module load proj/7.1.1
module load sqlite/3.33.0
module load curl/8.4.0
module load R/4.3.0
module load pandoc/2.17.1.1
cd project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses


###### Plink - make bed/bim files while filtering for LD######      === named "plinkBed.sh" ===
module load plink
plink --vcf meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Juglansfiltered
##### vcf 2 .raw(plink) #####  #for RDA input
plink --vcf meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Juglansfiltered.prune.in --make-bed --pca var-wts --recode A --out Juglansfiltered
##### vcf to LD trimmed vcf #####
plink --vcf meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Juglansfiltered.prune.in --recode vcf --out Juglansfiltered
#produced 596,513 SNPs #genotype rate 0.927 (after removing 0 ind =142 indv)


##remove Jcal.S.SBO.10 and Jcal.S.SBO.16
vcftools --gzvcf Juglans_jhin.repeatsOut.ef.vcf.gz --remove-indv Jcal.S.SBO.10 --remove-indv Jcal.S.SBO.16 --min-alleles 2 --max-alleles 2 --min-meanDP 5 --minDP 5 --max-missing 0.9 --maf 0.01 --recode --stdout | gzip -c > Jcal.S.SBO.10n16rem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz

###### Plink - make bed/bim files while filtering for LD######      === named "plinkBed.sh" ===
module load plink
plink --vcf Jcal.S.SBO.10n16rem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Juglansfiltered
##### vcf 2 .raw(plink) #####  #for RDA input
plink --vcf Jcal.S.SBO.10n16rem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Juglansfiltered.prune.in --make-bed --pca var-wts --recode A --out Juglansfiltered
##### vcf to LD trimmed vcf #####
plink --vcf Jcal.S.SBO.10n16rem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Juglansfiltered.prune.in --recode vcf --out Juglansfiltered
#produced 714,517 SNPs #genotype rate 0.933 (after removing 2 ind =140 indv)



module load plink
plink --vcf outliersrem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out outliersrem_Juglansfiltered
##### vcf 2 .raw(plink) #####  #for RDA input
plink --vcf outliersrem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract outliersrem_Juglansfiltered.prune.in --make-bed --pca var-wts --recode A --out outliersrem_Juglansfiltered
##### vcf to LD trimmed vcf #####
plink --vcf outliersrem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract outliersrem_Juglansfiltered.prune.in --recode vcf --out outliersrem_Juglansfiltered
#produced 417,463 SNPs #genotype rate 0.933 (after removing 19 ind = 123 indv)


module load plink
plink --vcf moreoutliersrem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out moreoutliersrem_Juglansfiltered
##### vcf 2 .raw(plink) #####  #for RDA input
plink --vcf moreoutliersrem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract moreoutliersrem_Juglansfiltered.prune.in --make-bed --pca var-wts --recode A --out moreoutliersrem_Juglansfiltered
##### vcf to LD trimmed vcf #####
plink --vcf moreoutliersrem_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract moreoutliersrem_Juglansfiltered.prune.in --recode vcf --out moreoutliersrem_Juglansfiltered
#produced 389,867 SNPs #genotype rate 0.938 (after removing 23 ind = 119 indv)


#####SPLIT BY SPECIES#######
module load plink
plink --vcf RIV3nSD33nA133remJhinONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out RIV3nSD33nA133rem_JhinONLY_filtered
##### vcf 2 .raw(plink) #####  #for RDA input
plink --vcf RIV3nSD33nA133remJhinONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract RIV3nSD33nA133rem_JhinONLY_filtered.prune.in --make-bed --pca var-wts --recode A --out RIV3nSD33nA133rem_JhinONLY_filtered
##### vcf to LD trimmed vcf #####
plink --vcf RIV3nSD33nA133remJhinONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract RIV3nSD33nA133rem_JhinONLY_filtered.prune.in --recode vcf --out RIV3nSD33nA133rem_JhinONLY_filtered
#produced 100,528 SNPs #genotype rate 0.938 (59 indv)

module load plink
plink --vcf SCA10nSON1nCON4nLA3n38nSD2n72remJcalONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered
##### vcf 2 .raw(plink) #####  #for RDA input
plink --vcf SCA10nSON1nCON4nLA3n38nSD2n72remJcalONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.prune.in --make-bed --pca var-wts --recode A --out SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered
##### vcf to LD trimmed vcf #####
plink --vcf SCA10nSON1nCON4nLA3n38nSD2n72remJcalONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.prune.in --recode vcf --out SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered
#produced 715,275 SNPs #genotype rate 0.977 (39 indv)



###### Rename chromosomes ########
sed -i.bak 's/CM084183\.1/1/g' Juglansfiltered.bim
sed -i.bak 's/CM084184\.1/2/g' Juglansfiltered.bim
sed -i.bak 's/CM084185\.1/3/g' Juglansfiltered.bim
sed -i.bak 's/CM084186\.1/4/g' Juglansfiltered.bim
sed -i.bak 's/CM084187\.1/5/g' Juglansfiltered.bim
sed -i.bak 's/CM084188\.1/6/g' Juglansfiltered.bim
sed -i.bak 's/CM084189\.1/7/g' Juglansfiltered.bim
sed -i.bak 's/CM084190\.1/8/g' Juglansfiltered.bim
sed -i.bak 's/CM084191\.1/9/g' Juglansfiltered.bim
sed -i.bak 's/CM084192\.1/10/g' Juglansfiltered.bim
sed -i.bak 's/CM084193\.1/11/g' Juglansfiltered.bim
sed -i.bak 's/CM084194\.1/12/g' Juglansfiltered.bim
sed -i.bak 's/CM084195\.1/13/g' Juglansfiltered.bim
sed -i.bak 's/CM084196\.1/14/g' Juglansfiltered.bim
sed -i.bak 's/CM084197\.1/15/g' Juglansfiltered.bim
sed -i.bak 's/CM084198\.1/16/g' Juglansfiltered.bim


sed -i.bak 's/CM084183\.1/1/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084184\.1/2/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084185\.1/3/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084186\.1/4/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084187\.1/5/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084188\.1/6/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084189\.1/7/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084190\.1/8/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084191\.1/9/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084192\.1/10/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084193\.1/11/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084194\.1/12/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084195\.1/13/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084196\.1/14/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084197\.1/15/g' outliersrem_Juglansfiltered.bim
sed -i.bak 's/CM084198\.1/16/g' outliersrem_Juglansfiltered.bim

sed -i.bak 's/CM084183\.1/1/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084184\.1/2/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084185\.1/3/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084186\.1/4/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084187\.1/5/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084188\.1/6/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084189\.1/7/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084190\.1/8/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084191\.1/9/g'  SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084192\.1/10/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084193\.1/11/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084194\.1/12/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084195\.1/13/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084196\.1/14/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084197\.1/15/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim
sed -i.bak 's/CM084198\.1/16/g' SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bim

###### Heterozygosity ######     === named "hetCalc.sh" ===
module load vcftools/0.1.16
vcftools --vcf manyrem_Qdoufiltered.vcf --out het_manyrem_Qdoufiltered --het 

#rename samples in vcf to remove duplicate naming
bcftools reheader -s Jcal_new_sample_names.txt -o SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered_renamed.vcf SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.vcf


##### vcf 2 .snp #####  #for GF input
#module load vcftools
vcftools --vcf manyrem_Qdoufiltered_renamed.vcf --012 --out manyrem_Qdoufiltered_renamed_snp
cut -f2- manyrem_Qdoufiltered_renamed_snp.012 | sed 's/-1/NA/g' >manyrem_Qdoufiltered_renamed_snp.temp
tr -d '\t' <manyrem_Qdoufiltered_renamed_snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - manyrem_Qdoufiltered_renamed_snp.012.indv) <(echo "" | cat header - manyrem_Qdoufiltered_renamed_snp.temp) > manyrem_Qdoufiltered_renamed_snp.forR
rm header manyrem_Qdoufiltered_renamed_snp.temp


####### admixture loop #######
for i in {1..5}; do admixture --cv ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/outliersrem_Juglansfiltered.bed $i > ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/outliersrem_Juglansfiltered.$i.log.out; done
grep -h CV ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/outliersrem_Juglansfiltered.*log.out


for i in {6..10}; do admixture --cv ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.bed $i > ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.$i.log.out; done
grep -h CV ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.*log.out




###PCA in R
library(ggplot2)
pca<-read.table("Juglansfiltered.eigenvec")
eigenval<-scan("Juglansfiltered.eigenval")
pca<-pca[,-1]
names(pca)[1]<-"ind"
names(pca)[2:ncol(pca)]<-paste0("PC",1:(ncol(pca)-1))
pve<-data.frame(PC=1:20,pve=eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
#PC2 + PC3
b <- ggplot(pca, aes(PC2, PC3)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

#with labels
b <- ggplot(pca, aes(PC1, PC2,label=ind)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+geom_text(hjust=1,vjust=0)

b <- ggplot(pca, aes(PC2, PC3,label=ind)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))+geom_text(hjust=0,vjust=0)









###PCA in R
library(ggplot2)
pca<-read.table("SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.eigenvec")
eigenval<-scan("SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.eigenval")
pca<-pca[,-1]
names(pca)[1]<-"ind"
names(pca)[2:ncol(pca)]<-paste0("PC",1:(ncol(pca)-1))
pve<-data.frame(PC=1:20,pve=eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
#PC2 + PC3
b <- ggplot(pca, aes(PC2, PC3)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

#with labels
b <- ggplot(pca, aes(PC1, PC2,label=ind)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+geom_text(hjust=1,vjust=0)

b <- ggplot(pca, aes(PC2, PC3,label=ind)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))+geom_text(hjust=0,vjust=0)

#color by latitude
library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
JCenv<-data.frame(JCenv)
colnames(JCenv) <- c("sample","lon", "lat")
pca <- merge(pca, JCenv[, c("sample", "lon")], by.x = "ind", by.y = "sample", all.x = TRUE)
b <- ggplot(pca, aes(PC1, PC2, colour = lon)) +
  geom_point(size = 3) +
  coord_equal() +
  theme_light()
b + 
  scale_colour_viridis_c(option = "turbo") +
  labs(colour = "Longitude") +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))









######## turn into raster
# Load the required packages
library(gstat)
library(sp)
library(raster)

library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/TNC_Qdou_coordsforGF_manyrem.csv")
geoQD<-subset(JCenv, select=c(x,y))
colnames(geoQD)<-c("Longitude","Latitude")
geoQD<-data.frame(geoQD)
pcaQD<-read.table("manyrem_Qdoufiltered.eigenvec")
eigenvalQD<-scan("manyrem_Qdoufiltered.eigenval")
pcaQD<-pcaQD[,-1]
names(pcaQD)[1]<-"ind"
names(pcaQD)[2:ncol(pcaQD)]<-paste0("PC",1:(ncol(pcaQD)-1))
pcaQD<-pcaQD[,-1]
pca_geo <- cbind(geoQD, pcaQD)
coordinates(pca_geo) <- ~Longitude + Latitude
##create a prediction grid
# Define the extent of the region you want to cover
x.range <- c(-123.86,-117.29)
y.range <- c(34.28,41.30)
# Create a grid of prediction points
grd <- expand.grid(Longitude = seq(from = x.range[1], to = x.range[2], by = 0.001),
                   Latitude = seq(from = y.range[1], to = y.range[2], by = 0.001))
# Convert to spatial points
coordinates(grd) <- ~ Longitude + Latitude
gridded(grd) <- TRUE  # Convert the grid to a gridded object
# Fit a variogram model (you can try different models such as "Sph", "Exp", etc.)
vgm_model <- variogram(PC1 ~ 1, data = pca_geo)
fit_model <- fit.variogram(vgm_model, model = vgm("Sph"))
# Perform kriging
kriging_result <- krige(PC1 ~ 1, pca_geo, grd, model = fit_model)
# Convert to raster for easy plotting
raster_kriged <- raster(kriging_result)
writeRaster(raster_kriged, "QDrangePCA.tif", "GTiff", overwrite=TRUE)


###############################################################################



###Het and Fis in excel

###Fst in vcftools
#make pop.txt files with one column containing the list of indivuals to keep for that pop
#use weighted print out
vcftools --vcf manyrem_Qdoufiltered_renamed.vcf --weir-fst-pop fst_DyeCreek.txt --weir-fst-pop fst_Randall.txt --out QDF_DCvR


########################################################  ANALYSES  ###############################################
qrsh -l h_rt=24:00:00,h_vmem=25G -pe shared 12
module load gcc/10.2.0
module load gdal/3.1.3
module load geos/3.11.1
module load proj/7.1.1
module load sqlite/3.33.0
module load curl/8.4.0
module load R/4.3.0
module load pandoc/2.17.1.1
module load cmake/3.30.0
cd project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONlY

###initial GF clim varibale importance on cluster

library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
JCenv<-data.frame(JCenv)
Coordinates<-subset(JCenv, select=c(x,y))
Coordinates<-data.frame(Coordinates)


library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA
remove_NAs_stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove_NAs_stack(ras_8110)

sample.coord.sp <- vect(Coordinates, geom = c("x", "y"), crs = crs(ras_8110))
clim.points_8110 <- extract(ras_8110, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
clim.points_8110<-data.frame(clim.points_8110)





###RDA to select climate variables and adaptive SNPs###
library(psych)
library(vegan)
library(readr)
library(adegenet)
library(parallel)

Env_scale <- scale(clim.points_8110, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
scale_cp5180 <- attr(Env_scale, 'scaled:scale')
center_cp5180 <- attr(Env_scale, 'scaled:center')
Env_scale<-as.data.frame(Env_scale)
row.names(Env_scale)<-c(JCenv$sample)
#colnames(Env_scale)<-c("aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
subEnv_scale<-subset(Env_scale,select=-c(ID))

colnames(Coordinates) <- c("Latitude", "Longitude")
Coordinates$LatxLon<-Coordinates$Latitude * Coordinates$Longitude


snp <- read.PLINK("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered.raw",parallel = TRUE)
library(parallel)
cl<-makeCluster(detectCores())
snp.imp <- parApply(cl,snp, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) #remove NAs and replace with most common genotype


pca.res <- prcomp(snp.imp, center = TRUE, scale. = TRUE)
eigenval <- (pca.res$sdev)^2
pve <- eigenval / sum(eigenval) * 100
pca <- data.frame(ind = rownames(snp.imp),
                     pca.res$x[,1:10])  # first 10 PCs

pve.df <- data.frame(PC = 1:length(pve), pve = pve)

library(ggplot2)
ggplot(pve.df[1:20,], aes(PC, pve)) +
  geom_bar(stat = "identity") +
  ylab("Percentage variance explained") +
  theme_light()

ggplot(pca, aes(PC1, PC2, label = ind)) +
  geom_point() +
  geom_text(size = 3, vjust = -0.5) +
  theme_classic()

PCs<-subset(pca,select=PC1)

Variables <- data.frame(PCs,Coordinates,subEnv_scale)








#pRDA
pRDAfull <- rda(snp.imp ~ PC1+Latitude+Longitude+LatxLon+aet+awc+cwd+depth+pct_clay+ph+ppt_djf+ppt_jja+terrain+tmn,  data=Variables,parallel=TRUE)
RsquareAdj(pRDAfull)
#anova(pRDAfull,parallel=TRUE)

## Pure climate model
pRDAclim <- rda(snp.imp ~ aet+awc+cwd+depth+pct_clay+ph+ppt_djf+ppt_jja+terrain+tmn + Condition(PC1+Latitude+Longitude+LatxLon),  data=Variables,parallel=TRUE)
RsquareAdj(pRDAclim)
#anova(pRDAclim,parallel=TRUE)

##Pure geography model 
pRDAgeog <- rda(snp.imp ~ Latitude+Longitude+LatxLon + Condition(PC1+aet+awc+cwd+depth+pct_clay+ph+ppt_djf+ppt_jja+terrain+tmn),  data=Variables,parallel=TRUE)
RsquareAdj(pRDAgeog)
#anova(pRDAgeog,parallel=TRUE)

##Pure population structure model
pRDAstruct <- rda(snp.imp ~ PC1 + Condition(aet+awc+cwd+depth+pct_clay+ph+ppt_djf+ppt_jja+terrain+tmn+Latitude+Longitude+LatxLon), data=Variables,parallel=TRUE)
RsquareAdj(pRDAstruct)
#anova(pRDAstruct,parallel=TRUE)


pRDAfull
pRDAclim
pRDAgeog
pRDAstruct

library(psych)
library(vegan)
library(readr)
library(adegenet)
library(parallel)
library(robust)
library(qvalue)
library(ggplot2)
library(dplyr)
RDA_env_constrained <- rda(snp.imp ~ aet+awc+cwd+depth+pct_clay+ph+ppt_djf+ppt_jja+terrain+tmn+Condition(PC1+Latitude+Longitude+LatxLon),  data=Variables,parallel=TRUE)
ggplot() +geom_line(aes(x=c(1:length(RDA_env_constrained$CCA$eig)), y=as.vector(RDA_env_constrained$CCA$eig)), linetype="dotted",size = 1.5, color="darkgrey") +geom_point(aes(x=c(1:length(RDA_env_constrained$CCA$eig)), y=as.vector(RDA_env_constrained$CCA$eig)), size = 3,color="darkgrey") +scale_x_discrete(name = "Ordination axes", limits=c(1:10)) +ylab("Inertia") +theme_bw()
screeplot(RDA_env_constrained, main="Eigenvalues of constrained axes")
#***************************************rdadapt function**************************************
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
#********************************************************************************************
rdadapt_env_constrained <- rdadapt(RDA_env_constrained, 2)
thres_env_unc <- 0.01/length(rdadapt_env_constrained$p.values)
outliers_constrained <- data.frame(Loci = colnames(snp.imp)[which(rdadapt_env_constrained$p.values<thres_env_unc)], p.value = rdadapt_env_constrained$p.values[which(rdadapt_env_constrained$p.values<thres_env_unc)], contig = unlist(lapply(strsplit(colnames(snp.imp)[which(rdadapt_env_constrained$p.values<thres_env_unc)], split = "_"), function(x) x[1])))
outliers_constrained <- outliers_constrained[order(outliers_constrained$contig, outliers_constrained$p.value),]
outliers_rdadapt_env_constrained <- as.character(outliers_constrained$Loci[!duplicated(outliers_constrained$contig)])
#10,819 SNPs
#--------------------------------------------------------choosing the top 208 loci (0.01%)-------------------------------------------------------------------
outliers_ordered_constrained<-outliers_constrained[order(outliers_constrained$p.value,decreasing=F),]
#########top 1% of outliers)##########
#Calculate the number of rows corresponding to the top 10% of the dataset
total_rows <- nrow(outliers_ordered_constrained)
top_percent <- ceiling(0.1 * total_rows)

# Select the top 10% of outliers based on p.value
top_outliers_constrained_1percent <- top_n(outliers_ordered_constrained, top_percent, -p.value)
#1083
###########################################

locus_scores <- scores(RDA_env_constrained, choices=c(1:4), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers_constrained$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%top_outliers_constrained_1percent$Loci] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env_constrained, choices=c(1,2), display="bp")) # pull the biplot scores
# Constrained eigenvalues
eig_constrained <- RDA_env_constrained$CCA$eig
var_exp <- eig_constrained / sum(eig_constrained) * 100
RDA1_perc <- round(var_exp[1], 2)
RDA2_perc <- round(var_exp[2], 2)

ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab(paste0("RDA 1 (", RDA1_perc, "%)")) +
  ylab(paste0("RDA 2 (", RDA2_perc, "%)")) +  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

Outliers <- rep("Neutral", length(colnames(snp.imp)))
Outliers[colnames(snp.imp)%in%outliers_constrained$Loci] <- "All outliers"
Outliers[colnames(snp.imp)%in%top_outliers_constrained_1percent$Loci] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(snp.imp)), 
                           pvalues = rdadapt_env_constrained$p.values, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env_unc), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))
  
write.csv(top_outliers_constrained_1percent$contig, "JcalONLY_RDA10p_RoseVariables_terra_PC1_LatxLon_outliers_con.csv")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
########plot samples onto snp RDA colored by locality
row.names(Variables)<-c(JCenv$sample)
Variables$Group <- gsub("QD_([A-Z]{2})_([A-Z]{3})_.*", "\\1_\\2", rownames(Variables))
table(Variables$Group)  # See count of each unique group
# Define colors for each group
group_colors <- c("AC_TMA" = "maroon1", "CA_KRN" = "blue", "DB_SLO" = "green4",
                  "EI_SLO" = "green4", "JH_FRZ" = "green", "JH_MON" = "powderblue",
                  "JH_MPS" = "orange", "JH_NAP" = "red4", "JH_PLC" = "red",
                  "SB_SBR" = "orchid4", "SH_TMA"="maroon1","ZP_KRN"="blue","QD_CA_KR_100"="blue")

# Plot RDA as usual
plot(RDA_env_unconstrained, scaling = 3, display = c("sites", "species","bp"))
sample_scores <- scores(RDA_env_unconstrained, choices = c(1,2), display = "sites", scaling = 3)
snp_scores <- scores(RDA_env_unconstrained, choices = c(1,2), display = "species", scaling = 3)
points(snp_scores[,1], snp_scores[,2], col = "grey", pch = 16, cex = 1.2)
points(sample_scores[,1], sample_scores[,2], col = group_colors[Variables$Group], pch = 16, cex = 1.2)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------



####### LFMM ########
##https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
library(psych)
library(vegan)
library(readr)
library(adegenet)
library(parallel)
library(lfmm)
library(qvalue)

pred.pca<-rda(clim.points_8110[,-1],scale=T)
screeplot(pred.pca, main = "Screeplot of Wolf Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")
#2 PCs

gen.pca <- rda(snp.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")


pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)
pred.PC2 <- scores(pred.pca, choices=2, display="sites", scaling=0)
#pred.PC3 <- scores(pred.pca, choices=3, display="sites", scaling=0)

pc1.lfmm<-lfmm_ridge(Y=snp.imp,X= pred.PC1 ,K=1)
pc1.pv<-lfmm_test(Y=snp.imp,X=pred.PC1,lfmm=pc1.lfmm,calibrate="gif")
pc1.pv$gif
pc1.qv<-qvalue(pc1.pv$pvalue)$qvalues
length(which(pc1.qv<0.1))
length(which(pc1.qv<0.01))
length(which(pc1.qv<0.001))

pc2.lfmm<-lfmm_ridge(Y=snp.imp,X= pred.PC2 ,K=1)
pc2.pv<-lfmm_test(Y=snp.imp,X=pred.PC2,lfmm=pc2.lfmm,calibrate="gif")
pc2.pv$gif
pc2.qv<-qvalue(pc2.pv$pvalue)$qvalues
length(which(pc2.qv<0.1))
length(which(pc2.qv<0.01))
length(which(pc2.qv<0.001))

#pc3.lfmm<-lfmm_ridge(Y=snp.imp,X= pred.PC3 ,K=1)
#pc3.pv<-lfmm_test(Y=snp.imp,X=pred.PC3,lfmm=pc3.lfmm,calibrate="gif")
#pc3.pv$gif
#pc3.qv<-qvalue(pc3.pv$pvalue)$qvalues
#length(which(pc3.qv<0.1))
#length(which(pc3.qv<0.01))
#length(which(pc3.qv<0.001))
#
#pc4.lfmm<-lfmm_ridge(Y=snp.imp,X= pred.PC4 ,K=1)
#pc4.pv<-lfmm_test(Y=snp.imp,X=pred.PC4,lfmm=pc4.lfmm,calibrate="gif")
#pc4.pv$gif
#pc4.qv<-qvalue(pc4.pv$pvalue)$qvalues
#length(which(pc4.qv<0.1))
#length(which(pc4.qv<0.01))
#length(which(pc4.qv<0.001))


#pc1.FDR.01 <- colnames(snp.imp)[which(pc1.qv < 0.01)]
#pc2.FDR.01 <- colnames(snp.imp)[which(pc2.qv < 0.01)]
#pc3.FDR.01 <- colnames(snp.imp)[which(pc3.qv < 0.01)]
#loci.all_FDR01<-c(pc1.FDR.01,pc2.FDR.01,pc3.FDR.01)
#unique_loci.all_FDR01 <- unique(loci.all_FDR01)

pc1.FDR.1 <- colnames(snp.imp)[which(pc1.qv < 0.1)]
pc2.FDR.1 <- colnames(snp.imp)[which(pc2.qv < 0.1)]
#pc3.FDR.1 <- colnames(snp.imp)[which(pc3.qv < 0.1)]
#pc4.FDR.1 <- colnames(snp.imp)[which(pc4.qv < 0.1)]
loci.all_FDR1<-c(pc1.FDR.1,pc2.FDR.1)
unique_loci.all_FDR1 <- unique(loci.all_FDR1)
#221

#write.csv(unique_loci.all_FDR1, "QD_LFMM1_norun_pc_outliers.csv")


#aet.env<-subset(clim.points_8110, select=aet_1951_1980)
#aet.lfmm<-lfmm_ridge(Y=snp.imp,X= aet.env ,K=1)
#aet.pv<-lfmm_test(Y=snp.imp,X=aet.env,lfmm=aet.lfmm,calibrate="gif")
#aet.pv$gif  #values should be near 1
#hist(aet.pv$pvalue[,1])
#hist(aet.pv$calibrated.pvalue[,1])
#aet.qv<-qvalue(aet.pv$pvalue)$qvalues
#length(which(aet.qv<0.1))
#length(which(aet.qv<0.01))
#length(which(aet.qv<0.001))
#aet.FDR.001 <- colnames(snp.imp)[which(aet.qv < 0.001)]
##368 SNPs
#aet.FDR.01 <- colnames(snp.imp)[which(aet.qv < 0.01)]
##1211 SNPs
#
#cwd.env<-subset(clim.points_8110, select=cwd_1951_1980)
#cwd.lfmm<-lfmm_ridge(Y=snp.imp,X= cwd.env ,K=1)
#cwd.pv<-lfmm_test(Y=snp.imp,X=cwd.env,lfmm=cwd.lfmm,calibrate="gif")
#cwd.pv$gif  #values should be near 1
#hist(cwd.pv$pvalue[,1])
#hist(cwd.pv$calibrated.pvalue[,1])
#cwd.qv<-qvalue(cwd.pv$pvalue)$qvalues
#length(which(cwd.qv<0.1))
#length(which(cwd.qv<0.01))
#length(which(cwd.qv<0.001))
#cwd.FDR.001 <- colnames(snp.imp)[which(cwd.qv < 0.001)]
##4182 SNPs
#cwd.FDR.01 <- colnames(snp.imp)[which(cwd.qv < 0.01)]
##1688 SNPs
#
#pet.env<-subset(clim.points_8110, select=pet_1951_1980)
#pet.lfmm<-lfmm_ridge(Y=snp.imp,X= pet.env ,K=1)
#pet.pv<-lfmm_test(Y=snp.imp,X=pet.env,lfmm=pet.lfmm,calibrate="gif")
#pet.pv$gif  #values should be near 1
#hist(pet.pv$pvalue[,1])
#hist(pet.pv$calibrated.pvalue[,1])
#pet.qv<-qvalue(pet.pv$pvalue)$qvalues
#length(which(pet.qv<0.1))
#length(which(pet.qv<0.01))
#length(which(pet.qv<0.001))
#pet.FDR.001 <- colnames(snp.imp)[which(pet.qv < 0.001)]
##631 SNPs
#pet.FDR.01 <- colnames(snp.imp)[which(pet.qv < 0.01)]
##3162 SNPs
#
#ppt.env<-subset(clim.points_8110, select=ppt_1951_1980)
#ppt.lfmm<-lfmm_ridge(Y=snp.imp,X= ppt.env ,K=1)
#ppt.pv<-lfmm_test(Y=snp.imp,X=ppt.env,lfmm=ppt.lfmm,calibrate="gif")
#ppt.pv$gif  #values should be near 1
#hist(ppt.pv$pvalue[,1])
#hist(ppt.pv$calibrated.pvalue[,1])
#ppt.qv<-qvalue(ppt.pv$pvalue)$qvalues
#length(which(ppt.qv<0.1))
#length(which(ppt.qv<0.01))
#length(which(ppt.qv<0.001))
#ppt.FDR.001 <- colnames(snp.imp)[which(ppt.qv < 0.001)]
##4568 SNPs
#ppt.FDR.1 <- colnames(snp.imp)[which(ppt.qv < 0.1)]
##7870 SNPs
#
#rch.env<-subset(clim.points_8110, select=rch_1951_1980)
#rch.lfmm<-lfmm_ridge(Y=snp.imp,X= rch.env ,K=1)
#rch.pv<-lfmm_test(Y=snp.imp,X=rch.env,lfmm=rch.lfmm,calibrate="gif")
#rch.pv$gif  #values should be near 1
#hist(rch.pv$pvalue[,1])
#hist(rch.pv$calibrated.pvalue[,1])
#rch.qv<-qvalue(rch.pv$pvalue)$qvalues
#length(which(rch.qv<0.1))
#length(which(rch.qv<0.01))
#length(which(rch.qv<0.001))
#rch.FDR.001 <- colnames(snp.imp)[which(rch.qv < 0.001)]
##423 SNPs
#rch.FDR.01 <- colnames(snp.imp)[which(rch.qv < 0.01)]
##1033 SNPs
#
#run.env<-subset(clim.points_8110, select=run_1951_1980)
#run.lfmm<-lfmm_ridge(Y=snp.imp,X= run.env ,K=1)
#run.pv<-lfmm_test(Y=snp.imp,X=run.env,lfmm=run.lfmm,calibrate="gif")
#run.pv$gif  #values should be near 1
#hist(run.pv$pvalue[,1])
#hist(run.pv$calibrated.pvalue[,1])
#run.qv<-qvalue(run.pv$pvalue)$qvalues
#length(which(run.qv<0.1))
#length(which(run.qv<0.01))
#length(which(run.qv<0.001))
#run.FDR.001 <- colnames(snp.imp)[which(run.qv < 0.001)]
##118 SNPs
#run.FDR.1 <- colnames(snp.imp)[which(run.qv < 0.1)]
##488 SNPs
#
#
#str.env<-subset(clim.points_8110, select=str_1951_1980)
#str.lfmm<-lfmm_ridge(Y=snp.imp,X= str.env ,K=1)
#str.pv<-lfmm_test(Y=snp.imp,X=str.env,lfmm=str.lfmm,calibrate="gif")
#str.pv$gif  #values should be near 1
#hist(str.pv$pvalue[,1])
#hist(str.pv$calibrated.pvalue[,1])
#str.qv<-qvalue(str.pv$pvalue)$qvalues
#length(which(str.qv<0.1))
#length(which(str.qv<0.01))
#length(which(str.qv<0.001))
#str.FDR.001 <- colnames(snp.imp)[which(str.qv < 0.001)]
#str.FDR.1 <- colnames(snp.imp)[which(str.qv < 0.1)]
#
#tmn.env<-subset(clim.points_8110, select=tmn_1951_1980)
#tmn.lfmm<-lfmm_ridge(Y=snp.imp,X= tmn.env ,K=1)
#tmn.pv<-lfmm_test(Y=snp.imp,X=tmn.env,lfmm=tmn.lfmm,calibrate="gif")
#tmn.pv$gif  #values should be near 1
#hist(tmn.pv$pvalue[,1])
#hist(tmn.pv$calibrated.pvalue[,1])
#tmn.qv<-qvalue(tmn.pv$pvalue)$qvalues
#length(which(tmn.qv<0.1))
#length(which(tmn.qv<0.01))
#length(which(tmn.qv<0.001))
#tmn.FDR.001 <- colnames(snp.imp)[which(tmn.qv < 0.001)]
##3 SNPs
#tmn.FDR.01 <- colnames(snp.imp)[which(tmn.qv < 0.01)]
##41 SNPs
#
#tmx.env<-subset(clim.points_8110, select=tmx_1951_1980)
#tmx.lfmm<-lfmm_ridge(Y=snp.imp,X= tmx.env ,K=1)
#tmx.pv<-lfmm_test(Y=snp.imp,X=tmx.env,lfmm=tmx.lfmm,calibrate="gif")
#tmx.pv$gif  #values should be near 1
#hist(tmx.pv$pvalue[,1])
#hist(tmx.pv$calibrated.pvalue[,1])
#tmx.qv<-qvalue(tmx.pv$pvalue)$qvalues
#length(which(tmx.qv<0.1))
#length(which(tmx.qv<0.01))
#length(which(tmx.qv<0.001))
#tmx.FDR.001 <- colnames(snp.imp)[which(tmx.qv < 0.001)]
##80 SNPs
#tmx.FDR.1 <- colnames(snp.imp)[which(tmx.qv < 0.1)]
##453 SNPs
#
#
#loci.all_FDR1<-c(ppt.FDR.1,str.FDR.1,tmx.FDR.1)
#unique_loci.all_FDR1 <- unique(loci.all_FDR1)
#
#loci.all_FDR001<-c(aet.FDR.001,pet.FDR.001,cwd.FDR.001,ppt.FDR.001,rch.FDR.001,run.FDR.001,str.FDR.001,tmn.FDR.001,tmx.FDR.001)
#unique_loci.all_FDR001 <- unique(loci.all_FDR001)
#
#write.csv(unique_loci.all_FDR001, "QA_LFMM001_outliers.csv")


#remove duplicates between LFMM and RDA
RDAloci<-top_outliers_unconstrained_10percent$Loci
combined_LFMMRDA_loci<-c(unique_loci.all_FDR1,RDAloci)
unique_combinedloci <- unique(combined_LFMMRDA_loci)
#958

write.csv(unique_combinedloci, "QD_RoseVariables_terra_RDAfour_LFMM1pca_outliers.csv")




###### extract specific adaptive loci ######
module load vcftools/0.1.16
vcftools --gzvcf SCA10nSON1nCON4nLA3n38nSD2n72remJcalONLY_meanDP5_genoDP5_MAF0.01_missing0.9.Juglans_jhin.repeatsOut.ef.vcf.gz --positions ~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY/JcalONLY_RDA10p_RoseVariables_terra_PC1_LatxLon_outliers_con.txt --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon.vcf.gz

vcftools --gzvcf JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon.vcf.gz --012 --out JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp
cut -f2- JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.012 | sed 's/-1/NA/g' >JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.temp
tr -d '\t' <JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.012.pos | tr '\n' '\t' | sed 's/[[:space:]]*$//' >header
paste <(echo "ID" | cat - JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.012.indv) <(echo "" | cat header - JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.temp) > JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.forR
rm header JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.temp



##back in R
SUBsnp <- read.table("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/JcalONLY_RDAonly10p_LatxLonPC1_outliers_uncon_snp.forR", header = T, row.names = 1)
library(parallel)
cl<-makeCluster(detectCores())
SUBsnp.imp <- parApply(cl,SUBsnp, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x)))))) #remove NAs and replace with most common genotype

library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
JCenv<-data.frame(JCenv)
Coordinates<-subset(JCenv, select=c(x,y))
Coordinates<-data.frame(Coordinates)


library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA
remove_NAs_stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove_NAs_stack(ras_8110)

sample.coord.sp <- vect(Coordinates, geom = c("x", "y"), crs = crs(ras_8110))
clim.points_8110 <- extract(ras_8110, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
clim.points_8110<-data.frame(clim.points_8110)


library(gradientForest)
env.gf <- cbind(clim.points_8110[,-1])
maxLevel <- log2(0.368*nrow(env.gf)/2)
gf <- gradientForest(cbind(env.gf, SUBsnp.imp), predictor.vars=colnames(env.gf), response.vars=colnames(SUBsnp.imp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)



library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp45_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C45_2040 <- rast(tif_files)
layer_names <- gsub("CNRMrcp45_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C45_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C45_2040 <- project(ras_C45_2040, new_crs, method = "bilinear")
ras_C45_2040[ras_C45_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp45_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C45_2070 <- rast(tif_files)
layer_names <- gsub("CNRMrcp45_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C45_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C45_2070 <- project(ras_C45_2070, new_crs, method = "bilinear")
ras_C45_2070[ras_C45_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp85_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C85_2040 <- rast(tif_files)
layer_names <- gsub("CNRMrcp85_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C85_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C85_2040 <- project(ras_C85_2040, new_crs, method = "bilinear")
ras_C85_2040[ras_C85_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/cnrm_rcp85_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_C85_2070 <- rast(tif_files)
layer_names <- gsub("CNRMrcp85_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_C85_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_C85_2070 <- project(ras_C85_2070, new_crs, method = "bilinear")
ras_C85_2070[ras_C85_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp45_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H45_2040 <- rast(tif_files)
layer_names <- gsub("HADES_rcp45_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H45_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H45_2040 <- project(ras_H45_2040, new_crs, method = "bilinear")
ras_H45_2040[ras_H45_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp45_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H45_2070 <- rast(tif_files)
layer_names <- gsub("HADES_rcp45_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H45_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H45_2070 <- project(ras_H45_2070, new_crs, method = "bilinear")
ras_H45_2070[ras_H45_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp85_2040_2069', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H85_2040 <- rast(tif_files)
layer_names <- gsub("HADES_rcp85_2040_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H85_2040) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H85_2040 <- project(ras_H85_2040, new_crs, method = "bilinear")
ras_H85_2040[ras_H85_2040 == -9999] <- NA

tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/hades_rcp85_2070_2099', 
                        pattern = "\\.tif$", full.names = TRUE)
ras_H85_2070 <- rast(tif_files)
layer_names <- gsub("HADES_rcp85_2070_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras_H85_2070) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_H85_2070 <- project(ras_H85_2070, new_crs, method = "bilinear")
ras_H85_2070[ras_H85_2070 == -9999] <- NA


tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA

remove.NAs.stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove.NAs.stack(ras_8110)
ras_H45_2040 <- remove.NAs.stack(ras_H45_2040)
ras_H45_2070 <- remove.NAs.stack(ras_H45_2070)
ras_H85_2040 <- remove.NAs.stack(ras_H85_2040)
ras_H85_2070 <- remove.NAs.stack(ras_H85_2070)
ras_C45_2040 <- remove.NAs.stack(ras_C45_2040)
ras_C45_2070 <- remove.NAs.stack(ras_C45_2070)
ras_C85_2040 <- remove.NAs.stack(ras_C85_2040)
ras_C85_2070 <- remove.NAs.stack(ras_C85_2070)



clim.land_8110 <- as.data.frame(ras_8110, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_8110[clim.land_8110 == -9999] <- NA
clim.land_8110 <- na.omit(clim.land_8110)
clim.land_H45_2040 <- as.data.frame(ras_H45_2040, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_H45_2040[clim.land_H45_2040 == -9999] <- NA
clim.land_H45_2040 <- na.omit(clim.land_H45_2040)
clim.land_H45_2070 <- as.data.frame(ras_H45_2070, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_H45_2070[clim.land_H45_2070 == -9999] <- NA
clim.land_H45_2070 <- na.omit(clim.land_H45_2070)
clim.land_H85_2040 <- as.data.frame(ras_H85_2040, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_H85_2040[clim.land_H85_2040 == -9999] <- NA
clim.land_H85_2040 <- na.omit(clim.land_H85_2040)
clim.land_H85_2070 <- as.data.frame(ras_H85_2070, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_H85_2070[clim.land_H85_2070 == -9999] <- NA
clim.land_H85_2070 <- na.omit(clim.land_H85_2070)
clim.land_C45_2040 <- as.data.frame(ras_C45_2040, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_C45_2040[clim.land_C45_2040 == -9999] <- NA
clim.land_C45_2040 <- na.omit(clim.land_C45_2040)
clim.land_C45_2070 <- as.data.frame(ras_C45_2070, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_C45_2070[clim.land_C45_2070 == -9999] <- NA
clim.land_C45_2070 <- na.omit(clim.land_C45_2070)
clim.land_C85_2040 <- as.data.frame(ras_C85_2040, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_C85_2040[clim.land_C85_2040 == -9999] <- NA
clim.land_C85_2040 <- na.omit(clim.land_C85_2040)
clim.land_C85_2070 <- as.data.frame(ras_C85_2070, cells = TRUE, xy = TRUE, na.rm = FALSE)
clim.land_C85_2070[clim.land_C85_2070 == -9999] <- NA
clim.land_C85_2070 <- na.omit(clim.land_C85_2070)

clim.land_H45_2040 <- clim.land_H45_2040[clim.land_H45_2040$cell %in% clim.land_8110$cell, ]
clim.land_H45_2070 <- clim.land_H45_2070[clim.land_H45_2070$cell %in% clim.land_8110$cell, ]
clim.land_H85_2040 <- clim.land_H85_2040[clim.land_H85_2040$cell %in% clim.land_8110$cell, ]
clim.land_H85_2070 <- clim.land_H85_2070[clim.land_H85_2070$cell %in% clim.land_8110$cell, ]
clim.land_C45_2040 <- clim.land_C45_2040[clim.land_C45_2040$cell %in% clim.land_8110$cell, ]
clim.land_C45_2070 <- clim.land_C45_2070[clim.land_C45_2070$cell %in% clim.land_8110$cell, ]
clim.land_C85_2040 <- clim.land_C85_2040[clim.land_C85_2040$cell %in% clim.land_8110$cell, ]
clim.land_C85_2070 <- clim.land_C85_2070[clim.land_C85_2070$cell %in% clim.land_8110$cell, ]
clim.land_8110 <- clim.land_8110[clim.land_8110$cell %in% clim.land_C85_2070$cell, ]







clim.land_8110 <- subset(clim.land_8110, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_H45_2040 <- subset(clim.land_H45_2040, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_H45_2070 <- subset(clim.land_H45_2070, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_H85_2040 <- subset(clim.land_H85_2040, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_H85_2070 <- subset(clim.land_H85_2070, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_C45_2040 <- subset(clim.land_C45_2040, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_C45_2070 <- subset(clim.land_C45_2070, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_C85_2040 <- subset(clim.land_C85_2040, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))
clim.land_C85_2070 <- subset(clim.land_C85_2070, select=c(cell,aet,awc,cwd,depth,pct_clay,ph,ppt_djf,ppt_jja,terrain,tmn))




colnames(clim.land_8110)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_H45_2040)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_H45_2070)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_H85_2040)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_H85_2070)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_C45_2040)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_C45_2070)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_C85_2040)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")
colnames(clim.land_C85_2070)<-c("ID","aet","awc","cwd","depth","pct_clay","ph","ppt_djf","ppt_jja","terrain","tmn")


library(gradientForest)
pred_8110 <- predict(gf, clim.land_8110[,-1])  #note the removal of the cell ID column with [,-1])
pred_H45_2040 <- predict(gf, clim.land_H45_2040[,-1])  #note the removal of the cell ID column with [,-1])
pred_H45_2070 <- predict(gf, clim.land_H45_2070[,-1])  #note the removal of the cell ID column with [,-1])
pred_H85_2040 <- predict(gf, clim.land_H85_2040[,-1])  #note the removal of the cell ID column with [,-1])
pred_H85_2070 <- predict(gf, clim.land_H85_2070[,-1])  #note the removal of the cell ID column with [,-1])
pred_C45_2040 <- predict(gf, clim.land_C45_2040[,-1])  #note the removal of the cell ID column with [,-1])
pred_C45_2070 <- predict(gf, clim.land_C45_2070[,-1])  #note the removal of the cell ID column with [,-1])
pred_C85_2040 <- predict(gf, clim.land_C85_2040[,-1])  #note the removal of the cell ID column with [,-1])
pred_C85_2070 <- predict(gf, clim.land_C85_2070[,-1])  #note the removal of the cell ID column with [,-1])


pcaToRaster <- function(snpPreds, rast, mapCells){
   
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
    
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  r_rast <- g_rast <- b_rast <- rast[[1]] * NA
  values(r_rast)[mapCells] <- scalR
  values(g_rast)[mapCells] <- scalG
  values(b_rast)[mapCells] <- scalB
  ##stacks color rasters
  outRast <- c(r_rast, g_rast, b_rast)
  names(outRast) <- c("red", "green", "blue")
  return(outRast)
}


#genomic turnover
RGBmap_8110 <- pcaToRaster(pred_8110, ras_8110$aet, clim.land_8110$ID)
#pdf("8110gf.pdf")
#plotRGB(RGBmap_8110)
#dev.off()
writeRaster(RGBmap_8110, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFturnover_8110.tif", filetype = "GTiff", overwrite = TRUE)

#local genomic offset

genOffset_8110_H452040 <- sqrt((pred_H45_2040[,1]-pred_8110[,1])^2+(pred_H45_2040[,2]-pred_8110[,2])^2
                    +(pred_H45_2040[,3]-pred_8110[,3])^2+(pred_H45_2040[,4]-pred_8110[,4])^2
					+(pred_H45_2040[,5]-pred_8110[,5])^2+(pred_H45_2040[,6]-pred_8110[,6])^2
					+(pred_H45_2040[,7]-pred_8110[,7])^2+(pred_H45_2040[,8]-pred_8110[,8])^2
					+(pred_H45_2040[,9]-pred_8110[,9])^2+(pred_H45_2040[,10]-pred_8110[,10])^2
                    )
mask<-ras_8110$aet
values(mask) <- NA
mask[clim.land_H45_2040$ID] <- genOffset_8110_H452040
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_H452040.tif", filetype = "GTiff", overwrite=TRUE)

genOffset_8110_H452070 <- sqrt((pred_H45_2070[,1]-pred_8110[,1])^2+(pred_H45_2070[,2]-pred_8110[,2])^2
                    +(pred_H45_2070[,3]-pred_8110[,3])^2+(pred_H45_2070[,4]-pred_8110[,4])^2
					+(pred_H45_2070[,5]-pred_8110[,5])^2+(pred_H45_2070[,6]-pred_8110[,6])^2
					+(pred_H45_2070[,7]-pred_8110[,7])^2+(pred_H45_2070[,8]-pred_8110[,8])^2
					+(pred_H45_2070[,9]-pred_8110[,9])^2+(pred_H45_2070[,10]-pred_8110[,10])^2
                    )
mask<-ras_8110$aet
values(mask) <- NA
mask[clim.land_H45_2070$ID] <- genOffset_8110_H452070
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_H452070.tif", filetype = "GTiff", overwrite=TRUE)



genOffset_8110_H852040 <- sqrt((pred_H85_2040[,1]-pred_8110[,1])^2+(pred_H85_2040[,2]-pred_8110[,2])^2
                    +(pred_H85_2040[,3]-pred_8110[,3])^2+(pred_H85_2040[,4]-pred_8110[,4])^2
					+(pred_H85_2040[,5]-pred_8110[,5])^2+(pred_H85_2040[,6]-pred_8110[,6])^2
					+(pred_H85_2040[,7]-pred_8110[,7])^2+(pred_H85_2040[,8]-pred_8110[,8])^2
					+(pred_H85_2040[,9]-pred_8110[,9])^2+(pred_H85_2040[,10]-pred_8110[,10])^2
                    )
mask<-ras_8110$aet
values(mask) <- NA
mask[clim.land_H85_2040$ID] <- genOffset_8110_H852040
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_H852040.tif", filetype = "GTiff", overwrite=TRUE)

genOffset_8110_H852070 <- sqrt((pred_H85_2070[,1]-pred_8110[,1])^2+(pred_H85_2070[,2]-pred_8110[,2])^2
                    +(pred_H85_2070[,3]-pred_8110[,3])^2+(pred_H85_2070[,4]-pred_8110[,4])^2
					+(pred_H85_2070[,5]-pred_8110[,5])^2+(pred_H85_2070[,6]-pred_8110[,6])^2
					+(pred_H85_2070[,7]-pred_8110[,7])^2+(pred_H85_2070[,8]-pred_8110[,8])^2
					+(pred_H85_2070[,9]-pred_8110[,9])^2+(pred_H85_2070[,10]-pred_8110[,10])^2
                    )
mask <- ras_8110$aet
mask[] <- NA  # Set initial raster values to NA
mask[clim.land_H85_2070$ID] <- genOffset_8110_H852070
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_H852070.tif", filetype = "GTiff", overwrite=TRUE)


genOffset_8110_C452040 <- sqrt((pred_C45_2040[,1]-pred_8110[,1])^2+(pred_C45_2040[,2]-pred_8110[,2])^2
                    +(pred_C45_2040[,3]-pred_8110[,3])^2+(pred_C45_2040[,4]-pred_8110[,4])^2
					+(pred_C45_2040[,5]-pred_8110[,5])^2+(pred_C45_2040[,6]-pred_8110[,6])^2
					+(pred_C45_2040[,7]-pred_8110[,7])^2+(pred_C45_2040[,8]-pred_8110[,8])^2
					+(pred_C45_2040[,9]-pred_8110[,9])^2+(pred_C45_2040[,10]-pred_8110[,10])^2
                    )
mask<-ras_8110$aet
values(mask) <- NA
mask[clim.land_C45_2040$ID] <- genOffset_8110_C452040
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_C452040.tif", filetype = "GTiff", overwrite=TRUE)

genOffset_8110_C452070 <- sqrt((pred_C45_2070[,1]-pred_8110[,1])^2+(pred_C45_2070[,2]-pred_8110[,2])^2
                    +(pred_C45_2070[,3]-pred_8110[,3])^2+(pred_C45_2070[,4]-pred_8110[,4])^2
					+(pred_C45_2070[,5]-pred_8110[,5])^2+(pred_C45_2070[,6]-pred_8110[,6])^2
					+(pred_C45_2070[,7]-pred_8110[,7])^2+(pred_C45_2070[,8]-pred_8110[,8])^2
					+(pred_C45_2070[,9]-pred_8110[,9])^2+(pred_C45_2070[,10]-pred_8110[,10])^2
                    )
mask<-ras_8110$aet
values(mask) <- NA
mask[clim.land_C45_2070$ID] <- genOffset_8110_C452070
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_C452070.tif", filetype = "GTiff", overwrite=TRUE)



genOffset_8110_C852040 <- sqrt((pred_C85_2040[,1]-pred_8110[,1])^2+(pred_C85_2040[,2]-pred_8110[,2])^2
                    +(pred_C85_2040[,3]-pred_8110[,3])^2+(pred_C85_2040[,4]-pred_8110[,4])^2
					+(pred_C85_2040[,5]-pred_8110[,5])^2+(pred_C85_2040[,6]-pred_8110[,6])^2
					+(pred_C85_2040[,7]-pred_8110[,7])^2+(pred_C85_2040[,8]-pred_8110[,8])^2
					+(pred_C85_2040[,9]-pred_8110[,9])^2+(pred_C85_2040[,10]-pred_8110[,10])^2
                    )
mask<-ras_8110$aet
values(mask) <- NA
mask[clim.land_C85_2040$ID] <- genOffset_8110_C852040
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_C852040.tif", filetype = "GTiff", overwrite=TRUE)

genOffset_8110_C852070 <- sqrt((pred_C85_2070[,1]-pred_8110[,1])^2+(pred_C85_2070[,2]-pred_8110[,2])^2
                    +(pred_C85_2070[,3]-pred_8110[,3])^2+(pred_C85_2070[,4]-pred_8110[,4])^2
					+(pred_C85_2070[,5]-pred_8110[,5])^2+(pred_C85_2070[,6]-pred_8110[,6])^2
					+(pred_C85_2070[,7]-pred_8110[,7])^2+(pred_C85_2070[,8]-pred_8110[,8])^2
					+(pred_C85_2070[,9]-pred_8110[,9])^2+(pred_C85_2070[,10]-pred_8110[,10])^2
                    )
mask <- ras_8110$aet
mask[] <- NA  # Set initial raster values to NA
mask[clim.land_C85_2070$ID] <- genOffset_8110_C852070
writeRaster(mask, "JcalONLY_RDAonly10p_LatxLonPC1_outliers_con_GFgenOffset_8110_C852070.tif", filetype = "GTiff", overwrite=TRUE)





############### wingen ######################
do_everything_for_me(liz_vcf, liz_coords, CA_env)
library(algatr)
library(wingen)
library(raster)
library(terra)
library(vcfR)
library(adegenet)
library(parallel)
library(Rcpp)

#++++++++++++++++++++++++++++++++++++
#HOW TO SAVE RASTERS FOR NEXT SESSION
s <- lapply(wgd_het, wrap)
saveRDS(s, 'wgd_het.rds')

wgd_het <- readRDS("wgd_het.rds")
wgd_het <- lapply(wgd_het, rast)
#++++++++++++++++++++++++++++++++++++



options(future.globals.maxSize = 8000 * 1024^2)

vcf<-read.vcfR('~/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/manyrem_Qdoufiltered.vcf')

library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/TNC_Qdou_coordsforGF_manyrem.csv")
Coordinates<-subset(JCenv, select=c(x,y))
Coordinates<-data.frame(Coordinates)

library(readxl)
Obs<-read_excel("/u/home/r/rcbuck/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/QDcalflora_LocationsforRaster.xlsx")
Observations<-subset(Obs,select=c(Longitude,Latitude))
colnames(Observations)<- c("x","y")
Observations<-as.data.frame(Observations)
options(future.globals.maxSize = 8000 * 1024^2)
qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)


preview_gd(qdou_lyr,Coordinates,wdim = 15,fact = 10)
wgd_het <- window_gd(vcf,Coordinates,qdou_lyr,stat = "Ho",wdim = 15,fact = 10,rarify=T)
#save SpatRast for next time
s <- lapply(wgd_het, wrap)
saveRDS(s, 'wgd_het.rds')
#load back in during next session
#wgd_het <- readRDS("wgd_het.rds")
#wgd_het <- lapply(wgd_het, rast)



load_algatr_example()
Cali_envlayer <- rast(CA_env$CA_rPCA1)


par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot map of pi
#pdf("wgd_het.pdf")
plot_gd(wgd_het, main = "Moving window Heterozygosity", bkg = Cali_envlayer)
#dev.off()

# Plot sample count map
pdf("wgd_het_samplecounts.pdf")
plot_count(wgd_het, main = "Sample counts")
dev.off()


kgd <- krig_gd(wgd_het, index = 1, grd= qdou_lyr)
summary(kgd)
#save SpatRast for next time
k <- lapply(kgd, wrap)
saveRDS(k, 'kgd.rds')
#load back in during next session
#kgd <- readRDS("kgd.rds")
#kgd <- lapply(kgd, rast)


par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot kriged map of pi
#pdf("Kriged_wgd_het.pdf")
plot_gd(kgd, main = "Kriged Heterozygosity")
#dev.off()

# Plot kriged sample count map
pdf("Kriged_het_samplecounts.pdf")
plot_count(kgd, main = "Kriged sample counts")
dev.off()

kgd_het_ras<-raster(kgd)
writeRaster(kgd_het_ras,"QD_krigedHo_biggerwindow.tif","GTiff",overwrite=T)

#do pi (nucleotide diversity)
wgd_pi <- window_gd(vcf,Coordinates,qdou_lyr,stat = "pi",wdim = 15,fact = 10,rarify=T)
kgd_pi <- krig_gd(wgd_pi, index = 1, grd= qdou_lyr)
kgd_pi_ras<-raster(kgd_pi)
writeRaster(kgd_pi_ras,"QD_krigedPi_biggerwindow.tif","GTiff",overwrite=T)



#mgd_1 <- mask_gd(kgd, kgd[["sample_count"]], minval = 1)
#mgd_2 <- mask_gd(kgd, kgd[["sample_count"]], minval = 2)
#
#par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
#pdf("Masked_het_minval1.pdf")
#plot_gd(mgd_1, main = "Kriged & masked het minval 1", bkg = Cali_envlayer)
#dev.off()
#
#pdf("Masked_het_minval2")
#plot_gd(mgd_2, main = "Kriged & masked het minval2", bkg = Cali_envlayer)
#dev.off()
#
## Resample qagr_lyr based on masked layer
#r <- terra::resample(Cali_envlayer, mgd_1)
## Perform masking
#mgd <- mask_gd(mgd_1, r)
## Plot masked map
#pdf("Cali_masked_het.pdf")
#plot_gd(mgd, main = "Kriged & masked pi", bkg = r)
#dev.off()

#$$$$$$$$$$$$$$$$$$$$$$$$$
kgd_pi_ras<-raster(kgd_pi)
writeRaster(kgd_pi_ras,"krigedPi_qlob.tif","GTiff",overwrite=T)
#$$$$$$$$$$$$$$$$$$$$$$$$$




library(algatr)
library(wingen)
library(raster)
library(terra)
library(vcfR)
library(adegenet)
library(parallel)
library(Rcpp)
library(readr)

####DO BY PRESERVE ###
library(readxl)
Obs<-read_excel("~/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/QDcalflora_LocationsforRaster.xlsx")
Observations<-subset(Obs,select=c(Longitude,Latitude))
colnames(Observations)<- c("x","y")
Observations<-as.data.frame(Observations)
options(future.globals.maxSize = 8000 * 1024^2)
qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)

load_algatr_example()
Cali_envlayer <- rast(CA_env$CA_rPCA1)


qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)
#DyeCreek
vcf_DyeCreek<-read.vcfR('~/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/DyeCreekonly_LDfiltered.vcf')
Denv <- read_csv("/u/home/r/rcbuck/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/TNC_Qdou_coordsDYECREEKonly.csv")
geoDyeCreek<-subset(Denv, select=c(x,y))
colnames(geoDyeCreek)<-c("Longitude","Latitude")
geoDyeCreek<-data.frame(geoDyeCreek)
#crop obs to DyeCreek
ext(qdou_lyr)<-ext(-122.14,-121.85,40.05,40.178)
preview_gd(qdou_lyr,geoDyeCreek,wdim = 15,fact = 10)
wgd_het_Dye <- window_gd(vcf_DyeCreek,geoDyeCreek,qdou_lyr,stat = "Ho",wdim = 10,fact = 15,rarify=T)
#save SpatRast for next time
sDye <- lapply(wgd_het_Dye, wrap)
saveRDS(sDye, 'wgd_het_Dye.rds')
#load back in during next session
#wgd_het_Dye <- readRDS("wgd_het_Dye.rds")
#wgd_het_Dye <- lapply(wgd_het_Dye, rast)
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot map of pi
plot_gd(wgd_het_Dye, main = "Moving window Heterozygosity", bkg = qdou_lyr)
# Plot sample count map
plot_count(wgd_het_Dye, main = "Sample counts")

kgd_Dye <- krig_gd(wgd_het_Dye, index = 1, grd= qdou_lyr)
summary(kgd_Dye)
#save SpatRast for next time
kDye <- lapply(kgd_Dye, wrap)
saveRDS(kDye, 'kgd_Dye.rds')
#load back in during next session
#kgd <- readRDS("kgd.rds")
#kgd <- lapply(kgd, rast)
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot kriged map of pi
plot_gd(kgd_Dye, main = "Kriged Heterozygosity")
# Plot kriged sample count map
plot_count(kgd_Dye, main = "Kriged sample counts")
kgd_het_ras_Dye<-raster(kgd_Dye)
writeRaster(kgd_het_ras_Dye,"krigedHo_DyeCreek_biggerwindow.tif","GTiff",overwrite=T)

#pi (nucleotide diversity)
wgd_pi_Dye <- window_gd(vcf_DyeCreek,geoDyeCreek,qdou_lyr,stat = "pi",wdim = 10,fact = 15,rarify=T)
kgd_pi_Dye <- krig_gd(wgd_pi_Dye, index = 1, grd= qdou_lyr)
kgd_pi_ras_Dye<-raster(kgd_pi_Dye)
writeRaster(kgd_pi_ras_Dye,"krigedPi_DyeCreek_biggerwindow.tif","GTiff",overwrite=T)



qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)
#Randall
vcf_Randall<-read.vcfR('~/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/Randallonly_LDfiltered.vcf')
Renv <- read_csv("/u/home/r/rcbuck/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/TNC_Qdou_coordsRANDALLonly.csv")
geoRandall<-subset(Renv, select=c(x,y))
colnames(geoRandall)<-c("Longitude","Latitude")
geoRandall<-data.frame(geoRandall)
#crop obs to Randall
qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)
ext(qdou_lyr)<-ext(-118.72,-118.37,35.13,35.42)
preview_gd(qdou_lyr,geoRandall,wdim = 15,fact = 10)
wgd_het_Rand <- window_gd(vcf_Randall,geoRandall,qdou_lyr,stat = "Ho",wdim = 15,fact = 10,rarify=T)
#save SpatRast for next time
sRand <- lapply(wgd_het_Rand, wrap)
saveRDS(sRand, 'wgd_het_Rand.rds')
#load back in during next session
#wgd_het_Dang <- readRDS("wgd_het_Dang.rds")
#wgd_het_Dang <- lapply(wgd_het_Dang, rast)
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot map of pi
plot_gd(wgd_het_Rand, main = "Moving window Heterozygosity", bkg = qdou_lyr)
# Plot sample count map
plot_count(wgd_het_Rand, main = "Sample counts")

kgd_Rand <- krig_gd(wgd_het_Rand, index = 1, grd= qdou_lyr)
summary(kgd)
#save SpatRast for next time
kRand <- lapply(kgd_Rand, wrap)
saveRDS(kRand, 'kgd_Rand.rds')
#load back in during next session
#kgd_Rand <- readRDS("kgd_Rand.rds")
#kgd_Rand <- lapply(kgd_Rand, rast)
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot kriged map of pi
plot_gd(kgd_Rand, main = "Kriged Heterozygosity")
# Plot kriged sample count map
plot_count(kgd_Rand, main = "Kriged sample counts")
kgd_het_ras_Rand<-raster(kgd_Rand)
writeRaster(kgd_het_ras_Rand,"krigedHo_Randall_biggerwindow.tif","GTiff",overwrite=T)

#pi (nucleotide diversity)
wgd_pi_Rand <- window_gd(vcf_Randall,geoRandall,qdou_lyr,stat = "pi",wdim = 15,fact = 10,rarify=T)
kgd_pi_Rand <- krig_gd(wgd_pi_Rand, index = 1, grd= qdou_lyr)
kgd_pi_ras_Rand<-raster(kgd_pi_Rand)
writeRaster(kgd_pi_ras_Rand,"krigedPi_Randall_biggerwindow.tif","GTiff",overwrite=T)



qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)
#Santa Rosa
vcf_LasPiletas<-read.vcfR('~/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/LasPiletasonly_LDfiltered.vcf')
LPenv <- read_csv("/u/home/r/rcbuck/project-vlsork/TNC_Qdouglasii/MarkedDuplicates/vcfs_bychr_REDO/TNC_Qdou_coordsLASPILETASonly.csv")
geoLasPiletas<-subset(LPenv, select=c(x,y))
colnames(geoLasPiletas)<-c("Longitude","Latitude")
geoLasPiletas<-data.frame(geoLasPiletas)
#crop obs to Santa Rosa
qdou_lyr <- coords_to_raster(Observations, res = 0.01, buffer = 9)
ext(qdou_lyr)<-ext(-120.17,-120.02,35.23,35.339)
preview_gd(qdou_lyr,geoLasPiletas,wdim = 15,fact = 10)
wgd_het_LP <- window_gd(vcf_LasPiletas,geoLasPiletas,qdou_lyr,stat = "Ho",wdim = 10,fact = 15,rarify=T)
#save SpatRast for next time
sLP <- lapply(wgd_het_LP, wrap)
saveRDS(sLP, 'wgd_het_LP.rds')
#load back in during next session
#wgd_het_LP <- readRDS("wgd_het_LP.rds")
#wgd_het_LP <- lapply(wgd_het_LP, rast)
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot map of pi
plot_gd(wgd_het_LP, main = "Moving window Heterozygosity", bkg = qdou_lyr)
# Plot sample count map
plot_count(wgd_het_LP, main = "Sample counts")

kgd_LP <- krig_gd(wgd_het_LP, index = 1, grd= qdou_lyr)
summary(kgd_LP)
#save SpatRast for next time
kLP <- lapply(kgd_LP, wrap)
saveRDS(kLP, 'kgd_LP.rds')
#load back in during next session
#kgd_LP <- readRDS("kgd_LP.rds")
#kgd_LP <- lapply(kgd_LP, rast)
par(mfrow = c(1, 2), oma = rep(1, 4), mar = rep(2, 4))
# Plot kriged map of pi
plot_gd(kgd_LP, main = "Kriged Heterozygosity")
# Plot kriged sample count map
plot_count(kgd_LP, main = "Kriged sample counts")
kgd_het_ras_LP<-raster(kgd_LP)
writeRaster(kgd_het_ras_LP,"krigedHo_LasPiletas_biggerwindow.tif","GTiff",overwrite=T)

#pi (nucleotide diversity)
wgd_pi_LP <- window_gd(vcf_LasPiletas,geoLasPiletas,qdou_lyr,stat = "pi",wdim = 10,fact = 15,rarify=T)
kgd_pi_LP <- krig_gd(wgd_pi_LP, index = 1, grd= qdou_lyr)
kgd_pi_ras_LP<-raster(kgd_pi_LP)
writeRaster(kgd_pi_ras_LP,"krigedPi_LasPiletas_biggerwindow.tif","GTiff",overwrite=T)







####Map differences in climate variables (current to future)
diff_5180_2070_C85_aet<-ras_M45_5180$aet_1951_1980-ras_C85_2070$aet_2070_2099
writeRaster(diff_5180_2070_C85_aet,"diff_5180_2070_C85_aet.tif","GTiff",overwrite=T)
diff_5180_2070_C85_pet<-ras_M45_5180$pet_1951_1980-ras_C85_2070$pet_2070_2099
writeRaster(diff_5180_2070_C85_pet,"diff_5180_2070_C85_pet.tif","GTiff",overwrite=T)
diff_5180_2070_C85_ppt<-ras_M45_5180$ppt_1951_1980-ras_C85_2070$ppt_2070_2099
writeRaster(diff_5180_2070_C85_ppt,"diff_5180_2070_C85_ppt.tif","GTiff",overwrite=T)
diff_5180_2070_C85_tmx<-ras_M45_5180$tmx_1951_1980-ras_C85_2070$tmx_2070_2099
writeRaster(diff_5180_2070_C85_tmx,"diff_5180_2070_C85_tmx.tif","GTiff",overwrite=T)
diff_5180_2070_C85_tmn<-ras_M45_5180$tmn_1951_1980-ras_C85_2070$tmn_2070_2099
writeRaster(diff_5180_2070_C85_tmn,"diff_5180_2070_C85_tmn.tif","GTiff",overwrite=T)
diff_5180_2070_C85_run<-ras_M45_5180$run_1951_1980-ras_C85_2070$run_2070_2099
writeRaster(diff_5180_2070_C85_run,"diff_5180_2070_C85_run.tif","GTiff",overwrite=T)
diff_5180_2070_C85_rch<-ras_M45_5180$rch_1951_1980-ras_C85_2070$rch_2070_2099
writeRaster(diff_5180_2070_C85_rch,"diff_5180_2070_C85_rch.tif","GTiff",overwrite=T)
diff_5180_2070_C85_str<-ras_M45_5180$str_1951_1980-ras_C85_2070$str_2070_2099
writeRaster(diff_5180_2070_C85_str,"diff_5180_2070_C85_str.tif","GTiff",overwrite=T)
diff_5180_2070_C85_cwd<-ras_M45_5180$cwd_1951_1980-ras_C85_2070$cwd_2070_2099
writeRaster(diff_5180_2070_C85_cwd,"diff_5180_2070_C85_cwd.tif","GTiff",overwrite=T)



#isolation by distance
#isolation by distance
#for individuals
library(vcfR)
library(adegenet)
library(geosphere)
library(vegan)
library(dartR)


library(dartR)
library(readr)
coords <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
colnames(coords) <- c("sample","lon", "lat")
#coords<-subset(coords,select=c(lon,lat))
coords<-data.frame(coords)

# Load your data as a genlight object
gl_data <- gl.read.vcf('/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/SCA10nSON1nCON4nLA3n38nSD2n72rem_JcalONLY_filtered_renamed.vcf')
gl_data@other$latlon <- coords

gen_dist <- dist(as.matrix(gl_data))
coords <- coords[match(indNames(gl_data), coords$sample), ]

geo_matrix <- distm(coords[, c("lon", "lat")], fun = distHaversine) / 1000
geo_dist <- as.dist(geo_matrix)

ibd_result <- mantel(gen_dist, geo_dist, method = "pearson", permutations = 999)
print(ibd_result)




##plotting
# Flatten the distance matrices
Dgeo_vec <- as.vector(geo_dist)
Dgen_vec <- as.vector(gen_dist)

# Fit the regression
# Fit linear model
lm_fit <- lm(Dgen_vec ~ Dgeo_vec)

# Create new x values for smooth line
x_new <- seq(min(Dgeo_vec), max(Dgeo_vec), length.out = 200)

# Get predicted values + confidence interval
pred <- predict(lm_fit, newdata = data.frame(Dgeo_vec = x_new), interval = "confidence", level = 0.95)

# Plot points (black)
plot(Dgeo_vec, Dgen_vec,
     xlab = "Geographic Distance (km)",
     ylab = "Genetic Differentiation (Euclidean)",
     main = "Isolation by Distance (IBD) among J. californica individuals",
     pch = 16, col = "black")

# Add confidence interval as a shaded polygon (grey)
polygon(c(x_new, rev(x_new)), c(pred[, "upr"], rev(pred[, "lwr"])),
        col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)

# Add regression line (blue)
lines(x_new, pred[, "fit"], col = "blue", lwd = 2)

# Add equation in bottom-right
slope <- formatC(coef(lm_fit)[2], format = "e", digits = 2)
intercept <- formatC(coef(lm_fit)[1], format = "e", digits = 2)
R2 <- round(summary(lm_fit)$r.squared, 3)
pval_mantel <- 0.001 #set from mantel test
eq <- bquote(italic(y) == .(intercept) + .(slope) %.% italic(x) * ", " ~ R^2 == .(R2) * ", " ~ italic(p) == .(pval_mantel))
usr <- par("usr")
xpos <- usr[2] - 0.05*(usr[2]-usr[1])
ypos <- usr[3] + 0.05*(usr[4]-usr[3])
text(xpos, ypos, labels = eq, adj = c(1,0), cex=0.9)





library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
JCenv<-data.frame(JCenv)
Coordinates<-subset(JCenv, select=c(x,y))
Coordinates<-data.frame(Coordinates)

library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA
remove_NAs_stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove_NAs_stack(ras_8110)

sample.coord.sp <- vect(Coordinates, geom = c("x", "y"), crs = crs(ras_8110))
clim.points_8110 <- extract(ras_8110, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
clim.points_8110<-data.frame(clim.points_8110)

clim.points_8110$sample_id <- JCenv$sample
subEnv<-clim.points_8110
row.names(subEnv) <- subEnv$sample_id
subEnv <- subEnv[match(indNames(gl_data), subEnv$sample_id), ]
env_scaled <- scale(subEnv[, sapply(subEnv, is.numeric)])
env_dist <- dist(env_scaled)

mantel_env <- mantel.partial(gen_dist, env_dist, geo_dist,
                             method = "pearson",
                             permutations = 9999)
mantel_env



# ---- Flatten matrices ----
# Flatten the distance matrices to vectors
Denv_vec <- as.vector(env_dist)      # environmental distance (Euclidean, scaled)
Dgen_vec <- as.vector(gen_dist)      # genetic distance (Euclidean)

# ---- Fit linear model (genetic ~ environmental distance) ----
lm_fit_env <- lm(Dgen_vec ~ Denv_vec)

# ---- Create x values for smooth regression line ----
x_new <- seq(min(Denv_vec), max(Denv_vec), length.out = 200)

# ---- Predicted values + confidence interval ----
pred_env <- predict(lm_fit_env, newdata = data.frame(Denv_vec = x_new),
                    interval = "confidence", level = 0.95)

# ---- Plot points ----
plot(Denv_vec, Dgen_vec,
	 xlab = "Environmental Distance (Euclidean, scaled)",
     ylab = "Genetic Distance (Euclidean)",
     main = "Isolation by Environment (IBE) among J. californica Individuals",
     pch = 16, col = "black")

# Add confidence interval (grey shaded polygon)
polygon(c(x_new, rev(x_new)), c(pred_env[, "upr"], rev(pred_env[, "lwr"])),
        col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)

# Add regression line (blue)
lines(x_new, pred_env[, "fit"], col = "blue", lwd = 2)

# Add regression stats to plot
slope <- formatC(coef(lm_fit_env)[2], format = "e", digits = 2)
intercept <- formatC(coef(lm_fit_env)[1], format = "e", digits = 2)
R2 <- round(summary(lm_fit_env)$r.squared, 3)
pval_mantel <- signif(mantel_env$signif, 2)  # from partial Mantel test

eq <- bquote(italic(y) == .(intercept) + .(slope) %.% italic(x) *
             ", " ~ R^2 == .(R2) * ", " ~ italic(p) == .(pval_mantel))

usr <- par("usr")
xpos <- usr[2] - 0.05 * (usr[2] - usr[1])
ypos <- usr[3] + 0.05 * (usr[4] - usr[3])
text(xpos, ypos, labels = eq, adj = c(1, 0), cex = 0.9)






######### EXTRACT SNPS AND FUNCTION POSITIONS ##############
module load conda
conda activate liftoff_env
liftoff \
  GCA_041380795.1_Jhindsii_Rawlins_hap2_genomic.fna \
  GCA_001411555.2_Walnut_2.0_genomic.fna \
  -g Jregia.gff \
  -o JregiaGenes_lifted_to_Jhindsii.gff \
  -p 8 \
  -u unmapped_Jregia_genes.gff \
  -copies \
  -polish



cp JcalONLY_RDA10p_RoseVariables_terra_PC1_LatxLon_outliers_con.txt snps.txt
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, "SNP_" NR}' snps.txt > snps.bed
nano snps.bed  # OPEN FILE AND REMOVE LAST line
# Sort the SNP file
sort -k1,1 -k2,2n snps.bed > snps.sorted.bed
head snps.sorted.bed



grep -P '\t(gene|lnc_RNA)\t' JregiaGenes_lifted_to_Jhindsii.gff > genes_tmp.gff3
awk -F'\t' '{ 
    match($9, /ID=([^;]+)/, arr); 
    print $1"\t"$4-1"\t"$5"\t"arr[1] 
}' genes_tmp.gff3 > genes.bed

# Sort BED
sort -k1,1 -k2,2n genes.bed > JregiaReference.genes.bed


module load bedtools
bedtools closest -a snps.sorted.bed -b JregiaReference.genes.bed -D a -t first > JCsnps_nearest_genes.tsv
















library(terra)
library(dplyr)
library(ggplot2)


library(readr)
JCenv <- read_csv("/u/home/r/rcbuck/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/JcalONLY_coordsforGF.csv")
JCenv<-data.frame(JCenv)
Coordinates<-subset(JCenv, select=c(x,y))
Coordinates<-data.frame(Coordinates)


library(terra)
tif_files <- list.files('/u/home/r/rcbuck/project-vlsork/ClimateVariables/RoseVariables/Current_1981_2010', 
                        pattern = "\\.tif$", full.names = TRUE)
ras <- rast(tif_files)
layer_names <- gsub("Current_1981_2010_", "", gsub("\\.tif$", "", basename(tif_files)))
names(ras) <- layer_names
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
ras_8110 <- project(ras, new_crs, method = "bilinear")
crs(ras_8110)
ras_8110[ras_8110 == -9999] <- NA
remove_NAs_stack <- function(r) {
  # Create a mask identifying pixels where at least one layer is NOT NA
  mask_layer <- sum(!is.na(r)) > 0  # TRUE where at least one layer has data
  
  # Apply mask to keep only those pixels
  r <- mask(r, mask_layer, maskvalue=0)  
  
  return(r)
}
ras_8110<- remove_NAs_stack(ras_8110)

sample.coord.sp <- vect(Coordinates, geom = c("x", "y"), crs = crs(ras_8110))
sample_env <- terra::extract(ras_8110, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
sample_env <- cbind(JCenv, sample_env)

# Mask the climate rasters to species range
JCrange <- vect("~/project-vlsork/Juglans/trimmedFastqs/MarkedDuplicates/vcfs_bychr/Analyses/ALL_Jcal_Merge_10mBuffer.shp")
new_crs <- "EPSG:4326"  # Terra uses EPSG codes directly
JCrange <- terra::project(JCrange, new_crs)
range_env <- mask(ras_8110, JCrange)

# Convert all cells in the range to a dataframe
range_df <- as.data.frame(range_env, xy = TRUE, na.rm = TRUE)


range_df$group <- "Range"
sample_env$group <- "Sampled"

range_vars <- names(ras_8110)

env_compare <- bind_rows(
  range_df[, c(range_vars, "group")],
  sample_env[, c(range_vars, "group")]
)


library(tidyr)
env_compare_long <- env_compare %>%
  pivot_longer(
    cols = all_of(range_vars), 
    names_to = "variable", 
    values_to = "value"
  )




ggplot(env_compare_long, aes(x = value, fill = group)) +
  geom_density(alpha = 0.4, linewidth = 0.3) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  scale_fill_manual(
    values = c("Range" = "#BCA4E3",   # soft purple
               "Sampled" = "#F8766D") # soft red
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey40"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right"
  ) +
  labs(
    x = "Value",
    y = "Density",
    fill = "Source"
  )  
