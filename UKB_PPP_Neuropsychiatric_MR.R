#IV selection#
# Criteria:
# (1) Selected SNPs at P < 5e-8
# (2) Excluded MHC region (chr6: 25.5-34.0Mb)
# (3) LD clumping (r2 < 0.01 within 1Mb)
# (4) Selected SNPs with a F-statistic > 10
# (5) Excluded SNPs associated with more than five proteins (highly pleiotropic)

list <- list.files("/public/home/heyu/UKB-PPP/MR/pGWAS_summary/Clump_1Mb_0.01 (new)",pattern = "*.clumped")
clump_sum <- data.frame()
for (i in 1:length(list)){
  clump <- as.data.frame(fread(paste0("/public/home/heyu/UKB-PPP/MR/pGWAS_summary/Clump_1Mb_0.01 (new)/",list[i]),header = T))
  clump <- as.data.frame(clump[,3])
  colnames(clump)[1] <- "rsid"
  data <- as.data.frame(fread(paste0("/public/home/heyu/UKB-PPP/MR/pGWAS_summary/MHC_excluded/",gsub(".clumped",".txt",list[i])),header = T))
  clump <- merge(data,clump,by="rsid")
  clump_sum <- rbind(clump_sum,clump)
}
clump_sum <- clump_sum[which(clump_sum$F_statistic > 10),]
count <- aggregate(clump_sum$rsid, by = list(clump_sum$rsid), FUN = length)
count <- count[which(count$x <= 5),]
count <- as.data.frame(count[,1])
colnames(count)[1] <- "rsid"
clump_sum <- merge(count,clump_sum,by="rsid")

#Mendelian randomization#
library(TwoSampleMR)
mr_ins <- format_data(clump_sum,type = "exposure",header = T,phenotype_col = "PRO",snp_col = "rsid",beta_col = "BETA",se_col = "SE",eaf_col = "A1FREQ",effect_allele_col = "ALLELE1",other_allele_col = "ALLELE0",pval_col = "Pvalue")
mr_out <- extract_outcome_data(snps = mr_ins$SNP,outcomes = "ieu-b-2",proxies = FALSE)
pheno <- as.data.frame(fread("/public/home/heyu/FinnGen/AD/finngen_R10_G6_AD_WIDE.gz",header = T))
pheno$pheno <- "AD"
mr_out <- format_data(pheno,type = "outcome",header = T,snps = mr_ins$SNP,phenotype_col = "pheno",snp_col = "rsids",beta_col = "beta",se_col = "sebeta",effect_allele_col ="alt",other_allele_col = "ref",eaf_col = "af_alt",pval_col = "pval")
mr_harmo <- harmonise_data(exposure_dat = mr_ins,outcome_dat = mr_out)
try(mr_res <- mr(mr_harmo,method_list = c("mr_wald_ratio","mr_ivw")))
write.csv(mr_res,"/public/home/heyu/FinnGen/AD/G6_AD_WIDE.csv",row.names = F)

#Colocalization#
library(coloc)
library(GWAS.utils)
library(locuscomparer)
pheno <- as.data.frame(fread("/public/home/heyu/UKB-PPP/MR/Outcome/ieu-b-2.txt",header = T))
pheno$P <- 10^-(pheno$LP)
colnames(pheno)[12] <- "se"
data <- clump_sum[which(clump_sum$PRO=="SNAP25"),]
table(data$CHROM)
chri <- data[which(data$CHROM==19),]
Exp <- as.data.frame(fread("/public/home3/UKB-PPP_pGWAS/SNAP25_P60880_OID30811_v1_Neurology_II/discovery_chr19_SNAP25:P60880:OID30811:v1:Neurology_II.gz",header = T))
list <- c(chri$rsid)
for (i in 1:length(list)) {
  p_min <- chri[which(chri$rsid==list[i]),]
  exp <- Exp[which(Exp$GENPOS>=p_min$GENPOS-500000&Exp$GENPOS<=p_min$GENPOS+500000),]
  exp <- exp[,-14]
  exp$GP37 <- as.numeric(gsub("[^0-9]","",regmatches(exp$ID,regexpr(":(.*?):",exp$ID))))
  exp$pval <- 10^-(exp$LOG10P)
  exp <- merge(exp,rsid,by="ID")
  out <- pheno[which(pheno$seqnames == 19&pheno$start>=p_min$GP37-500000&pheno$start<=p_min$GP37+500000),]
  data1 <- exp[,c(6,8,10,11,14:16)]
  data2 <- out[,c(2,11:12,21,23)]
  harmo <- merge(data1,data2,by="rsid")
  harmo <- harmo[!duplicated(harmo$rsid),]
  col_res <- coloc.abf(dataset1 = list(beta=harmo$BETA,varbeta=(harmo$SE)^2,snp=harmo$rsid,position=harmo$GP37,type="quant",MAF=eaf2maf(harmo$A1FREQ),pvalues=harmo$pval,N=harmo$N),dataset2 = list(beta=harmo$ES,varbeta=(harmo$se)^2,snp=harmo$rsid,position=harmo$start,type="cc",pvalues=harmo$P))
}
exp <- harmo[,c(1,7)]
out <- harmo[,c(1,11)]
colnames(out)[2] <- "pval"
fwrite(exp,"D:/Yu Lab/UKB-PPP/COLOC/pQTL.txt",sep="\t",quote = F)
fwrite(out,"D:/Yu Lab/UKB-PPP/COLOC/GWAS.txt",sep="\t",quote = F)
locuscompare(in_fn1 = "D:/Yu Lab/UKB-PPP/COLOC/pQTL.txt",in_fn2 = "D:/Yu Lab/UKB-PPP/COLOC/GWAS.txt",title1 = "pQTL",title2 = "GWAS")