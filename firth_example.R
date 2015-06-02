library(Firth)

# smaller data files
load("C:/data/snpinfo/SNPInfo_HumanExome-12v1_rev5_analysisCols.Rdata")
snpNameCol <- "Name"
aggregateByCol <- "SKATgene"
chrCol <- "Chr"

# phenotype
## binomial
phenotype.file <- "C:/data/phenotype/AA_EC2_phen_hypert05.csv"
phenotype.formula <- "hypert05~v1age01+age_sq+sex+bmi01+pc1+pc3+pc5+pc7+pc8"
phenotype.id <- "id"

p <- read.csv(phenotype.file, header=TRUE)
pheno <- na.omit(p[ , c(phenotype.id, colnames(get_all_vars(as.formula(phenotype.formula), data=p)))])
rownames(pheno) <- pheno[ , phenotype.id]
## end binomial


genotype.files <- list.files(path="C:/data/genotype", pattern="^AA", full.names=TRUE)[1:2]

load(genotype.files[1])

subjects <- intersect(rownames(pheno), rownames(GT))
snps <- intersect(snpinfo[ , snpNameCol], colnames(GT))

chrs <- unique(snpinfo[snpinfo[ , snpNameCol] %in% snps, chrCol])

Zs <-  GT[subjects, snps]
SNPI <- snpinfo[snpinfo[ , chrCol] %in% chrs, ]
PHENO <- pheno[subjects, ]

