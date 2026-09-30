library(ieugwasr)
api_status()

library(data.table)

dt <- fread("/Users/artemisbriasouli/metal_bun_meta_all1.TBL.map.annot.gc")
head(dt)



# Install LDlinkR if you haven't already
install.packages("LDlinkR")

# Load the LDlinkR library
library(LDlinkR)

library(data.table)
library(TwoSampleMR)

dt <- fread("metal_eGFR_meta1.TBL.map.annot.gc.gz")

# remove missing rsIDs
dt <- dt[!is.na(RSID)]

# convert to data.frame
dt <- as.data.frame(dt)

outcome_dat <- format_data(
        dt,
        type = "outcome",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr.GC",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col = "Freq1",
        pval_col = "P.value.GC"
)

library(TwoSampleMR)

# remove missing rsIDs if any
dt <- dt[!is.na(RSID), ]

# convert from data.table to data.frame
dt <- as.data.frame(dt)

outcome_dat <- format_data(
        dt,
        type = "outcome",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col = "Freq1",
        pval_col = "P-value"
)

harm <- harmonise_data(
        exposure_dat = bmi_dec,
        outcome_dat = outcome_dat)


loo <- mr_leaveoneout(harm)


dt <- fread(cmd = "gunzip -c /Users/artemisbriasouli/Desktop/CKDGen_eGFR-decline_overall_adjDM.txt")


mr_leaveoneout_plot(loo)

# Define your SNPs
snp1 <- "rs7540681"
snp2 <- "your_target_SNP"

# Retrieve LD information
ld_data <- ld_pair(snp1 = snp1, snp2 = snp2)

# View LD information
print(ld_data)

your_ld_function


install.packages("snpStats")
library(snpStats)

install.packages("SNPRelate")
library(SNPRelate)




# Install necessary packages
install.packages("tidyverse")
install.packages("MendelianRandomization")
install.packages("MRInstruments")

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("remotes")
library(TwoSampleMR)
library(ggplot2)
library(MendelianRandomization)
library(tidyverse)

library(readxl)
bmi_gwas <- read_excel("/Users/artemisbriasouli/Desktop/BRAIN_ADIPOSE_WEIGHTED.xlsx")
heart_failure_gwas <- read.table("/Users/artemisbriasouli/Desktop/HERMES_Jan2019_HeartFailure_summary_data.txt", header = TRUE, sep = "\t")
# Your previous steps for loading and preparing data remain the same

clump_data(
        bmi_wood_gwas,
        clump_kb = 10000,
        clump_r2 = 0.001,
        clump_p1 = 1,
        clump_p2 = 1,
        pop = "EUR",
        bfile = NULL,
        plink_bin = NULL
)

BMI_BIG <- read_exposure_data("/Users/artemisbriasouli/Not_iCloud/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt", clump = FALSE, sep = "\t", phenotype_col = "",
                              snp_col = "SNP",
                              beta_col = "BETA",
                              se_col = "SE",
                              eaf_col = "Freq_Tested_Allele_in_HRS",
                              effect_allele_col = "Tested_Allele",
                              other_allele_col = "Other_Allele",
                              pval_col = "P",
                              units_col = "",
                              ncase_col = "",
                              ncontrol_col = "",
                              samplesize_col = "N",
                              gene_col ="",
                              id_col = "",
                              min_pval = 1e-200,
                              log_pval = FALSE,
                              chr_col = "CHR",
                              pos_col = "POS"
)

CLUMP_BMI <- clump_data(
        BMI_BIG,
        clump_kb = 10000,
        clump_r2 = 0.01,
        clump_p1 = 5e-08,
        clump_p2 = 1,
        pop = "EUR",
        bfile = NULL,
        plink_bin = NULL
)

data_decline <- data.table::fread("/Users/artemisbriasouli/Downloads/CKDGen_eGFR-decline_overall_adjDM.txt", header = TRUE, sep = "\t")
nrow(data)
data <- data[grepl("^rs", data$SNP), ]

decline_data <- read_outcome_data(
        filename = "/Users/artemisbriasouli/Downloads/CKDGen_eGFR-decline_overall_adjDM.txt",
        snps = CLUMP_BMI$SNP,
        sep = "\t",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P-value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "N_total_sum"
)

str(dt)

detach("package:TwoSampleMR", unload = TRUE)
detach("package:plotly", unload = TRUE)




EGFR_bmi <- read_outcome_data("/Users/artemisbriasouli/Downloads/CKDGen_eGFR-decline_overall_adjDM.txt",
        snps = final_bmi_tissue_partitioned_data$SNP,
        sep = "\t",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P-value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "N_total_sum"
)

snps_adipose$R2_snp <-
        2 * snps_adipose$eaf.exposure *
        (1 - snps_adipose$eaf.exposure) *
        snps_adipose$beta.exposure^2

R2_adipose <- sum(snps_adipose$R2_snp, na.rm = TRUE)

R2_adipose

harmonise_BMI_EGFR <- harmonise_data(metabolic_data, EGFR_bmi, action = 1)

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')

if (is.null(bmi2018_exp_dat)) {
        stop("The result from extract_instruments is NULL.")
}
if (!is.data.frame(bmi2018_exp_dat)) {
        stop("The result from extract_instruments is not a data frame.")
}

bmi2018_exp_dat <- extract_instruments(outcomes = 'ieu-b-40',  
                                       p1 = 5e-08,
                                       clump = TRUE,
                                       p2 = 5e-08,
                                       r2 = 0.001,
                                       kb = 10000,
                                       opengwas_jwt = ieugwasr::get_opengwas_jwt(),
                                       force_server = FALSE)

print(bmi2018_exp_dat)
str(bmi2018_exp_dat)



bmi2018_exp_dat <- extract_instruments(outcomes = 'ieu-b-40',  p1 = 5e-08,
                                       clump = TRUE,
                                       p2 = 5e-08,
                                       r2 = 0.001,
                                       kb = 10000,
                                       opengwas_jwt = ieugwasr::get_opengwas_jwt(),
                                       force_server = FALSE)

bmi2018_exp_dat <- extract_instruments(outcomes = 'ieu-b-40',  p1 = 5e-08,
                                       clump = TRUE,
                                       p2 = 5e-08,
                                       r2 = 0.001,
                                       kb = 10000,
                                       opengwas_jwt = ieugwasr::get_opengwas_jwt(),
                                       force_server = FALSE)
str(bmi2018_exp_dat)
bmi2018_exp_dat_r2 <- clump_data(bmi2018_exp_dat, clump_r2 = 0.001)
get_r_from_bmaf(bmi2018_exp_dat_r2$beta.exposure, bmi2018_exp_dat_r2$eaf.exposure)

input <- read.csv(temp_file)
read_exposure_data(
        bmi_gwas,
        sep = " ",
        phenotype_col = "Phenotype",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval",
        units_col = "units",
        ncase_col = "ncase",
        ncontrol_col = "ncontrol",
        samplesize_col = "samplesize",
        gene_col = "gene",
        id_col = "id",
        min_pval = 1e-200,
        log_pval = FALSE,
        pval_threshold = 5e-08,
        plink_bin = NULL,
        bfile = NULL,
        clump_r2 = 0.001,
        clump_kb = 10000,
        clump_p1 = 5e-08,
        clump = TRUE,
        pop = "EUR",
        harmonise_strictness = 1
)


# Perform multivariable MR analysis
mr_results <- mr_iv(
        exposure = bmi ~ beta + se + samplesize,
        outcome = heart_failure ~ b + se + N,
        data = mr_data
)

# View results
summary(mr_results)

# Select relevant columns
bmi_data <- bmi_gwas %>%
        select(SNP, effect_allele, se, samplesize, beta_brain, beta_adipose, pval) %>%
        rename(Exposure_SNP = SNP) %>%
        mutate(new_column= BMI_brain)

heart_failure_data <- heart_failure_gwas %>%
          select(SNP, freq, A1, A2, b, se, N) %>%
           rename(Outcome_SNP = SNP)

# Merge datasets
mr_data <- merge(bmi_data, heart_failure_data, by = c("Exposure_SNP", "Outcome_SNP"), all = TRUE)

# Check column names
colnames(bmi_data)
colnames(heart_failure_data)

exposure_dat <- mv_extract_exposures_local(input_data, sep=" ", phenotype_col - "ID", snp_col = "SNP", beta_col_1 = "beta_brain")

# Check uniqueness of the merge columns
table(duplicated(bmi_data$Exposure_SNP))
table(duplicated(heart_failure_data$Outcome_SNP))

# Check data types
str(bmi_data$Exposure_SNP)
str(heart_failure_data$Outcome_SNP)

# Check for duplicates in Exposure_SNP and Outcome_SNP columns
any(duplicated(bmi_data$Exposure_SNP))
any(duplicated(heart_failure_data$Outcome_SNP))

library(dplyr)

# Merge datasets using dplyr's merge function

mr_results <- mr_ivw(
        exposures = list(bmi = bmi_data$beta,
                         se_bmi = bmi_data$se,
                         samplesize_bmi = bmi_data$samplesize),
        outcomes = list(heart_failure = heart_failure_data$b,
                        se_heart_failure = heart_failure_data$se,
                        N_heart_failure = heart_failure_data$N)
)

# View results
summary(mr_results)

mr_results <- mr_multivar(
        exposure = bmi ~ beta + se + samplesize,
        outcome = heart_failure ~ b + se + N,
        data = mr_data
)

bmi_instrument <- mr_instrument(beta = bmi_data$beta,
                                se = bmi_data$se,
                                n = bmi_data$samplesize)
heart_failure_instrument <- mr_instrument(beta = heart_failure_data$b,
                                          se = heart_failure_data$se,
                                          n = heart_failure_data$N)

multivar_mr <- mr_multivar(instrument = list(bmi_instrument, heart_failure_instrument))

# View results
summary(multivar_mr)

#examples
id_exposure <- c("ieu-a-299", "ieu-a-300", "ieu-a-302")
id_outcome <- "ieu-a-7"


new_data <- bmi_gwas %>%
        select(SNP, effect_allele, other_allele, eaf, beta_brain, se, samplesize, pval) %>%
        mutate(ID = "bmi_brain") 

new_data_1 <- bmi_gwas %>%
         select(SNP, effect_allele, other_allele, eaf, beta_adipose, se, samplesize, pval) %>%
        mutate(ID = "bmi_adipose") 

# Merge new_data and new_data_1 datasets
merged_data <- merge(new_data, new_data_1, all = TRUE)

# View the structure of the merged dataset
str(merged_data)

new_data <- new_data %>%
        rename(beta = beta_brain)

library(dplyr)

new_data_1 <- new_data_1 %>%
        rename(beta = beta_adipose)

combined_data <- merged_data %>%
        rename(SNP = SNP,
               beta = beta.outcome,
               se = se.outcome,
               effect_allele = effect_allele.outcome,
               other_allele = other_allele.outcome,
               eaf = eaf.outcome,
               ID = outcome,
               pval= pval.exposure)

column_mapping <- c(
        SNP = "SNP",
        beta = "beta.outcome",
        se = "se.outcome",
        effect_allele = "effect_allele.outcome",
        other_allele = "other_allele.outcome",
        eaf = "eaf.outcome",
        ID = "exposure_dat",
        beta_exposure = "beta.exposure",
        se_exposure = "se.exposure",
        effect_allele_exposure = "effect_allele.exposure",
        other_allele_exposure = "other_allele.exposure",
        eaf_exposure = "eaf.exposure"
)

combined_data <- combined_data %>%
        rename_with(~column_mapping, matches("^(SNP|beta\\.outcome|se\\.outcome|effect_allele\\.outcome|other_allele\\.outcome|eaf\\.outcome|outcome|exposure_dat|beta\\.exposure|se\\.exposure|effect_allele\\.exposure|other_allele\\.exposure|eaf\\.exposure)$"))

new_column_names_1 <- c(
        "SNP",
        "eaf.outcome",
        "effect_allele.outcome",
        "other_allele.outcome",
        "beta.outcome",
        "se.outcome",
        "size.outcome"
)

exposure_dat <- mv_extract_exposures_local ("/Users/artemisbriasouli/stacked_data.csv", sep = ",", phenotype_col = "exposure",
                                                                                   snp_col = "ID",
                                                                                       beta_col = "beta.exposure",
                                                                                       se_col = "se.exposure",
                                                                                       eaf_col = "eaf.exposure",
                                                                                       effect_allele_col = "effect_allele.exposure",
                                                                                       other_allele_col = "other_allele.exposure",
                                                                                       pval_col = "pval.exposure",
                                                                                       samplesize_col = "size.exposure",
                                                                                       min_pval = 1e-200,
                                                                                       log_pval = FALSE,
                                                                                       pval_threshold = 5e-08,
                                                                                       clump_r2 = 0.001,
                                                                                       clump_kb = 10000,
                                                                                       harmonise_strictness = 1)

exposure_dat <- mv_extract_exposures_local(input_data)

# Convert dataframe to character string with newline characters
combined_data_string <- capture.output(write.table(combined_data, row.names = FALSE, sep = ",", col.names = TRUE))

# Join the lines into a single character string
combined_data_string <- paste(combined_data_string, collapse = "\n")

# Write combined_data to a temporary CSV file
temp_file <- tempfile(fileext = ".csv")
write.csv(merged_data, file = temp_file, row.names = FALSE)

# Use the file name as input
input_data <- temp_file

# Write combined_data to a temporary CSV file
temp_file <- tempfile(fileext = ".csv")
write.csv(combined_data, file = temp_file, row.names = FALSE)

# Use the file name as input
input_data <- temp_file

new_column_names <- c(
     "SNP",
        "effect_allele.exposure",
           "other_allele.exposure",
            "eaf.exposure",
              "beta.exposure",
            "se.exposure",
             "size.exposure",
           "pval.exposure",
               "exposure"
 )


# Assuming your dataframe is named 'df' and the current column names are 'col1', 'col2', etc.

# Rename the columns individually
filtered_snps <- filtered_snps %>%
        rename(
               beta.exposure = beta,
               se.exposure = se,
               effect_allele.exposure = effect_allele,
)               other_allele.exposure = other_allele)

# Add the new column with "brain" for each row
df$id.exposure <- "brain"




outcome_dat <- read_outcome_data(
        "/Users/artemisbriasouli/pattaro_gwas.txt",
        snps = data$SNP,
        sep = ",",
        phenotype_col = "",
        snp_col = "rsID",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "freqA1",
        effect_allele_col = "allele1",
        other_allele_col = "allele2",
        pval_col = "pval",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "N")

# 1. Filter outcome SNPs to only those found in the exposure data
outcome_filtered <- decline_uni_outcome[decline_uni_outcome$SNP %in% wide_df$SNP, ]

# 2. Merge exposure and outcome data by SNP
merged_df <- wide_df %>%
        inner_join(outcome_filtered, by = "SNP")

merged_df <- wide_df %>%
        inner_join(outcome_dat, by = "SNP")


# 3. Extract matrices
bx <- as.matrix(merged_df[, c("beta.exposure.metabolic", "beta.exposure.mental")])
bxse <- as.matrix(merged_df[, c("se.exposure.metabolic", "se.exposure.mental")])
by <- merged_df$beta.outcome
byse <- merged_df$se.outcome
snps <- merged_df$SNP

# Optional: Allele info and EAFs
effect_allele <- merged_df[, c("effect_allele.exposure.metabolic", "effect_allele.exposure.mental")]
other_allele <- merged_df[, c("other_allele.exposure.metabolic", "other_allele.exposure.mental")]
eaf <- merged_df[, c("eaf.exposure.metabolic", "eaf.exposure.mental")]

# 4. Create MVMR input object
mvdat <- mr_mvinput(
        bx = bx,
        bxse = bxse,
        by = by,
        byse = byse,
        snps = snps,
        exposure = c("metabolic", "mental"),
        outcome = "your_outcome_name",
        effect_allele = effect_allele,
        other_allele = other_allele,
        eaf = eaf
)


decline_outcome <- read_outcome_data(
        "/Users/artemisbriasouli/Downloads/CKDGen_eGFR-decline_CKD.txt",
        snps = CLUMP_BMI$SNP,
        sep = "\t",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P-value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "N_total_sum")

Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum

wuttke_egfr_outcome_data <- read_outcome_data("/Users/artemisbriasouli/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt",
        snps = brain_snps$SNP,
        sep = ",",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P-value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "n_total_sum")

#Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum

mv_outcome_dat <- extract_outcome_data(brain_and_adipose$SNP, ieu-a-1120)


outcome_dat <- extract_outcome_data(exposure_dat$SNP, heart_failure_data_input,
proxies = 1, rsq = 0.6, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

Y <- mvdat_dec$outcome_beta
SE_Y <- mvdat_dec$outcome_se

X1 <- data.frame(mvdat_dec[["exposure_beta"]])$bmi_adipose
X2 <- data.frame(mvdat_dec[["exposure_beta"]])$bmi_brain

SE_X1 <- data.frame(mvdat_dec[["exposure_se"]])$bmi_adipose
SE_X2 <- data.frame(mvdat_dec[["exposure_se"]])$bmi_brain

SNP <- data.frame(mvdat_ckd_wuttke[["exposure_SE"]])$1

# Convert row names to a column in the dataframe
df_with_row_as_column <- cbind(RowName = rownames(mvdat$exposure_se), mvdat$exposure_se)

# View the dataframe
print(df_with_row_as_column)

# Now, the first column can be accessed as a normal column
first_column <- df_with_row_as_column[, "RowName"]
print(first_column)
# Extract row names
first_column <- rownames(mvdat$exposure_se)

# Print the result
print(first_column)

# Example data extraction (adjust these according to your actual data structure)
SNP <- rownames(df$exposure_se)        # Extract SNPs (assuming they are row names)
X1 <- df$exposure_se[, "bmi_adipose"]  # Extract first exposure
X2 <- df$exposure_se[, "bmi_brain"]    # Extract second exposure
Y <- df$outcome$Y                      # Extract outcome (adjust if needed)
SE_Y <- df$outcome$SE_Y                # Extract standard error for Y (adjust if needed)
SE_X1 <- df$exposure_se[, "bmi_adipose_se"]  # Extract standard error for X1 (adjust if needed)
SE_X2 <- df$exposure_se[, "bmi_brain_se"]    # Extract standard error for X2 (adjust if needed)

# Combine into a new dataframe
raw_data_cys <- data.frame(X1, X2, Y, SE_Y, SE_X1, SE_X2)

# View the new dataframe
print(raw_data)

F.data <- format_mvmr(BXGs = raw_data_cys[,c(2,3)],
                                        BYG = raw_data[,4],
                                        seBXGs = raw_data[,c(6,7)],
                                        seBYG = raw_data[,5],
                                        RSID = raw_data[,1])


mvmrcovmatrix<-matrix(c(1,-0.1,-0.05,-0.1,1,0.2,-0.05,0.2,1), nrow = 3, ncol = 3)
Xcovmat<-phenocov_mvmr(mvmrcovmatrix,F.data[,7:9])

sres2 <- strength_mvmr(r_input = F.data, gencov = Xcovmat)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)

res1 <- qhet_mvmr(F.data, mvmrcovmatrix, CI = F, iterations = 100)

#genevievecodeupdate

mvdat_dec <- mv_harmonise_data(final_bmi_tissue_partitioned_data, egfr_dec_outcome_dat, harmonise_strictness = 1)

F.data_cys <- format_mvmr(BXGs = raw_data_cys[,c(1,2)], BYG = raw_data_cys$Y, seBXGs = raw_data_cys[,c(5,6)], seBYG = raw_data_cys$SE_Y)



pheno_mat <- cor(raw_data_cys[,c(1,2)])

seBXGs=raw_data_cys[,c(5,6)]

pos_cys <- phenocov_mvmr(pheno_mat, seBXGs)

fstats_cys <- strength_mvmr(r_input=F.data_cys,gencov=pos_cys)[1,]


#genevieve_code
F.data<-format_mvmr(BXGs=mvdat$exposure_beta[,c(1,2)],
                    
                    BYG=mvdat$outcome_beta,
                    
                    seBXGs=mvdat$exposure_se[,c(1,2)],
                    
                    seBYG=mvdat$outcome_se,
                    
                    RSID=row.names(mvdat$exposure_beta))



pheno_mat <- cor(mvdat$exposure_beta[,c(1,2)])

seBXGs=mvdat$exposure_se[,c(1,2)]

pos <- phenocov_mvmr(pheno_mat, seBXGs)

fstats <- strength_mvmr(r_input=F.data,gencov=pos)[1,]




#need to convert the relevant columns of your data frames into matrices before passing them to mr_mvinput

bx_matrix <- as.matrix(cbind(X1, X2))
bxse_matrix <- as.matrix(cbind(SE_X1, SE_X2))

Y_vector <- as.vector(Y)
SE_Y_vector <- as.vector(SE_Y)

mrmvinput <- mr_mvinput(bx = bx_matrix, bxse = bxse_matrix, by = Y_vector, byse = SE_Y_vector,  exposure = "exposure",
                        outcome = "outcome",
                        snps = "snp",
                        effect_allele = NA,
                        other_allele = NA,
                        eaf = NA)

mr_mvinput(
        bx = matrix(),
        bxse = matrix(),
        by = 0,
        byse = 0,
        correlation = matrix(),
        exposure = "exposure",
        outcome = "outcome",
        snps = "snp",
        effect_allele = NA,
        other_allele = NA,
        eaf = NA
)

# Load required libraries
library(tidyr)
library(dplyr)

# Load your data
hea
# Reshape to wide format: one row per SNP, with separate columns for each exposure
wide_df <- df %>%
        pivot_wider(
                id_cols = SNP,
                names_from = exposure,
                values_from = c(beta.exposure, se.exposure),
                names_glue = "{.value}.{exposure}"
        )

# Extract matrices for mr_mvinput
bx <- as.matrix(wide_df[, c("beta.exposure.metabolic", "beta.exposure.mental")])
bxse <- as.matrix(wide_df[, c("se.exposure.metabolic", "se.exposure.mental")])
snps <- wide_df$SNP

# Optional: Check number of SNPs and preview data
print(length(snps))
head(bx)
head(bxse)




mr_egger_result <- mr_mvegger(
        mrmvinput,
        orientate = 1,
        correl = FALSE,
        distribution = "normal",
        alpha = 0.05
)

# Estimates and standard errors
estimates <- c(0.126, 0.194, 0.000)  # The intercept is typically not plotted, but included here for reference
std_errors <- c(0.260, 0.278, 0.004)

# Calculate 95% confidence intervals
lower_ci <- estimates - 1.96 * std_errors
upper_ci <- estimates + 1.96 * std_errors

# Labels for the exposures
exposures <- c("Exposure 1", "Exposure 2", "Intercept")

# Set up the plot area
plot(estimates, 1:length(estimates), 
     xlim = c(min(lower_ci), max(upper_ci)), 
     pch = 19, xlab = "Effect Size", 
     ylab = "Exposures", 
     main = "Multivariable MR-Egger Forest Plot", 
     yaxt = "n")  # Suppress y-axis labels for custom labels

# Add custom y-axis labels
axis(2, at = 1:length(exposures), labels = exposures, las = 1)

# Adding error bars for confidence intervals
arrows(x0 = lower_ci, y0 = 1:length(estimates), 
       x1 = upper_ci, y1 = 1:length(estimates), 
       code = 3, angle = 90, length = 0.05)

# Optionally, add a grid
abline(v = 0, lty = 2, col = "grey")

# Estimates and standard errors for each exposure
estimates <- c(0.126, 0.194)  # Effect sizes for Exposure 1 and Exposure 2
std_errors <- c(0.260, 0.278)  # Corresponding standard errors

# Calculate 95% confidence intervals
lower_ci <- estimates - 1.96 * std_errors
upper_ci <- estimates + 1.96 * std_errors

# Labels for the exposures
exposures <- c("bmi_adipose", "bmi_brain")

# Set up the plot area
plot(estimates, 1:length(estimates), 
     xlim = c(min(lower_ci), max(upper_ci)), 
     pch = 19, xlab = "Causal Estimate (Effect Size)", 
     ylab = "Exposures", 
     main = "Multivariable MR Forest Plot", 
     yaxt = "n")  # Suppress y-axis labels for custom labels

# Add custom y-axis labels
axis(2, at = 1:length(exposures), labels = exposures, las = 1)

# Adding error bars for confidence intervals
arrows(x0 = lower_ci, y0 = 1:length(estimates), 
       x1 = upper_ci, y1 = 1:length(estimates), 
       code = 3, angle = 90, length = 0.05)

# Optionally, add a reference line at zero effect
abline(v = 0, lty = 2, col = "grey")



# Example using base R to construct a forest plot
plot(mr_egger_result@Exposure1$Estimate, # X-axis: Effect size of exposure 1
     1:length(mr_egger_result@Exposure1$Estimate), # Y-axis: SNPs
     xlim = c(-1, 1), # Adjust this range based on your data
     pch = 19, # Use solid circles for points
     xlab = "Effect Size (Exposure 1)", 
     ylab = "SNPs", 
     main = "MVMR Forest Plot: Exposure 1")

# Add error bars for confidence intervals
arrows(x0 = mr_egger_result@Exposure1$CI.Lower,
       y0 = 1:length(mr_egger_result@Exposure1$Estimate),
       x1 = mr_egger_result@Exposure1$CI.Upper,
       y1 = 1:length(mr_egger_result@Exposure1$Estimate),
       code = 3, angle = 90, length = 0.05)




ivreg_results <- ivreg(Y ~ X[,1] + X[,2] - 1, weights = 1 / SE_Y^2)
summary(ivreg_results)

exposure1 <- data.frame(mvdat[["exposure_beta"]])$bmi_adipose
exposure2 <- data.frame(mvdat[["exposure_beta"]])$bmi_brain

install.packages("VariantAnnotation")
library(VariantAnnotation)

install.packages("vcFR")
library(VCFR)

install.packages("VariantTools")
library(VariantTools)

library(VariantAnnotation)

vcf_file <- "/Users/artemisbriasouli/Desktop/ebi-a-GCST003374.vcf"
vcf_data <- readVcf(vcf_file)

# Load necessary packages
library(VariantAnnotation)

# Read VCF file into R
vcf_file <- "your_file.vcf"
vcf_data <- readVcf(vcf_file)

# View summary of vcf_data
summary(vcf_data)

# Extract SNP identifiers
snps <- rownames(vcf_data)

# Extract genotype data (if needed)
genotypes <- geno(vcf_data)

# Extract other relevant information (e.g., effect sizes, standard errors, allele frequencies)
# Example: Extract INFO field
info_field <- info(vcf_data)

# Subset data if needed (e.g., to include only specific SNPs)
# Example: Subset data to include only SNPs with minor allele frequency > 0.05
maf_threshold <- 0.05
snps_filtered <- snps[info_field$AF > maf_threshold]
vcf_data_filtered <- vcf_data[snps_filtered]

# Merge with exposure data (if needed)
# Example: Merge with exposure data from a separate source using SNP identifiers
merged_data <- merge(vcf_data, exposure_dat, by = "SNP", all = TRUE)

# Perform further analysis (e.g., Mendelian randomization)
# Example: Perform MR analysis using appropriate methods and packages

# Evaluate results and interpret findings

# Clean up if necessary
rm(vcf_data)  # Remove vcf_data object from memory

# Load necessary packages
library(VariantAnnotation)

# Read VCF file into R
vcf_file <- "your_file.vcf"
vcf_data <- readVcf(vcf_file)

# Extract INFO field
info_field <- info(vcf_data)

# Specify the path to your VCF file
vcf_file <- "/Users/artemisbriasouli/Desktop/ebi-a-GCST003374.vcf"

# Read the VCF file
vcf <- read.vcfR(vcf_file)

# View the VCF object
print(vcf)

# Summarize the VCF object
summary(vcf)

# View the first few rows of the VCF object
head(vcf@fix)
head(vcf@gt)

# Extract the fixed fields and the genotype matrix as data frames
fixed_df <- as.data.frame(vcf@fix)
gt_df <- as.data.frame(vcf@gt)

# Combine fixed fields and genotype matrix into one data frame
vcf_df <- cbind(fixed_df, gt_df)

# Specify the path for the output TXT file
output_txt_file <- "/Users/artemisbriasouli/Desktop/ebi-a-GCST003374.txt"

# Write the combined data frame to a TXT file
write.table(vcf_df, file = output_txt_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Sample data: matrix format
data <- matrix(c("ES:SE:LP:AF:ID", "0.059:0.064:0.455932:0.032:rs12565286", 
                 "ES:SE:LP:AF:ID", "0.075:0.066:0.60206:0.034:rs11804171",  
                 "ES:SE:LP:AF:ID", "-0.067:0.068:0.49485:0.934:rs2977670",  
                 "ES:SE:LP:AF:ID", "-0.11:0.083:0.721246:0.062:rs12138618", 
                 "ES:SE:LP:AF:ID", "-0.02:0.028:0.327902:0.835:rs3094315",  
                 "ES:SE:LP:AF:ID", "0.0081:0.036:0.0861861:0.885:rs3131968"), 
               ncol = 2, byrow = TRUE)

# Extract the column names
column_names <- strsplit(data[1, 1], ":")[[1]]

# Split the data column
split_data <- t(sapply(data[, 2], function(x) strsplit(x, ":")[[1]]))

# Convert to a data frame
df <- data.frame(split_data, stringsAsFactors = FALSE)
colnames(df) <- column_names

# Specify the path for the output TXT file
output_txt_file <- "/Users/artemisbriasouli/Desktop/ebi-a-GCST003374_1.txt"

# Write the data frame to a TXT file
write.table(df, file = output_txt_file, sep = "\t", row.names = FALSE, quote = FALSE)

# Print the data frame to the console
print(df)



# Load necessary libraries
library(forestplot)

# Example data (replace with your own data)
mr_results <- data.frame(
        SNP = c("rs123456", "rs654321", "rs789012"),
        Trait = c("Trait A", "Trait B", "Trait C"),
        Effect = c(0.2, -0.1, 0.3),  # Effect estimates
        CI_lower = c(0.1, -0.2, 0.2),  # Lower bound of confidence interval
        CI_upper = c(0.3, 0, 0.4)      # Upper bound of confidence interval
)

# Create a forest plot
forestplot(
        mean = mr_results$Effect,
        lower = mr_results$CI_lower,
        upper = mr_results$CI_upper,
        labels = mr_results$SNP,
        title = "Multivariable MR Results",
        xlab = "Effect Size",
        is.summary = c(TRUE, FALSE, TRUE),
        col = fpColors(box = "royalblue", lines = "black", zero = "black"),
        txt_gp = fpTxtGp(label = gpar(fontsize = 10)),
        xticks = c(-0.5, 0, 0.5),
        new_page = TRUE
)

# Extract necessary information (e.g., effect sizes, standard errors, allele frequencies)
# Assuming your VCF file contains these annotations in the INFO field

# Example: Extract effect sizes (beta), standard errors (se), and allele frequencies (af)
effect_sizes <- info_field$BETA
standard_errors <- info_field$SE
allele_frequencies <- info_field$AF

# View extracted information
head(effect_sizes)
head(standard_errors)
head(allele_frequencies)

mvdat_ckd_wuttke <- mv_harmonise_data(final_bmi_tissue_partitioned_data, wuttke_outcome_data, action=1)
mvdat_ckd_wuttke <- mv_harmonise_data(final_bmi_tissue_partitioned_data, wuttke_outcome_data, harmonise_strictness=1)
mvdat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
View(mvdat_ckd_wuttke)
res_wuttke <- mv_multiple(mvdat_ckd_wuttke)
View(res_wuttke)


mv_egfr_dec <- mv_harmonise_data(data, ckd_new_outcome_dat, action = 1)

mv_egfr_new <- mv_harmonise_data(final_bmi_tissue_partitioned_data, wuutke_egfr_outcome_data, harmonise_strictness = 1)

# Example beta and standard error values
beta <- -0.05122421	
	
bmi_dec$exposure <- "bmi"
bmi_dec$id.exposure <- "bmi"	

 # Beta estimate from Mendelian randomization
se <- 0.002241317 # Standard error of the beta

# Step 1: Calculate the 95% Confidence Interval
lower_CI <- beta - 1.96 * se
upper_CI <- beta + 1.96 * se

# Step 2: Display the results
lower_CI
upper_CI

Chr Pos_b37 RSID Allele1 Allele2 Freq1 Effect StdErr P-value n_total_sum
# Beta and standard error from MR results
beta <- 0.15106545 						
	

	
	



 # Example beta value
se <- 0.02898718

# Calculate the odds ratio (OR)
OR <- exp(beta)

# Calculate the 95% confidence intervals (CI)
lower_CI <- exp(beta - 1.96 * se)
upper_CI <- exp(beta + 1.96 * se)

# Output the results
cat("Odds Ratio (OR):", OR, "\n")
cat("95% Confidence Interval: (", lower_CI, ",", upper_CI, ")\n")

# Beta and standard error from MR results
beta <- 0.15106545 						

    # Example beta value
se <- 0.02898718
# Calculate the 95% confidence intervals (CI)
lower_CI_beta <- beta - 1.96 * se
upper_CI_beta <- beta + 1.96 * se

# Output the results
cat("95% Confidence Interval for Beta: (", lower_CI_beta, ",", upper_CI_beta, ")\n")

ieugwasr::ld_clump()

install.packages("dplyr")
library(openxlsx)
library(dplyr)
library(openxlsx)
metabolic_data <- read.xlsx("/Users/artemisbriasouli/Not_iCloud/BOOK3-MR.xlsx")

filtered_snps_brain <- gwas_data %>%
        filter(PPA4_brain > 0.8)

snps_adipose <- read_excel("/Users/artemisbriasouli/filtered_snps.xlsx", sheet = 1)  # Read the first sheet


write.xlsx(filtered_snps_brain, 
           library(readxl)
)

library(TwoSampleMR)
library(ggplot2)

above_threshold_file <- paste0("file://", file.path(getwd(), "above_threshold.csv"))
below_threshold_file <- paste0("file://", file.path(getwd(), "below_threshold.csv"))


brain_exp <- read_exposure_data(
        "below_threshold.csv",
        clump = FALSE,
        sep = ",",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval",
        samplesize_col = "samplesize",
        min_pval = 1e-200,
        log_pval = FALSE,
        chr_col = "chr"
)


"/Users/artemisbriasouli/Desktop/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt",
snps = final_bmi_tissue_partitioned_data$SNP,
sep = ",",
phenotype_col = "",
snp_col = "RSID",
beta_col = "Effect",
se_col = "StdErr",
eaf_col = "Freq1",
effect_allele_col = "Allele1",
other_allele_col = "Allele2",
pval_col = "P-value",
units_col = "",
ncase_col = "",
ncontrol_col = "",
samplesize_col = "n_total_sum"

chd_out_dat_brain <- extract_outcome_data(snps = brain_exp$SNP, outcomes = 'ieu-a-7')

chd_out_dat_adipose <- extract_outcome_data(snps = adipose_exp$SNP, outcomes = 'ieu-a-7')

dat_metab <- harmonise_data(decline_outcome, CLUMP_BMI, action = 1)
dat_brain <- harmonise_data(brain_exp, chd_out_dat_brain, action = 1)

res <- mr(dat_metab, parameters = default_parameters(),
          method_list = subset(mr_method_list(), use_by_default)$obj
)

res_brain <- mr(dat_brain, parameters = default_parameters(),
                method_list = subset(mr_method_list(), use_by_default)$obj
)

afib <- extract_outcome_data(snps = brain_exp$SNP, outcomes = 'ebi-a-GCST006414')
afib_adipose <- extract_outcome_data(snps = adipose_exp$SNP, outcomes = 'ebi-a-GCST006414')
afib_dat <- harmonise_data(brain_exp, afib, action = 2)
afib_adipose_dat <- harmonise_data(adipose_exp, afib_adipose, action = 2)
res_afib_adipose <- mr(afib_adipose_dat)
res_afib <- mr(afib_dat)
write.csv(above_threshold, above_threshold_file, row.names = FALSE)
write.csv(below_threshold, below_threshold_file, row.names = FALSE)
filtered_adipose_snps_filename <- "filtered_adipose_snps.xlsx"

convert_to_filename <- function(string) {
        # Remove invalid characters
        filename <- gsub("[^A-Za-z0-9._-]", "_", string)
        # Remove leading and trailing underscores
        filename <- gsub("^_+|_+$", "", filename)
        return(filename)
}


# Assuming your data frame is called outcome_df and the OR column is named "OR"
decline_uni_outcome$beta <- log(decline_uni_outcome$OR)

# Step 1: Check if OR column has any non-numeric values
# Convert OR to numeric, coercing non-numeric values to NA
decline_uni_outcome$OR <- as.numeric(decline_uni_outcome$OR)

# Step 2: Remove or handle NA values
# Option A: Remove rows where OR is NA
outcome_df <- outcome_df[!is.na(outcome_df$OR), ]

# Option B: Handle NA values (replace them, if needed, or keep them if it’s acceptable)

# Step 3: Calculate beta from OR
outcome_df$beta <- log(outcome_df$beta.outcome)

library(dplyr)

# Assuming your data frame is named df
decline_uni_outcome <- decline_uni_outcome %>% select(-beta.outcome)



decline_uni_outcome_1 <- read_outcome_data(
        "/Users/artemisbriasouli/Desktop/cys.txt",
        snps = bmi_dec$SNP,
        sep = "\t",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P.value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "n"
)

decline_uni_outcome_21 <- read_outcome_data(
        "/Users/artemisbriasouli/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt",
        snps = snps_brain$SNP,
        sep = " ",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P-value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "n_total_sum"
)


decline_uni_outcome_3 <- read_outcome_data(
        "/Users/artemisbriasouli/pattaro_gwas.txt",
        snps = bmi_dec$SNP,
        sep = ",",
        phenotype_col = "",
        snp_col = "rsID",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "freqA1",
        effect_allele_col = "allele1",
        other_allele_col = "allele2",
        pval_col = "pval",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "N"
)

decline_uni_outcome_4 <- read_outcome_data(
        "/Users/artemisbriasouli/bun.txt",
        snps = bmi_dec$SNP,
        sep = "\t",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P.value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "n"
)








sep = ",",
phenotype_col = "",
snp_col = "RSID",
beta_col = "Effect",
se_col = "StdErr",
eaf_col = "Freq1",
effect_allele_col = "Allele1",
other_allele_col = "Allele2",
pval_col = "P-value",
units_col = "",
ncase_col = "",
ncontrol_col = "",
samplesize_col = "n_total_sum"

input_string <- 'above_threshold.csv'
filename <- convert_to_filename(input_string)
print(filename)


above_threshold_filename <- create_valid_filename("above_threshold.csv")
below_threshold_filename <- create_valid_filename("below_threshold.csv")


file_path <- "/Users/artemisbriasouli/Not_iCloud/filtered_adipose_snps.xlsx"  # Adjust the file path accordingly
file.exists(file_path)

filename <- gsub('.*"', '', file_path)  # Remove everything before the last double quote
filename <- gsub('"', '', filename)   


filename <- gsub('.*"', '', file_path)  # Remove everything before the last double quote
filename <- gsub('"', '', filename)     # Remove the double quotes

# Print the extracted filename
print(filename)

filename <- gsub('.*"', '', file_path)  # Remove everything before the last double quote
filename <- gsub('"', '', filename) 

pattaro_data_2025 <- read.csv("/Users/artemisbriasouli/formatted_round3meta_CKD_overall_IV_2GC_b36_MAFget005_Nget50_20120725_b37.csv")

final_bmi_tissue_partitioned_data <- read.csv("/Users/artemisbriasouli/my_data.csv")

bmi_tissue_data <- read.csv("/Users/artemisbriasouli/Desktop/my_data.csv")

# Example values for beta and standard error
b <- -0.03 	
	
		     # Beta coefficient
se <- 0.01 	
 # Standard error

# Calculate the Odds Ratio (OR)
OR <- exp(b)

# Calculate the 95% Confidence Interval for the beta
CI_beta_lower <- b - 1.96 * se
CI_beta_upper <- b + 1.96 * se

# Convert the beta confidence interval to OR confidence interval
CI_OR_lower <- exp(CI_beta_lower)
CI_OR_upper <- exp(CI_beta_upper)

# Display the results
cat("Odds Ratio (OR):", OR, "\n")
cat("95% Confidence Interval for OR: [", CI_OR_lower, ",", CI_OR_upper, "]\n")

# Extracting beta and standard error from the MR result
beta <- -0.02244497   # effect estimate (beta)
se <- 0.01660258   # standard error of the estimate

# Calculating the 95% confidence interval
z_value <- 1.96  # for 95% confidence
lower_ci <- beta - z_value * se
upper_ci <- beta + z_value * se

# Displaying the results
cat("Beta (Causal Effect Estimate):", beta, "\n")
cat("95% Confidence Interval: (", lower_ci, ",", upper_ci, ")\n")




filtered_snps <- gwas_data %>%
        filter(PPA4_adipose > 0.8)

write.csv(above_threshold, "above_threshold.csv", row.names = FALSE)
write.csv(above_thereshold_1, "above_threshold_1.csv", row.names = FALSE)

above_threshold <- read.csv("above_threshold.csv")
above_threshold_1<- read.csv("above_threshold_1.csv")

columns_to_remove <- c("PPA4_brain", "PPA4_adipose", "Gene_brain", "Gene_adipose")
above_threshold <- above_threshold[, !(names(above_threshold) %in% columns_to_remove)]
below_threshold <- above_threshold_1[, !(names(above_threshold_1) %in% columns_to_remove)]

# Save the updated data into two separate CSV files
write.csv(above_threshold, "above_threshold.csv", row.names = FALSE)
write.csv(below_threshold, "below_threshold.csv", row.names = FALSE)


common_names <- intersect(names(new_data), names(new_data_1))
new_data_1 <- new_data_1[, names(new_data_1) %in% common_names]
stacked_data <- rbind(new_data, new_data_1)

# Check the number of columns in each dataset
ncol_dataset1 <- ncol(new_data)
ncol_dataset2 <- ncol(new_data_1)

# Print the number of columns
print(ncol_dataset1)
print(ncol_dataset2)

library(dplyr)

# Assuming your stacked_data dataframe contains the combined dataset
# Merge the two last columns, replacing NA values
stacked_data <- input %>%
        mutate(beta = coalesce(beta_adipose, beta_brain)) %>%
        select(-beta_adipose, -beta_brain)



# Write the data_wood_mr dataframe to a text file
write.table(pattaro_data_2025, file = "/Users/artemisbriasouli/pattaro_data_2025.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



bmi_wood_final <- read_exposure_data(
        "/Users/artemisbriasouli/Not_iCloud/data_wood_mr.txt",
        clump = FALSE,
        sep = "\t",
        phenotype_col = "",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        eaf_col = "Freq_Tested_Allele_in_HRS",
        effect_allele_col = "Tested_Allele",
        other_allele_col = "Other_Allele",
        pval_col = "P",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "N",
        gene_col = "",
        id_col = "id.exposure",
        min_pval = 1e-200,
        log_pval = FALSE,
        chr_col = "chr",
        pos_col = "pos"
)


decline_outcome_brain <- read_outcome_data(
        "/Users/artemisbriasouli/data_decline_rapid.txt",
                snps = filtered_snps_brain$SNP,
                sep = "\t",
                phenotype_col = "",
               snp_col = "RSID",
        beta_col = "OR",
               se_col = "StdErr",
               eaf_col = "Freq1",
               effect_allele_col = "Allele1",
               other_allele_col = "Allele2",
               pval_col = "P.value",
       units_col = "",
         ncase_col = "",
                ncontrol_col = "",
                samplesize_col = "N_total_sum",
                gene_col = "",
                id_col = "",
                min_pval = 1e-200,
                log_pval = FALSE,
                chr_col = "",
                pos_col = ""
        )

write.csv(stacked_data, file = "my_data.csv", row.names = FALSE)

mvmr_exposures <- mv_extract_exposures_local(
        data,
        sep = ",",
        phenotype_col = "Phenotype",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval",
        id_col = "ID",
        min_pval = 1e-200,
        log_pval = FALSE,
        harmonise_strictness = 1
)

id_exposure <- c("ieu-a-299", "ieu-a-300", "ieu-a-302")
id_outcome <- "ieu-a-7"

id_df <- data.frame(ID = stacked_data$ID)

mv_exposure_dat <- mv_extract_exposures(ID)

mv_outcome_dat <- extract_outcome_data(mvmr_exposures$SNP, 'ieu-a-7')
mvdat <- mv_harmonise_data(mvmr_exposures, mv_outcome_dat, harmonise_strictness = 1)

stacked_data$Phenotype <- rep(c("adipose", "brain"), length.out = nrow(stacked_data))


id_exposure <- c("ieu-a-299", "ieu-a-300", "ieu-a-302")
id_outcome <- "ieu-a-7"

stacked_data <- stacked_data %>%
        select(-samplesize)

id_outcome <- "ieu-a-7"

mv_outcome_dat <- extract_outcome_data(stacked_data$SNP, id_outcome)

write.csv(stacked_data, file = "/Users/artemisbriasouli/stacked_data.csv", row.names = FALSE)

write.csv(kidney_data, file = "/Users/artemisbriasouli/kidney_data.csv", row.names = FALSE)

kidney_data <- read.table("/Users/artemisbriasouli/Downloads/metal_eGFR_meta_ea_6_7_NK_20230328_1.TBL.map.annot.gc", header = TRUE, sep = "\t")

colnames(kidney_data)[colnames(kidney_data) == "Allele1"] <- "effect_allele.outcome"

male_data <- fread("/Users/artemisbriasouli/Downloads/egfr_creat_male_REGENIE.imputed.finalResults.txt", header = TRUE, sep = "\t")



# Define the specific BETA and SE values
beta <- -0.03 	

	
se <- 	0.01 

# Calculate OR
OR <- exp(BETA)

# Calculate 95% confidence interval for OR
CI_lower <- exp(BETA - 1.96 * SE)
CI_upper <- exp(BETA + 1.96 * SE)

# Print results
cat("OR:", OR, "\n")
cat("95% CI:", CI_lower, "-", CI_upper, "\n")

# Define specific BETA and SE values
BETA <- -0.05179975	
SE <- 0.007331999

# Calculate 95% CI for BETA
CI_lower <- BETA - 1.96 * SE
CI_upper <- BETA + 1.96 * SE

# Print results
cat("95% CI for BETA:", CI_lower, "-", CI_upper, "\n")




options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
dis_dat <- extract_outcome_data(snps = data$SNP, outcomes = kidney_data)
dis_dat <- extract_outcome_data(data$SNP, kidney_data)

colnames(kidney_data)[colnames(kidney_data) == "RSID"] <- "SNP"

egfr_dec_data <- read.table("/Users/artemisbriasouli/Desktop/CKDGen_eGFR-decline_overall_adjDM.txt", header = TRUE, sep = "\t")
kidney_data <- subset(kidney_data, select = -c(1, 9, 10, 13, 14))

head(outcome_dat)

mv_bun_dat <- mv_harmonise_data(final_bmi_tissue_partitioned_data, ckd_new_outcome_dat, harmonise_strictness = 1)
VIEWres_cys


res_single <- mr_singlesnp(mv_bun_dat)
p2 <- mr_forest_plot(res_single)


decline_data <- read.table(gzfile("/Users/artemisbriasouli/Downloads/CKDGen_rapid3_overall.txt.gz"), header = TRUE, sep = "\t")
/Users/artemisbriasouli/Downloads/CKDGen_rapid3_overall.txt.gz


mvdat_dec <- mv_harmonise_data(final_bmi_tissue_partitioned_data, cys_outcome_dat, harmonise_strictness = 1)
)                       )
res_bun

res_new_bun <- mv_basic(mvdat_regenie)
res_new_bun

res_new_cys <- mv_residual(mv_cys_dat)
res_new_cys


forest_plot_1_to_many(
        mr_res = "res_new_cys",
        b = "b",
        se = "se",
        TraitM = "outcome",
        col1_width = 1,
        col1_title = "",
        exponentiate = FALSE,
        trans = "identity",
        ao_slc = TRUE,
        lo = NULL,
        up = NULL,
        by = NULL,
        xlab = "Effect (95% confidence interval)",
        addcols = NULL,
        addcol_widths = NULL,
        addcol_titles = "",
        subheading_size = 6,
        shape_points = 15,
        colour_scheme = "black",
        col_text_size = 5,
        weight = NULL
)

mvdat_regenie <- mv_harmonise_data(final_bmi_tissue_partitioned_data, female_regienie_outcome_dat, harmonise_strictness = 1)
generate_odds_ratios(res)
kidney_data <- mutate(kidney_data, outcome = "kidneydecline")




outcon
#Error in r$status_code : $ operator is invalid for atomic vectors

outcome_dat <- extract_outcome_data(data$SNP, outcomes= c('/Users/artemisbriasouli/kidney_data.txt'))

outcome_dat <- read_outcome_data(
                 filename='/Users/artemisbriasouli/kidney_data.txt',
                 snps = data$SNP,
                 sep = "\t",
                 phenotype_col = "",
                 snp_col = "RSID",
                 beta_col = "Effect",
                 se_col = "StdErr",
                 eaf_col = "Freq1",
                 effect_allele_col = "Allele1",
                 other_allele_col = "Allele2",
                 pval_col = "P.value",
                 units_col = "",
                 ncase_col = "",
                 ncontrol_col = "",
                 samplesize_col = "n",
                 gene_col = "",
                 id_col = "",
                 min_pval = 1e-200,
                 log_pval = FALSE,
                 chr_col = "chr",
                 pos_col = "pos")

kidney_data <- read.csv("/Users/artemisbriasouli/kidney_data.csv", stringsAsFactors = FALSE)
write.table(kidney_data, file = "/Users/artemisbriasouli/kidney_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(data, file = "/Users/artemisbriasouli/stacked_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(data, file = "stacked_data.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
data <- read.csv("/Users/artemisbriasouli/stacked_data.csv", stringsAsFactors = FALSE)

target_snp <- "rs7540681"
snp_exists <- target_snp %in% kidney_data$SNP

if (snp_exists) {
            print(paste("SNP", target_snp, "exists in the kidney dataset."))
        } else {
             print(paste("SNP", target_snp, "does not exist in the kidney dataset."))
        }
#[1] "SNP rs7540681 does not exist in the kidney dataset."

cyss_data <- read.table("/Users/artemisbriasouli/metal_bun_meta_all1.TBL.map.annot.gc", header = TRUE, sep = "\t")
bunn_data <- read.table("/Users/artemisbriasouli/metal_eGFRcys_meta1.TBL.map.annot.gc", header = TRUE, sep = "\t")
head(bun_data)

write.table(bun_data, "/Users/artemisbriasouli/bun.txt", sep = "\t", quote = FALSE, row.names = FALSE)

what_data <- read.table("/Users/artemisbriasouli/updated_bmi_tissue_partitioned.txt", header = TRUE, sep = "\t")
bmi_dec <- read.table("/Users/artemisbriasouli/bmi.effs.txt", sep = "\t", header = TRUE)
female_outcome_dat_twosample <- read_outcome_data(
        filename= "/Users/artemisbriasouli/Downloads/output_file_artemis_femme.txt",
        snps = bmi_dec$SNP,
        sep = "\t",
        phenotype_col = "ID",
        snp_col = "RSID",
        beta_col = "BETA",
        se_col = "SEBETA_GC",
        eaf_col = "MAF",
        effect_allele_col = "REF",
        other_allele_col = "ALT",
        pval_col = "PVALUE_GC",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "",
        gene_col = "",
        id_col = "",
        min_pval = 1e-200,
        log_pval = FALSE,
        chr_col = "CHROM",
        pos_col = "POS")

library(ggplot2)
library(dplyr)
library(forcats) 

# Load required libraries
library(ggforestplot)
library(tidyverse)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(forcats)  # For ordering factors

# Define color mapping for exposures
exposure_colors <- c("Brain Tissue" = "#1f78b4",  # Blue
                     "Adipose Tissue" = "#33a02c", # Green
                     "BMI Total" = "#e31a1c")  # Red (if used)

# Function to create a forest plot
plot_forest <- function(df, title, log_scale = FALSE) {
        ggplot(df, aes(x = fct_reorder(Outcome, Beta_OR), 
                       y = Beta_OR, 
                       ymin = Lower_CI, 
                       ymax = Upper_CI, 
                       color = Category)) +
                geom_pointrange(size = 1.2) +  # Point estimate and CI
                geom_hline(yintercept = ifelse(log_scale, 1, 0), 
                           linetype = "dashed", color = "black") +
                coord_flip() +  # Flip to make it a vertical forest plot
                scale_color_manual(values = exposure_colors) +
                labs(title = title,
                     x = "Outcome",
                     y = ifelse(log_scale, "Odds Ratio (95% CI)", "Beta (95% CI)"),
                     color = "Exposure") +
                theme_minimal(base_size = 14)
}

# Univariable MR Data
df_univariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "BMI_brain", 1.25, 1.02, 1.54, "Brain Tissue", TRUE,  # OR
        "CKD (Pattaro)", "BMI_adipose", 1.26, 0.96, 1.66, "Adipose Tissue", TRUE,  # OR
        "CKD (Wuttke)", "BMI_brain", 1.00, 0.99, 1.02, "Brain Tissue", TRUE,  # OR
        "CKD (Wuttke)", "BMI_adipose", 0.99, 0.98, 1.01, "Adipose Tissue", TRUE,  # OR
        "eGFR(decline)", "BMI_brain", 1.16, 1.02, 1.32, "Brain Tissue", FALSE,
        "eGFR(decline)", "BMI_adipose", 1.37, 1.16, 1.62, "Adipose Tissue", FALSE,
        "eGFRcys", "BMI_brain", -0.05, -0.06, -0.04, "Brain Tissue", FALSE,
        "eGFRcys", "BMI_adipose", -0.05, -0.07, -0.03, "Adipose Tissue", FALSE
)

# Multivariable MR Data
df_multivariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "BMI_brain", 1.19, 0.81, 1.73, "Brain Tissue", TRUE,  # OR
        "CKD (Pattaro)", "BMI_adipose", 1.12, 0.71, 1.75, "Adipose Tissue", TRUE,  # OR
        "CKD (Wuttke)", "BMI_brain", 1.18, 0.90, 1.55, "Brain Tissue", TRUE,  # OR
        "CKD (Wuttke)", "BMI_adipose", 1.20, 0.87, 1.70, "Adipose Tissue", TRUE,  # OR
        "eGFR(decline)", "BMI_brain", 1.03, 0.93, 1.2, "Brain Tissue", FALSE,
        "eGFR(decline)", "BMI_adipose", 1.07, 0.94, 1.2, "Adipose Tissue", FALSE,
        "eGFRcys", "BMI_brain", -0.04, -0.07, -0.01, "Brain Tissue", FALSE,
        "eGFRcys", "BMI_adipose", -0.02, -0.06, 0.01, "Adipose Tissue", FALSE
)

# Generate the Univariable MR Forest Plot (Beta)
plot_forest(
        df_univariable %>% filter(Log_Scale == FALSE),
        title = "Univariable MR: Beta Estimates",
        log_scale = FALSE
)

# Generate the Univariable MR Forest Plot (OR)
plot_forest(
        df_univariable %>% filter(Log_Scale == TRUE),
        title = "Univariable MR: Odds Ratios",
        log_scale = TRUE
)

# Generate the Multivariable MR Forest Plot (Beta)
plot_forest(
        df_multivariable %>% filter(Log_Scale == FALSE),
        title = "Multivariable MR: Beta Estimates",
        log_scale = FALSE
)

# Generate the Multivariable MR Forest Plot (OR)
plot_forest(
        df_multivariable %>% filter(Log_Scale == TRUE),
        title = "Multivariable MR: Odds Ratios",
        log_scale = TRUE
)


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define color mapping for exposures
exposure_colors <- c("Brain BMI" = "#1f78b4",  # Blue
                     "Adipose BMI" = "#33a02c", # Green
                     "Total BMI" = "#e31a1c")  # Red

# Univariable MR Data (Example - Fill in with your full dataset)
df_univariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", 1.25, 1.02, 1.54, "Brain BMI", TRUE,  # OR
        "CKD (Pattaro)", "Adipose BMI", 1.26, 0.96, 1.66, "Adipose BMI", TRUE,  # OR
        "CKD (Pattaro)", "Total BMI", 1.24, 1.2, 1.3, "Total BMI", TRUE,  # OR
        "eGFRcys", "Brain BMI", -0.05, -0.06, -0.04, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.05, -0.07, -0.03, "Adipose BMI", FALSE,
        "eGFRcys", "Total BMI", -0.05, -0.06, -0.046, "Total BMI", FALSE
)

# Multivariable MR Data (Example - Fill in with your full dataset)
df_multivariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", 1.19, 0.81, 1.73, "Brain BMI", TRUE,  # OR
        "CKD (Pattaro)", "Adipose BMI", 1.12, 0.71, 1.75, "Adipose BMI", TRUE,  # OR
        "CKD (Pattaro)", "Total BMI", 1.17, 1.1, 1.23, "Total BMI", TRUE,  # OR
        "eGFRcys", "Brain BMI", -0.04, -0.07, -0.01, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.02, -0.06, 0.01, "Adipose BMI", FALSE,
        "eGFRcys", "Total BMI", -0.05, -0.06, -0.046, "Total BMI", FALSE
)

# Function to create a single forest plot
plot_forest <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Order factors for better visualization
        
        ggplot(df, aes(x = Outcome, y = Beta_OR, ymin = Lower_CI, ymax = Upper_CI, color = Category, group = Exposure)) +
                geom_pointrange(size = 1.2, position = position_dodge(width = 0.6)) +  # Point estimate and CI
                geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Null effect for ORs
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors
                labs(title = title,
                     x = "Outcome",
                     y = "Effect Estimate (95% CI)",
                     color = "Exposure") +
                theme_minimal(base_size = 14)
}

plot_forest(df_univariable, "Univariable MR: BMI Pathways and Kidney Function")

# Generate the **Multivariable MR Forest Plot**
plot_forest(df_multivariable, "Multivariable MR: BMI Pathways and Kidney Function")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define color mapping for exposures
exposure_colors <- c("Brain BMI" = "#1f78b4",  # Blue
                     "Adipose BMI" = "#33a02c", # Green
                     "Total BMI" = "#e31a1c")  # Red

# Univariable MR Data (Example - Fill in with your full dataset)
df_univariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", 1.25, 1.02, 1.54, "Brain BMI", TRUE,  # OR
        



dat_femme_sex_twosample <- harmonise_data(bmi_dec, female_outcome_dat_twosample, action = 1)

# Assuming your DataFrame is called df
bmi_dec$exposure <- "femme"
bmi_dec$id.exposure <- "bmifemme"



outcome_pattaro <- read_outcome_data(
        filename = "/Users/artemisbriasouli/metal_eGFR_meta_ea1.TBL.map.annot.gc",
        snps = snps_brain$SNP,
        sep = "\t",
        phenotype_col = "",
        snp_col = "RSID",
        beta_col = "Effect",
        se_col = "StdErr",
        eaf_col = "Freq1",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        pval_col = "P.value",
        units_col = "",
        ncase_col = "",
        ncontrol_col = "",
        samplesize_col = "n"
        )

dat_regenie <- harmonise_data(bmi_dec, outcome_pattaro, action = 1)

mr <- mr(dat_regenie)

heterogeneity <- mr_heterogeneity(dat_regenie)
print(heterogeneity)

egger_intercept <- mr_pleiotropy_test(dat_regenie)
print(egger_intercept)


dat_bmi_kidney <- harmonise_data(bmi_height, outcome_transethnic_bmi_height, action = 1)

pattaro_file <- read.csv("/Users/artemisbriasouli/pattaro_gwas.txt")

# For binary outcomes (log OR scale)
beta <- 0.005967454																			
se <-  0.005035531
        
# Example for beta outcome
beta <- 0																							
se <- 0.001402239

lower_CI <- beta - 1.96 * se
upper_CI <- beta + 1.96 * se

cat("Beta:", round(beta, 3), "\n")
cat("95% CI:", round(lower_CI, 3), "-", round(upper_CI, 3), "\n")




      calculate_or_ci <- function(beta, se) {
        OR <- exp(beta)
        CI_lower <- exp(beta - 1.96 * se)
        CI_upper <- exp(beta + 1.96 * se)
        return(data.frame(OR = OR, CI_Lower = CI_lower, CI_Upper = CI_upper))
}

# For continuous outcomes
calculate_beta_ci <- function(beta, se) {
        CI_lower <- beta - 1.96 * se
        CI_upper <- beta + 1.96 * se
        return(data.frame(Beta = beta, CI_Lower = CI_lower, CI_Upper = CI_upper))
}




#  Read dataset
df <- fread("/Users/artemisbriasouli/Desktop/processed_CKD_NEW_data.txt", header = TRUE)

# Check if data is read correctly
if(nrow(df) == 0) {
        stop("The dataset is empty or not read correctly. Please check the file path and content.")
}

# Remove first row
df <- df[-1, ]

# Check if there is still data left after removing the first row
if(nrow(df) == 0) {
        stop("No data left after removing the first row. Please check the dataset.")
}

#  Remove the last column
df <- df[, -ncol(df)]

# Save cleaned dataset
write.table(df, "/Users/artemisbriasouli/Desktop/cleaned_dataset.csv", row.names = FALSE)

library(data.table)

# Step 1: Read the dataset
df <- fread("/Users/artemisbriasouli/Desktop/ckd_new.txt", header = TRUE, fill = TRUE)

# Check if data is read correctly
if(nrow(df) == 0 || ncol(df) == 0) {
        stop("The dataset is empty or not read correctly. Please check the file path and content.")
}

# Step 2: Verify Column Names
if(!"outcome" %in% colnames(df)) {
        stop("Column 'outcome' does not exist in the dataset.")
}


library(readxl)

# Load Excel file
df <- read_excel("/Users/artemisbriasouli/METABOLIC_AND_MENTAL_GENE_ONTOLOGY.xlsx")

# Keep rows where either column 8 or 9 is TRUE
filtered_df <- df[df[[8]] == TRUE | df[[9]] == TRUE, ]

# Write to new CSV
write.csv(filtered_df, "/Users/artemisbriasouli/filtered_metabolic_and_mental.csv", row.names = FALSE)


# Step 3: Review Replacement Process (if applicable)
# Replace the "outcome" column with a single value
# For example:
df$outcome <- "your_value"

# Save the cleaned dataset
write.table(df, "cleaned_dataset.txt", row.names = FALSE, sep = "\t", quote = FALSE)

View(bmi_cojo)

bmi_pattaro <- read_exposure_data(filename = "/Users/artemisbriasouli/Desktop/Meta-analysis_Wood_et_al+UKBiobank_2018.txt",  clump = FALSE,
                               sep = "\t",
                               phenotype_col = "",
                               snp_col = "SNP",
                               beta_col = "BETA",
                               se_col = "SE",
                               eaf_col = "Freq_Tested_Allele_in_HRS",
                               effect_allele_col = "Tested_Allele",
                               other_allele_col = "Other_Allele",
                               pval_col = "P",
                               units_col = "",
                               ncase_col = "",
                               ncontrol_col = "",
                               samplesize_col = "N",
                               gene_col = "",
                               id_col = "",
                               min_pval = 1e-200,
                               log_pval = FALSE,
                               chr_col = "CHR",
                               pos_col = "POS"
)

# Calculate the F-statistic for each SNP
dat_pattaro$F_statistic <- (dat_pattaro$beta.exposure / dat_pattaro$se.exposure)^2

# Check the first few F-statistics to ensure correctness
head(dat_pattaro$F_statistic)

# Count the number of strong instruments (F-statistic > 10)
num_strong_instruments <- sum(dat_pattaro$F_statistic > 10)

# Print the result
print(paste("Number of strong instruments (F-statistic > 10):", num_strong_instruments))


mr_pattaro_hetero <- mr_heterogeneity(dat_pattaro)
View(mr_pattaro_hetero)
mr_pattaro_pleiotropy <- mr_pleiotropy_test(dat_pattaro)
View(mr_pattaro_pleiotropy)

#Steiger
# Assume your data is in a dataframe called dat_pattaro
# and it contains columns: beta_exposure, beta_outcome, se_exposure, se_outcome

# Calculate R^2 for exposure and outcome
dat_pattaro$R2_exposure <- (dat_pattaro$beta.exposure^2) / 
        (dat_pattaro$beta.exposure^2 + dat_pattaro$se.exposure^2)

dat_pattaro$R2_outcome <- (dat_pattaro$beta.outcome^2) / 
        (dat_pattaro$beta.outcome^2 + dat_pattaro$se.outcome^2)

# Apply Steiger filtering: retain SNPs where R2 for exposure > R2 for outcome
dat_pattaro_steiger_filtered <- subset(dat_pattaro, R2_exposure > R2_outcome)

# Check how many SNPs are retained after filtering
num_retained_snps <- nrow(dat_pattaro_steiger_filtered)
print(paste("Number of SNPs retained after Steiger filtering:", num_retained_snps))

# Optionally, save the filtered dataset to a new file
write.table(dat_pattaro_steiger_filtered, file = "dat_pattaro_steiger_filtered.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# Install and load the mr.raps package if you haven't already
install.packages("mr.raps")
library(mr.raps)

heterogeneity <- mr_heterogeneity(dat_metab)
print(heterogeneity)


# Use the relevant columns from your dat_pattaro dataframe
beta_exposure <- dat_pattaro$beta.exposure
beta_outcome <- dat_pattaro$beta.outcome
se_exposure <- dat_pattaro$se.exposure
se_outcome <- dat_pattaro$se.outcome

# Perform MR-RAPS analysis using your dataframe
mr_raps_result <- mr.raps(beta.1 = beta_exposure, 
                          beta.2 = beta_outcome, 
                          se_1 = se_exposure, 
                          se.2 = se_outcome)

# Print the results
print(mr_raps_result)

# Assuming your dataframe is named dat_pattaro and has columns 'PVAL.EXPOSURE' and 'PVAL.OUTCOME'

# Filter the dataframe based on the p-value threshold for both exposure and outcome
filtered_dat_pattaro <- dat_pattaro[dat_pattaro$pval.exposure < 5e-8 & dat_pattaro$pval.outcome < 5e-8, ]

# View the filtered dataframe
print(filtered_dat_pattaro)


write.table(wood, file = "/Users/artemisbriasouli/Desktop/Meta-analysis_Wood_et_al+UKBiobank_2018.txt", sep = "\t", quote = FALSE, row.names = FALSE)

bmi_wood_gwas <- read.table("/Users/artemisbriasouli/Desktop/Meta-analysis_Wood_et_al+UKBiobank_2018.txt", header = TRUE, sep = "\t")

pattaro_data <- read.csv("/Users/artemisbriasouli/Desktop/formatted_round3meta_CKD_overall_IV_2GC_b36_MAFget005_Nget50_20120725_b37.txt", header = TRUE)

output_filename <- "processed_CKD_data.txt"
write.table(ckd_data, file = "/Users/artemisbriasouli/Desktop/processed_CKD_NEW_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)


library(dplyr)
library(tidyr)


summary_stats <- read.table("/Users/artemisbriasouli/Downloads/output_file_femme.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE)
library(data.table)

ref_1kg <- read.table("/Users/artemisbriasouli/Downloads/1kg.h234.rsid.txt", header = FALSE, sep = " ", stringsAsFactors = FALSE, col.names = c("CHR_POS_ALLELES", "RSID"))

ref_1kg <- ref_1kg %>%
        separate(CHR_POS_ALLELES, into = c("CHR_POS", "ALLELES"), sep = "_") %>%
        separate(CHR_POS, into = c("CHR", "POS"), sep = ":", convert = TRUE) %>%
        select(CHR, POS, RSID)  # Keep only relevant columns

ref_1kg <- fread("/Users/artemisbriasouli/Downloads/1kg.h234.rsid.txt")
colnames(summary_stats)[colnames(summary_stats) %in% c("chromosome", "position")] <- c("CHR", "POS")

#Merge summary statistics with 1KG RSID reference
merged_data <- left_join(summary_stats, ref_1kg, by = c("CHR", "POS"))

write.table(merged_data, "merged_summary_stats_with_rsid.txt", row.names = FALSE, sep = "\t", quote = FALSE)



# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define color mapping for exposures
exposure_colors <- c("Brain BMI" = "#1f78b4",  # Blue
                     "Adipose BMI" = "#33a02c", # Green
                     "Total BMI" = "#e31a1c")  # Red

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define color mapping for exposures
exposure_colors <- c("Brain BMI" = "#1f78b4",  # Blue
                     "Adipose BMI" = "#33a02c", # Green
                     "Total BMI" = "#e31a1c")  # Red (only in Univariable MR)

# Univariable MR Data (Includes **Total BMI**)
df_univariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", log(1.25), log(1.02), log(1.54), "Brain BMI", TRUE,  # Log(OR)
        "CKD (Pattaro)", "Adipose BMI", log(1.26), log(0.96), log(1.66), "Adipose BMI", TRUE,  # Log(OR)
        "CKD (Pattaro)", "Total BMI", log(1.24), log(1.2), log(1.3), "Total BMI", TRUE,  # Log(OR)
        "CKD (Wuttke)", "Brain BMI", log(1.00), log(0.99), log(1.02), "Brain BMI", TRUE,
        "CKD (Wuttke)", "Adipose BMI", log(0.99), log(0.98), log(1.01), "Adipose BMI", TRUE,
        "CKD (Wuttke)", "Total BMI", log(1.17), log(1.1), log(1.23), "Total BMI", TRUE,
        "eGFRcys", "Brain BMI", -0.05, -0.06, -0.04, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.05, -0.07, -0.03, "Adipose BMI", FALSE,
        "eGFRcys", "Total BMI", -0.05, -0.06, -0.046, "Total BMI", FALSE,
        "eGFRcrea", "Brain BMI", -0.006, -0.02, 0.006, "Brain BMI", FALSE,
        "eGFRcrea", "Adipose BMI", -0.008, -0.02, 0.008, "Adipose BMI", FALSE,
        "eGFRcrea", "Total BMI", -0.006, -0.009, -0.002, "Total BMI", FALSE,
        "BUN", "Brain BMI", 0.005, -0.006, 0.02, "Brain BMI", FALSE,
        "BUN", "Adipose BMI", 0.02, 0.002, 0.03, "Adipose BMI", FALSE,
        "BUN", "Total BMI", 0.02, 0.01, 0.021, "Total BMI", FALSE,
        "eGFR(decline)", "Brain BMI", 1.16, 1.02, 1.32, "Brain BMI", FALSE,
        "eGFR(decline)", "Adipose BMI", 1.37, 1.16, 1.62, "Adipose BMI", FALSE,
        "eGFR(decline)", "Total BMI", 1.24, 1.1, 1.38, "Total BMI", FALSE
)

# Multivariable MR Data (EXCLUDES **Total BMI**)
df_multivariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~Lower_CI, ~Upper_CI, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", log(1.19), log(0.81), log(1.73), "Brain BMI", TRUE,
        "CKD (Pattaro)", "Adipose BMI", log(1.12), log(0.71), log(1.75), "Adipose BMI", TRUE,
        "CKD (Wuttke)", "Brain BMI", log(1.18), log(0.90), log(1.55), "Brain BMI", TRUE,
        "CKD (Wuttke)", "Adipose BMI", log(1.20), log(0.87), log(1.70), "Adipose BMI", TRUE,
        "eGFRcys", "Brain BMI", -0.04, -0.07, -0.01, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.02, -0.06, 0.01, "Adipose BMI", FALSE,
        "eGFRcrea", "Brain BMI", 0.001, -0.02, 0.02, "Brain BMI", FALSE,
        "eGFRcrea", "Adipose BMI", -0.01, -0.04, 0.01, "Adipose BMI", FALSE,
        "BUN", "Brain BMI", -0.005, -0.03, 0.01, "Brain BMI", FALSE,
        "BUN", "Adipose BMI", 0.02, -0.004, 0.04, "Adipose BMI", FALSE,
        "eGFR(decline)", "Brain BMI", 1.03, 0.93, 1.2, "Brain BMI", FALSE,
        "eGFR(decline)", "Adipose BMI", 1.07, 0.94, 1.2, "Adipose BMI", FALSE
)

# Function to create the forest plot
plot_forest <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Ensure correct order
        
        ggplot(df, aes(x = Outcome, y = Beta_OR, ymin = Lower_CI, ymax = Upper_CI, color = Category, group = Exposure)) +
                geom_pointrange(size = 0.8, position = position_dodge(width = 0.6)) +  # Make dots smaller & thinner
                geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.1, position = position_dodge(0.6)) +  # Ensure CIs are visible
                geom_hline(data = df %>% filter(Log_Scale == TRUE), aes(yintercept = 0), linetype = "dashed", color = "black") +  # Log(OR) null line
                geom_hline(data = df %>% filter(Log_Scale == FALSE), aes(yintercept = 0), linetype = "dashed", color = "black") +  # Beta null line
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors
                labs(title = title,
                     x = "Outcome",
                     y = "Effect Estimate (95% CI) [Log Odds for CKD, Beta for Others]",
                     color = "Exposure") +
                theme_minimal(base_size = 14)
}

# Generate **Univariable MR Forest Plot** (with Total BMI, log(OR) for CKD)
plot_forest(df_univariable, "Univariable MR: BMI Pathways and Kidney Function")

# Generate **Multivariable MR Forest Plot** (without Total BMI)
plot_forest(df_multivariable, "Multivariable MR: Kidney Function markers")

# Given values
beta <- 2.254274
	
      # Example beta estimate
se <- 	2.2571721  # Example standard error

# Compute 95% confidence interval
lower_CI <- beta - 1.96 * se
upper_CI <- beta + 1.96 * se

# Print result
c(lower_CI, upper_CI)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create dataframe for Univariable MR
univariable_df <- data.frame(
        Outcome = c("CKD (Pattaro)", "CKD (Wuttke)", "CKD new", "eGFR decline", "eGFRcys", "eGFRcrea", "BUN"),
        Exposure = c("BMI", "BMI_brain", "BMI_adipose", "BMI", "BMI_brain", "BMI_adipose", "BMI"),
        Beta_OR = c(1.24, 1.25, 1.26, 0.15, -0.05, -0.006, 0.02),
        SE = c(0.04, NA, NA, 0.07, 0.007, 0.002, 0.002),
        Lower_CI = c(1.2, 1.02, 0.96, 0.02, -0.06, -0.009, 0.01),
        Upper_CI = c(1.3, 1.54, 1.66, 0.29, -0.04, -0.002, 0.021)
)

# Convert CKD outcomes to log(OR)
univariable_df <- univariable_df %>%
        mutate(Beta_OR = ifelse(Outcome %in% c("CKD (Pattaro)", "CKD (Wuttke)", "CKD new"), log(Beta_OR), Beta_OR),
               Lower_CI = ifelse(Outcome %in% c("CKD (Pattaro)", "CKD (Wuttke)", "CKD new"), log(Lower_CI), Lower_CI),
               Upper_CI = ifelse(Outcome %in% c("CKD (Pattaro)", "CKD (Wuttke)", "CKD new"), log(Upper_CI), Upper_CI))

# Plot Univariable MR
ggplot(univariable_df, aes(x = Exposure, y = Beta_OR, ymin = Lower_CI, ymax = Upper_CI)) +
        geom_pointrange() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_flip() +
        theme_minimal() +
        labs(title = "Univariable MR Forest Plot", y = "Beta (or log OR)", x = "Exposure")

# Create dataframe for Multivariable MR
multivariable_df <- data.frame(
        Outcome = c("CKD (Pattaro)", "CKD (Wuttke)", "eGFR decline", "eGFRcys", "eGFRcrea", "BUN"),
        Exposure = c("BMI_brain", "BMI_adipose", "BMI_brain", "BMI_adipose", "BMI_brain", "BMI_adipose"),
        Beta_OR = c(1.19, 1.12, 0.03, -0.04, 0.001, -0.005),
        SE = c(0.23, 0.19, 0.06, 0.01, 0.01, 0.01),
        Lower_CI = c(0.81, 0.71, -0.08, -0.07, -0.02, -0.03),
        Upper_CI = c(1.73, 1.75, 0.14, -0.01, 0.02, 0.01)
)

# Convert CKD outcomes to log(OR)
multivariable_df <- multivariable_df %>%
        mutate(Beta_OR = ifelse(Outcome %in% c("CKD (Pattaro)", "CKD (Wuttke)"), log(Beta_OR), Beta_OR),
               Lower_CI = ifelse(Outcome %in% c("CKD (Pattaro)", "CKD (Wuttke)"), log(Lower_CI), Lower_CI),
               Upper_CI = ifelse(Outcome %in% c("CKD (Pattaro)", "CKD (Wuttke)"), log(Upper_CI), Upper_CI))

# Plot Multivariable MR
ggplot(multivariable_df, aes(x = Exposure, y = Beta_OR, ymin = Lower_CI, ymax = Upper_CI)) +
        geom_pointrange() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_flip() +
        theme_minimal() +
        labs(title = "Multivariable MR Forest Plot", y = "Beta (or log OR)", x = "Exposure")

calculate_CI <- function(beta, se) {
        lower_CI <- beta - (1.96 * se)
        upper_CI <- beta + (1.96 * se)
        return(c(lower_CI, upper_CI))
}

# Univariable MR Data (Includes **Total BMI**)
df_univariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~SE, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", log(1.25), NA, "Brain BMI", TRUE,  # Log(OR)
        "CKD (Pattaro)", "Adipose BMI", log(1.26), NA, "Adipose BMI", TRUE,  # Log(OR)
        "CKD (Pattaro)", "Total BMI", log(1.24), NA, "Total BMI", TRUE,  # Log(OR)
        "CKD (Wuttke)", "Brain BMI", log(1.00), NA, "Brain BMI", TRUE,
        "CKD (Wuttke)", "Adipose BMI", log(0.99), NA, "Adipose BMI", TRUE,
        "CKD (Wuttke)", "Total BMI", log(1.17), NA, "Total BMI", TRUE,
        "eGFRcys", "Brain BMI", -0.05, 0.007, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.05, 0.01, "Adipose BMI", FALSE,
        "eGFRcys", "Total BMI", -0.05, 0.002, "Total BMI", FALSE,
        "eGFRcrea", "Brain BMI", -0.005, 0.006, "Brain BMI", FALSE,
        "eGFRcrea", "Adipose BMI", -0.008, 0.008, "Adipose BMI", FALSE,
        "eGFRcrea", "Total BMI", -0.006, 0.002, "Total BMI", FALSE,
        "BUN", "Brain BMI", 0.005, 0.006, "Brain BMI", FALSE,
        "BUN", "Adipose BMI", 0.02, 0.008, "Adipose BMI", FALSE,
        "BUN", "Total BMI", 0.02, 0.002, "Total BMI", FALSE,
        "eGFR(decline)", "Brain BMI", 0.15, 0.07, "Brain BMI", FALSE,
        "eGFR(decline)", "Adipose BMI", 0.32, 0.08, "Adipose BMI", FALSE,
        "eGFR(decline)", "Total BMI", 2.25, 2.26, "Total BMI", FALSE
) %>%
        rowwise() %>%
        mutate(CI = list(calculate_CI(Beta_OR, SE))) %>%
        mutate(Lower_CI = CI[1], Upper_CI = CI[2]) %>%
        select(-CI)

# Multivariable MR Data (EXCLUDES **Total BMI**)
df_multivariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~SE, ~Category, ~Log_Scale,
        "CKD (Pattaro)", "Brain BMI", log(1.19), NA, "Brain BMI", TRUE,
        "CKD (Pattaro)", "Adipose BMI", log(1.12), NA, "Adipose BMI", TRUE,
        "CKD (Wuttke)", "Brain BMI", log(1.18), NA, "Brain BMI", TRUE,
        "CKD (Wuttke)", "Adipose BMI", log(1.20), NA, "Adipose BMI", TRUE,
        "eGFRcys", "Brain BMI", -0.04, 0.01, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.02, 0.02, "Adipose BMI", FALSE,
        "eGFRcrea", "Brain BMI", 0.001, 0.01, "Brain BMI", FALSE,
        "eGFRcrea", "Adipose BMI", -0.01, 0.01, "Adipose BMI", FALSE,
        "BUN", "Brain BMI", -0.005, 0.01, "Brain BMI", FALSE,
        "BUN", "Adipose BMI", 0.02, 0.01, "Adipose BMI", FALSE,
        "eGFR(decline)", "Brain BMI", 1.03, 0.93, 1.2, "Brain BMI", FALSE,
        "eGFR(decline)", "Adipose BMI", 1.07, 0.94, 1.2, "Adipose BMI", FALSE
) %>%
        rowwise() %>%
        mutate(CI = list(calculate_CI(Beta_OR, SE))) %>%
        mutate(Lower_CI = CI[1], Upper_CI = CI[2]) %>%
        select(-CI)

# Function to create the forest plot
plot_forest <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Ensure correct order
        
        ggplot(df, aes(x = Outcome, y = Beta_OR, ymin = Lower_CI, ymax = Upper_CI, color = Category, group = Exposure)) +
                geom_pointrange(size = 0.8, position = position_dodge(width = 0.6)) +  # Make dots smaller
                geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.1, position = position_dodge(0.6)) +  # Ensure CIs are visible
                geom_hline(data = df %>% filter(Log_Scale == TRUE), aes(yintercept = 0), linetype = "dashed", color = "black") +  # Log(OR) null line
                geom_hline(data = df %>% filter(Log_Scale == FALSE), aes(yintercept = 0), linetype = "dashed", color = "black") +  # Beta null line
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors
                labs(title = title,
                     x = "Outcome",
                     y = "Effect Estimate (95% CI) [Log Odds for CKD, Beta for Others]",
                     color = "Exposure") +
                theme_minimal(base_size = 14)
}

# Generate **Univariable MR Forest Plot** (with Total BMI, log(OR) for CKD)
plot_forest(df_univariable, "Univariable MR: tissue-stratified BMI and Kidney Function")

# Generate **Multivariable MR Forest Plot** (without Total BMI)
plot_forest(df_multivariable, "Multivariable MR: tissue-stratified BMI and Kidney Function")

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define color mapping for exposures
exposure_colors <- c("Brain BMI" = "#1f78b4",  # Blue
                     "Adipose BMI" = "#33a02c", # Green
                     "Total BMI" = "#e31a1c")  # Red (Only in Univariable MR)

# Function to calculate 95% confidence intervals
calculate_CI <- function(beta, se) {
        lower_CI <- beta - (1.96 * se)
        upper_CI <- beta + (1.96 * se)
        return(c(lower_CI, upper_CI))
}

# Univariable MR Data (Includes **Total BMI**, Omits eGFR decline)
df_univariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~SE, ~Category, ~Log_Scale,
        "eGFRcys", "Brain BMI", -0.05, 0.007, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.05, 0.01, "Adipose BMI", FALSE,
        "eGFRcys", "Total BMI", -0.05, 0.002, "Total BMI", FALSE,
        "eGFRcrea", "Brain BMI", -0.005, 0.006, "Brain BMI", FALSE,
        "eGFRcrea", "Adipose BMI", -0.008, 0.008, "Adipose BMI", FALSE,
        "eGFRcrea", "Total BMI", -0.006, 0.002, "Total BMI", FALSE,
        "BUN", "Brain BMI", 0.005, 0.006, "Brain BMI", FALSE,
        "BUN", "Adipose BMI", 0.02, 0.008, "Adipose BMI", FALSE,
        "BUN", "Total BMI", 0.02, 0.002, "Total BMI", FALSE,
        "eGFRdecline", "Brain BMI", 0.07, 0.03, "Brain BMI", FALSE,
        "eGFRdecline", "Adipose BMI", 0.08, 0.04, "Adipose BMI", FALSE,
        "eGFRdecline", "Total BMI",0.03, 0.01, "Total BMI", FALSE,
        
) %>%
        rowwise() %>%
        mutate(CI = list(calculate_CI(Beta_OR, SE))) %>%
        mutate(Lower_CI = CI[1], Upper_CI = CI[2]) %>%
        select(-CI)

# Multivariable MR Data (EXCLUDES **Total BMI**, Omits eGFR decline)
df_multivariable <- tribble(
        ~Outcome, ~Exposure, ~Beta_OR, ~SE, ~Category, ~Log_Scale,
        "eGFRcys", "Brain BMI", -0.04, 0.01, "Brain BMI", FALSE,
        "eGFRcys", "Adipose BMI", -0.02, 0.02, "Adipose BMI", FALSE,
        "eGFRcys", "Total BMI", -0.05, 0.002, "Total BMI", FALSE,
        "eGFRcrea", "Brain BMI", 0.001, 0.01, "Brain BMI", FALSE,
        "eGFRcrea", "Adipose BMI", -0.01, 0.01, "Adipose BMI", FALSE,
        "eGFRcrea", "Total BMI", -0.006, 0.002, "Total BMI", FALSE,
        "BUN", "Brain BMI", -0.005, 0.01, "Brain BMI", FALSE,
        "BUN", "Adipose BMI", 0.02, -0.004, "Adipose BMI", FALSE,
        "BUN", "Total BMI", 0.02, 0.002, "Total BMI", FALSE,
        "eGFRdecline", "Brain BMI", 0.03, 0.06, "Brain BMI", FALSE,
        "eGFRdecline", "Adipose BMI", 0.06, 0.07, "Adipose BMI", FALSE,
        "eGFRdecline", "Total BMI", 0.03, 0.01, "Total BMI", FALSE
        ) %>%
        rowwise() %>%
        mutate(CI = list(calculate_CI(Beta_OR, SE))) %>%
        mutate(Lower_CI = CI[1], Upper_CI = CI[2]) %>%
        select(-CI)

# Function to create the forest plot
plot_forest <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Ensure correct order
        
        ggplot(df, aes(x = Outcome, y = Beta_OR, ymin = Lower_CI, ymax = Upper_CI, color = Category, group = Exposure)) +
                geom_pointrange(size = 0.8, position = position_dodge(width = 0.6)) +  # Make dots smaller
                geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.1, position = position_dodge(0.6)) +  # Ensure CIs are visible
                geom_hline(data = df %>% filter(Log_Scale == TRUE), aes(yintercept = 0), linetype = "dashed", color = "black") +  # Log(OR) null line
                geom_hline(data = df %>% filter(Log_Scale == FALSE), aes(yintercept = 0), linetype = "dashed", color = "black") +  # Beta null line
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors
                labs(title = title,
                     x = "Outcome",
                     y = "Effect Estimate (95% CI)",
                     color = "Exposure") +
                theme_minimal(base_size = 14)
}

# Generate **Univariable MR Forest Plot** (with Total BMI, log(OR) for CKD)


# Generate **Multivariable MR Forest Plot** (without Total BMI)
plot_forest(df_multivariable, "Multivariable MR: Kidney Function markers")

"CKD (Pattaro)", "Brain BMI", log(1.25), NA, "Brain BMI", TRUE,
"CKD (Pattaro)", "Adipose BMI", log(1.26), NA, "Adipose BMI", TRUE,
"CKD (Pattaro)", "Total BMI", log(1.24), NA, "Total BMI", TRUE,
"CKD (Wuttke)", "Brain BMI", log(1.00), NA, "Brain BMI", TRUE,
"CKD (Wuttke)", "Adipose BMI", log(0.99), NA, "Adipose BMI", TRUE,
"CKD (Wuttke)", "Total BMI", log(1.17), NA, "Total BMI", TRUE,

"CKD (Pattaro)", "Brain BMI", log(1.19), NA, "Brain BMI", TRUE,
"CKD (Pattaro)", "Adipose BMI", log(1.12), NA, "Adipose BMI", TRUE,
"CKD (Wuttke)", "Brain BMI", log(1.18), NA, "Brain BMI", TRUE,
"CKD (Wuttke)", "Adipose BMI", log(1.20), NA, "Adipose BMI", TRUE,

# CKD and CKDnew Data (Odds Ratios and 95% CI for BMI, BMI_brain, and BMI_adipose)
df_ckd <- tribble(
        ~Outcome, ~Exposure, ~OR, ~CI_Lower, ~CI_Upper, ~Category, ~Log_Scale,
        "CKD", "BMI", 1.24, 1.2, 1.3, "BMI", FALSE,
        "CKD", "BMI_brain", 1.25, 1.02, 1.54, "BMI_brain", FALSE,
        "CKD", "BMI_adipose", 1.26, 0.96, 1.66, "BMI_adipose", FALSE,
        "CKDnew", "BMI", 1.17, 1.1, 1.23, "BMI", FALSE,
        "CKDnew", "BMI_brain", 1.00, 0.99, 1.02, "BMI_brain", FALSE,
        "CKDnew", "BMI_adipose", 0.99, 0.98, 1.01, "BMI_adipose", FALSE
)

# Function to create the forest plot for CKD and CKDnew
plot_forest_ckd <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Ensure correct order
        
        ggplot(df, aes(x = Outcome, y = OR, ymin = CI_Lower, ymax = CI_Upper, color = Category, group = Exposure)) +
                geom_pointrange(size = 0.8, position = position_dodge(width = 0.6)) +  # Make dots smaller
                geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, position = position_dodge(0.6)) +  # Ensure CIs are visible
                geom_hline(data = df %>% filter(Log_Scale == TRUE), aes(yintercept = 1), linetype = "dashed", color = "black") +  # OR null line
                geom_hline(data = df %>% filter(Log_Scale == FALSE), aes(yintercept = 1), linetype = "dashed", color = "black") +  # Beta null line
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors
                labs(title = title,
                     x = "Outcome",
                     y = "Effect Estimate (95% CI) [Odds Ratio]",
                     color = "Exposure") +
                theme_minimal(base_size = 14)
}

# Generate Forest Plot for CKD and CKDnew (based on Odds Ratios)
plot_forest_ckd(df_ckd, "Univariable MR: tissue-stratified BMI and CKD")

# CKD and CKDnew Data (Odds Ratios and 95% CI for BMI, BMI_brain, and BMI_adipose)
df_ckd <- tribble(
        ~Outcome, ~Exposure, ~OR, ~CI_Lower, ~CI_Upper, ~Category, ~Log_Scale,
        "CKDpattaro", "BMI", 1.24, 1.2, 1.3, "Total BMI", FALSE,
        "CKDpattaro", "BMI_brain", 1.25, 1.02, 1.54, "Brain BMI", FALSE,
        "CKDpattaro", "BMI_adipose", 1.26, 0.96, 1.66, "Adipose BMI", FALSE,
        "CKDwuttke", "BMI", 1.16, 1.1099, 1.23, "Total BMI", FALSE,
        "CKDwuttke", "BMI_brain", 1.29, 1.1, 1.513, "Brain BMI", FALSE,
        "CKDwuttke", "BMI_adipose", 1.313, 1.063, 1.62, "Adipose BMI", FALSE
)

# Custom colors for each exposure category
exposure_colors <- c("Total BMI" = "blue", "Brain BMI" = "green", "Adipose BMI" = "orange")

# Function to create the forest plot for CKD and CKDnew
plot_forest_ckd <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Ensure correct order
        
        ggplot(df, aes(x = Outcome, y = OR, ymin = CI_Lower, ymax = CI_Upper, color = Category, group = Exposure)) +
                geom_pointrange(size = 1.2, position = position_dodge(width = 0.6)) +  # Make points more visible
                geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, position = position_dodge(0.6)) +  # Ensure CIs are visible
                geom_hline(data = df %>% filter(Log_Scale == TRUE), aes(yintercept = 1), linetype = "dashed", color = "black") +  # OR null line
                geom_hline(data = df %>% filter(Log_Scale == FALSE), aes(yintercept = 1), linetype = "dashed", color = "black") +  # Beta null line
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors for exposures
                labs(title = title,
                     x = "Outcome",
                     y = "Odds Ratio",
                     color = "Exposure") +
                theme_minimal(base_size = 14) +
                theme(axis.text.y = element_text(size = 12),  # Adjust axis text size
                      legend.position = "right")  # Move the legend to the right
}


# Multivariable MR Data (CKD (Pattaro et al.), CKDnew (Wuttke et al.))
df_multivariable_ckd <- tribble(
        ~Outcome, ~Exposure, ~OR, ~CI_Lower, ~CI_Upper, ~Category, ~Log_Scale,
        "CKDpattaro", "BMI_brain", 1.19, 0.81, 1.73, "Brain BMI", FALSE,
        "CKDpattaro", "BMI_adipose", 1.12, 0.71, 1.75, "Adipose BMI", FALSE,
        "CKDwuttke", "BMI_brain", 1.18, 0.9, 1.55, "Brain BMI", FALSE,
        "CKDwuttke", "BMI_adipose", 1.20, 0.87, 1.7, "Adipose BMI", FALSE,
        "CKDpattaro", "Total BMI", 1.17, 1.1, 1.23, "Total BMI", FALSE,
        "CKDwuttke", "Total BMI", 1.24, 1.2, 1.3, "Total BMI", FALSE
)

# Custom colors for each exposure category
exposure_colors <- c("Total BMI" = "blue", "Brain BMI" = "green", "Adipose BMI" = "orange")

# Function to create the forest plot for CKD and CKDnew (Multivariable MR)
plot_forest_multivariable <- function(df, title) {
        df <- df %>%
                mutate(Outcome = factor(Outcome, levels = rev(unique(Outcome))))  # Ensure correct order
        
        ggplot(df, aes(x = Outcome, y = OR, ymin = CI_Lower, ymax = CI_Upper, color = Category, group = Exposure)) +
                geom_pointrange(size = 1.2, position = position_dodge(width = 0.6)) +  # Make points more visible
                geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, position = position_dodge(0.6)) +  # Ensure CIs are visible
                geom_hline(data = df %>% filter(Log_Scale == TRUE), aes(yintercept = 1), linetype = "dashed", color = "black") +  # OR null line
                geom_hline(data = df %>% filter(Log_Scale == FALSE), aes(yintercept = 1), linetype = "dashed", color = "black") +  # Beta null line
                coord_flip() +  # Flip axes for readability
                scale_color_manual(values = exposure_colors) +  # Custom colors for exposures
                labs(title = title,
                     x = "Outcome",
                     y = "Odds Ratio",
                     color = "Exposure") +
                theme_minimal(base_size = 14) +
                theme(axis.text.y = element_text(size = 12),  # Adjust axis text size
                      legend.position = "right")  # Move the legend to the right
}

# Generate Forest Plot for Multivariable MR: CKD and CKDnew (based on Odds Ratios)

bx_matrix <- as.matrix(mvdat_dec$exposure_beta)

bxse_matrix <- as.matrix(mvdat_dec$exposure_se)

Y_vector <- as.vector(mvdat_dec$outcome_beta)

SE_Y_vector <- as.vector(mvdat_dec$outcome_se)

mrmvinput <- mr_mvinput(
        bx = bx_matrix,
        bxse = bxse_matrix,
        by = Y_vector,
        byse = SE_Y_vector,
        exposure = "exposure",
        outcome = "eGFRcys",
        snps = "snp",
        effect_allele = NA,
        other_allele = NA,
        eaf = NA
)

mvivw_cys <- mr_mvivw(
        mrmvinput,
        model = "random"
)

mvivw_cys

class(mvivw_cys)
slotNames(mvivw_cys)

mvivw_cys@Estimate
mvivw_cys@StdError

mvivw_cys@Correlation

mvivw_cys@Heter.Stat

mvegger_cys <- mr_mvegger(mrmvinput)
mvegger_cys

adipose_snps <- unique(snps_adipose$SNP)
brain_snps   <- unique(snps_brain$SNP)
bmi_snps <- unique(bmi_dec$SNP)

wuttke_adipose <- pattaro_file[
        pattaro_file$rsID %in% adipose_snps,
]

wuttke_brain <- pattaro_file[
        pattaro_file$rsID %in% brain_snps,
]

wuttke_bmi <- pattaro_file[
        pattaro_file$rsID %in% bmi_snps,
]

wuttke1_adipose <- bunn_data[
        bunn_data$RSID %in% adipose_snps,
]

wuttke1_brain <- bunn_data[
        bunn_data$RSID %in% brain_snps,
]

wuttke1_bmi <- bunn_data[
        bunn_data$RSID %in% bmi_snps,
]

wuttke2_adipose <- egfr_dec_data[
        egfr_dec_data$RSID %in% adipose_snps,
]

wuttke2_brain <- egfr_dec_data[
        egfr_dec_data$RSID %in% brain_snps,
]

wuttke2_bmi <- egfr_dec_data[
        egfr_dec_data$RSID %in% bmi_snps,
]

wuttke_adipose_outcome <- format_data(
        wuttke_adipose,
        type = "outcome",
        snp_col = "rsID",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "allele1",
        other_allele_col = "allele2",
        eaf_col = "freqA1",
        pval_col = "pval",
        samplesize_col = "N",
        phenotype_col = "pattaro"
)

dat_wuttke_adipose <- harmonise_data(
        exposure_dat = snps_adipose,
        outcome_dat = wuttke_adipose
)

dat_wuttke_bmi <- harmonise_data(
        exposure_dat = bmi_dec,
        outcome_dat = wuttke_adipose_outcome
)

dat_wuttke_brain <- harmonise_data(
        exposure_dat = snps_brain,
        outcome_dat = wuttke_egfr_outcome_data
)

dat_wuttke_adipose1 <- harmonise_data(
        exposure_dat = adipose_exposure,
        outcome_dat = wuttke_adipose_outcome
)

dat_wuttke_bmi1 <- harmonise_data(
        exposure_dat = adipose_exposure,
        outcome_dat = wuttke_adipose_outcome
)

dat_wuttke_brain1 <- harmonise_data(
        exposure_dat = brain_exposure,
        outcome_dat = wuttke_brain_outcome
)