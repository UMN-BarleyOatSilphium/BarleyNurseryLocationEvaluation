## ASREML modeling

library(asreml)
library(tibble)
library(magrittr)

# Set the working directory
setwd("G:/AGRO/BARLEY_LAB/Jeff/Projects/BarleyNurseryAnalysis/")
proj_dir <- getwd()

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


# Read or load data
# Phenotypic data to model
load(file.path(result_dir, "phenotype_data_tomodel.RData"))

# Trial metadata
trial_metadata <- read.csv(file = file.path(data_dir, "nursery_trial_metadata_use1.csv"))

# Line metadata
line_metadata <- read.csv(file = file.path(data_dir, "nursery_entry_metadata_use.csv"))

# Pedigree relationship matrix
load(file = file.path(data_dir, "nursery_pedigree_relmat.RData"))



## Convert the pedigree mat to A
A <- as.matrix(pedigree_Amat[genotypes, genotypes]); rm(pedigree_Amat)

## Subset more relevant traits
traits_tokeep <- c("GrainProtein", "MaltExtract", "GrainYield", "HeadingDate", "PlantHeight")


## Iterate over rows in this data.frame
for (i in seq_len(nrow(pheno_to_model))) {

  df1 <- droplevels(pheno_to_model$data[[i]])
    # Create interaction
  df1$ge <- interaction(df1$line_name, df1$environment, sep = ":")
  
  # ## Subset data to demo
  # df1 <- subset(df1, year1 < 2008) %>%
  #   droplevels()
  
  # Subset the A matrix
  A_use <- subset(A, row.names(A) %in% levels(df1$line_name), colnames(A) %in% levels(df1$line_name))
  
  # Calculate the inverse
  Ainv <- solve(Matrix(A_use, sparse = TRUE))
  Ainv <- structure(as.matrix(Ainv), dimnames = dimnames(A_use), "INVERSE" = TRUE)
  
  ## Build a relationship matrix for GxE
  # E matrix for environments
  # E <- diag(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  E <- Diagonal(nlevels(df1$environment)); dimnames(E) <- replicate(2, levels(df1$environment), simplify = F)
  
  # E <- select(df1, line_name, environment, value) %>%
  #   spread(environment, value) %>%
  #   select(levels(df1$environment)) %>%
  #   cor(., use = "pairwise.complete.obs")
  
  # Inverse of AE is equal to kronecker(Ainv, E)
  AEinv <- kronecker(Ainv, E, make.dimnames = TRUE)
  dimnames(AE) <- replicate(2, paste(rep(row.names(A_use), each = ncol(E)), row.names(E), sep = ":"), simplify = F)
  
  
  
  
  ## Fit models
  # Fixed environment, random genotype, CS variance
  model <- asreml(fixed = value ~ 1 + environment, 
                   random = ~ line_name,
                   workspace = "1gb",
                   data = df1)
  
  # Fixed environment, random genotype (regular Amat), diag variance
  model <- asreml(fixed = value ~ 1 + environment, 
                  random = ~ vm(line_name, A_use),
                  residual = ~ dsum(~ units | environment),
                  fail = "soft",
                  workspace = "1gb",
                  data = df1)
  
  # Fixed environment, random genotype (Amat inverse), diag variance
  model <- asreml(fixed = value ~ 1 + environment, 
                  random = ~ vm(line_name, Ainv),
                  residual = ~ dsum(~ units | environment),
                  fail = "soft",
                  workspace = "1gb",
                  data = df1)
  
  # Fixed environment, random genotype + GxE, diag variance
  model <- asreml(fixed = value ~ 1 + environment, 
                  random = ~ vm(line_name, A_use) + vm(ge, AEinv),
                  # residual = ~ dsum(~ units | environment),
                  
                  workspace = "4gb",
                  data = df1)
  
  # Fixed environment and genotype, random GxE with FA covariance, diag variance
  model <- asreml(fixed = value ~ 1 + environment + line_name, 
                  random = ~ fa(ge, 1),
                  # residual = ~ dsum(~ units | environment),
                  workspace = "4gb",
                  data = df1)
  
  
  
  
  # Fixed environment, random genotype + GxE, diag variance
  model <- asreml(fixed = value ~ 1 + environment, 
                  random = ~ vm(line_name, A_use) + id(vm(line_name, A_use)):diag(environment),
                  residual = ~ dsum(~ units | environment),
                  workspace = "4gb",
                  data = df1)
  
  # View estimate variance components
  summary(model)$varcomp
  
  # Fixed environment, random genotype and GxE, CS variance
  model <- asreml(fixed = value ~ 1 + environment, 
                   random = ~ line_name + line_name:environment,
                   residual = ~ idv(units),
                   workspace = 5e8,
                   data = df1)
  
  # Fixed environment, random genotype (with relmat) and GxE, CS variance
  model <- asreml(fixed = value ~ 1 + environment, 
                  random = ~ vm(line_name, A) + line_name:environment,
                  residual = ~ units,
                  workspace = 5e8,
                  data = df1)
  
  # Fixed environment, random genotype (with relmat) and GxE, CS variance
  model <- asreml(fixed = value ~ 1 + environment, 
                  random = ~ vm(line_name, A) + line_name:environment,
                  residual = ~ diag(environment),
                  workspace = 5e8,
                  data = df1)
  
  # View estimate variance components
  summary(model)$varcomp
  
  # Calculate heritability
  vpredict(object = model, xform = h2 ~ V1 / (V1 + V3))
  
  # Retried BLUPs
  summary(model, coef = T)$coef.random
  
  # Retrieve blues
  summary(model, coef = T)$coef.fixed
  
  
  
  
  # Fit the model
  fit <- mmer(fixed = value ~ 1 + environment,
              random = ~vs(line_name, Gu = A_use) + vs(line_name:environment, Gu = AE),
              rcov = ~vs(ds(environment), units),
              data = df1, date.warning = FALSE)


























