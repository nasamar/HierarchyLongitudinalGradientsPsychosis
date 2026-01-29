################################################################################
# Script to apply COMBATLS on MIND networks 
################################################################################

#  Copyright (C) 2026 University of Seville
# 
#  Written by Natalia García San Martín (ngarcia1@us.es)
# 
#  This file is part of Hierarchy Longitudinal MIND Gradients Psychosis toolkit.
#
#  Hierarchy Longitudinal MIND Gradients Psychosis toolkit is free software: 
#  you can redistribute it and/or modify it under the terms of the 
#  GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
# 
#  Hierarchy Longitudinal MIND Gradients Psychosis toolkit is distributed in the hope that 
#  it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
#  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Hierarchy Longitudinal MIND Gradients Psychosis toolkit. If not, see 
#  <https://www.gnu.org/licenses/>.

rm(list=ls()) # Previous data cleaning in memory
gc()

library(ComBatFamily)
library(R.matlab)
library(stringr)
library(gamlss)

library(future.apply)

plan(multisession)  # Usa múltiples núcleos (multisession funciona en todos los SO)


location <- 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/'


parcellation <- 'aparc_500_sym_longitudinal'
# parcellation <- 'subcortical_longitudinal'


types <- c('COMBATLS','COMBATLS_covars')
types <- c('COMBATLS_covars')

read <- paste0('datasets/MIND/MIND_networks_PAFIP/',parcellation)

setwd(paste0(location,read))

# OPEN MIND
subjects_csv <- list.files(getwd(),pattern = "\\.csv$")
# subjects_csv <- gsub('_500.*','',subjects_csv)

n_regions <- ifelse(grepl('aparc_500',parcellation),318,
                    ifelse(parcellation == 'subcortical_longitudinal',14,
                      ifelse(parcellation == 'Yeo2011_7Networks_longitudinal',14,
                           ifelse(parcellation == 'vonEconomo_longitudinal',86,
                                  ifelse(grepl('aparc',parcellation),68)))))
              
# # MIND_vectors <- as.data.frame(matrix(NA, nrow = length(subjects_csv), ncol = n_regions*n_regions)) 
MIND_vectors_list <- vector("list", length(subjects_csv))
rownames_list <- character(length(subjects_csv))

for (i in seq_along(subjects_csv)) {
  subject_csv <- subjects_csv[i]
  MIND <- as.matrix(read.csv(subject_csv, row.names = 1))

  if (i==1){
    regions <- colnames(MIND)
  }
  if (ncol(MIND) < n_regions){
    lacking_region <- setdiff(regions,colnames(MIND))

    index_lacking_region <- which(regions==lacking_region)
    index_symmetric_region <-  which(grepl(gsub('^[^_]*_','',lacking_region),colnames(MIND)))
    MIND <- cbind(MIND[,1:index_lacking_region-1],MIND[,index_symmetric_region],MIND[,(index_lacking_region+1):n_regions-1])
    MIND <- rbind(MIND[1:index_lacking_region-1,],MIND[index_symmetric_region,],MIND[(index_lacking_region+1):n_regions-1,])
  }
  # MIND_vector <- as.vector(as.matrix(MIND))
  # MIND_vectors[i,] <- MIND_vector
  MIND_vectors_list[[i]] <- as.vector(MIND)

  if (grepl('longitudinal',parcellation) | parcellation=='subcortical_longitudinal'){
    rownames_list[i] <- paste0(str_extract(subject_csv, "(?<=sub-)[0-9]+"), "_", str_extract(subject_csv, "(?<=ses-)[0-9]+"))
  }else{
    rownames_list[i] <- gsub('_500.*','',subject_csv,gsub('.csv','',subject_csv))
  }

  if (i %% 100 == 0) gc()
}

MIND_vectors <- as.data.frame(do.call(rbind, MIND_vectors_list))
rownames(MIND_vectors) <- rownames_list
rm(MIND_vectors_list); gc()

covariates <- read.csv(paste0(location,'datasets/PAFIP/Covariates_complete.csv'),row.names=1)
covariates <- covariates[rownames(covariates)%in%rownames_list,]
if (parcellation == 'subcortical_longitudinal'){
  excluded_subjects <- read.csv(paste0(location,'datasets/PAFIP/ENIGMA_Shape/QA_Status.csv'),row.names=1)
  excluded_subjects <- rownames(excluded_subjects)[apply(excluded_subjects, 1, function(x) any(x == 0))]
  excluded_subjects <- paste0(str_extract(excluded_subjects, "(?<=sub-)[0-9]+"), "_", str_extract(excluded_subjects, "(?<=ses-)[0-9]+"))
  covariates <- covariates[rownames(covariates)[!rownames(covariates)%in%excluded_subjects],]
  
}

# covariates <- covariates[covariates$Assessment==1,]

Euler <- read.csv(paste0(location,'datasets/PAFIP/euler_longitudinal.csv'),row.names=1)
rownames(Euler) <- gsub('ses-','',gsub('sub-','',rownames(Euler)))
Euler <- Euler[intersect(rownames(covariates), rownames(Euler)),]

covariates <- covariates[intersect(rownames(covariates), rownames(Euler)),]

covariates$Machine_Teslas <- as.factor(covariates$Machine_Teslas)


if (grepl('longitudinal',parcellation) & !grepl('subcortical',parcellation)){
  global <- read.csv(paste0(location,'datasets/PAFIP/aparc_formatted_1293s_fs_v7.4.1_long_FUP_13-01-2025/global_aseg_passedqc.csv'),row.names=1) 
  rownames(global) <- gsub('ses-','',gsub('sub-','',row.names(global)))
  
}


MIND_vectors <- MIND_vectors[rownames(MIND_vectors)%in%rownames(covariates),]

index_zeros = which(apply(MIND_vectors, 2, function(col) all(col == 0))) # quitar columnas == 0

covariates <- covariates[rownames(MIND_vectors),]
if (parcellation=='aparc_longitudinal'){
  global <- global[rownames(covariates),]
}

for (type in types){ 
  
  formula <- y ~ Sex+Case+Sex*Case+Age_MRI
  # formula <- y ~ Sex+Case+Age_MRI

  formula_eTIV <- y ~ Sex
  
  # MIND
  if (grepl('500',parcellation)){

    MIND_vectors_split <- MIND_vectors[,-index_zeros]
    MIND_vectors_split_1 <- MIND_vectors_split[,c(1:(ncol(MIND_vectors_split)/3))]
    MIND_vectors_split_2 <- MIND_vectors_split[,c((ncol(MIND_vectors_split)/3+1):(2*ncol(MIND_vectors_split)/3))]
    MIND_vectors_split_3 <- MIND_vectors_split[,c((2*ncol(MIND_vectors_split)/3+1):(ncol(MIND_vectors_split)))]

    if (type=='COMBATLS_covars'){
      MIND_combat_results_1 <-combatls(as.matrix(MIND_vectors_split_1),covariates$Machine_Teslas,covariates[,c('Sex','Case','Age_MRI')],formula)
      cat("\014")
      print('1')
      MIND_combat_results_2 <-combatls(as.matrix(MIND_vectors_split_2),covariates$Machine_Teslas,covariates[,c('Sex','Case','Age_MRI')],formula)
      cat("\014")
      print('2')
      MIND_combat_results_3 <-combatls(as.matrix(MIND_vectors_split_3),covariates$Machine_Teslas,covariates[,c('Sex','Case','Age_MRI')],formula)
      cat("\014")
      print('3')

    }else{
      MIND_combat_results_1 <- combatls(as.matrix(MIND_vectors_split_1),covariates$Machine_Teslas)
      MIND_combat_results_2 <- combatls(as.matrix(MIND_vectors_split_2),covariates$Machine_Teslas)
      MIND_combat_results_3 <- combatls(as.matrix(MIND_vectors_split_3),covariates$Machine_Teslas)

    }

    MIND_vectors_combat <- as.matrix(MIND_vectors)
    MIND_vectors_combat[,setdiff(c(1:ncol(MIND_vectors)),index_zeros)] <- cbind(MIND_combat_results_1$dat.combat,MIND_combat_results_2$dat.combat,MIND_combat_results_3$dat.combat)

  }else{

    if (type=='COMBATLS_covars'){
      MIND_combat_results <-combatls(as.matrix(MIND_vectors[,-index_zeros]),covariates$Machine_Teslas,covariates[,c('Sex','Case','Age_MRI')],formula)

    }else{
      MIND_combat_results <- combatls(as.matrix(MIND_vectors[,-index_zeros]),covariates$Machine_Teslas)

    }

    MIND_vectors_combat <- as.matrix(MIND_vectors)
    MIND_vectors_combat[,setdiff(c(1:ncol(MIND_vectors)),index_zeros)] <- MIND_combat_results$dat.combat
  }

  gc()



  # SAVE RESULTS
  MIND_combat <- as.list(1:nrow(MIND_vectors_combat),)
  i = 1
  for (subject in rownames(MIND_vectors_combat)){

    MIND_combat[[i]] <- matrix(unlist(MIND_vectors_combat[subject,]),nrow=n_regions,ncol=n_regions)
    colnames(MIND_combat[[i]]) <- colnames(MIND)
    rownames(MIND_combat[[i]]) <- colnames(MIND)

    write.csv(MIND_combat[[i]],file=paste0(location,'Code/R/Connectivity/',parcellation,'/',type,'/',subject,'.csv'),row.names = TRUE)

    i = i + 1
  }

  
  # Estimated Total Intracranial Volume (eTIV)
  if (grepl('aparc',parcellation)&!grepl('500_sym',parcellation)){
    if (type=='COMBATLS_covars'){
      
      
      
      etiv_combat_results <-combatls(global[,'EstimatedTotalIntraCranialVol'],covariates$Machine_Teslas,covariates[,'Sex',drop=FALSE],formula_eTIV)
      etivs_combat_results <-combatls(global[,c('BrainSegVol.to.eTIV','MaskVol.to.eTIV','EstimatedTotalIntraCranialVol')],covariates$Machine_Teslas,covariates[,'Sex',drop=FALSE],formula_eTIV,eb=FALSE)
      
      
    }
    
    
    global_combat <- global_combat_results$dat.combat
    
    global_combat[,c('BrainSegVol.to.eTIV','MaskVol.to.eTIV','EstimatedTotalIntraCranialVol')] <- etivs_combat_results$dat.combat
    
    

    # SAVE RESULTS
    if (grepl('longitudinal',parcellation)){
      write.csv(global_combat,file=paste0(location,'datasets/PAFIP/aparc_formatted_1293s_fs_v7.4.1_long_FUP_13-01-2025/COMBATLS',ifelse(all(covariates$Assessment==1),'_cross_sectional',''),'/global_aparc_',type,ifelse(all(covariates$Assessment==1),'_cross_sectional',''),'.csv'),row.names = TRUE)
    }
    
  }
  
}





