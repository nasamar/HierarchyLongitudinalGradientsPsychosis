################################################################################
# Script to perform multiple regression (lm) and linear mixed modeling (lmer) 
# on MIND networks
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

######################  Local structural Longitudinal Analisis ###############################################################


# 1. Initial preparation
####################################################################################

# 1.1. Clean Memory and load libraries
rm(list=ls())

library(lmerTest)

location <- 'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/R/'
source(paste0(location,'regional_brainmap_representation_borders.R'))
source(paste0(location,'regional_brainmap_representation_subcortical.R'))



combat <- 'COMBATLS_covars'
# combat <- 'COMBATLS'
# combat <- 'COMBAT'
# combat <- ''

parcellation <- 'aparc_500_sym'
parcellation <- 'subcortical'


y_var <- 'degrees'
y_var <- 'gradients'


gradient <- '_G1'
# gradient <- ''
y_var_G <- paste0(y_var,gradient)

normalization <- ''
normalization <- '/FEP+CN'
# normalization <- '/FEP_CN'
# normalization <- '/CN'

raw_or_residuals <- 'raw'
# raw_or_residuals <- 'residuals'



sessions <- c(1,2,3,5,10,15,20)
# sessions <- c(1,2,3,5,10,15)


GeneralizedLinearModelGauss <- function(D,type){
  # Zones in Desikan-Killiany atlas
  if (type=='mixed'){
    names_Objetives <- colnames(D)[c((which(colnames(D)=='mean_euler')+1):ncol(D))]
    
  }else if (grepl('global',type)){
    names_Objetives <- paste0('global_',y_var)
  }else if (type=='no mixed regional'){
    names_Objetives <- colnames(D)[c((which(colnames(D)=='mean_euler')+1):(ncol(D)-1))] 
  }

  if (!grepl('no mixed',type)){  
    # Names of variables used in interpolation. Used as names of columns in interpolationResults table
   v<-c('dx1','Age_inclusion','Sex1','eTIV','mean_euler','CPZ_equivalent','Treatment_Time' ,'dx1:Treatment_Time' ,'CPZ_equivalent:Treatment_Time') 
    
    if (!all(D[,'protocolChange15']==0)){
      v<-c(v,'protocol_Change15')
    }
    
  }else{
    v<-c('dx1','Age_inclusion','Sex1','eTIV','mean_euler')
    
  }

  # construction of Table for saving beta and p values
  res <- data.frame(matrix(0,ncol=length(names_Objetives),nrow=length(v)))
  names(res)<-names_Objetives
  row.names(res)<-v
  
  resPValues <- res
  
  
  n<-length(v)+1
  
  if (!grepl('global',type)){
    for(i in 1:length(names_Objetives)){
      
      y<-names_Objetives[[i]]
      D[y]<-as.vector(scale(D[y]))
      
      if (type=='mixed'){
        formulaText<-paste0(y,'~1+Case+Age_inclusion+Sex+eTIV+mean_euler+CPZ_equivalent+Treatment_Time+Case:Treatment_Time+CPZ_equivalent:Treatment_Time',ifelse(!all(D[,'protocolChange15']==0),'+protocolChange15',''),'+(1|Subject)',sep='')

        # print(formulaText)
        
        Model <- lmer(formulaText,data=D)

        
      }else if (type=='no mixed regional'){
        formulaText<-paste0(y,'~1+Case+Age_inclusion+Sex+eTIV+mean_euler',sep='')
        
        Model <- lm(formulaText,data=D)
          
      }
      # print(summary(Model))
      res[,i]<-(as.numeric((coef(summary(Model)))[2:n,ncol(coef(summary(Model)))-1]))
      resPValues[,i]<-(as.numeric((coef(summary(Model)))[2:n,ncol(coef(summary(Model)))]))
      
    }
    

  }else{
    D[paste0('global_',y_var)]<-scale(D[paste0('global_',y_var)])
    
    if (type=='mixed global'){
      formulaText<-paste0('global_',y_var,'~1+Case+Age_inclusion+Sex+eTIV+mean_euler+CPZ_equivalent+Treatment_Time+Case:Treatment_Time+CPZ_equivalent:Treatment_Time',ifelse(!all(D[,'protocolChange15']==0),'+protocolChange15',''),'+(1|Subject)',sep='')
      
      print(formulaText)
      Model <- lmer(formulaText,data=D)
      
      
    }else if (type=='no mixed global'){
      formulaText<-paste0('global_',y_var,'~1+Case+Age_inclusion+Sex+eTIV+mean_euler',sep='')
      
      Model <- lm(formulaText,data=D)
        
    }
    
    print(summary(Model))
    
    res[,1]<-(as.numeric((coef(summary(Model)))[2:n,ncol(coef(summary(Model)))-1]))
    resPValues[,1]<-(as.numeric((coef(summary(Model)))[2:n,ncol(coef(summary(Model)))]))
    
  }
  

  # transpose results tables for correction of p-values
  res<-t(res)
  res_complete <- res
  resPValues<-t(resPValues)
  
  # FRD correction of P values and selection of significant betas
  for(z in 1:ncol(res)){
    resPValues[,z]<- p.adjust(resPValues[,z], method="fdr")
  }
  
  for(i in 1:nrow(res)){
    for(j in 1:ncol(res)){
      if(resPValues[i,j]>0.05){
        res[i,j]<-NaN  # we set to NULL the beta parameters of the variables with p values < 0.05
      }
    }
    
  }
  
  return(list(res = res, res_complete = res_complete))
}


# 1.3. DATA LOADING
Data_Base_covariates <-'C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/Covariates_complete.csv'
covariates <-read.csv(Data_Base_covariates,row.names = 1)

Euler <- read.csv('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/euler_longitudinal.csv',row.names = 1)
rownames(Euler) <- gsub('ses-','',gsub('sub-','',rownames(Euler)))

eTIV <- read.csv(paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/datasets/PAFIP/aparc_formatted_1293s_fs_v7.4.1_long_FUP_13-01-2025/COMBATLS/global_aparc_',combat,'.csv'),row.names = 1)
eTIV <- eTIV[,'EstimatedTotalIntraCranialVol',drop=FALSE]


# Degrees
Data_Base_FEP <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/3. MIND_long/Data/degree/',parcellation,'/',combat,'/degree_68_FEP',ifelse(raw_or_residuals=='residuals','_residuals',''),'.csv')
Data_Base_CN <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/3. MIND_long/Data/degree/',parcellation,'/',combat,'/degree_68_CN',ifelse(raw_or_residuals=='residuals','_residuals',''),'.csv')
  

# Gradients
if (y_var == 'gradients'){
  Data_Base_FEP <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients_normalized/',parcellation,'/',combat,normalization,'/',raw_or_residuals,'/',y_var_G,'_FEP.csv')
  Data_Base_CN <-paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients_normalized/',parcellation,'/',combat,normalization,'/',raw_or_residuals,'/',y_var_G,'_CN.csv')
}

y_var_FEP <-read.csv(Data_Base_FEP, row.names = 1)
y_var_CN <-read.csv(Data_Base_CN, row.names = 1)
morphometry_var <- rbind(y_var_CN,y_var_FEP)

covariates <- covariates[rownames(morphometry_var),]
Euler <- Euler[rownames(covariates),]
eTIV <- eTIV[rownames(covariates),]

covariates$Treatment_Time <- covariates$Age_MRI-covariates$Age_inclusion
covariates$eTIV <- eTIV


morphometry_var_longitudinal <- data.frame(morphometry_var[rownames(covariates),],row.names = rownames(covariates))

covariates_longitudinal <- covariates[,c(1:9,which(colnames(covariates)%in%c('eTIV','CPZ_equivalent','Treatment_Time','protocolChange15')))]


Data_Base_Global <- merge(merge(covariates_longitudinal,Euler,by=0),morphometry_var_longitudinal,by.x = 1,by.y = 0)

Data_Base_Global <- Data_Base_Global[Data_Base_Global$Assessment %in% sessions,]
rownames(Data_Base_Global) <- Data_Base_Global[,1]



# LONGITUDINAL
{
# for participants in base line we set medication to zero
Data_Base_Global$CPZ_equivalent<-ifelse(Data_Base_Global$Assessment==1,0,Data_Base_Global$CPZ_equivalent)
# Creation of a binary variable dx for diagnosis (dx1 psychotic, 0 control)
# To Controls participants we assigned medication  Zero
Data_Base_Global$CPZ_equivalent<-ifelse(Data_Base_Global$Case==0,0.0,Data_Base_Global$CPZ_equivalent)
# fixed dx variable as an R factor
Data_Base_Global$Case<-as.factor(Data_Base_Global$Case)
# Normalize continuous relevant variables
Data_Base_Global$Age_inclusion<-as.vector(scale(Data_Base_Global$Age_inclusion))
Data_Base_Global$Treatment_Time<-as.vector(scale(Data_Base_Global$Treatment_Time))
Data_Base_Global$CPZ_equivalent<-as.vector(scale(Data_Base_Global$CPZ_equivalent))
Data_Base_Global$eTIV<-as.vector(scale(Data_Base_Global$eTIV))
Data_Base_Global$mean_euler <- scale(Data_Base_Global$mean_euler)
# fixed Sex variable as an R factor
Data_Base_Global$Sex<-as.factor(Data_Base_Global$Sex)


Data_Base_Global_long <- Data_Base_Global[rownames(covariates_longitudinal),]

InterpolationResults <- GeneralizedLinearModelGauss(Data_Base_Global,'mixed')
print(InterpolationResults)
if (y_var=='degrees'){
  write.csv(InterpolationResults,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/degrees/',parcellation,'/InterpolationResults_',y_var,'.csv'),row.names = TRUE)
}else{
  write.csv(InterpolationResults,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients',ifelse(grep('CN',normalization),'_normalized',''),'/',parcellation,'/',combat,normalization,'/',raw_or_residuals,'/InterpolationResults_',y_var,ifelse(y_var%in%c('gradients','range','std'),gradient,''),'.csv'),row.names = TRUE)
  
}

Data_Base_Global[,paste0('global_',y_var)] = apply(as.matrix(Data_Base_Global[,c((which(colnames(Data_Base_Global)=='mean_euler')+1):ncol(Data_Base_Global))]), 1, mean, na.rm = TRUE)
InterpolationResults_global <- GeneralizedLinearModelGauss(Data_Base_Global,'mixed global')
print(InterpolationResults_global)

}

# BASELINE
Data_Base_Global_baseline = Data_Base_Global[Data_Base_Global$Assessment==1,]
Data_Base_Global_baseline[,paste0('global_',y_var)] = apply(as.matrix(Data_Base_Global_baseline[,c((which(colnames(Data_Base_Global)=='mean_euler')+1):ncol(Data_Base_Global_baseline))]), 1, mean, na.rm = TRUE)
InterpolationResults_baseline_regional <- GeneralizedLinearModelGauss(Data_Base_Global_baseline,'no mixed regional')
if (y_var=='degrees'){
  write.csv(InterpolationResults_baseline_regional,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/degrees/',parcellation,'/InterpolationResults_',y_var,'_baseline.csv'),row.names = TRUE)
}else{
  write.csv(InterpolationResults_baseline_regional,file=paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Code/MATLAB/Connectivity/Longitudinal/gradients',ifelse(grep('CN',normalization),'_normalized',''),'/',parcellation,'/',combat,normalization,'/',raw_or_residuals,'/InterpolationResults_',y_var,ifelse(y_var%in%c('gradients','range','std'),gradient,''),'_baseline.csv'),row.names = TRUE)

}

InterpolationResults_baseline_global <- GeneralizedLinearModelGauss(Data_Base_Global_baseline,'no mixed global')
print(InterpolationResults_baseline_global)


print('End')


##############################################################################################################
# 3.1 . Representation of data
if (parcellation == 'subcortical'){
  variables<- c('dx1','Age_inclusion','Sex1','eTIV')
  for (i in variables){
    data<-data.frame(values=InterpolationResults_baseline_regional$res_complete[,i],sig=!is.na(InterpolationResults_baseline_regional$res[,i]))
    map<-eval(parse(text = paste0('regional_brainmap_representation_',parcellation)))
    
    print(map(data,paste0('Baseline ', y_var,'-> Fixed effect: ',i),max(InterpolationResults_baseline_regional$res_complete,na.rm = TRUE),min(InterpolationResults_baseline_regional$res_complete,na.rm = TRUE),0,'horizontal'))
    ggsave(paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Molecular/plots/Connectivity_longitudinal/Baseline_', y_var,'_Fixed_effect_',gsub(':','_',i),'.png'), width = 8, height = 5, dpi = 300)
  
  }
  
  variables<- c('dx1','Age_inclusion','Sex1','eTIV','mean_euler','CPZ_equivalent','Treatment_Time','dx1:Treatment_Time','CPZ_equivalent:Treatment_Time','protocolChange15' )

  for (i in variables){
  
    data<-data.frame(values=InterpolationResults$res_complete[,i],sig=!is.na(InterpolationResults$res[,i]))

    map<-eval(parse(text = paste0('regional_brainmap_representation_',parcellation)))
      
    print(map(data,paste0('Local ', y_var,'-> Fixed effect: ',i),max(InterpolationResults$res_complete,na.rm = TRUE),min(InterpolationResults$res_complete,na.rm = TRUE),0,'horizontal'))
    ggsave(paste0('C:/Users/usuario/OneDrive - UNIVERSIDAD DE SEVILLA/Natalia/Molecular/plots/Connectivity_longitudinal/Local_', y_var,'_Fixed_effect_',gsub(':','_',i),'.png'), width = 8, height = 5, dpi = 300)
    
  }

}
