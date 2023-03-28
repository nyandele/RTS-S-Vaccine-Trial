
######################################### RTS,S Dataset Analysis ##############################################
##########################################################################################################

########### Demograpy Dataset
rtss_dem = read.table("demography.csv", header = TRUE, sep = ',')
str(rtss_dem)
head(rtss_dem)
names(rtss_dem)

#omitting unneeded columns (and rows) in Demography dataset
rtss_dem = rtss_dem[, -c(1,3,5:10,12:15,17:22,25:30)]
rtss_dem = rtss_dem[, -c(12,13,16:18,20, 22,24,25,26,28,29,32:46)]
rtss_dem = rtss_dem[, -c(24,25,26,30,31,32,36)]

#Rename some variables and omit unneeded columns
colnames(rtss_dem)
names(rtss_dem)[6] = "AGE_M"
names(rtss_dem)[13] = "AGE_W"
rtss_dem = rtss_dem[, -c(2,12)]
names(rtss_dem)[14] = "AGE3_M"
names(rtss_dem)[17] = "AGE4_W"
rtss_dem = rtss_dem[, -c(13,16)]
rtss_dem = rtss_dem[, -c(4)]
rtss_dem = rtss_dem[, -c(1,9:14,16:18,23)]

#Rearrange columns in suitable order
colnames(rtss_dem)
rtss_dem = rtss_dem[,c(1,12,13,2,3,8,4,5,6,7,9,10,11)]

##Correcting Date formulas
str(rtss_dem)
rtss_dem$STRTITT = strptime(rtss_dem$STRTITT, format="%d-%b-%y")
rtss_dem$STRTITT = as.Date(rtss_dem$STRTITT, format = "%Y-%m-%d")
rtss_dem$ENDITT_SE = strptime(rtss_dem$ENDITT_SE, format="%d-%b-%y")
rtss_dem$ENDITT_SE = as.Date(rtss_dem$ENDITT_SE, format = "%Y-%m-%d")
str(rtss_dem)

#Rearrange columns and delete unneeded ones
names(rtss_dem)
rtss_dem = rtss_dem[,-c(2,13)]
rtss_dem = rtss_dem[,c(1,3,4,6,7,8,2,9,5,10,11)]

########### Hematology Dataset
library(readxl)
rtss_hem = read_excel("Mal05.Heamatology_ummi.xlsx")
str(rtss_hem)
#omitting unneeded columns in Hematology dataset
names(rtss_hem)
rtss_hem = rtss_hem[, -c(1,4,5:14,16:26,28:34)]
names(rtss_hem)[5] = "dat_run"
rtss_hem = rtss_hem[,c(1,2,5,3,4)]

#Sorting rows by PID and dat_run
str(rtss_hem)
rtss_hem$dat_run = as.Date(rtss_hem$dat_run)
str(rtss_hem)
rtss_hem = with(rtss_hem, rtss_hem[order(PID, dat_run),])
rtss_hem = rtss_hem[!duplicated(rtss_hem$PID),]
table(duplicated(rtss_hem))

############ Merge Demography to Hematology
rtss_dem_hem = merge(rtss_dem, rtss_hem, by = c("PID"))

#if you want to export rtss_dem_hem into your directory, use this code
write.csv(rtss_dem_hem, file = "rtss_dem_hem3.csv")

#Check for 0's and NAs wrt monocytes and lymphocytes
summary(rtss_dem_hem$mono_n==0.00)
summary(rtss_dem_hem$lym_n==0.00)
summary(is.na(rtss_dem_hem$mono_n))
summary(is.na(rtss_dem_hem$lym_n))
#Omit rows with 0's and NA's wrt mono and lym
rtss_dem_hem = rtss_dem_hem[!(rtss_dem_hem$mono_n==0.00),]
rtss_dem_hem = rtss_dem_hem[!(rtss_dem_hem$lym_n==0.00),]

########### Treatment allocation dataset

### Import dataset with treatment allocation
treat_allocation = read.table("Treatment allocation.csv", header = TRUE, sep = ',')
names(treat_allocation)
treat_allocation = treat_allocation[, -c(2:70,72:98)]
str(treat_allocation)
names(treat_allocation)
names(treat_allocation)[2] = "Treat_allocation"

### Merge demography, hematology, and treatment allocation information
rtss_dem_hem_alloc = merge(rtss_dem_hem, treat_allocation, by = c("PID"))
names(rtss_dem_hem_alloc)
str(rtss_dem_hem_alloc)
names(rtss_dem_hem_alloc)[2] = "Vaccination_date"
names(rtss_dem_hem_alloc)[3] = "End_fu_date"
names(rtss_dem_hem_alloc)
rtss_dem_hem_alloc = rtss_dem_hem_alloc[,c(1,2,3,16,4,5,6,7,8,9,10,11,12,13,14,15)]

########### Parasitology dataset
library(readxl)
rtss_par = read_excel("pr-mal055-try.xlsx")
str(rtss_par)
names(rtss_par)

#omitting unneeded columns in Parasitology dataset
rtss_par = rtss_par[, -c(1,5:12,16:28,30,31)]
names(rtss_par)

#Making a subset of entries from the first reading only
rtss_par = subset(rtss_par, rtss_par$reading==1)

#omitting more unneeded columns
names(rtss_par)
rtss_par = rtss_par[, -c(2,4,5,7)]
names(rtss_par)[2] = "diagnosis_date"

#Changing to date format
str(rtss_par)
rtss_par$diagnosis_date = as.Date(rtss_par$diagnosis_date)

#Sorting rows by studyno and date_taken
rtss_par = with(rtss_par, rtss_par[order(studyno, diagnosis_date),])
rtss_par = rtss_par[!duplicated(rtss_par$studyno),]

#### Merge rtss_par with rtss_dem_hem_alloc
rtss_all = merge(rtss_dem_hem_alloc, rtss_par, by = c("studyno"))



##################################### XXXXXXXXXXXXXXX Cox Regression  XXXXXXXXXXXXXXXXXXX #######################

#Convert age category (agec) from character to factor
class(rtss_initial_01$agec)
rtss_initial_01$agec = factor(rtss_initial_01$agec)

######Testing model assumptions

### Proportional Hazards Assumption  --> test whether the relative difference in hazards between the groups stays constant over time
#   --> this is tested in all covariates of the model

# Testing for proportional hazards assumption
library(survival)
load(survminer)
library(survminer)
residual_cox = coxph(Surv(time, malaria) ~ treat_allocation_cat + agec + sex + outp_distance + bednt_ + ML, data = rtss_initial_01)
residual_cox
test_ph = cox.zph(residual_cox)
test_ph

final_model = coxph(formula = Surv(time, malaria) ~ outp_distance + treat_allocation_cat + 
                      agec + sex + bednt_ + I(((ML + 0.1)/0.1)^1)*treat_allocation_cat, data = rtss_initial_01)
final_model

### Strata-specific adjusted hazard rates aHR

fit <- survreg(Surv(time, malaria) ~ outp_distance + treat_allocation_cat + 
                 agec + sex + bednt_ + I(((ML + 0.1)/0.1)^1)*treat_allocation_cat, data = rtss_initial_01,
               dist = "weibull")
summary(fit)

fit2 <- survreg(Surv(time, malaria) ~ outp_distance + treat_allocation_cat + ML +
                  agec + sex + bednt_ + ML*treat_allocation_cat, data = rtss_initial_01,
                dist = "weibull")

# Calculate the strata-specific adjusted hazard rates
strat <- strata(fit2, transform = "weibull")


## another try
datadist(rtss_initial_01)
label(outp_distance) <- "outp_distance"
label(treat_allocation_cat) <- "treat_allocation_cat"
label(ML) <- "ML"
label(agec) <- "agec"
label(sex) <- "sex"
label(bednt_) <- "bednt_"

fit <- cph(Surv(time, malaria) ~ outp_distance + treat_allocation_cat + ML +
             agec + sex + bednt_ + ML*treat_allocation_cat, x=TRUE, y=TRUE,
           surv = TRUE, strata =treat_allocation_cat, data = rtss_initial_01)
summary(fit)
str(cph)
names(rtss_initial_01)
