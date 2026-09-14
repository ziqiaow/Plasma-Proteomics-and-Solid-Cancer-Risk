library(readstata13)
library(data.table)

#replace cancer name with each target site

in.dir = "//win.ad.jhu.edu/cloud/welchcenter$/ARIC/users/zwang389/cancer/liver/"
visit_data_dir = "//win.ad.jhu.edu/cloud/aric-cancer$/STATA NP/"

#liver cancer data
f.aric_derive2_liver <- paste(visit_data_dir,"liver_cancer_2015_np.dta",sep="") 
d.aric_derive2_liver <- read.dta13(f.aric_derive2_liver) 
dim(d.aric_derive2_liver) 
#[1] 15653   17

#overview of the data
table(d.aric_derive2_liver$prvcancr)
#   0    1 
#14743   910 
table(d.aric_derive2_liver$incidence) #total incident cases
#   0    1 
#14709    34 
table(d.aric_derive2_liver$prvcancr,d.aric_derive2_liver$incidence, useNA = "ifany")
#        0    1 <NA>
# 0 14709    34     0
#1     0     0   910


length(which(is.na(d.aric_derive2_liver$incidence_t1))) #prevalent cases at baseline before V1
#[1] 910
length(which(is.na(d.aric_derive2_liver$incidence_t2) & (d.aric_derive2_liver$prvcancr == 1 | d.aric_derive2_liver$incidence==1))) #prevalent cases or missing data at V2
#[1] 915

length(which(is.na(d.aric_derive2_liver$incidence_t2) & (d.aric_derive2_liver$incidence==1))) #incident at V1
#[1] 5

length(which(!is.na(d.aric_derive2_liver$incidence_t2) & (d.aric_derive2_liver$incidence==1))) #incident at or after v2
#[1] 29

length(which(d.aric_derive2_liver$incidence_t2 < 5 & d.aric_derive2_liver$incidence == 1)) #incident cases from V2 within 5 yrs
#[1]1
length(which(d.aric_derive2_liver$incidence_t2 > 5 & d.aric_derive2_liver$incidence == 1)) #incident cases from V2 after 5 yrs
#[1] 28

d.aric_derive2_liver$incidence_after_v2 = 0
d.aric_derive2_liver$incidence_after_v2[which(!is.na(d.aric_derive2_liver$incidence_t2) & (d.aric_derive2_liver$incidence==1))] = 1
d.aric_derive2_liver$incidence_after_v2[which(is.na(d.aric_derive2_liver$incidence_t2) & (d.aric_derive2_liver$prvcancr == 1 | d.aric_derive2_liver$incidence==1))] = NA #prevalent cases or missing data at V2
table(d.aric_derive2_liver$incidence_after_v2, useNA = "ifany")
#   0    1 <NA> 
#14526   190   937 
save(d.aric_derive2_liver,file=paste0(in.dir,"liver_incident_v2.RData"))

#load clinical data
load("//win.ad.jhu.edu/cloud/welchcenter$/ARIC/users/zwang389/cancer/colon/clinical_keep_v2.RData")
colnames(d.aric_derive2_liver)[1]="ID"
liver_final = merge(d.aric_derive2_liver,clinical_keep,by="ID")
save(liver_final,file=paste0(in.dir,"liver_incident_v2.RData"))

