
#### Puja 2024 try with cretaceous data ######


# Calling packges

install.packages("readxl")  # Run this ONCE
library(readxl)             # Load the package EVERY time

library(divDyn)
library(icosa)
library(betapart)
library(ggplot2)
library(vegan) ### needed for jacBeta function

# Insert original data from the PBDB 
# Delete the lines above the data in the table before insert

setwd("D:\\Thesis")
bivdata <- read_excel("D:/Thesis/AmalaData.xlsx")



# Extract the data needed for analysis in this paper from the original data

needcol = c("order",
            "family",
            "genus...25",
            "paleolng",
            "paleolat",
            "time_bins",
            "early_interval",
            "late_interval",
            "max_ma",
            "min_ma")
bivdata = bivdata[,needcol]

need_stages = c("Danian",
                "Selandian",
                "Thanetian",
                "Ypresian",
                "Lutetian",
                "Bartonian",
                "Priabonian",
                "Rupelian",
                "Chattian",
                "Aquitanian",
                "Burdigalian",
                "Langhian",
                "Serravallian",
                "Tortonian",
                "Messinian",
                "Zanclean",
                "Piacenzian",
                "Gelasian",
                "Clabrian",
                "Chibanian",
                "Late Pleistocene") 

# Two more stages at the begining and ending of studied interval were necessary for GF Ext/Org rates

bivdata = bivdata[which(bivdata$time_bins %in% need_stages),]


# Assign a stage number for each occurrence

stgno = rep(NA,nrow(bivdata))
for(i in 1:21){
  stgno[which(bivdata$time_bins == need_stages[i])] = i
}
bivdata = cbind.data.frame(bivdata,stgno)


# Removing occurrence without generic information
#bivdata5 = bivdata4[-which(bivdata4$genus == ""),]  

# Removing occurrence without paleocoordinate information

bivdata$paleolat = as.numeric(bivdata$paleolat)
bivdata$paleolng = as.numeric(bivdata$paleolng)

needcoor = intersect(which(-90 < bivdata$paleolat & bivdata$paleolat < 90),
                     which(-180 < bivdata$paleolng & bivdata$paleolng < 180))

bivdata = bivdata[needcoor,]  

nrow(bivdata) #finally 36058 occurrences were used in subsequent analysis.

# Export the data if needed
# write.csv (bivdata, file ="C:\\Users\\YJ\\Documents\\Bivdata-used.csv")


# Input the verified data here

#setwd("C:\\Users\\YJ\\Documents")
#bivdata = read.csv("Supplementary Tab. S1.csv",header = TRUE)


# Assign geogrids for each occurrence

lng = bivdata$paleolng
lat = bivdata$paleolat

hLow = hexagrid(c(5,5))
bivdata$geocell = locate(hLow, cbind(lng,lat), randomborder = F)


# Calculate the generic diversity& SQS diversity& origination rate & extinction rate

raw_result = divDyn(bivdata,tax = "genus...25", bin = "stgno")

plot(raw_result$divSIB,ylab = " raw richness",xlab="stgno",type="b",xaxt="n")

sqs0.8 = subsample(bivdata, iter=50, q=0.8, tax="genus...25", bin="stgno", type="sqs")

plot(sqs0.8$divSIB,ylab = "SQS richness",xlab="binno",type="b",xaxt="n")


# Export diversification results

write.csv(raw_result, file = "raw_result_Creta.csv")
write.csv(sqs0.8, file = "SQS0.8_Creta.csv")



# Calculate & Resample BC & beta-diversity for each stage

bivstagebcresult = matrix(data = NA, 21, 11)  
colnames(bivstagebcresult) = c("meanBC","SDBC","lowCIBC","upCIBC","genus","geocell","occurrences",
                               "jacbeta$beta.JAC","meanjbeta","lowjbeta","upjbeta")         
for (t in 1:21){
  bivstagebcresult[t,] = BC(t)
  print(t)
}   

# Export BC & beta-diversity results

write.csv(bivstagebcresult,file = "bivstagebcresul_Creta.csv")



##################### Don't need the rest ###########################################

########################################################################################


# Calculate BC of survivors and newcomers for each stage

#Sur_New_result = matrix(data = NA, 12, 12)
#colnames(Sur_New_result) = c("surmeanBC","SDBC","lowCIBC","upCIBC","surgenus","surgeocell",
                             "newmeanBC","SDBC","lowCIBC","upCIBC","newgenus","newgeocell")
#for (t in 3:14){
  #Sur_New_result[t-2,] = SNBC(t)
  #print(t)
#}


# Export BC of survivors and newcomers
#write.csv(Sur_New_result,file = "Sur_New_result_Creta.csv")                                



# Calculate BC of survivors in previous stage 


#prestagesurbivbcresult = matrix(data = NA, 12, 6)
#colnames(prestagesurbivbcresult) = c("K1meanBC","K1SDBC","K1lowCIBC","K1upCIBC","K1genus","K1geocell")

#for (t in 3:14){
  #prestagesurbivbcresult[t-2,] = SPBC(t)
  #print(t)
#}
#prestagesurbivbcresult


# Export BC of survivors in previous stage 

#write.csv(prestagesurbivbcresult,file = "prestagesurbivbcresult_Creta.csv")


############################################################################################

# Subsample BC & beta-diversity  of Maastrichtian, Danian

occdata2 = bivdata

subbcresult = matrix(data = NA, 20000, 2)
colnames(subbcresult) = c("BC","stage")
subbcresult = as.data.frame(subbcresult)

sub_beta_result = matrix(data = NA, 20000, 2)
colnames(sub_beta_result) = c("beta","stage")
sub_beta_result = as.data.frame(sub_beta_result)


# Maastrichtian
subbcresult$BC[1:10000] = subBC(13,40) ### 13 is the Maastrichtian stage no., 40 is the resampled geocell no. 
subbcresult$stage[1:10000] = "Maastrichtian"
sub_beta_result$beta[1:10000] = MSBD(13,40)
sub_beta_result$stage[1:10000] = "Maastrichtian"

# Danian
subbcresult$BC[10001:20000] = subBC(14,40)
subbcresult$stage[10001:20000] = "Danian"
sub_beta_result$beta[10001:20000] = MSBD(14,40)
sub_beta_result$stage[10001:20000] = "Danian"


##########################################################################################

# Plot the resample results

ggplot(subbcresult, aes(BC)) + geom_histogram()

ggplot(subbcresult, aes(BC, stage,fill = stage,colour = stage)) +
  geom_violin(alpha = 0.5)

ggplot(sub_beta_result, aes(beta)) + geom_histogram()

ggplot(sub_beta_result, aes(beta, stage,fill = stage,colour = stage)) +
  geom_violin(alpha = 0.5)



# Calculate the 95% confidence interval of the subsample results

BC_dat = subbcresult$BC[which(subbcresult$stage == "Maastrichtian")]
c(mean(BC_dat),sd(BC_dat),
  quantile(BC_dat,0.025,na.rm = TRUE),quantile(BC_dat,0.975,na.rm = TRUE))

BC_dat = subbcresult$BC[which(subbcresult$stage == "Danian")]
c(mean(BC_dat),sd(BC_dat),
  quantile(BC_dat,0.025,na.rm = TRUE),quantile(BC_dat,0.975,na.rm = TRUE))


beta_dat = sub_beta_result$beta[which(sub_beta_result$stage == "Maastrichtian")]
c(mean(beta_dat),sd(beta_dat),
  quantile(beta_dat,0.025,na.rm = TRUE),quantile(beta_dat,0.975,na.rm = TRUE))

beta_dat = sub_beta_result$beta[which(sub_beta_result$stage == "Danian")]
c(mean(beta_dat),sd(beta_dat),
  quantile(beta_dat,0.025,na.rm = TRUE),quantile(beta_dat,0.975,na.rm = TRUE))




# Export BC & beta-diversity subsample results by equal cells

write.csv(subbcresult,file = "subbcresult_Creta.csv")

write.csv(sub_beta_result,file = "sub_beta_result_Creta.csv")



#########################################################################################




#Functions



BC = function(stageno){
  stgocc = bivdata[which(bivdata$stgno == stageno),]
  n.row = length(levels(factor(stgocc$genus))) #number of genus in a stage
  n.col = length(levels(factor(stgocc$geocell))) #number of locations in a stage
  if (n.row < 10 | n.col<5){
    return(rr <<- 0)}
  else{rr <<- 1}  
  
  # Create geocell-taxa metrix for a stage
  
  a = nrow(stgocc) 
  levels(factor(stgocc$genus))  
  n.row = length(levels(factor(stgocc$genus))) #number of genus in a stage 
  n.col = length(levels(factor(stgocc$geocell))) #number of locations in a stage 
  initial.matrix = matrix (data = 0, n.row,n.col) 
  rownames(initial.matrix) = levels(factor(stgocc$genus)) 
  colnames(initial.matrix) = levels(factor(stgocc$geocell)) 
  
  
  for (i in 1:a){
    gname = levels(factor(stgocc$genus[i]))    
    cname = levels(factor(stgocc$geocell[i])) 
    initial.matrix[which(rownames(initial.matrix) == gname),which(colnames(initial.matrix) == cname)] = 1  
  }
  
  bcdata = initial.matrix
  mbetamartix = t(bcdata)   
  
  jacbeta = beta.multi(mbetamartix,index.family = "jaccard")        # calculate jacbeta directly when resampling is not required 
  
  
  #bootstrap??sample size equals to generic richness??replace = TRUE
  
  #subsample??sample size = 90%??80%??50% generic richness??replace = FALSE 
  
  
  BC.rs = rep(NA,10000)
  Rejbeta = rep(NA,10000)
  b = rep(NA,n.col)
  for(i in 1:10000){
    ssample<-sample(n.row,size=round(0.9*nrow(bcdata)),replace=FALSE)  # subsample by 90% genera
    bootsp<-bcdata[ssample,]
    subsamcolsum = as.matrix(colSums(bootsp))
    if(sum(which(subsamcolsum == 0))>0)      # remove taxa without occurrences in subsamples
    {bootsp = bootsp[,-which(subsamcolsum == 0)]}
    o<-sum(bootsp)
    BC.rs[i] = (o-nrow(bootsp))/(ncol(bootsp)*nrow(bootsp)-nrow(bootsp))
    tbootsp=t(bootsp)
    jbeta= beta.multi(tbootsp,index.family = "jaccard") 
    Rejbeta[i] = jbeta$beta.JAC
  }
  hist(BC.rs)
  hist(Rejbeta)
  
  #write.csv(BC.rs,file = paste(as.character(t),".csv"))  #Export subsample results of each stage if necessary
  #write.csv(Rejbeta,file = paste(as.character(t),".csv"))
  
  meanBC<-mean(BC.rs,na.rm=TRUE)
  SDBC<-sd(BC.rs,na.rm=TRUE)
  lowCIBC<-quantile(BC.rs,0.025,na.rm=TRUE)
  upCIBC<-quantile(BC.rs,0.975,na.rm=TRUE)
  meanjbeta<-mean(Rejbeta,na.rm=TRUE)
  lowjbeta<-quantile(Rejbeta,0.025,na.rm=TRUE)
  upjbeta<-quantile(Rejbeta,0.975,na.rm=TRUE)
  
  sur.resl<-c(meanBC,SDBC,lowCIBC,upCIBC,nrow(bcdata),ncol(bcdata),a,
              jacbeta$beta.JAC,meanjbeta,lowjbeta,upjbeta)
  sur.resl
  
} # Whole stage BC&beta function

SNBC = function(k){
  stagedata = bivdata[which(bivdata$stgno == k),]     #Stagedata of the k stage
  stagegenus = levels(factor(stagedata$genus))       
  genusno = length(stagegenus)              
  genustype = rep(NA,genusno)                 
  pre_one_stagedata = bivdata[which(bivdata$stgno == k-1),] 
  pre_stagedata = bivdata[which(bivdata$stgno < k),]  #pre-stagedata= data before k stage
  
  genustype[which(stagegenus %in% pre_one_stagedata$genus)] = "survivor"     #define survivor
  genustype[-which(stagegenus %in% pre_stagedata$genus)] = "newcomer"        #define newcomer
  
  genustypedata = cbind.data.frame(stagegenus,genustype)                      
  survivorgenus = stagegenus[which(genustypedata$genustype == "survivor")]    
  newcomergenus = stagegenus[which(genustypedata$genustype == "newcomer")]    
  
  # Survivor matrix
  sn.row = length(survivorgenus)                             
  sn.col = length(levels(factor(stagedata$geocell)))          
  sinitial.matrix = matrix (data = 0, sn.row,sn.col)          
  rownames(sinitial.matrix) = survivorgenus 
  colnames(sinitial.matrix) = levels(factor(stagedata$geocell)) 
  
  # Newcomer matrix
  nn.row = length(newcomergenus)
  nn.col = length(levels(factor(stagedata$geocell)))
  ninitial.matrix = matrix (data = 0, nn.row,nn.col) 
  rownames(ninitial.matrix) = newcomergenus 
  colnames(ninitial.matrix) = levels(factor(stagedata$geocell))       
  
  # calculating survivor matrix
  for (i in 1:nrow(stagedata)){
    gname = levels(factor(stagedata$genus[i]))
    cname = levels(factor(stagedata$geocell[i]))
    sinitial.matrix[which(rownames(sinitial.matrix) == gname),which(colnames(sinitial.matrix) == cname)] = 1
  }
  survivormatrix = sinitial.matrix
  scolsum = as.matrix(colSums(survivormatrix))
  if(sum(which(scolsum == 0))>0)
  {survivormatrix = survivormatrix[,-which(scolsum == 0)]}             
  
  # calculating newcomer matrix
  for (i in 1:nrow(stagedata)){
    gname = levels(factor(stagedata$genus[i]))
    cname = levels(factor(stagedata$geocell[i]))
    ninitial.matrix[which(rownames(ninitial.matrix) == gname),which(colnames(ninitial.matrix) == cname)] = 1
  }
  newcomermatrix = ninitial.matrix
  ncolsum = as.matrix(colSums(newcomermatrix))
  if(sum(which(ncolsum == 0))>0)
  {newcomermatrix = newcomermatrix[,-which(ncolsum == 0)]}          
  
  # Calculating BC
  bootbc = function(n){
    bcdata = n
    BC.rs = rep(NA,10000)
    for(i in 1:10000){                                                   
      ssample<-sample(nrow(bcdata),size=nrow(bcdata),replace=TRUE)
      bootsp<-bcdata[ssample,]
      subsamcolsum = as.matrix(colSums(bootsp))
      if(sum(which(subsamcolsum == 0))>0)
      {bootsp = bootsp[,-which(subsamcolsum == 0)]}                    
      o<-sum(bootsp)
      BC.rs[i] = (o-nrow(bootsp))/(ncol(bootsp)*nrow(bootsp)-nrow(bootsp))     
    }
    hist(BC.rs)                                                              # hist the BC results      
    meanBC<-mean(BC.rs,na.rm=TRUE)                                           # calculate mean value of subsample       
    SDBC<-sd(BC.rs,na.rm=TRUE)
    lowCIBC<-quantile(BC.rs,0.025,na.rm=TRUE)
    upCIBC<-quantile(BC.rs,0.975,na.rm=TRUE)                                        
    sur.resl<-c(meanBC,SDBC,lowCIBC,upCIBC,nrow(bcdata),ncol(bcdata))
    sur.resl
  }
  c(bootbc(survivormatrix),bootbc(newcomermatrix))
}     # BC function of survivors & newcomers

SPBC = function(k){
  stagedata = bivdata[which(bivdata$stgno == k),]
  stagegenus = levels(factor(stagedata$genus))
  genusno = length(stagegenus)
  genustype = rep(NA,genusno)
  prestagedata = bivdata[which(bivdata$stgno < k),]
  pre1stagedata = bivdata[which(bivdata$stgno == k-1),]
  
  for(i in 1:genusno){
    if (sum(which(prestagedata$genus == stagegenus[i]))>0){genustype[i] = "survivor"}
    else{genustype[i] = "newcomer"}
  }
  genustypedata = cbind.data.frame(stagegenus,genustype)
  survivorgenus = stagegenus[which(genustypedata$genustype == "survivor")]
  
  survivorno = length(survivorgenus)
  survivortype = rep(NA,survivorno)
  for(i in 1:survivorno){
    if (sum(which(pre1stagedata$genus == survivorgenus[i]))>0){survivortype[i] = "RT"}
    else{genustype[i] = "NRT"}
  }
  RTsurvivor = survivorgenus[which(survivortype == "RT")]
  
  # k-1 stage RTsurvivor matrix
  sn.row = length(RTsurvivor)
  sn.col = length(levels(factor(pre1stagedata$geocell)))
  sinitial.matrix = matrix (data = 0, sn.row,sn.col) 
  rownames(sinitial.matrix) = RTsurvivor 
  colnames(sinitial.matrix) = levels(factor(pre1stagedata$geocell))
  
  # k stage RTsurvivor matrix
  kn.row = length(RTsurvivor)
  kn.col = length(levels(factor(stagedata$geocell)))
  kinitial.matrix = matrix (data = 0, kn.row,kn.col) 
  rownames(kinitial.matrix) = RTsurvivor 
  colnames(kinitial.matrix) = levels(factor(stagedata$geocell))
  
  # Calculating k-1 stage RTsurvivor matrix
  for (i in 1:nrow(pre1stagedata)){
    gname = levels(factor(pre1stagedata$genus[i]))
    cname = levels(factor(pre1stagedata$geocell[i]))
    sinitial.matrix[which(rownames(sinitial.matrix) == gname),which(colnames(sinitial.matrix) == cname)] = 1
  }
  K1RTsurvivormatrix = sinitial.matrix
  ncol(K1RTsurvivormatrix)
  scolsum = as.matrix(colSums(K1RTsurvivormatrix))
  if(sum(which(scolsum == 0))>0)
  {K1RTsurvivormatrix = K1RTsurvivormatrix[,-which(scolsum == 0)]}
  
  # Calculating k stage RTsurvivor matrix
  for (i in 1:nrow(stagedata)){
    gname = levels(factor(stagedata$genus[i]))
    cname = levels(factor(stagedata$geocell[i]))
    kinitial.matrix[which(rownames(kinitial.matrix) == gname),which(colnames(kinitial.matrix) == cname)] = 1
  }
  KRTsurvivormatrix = kinitial.matrix
  scolsum = as.matrix(colSums(KRTsurvivormatrix))
  if(sum(which(scolsum == 0))>0)
  {KRTsurvivormatrix = KRTsurvivormatrix[,-which(scolsum == 0)]}
  
  
  # Calculating BC
  bootbc = function(n){
    bcdata = n
    BC.rs = rep(NA,10000)
    for(i in 1:10000){
      ssample<-sample(nrow(bcdata),size=nrow(bcdata),replace=TRUE)
      bootsp<-bcdata[ssample,]
      subsamcolsum = as.matrix(colSums(bootsp))
      if(sum(which(subsamcolsum == 0))>0)
      {bootsp = bootsp[,-which(subsamcolsum == 0)]}
      o<-sum(bootsp)
      BC.rs[i] = (o-nrow(bootsp))/(ncol(bootsp)*nrow(bootsp)-nrow(bootsp))
    }
    hist(BC.rs)
    meanBC<-mean(BC.rs,na.rm=TRUE)
    SDBC<-sd(BC.rs,na.rm=TRUE)
    lowCIBC<-quantile(BC.rs,0.025,na.rm=TRUE)
    upCIBC<-quantile(BC.rs,0.975,na.rm=TRUE)
    sur.resl<-c(meanBC,SDBC,lowCIBC,upCIBC,nrow(bcdata),ncol(bcdata))
    sur.resl
  }
  bootbc(K1RTsurvivormatrix)
}     # BC function of survivor in previous stage

subBC = function(p,k){
  stgocc = occdata2[which(occdata2$stgno == p),]
  n.row = length(levels(factor(stgocc$genus))) #number of genus in a stage
  n.col = length(levels(factor(stgocc$geocell))) #number of locations in a stage
  if (n.row < 10 | n.col<5){
    return(rr <<- 0)}
  else{rr <<- 1}
  
  # Stage matrix
  
  a = nrow(stgocc) 
  levels(factor(stgocc$genus))  
  n.row = length(levels(factor(stgocc$genus))) 
  levels(factor(stgocc$geocell)) 
  n.col = length(levels(factor(stgocc$geocell))) 
  initial.matrix = matrix (data = 0, n.row,n.col) 
  rownames(initial.matrix) = levels(factor(stgocc$genus)) 
  colnames(initial.matrix) = levels(factor(stgocc$geocell)) 
  initial.matrix
  
  # Calculating matrix
  for (i in 1:a){
    gname = levels(factor(stgocc$genus[i]))
    cname = levels(factor(stgocc$geocell[i]))
    initial.matrix[which(rownames(initial.matrix) == gname),which(colnames(initial.matrix) == cname)] = 1
  }
  
  bcdata = initial.matrix
  
  
  # Subsample
  
  BC.rs = rep(NA,10000)
  b = rep(NA,n.col)
  for(i in 1:10000){
    ssample<-sample(n.col,size=k,replace=FALSE)
    bootsp<-bcdata[,ssample]
    subsamrowsum = as.matrix(rowSums(bootsp))
    if(sum(which(subsamrowsum == 0))>0)
    {bootsp = bootsp[-which(subsamrowsum == 0),]}
    o<-sum(bootsp)
    BC.rs[i] = (o-nrow(bootsp))/(ncol(bootsp)*nrow(bootsp)-nrow(bootsp))
  }
  hist(BC.rs)
  BC.rs
}  # Subsample BC by equal cells function

MSBD = function(p,k){
  stgocc = occdata2[which(occdata2$stgno == p),]
  n.row = length(levels(factor(stgocc$genus))) #number of genus in a stage
  n.col = length(levels(factor(stgocc$geocell))) #number of locations in a stage
  if (n.row < 10 | n.col<5){
    return(rr <<- 0)}
  else{rr <<- 1}  
  
  # Creating geocell-taxa metrix for a stage
  
  a = nrow(stgocc) 
  levels(factor(stgocc$genus))  #genera in a stage
  n.row = length(levels(factor(stgocc$genus))) #number of genus in a stage
  n.col = length(levels(factor(stgocc$geocell))) #number of locations in a stage
  initial.matrix = matrix (data = 0, n.row,n.col) 
  rownames(initial.matrix) = levels(factor(stgocc$genus)) 
  colnames(initial.matrix) = levels(factor(stgocc$geocell)) 
  
  
  for (i in 1:a){
    gname = levels(factor(stgocc$genus[i]))
    cname = levels(factor(stgocc$geocell[i]))
    initial.matrix[which(rownames(initial.matrix) == gname),which(colnames(initial.matrix) == cname)] = 1
  }
  
  # Calculating miltiple-site beta diversity
  
  bcdata = initial.matrix 
  mbetamartix = t(bcdata)
  jbeta = rep(NA,10000)
  for(i in 1:10000){
    ssample<-sample(nrow(mbetamartix),size=k,replace=FALSE)
    bootsp<-mbetamartix[ssample,]
    subsamcolsum = as.matrix(colSums(bootsp))
    if(sum(which(subsamcolsum == 0))>0)      #remove taxa without occurrences in subsamples
    {bootsp = bootsp[,-which(subsamcolsum == 0)]}
    
    else{bootsp = bootsp}
    
    jacbeta = beta.multi(bootsp,index.family = "jaccard")
    jbeta[i] = jacbeta$beta.JAC
  }
  jbeta
}   # Subsample beta by equal cells function






