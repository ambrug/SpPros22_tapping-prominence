# anna.bruggeman@uni-bielefeld.de

## SUMMARY of script:
# 1) data cleaning/preprocessing
# 2) compute inter-participant agreement on RPT task
# 3)  - RandomForests on taping data by participant, then for full dataset, with function cforest()
#     - varimp plots
# 4)  - RandomForests on mean RPT, and on mean tapping data, with function ranger()
#     - varimp plots

## INPUT files
# alldf_groupUI.csv (created in tap_vs_RPT.R)
# wordinfo_manual.csv
# SUBTLEX-DE_cleaned_with_Google00.txt. Retrievable here: http://crr.ugent.be/SUBTLEX-DE/SUBTLEX-DE_txt_cleaned_with_Google00.zip

## OUTPUT files
# n/a

rm(list=ls())

library(lme4)
library(irr)
library(party)
library(ranger)
library(ggplot2)
library(randomForest)
library(tidyverse)
library(reshape2)

## add two manual functions
`%nin%` = Negate(`%in%`)
RSQUARE = function(y_actual,y_predict){cor(y_actual,y_predict)^2}

########## LOAD DATA, specify data types
alldf <- read.csv("alldf_groupUI.csv")[,-1]

# check class of columns
sapply(alldf, class)
colnames(alldf)

# to become character/factors
char2fac <- colnames(alldf)[c(1,2,4,5,6,7,14,16,17,25,37:44)]
alldf[char2fac] <- sapply(alldf[char2fac],as.factor)
alldf$p <- as.factor(alldf$p)
alldf$PA <- as.factor(alldf$PA)

# to become numeric
char2num <- colnames(alldf)[c(3,8:13,15,18:24,26:36,45:48)]
alldf[char2num] <- sapply(alldf[char2num],as.numeric) # warnings --> NAs introduced

##########################################
##############  PREPROCES acoustics: remove outliers
alldf[which(alldf$intens_v < 50),]$intens_v <- NA
## F0
alldf[which(alldf$filename == "M-A_left_0613" & alldf$syll == "C"),]$f0meanst <- NA
alldf[which(alldf$filename == "M-A_left_0613" & alldf$syll == "C"),]$f0maxst <- NA
alldf[which(alldf$filename == "M-A_left_0613" & alldf$syll == "C"),]$f0rangest <- NA
alldf[which(alldf$filename == "M-A_left_0613" & alldf$syll == "C"),]$f0meanhz <- NA
alldf[which(alldf$filename == "M-A_left_0613" & alldf$syll == "C"),]$f0maxhz <- NA
alldf[which(alldf$filename == "M-A_left_0613" & alldf$syll == "C"),]$f0rangehz <- NA
## CPPS - OK
boxplot(alldf$cpps)
## H1H2 - OK
boxplot(alldf$h1h2)
alldf[which(alldf$h1h2 > 30),]$h1h2 <- NA
## SLOPE - OK
boxplot(alldf$slope)

########################################################################
######### Inter-participant agreement in RPT-underlining 
########################################################################
# long to wide
pdata <- alldf[,c(2,3:6)]
levels(pdata$p) <- c(0,1)
pdata$p <- as.numeric(as.character(pdata$p))
pdata$tapper <- as.factor(pdata$tapper)
pdata$filename <- as.factor(pdata$filename)
pdata$syll <- as.factor(pdata$syll)
pdatawide <- dcast(pdata, filename+syll+syllno ~ tapper, value.var="p")
kappam.fleiss(pdatawide[4:13])# 0.505 for VP1-12
#########################
#########################

#########################
######### PA vs p correspondence for individuals: X2 and Cramér's V
#########################
alldf$tapper<- as.factor(alldf$tapper)
# create empty df
chisqtable <- data.frame(matrix(vector(), 0, 0))

# create an empty list 
dfs_participant <- list()

# add individual participant DFs to this list
for (i in 1:length(levels(alldf$tapper))){
  partdf <- alldf[alldf$tapper==levels(as.factor(alldf$tapper))[i],]
  dfs_participant[[i]] <- partdf
}

# run through list and select each participant's dataframe in turn
for (i in 1:length(dfs_participant)){
  df <- dfs_participant[[i]]
  tapper <- paste0("P",i)
  
  # calc values of Chisq table cells
  prom_PA <- sum(df$p == "prom" & df$PA == "pitch accented")
  prom_noPA <- sum(df$p == "prom" & df$PA != "pitch accented")
  noprom_PA <- sum(df$p != "prom" & df$PA == "pitch accented")
  noprom_noPA <- sum(df$p != "prom" & df$PA != "pitch accented")
  
  # create 2x2 matrix
  tbl = matrix(data=c(prom_PA, prom_noPA, noprom_PA, noprom_noPA), nrow=2, ncol=2, byrow=T)
  dimnames(tbl) = list(Prominence=c('prom', 'not prom'), Accentcategory=c('PA', 'no PA'))
  
  chi2 = chisq.test(tbl, correct=F)
  pvalue <- chi2$p.value # X2 and p value
  Vvalue <- sqrt(chi2$statistic / sum(tbl)) # Cramer's V X2)
  chisqtable <- rbind(chisqtable, c(Vvalue,pvalue,tapper))
}
colnames(chisqtable) <- c("cramerV", "pvalue", "participant")

range(chisqtable$cramerV) #0.38 - 0.75
mean(as.numeric(chisqtable$cramerV)) #0.58
#########################
#########################


###########################################################################
#### TRANSFORM / Z-SCORE acoustic variables by sentence
alldf <- alldf %>%
  group_by(filename) %>%
  mutate(sylldur_z = scale(sylldur),
         segdur_z = scale(segdur),
         intens_v_z = scale(intens_v),
         intens_z = scale(intens),
         f0meanst_v_z = scale(f0meanst_v),
         f0maxst_v_z = scale(f0maxst_v),
         f0rangest_v_z = scale(f0rangest_v),
         f0meanst_z = scale(f0meanst),
         f0maxst_z = scale(f0maxst),
         f0rangest_z = scale(f0rangest),
         cpps_z = scale(cpps),
         h1h2_z = scale(h1h2),
         slope_z = scale(slope),
  )

###############################################################
############## Random Forests
###############################################################

## add info à la Baumann & Winter 2018: POS, WordFreq, nsyll of word
wordinfo <- read.csv("wordinfo_manual.csv", 
                     sep=";", header=T)
wordinfo$word <- as.factor(wordinfo$word)
wordinfo$phrase.pos <- as.factor(wordinfo$phrase.pos)

subtlex <- read.table("SUBTLEX-DE_cleaned_with_Google00.txt", 
                      sep="\t", header=T, fileEncoding="iso-8859-1", quote="\"")
subtlex <- subtlex[,c(1,6)]

subtlex <- subtlex %>% 
  mutate_all(funs(str_replace(., ",", ".")))

# merge subtlex with wordinfo (words are manually changed to match with subtlex)
colnames(subtlex)[1] <- "word"
wordinfo <- merge(x = wordinfo, y = subtlex, by = "word", all.x = TRUE)

# keep only filename and syllno - thats will allow for correct merge
wordinfo$UI <- paste(wordinfo$filename, wordinfo$syllno, wordinfo$syll, sep="-")
wordinfo$filenamesyll <- as.factor(paste(wordinfo$filename, wordinfo$syllno, sep="-"))

# subset necessary columns
wordinfo <- wordinfo[,c(1,2,4:10)]

# nrows identical between alldf and alldf_merged. OK
alldf_merged <- merge(alldf, wordinfo, by=c("filename", "syllno", "syll"), all.x=T)
#create backup
alldf_premerge <- alldf
# merge wordinfo with existing alldf
alldf <- alldf_merged
unique(as.factor(alldf[which(is.na(alldf$lgSUBTLEX)),]$word)) # 13 words have no logFreq 

# further proc
alldf$POS <- as.factor(alldf$POS)
alldf$lgSUBTLEX <- as.factor(alldf$lgSUBTLEX)
alldf$stressable <- as.factor(alldf$stressable)

# save memory
rm(subtlex)


####################################################################################
############################ RANDOM FORESTS with CFOREST - TAPPING data only
####################################################################################
# create new response variable with either TapDur or TapForce

# P3/ZE and P11/RU get duration, the rest velocity=force
alldf$predictee <- alldf$DrumVeloc.log.sc
summary(alldf$predictee) #1832 NAs

alldf[alldf$tapper %in% c("ZE", "RU"),]$predictee <- NA
summary(alldf$predictee) #1475 NAs

for (z in 1:nrow(alldf)) {
  if (is.na(alldf$predictee[z])){
    alldf$predictee[z] <- alldf$DrumDur.log.sc[z]
  }
}
summary(alldf$predictee) #1832 NAs OK!
# create backup
alldf_full <- alldf
alldf <- alldf[!is.na(alldf$predictee),]

#################################
## set up list with participants for iteration
alldf$tapper <- as.factor(alldf$tapper)

# create an empty list
dfs_participant <- list()

# add individual participant DFs to this list
for (i in 1:length(levels(alldf$tapper))){
  partdf <- alldf[alldf$tapper==levels(as.factor(alldf$tapper))[i],]
  # add item to list which is this particular iteration's df
  dfs_participant[[i]] <- partdf
}

##################################################################################################################
###################################### RF by participant: ACOUSTIC + MIXED. Takes a very long time.
# create empty final df
varimp_cf1 <- data.frame(matrix(vector(), 0, 12))

#\dontrun{
for (i in 1:length(dfs_participant)){
  
  df <- dfs_participant[[i]]
  tapper_anon <- paste0("P",i)
  
  # also remove f0rangest <- if not in there then some other values also wont be
  df <- df[!(is.na(df$f0rangest_v)),]
  
  # need to remove NAs from predicted value
  df <- df[!(is.na(df$predictee)),]
  
  # empty temp df
  rf_temp <- data.frame(matrix(vector(), 0, 0))
  
  # run 100 iterations
  for (x in 1:100) {
    set.seed(x)
    
    # creation of random test/validation set depending on fixed seed
    trainnum <- sample(1:nrow(df), 0.7*nrow(df), replace=F) # training 70% of data
    df_train <- df[c(trainnum),]
    df_test <- df[-c(trainnum),]
    
    ####### ACOUSTIC ONLY
    #cforest_acoustic <- cforest(predictee ~ sylldur_z + segdur_z +
    #                                 f0meanst_v_z + f0maxst_v_z + 
    #                                 f0rangest_v_z + 
    #                                 intens_v_z +
    ##                                 h1h2_z + cpps_z + slope_z, 
    #                              data = df_train,
    #                               controls = cforest_unbiased(),
    #                               xtrafo = ptrafo, ytrafo = ptrafo, #defaults
    #                             scores = NULL, weights = NULL)
    #   pred <- predict(cforest_acoustic, newdata = df_test, type="response", OOB = TRUE)
    #   rsq <- RSQUARE(df_test$predictee, pred)
    #   rmsevalue <- rmse(df_test$predictee, pred)
    #   varimptemp1 <- unname(varimp(cforest_acoustic))
    #   columnames <- c(names(varimp(cforest_acoustic)), "r2", "rmse")
    
    ######## MIXED
    cforest_mixed <- cforest(predictee ~ sylldur_z + segdur_z +
                               slope_z +
                               cpps_z +
                               PA + POS + 
                               lgSUBTLEX + nsyllword +
                               phrase.pos,
                             data=df_train,
                             controls = cforest_unbiased(),
                             xtrafo = ptrafo, ytrafo = ptrafo, #defaults
                             scores = NULL, weights = NULL)
    pred <- predict(cforest_mixed, newdata = df_test, type="response", OOB = TRUE)
    rsq <- RSQUARE(df_test$predictee, pred)
    rmsevalue <- rmse(df_test$predictee, pred)
    varimptemp1 <- unname(varimp(cforest_mixed))
    
    # add 1 line to tempdf for this forest's var. rankings and rmse, r2
    rf_temp[x, 1] <- varimptemp1[1]
    rf_temp[x, 2] <- varimptemp1[2] 
    rf_temp[x, 3] <- varimptemp1[3]
    rf_temp[x, 4] <- varimptemp1[4] 
    rf_temp[x, 5] <- varimptemp1[5] 
    rf_temp[x, 6] <- varimptemp1[6]
    rf_temp[x, 7] <- varimptemp1[7] 
    rf_temp[x, 8] <- varimptemp1[8] 
    rf_temp[x, 9] <- varimptemp1[9]
    rf_temp[x, 10] <- rsq
    rf_temp[x, 11] <- rmsevalue
    rf_temp[x, 12] <- paste0("P", i)
    
    # for x in 100
  }
  
  # append mean values for each participant to new dataframe
  varimp_cf1[i,1] <- mean(rf_temp[,1])
  varimp_cf1[i,2] <- mean(rf_temp[,2])
  varimp_cf1[i,3] <- mean(rf_temp[,3])
  varimp_cf1[i,4] <- mean(rf_temp[,4])
  varimp_cf1[i,5] <- mean(rf_temp[,5])
  varimp_cf1[i,6] <- mean(rf_temp[,6])
  varimp_cf1[i,7] <- mean(rf_temp[,7])
  varimp_cf1[i,8] <- mean(rf_temp[,8])
  varimp_cf1[i,9] <- mean(rf_temp[,9])
  varimp_cf1[i,10] <- mean(rf_temp[,10])
  varimp_cf1[i,11] <- mean(rf_temp[,11])
  varimp_cf1[i,12] <- paste0("P", i)
  # for dfs/participant
  
}
# dontrun
#}

################################################
############ post-proc ACOUSTIC forests
################################################
# set colnames of empty dataframe
colnames(varimp_cf1) <- c(names(varimp(cforest_acoustic)), "r2", "rmse", "participant")

varimpVP <- varimp_cf1
mean(varimpVP$r2) 
range(varimpVP$r2) 

#wide to long
varimpVP <- gather(varimpVP, predictor, importance, sylldur_z:slope_z, factor_key=TRUE)
# arrange highest to lower importance
varimpVP <- varimpVP[order(varimpVP$importance, decreasing=F),]
colnames(varimpVP)[colnames(varimpVP) %in% "speaker"] <- "participant"

# NEED TO GROUP BY particiapnt
varimpVP$participant <- as.factor(varimpVP$participant)
# get top 3 variables pp
varimpVP_summary <- varimpVP %>%                         
  arrange(desc(importance)) %>% 
  group_by(participant) %>%
  slice(1:3)

################################################
############ post-proc MIXED forests
################################################
# set colnames of empty dataframe
colnames(varimp_cf1) <- c(names(varimp(cforest_mixed)), "r2", "rmse", "participant")
colnames(varimpVPmixed) <- c(names(varimp(cforest_mixed)), "r2", "rmse", "participant")

#wide to long
varimpVPmixed <- gather(varimpVPmixed, predictor, importance, sylldur_z:phrase.pos, factor_key=TRUE)
# arrange highest to lower importance
varimpVPmixed <- varimpVPmixed[order(varimpVPmixed$importance, decreasing=F),]

range(varimpVPmixed$r2) # 0.01 to 0.30
mean(varimpVPmixed$r2) # 0.11
# GROUP BY particiapnt
varimpVPmixed %>% 
  group_by(participant) %>%
  arrange(order(importance))

#################################################################### 
############ RFs ON FULL DATASET, ALL PARTICIPANTS TOGETHER
####################################################################
alldf_complete <- alldf[!(is.na(alldf$f0rangest_v)),] # N=5093 left. still get NA messages
sum(is.na(alldf_complete$predictee)) # f0 mean max, durs, VQs
alldf_complete <- alldf_complete[!(is.na(alldf_complete$intens_v)),] # N=3512

####################################################################
############################ ACOUSTIC

# empty final df # acoustic
rf_temp <- data.frame(matrix(vector(), 0, 0))

#\dontrun{
for (x in 1:100) {
  set.seed(x)
  
  trainnum <- sample(1:nrow(alldf_complete), 0.7*nrow(alldf_complete), replace=F) # training 70% of data
  df_train <- alldf_complete[c(trainnum),]
  df_test <- alldf_complete[-c(trainnum),]
  
  # RF on DrumForce with just acoustic predictors
  cforest_acoustic <- cforest(predictee ~ sylldur_z + segdur_z +
                                f0meanst_v_z + f0maxst_v_z + 
                                f0rangest_v_z + 
                                intens_v_z +
                                h1h2_z + cpps_z + slope_z, 
                              data = df_train,
                              controls = cforest_unbiased(),
                              xtrafo = ptrafo, ytrafo = ptrafo, #defaults
                              scores = NULL, weights = NULL)
  
  pred <- predict(cforest_acoustic, newdata = df_test, type="response", OOB = TRUE)
  rsq <- RSQUARE(df_test$predictee, pred)
  rmsevalue <- rmse(df_test$predictee, pred)
  varimptemp1 <- unname(varimp(cforest_acoustic))
  
  # add to temp dataframe for this x iteration
  rf_temp[x, 1] <- varimptemp1[1]
  rf_temp[x, 2] <- varimptemp1[2]
  rf_temp[x, 3] <- varimptemp1[3]
  rf_temp[x, 4] <- varimptemp1[4]
  rf_temp[x, 5] <- varimptemp1[5]
  rf_temp[x, 6] <- varimptemp1[6]
  rf_temp[x, 7] <- varimptemp1[7]
  rf_temp[x, 8] <- varimptemp1[8]
  rf_temp[x, 9] <- varimptemp1[9]
  rf_temp[x, 10] <- rsq
  rf_temp[x, 11] <- rmsevalue
  
  # for x in 100
}

varimp_all <- data.frame(matrix(vector(), 0, 12))
######################
### post-proc ACOUSTIC
######################
# append mean values for to new dataframe with 1 line
varimp_all[1,1] <- mean(rf_temp[,1])
varimp_all[1,2] <- mean(rf_temp[,2])
varimp_all[1,3] <- mean(rf_temp[,3])
varimp_all[1,4] <- mean(rf_temp[,4])
varimp_all[1,5] <- mean(rf_temp[,5])
varimp_all[1,6] <- mean(rf_temp[,6])
varimp_all[1,7] <- mean(rf_temp[,7])
varimp_all[1,8] <- mean(rf_temp[,8])
varimp_all[1,9] <- mean(rf_temp[,9]) #last predictor cpps
varimp_all[1,10] <- mean(rf_temp[,10]) #r2
varimp_all[1,11] <- mean(rf_temp[,11]) #rmse
varimp_all[1,12] <- "acoustic"
colnames(varimp_all) <- c(names(varimp(cforest_acoustic)), "r2", "rmse", "type")

####################################################################
############################ MIXED
# empty final df # mixed 
rf_temp2 <- data.frame(matrix(vector(), 0, 0))

for (x in 1:100) {
  set.seed(x)
  
  trainnum <- sample(1:nrow(alldf_complete), 0.7*nrow(alldf_complete), replace=F) # training 70% of data
  df_train <- alldf_complete[c(trainnum),]
  df_test <- alldf_complete[-c(trainnum),]
  
  # RF on DrumForce with just acoustic predictors
  cforest_mixed <- cforest(predictee ~ sylldur_z + segdur_z +
                             slope_z +
                             cpps_z +
                             PA + POS + 
                             lgSUBTLEX + nsyllword +
                             phrase.pos,
                           data = df_train,
                           controls = cforest_unbiased(),
                           xtrafo = ptrafo, ytrafo = ptrafo, #defaults
                           scores = NULL, weights = NULL)
  
  pred2 <- predict(cforest_mixed, newdata = df_test, type="response", OOB = TRUE)
  rsq2 <- RSQUARE(df_test$predictee, pred2)
  rmsevalue2 <- rmse(df_test$predictee, pred2)
  varimptemp2 <- unname(varimp(cforest_mixed))
  # add to temp dataframe for this x iteration
  rf_temp2[x, 1] <- varimptemp2[1]
  rf_temp2[x, 2] <- varimptemp2[2]
  rf_temp2[x, 3] <- varimptemp2[3]
  rf_temp2[x, 4] <- varimptemp2[4]
  rf_temp2[x, 5] <- varimptemp2[5]
  rf_temp2[x, 6] <- varimptemp2[6]
  rf_temp2[x, 7] <- varimptemp2[7]
  rf_temp2[x, 8] <- varimptemp2[8]
  rf_temp2[x, 9] <- varimptemp2[9]
  rf_temp2[x, 10] <- rsq2
  rf_temp2[x, 11] <- rmsevalue2
  colnames(rf_temp2) <- c(names(varimp(cforest_mixed)), "r2", "rmse")
}
mean(rf_temp2$r2)


######################
### post-proc MIXED
######################
varimp_mixed <- data.frame(matrix(vector(), 0, 12))

# append mean values for to new MIXED dataframe with 1 line
varimp_mixed[1,1] <- mean(rf_temp2[,1])
varimp_mixed[1,2] <- mean(rf_temp2[,2])
varimp_mixed[1,3] <- mean(rf_temp2[,3])
varimp_mixed[1,4] <- mean(rf_temp2[,4])
varimp_mixed[1,5] <- mean(rf_temp2[,5])
varimp_mixed[1,6] <- mean(rf_temp2[,6])
varimp_mixed[1,7] <- mean(rf_temp2[,7])
varimp_mixed[1,8] <- mean(rf_temp2[,8])
varimp_mixed[1,9] <- mean(rf_temp2[,9]) 
varimp_mixed[1,10] <- mean(rf_temp2[,10]) #r2
varimp_mixed[1,11] <- mean(rf_temp2[,11]) #rmse
varimp_mixed[1,12] <- "mixed"
colnames(varimp_mixed) <- c(names(varimp(cforest_mixed)), "r2", "rmse", "type")

## SUMMARY: means for entire POS dataframe: 0.07 r2 and range between participants
#mean(varimp_cf1$r2) #0.10, 
#mean(varimp_cf1$r2)range 0.00 to 0.27
########################################## 
########################################## END OF RFs on FULL DATASET


########################################## 
########################################## PLOT VARIMP FULL DATASET
## PLOT 
varimp_all <- varimp_all[,-1]
varimp_all <- varimp_all[,-12]
varimp_mixed <- varimp_mixed[,-1]
varimp_mixed <- varimp_mixed[,-12]
#wide to long
varimp_all <- gather(varimp_all, predictor, importance, sylldur_z:slope_z, factor_key=TRUE)
varimp_mixed <- gather(varimp_mixed, predictor, importance, sylldur_z:phrase.pos, factor_key=TRUE)
# arrange highest to lower importance
varimp_all <- varimp_all[order(varimp_all$importance, decreasing=F),]
varimp_mixed <- varimp_mixed[order(varimp_mixed$importance, decreasing=F),]
# reorder predictor levels: a la factor(sizes, levels = c("small", "medium", "large"))
varimp_all$predictor <- factor(varimp_all$predictor, levels=varimp_all$predictor)
varimp_mixed$predictor <- factor(varimp_mixed$predictor, levels=varimp_mixed$predictor)

varimp <- rbind(varimp_all, varimp_mixed)
# all lower case
varimp$type <- as.factor(varimp$type)
levels(varimp_all$predictor) <- c("f0 mean", "intens", "f0 max", "f0 range", "h1h2", "cpps", "alpha", "dur v", "dur syll")
####################################
## plot variable importance ACOUSTIC
####################################
ggplot() +
  geom_point(data=varimp_all,
             aes(x=importance, y=predictor), size=3) + 
  scale_x_continuous("Relative importance") +
  scale_y_discrete("") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20),
        legend.position="none")

####################################
## plot variable importance MIXED
####################################
# UPPER CASE for non-acoustic predictors
levels(varimp_mixed$predictor) <- c("N SYLL", "POS", "PHRASE-POS", "alpha", "cpps", "dur v", "FREQ", "dur syll", "PA")
ggplot() +
  geom_point(data=varimp_mixed,
             aes(x=importance, y=predictor), size=3) + 
  scale_x_continuous("Relative importance") +
  scale_y_discrete("") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20),
        legend.position="none")
####################################################################################


####################################################################################
############################ RANDOM FORESTS with RANGER: on RPT data
####################################################################################

############## remove incomplete observations
alldf_completep <- alldf[!is.na(alldf$slope_z),] 
alldf_completep <- alldf_completep[!is.na(alldf_completep$POS),] 
alldf_completep <- alldf_completep[!is.na(alldf_completep$lgSUBTLEX),] ##3320

#################################
################################# preproces 

## need to get mean p-scores (each syll then will have only 1 value)
alldf_completep$p <- as.character(alldf_completep$p)
alldf_completep[alldf_completep$p %in% "not prom",]$p <- 0
alldf_completep[alldf_completep$p == "prom",]$p <- 1
alldf_completep$p <- as.numeric(alldf_completep$p)

# new dataframe with N=403 mean p-scores
pdf <- alldf_completep %>%
  group_by(filename, syllno) %>%
  summarise(meanp = mean(p))

# create subset qith relevant columns for RF
ydf <- alldf_completep[,colnames(alldf_completep) %in% c("filename", "syllno", "sylldur_z", "segdur_z", "slope_z", "cpps_z", "h1h2_z", "f0meanst_v_z", "f0maxst_v_z", "f0rangest_v_z", "intens_z", "PA", "POS", "lgSUBTLEX", "nsyllword", "phrase.pos")]
# remove duplicate rows
ydf <- unique(ydf) #403 same as pdf. 
# merge syll values with mean p values
pdf <- merge(pdf, ydf, by=c("filename", "syllno"))

pdf <- pdf[!is.na(pdf$f0rangest_v_z),]
pdf <- pdf[!is.na(pdf$f0meanst_v_z),]
pdf <- pdf[!is.na(pdf$f0maxst_v_z),]
pdf <- pdf[!is.na(pdf$h1h2_z),]
# N=386 left


# ALL ACOUSTIC + NON-ACOUSTIC
myFormula_mixed <- as.formula('meanp ~ sylldur_z + segdur_z +
                             slope_z +
                             cpps_z +
                             POS + 
                             PA +
                             lgSUBTLEX + nsyllword +
                             phrase.pos')

# ALL ACOUSTIC + NON-ACOUSTIC - PA
#myFormula <- as.formula('meanp ~ sylldur_z + segdur_z +
#                             slope_z +
#                             cpps_z +
#                             f0meanst_v_z +
#                             f0rangest_v_z +
#                             f0maxst_v_z +
#                            intens_z +
#                            POS + 
#                             lgSUBTLEX + nsyllword +
#                             phrase.pos')

#ACOUSTIC
myFormula_ac <- as.formula('meanp ~ sylldur_z + segdur_z +
                             slope_z +
                             cpps_z +
                              h1h2_z +
                             f0meanst_v_z +
                             f0rangest_v_z +
                             f0maxst_v_z +
                             intens_z')

# QUICK RUN FOR BALLPARK NUMBERS
# MIXED predictors 0.58
myforest_red <- ranger(myFormula_mixed, data = pdf,
                       importance = 'permutation', mtry = 3, num.trees = 2000)
myforest_red$r.squared
# acoustic only 0.23
myforest_red <- ranger(myFormula_ac, data = pdf,
                       importance = 'permutation', mtry = 3, num.trees = 2000)
myforest_red$r.squared 

#############################################################
############# - MEAN RPT scores - based on MIXED PREDICTORS
#############################################################
# set up empty df
rangerdf_meanp_mixed <- data.frame(matrix(vector(), 0, 2))
for (x in 1:100) {
  set.seed(x)
  # create new random training and testing data
  trainnump <- sample(1:nrow(pdf), 0.7*nrow(pdf), replace=F) # training 70% of data
  df_train_p <- pdf[c(trainnump),]
  df_test_p <- pdf[-c(trainnump),]
  
  # train
  myforest_red_train_mixed <- ranger(myFormula_mixed, data = df_train_p,
                                     importance = 'permutation', mtry = 3, num.trees = 2000)
  # get predictions on test data
  myforest_red_preds_mixed <- predict(myforest_red_train_mixed,
                                      data = df_test_p)
  #r2 values: 
  multr2 <- summary(lm(scale(myforest_red_preds_mixed$predictions) ~
                         -1 + scale(df_test_p$meanp)))$r.squared
  rangerdf_meanp_mixed <- rbind(rangerdf_meanp_mixed, c(multr2, x))
}

#############################################################
############# - MEAN RPT scores - based on ACOUSTIC PREDICTORS
#############################################################
# set up empty df
rangerdf_meanp_ac <- data.frame(matrix(vector(), 0, 2))
for (x in 1:100) {
  set.seed(x)
  # create new random training and testing data
  trainnump <- sample(1:nrow(pdf), 0.7*nrow(pdf), replace=F) # training 70% of data
  df_train_p <- pdf[c(trainnump),]
  df_test_p <- pdf[-c(trainnump),]
  
  # train
  myforest_red_train_ac <- ranger(myFormula_ac, data = df_train_p,
                                  importance = 'permutation', mtry = 3, num.trees = 2000)
  # get predictions on test data
  myforest_red_preds_ac <- predict(myforest_red_train_ac,
                                   data = df_test_p)
  #r2 values: 
  multr2 <- summary(lm(scale(myforest_red_preds_ac$predictions) ~
                         -1 + scale(df_test_p$meanp)))$r.squared
  rangerdf_meanp_ac <- rbind(rangerdf_meanp_ac, c(multr2, x))
}

#####################
##### RESULTS
#####################
## MIXED FORESTS # R2=0.57
colnames(rangerdf_meanp_mixed) <- c("R2", "number_run")
mean(rangerdf_meanp_mixed$R2) 

## ACOUSTIC FORESTS # R2=0.23
colnames(rangerdf_meanp_ac) <- c("R2", "number_run")
mean(rangerdf_meanp_ac$R2) 



######################################################
##### PLOTTING VARIMP mean RPT based on Ranger forests
######################################################
### sort byimportance
varimp_BW <- as.data.frame(sort(myforest_red_train_mixed$variable.importance)) 
varimp_BW_ac <- as.data.frame(sort(myforest_red_train_ac$variable.importance))

# order DF based on importance
varimp_BW$var <- rownames(varimp_BW)
varimp_BW_ac$var <- rownames(varimp_BW_ac)

colnames(varimp_BW)[1] <- "importance"
colnames(varimp_BW_ac)[1] <- "importance"

varimp_BW_ordered <- varimp_BW %>%
  arrange(importance)
varimp_BW_ac_ordered <- varimp_BW_ac %>%
  arrange(importance)

varimp_BW_ordered$var <- factor(varimp_BW_ordered$var, levels = varimp_BW_ordered$var)
varimp_BW_ac_ordered$var <- factor(varimp_BW_ac_ordered$var, levels = varimp_BW_ac_ordered$var)

# rename levels
levels(varimp_BW_ordered$var) <- c("PHRASE-POS", "N SYLL", "POS", "cpps", "FREQ", "alpha", "dur syll", "dur v", "PA")
levels(varimp_BW_ac_ordered$var) <- c("f0 range", "h1h2", "intens", "f0 max", "f0 mean", "cpps", "alpha", "dur syll", "dur v")

###########################
##### PLOT MIXED predictors
########################### 
ggplot(data=varimp_BW_ordered, aes(x=importance, y=var)) +
  geom_point(size=3) +
  scale_y_discrete("") +
  scale_x_continuous("") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20),
        legend.position="none")

###########################
##### PLOT ACOUSTIC predictors
########################### 
ggplot(data=varimp_BW_ac_ordered, aes(x=importance, y=var)) +
  geom_point(size=3) +
  scale_y_discrete("") +
  scale_x_continuous("") +
  theme_bw() +
  theme(axis.text.y = element_text(size=20),
        legend.position="none")


####################################################################################
############################ RANDOM FORESTS with RANGER: on mean TAPPING data
####################################################################################

## ACOUSTIC + NON-ACOUSTIC predictors 
myFormula2_mixed <- as.formula('meanpredictee ~ sylldur_z + segdur_z +
                             slope_z +
                             cpps_z +
                             PA + POS + 
                             lgSUBTLEX + nsyllword +
                             phrase.pos')


## ACOUSTIC ONLY predictors 9x
myFormula2_ac <- as.formula('meanpredictee ~ sylldur_z + segdur_z +
                             slope_z +
                             cpps_z + h1h2_z +
                             f0meanst_v_z +
                             f0rangest_v_z +
                             f0maxst_v_z +
                             intens_v_z')


#without PA!
myFormula2 <- as.formula('meanpredictee ~ sylldur_z + segdur_z +
                             slope_z +
                             cpps_z +
                             POS + 
                             lgSUBTLEX + nsyllword +
                             phrase.pos')

## clean up
alldf_completepredictee <- alldf[!is.na(alldf$predictee),]
alldf_completepredictee <- alldf_completepredictee[!is.na(alldf_completepredictee$slope_z),] 
alldf_completepredictee <- alldf_completepredictee[!is.na(alldf_completepredictee$h1h2_z),] 
alldf_completepredictee <- alldf_completepredictee[!is.na(alldf_completepredictee$intens_v_z),] 
alldf_completepredictee <- alldf_completepredictee[!is.na(alldf_completepredictee$f0maxst_v_z),]
alldf_completepredictee<- alldf_completepredictee[!is.na(alldf_completepredictee$POS),] 
alldf_completepredictee <- alldf_completepredictee[!is.na(alldf_completepredictee$lgSUBTLEX),] #N=3177

## new dataframe with N=413 mean tpaping predictee values
pdf <- alldf_completepredictee %>%
  group_by(filename, syllno) %>%
  summarise(meanpredictee = mean(predictee))

# create subset qith relevant columns for RF
ydf <- alldf_completepredictee[,colnames(alldf_completepredictee) %in% c("filename", "syllno", "intens_v_z", "f0meanst_v_z", "f0rangest_v_z", "f0maxst_v_z", "sylldur_z", "segdur_z", "slope_z", "h1h2_z", "cpps_z", "PA", "POS", "lgSUBTLEX", "nsyllword", "phrase.pos")]
# remove duplicate rows
ydf <- unique(ydf) #385 same as datafram pdf
# merge syll values with mean p values
pdf <- merge(pdf, ydf, by=c("filename", "syllno"))

# QUICK RUN FOR BALLPARK NUMBERS
# MIXED predictors R=0.21
myforest_red <- ranger(myFormula2_mixed, data = pdf,
                       importance = 'permutation', mtry = 3, num.trees = 2000)
myforest_red$r.squared
# acoustic R=0.11
myforest_red <- ranger(myFormula2_ac, data = pdf,
                       importance = 'permutation', mtry = 3, num.trees = 2000)
myforest_red$r.squared 


#####################################################################
############# - MEAN TAP DUR/FORCE scores - based on MIXED PREDICTORS
#####################################################################
# set up empty df
rangerdf_meanpredictee_mixed <- data.frame(matrix(vector(), 0, 2))
for (x in 1:100) {
  set.seed(x)
  # create new random trainging and testing data
  trainnump <- sample(1:nrow(pdf), 0.7*nrow(pdf), replace=F) # training 70% of data
  df_train_p <- pdf[c(trainnump),]
  df_test_p <- pdf[-c(trainnump),]
  
  # train
  myforest_red_train_mixed <- ranger(myFormula2_mixed, data = df_train_p,
                                     importance = 'permutation', mtry = 3, num.trees = 2000)
  # get predictions on test data
  myforest_red_preds_mixed <- predict(myforest_red_train_mixed,
                                      data = df_test_p)
  #r2 values: 
  multr2 <- summary(lm(scale(myforest_red_preds_mixed$predictions) ~
                         -1 + scale(df_test_p$meanpredictee)))$r.squared
  rangerdf_meanpredictee_mixed <- rbind(rangerdf_meanpredictee_mixed, c(multr2, x))
}

#####################################################################
############# - MEAN TAP DUR/FORCE scores - based on ACOUSTIC PREDICTORS
###################################################################### set up empty df
rangerdf_meanpredictee_ac <- data.frame(matrix(vector(), 0, 2))
for (x in 1:100) {
  set.seed(x)
  # create new random trainging and testing data
  trainnump <- sample(1:nrow(pdf), 0.7*nrow(pdf), replace=F) # training 70% of data
  df_train_p <- pdf[c(trainnump),]
  df_test_p <- pdf[-c(trainnump),]
  
  # train
  myforest_red_train_ac <- ranger(myFormula2_ac, data = df_train_p,
                                  importance = 'permutation', mtry = 3, num.trees = 2000)
  # get predictions on test data
  myforest_red_preds_ac <- predict(myforest_red_train_ac,
                                   data = df_test_p)
  #r2 values: 
  multr2 <- summary(lm(scale(myforest_red_preds_ac$predictions) ~
                         -1 + scale(df_test_p$meanpredictee)))$r.squared
  rangerdf_meanpredictee_ac <- rbind(rangerdf_meanpredictee_ac, c(multr2, x))
}

#################
#### RESULTS
#################
# mixed R2= 0.22
colnames(rangerdf_meanpredictee_mixed) <- c("R2", "number_run")
mean(rangerdf_meanpredictee_mixed$R2) 

# acoustic R2= 0.12
colnames(rangerdf_meanpredictee_ac) <- c("R2", "number_run")
mean(rangerdf_meanpredictee_ac$R2) 


######################################################
##### PLOTTING VARIMP mean RPT based on Ranger forests
######################################################
varimp_BW <- as.data.frame(sort(myforest_red_train_mixed$variable.importance))
varimp_BW_ac <- as.data.frame(sort(myforest_red_train_ac$variable.importance))

# order DF based on importance
varimp_BW$var <- rownames(varimp_BW)
varimp_BW_ac$var <- rownames(varimp_BW_ac)

colnames(varimp_BW)[1] <- "importance"
colnames(varimp_BW_ac)[1] <- "importance"

varimp_BW_ordered <- varimp_BW %>%
  arrange(importance)
varimp_BW_ac_ordered <- varimp_BW_ac %>%
  arrange(importance)

varimp_BW_ordered$var <- factor(varimp_BW_ordered$var, levels = varimp_BW_ordered$var)
varimp_BW_ac_ordered$var <- factor(varimp_BW_ac_ordered$var, levels = varimp_BW_ac_ordered$var)

# rename level
levels(varimp_BW_ordered$var) <- c("cpps", "N SYLL", "POS", "dur v", "PHRASE-POS", "FREQ", "alpha", "dur syll", "PA")
levels(varimp_BW_ac_ordered$var) <- c("h1h2", "intens", "f0 range", "cpps", "f0 mean", "alpha", "dur v", "f0 max", "dur syll")


###########################
##### PLOT MIXED predictors
########################### 
ggplot(data=varimp_BW_ordered, aes(x=importance, y=var)) +
  geom_point(size=3) +
  scale_y_discrete("") +
  scale_x_continuous("") +
  theme_bw() +
  theme(axis.text.y = element_text(size=16),
        legend.position="none")

###########################
##### PLOT ACOUSTIC predictors
########################### 
ggplot(data=varimp_BW_ac_ordered, aes(x=importance, y=var)) +
  geom_point(size=3) +
  scale_y_discrete("") +
  scale_x_continuous("") +
  theme_bw() +
  theme(axis.text.y = element_text(size=16),
        legend.position="none")

