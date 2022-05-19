# anna.bruggeman@uni-bielefeld.de

## SUMMARY of script:
# 1) merge acoustic data with tapping responses, normalise per utterance
#   1b) save output dataframe for RandomForest analysis other script
# 2) get descriptive stats on tapping success
# 3) correlate RPT data with tapping, per participant

## INPUT files
# tapcorrect.txt
# tapwrong.txt
# pveloc.txt
# extracted_acoustic.txt
# GECO_stimuli_overview.csv

## OUTPUT files
# alldf_groupUI.csv -- input for RandomForest script

library(tidyverse)
library(lme4)

'%nin%' = Negate('%in%')

## READ IN DATA
# Lists of incorrect/correct tapped items
dc <- read.csv("tapcorrect.txt", header = F)
dw <- read.csv("tapwrong.txt", header = F)
# READ IN TAP + PROM Dataframe
df <- read.csv("pveloc.csv", sep=",")[,-1]
## ACOUSTICS info
ac <- read.csv("extracted_acoustic.txt")
## INFO about ambiguous syllabification, phrasal position etc.
stinfo <- read.csv("GECO_stimuli_overview.csv")[,-c(1,7)]

##### DOUBLECHK P and AC have same levels- OK
filenamesyllnosyll_p <- as.factor(paste(df$filename, df$syllno, df$syll, sep="-"))
filenamesyllnosyll_ac <- as.factor(paste(ac$filename, ac$syllno, ac$syll, sep="-"))
setdiff(unique(filenamesyllnosyll_p), filenamesyllnosyll_ac)

############################################################################################
############################################## DESCRIPTIVE STATS ON TAPPING SUCCESS
############################################################################################

## READ IN CORRECT/INCORRECT ITEMS
## crosscheck: correct and incorrectly tapped items ones are indeed exclusive
dc$UI <- paste(dc$V1, dc$V2, sep="-")
dw$UI <- paste(dw$V1, dw$V2, sep="-")
# should be 0
sum(dc$UI %in% dw$UI)
# single dataframe
dc$correct <- "y"
dw$correct <- "n"
d <- rbind(dc,dw)
d$total <- 53
colnames(d)[1] <- "filename"
d$filename <- substr(d$filename, 6, nchar(d$filename))
# 3 participants' totals are lower due to tech issue during task
d[d$V2 == "JV",]$total <- 52 
d[d$V2 == "MS",]$total <- 52 
d[d$V2 == "ZL",]$total <- 52 
## tapping success and difficult sentences
newd <- merge(d, stinfo, by="filename")

### CALCULATE TAPPING SUCCESS (MEANS/RANGE)
# overall means
totals <- data.frame(xtabs(data=d, ~V2+correct))
# add totals and correct for 3 participants
totals$total <- 53
totals[totals$V2 %in% c("MS","JV", "ZL"),]$total <- 52
totals$perc <- round(totals$Freq/totals$total*100, 0)

tapply(totals$perc, totals$correct, mean) # correct 68.9%
tapply(totals$perc, totals$correct, range) # correct: 60-85
tapply(totals$Freq, totals$correct, sum) # 436 total correct, 197 not.

# check means correct for sentences with ambiguous syllabification
ambsent <- subset(data.frame(xtabs(data=newd, ~V2+correct+ambsyll)), ambsyll == "y")
ambsent$perc <- ambsent$Freq / 12
tapply(ambsent$perc, ambsent$correct, mean) # 60.4% correct
tapply(ambsent$perc, ambsent$correct, range) # correct 33-83%
# check means correct for sentences with hesitation, repetitions
hessent <- subset(data.frame(xtabs(data=newd, ~V2+correct+nonlex)), nonlex == "y")
hessent$perc <- hessent$Freq / 7
tapply(hessent$perc, hessent$correct, mean) # 67% correct
tapply(hessent$perc, hessent$correct, range) # correct 14-100%


############################################################################################
############################################## CORRELATING TAP + RPT DATA
############################################################################################
df$p <- as.factor(df$p)
df$PA <- as.factor(df$PA)
df$DrumVeloc <- as.numeric(as.character(df$DrumVeloc))
df$DrumDur <- as.numeric(as.character(df$DrumDur))*1000

# merge acoustic data with tapping info dataframe 
fulldf <- merge(df, ac, by= c("filename", "syllno", "speaker", "syll")) # still same number of rows as df --> OK
# delete mono_ prefix from filename
fulldf$filename <- substr(fulldf$filename, 6, nchar(fulldf$filename))

# add info about ambiguous syllabification etc.
alldf <- merge(fulldf, stinfo, by=c("speaker","filename")) # still same number of obs 5654, OK
alldf$UI <- as.factor(paste(alldf$filename, alldf$tapper))

### DISTRIBUTION OF TAP VALUES - far from normal
hist(df$DrumVeloc, breaks=25)
hist(df$DrumDur, breaks=35)

# log trans and normalise DrumVeloc and DrumDur by UniqueIdentifier (i.e. each file+tapper combination)
alldf_groupUI <- alldf %>% 
  group_by(UI) %>%
  mutate(
    DrumVeloc.log.sc = scale(log(DrumVeloc)),
    DrumDur.log.sc = scale(log(DrumDur)),
    DrumVeloc.sc = scale(DrumVeloc),
    DrumDur.sc =scale(DrumDur))

###### TAP Distributions log+normed still not normal
shapiro.test(alldf_groupUI$DrumDur.log.sc) # sign diff from normal
shapiro.test(alldf_groupUI$DrumVeloc.log.sc) # sign diff from normal

## ALLDF for acoustics # still UI is correct
summary(alldf_groupUI)
alldf$tapper <- as.factor(alldf$tapper)

#################################################################
# SAVE df as input file for RandomForests
write.csv(alldf_groupUI, file="alldf_groupUI.csv")
#################################################################

############ MODEL + GRAPH correlating tap values with RPT; GLMs by participant

# create an empty list
dfs_participant <- list()
for (i in 1:length(levels(alldf$tapper))){
  partdf <- alldf_groupUI[alldf_groupUI$tapper==levels(alldf$tapper)[i],]
  dfs_participant[[i]] <- partdf
}

# create empty df
modeloutput <- data.frame(tapper=character(),
                          drummeasure=character(),
                          modelestimate=character(),
                          estimate=numeric(),
                          stderror=numeric(),
                          pvalue=numeric(),
                          stringsAsFactors=FALSE)

# run through list with all participants' subdataframes
for (i in 1:length(dfs_participant)){
  parti <- dfs_participant[[i]]
  tapper_anon <- paste0("P",i)
  
  # create DrumDur glm
  durmodel <- glm(formula = p ~ DrumDur.log.sc,
                  family = binomial, 
                  data = parti)
  dur_inter <- coef(summary(durmodel))[1,1]
  dur_inter_err <- coef(summary(durmodel))[1,2]
  # predicted value for "prominent" is intercept + estimate
  dur_est <- coef(summary(durmodel))[2,1] + dur_inter
  dur_est_err <- coef(summary(durmodel))[2,2]
  dur_pvalue <- coef(summary(durmodel))[2,4]
  modeloutput <- rbind(modeloutput, c(tapper_anon,"tap duration","intercept",dur_inter,dur_inter_err,dur_pvalue))
  modeloutput <- rbind(modeloutput, c(tapper_anon,"tap duration","estimate",dur_est,dur_est_err,dur_pvalue))
  
  # create DrumVeloc glm
  velocmodel <- glm(formula = p ~ DrumVeloc.log.sc,
                    family = binomial, 
                    data = parti)
  veloc_inter <- coef(summary(velocmodel))[1,1]
  veloc_inter_err <- coef(summary(velocmodel))[1,2]
  # predicted value for "prominent" is intercept + estimate
  veloc_est <- coef(summary(velocmodel))[2,1] + veloc_inter
  veloc_est_err <- coef(summary(velocmodel))[2,2]
  veloc_pvalue <- coef(summary(velocmodel))[2,4] 
  modeloutput <- rbind(modeloutput, c(tapper_anon,"tap force","intercept",veloc_inter,veloc_inter_err,veloc_pvalue))
  modeloutput <- rbind(modeloutput, c(tapper_anon,"tap force","estimate",veloc_est,veloc_est_err,veloc_pvalue))
}

modeloutput$predval <- modeloutput

# Change column names to work with easier and add column for prominence (intercept = "not prom")
colnames(modeloutput) <- c("tapper", "drummeasure", "modelpar", "estimate", "sterror", "pvalue")
modeloutput$prom <- as.factor(rep(c("n","y")))
summary(modeloutput)
modeloutput$pvalue <- as.numeric(modeloutput$pvalue)
modeloutput$estimate <- as.numeric(modeloutput$estimate)
modeloutput$sterror <- as.numeric(modeloutput$sterror)
modeloutput$drummeasure <- as.factor(modeloutput$drummeasure)
modeloutput$tapper <- as.factor(modeloutput$tapper)
# rename to P1-P12 for graph purposes
modeloutput$tapper <- factor(modeloutput$tapper, levels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "P11", "P12"))
modeloutput$typeprom <- as.factor(paste0(modeloutput$drummeasure,modeloutput$prom))
# add default value for model sign diff
modeloutput$diff <- "n"
# change default to "yes" if p < 0.05
modeloutput[modeloutput$pvalue < 0.05,]$diff <- "y"

# PLOT MODEL PREDICTIONS 
ggplot(data=modeloutput, aes(x=prom, y=estimate)) +
  geom_line(aes(group = drummeasure, colour = diff, size=diff), alpha = 0.5) + 
  geom_point(aes(fill=drummeasure, shape=drummeasure),size = 3, stroke = 1) +
  geom_errorbar(aes(ymin = estimate-sterror, ymax = estimate+sterror), width = 0.2) +
  scale_fill_manual(name = "",
                    breaks = c("tap duration", "tap force"),
                    labels = c("tap duration", "tap force"),
                    values = c("white", "gray")) +
  scale_shape_manual(name = "",
                     breaks = c("tap duration", "tap force"),
                     labels = c("tap duration", "tap force"),
                     values = c(21, 24)) +
  scale_colour_manual(values=c("gray", "chartreuse4"),
                      guide = "none") +
  scale_size_manual(values=c(1,1.5),
                    guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = c(21,24)))) +
  scale_y_continuous("Logistic regression model estimate\n(based on log-scaled measure)") +
  scale_x_discrete(labels=c("not\nprom", "prom")) +
  facet_wrap(~tapper, nrow=2) +
  theme_bw() +
  theme(strip.text.x = element_text(size=15, face="bold"),
        strip.text.y = element_text(size=15, face="bold"),
        strip.background = element_rect(colour="#999999", size=1, fill="white"),
        axis.title.y =element_text(size=15),
        axis.text.y = element_text(size=10),
        axis.text.x= element_text(size=10),
        axis.title.x=element_blank(),
        legend.key.size=unit(2, "line"),
        legend.direction="horizontal",
        legend.key = element_rect(colour = 'white', fill = "white"),
        legend.background = element_rect(colour = 'white', fill = "white"),
        legend.position="top",
        legend.box = "horizontal",
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))