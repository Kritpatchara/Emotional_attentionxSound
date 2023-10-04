# URL ----------------------------------------------------------
# fileURL = address
fileURL = "D:/EmotionalAttention"
# fileURL = "/Users/Pang/Documents/EmotionalAttention/"
#----------------------------------------------------------------
fileURL_txtdata = paste0(fileURL, "/data_envSound_csv")

# Read
files = dir(fileURL_txtdata, pattern = '.csv')
d = data.frame()
for (i in 1:length(files)) {
  temp = read.csv(file.path(fileURL_txtdata, files[i]), header = TRUE)
  temp$runNum = 1:nrow(temp)
  d = rbind(d, temp)
}

# Change column content
d$emo_HA_AN_NE[d$emo_HA_AN_NE == 1] = "HA"
d$emo_HA_AN_NE[d$emo_HA_AN_NE == 2] = "AN"
d$emo_HA_AN_NE[d$emo_HA_AN_NE == 3] = "NE"
d$csi_125_500[d$csi_125_500 == 1] = "125"
d$csi_125_500[d$csi_125_500 == 2] = "500"
d$coninc[d$coninc == 1] = "congruent"
d$coninc[d$coninc == 2] = "incongruent"
d$shape_diamondglass[d$shape_diamondglass == 1] = "diamond"
d$shape_diamondglass[d$shape_diamondglass == 2] = "hourglass"
d$resp[d$resp == 1] = "diamond"
d$resp[d$resp == 2] = "hourglass"

# Change column names
colnames(d) = c("subNum", "runNum",
                "gender", "hand", "age",
                "shape_diamondglass", "cueloc", "coninc",
                "emo_HA_AN_NE", "csi_125_500", "exshuffle",
                "resp", "RT", "hit",
                "isi", "iti", "hit_notslow", "slow",
                "stimonset")

# Actual subjects
actual_subs = c("PN1608")

# Subset
d = subset(d, subNum %in% actual_subs)

# Create a new column for sound condition (PN or NP)
actual_subs = c("PN2402", "PN2403", "PN1504", "PN1505",
                "PN2506", "PN2507", "PN1608", "PN1609",
                "NP1411", "NP2412", "NP2513", "NP1514",
                "NP1515", "NP2516", "NP2517", "NP1618",
                "NP1619", "NP2620")

# Create a new column for sound condition (P or N)
d$soundCond_P_N[d$soundCond == "PN" & d$runNum <= 144] <- "P"
d$soundCond_P_N[d$soundCond == "PN" & d$runNum > 144] <- "N"
d$soundCond_P_N[d$soundCond == "NP" & d$runNum <= 144] <- "N"
d$soundCond_P_N[d$soundCond == "NP" & d$runNum > 144] <- "P"

#---------------------------------------------------------------------
# bar plot for correctRT of different emotion (HA, AN, NE) x congruency (con x incon). 
d.correctRT = subset(d, hit == 1)

# use this when there are multiple subs
# d.agg = aggregate(data = d.correctRT, respRT ~ faceEmo + stimCongruency + subNum, mean)

# Source function ##############################
source("D:/EmotionalAttention/R/summarySE.R")
# source("flanker_emotionalAttention_supFunc.R")
################################################


createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "P")
createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "N")

createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "P")
createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "N")

