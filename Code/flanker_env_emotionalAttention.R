# URL ----------------------------------------------------------
# fileURL = address
fileURL = "~/Documents/EmotionalAttention"
# fileURL = "/Users/Pang/Documents/EmotionalAttention/"
#----------------------------------------------------------------
fileURL_txtdata = paste0(fileURL, "/data_envSound_csv")

# Read
files = dir(fileURL_txtdata, pattern = '.csv')
d = data.frame()
for (i in 1:length(files)) {
  temp = read.csv(file.path(fileURL_txtdata, files[i]), header = TRUE)
  temp$runNum = 1:nrow(temp)

    
  # Get aud type
  temp2 = strsplit(files[i],"_")[[1]][2]
  temp3 = strsplit(temp2, ".csv")[[1]][1]
  temp$audType = temp3

  # Get file type
  if(grepl("mock", files[i], ignore.case = T)){
    temp$fileType = "mock"
  } else {
    temp$fileType = "actual"
  }
  
  temp$fileName = files[i]
  # Bind data  
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
d$soundCond_P_N[d$audType == "audType1"] <- "P"
d$soundCond_P_N[d$audType == "audType2"] <- "N"

# Change column names
colnames(d) = c("subNum", "runNum",
                "gender", "hand", "age",
                "shape_diamondglass", "cueloc", "coninc",
                "emo_HA_AN_NE", "csi_125_500", "exshuffle",
                "resp", "RT", "hit",
                "isi", "iti", "hit_notslow", "slow",
                "stimonset", "audType", "fileType", "fileName",
                "soundCond_P_N")

# Subset
d = subset(d, fileType == "actual")

# Create a new column for sound condition (PN or NP)
d$soundCond = substr(d$subNum, 1,2)

# Sanity check
table(d$coninc, d$emo_HA_AN_NE)
table(d$coninc, d$csi_125_500)
table(d$coninc, d$cueloc)
table(d$coninc, d$soundCond_P_N)

#---------------------------------------------------------------------
# bar plot for correctRT of different emotion (HA, AN, NE) x congruency (con x incon). 
d.correctRT = subset(d, hit == 1)

# use this when there are multiple subs
# d.agg = aggregate(data = d.correctRT, respRT ~ faceEmo + stimCongruency + subNum, mean)

# Source function ##############################
source("D:/EmotionalAttention/R/summarySE.R")
# source("flanker_emotionalAttention_supFunc.R")
################################################

createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "P", CSI = 125)
createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "N", CSI = 125)
createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "P", CSI = 500)

createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "N", CSI = 500)

# createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "P", CSI = 125)
# createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "N", CSI = 125)
# createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "P", CSI = 500)
# createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "N", CSI = 500)

createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "P")
createEmoAtten_linePlot_soundEnv(d = d,plotType = "hitRate",soundEnv_P_N = "N")

createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "P")
createEmoAtten_linePlot_soundEnv(d = d,plotType = "correctRT",soundEnv_P_N = "N")

