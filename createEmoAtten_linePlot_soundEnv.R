createEmoAtten_linePlot_soundEnv = function(d, 
                                   plotType = "hitRate", 
                                   soundEnv_P_N = "", 
                                   CSI = "",
                                   subNumber = "") {
  # library
  library(Rmisc)
  library(ggplot2)
  library(rstatix)
  library(dplyr)
  library(magrittr)
  
  # subset by plotType
  if (plotType == "correctRT") {
    d.plotType = subset(d, hit == 1); 
    depVar = "RT"
    ymin = 0.00
    ymax = 1.00
  } else if (plotType == "hitRate") {
    d.plotType = d; 
    depVar = "hit"
    ymax = 1.00
    ymin = 0.00
  }
  
  # subset by subject
  if (subNumber != "") {
    d.plotType = subset(d.plotType, subNum == subNumber);
  }
  
  print("-------------------------------------------------------------")
  
  if (soundEnv_P_N == "" & CSI != "") {
    # subset by csi only
    d.bycsi = subset(d.plotType, csi_125_250 == CSI)
  } else if (soundEnv_P_N != "" & CSI == "") {
    # subset by soundEnv_P_N only
    d.bycsi = subset(d.plotType, soundCond_P_N == soundEnv_P_N)
  } else if (soundEnv_P_N == "" & CSI == "") {
    # All data
    d.bycsi = d.plotType
  } else {
    # subset by csi and faceOrient
    d.bycsi = subset(d.plotType, csi == CSI & soundCond_P_N == soundEnv_P_N)
  }
  
  # create summary of data with error-bar information ##########################
  if (plotType == "correctRT") {
    d.agg = d.bycsi %>% 
      aggregate(RT ~ emo_HA_AN_NE + coninc + subNum, mean) %>%
      convert_as_factor(emo_HA_AN_NE, coninc)
  } else if (plotType == "hitRate") {
    d.agg = d.bycsi %>% 
      aggregate(hit ~ emo_HA_AN_NE + coninc + subNum, mean) %>%
      convert_as_factor(emo_HA_AN_NE, coninc)
  }
  
  
  # Change column names
  colnames(d.agg) = c("faceEmo", "stimCongruency", "subNum", "depVar")
  
  # d.sum = summarySE(data = d.agg, 
  #                         measurevar = "depVar", 
  # groupvars = c("faceEmo","stimCongruency"))
  d.sum = summarySEwithin(data = d.agg,
                          measurevar = "depVar",
                          withinvars = c("faceEmo","stimCongruency"),
                          idvar = "subNum")
  colnames(d.sum) <- c("faceEmo","stimCongruency","N","depVar","sd","se","ci")
  print(d.sum)
  
  
  # Create Plot ################################################################
  p1 = ggplot(data = d.sum, mapping = aes(faceEmo, depVar, group = stimCongruency)) +
    geom_line(stat= "identity", aes(linetype = stimCongruency)) +
    geom_point() +
    geom_errorbar(aes(ymin = depVar-se, ymax = depVar+se),
                  width = 0.2, position= position_dodge(0.05)) +
    theme_classic() +
    scale_y_continuous(limits = c(ymin, ymax)) +
    ylab(plotType) + xlab("Face emotion") +
    ggtitle(paste("sound=", soundEnv_P_N, " , csi=", CSI))
  
  print(p1)
  
  # 
  # # Save PDF
  # if (saveFile) {
  #   fileURL = "~/Documents/EmotionalAttention/R/"
  #   if (subNumber =="") {
  #     path = paste0(fileURL,"Rplot/",plotType,"_",faceOri,"_",CSI,".pdf");
  #   } else {
  #     dir.create(paste0(fileURL,"Rplot/", "sub_", subNumber))
  #     path = paste0(fileURL,"Rplot/", "sub_", subNumber, "/", plotType,"_",faceOri,"_",CSI,".pdf")
  #   }
  #   
  #   # Save plot as pdf
  #   pdf(path, width = 5, height = 3)
  #   print(p1)
  #   dev.off()
  # } else {
  #   print(p1)
  # }
  # 
  # 
  # # Perform anova ##############################################################
  # if (plotType == "hitRate") {
  #   EmoAtten_ANOVA_hit(d, faceOri = faceOri, CSI = CSI, subNumber = subNumber);
  # } else if (plotType == "correctRT") {
  #   EmoAtten_ANOVA_correctRT(d, faceOri = faceOri, CSI = CSI, subNumber = subNumber)
  # }
  
}

