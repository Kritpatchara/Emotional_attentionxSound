EmoAtten_ANOVA_hit = function(d, 
                              soundEnv_P_N = "", 
                              csi_125_500 = "",
                              subNumber = "") {
  
  
  if (soundEnv_P_N == "" & csi_125_500 != "") {
    
    # ANOVA for csi_125_500 = csi_125_500
    d %>% subset(csi_125_500 == csi_125_500) %>%
      aggregate(hit ~ subNum + coninc + emo_HA_AN_NE + soundCond_P_N, mean) %>%
      convert_as_factor(subNum, coninc, emo_HA_AN_NE, soundCond_P_N) %>%
      anova_test(dv = hit,
                 wid = subNum,
                 within = c(emo_HA_AN_NE, soundCond_P_N, coninc)) %>%
      get_anova_table() %>%
      print()
    
  } else if (soundEnv_P_N != "" & csi_125_500 == "") {
    
    # ANOVA for soundEnv_P_N = soundCond_P_N
    d %>% subset(soundCond_P_N == soundEnv_P_N) %>%
      aggregate(hit ~ subNum + coninc + emo_HA_AN_NE + csi_125_500, mean) %>%
      convert_as_factor(subNum, coninc, emo_HA_AN_NE, csi_125_500) %>%
      anova_test(dv = hit,
                 wid = subNum,
                 within = c(emo_HA_AN_NE, csi_125_500, coninc)) %>%
      get_anova_table() %>%
      print()
    
  } else if (soundEnv_P_N == "" & csi_125_500 == "") {
    
    # ANOVA for all 
    d %>% aggregate(hit ~ subNum + coninc + emo_HA_AN_NE + soundCond_P_N + csi_125_500, mean) %>%
      convert_as_factor(subNum, coninc, emo_HA_AN_NE, soundCond_P_N, csi_125_500) %>%
      anova_test(dv = hit,
                 wid = subNum,
                 within = c(emo_HA_AN_NE, soundCond_P_N, csi_125_500, coninc)) %>%
      get_anova_table() %>%
      print()
    
  } else if (soundEnv_P_N != "" & csi_125_500 != ""){
    # ANOVA for csi_125_500 = csi_125_500 and soundEnv_P_N = soundCond_P_N
    d %>% subset(csi_125_500 == csi_125_500 & soundCond_P_N == soundEnv_P_N) %>%
      aggregate(hit ~ subNum + coninc + emo_HA_AN_NE, mean) %>%
      convert_as_factor(subNum, coninc, emo_HA_AN_NE) %>%
      anova_test(dv = hit,
                 wid = subNum,
                 within = c(emo_HA_AN_NE, coninc)) %>%
      get_anova_table() %>%
      print()
  }
  
  
}