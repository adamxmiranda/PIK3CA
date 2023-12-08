library(tidyverse)

wt_ex9_res <- read_tsv("qPCR12-4_WTex9.txt")
ex20_res <- read_tsv("qPCR12-4_ex20.txt")

wt_res <- wt_ex9_res %>% head(n = 48)
wt_res$Genotype <- c(rep("WT", 48))
wt_res_bytech <- wt_res %>% head(12)
ACTBrepCT <- wt_res %>% slice(13:24) %>% select(c(6)) %>% as.vector()
AREGrep1CT <- wt_res %>% slice(25:36) %>% select(c(6)) %>% as.vector()
AREGrep2CT <- wt_res %>% slice(37:48) %>% select(c(6)) %>% as.vector()
wt_res_bytech <- cbind(wt_res_bytech, ACTBrepCT, AREGrep1CT, AREGrep2CT)
names(wt_res_bytech) <- c("sample_ID", "rem1", "rem2", "rem3", "rem4", 
                          "ACTBrep1", "Genotype", "ACTBrep2", "AREGrep1", "AREGrep2")
wt_res_bytech$Treatment <- c(rep("Negative", 4), rep("AREG", 4), rep("none", 4))
wt_res_bytech <- wt_res_bytech %>% select(c(1,7, 11, 6, 8:10))
wt_res_bytech

ex9_res <- wt_ex9_res %>% slice(49:96)
ex9_res$Genotype <- c(rep("Ex9", 48))
ex9_res_bytech <- ex9_res %>% head(12)
ACTBrepCT <- ex9_res %>% slice(13:24) %>% select(c(6)) %>% as.vector()
AREGrep1CT <- ex9_res %>% slice(25:36) %>% select(c(6)) %>% as.vector()
AREGrep2CT <- ex9_res %>% slice(37:48) %>% select(c(6)) %>% as.vector()
ex9_res_bytech <- cbind(ex9_res_bytech, ACTBrepCT, AREGrep1CT, AREGrep2CT)
names(ex9_res_bytech) <- c("sample_ID", "rem1", "rem2", "rem3", "rem4", 
                           "ACTBrep1", "Genotype", "ACTBrep2", "AREGrep1", "AREGrep2")
ex9_res_bytech$Treatment <- c(rep("Negative", 4), rep("AREG", 4), rep("none", 4))
ex9_res_bytech <- ex9_res_bytech %>% select(c(1,7, 11, 6, 8:10))
ex9_res_bytech

ex20_res$Genotype <- c(rep("Ex20", 48))
ex20_res_bytech <- ex20_res %>% head(12)
ACTBrepCT <- ex20_res %>% slice(13:24) %>% select(c(6)) %>% as.vector()
AREGrep1CT <- ex20_res %>% slice(25:36) %>% select(c(6)) %>% as.vector()
AREGrep2CT <- ex20_res %>% slice(37:48) %>% select(c(6)) %>% as.vector()
ex20_res_bytech <- cbind(ex20_res_bytech, ACTBrepCT, AREGrep1CT, AREGrep2CT)
names(ex20_res_bytech) <- c("sample_ID", "rem1", "rem2", "rem3", "rem4", 
                            "ACTBrep1", "Genotype", "ACTBrep2", "AREGrep1", "AREGrep2")
ex20_res_bytech$Treatment <- c(rep("Negative", 4), rep("AREG", 4), rep("none", 4))
ex20_res_bytech <- ex20_res_bytech %>% select(c(1,7, 11, 6, 8:10))
ex20_res_bytech[ex20_res_bytech=="Undetermined"] <- 40
ex20_res_bytech$ACTBrep1 <- as.numeric(ex20_res_bytech$ACTBrep1)
ex20_res_bytech$ACTBrep2 <- as.numeric(ex20_res_bytech$ACTBrep2)
ex20_res_bytech$AREGrep1 <- as.numeric(ex20_res_bytech$AREGrep1)
ex20_res_bytech$AREGrep2 <- as.numeric(ex20_res_bytech$AREGrep2)
ex20_res_bytech


wt_res_bytech <- wt_res_bytech %>% mutate(ACTB_av = (ACTBrep1 + ACTBrep2)/2)
ex9_res_bytech <- ex9_res_bytech %>% mutate(ACTB_av = (ACTBrep1 + ACTBrep2)/2)
ex20_res_bytech <- ex20_res_bytech %>% mutate(ACTB_av = (ACTBrep1 + ACTBrep2)/2)
wt_res_bytech <- wt_res_bytech %>% mutate(AREG_av = (AREGrep1 + AREGrep2)/2)
ex9_res_bytech <- ex9_res_bytech %>% mutate(AREG_av = (AREGrep1 + AREGrep2)/2)
ex20_res_bytech <- ex20_res_bytech %>% mutate(AREG_av = (AREGrep1 + AREGrep2)/2)

wt_res_bytech <- wt_res_bytech %>% mutate(DeltCT = AREG_av - ACTB_av)
ex9_res_bytech <- ex9_res_bytech %>% mutate(DeltCT = AREG_av - ACTB_av)
ex20_res_bytech <- ex20_res_bytech %>% mutate(DeltCT = AREG_av - ACTB_av)

wt_none_av <- wt_res_bytech %>% filter(Treatment == "none") %>% select(c(10))
wt_none_av <- mean(wt_none_av$DeltCT)
ex9_none_av <- ex9_res_bytech %>% filter(Treatment == "none") %>% select(c(10))
ex9_none_av <- mean(ex9_none_av$DeltCT)
ex20_none_av <- ex20_res_bytech %>% filter(Treatment == "none") %>% select(c(10))
ex20_none_av <- mean(ex20_none_av$DeltCT)

wt_res_bytech <- wt_res_bytech %>% mutate(DeltDeltCT = DeltCT - wt_none_av)
ex9_res_bytech <- ex9_res_bytech %>% mutate(DeltDeltCT = DeltCT - ex9_none_av)
ex20_res_bytech <- ex20_res_bytech %>% mutate(DeltDeltCT = DeltCT - ex20_none_av)
wt_res_bytech <- wt_res_bytech %>% mutate(foldEXP = 2^(-DeltDeltCT))
ex9_res_bytech <- ex9_res_bytech %>% mutate(foldEXP = 2^(-DeltDeltCT))
ex20_res_bytech <- ex20_res_bytech %>% mutate(foldEXP = 2^(-DeltDeltCT))
total_df <- rbind(wt_res_bytech, ex9_res_bytech, ex20_res_bytech)
total_df

total_df$relexp <- (40- total_df$AREG_av)/(40 - total_df$ACTB_av)
total_df

WTnull <- total_df %>% filter(Genotype == "WT" & Treatment == "none")
E545Knull <- total_df %>% filter(Genotype == "Ex9" & Treatment == "none")
H1047Rnull <- total_df %>% filter(Genotype == "Ex20" & Treatment == "none")
WTrelnorm <- mean(WTnull$relexp)
E545Krelnorm <- mean(E545Knull$relexp)
H1047Rrelnorm <- mean(H1047Rnull$relexp)

WTall <- total_df %>% filter(Genotype == "WT")
E545Kall <- total_df %>% filter(Genotype == "Ex9")
H1047Rall <- total_df %>% filter(Genotype == "Ex20")

WTall$relnorm <- WTall$relexp / WTrelnorm
E545Kall$relnorm <- E545Kall$relexp / E545Krelnorm
H1047Rall$relnorm <- H1047Rall$relexp / H1047Rrelnorm

normedtotal_df <- rbind(WTall, E545Kall, H1047Rall)




