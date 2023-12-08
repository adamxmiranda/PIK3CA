firstmolecule <- read_tsv(file = "firstmolecule_exp.txt")
secondmolecule <- read_tsv(file = "secondmolecule_exp.txt")
thirdmolecule <- total_df
firstmolecule[firstmolecule=="AREG"] <- "AREG_1"
secondmolecule[secondmolecule=="AREG"] <- "AREG_2"
thirdmolecule[thirdmolecule=="AREG"] <- "AREG_3"


combined_df <- rbind(firstmolecule, secondmolecule, thirdmolecule)

combined_df_negative <- combined_df %>% filter(Treatment == "Negative")
combined_df_negativeWT <- combined_df_negative %>% filter(Genotype == "WT")
combined_df_negativeWT$foldEXP
combined_df_negativeE545K <- combined_df_negative %>% filter(Genotype == "Ex9")
combined_df_negativeE545K$foldEXP
combined_df_negativeH1047R <- combined_df_negative %>% filter(Genotype == "Ex20")
combined_df_negativeH1047R$foldEXP


combined_df_null <- combined_df %>% filter(Treatment == "none")
combined_df_nullWT <- combined_df_null %>% filter(Genotype == "WT")
combined_df_nullWT$foldEXP
combined_df_nullE545K <- combined_df_null %>% filter(Genotype == "Ex9")
combined_df_nullE545K$foldEXP
combined_df_nullH1047R <- combined_df_null %>% filter(Genotype == "Ex20")
combined_df_nullH1047R$foldEXP

combined_df_AREG1 <- combined_df %>% filter(Treatment == "AREG_1")
combined_df_AREG1WT <- combined_df_AREG1 %>% filter(Genotype == "WT")
combined_df_AREG1WT$foldEXP
combined_df_AREG1E545K <- combined_df_AREG1 %>% filter(Genotype == "Ex9")
combined_df_AREG1E545K$foldEXP
combined_df_AREG1H1047R <- combined_df_AREG1 %>% filter(Genotype == "Ex20")
combined_df_AREG1H1047R$foldEXP

combined_df_AREG2 <- combined_df %>% filter(Treatment == "AREG_2")
combined_df_AREG2WT <- combined_df_AREG2 %>% filter(Genotype == "WT")
combined_df_AREG2WT$foldEXP
combined_df_AREG2E545K <- combined_df_AREG2 %>% filter(Genotype == "Ex9")
combined_df_AREG2E545K$foldEXP
combined_df_AREG2H1047R <- combined_df_AREG2 %>% filter(Genotype == "Ex20")
combined_df_AREG2H1047R$foldEXP

combined_df_AREG3 <- combined_df %>% filter(Treatment == "AREG_3")
combined_df_AREG3WT <- combined_df_AREG3 %>% filter(Genotype == "WT")
combined_df_AREG3WT$foldEXP
combined_df_AREG3E545K <- combined_df_AREG3 %>% filter(Genotype == "Ex9")
combined_df_AREG3E545K$foldEXP
combined_df_AREG3H1047R <- combined_df_AREG3 %>% filter(Genotype == "Ex20")
combined_df_AREG3H1047R$foldEXP
