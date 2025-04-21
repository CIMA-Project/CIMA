library(TeachingDemos)

#画IKZF4
source("/CIMA/Script/plot_OmicsSMR_xQTL_with_legend.r")

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Treg_FOXP3_As_OPERA.IKZF4.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Treg_FOXP3_As_OPERA.IKZF4.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.00028,msmr_thresh=0.001,esmr_heidi = 0.1,msmr_heidi = 0.1,trait_name="As",window= 500,highlight = 'chr12_56050848')
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Treg_FOXP3_level_of_interleukin_12_subunit_beta_in_blood_OPERA.IKZF4.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Treg_FOXP3_level_of_interleukin_12_subunit_beta_in_blood_OPERA.IKZF4.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.000067,msmr_thresh=0.001,esmr_heidi = 0.05,msmr_heidi = 0.05,trait_name="level_of_interleukin_12_subunit_beta_in_blood",window= 500,highlight = 'chr12_56050848')
dev.off()

source("/CIMA/Script/plot_SMR.r")

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Treg_FOXP3_As_SMR.IKZF4.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Treg_FOXP3_As_gene_GWAS_SMR_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name="As") 
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Treg_FOXP3_level_of_interleukin_12_subunit_beta_in_blood_SMR.IKZF4.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Treg_FOXP3_level_of_interleukin_12_subunit_beta_in_blood_gene_GWAS_SMR_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name="level_of_interleukin_12_subunit_beta_in_blood") 
dev.off()

#画CCR6
source("/CIMA/Script/plot_OmicsSMR_xQTL_with_legend.r")
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_Bn_TCL1A_RA_OPERA.CCR6.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/Bn_TCL1A_RA_OPERA.CCR6.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=9.9e-10,msmr_thresh=6.5e-08,esmr_heidi = 0.2,msmr_heidi = 0.3,trait_name="RA",window= 500,highlight = 'chr6_167125409')
dev.off()

pdf("/CIMA/Result/plot/SMR_OPERA_plot/Bn_TCL1A_RA_OPERA.CCR6_effect.pdf", width = 8, height = 8)
omicSMREffectPlot(data =SMRData,exposure_probe = "chr6:167122097-167122598")
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cDC2_CD1C_RA_OPERA.CCR6.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cDC2_CD1C_RA_OPERA.CCR6.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=1.3e-12,msmr_thresh=1.3e-12,esmr_heidi = 0.04,msmr_heidi = 0.04,trait_name="RA",window= 500,highlight = 'chr6_167125409')
dev.off()


SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_Transitional_B_SOX4_RA_OPERA.CCR6.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/Transitional_B_SOX4_RA_OPERA.CCR6.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=1.7e-6,msmr_thresh=1.7e-6,esmr_heidi = 0.2,msmr_heidi = 0.2,trait_name="RA",window= 500,highlight = 'chr6_167125409')
dev.off()



source("/CIMA/Script/plot_SMR.r")
SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_Bn_TCL1A_RA_SMR.CCR6.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/Bn_TCL1A_RA_SMR.CCR6_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "RA")
dev.off()


SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_Transitional_B_SOX4_RA_SMR.CCR6.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/Transitional_B_SOX4_RA_SMR.CCR6_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "RA") 
dev.off()


SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cDC2_CD1C_RA_SMR.CCR6.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cDC2_CD1C_RA_SMR.CCR6_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "RA")
dev.off()



#画PADI2
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_CD14_RA_OPERA.PADI2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_RA_OPERA.PADI2.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.000002,msmr_thresh=0.000002,esmr_heidi = 0.05,msmr_heidi = 0.05,trait_name="RA",window= 210,highlight = 'chr1_17089967')
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_CD14_RA_SMR.PADI2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_RA_SMR.PADI2_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "RA")
dev.off()

#画TLR1
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tfh-like_CXCR5_macrophage_inflammatory_protein_1a_measurement_OPERA.TLR1.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tfh-like_CXCR5_macrophage_inflammatory_protein_1a_measurement_OPERA.TLR1.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.00003,msmr_thresh=0.00003,esmr_heidi = 0.05,msmr_heidi = 0.05,trait_name="RA",window= 500,highlight = 'chr4_38790903')
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tfh-like_CXCR5_macrophage_inflammatory_protein_1a_measurement_SMR.TLR1.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_RA_SMR.PADI2_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "RA")
dev.off()
SMRLocusPlot(data=SMRData,smr_thresh=0.00003,heidi_thresh = 0.05,plotWindow = 500)

#画S100A12
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_CD14_protein_S100_A12_measurement_OPERA.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_protein_S100_A12_measurement_OPERA.S100A12.pdf", width = 6, height = 10)
omicSMRLocusPlot(data=SMRData,esmr_thresh=1.2e-06,msmr_thresh=5.4e-9,esmr_heidi = 0.1,msmr_heidi = 0.1,trait_name="protein_S100_A12_measurement",window= 500,highlight = 'chr1_153366838')
dev.off()

pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_protein_S100_A12_measurement_OPERA.S100A12_effect.pdf", width = 8, height = 8)
omicSMREffectPlot(data =SMRData,exposure_probe = "chr1:153367064-153367565")
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_IL1B_protein_S100_A12_measurement_OPERA.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_IL1B_protein_S100_A12_measurement_OPERA.S100A12.pdf", width = 6, height = 10)
omicSMRLocusPlot(data=SMRData,esmr_thresh=4.2e-06,msmr_thresh=4.2e-06,esmr_heidi = 0.1,msmr_heidi = 0.1,trait_name="protein_S100_A12_measurement",window= 500,highlight = 'chr1_153366838')
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_IFI44L_protein_S100_A12_measurement_OPERA.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_IFI44L_protein_S100_A12_measurement_OPERA.S100A12.pdf", width = 6, height = 10)
omicSMRLocusPlot(data=SMRData,esmr_thresh=2.3e-07	,msmr_thresh=2.3e-07,esmr_heidi = 0.1,msmr_heidi = 0.1,trait_name="protein_S100_A12_measurement",window= 500,highlight = 'chr1_153366838')
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tn_CCR7_protein_S100_A12_measurement_OPERA.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tn_CCR7_protein_S100_A12_measurement_OPERA.S100A12.pdf", width = 6, height = 8)
omicSMRLocusPlot(data=SMRData,esmr_thresh=1.3e-07,msmr_thresh=1.3e-07,esmr_heidi = 0.1,msmr_heidi = 0.1,trait_name="protein_S100_A12_measurement",window= 500,highlight = 'chr1_153366838')
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_Bn_TCL1A_protein_S100_A12_measurement_OPERA.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/Bn_TCL1A_protein_S100_A12_measurement_OPERA.S100A12.pdf", width = 6, height = 8)
omicSMRLocusPlot(data=SMRData,esmr_thresh=4.6e-06,msmr_thresh=4.6e-06,esmr_heidi = 0.05,msmr_heidi = 0.05,trait_name="protein_S100_A12_measurement",window= 500,highlight = 'chr1_153366838')
dev.off()


source("/CIMA/Script/plot_SMR.r")
SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_CD14_protein_S100_A12_measurement_SMR.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/myplot_cMono_CD14_protein_S100_A12_measurement_SMR.S100A12_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "protein_S100_A12_measurement")
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_IL1B_protein_S100_A12_measurement_SMR.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/myplot_cMono_IL1B_protein_S100_A12_measurement_SMR.S100A12_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "protein_S100_A12_measurement")
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tn_CCR7_protein_S100_A12_measurement_SMR.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/myplot_CD4_Tn_CCR7_protein_S100_A12_measurement_SMR.S100A12_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "protein_S100_A12_measurement")
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_Bn_TCL1A_protein_S100_A12_measurement_SMR.S100A12.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/myplot_Bn_TCL1A_protein_S100_A12_measurement_SMR.S100A12_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "protein_S100_A12_measurement")
dev.off()

#画TMEM258和FADS2

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tfh-like_CXCR5_LPC_20_4_OPERA.TMEM258.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tfh-like_CXCR5_LPC_20_4_OPERA.TMEM258.pdf", width = 6, height = 10)
omicSMRLocusPlot(data=SMRData,esmr_thresh=1.463644e-06,msmr_thresh=1.463644e-06,esmr_heidi = 0.35,msmr_heidi = 0.35,trait_name="LPC_20:4",window= 500,highlight = 'chr11_61781087')
dev.off()

SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD8_CTL_GZMB_stem_Cell_Factor_measurement_OPERA.FADS2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD8_CTL_GZMB_stem_Cell_Factor_measurement_OPERA.FADS2.pdf", width = 6, height = 10)
omicSMRLocusPlot(data=SMRData,esmr_thresh=3.3e-07,msmr_thresh=3.3e-07,esmr_heidi = 0.35,msmr_heidi = 0.35,trait_name="stem_Cell_Factor_measurement",window= 500,highlight = 'chr11_61811991')
dev.off()

source("/CIMA/Script/plot_SMR.r")
SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tfh-like_CXCR5_LPC_20_4_smr.TMEM258.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tfh-like_CXCR5_LPC_20_4_smr.TMEM258_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "LPC20:4")
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD8_CTL_GZMB_stem_Cell_Factor_measurement_smr.FADS2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD8_CTL_GZMB_stem_Cell_Factor_measurement_smr.FADS2_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "stem_Cell_Factor_measurement")
dev.off()

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD8_CTL_GZMB_TNF_related_apoptosis_inducing_ligand_measurement_smr.FADS2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD8_CTL_GZMB_TNF_related_apoptosis_inducing_ligand_measurement_smr.FADS2_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name = "TNF_related_apoptosis_inducing_ligand_measurement")
dev.off()

#画FADS2
detected_list =readRDS('/CIMA/Result/summary/20250212_gene_detected_list_by_CT.rds')
selected_names <- names(detected_list)[sapply(detected_list, function(x) any(grepl("PADI2", x)))]
selected_names

#画SLC16A11
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tn_CCR7_T2D_OPERA.SLC16A11.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tn_CCR7_T2D_OPERA.SLC16A11.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.0000000701,msmr_thresh=0.00000127,esmr_heidi = 0.2,msmr_heidi = 0.2,trait_name="T2D",window= 500,highlight = 'chr17_7037074')
dev.off()

pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tn_CCR7_T2D_OPERA.SLC16A11_effect.pdf", width = 8, height = 8)
omicSMREffectPlot(data =SMRData,exposure_probe = "chr17:7042106-7042607")
dev.off()

source("/CIMA/Script/plot_SMR.r")

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_CD4_Tn_CCR7_T2D_SMR.SLC16A11.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/CD4_Tn_CCR7_T2D_SMR.SLC16A11_gene_GWAS_SMR_effect.pdf", width = 8, height = 8)
SMREffectPlot(data=SMRData, trait_name="T2D") 
dev.off()

#画HEMGN
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_MK_GP9_C_X_C_motif_chemokine_5_measurement_OPERA.HEMGN.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/MK_GP9_C_X_C_motif_chemokine_5_measurement_OPERA.HEMGN.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.00003,msmr_thresh=0.00003,esmr_heidi = 0.001,msmr_heidi = 0.001,trait_name="test",window= 500,highlight = 'chr9_97922475')
dev.off()

source("/CIMA/Script/plot_SMR.r")

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_MK_GP9_C_X_C_motif_chemokine_5_measurement_SMR.HEMGN.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/MK_GP9_C_X_C_motif_chemokine_5_measurement_SMR.HEMGN_gene_GWAS_SMR_effect.pdf", width = 8, height = 8)
SMREffectPlot(SMRData)
dev.off()

#画TMEM199
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_CD14_TNF_related_activation_induced_cytokine_measurement_OPERA.TMEM199.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_TNF_related_activation_induced_cytokine_measurement_OPERA.TMEM199.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.00003,msmr_thresh=0.01,esmr_heidi = 0.02,msmr_heidi = 0.001,trait_name="test",window= 500,highlight = 'chr17_28395709')
dev.off()

source("/CIMA/Script/plot_SMR.r")

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_cMono_CD14_TNF_related_activation_induced_cytokine_measurement_SMR.TMEM199.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/cMono_CD14_TNF_related_activation_induced_cytokine_measurement_SMR.TMEM199_gene_GWAS_SMR_effect.pdf", width = 8, height = 8)
SMREffectPlot(SMRData)
dev.off()

#画FADS2
SMRData = ReadomicSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_NKT_NCR1_TNF_related_apoptosis_inducing_ligand_measurement_OPERA.FADS2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/NKT_NCR1_TNF_related_apoptosis_inducing_ligand_measurement_OPERA.FADS2.pdf", width = 6, height = 12)
omicSMRLocusPlot(data=SMRData,esmr_thresh=0.000000404213,msmr_thresh=0.000000404213,esmr_heidi = 0,msmr_heidi = 0,trait_name="test",window= 500,highlight = 'chr11_61820833')
dev.off()

source("/CIMA/Script/plot_SMR.r")

SMRData = ReadSMRData("/CIMA/Result/plot/SMR_OPERA_plot/plot/myplot_NKT_NCR1_TNF_related_apoptosis_inducing_ligand_measurement_SMR.FADS2.txt")
pdf("/CIMA/Result/plot/SMR_OPERA_plot/NKT_NCR1_TNF_related_apoptosis_inducing_ligand_measurement_SMR.FADS2_gene_GWAS_SMR_effect.pdf", width = 8, height = 8)
SMREffectPlot(SMRData)
dev.off()
