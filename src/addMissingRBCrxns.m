function model = addMissingRBCrxns(model)
% add missing RBC reactions to WBM models
% What was the reference again for that? Some published RBC model.
% IT Dec 2022

RBC_reactions={
'RBC_2KMBTA'	'2-keto-4-methylthiobutyrate transaminase'	'RBC_2kmb[c] + RBC_glu_L[c] <=> RBC_akg[c] + RBC_met_L[c]'	'Methionine Salvage'	'Tat.1'	'2.6.1.5'	'YES'	'NO'	''''
'RBC_AP4AH1'	'Ap4A hydrolase, asymmetrically'	'RBC_ap4a[c] + RBC_h2o[c] -> RBC_amp[c] + RBC_atp[c] + 2 RBC_h[c]'	'Nucleotides'	'NudT2.2, NudT2.1, NudT2.3'	'EC-3.6.1.17'	'YES'	'NO'	'Hankin et al. 1995 Int J of Biochem and Cell Bio 27(2):201-206'
'RBC_ARD1_hs'	'Acireductone Dixoygenase 1 (homo sapiens)'	'RBC_dhmtp[c] + RBC_o2[c] -> RBC_2kmb[c] + RBC_for[c] + 2 RBC_h[c]'	'Methionine Salvage'	'Adi1'	'1.13.11.53'	'NO'	'NO'	''''
'RBC_ARGt5r'	'L-arginine transport via diffusion'	'RBC_arg_L[e] <=> RBC_arg_L[c]'	'Transport, Extracellular'	''''	''''	''''	'NO'	'Mendes et al. Clinical Science (1997) 93(1):57-64'
'RBC_ARS_hs'	'acireductone synthase (homo sapiens)'	'RBC_dkmpp[c] + RBC_h2o[c] -> RBC_dhmtp[c] + RBC_pi[c]'	'Methionine Salvage'	'Enoph1'	'3.1.3.77'	'NO'	'NO'	''''
'RBC_BILIRBU'	'bilirubin UDP-glucuronosyltransferase'	'RBC_bilirub[c] + 2 RBC_h[c] + RBC_udpglcur[c] -> RBC_bilglcur[c] + RBC_udp[c]'	'Heme Degradation'	'Ugt1A1.1, Ugt1A4.1'	'EC-2.4.1.17'	'NO'	'YES'	'Anderson et al. Am J Med 63(3):359-364 (1977)'
'RBC_C160CPT2rbc'	'C160 reverse'	'RBC_coa[c] + RBC_pmtcrn[c] -> RBC_crn[c] + RBC_pmtcoa[c]'	'Carnitine shuttle'	'Cpt2.1'	''''	'NO'	'YES'	'Arduini et al. JBC 1992 267(18):12673-12681'
'RBC_C181CPT2rbc'	'carnitine octadecenoyl transferase'	'RBC_coa[c] + RBC_odecrn[c] -> RBC_crn[c] + RBC_odecoa[c]'	'Carnitine shuttle'	'Cpt2.1'	''''	'NO'	'YES'	'Arduini et al. JBC 1992 267(18):12673-12681'
'RBC_CAT'	'catalase'	'2 RBC_h2o2[c] -> 2 RBC_h2o[c] + RBC_o2[c]'	'ROS Detoxification'	'Cat.1'	'EC-1.11.1.6'	'YES'	'NO'	'Puget et al. 1977 Biochem Biophys Res Comm 77(4):1525'
'RBC_DM_nadh'	'demand NADH (similar to methemoglobin redox)'	'RBC_nadh[c] -> RBC_h[c] + RBC_nad[c]'	''''	''''	''''	''''	''''	''''
'RBC_EX_ascb_L(e)_[bc]'	'L-Ascorbate exchange'	'RBC_ascb_L[e] <=> ascb_L[bc]'	''''	''''	''''	''''	''''	''''
'RBC_EX_hcyst_L(e)_[bc]'	'L-homocysteine exchange'	'RBC_hcys_L[e] <=> hcys_L[bc]'	''''	''''	''''	''''	''''	''''
'RBC_EX_ncam(e)_[bc]'	'Nicotinamide exchange'	'RBC_ncam[e] <=> ncam[bc]'	''''	''''	''''	''''	''''	''''
'RBC_EX_normete_L(e)_[bc]'	'L-Normetanephrine exchange'	'RBC_normete_L[e] <=> normete_L[bc]'	''''	''''	''''	''''	''''	''''
'RBC_EX_spmd(e)_[bc]'	'Spermidine exchange'	'RBC_spmd[e] <=> spmd[bc]'	''''	''''	''''	''''	''''	''''
'RBC_FACOAL160'	'fatty-acid--CoA ligase (hexadecanoate)'	'RBC_atp[c] + RBC_coa[c] + RBC_hdca[c] <=> RBC_amp[c] + RBC_pmtcoa[c] + RBC_ppi[c]'	'Fatty acid activation'	'Acsl1.1'	'EC-6.2.1.3'	'YES'	'YES'	'Arduini et al. JBC 1992 267(18):12673-12681'
'RBC_FACOAL181'	'fatty-acid--CoA ligase (octadecenoate)'	'RBC_atp[c] + RBC_coa[c] + RBC_ocdcea[c] <=> RBC_amp[c] + RBC_odecoa[c] + RBC_ppi[c]'	'Fatty acid activation'	'Acsl1.1'	'EC-6.2.1.3'	'YES'	'YES'	'Arduini et al. JBC 1992 267(18):12673-12681'
'RBC_FCLT'	'Ferrochelatase'	'RBC_fe2[c] + RBC_ppp9[c] -> 2 RBC_h[c] + RBC_pheme[c]'	'Heme Biosynthesis'	'Fech.1, Fech.2'	'EC-4.99.1.1'	'NO'	'YES'	'Anderson et al. Am J Med 63(3):359-364 (1977)'
'RBC_FE2t1'	'iron (II) transport'	'RBC_fe2[e] <=> RBC_fe2[c]'	'Transport, Extracellular'	''''	''''	''''	'NO'	''''
'RBC_GLYCt'	'glycerol transport via channel'	'RBC_glyc[e] <=> RBC_glyc[c]'	'Transport, Extracellular'	''''	''''	''''	'NO'	'Bruckdorfer et al. 1969 BBA - Biomembranes 183(2):334-345'
'RBC_GLYt7_211_r'	'glycine reversible transport via sodium and chloride symport (2:1:1)'	'RBC_cl[e] + RBC_gly[e] + 2 RBC_na1[e] <=> RBC_cl[c] + RBC_gly[c] + 2 RBC_na1[c]'	'Transport, Extracellular'	'Slc6a9.1, Slc6a9.2'	''''	'NO'	'NO'	'Tunnicliff Comp Biochem Physio 1994 108A(4):471-8'
'RBC_GMPR'	'GMP reductase'	'RBC_gmp[c] + 2 RBC_h[c] + RBC_nadph[c] -> RBC_imp[c] + RBC_nadp[c] + RBC_nh4[c]'	'Nucleotides'	'GmpR.1'	'EC-1.7.1.7'	'YES'	'NO'	'Bontemps et al. Biochem J 1988 250(3):687-96'
'RBC_GULND'	'gulonate dehydrogenase'	'RBC_glcur[c] + RBC_h[c] + RBC_nadph[c] <=> RBC_guln[c] + RBC_nadp[c]'	'Pentose and Glucuronate Interconversions'	''''	'EC-1.1.1.19'	''''	'NO'	'Surgenor The RBC 1974'
'RBC_HCO3E'	'Carbonic anhydrase'	'RBC_co2[c] + RBC_h2o[c] <=> RBC_h[c] + RBC_hco3[c]'	'Miscellaneous'	'Ca1.1, Ca2.1, Ca3.1'	'EC-4.2.1.1'	'YES'	'NO'	'Surgenor The RBC 1974'
'RBC_Ht'	'proton diffusion'	'RBC_h[e] <=> RBC_h[c]'	'Transport, Extracellular'	''''	''''	''''	''''	''''
'RBC_L_LACt2r'	'L_Lactate reversible transport via proton symport'	'RBC_h[e] + RBC_lac_L[e] <=> RBC_h[c] + RBC_lac_L[c]'	'Transport, Extracellular'	'Slc16a1.1'	''''	'YES'	'NO'	'McLellan et al. 1988 Mechanisms of Ageing and Dev 48(1):63-71, Thornalley Biochem J 1988;254:751-755.'
'RBC_LNLCCPT2rbc'	'carnitine transferase'	'RBC_coa[c] + RBC_lnlccrn[c] -> RBC_crn[c] + RBC_lnlccoa[c]'	'Fatty acid activation'	'Acsl1.1'	'EC-6.2.1.3'	'YES'	'YES'	'Arduini et al. JBC 1992 267(18):12673-12681'
'RBC_MALt'	'Malate exchange, diffusion'	'RBC_mal_L[e] <=> RBC_mal_L[c]'	'Transport, Extracellular'	''''	'EC-3.2.1.20'	''''	'NO'	'Simpson et al. Biochimica et Biophysica Acta 1982 707:191-200'
'RBC_NACt'	'Nicotinic acid uptake'	'RBC_nac[e] <=> RBC_nac[c]'	'Transport, Extracellular'	''''	''''	''''	'NO'	'Sofue et al. Biochem J (1992) 288:669-674, Totskii et al. Ukr Biokhim Zh (1979) 51(3):241-5'
'RBC_NCAMt'	'Nicotinamide acid uptake'	'RBC_ncam[e] <=> RBC_ncam[c]'	'Transport, Extracellular'	''''	''''	''''	'NO'	'Sofue et al. Biochem J (1992) 288:669-674, Totskii et al. Ukr Biokhim Zh (1979) 51(3):241-5'
'RBC_OROATP'	'Orotate transporter'	'RBC_h[e] + RBC_orot[e] -> RBC_h[c] + RBC_orot[c]'	'Transport, Extracellular'	''''	''''	''''	'NO'	'Webster et al. 2001 The Metabolic and Molecular Bases of Inherited Diseases, "Hereditary orotic aciduria and other disorders of pyrimidine metabolism"'
'RBC_PIt'	'Inorganic phosphate exchange, diffusion'	'RBC_pi[c] <=> RBC_pi[e]'	'Transport, Extracellular'	''''	''''	''''	''''	''''
'RBC_PMANM'	'phosphomannomutase'	'RBC_man1p[c] <=> RBC_man6p[c]'	'Fructose and Mannose Metabolism'	'Pgm2, Pmm2.1'	'EC-5.4.2.8'	'YES'	'NO'	'Beutler et al. 1969 J Clin Invest 48(3):461-6'
'RBC_PPPGO'	'protoporphyrinogen oxidase (aerobic)'	'RBC_pppg9[c] + 1.5 RBC_o2[c] -> RBC_ppp9[c] + 3 RBC_h2o[c]'	'Heme Biosynthesis'	'Ppox.1'	'EC-1.3.3.4'	'NO'	'NO'	'Anderson et al. Am J Med 63(3):359-364 (1977)'
'RBC_PTRCt'	'puterscine transport (faciliated)'	'RBC_ptrc[e] <=> RBC_ptrc[c]'	'Transport, Extracellular'	''''	''''	'NO'	'NO'	'Fukumoto et al. Biochimica et Biophysica Acta 1996 1282(1):48-56'
'RBC_PYAM5PO'	'pyridoxamine 5-phosphate oxidase'	'RBC_h2o[c] + RBC_o2[c] + RBC_pyam5p[c] -> RBC_h2o2[c] + RBC_nh4[c] + RBC_pydx5p[c]'	'Vitamin B6 Metabolism'	'PnpO.1'	'EC-1.4.3.5'	'YES'	'YES'	'Giraud et al. Am J Clin Nutr 1995 62:104-9'
'RBC_RIBFLVt3o'	'riboflavin transport out (ATP dependent)'	'RBC_atp[c] + RBC_h2o[c] + RBC_ribflv[c] -> RBC_adp[c] + RBC_h[c] + RBC_pi[c] + RBC_ribflv[e]'	'Transport, Extracellular'	'Abcg2'	''''	'NO'	'YES'	'van Herwaarden et al. 2007 Mol + Cell Biol 27(4):1247-1253'
'RBC_SPMDt'	'spermidine transport (faciliated)'	'RBC_spmd[e] <=> RBC_spmd[c]'	'Transport, Extracellular'	''''	''''	'NO'	'NO'	'Fukumoto et al. Biochimica et Biophysica Acta 1996 1282(1):48-56'    

};
gprs(1:size(RBC_reactions,1)) = {''};
model = addReactionsHH(model,RBC_reactions(:,1),RBC_reactions(:,2),RBC_reactions(:,3),gprs,RBC_reactions(:,4));