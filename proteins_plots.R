# ECM = CC	GO:0031012	extracellular matrix
# cardiovascular = BP	GO:0007507	heart development
# somatotroph_axis = BP	GO:0070848	response to growth factor
# BMP = BP	GO:0071772	response to BMP

top10LAKI <- dplyr::arrange(tab, P.value.adj.LAKI_KOvWT)[1:10, "AptName"]
top20LAKI <- dplyr::arrange(tab, P.value.adj.LAKI_KOvWT)[1:20, "AptName"]
top10FACE <- dplyr::arrange(tab, P.value.adj.ZMPSTE_KOvWT)[1:10, "AptName"]
top20FACE <- dplyr::arrange(tab, P.value.adj.ZMPSTE_KOvWT)[1:20, "AptName"]
top10KOvsWT <- dplyr::arrange(tab, P.value.adj.KOvWT)[1:10, "AptName"]
top20KOvsWT <- dplyr::arrange(tab, P.value.adj.KOvWT)[1:20, "AptName"]

top20LAKI_up <- dplyr::filter(tab, Coef.LAKI_KOvWT > 0) %>% 
                dplyr::arrange(P.value.adj.LAKI_KOvWT) %>% 
                head(20) %>% 
                dplyr::pull(AptName)
top20LAKI_down <- dplyr::filter(tab, Coef.LAKI_KOvWT < 0) %>% 
                  dplyr::arrange(P.value.adj.LAKI_KOvWT) %>% 
                  head(20) %>% 
                  dplyr::pull(AptName)

top20FACE_up <- dplyr::filter(tab, Coef.ZMPSTE_KOvWT > 0) %>% 
                dplyr::arrange(P.value.adj.ZMPSTE_KOvWT) %>% 
                head(20) %>% 
                dplyr::pull(AptName)
top20FACE_down <- dplyr::filter(tab, Coef.ZMPSTE_KOvWT < 0) %>% 
                  dplyr::arrange(P.value.adj.ZMPSTE_KOvWT) %>% 
                  head(20) %>% 
                  dplyr::pull(AptName)

top20KOvsWT_up <- dplyr::filter(tab, Coef.KOvWT > 0) %>% 
                  dplyr::arrange(P.value.adj.KOvWT) %>% 
                  head(20) %>% 
                  dplyr::pull(AptName)
top20KOvsWT_down <- dplyr::filter(tab, Coef.KOvWT < 0) %>% 
                    dplyr::arrange(P.value.adj.KOvWT) %>% 
                    head(20) %>% 
                    dplyr::pull(AptName)


ECM_LAKI_list <- unlist(stringr::str_split("MFAP4/FMOD/DPT/COL2A1/DAG1/TNXB/FLRT2/ASPN/SERPINF1/PKM/EFEMP1/CD248/CCDC80/CLEC3B/SDC3/COL1A1/IGFALS/COL3A1/FBLN5/TGFBI/ANOS1/HAPLN1/FBLN7/NDNF/HTRA1/WNT5B/SMOC2/VTN/SOST/COL9A3/CASK/MFAP2/LRRC32/TIMP4/SFRP2/F2/CTSD/MATN2/ANGPTL1/NPNT/WNT7A/IHH/CILP/TNC/ENTPD2/APOH/F9/NPPA/A2M/AGRN",
                                      pattern = "/"))
ECM_FACE_list <- unlist(stringr::str_split("FMOD/DPT/DAG1/CD248/ASPN/THBS3/TNXB/CCDC80/CLEC3B/NID2/EFEMP1/VTN/COCH",
                                      pattern = "/"))
cardiovascular_LAKI_list <- unlist(stringr::str_split("IGF1/COL2A1/FLRT2/SNAI2/NOTCH1/DSG2/GSK3A/RBP4/MAPK14/MAPK3/JAG1/GRK2/COL3A1/NFATC4/EGFR/NRG1/FGFRL1/NDST1/SFRP2/PDLIM3/FLNA/BMPR1A/DLL1/CERT1/SRI/PKD2/POU6F1/ROBO2/FOXC2/STK4/IHH/FGF8/NOG/FKBP1A/ANKRD1/NPPA/MTHFD1",
                                           pattern = "/"))
somatotroph_axis_LAKI_list <- unlist(stringr::str_split("FMOD/CBL/COL2A1/FLRT2/SNAI2/NOTCH1/ASPN/RGMA/CDH5/TWSG1/CLEC3B/COL1A1/FLT4/MSTN/GREM2/MAPK14/HJV/MAPK3/FGFR1/AKT1/COL3A1/SRC/PDGFD/NDNF/HTRA1/PTP4A3/EHD1/MAPKAPK2/GAS1/SMOC2/EGFR/VTN/FGFRL1/SOST/VEGFD/CASK/ZFYVE27/ROR2/NDST1/FSTL3/LRRC32/SFRP2/BMPR1A/DLL1/ADAM17/PRKCB/NPNT/WNT7A/CILP/MICOS10-NBL1/FGF8/NOG/FKBP1A/ANKRD1/GIPC1/BMP15/PELO/GDF2",
                                           pattern = "/"))
BMP_LAKI_list <- unlist(stringr::str_split("COL2A1/NOTCH1/RGMA/CDH5/TWSG1/GREM2/HJV/MAPK3/HTRA1/SOST/ROR2/FSTL3/SFRP2/BMPR1A/MICOS10-NBL1/NOG/BMP15/PELO/GDF2",
                                           pattern = "/"))

# proteins_plots <- data.frame(AptName = names(mice_uniovi)[!names(mice_uniovi) %in% sample_info]) %>% 
#                   dplyr::left_join(analytes[, c("AptName", "EntrezGeneID", "EntrezGeneSymbol")]) %>% 
#                   dplyr::mutate(top10LAKI = ifelse(AptName %in% top10LAKI, "X", ""),
#                                 top20LAKI = ifelse(AptName %in% top20LAKI, "X", ""),
#                                 top20LAKI_up = ifelse(AptName %in% top20LAKI_up, "X", ""),
#                                 top20LAKI_down = ifelse(AptName %in% top20LAKI_down, "X", ""),
#                               
#                                 top10FACE = ifelse(AptName %in% top10FACE, "X", ""),
#                                 top20FACE = ifelse(AptName %in% top20FACE, "X", ""),
#                                 top20FACE_up = ifelse(AptName %in% top20FACE_up, "X", ""),
#                                 top20FACE_down = ifelse(AptName %in% top20FACE_down, "X", ""),
#                                 
#                                 top10KOvsWT = ifelse(AptName %in% top10KOvsWT, "X", ""),
#                                 top20KOvsWT = ifelse(AptName %in% top20KOvsWT, "X", ""),
#                                 top20KOvsWT_up = ifelse(AptName %in% top20KOvsWT_up, "X", ""),
#                                 top20KOvsWT_down = ifelse(AptName %in% top20KOvsWT_down, "X", ""),
#                                 
#                                 ECM_LAKI = ifelse(EntrezGeneSymbol %in% ECM_LAKI_list, "X", ""),
#                                 ECM_FACE = ifelse(EntrezGeneSymbol %in% ECM_FACE_list, "X", ""),
#                                 cardiovascular_LAKI = ifelse(EntrezGeneSymbol %in% cardiovascular_LAKI_list, "X", ""),
#                                 somatotroph_axis_LAKI = ifelse(EntrezGeneSymbol %in% somatotroph_axis_LAKI_list, "X", ""),
#                                 BMP_LAKI = ifelse(EntrezGeneSymbol %in% BMP_LAKI_list, "X", "")
#                                 ) %>% 
#                                 dplyr::filter(if_any(!all_of(c("AptName", "EntrezGeneID", "EntrezGeneSymbol")), ~ grepl("X",.))) %>% 
#                                 dplyr::left_join(tab[, c()], by = "AptName",)

proteins_plots2 <-tab %>% 
                    dplyr::select(c("AptName", "EntrezGeneID", "EntrezGeneSymbol", "Coef.KOvWT", "Coef.LAKI_KOvWT", "Coef.ZMPSTE_KOvWT",
                                    "P.value.adj.KOvWT", "P.value.adj.LAKI_KOvWT", "P.value.adj.ZMPSTE_KOvWT")) %>% 
                    dplyr::mutate(top10LAKI = ifelse(AptName %in% top10LAKI, "X", ""),
                                  top20LAKI = ifelse(AptName %in% top20LAKI, "X", ""),
                                  top20LAKI_up = ifelse(AptName %in% top20LAKI_up, "X", ""),
                                  top20LAKI_down = ifelse(AptName %in% top20LAKI_down, "X", ""),
                                  
                                  top10FACE = ifelse(AptName %in% top10FACE, "X", ""),
                                  top20FACE = ifelse(AptName %in% top20FACE, "X", ""),
                                  top20FACE_up = ifelse(AptName %in% top20FACE_up, "X", ""),
                                  top20FACE_down = ifelse(AptName %in% top20FACE_down, "X", ""),
                                  
                                  top10KOvsWT = ifelse(AptName %in% top10KOvsWT, "X", ""),
                                  top20KOvsWT = ifelse(AptName %in% top20KOvsWT, "X", ""),
                                  top20KOvsWT_up = ifelse(AptName %in% top20KOvsWT_up, "X", ""),
                                  top20KOvsWT_down = ifelse(AptName %in% top20KOvsWT_down, "X", ""),
                                  
                                  ECM_LAKI = ifelse(EntrezGeneSymbol %in% ECM_LAKI_list, "X", ""),
                                  ECM_FACE = ifelse(EntrezGeneSymbol %in% ECM_FACE_list, "X", ""),
                                  cardiovascular_LAKI = ifelse(EntrezGeneSymbol %in% cardiovascular_LAKI_list, "X", ""),
                                  somatotroph_axis_LAKI = ifelse(EntrezGeneSymbol %in% somatotroph_axis_LAKI_list, "X", ""),
                                  BMP_LAKI = ifelse(EntrezGeneSymbol %in% BMP_LAKI_list, "X", "")
                    ) %>% 
                    dplyr::filter(if_any(!all_of(c("AptName", "EntrezGeneID", "EntrezGeneSymbol")), ~ grepl("X",.)))
                    
                

write.xlsx(proteins_plots2, "proteins_plots.xlsx")
