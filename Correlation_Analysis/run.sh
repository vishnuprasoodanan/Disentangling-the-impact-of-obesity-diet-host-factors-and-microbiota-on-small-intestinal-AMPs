R CMD BATCH Corr_Exp_Wise.R
mkdir Exp_Wise_Correlation Diet_Wise_Correlation Source_Wise_Correlation
mv Si8_FP0* Exp_Wise_Correlation
mv WSD_* Diet_Wise_Correlation
mv Chow_* Diet_Wise_Correlation
mv HFD_* Diet_Wise_Correlation 
mv Content* Source_Wise_Correlation
mv MAB* Source_Wise_Correlation
