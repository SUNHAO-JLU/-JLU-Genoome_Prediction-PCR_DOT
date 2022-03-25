# Genome_Prediction-PCR_DOT
We developed a method PCR-DOT to perform principal components regression based on the dot product theory. We compared the predictive ability of different approaches based on the simulated and real genomic data, and in general, the prediction ability of PCR-dot was better than the previous PCR-eigen and PCR-SS methods.
# Usage: R code:
# Reference data : PCR model
source("PCR_HAT_dot_return_model.r")
PCR_HAT_dot=PCR_HAT_dot_input_rx_ry(reference_data_x,reference_data_y)
  
v=PCR_HAT_dot[[1]]      #   x=udv ; p=xv
pos=PCR_HAT_dot[[2]]    #  the postion of selected PCs
bata=PCR_HAT_dot[[3]]   #  the coefficients

# Prediction  
tx_pc=test_data_x%*%v
predict_y=bata[1]+tx_pc[,pos]%*%as.matrix(bata[-1])
