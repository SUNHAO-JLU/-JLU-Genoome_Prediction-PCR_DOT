# Genome_Prediction-PCR_DOT

Usage:
R code:
  
  # Reference data 
  source("PCR_HAT_dot_return_model.r")
  PCR_HAT_dot=PCR_HAT_dot_input_rx_ry(reference_data_x,reference_data_y)
  
  v=PCR_HAT_dot[[1]]      #   x=udv ; p=xv
  pos=PCR_HAT_dot[[2]]    #  the postion of selected PCs
  bata=PCR_HAT_dot[[3]]   #  the coefficients

  # Prediction  
  tx_pc=test_data_x%*%v
  predict_y=bata[1]+tx_pc[,pos]%*%as.matrix(bata[-1])
