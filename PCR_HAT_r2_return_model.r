###############################################################################################
###############################################################################################
## refX dimensionality: n*m ;  n denotes sample size; m denotes variables (SNPs)
given_refX_get_pc_D_v=function(refX){
                                      X=refX
									  sample_size=nrow(X)
									  
                                      e=eigen(X%*%t(X))
                                      pc=e$vectors[,-sample_size]%*%(diag(sqrt(e$values[-sample_size])))
                                      
									  #pc[1:3,1:3]                                      
									  D=e$values[-sample_size]                                      
									  V=t(X)%*%e$vectors[,-sample_size]%*%diag(sqrt(1/D))
return(list(pc,
            D,
			V
            ))
}

###############################################################################################
###############################################################################################
## LOO cross-validation
LOO_hat_PDP=function(pc,D,ry){                           
						      						   
						      solve_d=diag(1/D,length(D),length(D))
						   
						      H=pc%*%solve_d%*%t(pc)
						   
						      error=H%*%ry-ry
						   
						  					 
                              press=(error/(1-diag(H)))^2
                              PRESS<-sum(press)
                    
					        return(PRESS)
}



###############################################################################################
###############################################################################################
##  provide rx ry, return the best PCs and model
PCR_HAT_r2_input_rx_ry=function(rx,ry){

    eigen.rx=given_refX_get_pc_D_v(rx)
    
	pc.rx=eigen.rx[[1]]
	D.rx=eigen.rx[[2]]
	V.rx=eigen.rx[[3]]
	
	
	rx=pc.rx


	
##  0.90 variance
##	
rxv=apply(rx,2,var)
	
	for(i in 1: ncol(rx) ){ vv=sum(rxv[1:i])/sum(rxv)
                            if(0.90<= vv ){break}
    					  }
	print(i) ## the max PC numbers
	v90=i
	
	rx=rx[,1:v90]


##  PCR_LOO_r2	
##	     	
r2=apply(rx,2,function(ppp){ (cor(ppp,ry))^2  }    )	   
			
				
##  rank: r2   	 
##	     	
rankrfx=rank(-r2)
   
	   
##  the best PCs	 
##	   	   
press=NULL 
   
    for(i in 1:ncol(rx)){
                      pos=which(rankrfx<=i)
                      press[i]=LOO_hat_PDP(rx[,pos],D.rx[pos],ry)
	                }
     
	best_pcs_press=which.min(press)
    #min(press)

    pos=which(rankrfx<=best_pcs_press)
    length(pos)
    #pos 
	
##  the PCR_r2_model	 
##	  
pc_bata=summary(lm(ry~rx[,pos]))$coefficients[,1]
	    
return(list(            
            V.rx,
			pos,
			pc_bata	   
		    ))
}

##  prediction	 
##
##  tx_pc=tx%*%V.rx	  
##  pre_y=bata[1]+tx_pc[,pos]%*%as.matrix(bata[-1])