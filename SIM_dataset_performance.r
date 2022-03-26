
source("PCR_HAT_dot_return_model.r")
source("PCR_HAT_eigen_return_model.r")
source("PCR_HAT_r2_return_model.r")

###########################################################	
###########################################################
## SIM dataset
phe <- read.table(file="sim.phe.txt",header=TRUE)                    
gen <- read.csv(file="sim.gen.txt",header=TRUE)

yused=1
#yused=2
#yused=3

y=as.numeric(as.character(phe[,yused]))
sample_size=length(y)

G=as.matrix(t(gen))
G[G==1]=2
G[G==0]=1
G[G==-1]=0

X=apply(G,2,function(ggg){ f=mean(ggg)/2
                           (ggg-f*2)/sqrt(2*f*(1-f))
                         })

print (yused)
###########################################################	
###########################################################
yx=cbind(y,X)  
###########################################################	
###########################################################
k = 5    ##cross validation fold
cv_fold = rep(1:k,ceiling(sample_size/k))[1:sample_size] 

###########################################################	
###########################################################
for(cvseed in 1:10){

set.seed(cvseed)
ff = sample(cv_fold,sample_size)
#write.csv(ff, paste("cvseed_",cvseed,"cv_fold.csv"))	##write seed cv

cvlist <- lapply( 1:k, function(cvcvcv){which(ff==cvcvcv)} )     


    R2=lapply(c(1:k),function(parts_parts){                          
		                    train_yx = yx[-(cvlist[[parts_parts]]),]   ##train or Ref data
				            test_yx = yx[cvlist[[parts_parts]],]      ##test data
			  				  
	    ry=train_yx[,1]
	    rx=train_yx[,-1]
	    tx=test_yx[,-1]
	    ty=test_yx[,1]
###########################################################	
###########################################################
##      PCR-DOT
###########################################################		
	    PCR_HAT_dot=PCR_HAT_dot_input_rx_ry(rx,ry)
	   
	           v=PCR_HAT_dot[[1]]
	         pos=PCR_HAT_dot[[2]]
	        bata=PCR_HAT_dot[[3]]
	
	    tx_pc=tx%*%v
        pre_y=bata[1]+tx_pc[,pos]%*%as.matrix(bata[-1])
	
	    ability_dot=cor(pre_y,ty)
	     pcnums_dot=length(pos)
###########################################################
###########################################################
##      PCR-eigen
###########################################################		
	    PCR_HAT_eigen=PCR_HAT_eigen_input_rx_ry(rx,ry)
	   
	           v=PCR_HAT_eigen[[1]]
	         pos=PCR_HAT_eigen[[2]]
	        bata=PCR_HAT_eigen[[3]]
	
	    tx_pc=tx%*%v
        pre_y=bata[1]+tx_pc[,pos]%*%as.matrix(bata[-1])
	
	    ability_eigen=cor(pre_y,ty)
	     pcnums_eigen=length(pos)
###########################################################
###########################################################
##      PCR-SS
###########################################################		
	    PCR_HAT_r2=PCR_HAT_r2_input_rx_ry(rx,ry)
	   
	           v=PCR_HAT_r2[[1]]
	         pos=PCR_HAT_r2[[2]]
	        bata=PCR_HAT_r2[[3]]
	
	    tx_pc=tx%*%v
        pre_y=bata[1]+tx_pc[,pos]%*%as.matrix(bata[-1])
	
	    ability_r2=cor(pre_y,ty)
	     pcnums_r2=length(pos)
###########################################################
###########################################################
###########################################################
	result=data.frame(	                  
					  ability_dot,  
					   pcnums_dot,
					
					  ability_eigen,
					   pcnums_eigen,
					  
					  ability_r2,
					   pcnums_r2				                   
					)
    return(result)
    })

    #R2
OUT= matrix( unlist(R2) ,6,5)
print(OUT)

if(cvseed == 1){OUTS= matrix(c("ability_dot","pcnums_dot","ability_eigen","pcnums_eigen","ability_r2","pcnums_r2"),6,1)}
OUTS=cbind(OUTS,OUT)

}
results=OUTS
write.csv(results,row.names=FALSE,paste("pcr_SIM_r2_y",yused,".csv",sep="_"))
