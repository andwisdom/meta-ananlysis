study.names<-c("GSE1297", "GSE5281") #a vector of study`names
n<-length(study.names) #n is the number of studies

library("MetaDE")
AD.raw<- MetaDE.Read(study.names, skip=rep(1,n), via = "txt", matched=T, log=F )
AD.merged<- MetaDE.merge(AD.raw)
w<- MetaDE.filter(AD.merged,c(0.0001,0.0001))
meta.res.FEM<-MetaDE.rawdata(x=w,ind.method = rep('modt',n),meta.method = "FEM",nperm = 1000,ind.tail="high",paired=rep(FALSE,n))
meta.res.REM<-MetaDE.rawdata(x=w,ind.method = rep('modt',n),meta.method = "REM",nperm = 1000,ind.tail="high",paired=rep(FALSE,n))


tau2<-meta.res.REM$meta.analysis$tau2
sigma2<-meta.res.FEM$ind.Var
Q<-meta.res.REM$meta.analysis$Qval
ind.ES<-meta.res.FEM$ind.ES
ES_FEM<-meta.res.FEM$meta.analysis$mu.hat
ES_REM<-meta.res.REM$meta.analysis$mu.hat

tau2<-data.matrix(tau2)
sigma2<-data.matrix(sigma2)
Q<-data.matrix(Q)
ind.ES<-data.matrix(ind.ES)
ES_FEM<-data.matrix(ES_FEM)
ES_REM<-data.matrix(ES_REM)


D2<-c(0) #Set initial value 

n<-nrow(ind.ES)  #Calculate the number of genes 
m<-ncol(ind.ES)   #Calculate the number of study

for (i in 1:n){
	tau2_e<-tau2[i]
	sigma2_e<- sigma2[i,]
	Q_e<-Q[i]
	ind.ES_e<-ind.ES[i,]
	ES_FEM_e<-ES_FEM[i]
	ES_REM_e<-ES_REM[i]
	#calcute S_MM,S_1,S_2,S_3,_S_4
	S_MM_e<-0
	S_1<-0
	S_2<-0
	S_3<-0
	S_4<-0
	for (j in 1:m){
		wt_1<-1/(sigma2_e[j])
		wt_2<-1/(sigma2_e[j] + tau2_e)
		S_MM_e<-S_MM_e + (wt_2)*((ind.ES_e[j]-ES_REM_e)**2)
		S_1<-S_1 + wt_1
		S_2<-S_2 + (wt_1)**2
	}	
    	#calcute D2
	D2_e<-(Q_e-S_MM_e)/((S_1)-(S_2/S_1))
	if (tau2_e==0) D2_e<-0
	if (D2_e>0) D2_e<-D2_e
	if (D2_e<0) D2_e<-0
	D2<-c(D2,D2_e)
	##
	#calcute K2	
}
D2<-D2[-1]




D2<-as.vector(D2)

sigma2_add_D2<-sigma2 + D2



sigma2_add_D2<-as.matrix(sigma2_add_D2)

ES<-meta.res.FEM$ind.ES

D2<-list(ES=ES,Var=sigma2_add_D2)
D2_res<-MetaDE.ES(x=D2,meta.method='FEM')
p_value<-D2_res$pval
View(p_value)
class(p_value)
p_value2<-p_value[order(p_value,decreasing = FALSE)]
p_value<-as.matrix(p_value2)
View(p_value)
D2_res$pval<-p_value
View(D2_res$pval)
FDR<-D2_res$FDR


FDR2<-as.vector(FDR)
FDR2<-cbind(rownames(FDR),FDR)
FDR<-as.matrix(FDR)
FDR3<-FDR2[order(FDR2[,2], decreasing = FALSE)]
FDR3<-as.matrix(FDR3)
FDR2<-FDR[order(FDR[,1], decreasing = FALSE)]
FDR2<-as.matrix(FDR2)
FDR<-cbind(FDR3, FDR2)
D2_res$FDR<-FDR

