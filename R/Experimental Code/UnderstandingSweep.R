data<-matrix(seq(1,12),ncol=4,nrow=3,byrow=TRUE)
data

colMean <- apply(data,2,mean)
colMean

rowMean<-apply(data,1,mean)

myfunc<-function(x){
    t=1
    for (i in x){t=t*i}
    return(t)
}

multCol<-apply(data,2,myfunc)
multCol

data
sweep(data,2,c(3,4,5,6),"-")
