AddUpGenes<-function(dataup,tree,mean_up,sd_up,norm,purity,sim_num,filename){
tree <- read.csv(tree, header=T, sep='\t')
data<-read.csv(dataup,sep="\t",row.names = 1)
gene_num<-length(data[,1])
node.final<-matrix(NA,sim_num, gene_num)
rate.final<-matrix(NA,sim_num, gene_num)
node_num<-10
prop.vec<-rep(NA,node_num)
for(i in 1:node_num){
  level1<-tree$gene[i]
  b<-as.character(tree$phi[i])
  prop<-as.numeric(substr(b,2,nchar(b)-2))
  prop.vec[i]<-prop
}
prop.vec<-prop.vec[!is.na(prop.vec)]
prop.vec<-prop.vec[prop.vec>=0.1]
likely.total<-rep(NA,sim_num)
rate1<-rnorm(gene_num,mean_up,sd_up)
rate<-exp(rate1)
for (c in 1:sim_num){
  rate.final[c,]<-rate #original set a rate
  node<-rep(NA,gene_num)
  # print(c(c,rate))
  for (k in 1:gene_num){   # for each gene get its node
    likely<-rep(NA,length(prop.vec))
    lambda1<-((1-prop.vec*purity)+rate[k]*prop.vec*purity)*norm
    for(j in 1:length(prop.vec)){
      mu<- data[k,1]*lambda1[j]
      logfaca <- 0
      for (i in 1:data[k,2])
      {
        logfaca <-logfaca+log(i)  
      }
      likely[j]<--mu+data[k,2]*log(mu)-logfaca
    }
    #   print(c(c,likely))
    node[k]<-which(likely == max(likely), arr.ind = TRUE)[1]
  }
  node.final[c,]<-node
  rate1.new<-rnorm(gene_num,log(rate),0.8)
  rate.new<-exp(rate1.new)
  prop.node<-prop.vec[node]
  likely.new<-likely.old<-likely.llh<-rep(NA,gene_num)  # for each gene update its rate
  for (k in 1:gene_num){
    lambda1<-((1-prop.node[k]*purity)+rate[k]*prop.node[k]*purity)*norm
    lambda1.new<-((1-prop.node[k]*purity)+rate.new[k]*prop.node[k]*purity)*norm
    mu.old<- data[k,1]*lambda1
    mu.new<- data[k,1]*lambda1.new
    logfaca <- 0
    for (i in 1:data[k,2])
    {
      logfaca <-logfaca+log(i)  
    }
    likely.old[k]<--mu.old+data[k,2]*log(mu.old)-logfaca+log(dnorm(log(rate[k]),mean_up,sd_up))
    likely.new[k]<--mu.new+data[k,2]*log(mu.new)-logfaca+log(dnorm(log(rate.new[k]),mean_up,sd_up)) 
    
    if ((likely.new[k]-likely.old[k])>log(runif(1))){
      rate[k]<-rate.new[k]
    }
  }
  for (k in 1:gene_num){
    lambda.llh<-((1-prop.node[k]*purity)+rate[k]*prop.node[k]*purity)*norm
    mu.llh<- data[k,1]*lambda.llh
    
    logfaca <- 0
    for (i in 1:data[k,2])
    {
      logfaca <-logfaca+log(i)  
    }
    likely.llh[k]<--mu.llh+data[k,2]*log(mu.llh)-logfaca
  }
  likely.total[c]<-sum(likely.llh)
}
#variable<- apply(likely, 1, which.max)
#node[m,l]<-as.numeric(names(sort(table(variable),decreasing=TRUE)[1]))
print(which(likely.total == max(likely.total), arr.ind = TRUE)[1])
node.idx<-  which(likely.total == max(likely.total), arr.ind = TRUE)[1]
#write.csv(cbind(node.final,likely.total),file=paste("uptree/node.final",name, id,".csv",sep=""))
#write.csv(rate.final,file="rate.final.csv",sep=""))
save(node.final,node.idx,likely.total,rate.final,file=paste(filename,"_result.rdata",sep=""))
}
AddUpGenes("dataup.txt","tree.txt",2,0.8,1,1,50,"data")
