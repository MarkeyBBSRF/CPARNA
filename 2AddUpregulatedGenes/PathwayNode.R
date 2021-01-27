pathway<-function(filename,tree,path){
load(file=paste(filename,"_result.rdata",sep=""))
load(file=path)
tree <- read.csv(tree, header=T, sep='\t')
node.llh<-node.final[node.idx,]
gene.up<-rownames(data)
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
node_num<-length(prop.vec)
match.up1<-matrix(NA,length(path),node_num)
GeneinNodeup<-rep(NA,node_num)  
for(i in 1:node_num){
  for (j in 1:length(path)){
  genenumber<-gene.up[which(node.llh == i)]
  GeneinNodeup[i]<-length(genenumber)
  match.up1[j,i]<-sum(genenumber %in% path[[j]])
  }
}
match.n<-matrix(NA,length(path),node_num)
GeneinNodedn<-rep(NA,node_num) 

for(i in 1:node_num){
  for (j in 1:length(path)){
  level1<-tree$gene[i]
  genenumber<-as.character(strsplit(as.character(level1), "; ")[[1]]) 
  GeneinNodedn[i]<-length(genenumber)
  match.n[j,i]<-sum(genenumber %in% path[[j]])
  }
}
#save(match.n,file=paste("matrixdown",name,'.rdata'))
GeneinNode<-GeneinNodedn+GeneinNodeup
match.t<-match.n+match.up1
res.p<-matrix(NA,length(path),node_num)
################  
for(m in 1:length(path)){
  for (i in 1:node_num){
    a11<-match.t[m,i]
    a12<-GeneinNode[i]-a11
    a21<-length(path[[m]])-a11
    a22<-length(unlist(path))-a11-a12-a21
    res.p[m,i]<-fisher.test(matrix(c(a11,a21,a12,a22),ncol=2),alternative="greater")$p.value
  }
}
round(res.p,4)
res<-rbind(res.p,match.n,match.n/rowSums(match.n))
save(res.p,match.n,file=paste(filename,"_pathwayRes.rdata",sep=""))
}

pathway("data","tree.txt","PathBreast.rdata")
