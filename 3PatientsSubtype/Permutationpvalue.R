permutaionSinglePT<-function(dataup,filename,tree){
  tree <- read.csv(tree, header=T, sep='\t')
  data<-read.csv(dataup,sep="\t",row.names = 1)
  load(file=paste(filename,"_result.rdata",sep=""))
  node.llh<-node.final[node.idx,]
gene.up<-rownames(data)
node_num<-15
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
geneinnode<-list()
for(i in 1:node_num){
  genenumber<-gene.up[which(node.llh == i)]
  level1<-tree$gene[i]
  genenumber2<-as.character(strsplit(as.character(level1), "; ")[[1]])
  geneinnode[[i]]<-c(genenumber,genenumber2)
}
order<-list()
order[[1]]<-0
if(node_num>1){
  for(i in 2:node_num){
    num<-tree$parent[i]
    m<-as.numeric(as.character(num))+1  
    aaa<-m
    while (m>1){
      m<-as.numeric(as.character(tree$parent[m]))+1
      aaa<-c(m,aaa)
    }
    order[[i]]<-aaa
  }
}
save(geneinnode,order,file=paste("geneinnode",filename,'.rdata'))
}

temp<-function(gene1,gene2,geneinnode,order){
  geneontree<-unlist(geneinnode)
  if (is.element(gene1, geneontree) & sum(gene2==geneontree)==0){
    prob<-1 
  }else if(is.element(gene2, geneontree) & sum(gene1==geneontree)==0){
    prob<- -1 
  }else if(is.element(gene1, geneontree) & is.element(gene2, geneontree)){
    for (i in 1:length(geneinnode)){
      if(is.element(gene1,geneinnode[[i]])){
        node_pos1<-i
      }
    } 
    for (i in 1:length(geneinnode)){
      if(is.element(gene2,geneinnode[[i]])){
        node_pos2<-i
      }
    }
    if(is.element(node_pos2,order[[node_pos1]])){
      prob<- -1
    }else if(is.element(node_pos1,order[[node_pos2]])){  ####gene1 is the root of gene2
      prob<- 1
    }else{
      prob<-0
    }
  }else (prob<-0)
}
permutation<-function(path1,path2,geneinnode,order,NoofPerm,seed){
len_path1<-length(path1)
len_path2<-length(path2)
permutation<-rep(NA,NoofPerm)
set.seed(seed)
for (p in 1:NoofPerm){
path1_perm<-sample(c(path1,path2),length(path1),replace=F)
path2_perm<-setdiff(c(path1,path2),path1_perm)
prop<-rep(NA,length(path1_perm)*length(path2_perm))
k<-1
for (i in 1:length(path1_perm)){
  for (j in 1:length(path2_perm)){
    prop[k]<-temp(path1_perm[i],path2_perm[j],geneinnode,order)
    k<-k+1 
  }
}
permutation[p]<-sum(prop)
}
return(permutation)
}

stat_score<-function(path1,path2,geneinnode,order){
  len_path1<-length(path1)
  len_path2<-length(path2)
  prop<-rep(NA,len_path1*len_path2)
  k<-1
  for (i in 1:length(path1)){
    for (j in 1:length(path2)){
      prop[k]<-temp(path1[i],path2[j],geneinnode,order)
      k<-k+1 
    }
  }
 return(sum(prop))
}

load(paste("geneinnode",filename,'.rdata'))
load("PathBreast.rdata")
l<-1
per<-NULL
score<-res_per<-rep(NA)
NoofPerm<-10
for (m in 1:(length(path)-1)){
  for (n in (m+1):length(path)){
      per[[l]]<-permutation(path[[m]],path[[n]],geneinnode,order,NoofPerm=NoofPerm,seed=1234)
      score[l]<-stat_score(path[[m]],path[[n]],geneinnode,order)
      res_per[l]<-sum(per[[l]]>score[l])/NoofPerm
      l<-l+1
      print(l)
  }
}
save(per,score,res_per,file=paste("permutation",filename,'.rdata'))
