readComp = function(f){
  r = read.table(f,skip = 5,header = T)
  l = as.numeric(sapply(strsplit(readLines(f)[2:3],' '),'[',3))
  r$ncells1=l[1]
  r$ncells2=l[2]
  r
}

castXYtable = function(x,y,i){
  ys = sort(unique(y))
  xs = sort(unique(x))
  m = matrix(NA,nrow=length(xs),ncol=length(ys),dimnames=list(xs,ys))
  x = as.character(x)
  y = as.character(y)
  for(j in 1:length(x))
    m[x[j],y[j]]= i[j]
  m
}


loadAllComps = function(names,path_shared,path_soc){
  res = NULL
  for(n1 in names){
    for(n2 in names){
      r = readComp(paste0(path_shared,'/map.',n1,'-',n2))
      r$name1=n1
      r$name2=n2
      res = rbind(res,r)
    }
  }
  res$name1_ = paste0(res$name1,'|',res$experiment1_cluster)
  res$name2_ = paste0(res$name2,'|',res$experiment2_cluster)
  mtx=castXYtable(res$name1_,res$name2_,res$loss)
  clsizes = NULL
  for(n in names){
    cl = read.table(paste0(path_soc,'/',n,'/clusters.tsv'),header = TRUE)
    clsizes = rbind(clsizes,data.frame(name=n,group='doublet',count=sum(cl$status=='doublet')))
    clsizes = rbind(clsizes,data.frame(name=n,group='unassigned',count=sum(cl$status=='unassigned')))
    for(cln in sort(unique(cl$assignment[cl$status=='singlet']))){
      clsizes = rbind(clsizes,data.frame(name=n,group=cln,count=sum(cl$status=='singlet' & cl$assignment==cln)))
    }
    
  }
  #clsizes = reshape::cast(clsizes,name ~ group,value = 'count')
  list(res=res,mtx=as.matrix(mtx),clsizes=clsizes)
}

args = commandArgs(trailingOnly = TRUE)

if(!(length(args) %in% c(3,4)))
  stop("Please provide free four commandline arguments:
        path_to_shared_sample_folder path_to_souporcell_folder paths_to_samplefile [expected_number_of_genotypes]")

sids = read.table(args[3])
cmp = loadAllComps(sids$V1,args[1],args[2])
rownames(cmp$clsizes) = paste0(cmp$clsizes$name,'|',cmp$clsizes$group)
if(length(args)==3){
  library(cluster)
  gap.stat = clusGap(cmp$mtx, FUNcluster = function(d,k){
    hcl = hclust(as.dist(d))
    list(cluster=cutree(hcl,k))
  }, K.max = nrow(cmp$mtx)-1)
  gap = gap.stat$Tab[,'gap']
  # look for first maximum
  ngenotypes = NA
  if(length(gap)>2){
    l = length(gap)
    inx = which(gap[-c(1,l)] > gap[-c(l-1,l)] & gap[-c(1,l)] > gap[-1:-2])
    if(length(inx)>0)
      ngenotypes = min(inx)
  }
  if(is.na(ngenotypes))
    ngenotypes = which.max(gap)
}else
  ngenotypes = as.integer(args[4])
hcl = hclust(as.dist(cmp$mtx))
cl = cutree(hcl,k=ngenotypes)
cmp$clsizes$genotype_cluster = cl[rownames(cmp$clsizes)]

set.seed(1)
cols = sample(colors(distinct=TRUE))
ann = cmp$clsizes[colnames(cmp$mtx),]

# save results
dir.create('group_samples',recursive = TRUE,showWarnings = FALSE)
size = 5 + nrow(cmp$mtx)*0.5
pdf('group_samples/shared_samples_loss_heatmap.pdf',w=10,h=10)
heatmap(cmp$mtx,symm = T,distfun = function(x)as.dist(x),margins = c(20,20),col=rev(heat.colors(100)),
        ColSideColors = cols[as.integer(factor(ann$name))],
        RowSideColors = cols[as.integer(factor(ann$genotype_cluster))])
dev.off()


write.csv(cmp$mtx,'group_samples/shared_samples_loss.csv')
write.csv(cmp$clsizes,'group_samples/shared_samples_clusters.csv')


