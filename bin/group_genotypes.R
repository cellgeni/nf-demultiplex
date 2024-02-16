readComp = function(f){
  r = read.table(f,skip = 5,header = T)
  l = as.numeric(sapply(strsplit(readLines(f)[2:3],' '),'[',3))
  r$ncells1=l[1]
  r$ncells2=l[2]
  r
}

loadAllComps = function(names,path){
  res = NULL
  for(n1 in names){
    for(n2 in names){
      r = readComp(paste0(path,'/shared_samples/map.',n1,'-',n2))
      r$name1=n1
      r$name2=n2
      res = rbind(res,r)
    }
  }
  res$name1_ = paste0(res$name1,'|',res$experiment1_cluster)
  res$name2_ = paste0(res$name2,'|',res$experiment2_cluster)
  mtx=reshape::cast(res,name1_~name2_,value='loss')
  clsizes = NULL
  for(n in names){
    cl = read.table(paste0(path,'/',n,'/clusters.tsv'),header = TRUE)
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

if(length(args)!=3 || !file.exists(args[1]))
  stop("Please provide three commandline arguments:
       paths_to_samplefile path_to_souporcell_folder expected_number_of_genotypes")

sids = read.table(args[1])
cmp = loadAllComps(sids$V1,paste0(args[2]))
rownames(cmp$clsizes) = paste0(cmp$clsizes$name,'|',cmp$clsizes$group)

hcl = hclust(as.dist(cmp$mtx))
cl = cutree(hcl,k=as.integer(args[3]))
cmp$clsizes$genotype_cluster = cl[rownames(cmp$clsizes)]

set.seed(1)
cols = sample(colors(distinct=TRUE))
ann = cmp$clsizes[colnames(cmp$mtx),]
pdf('shared_samples_loss_heatmap.pdf',w=10,h=10)
heatmap(cmp$mtx,symm = T,distfun = function(x)as.dist(x),margins = c(10,10),col=rev(heat.colors(100)),
        ColSideColors = cols[as.integer(factor(ann$name))],
        RowSideColors = cols[as.integer(factor(ann$genotype_cluster))])
dev.off()

# save results
write.csv(cmp$mtx,'shared_samples_loss.csv')
write.csv(cmp$clsizes,'shared_samples_clusters.csv')


