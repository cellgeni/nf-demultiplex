library(visutils)
samplefile = commandArgs(trailingOnly=TRUE)
samples = read.table(samplefile,header = F)
colnames(samples) = c('sid','bam_path','barcodes_path','n')


soc = lapply(seq_len(nrow(samples)),function(i)read.table(paste0('souporcell/',samples$sid[i],'/clusters.tsv'),header = T,row.names = 1))

# vireo cannot be run with k=1, so lets put NULL on these places
vir = lapply(seq_len(nrow(samples)),function(i){
  fn = paste0('vireo/',samples$sid[i],'-vireo/donor_ids.tsv')
  if(!file.exists(fn))
    return(NULL)
  read.table(fn,header = T,row.names = 1)
})

names(soc) = names(vir) = samples$sid
# make sure barcode order is the same
for(i in seq_len(nrow(samples))){
  if(!is.null(vir[[i]]))
    vir[[i]] = vir[[i]][rownames(soc[[i]]),]
  soc[[i]]$donor = soc[[i]]$assignment
  f = soc[[i]]$status!='singlet'
  soc[[i]]$donor[f] = soc[[i]]$status[f]
}

# _comp soc vs vireo ###########
soc2vir = lapply(samples$sid,function(s){
  if(is.null(vir[[s]]))
    return(NULL)
  table(soc[[s]]$donor,vir[[s]]$donor_id)
})
names(soc2vir) = samples$sid

pdf('soc2vir_hm.pdf',w=4*4,h=3*4)
par(mfrow=c(4,4),mar=c(6,5,1,1),oma=c(0,0,1,0),bty='n')
for(i in seq_len(nrow(samples))){
  s = samples$sid[i]
  if(!is.null(soc2vir[[s]]))
    imageWithText(log1p(soc2vir[[s]]),soc2vir[[s]],main=s,xlab='souporcell',ylab='vireo')
}
mtext(samples$tic[i],3,outer = TRUE)
dev.off()

# find best match
pairs = NULL
for(i in which(samples$n>1)){
  x = soc2vir[[i]][as.character(0:(samples$n[i]-1)),paste0('donor',0:(samples$n[i]-1))]
  p = RcppHungarian::HungarianSolver(-x)$pairs 
  p = data.frame(sample_id = samples$sid[i],
             souporcell = rownames(x)[p[,1]],
             vireo = colnames(x)[p[,2]])
  p$intersect2union = NA
  
  x = soc2vir[[i]]
  for(j in seq_len(nrow(p)))
    p$intersect2union[j] = x[p$souporcell[j],p$vireo[j]] / (sum(x[p$souporcell[j],]) + sum(x[,p$vireo[j]]) - x[p$souporcell[j],p$vireo[j]]) 
  pairs = rbind(pairs,p)
}

write.csv(pairs,'soc2vir_best_match.csv')

