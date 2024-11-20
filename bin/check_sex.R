args = commandArgs(trailingOnly = TRUE)
samples = read.table(args[1])[,1]

dir.create('sex',recursive = TRUE,showWarnings = FALSE)

# check it carefully 
gex = lapply(samples,function(s)Seurat::Read10X(paste0('gex/',s)))
soc = lapply(samples,function(s)read.table(paste0('souporcell/',s,'/clusters.tsv'),header = T,row.names = 1))

# vireo cannot berun with k=1, so lets put NULL on these places
vir = lapply(samples,function(s){
  fn = paste0('vireo/',s,'-vireo/donor_ids.tsv')
  if(!file.exists(fn))
    return(NULL)
  read.table(fn,header = T,row.names = 1)
})

names(gex) = names(soc) = names(vir) = samples


# make sure barcode order is the same
for(i in 1:length(samples)){
  if(!is.null(vir[[i]]))
    vir[[i]] = vir[[i]][colnames(gex[[i]]),]
  soc[[i]] = soc[[i]][colnames(gex[[i]]),]
  soc[[i]]$donor = soc[[i]]$assignment
  f = soc[[i]]$status!='singlet'
  soc[[i]]$donor[f] = soc[[i]]$status[f]
}

# _looks on gender ############
# sex_gids = c('F'='XIST','M'='RPS4Y1')
calSexGeneExpression = function(gex,groups,thr=2){
  result = NULL
  for(s in names(gex)){
    if(is.null(groups[[s]])) next
    
    
    pb = as.matrix(visutils::calcColSums(gex[[s]],groups[[s]]$donor))
    pb = sweep(pb,2,colSums(pb),'/')*1e4
    
    nc = as.matrix(visutils::calcColSums(gex[[s]]>0,groups[[s]]$donor))
    nc = sweep(nc,2,colSums(nc),'/')
    
    donors =  setdiff(colnames(pb),c('doublet','unassigned'))
    for(d in donors){
      t = data.frame(sample = s,
                     genotype = d,
                     XIST_cpm = pb['XIST',d],
                     XIST_fcells = nc['XIST',d],
                     RPS4Y1_cpm = pb['RPS4Y1',d],
                     RPS4Y1_fcells = nc['RPS4Y1',d])
      result = rbind(result,t)
    }
  } 
  result$sex = 'undef'
  result$sex[result$XIST_cpm     > result$RPS4Y1_cpm*thr] = 'F'
  result$sex[result$XIST_cpm*thr < result$RPS4Y1_cpm] = 'M'
  result
}
thr=2
soc_sex = calSexGeneExpression(gex,soc,thr=thr)
vir_sex = calSexGeneExpression(gex,vir,thr=thr)

plotSexGeneExp = function(d,thr,cols=c('undef'='gray','F'='red','M'='blue'),...){
  pc = c(d$XIST_cpm,d$RPS4Y1_cpm)
  if(any(pc>0)){
    pc = min(pc[pc>0])/2
  }else{
    pc = 1
  }
  
  x = log10(ifelse(d$XIST_cpm>0,d$XIST_cpm,pc))
  r = log10(ifelse(d$RPS4Y1_cpm>0,d$RPS4Y1_cpm,pc))
  lim = range(x,r)
  plot(x,r,xlim = lim,ylim=lim,xlab='XIST l10cpm (Female)',ylab='RPS4Y1 l10cpm (Male)',pch=19,col=cols[d$sex],...)
  abline(a=0,b=1)
  abline(a=-log10(thr),b=1,lty=2)
  abline(a= log10(thr),b=1,lty=2)
}

plotSexGeneExpNCell = function(d,cols=c('undef'='gray','F'='red','M'='blue'),...){
  x = d$XIST_fcells*100
  r = d$RPS4Y1_fcells*100
  lim = range(x,r)
  plot(x,r,xlim = lim,ylim=lim,xlab='XIST %cells (Female)',ylab='RPS4Y1 %cells (Male)',pch=19,col=cols[d$sex],...)
  abline(a=0,b=1)
}


# save results ########
write.csv(soc_sex,'sex/souporcell_sex_assignment.csv')
write.csv(vir_sex,'sex/vireo_sex_assignment.csv')

pdf('sex/donor_XIST_RPS4Y1_expression.pdf',height = 7,width = 7)
par(mfrow=c(2,2),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
plotSexGeneExp(soc_sex,thr=thr,main='Souporcell')
plotSexGeneExpNCell(soc_sex,main='Souporcell')

plotSexGeneExp(vir_sex,thr=thr,main='Vireo')
plotSexGeneExpNCell(vir_sex,main='Vireo')
dev.off()

gids = c('F'='XIST','M'='RPS4Y1') 

nrow = floor(sqrt(length(samples)))
ncol = ceiling(length(samples)/nrow)

dcol = visutils::char2col(unlist(lapply(soc,function(x)x$donor)))
png('sex/XIST_RPS4Y1_in_cells_souporcell.png',units = 'in',res = 400,height = nrow*3,width = ncol*3)
par(mfrow=c(nrow,ncol),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
for(sid in names(soc)){
  exp = sweep(t(as.matrix(gex[[sid]][gids,])),1,Matrix::colSums(gex[[sid]]),'/')*1e4
  plot(exp+1+runif(length(exp),-0.2,0.2),col=dcol[soc[[sid]]$donor],log='xy',pch=16,cex=0.3,main=sid,xlab='CP10K(XIST) + 1',ylab='CP10K(RPS4Y1) + 1')
  f = names(dcol) %in% soc[[sid]]$donor
  legend('topright',col=dcol[f],legend=names(dcol)[f],pch=16)
}
dev.off()


dcol = visutils::char2col(unlist(lapply(vir,function(x)x$donor)))
png('sex/XIST_RPS4Y1_in_cells_vireo.png',units = 'in',res = 400,nrow*3,width = ncol*3)
par(mfrow=c(nrow,ncol),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1))
for(sid in names(vir)){
  if(is.null(vir[[sid]])){
    plot(1,t='n',main=sid)
    next
  }
  exp = sweep(t(as.matrix(gex[[sid]][gids,])),1,Matrix::colSums(gex[[sid]]),'/')*1e4
  plot(exp+1+runif(length(exp),-0.2,0.2),col=dcol[vir[[sid]]$donor],log='xy',pch=16,cex=0.3,main=sid,xlab='CP10K(XIST) + 1',ylab='CP10K(RPS4Y1) + 1')
  f = names(dcol) %in% vir[[sid]]$donor
  legend('topright',col=dcol[f],legend=names(dcol)[f],pch=16)
}
dev.off()


