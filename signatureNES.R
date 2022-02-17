#!/usr/bin/Rscript

library(Seurat)
#
# NJ01, multifocal
rd001 = readRDS('RD-20180817-001-SR18271.seurat_final.rds')
rd002 = readRDS('RD-20180817-002-SR18271.seurat_final.rds')
rd001 = UpdateSeuratObject(rd001) 
rd002 = UpdateSeuratObject(rd002) 
rd0012.anchors <- FindIntegrationAnchors(object.list = c(rd001,rd002), dims = 1:20)
rd0012 <- IntegrateData(anchorset=rd0012.anchors,dims=1:20)
rd001t = subset(rd0012,cells=which(rd0012$putativeTumor2==1 & rd0012$orig.ident=='RD-20180817-001-SR18271'))
rd002t = subset(rd0012,cells=which(rd0012$putativeTumor2==1 & rd0012$orig.ident=='RD-20180817-002-SR18271' & rd0012$res.2!=22 & rd0012$res.2!=20)) #  & 
rd0012$res.2!=20
mat001t = as.matrix(rd001t[['RNA']]@data)
mat002t = as.matrix(rd002t[['RNA']]@data)
rd001t.f = subset(rd001t,cells=which(unlist(mat001t['PTPRC',])==0))
rd002t.f = subset(rd002t,cells=which(unlist(mat002t['PTPRC',])==0))
mat001t.f = as.matrix(rd001t.f[['RNA']]@data)
mat002t.f = as.matrix(rd002t.f[['RNA']]@data)
exclude.idx = unique(c(
  which(apply(mat001t.f,1,sum)==0),
  which(apply(mat002t.f,1,sum)==0)
))
mat001t.ff = mat001t.f[-exclude.idx,]
mat002t.ff = mat002t.f[-exclude.idx,]
stat001 = as.data.frame(t(apply(mat001t.ff,1,function(x){
  c(
    length(which(x>0))/length(x),
    sd(x)/mean(x)
  )
})))
colnames(stat001) = c('ratio','cv')
stat002 = as.data.frame(t(apply(mat002t.ff,1,function(x){
  c(
    length(which(x>0))/length(x),
    sd(x)/mean(x)
  )
})))
colnames(stat002) = c('ratio','cv')
stat012 = as.data.frame(t(apply(cbind(mat001t.ff,mat002t.ff),1,function(x,i){
  c(
    mean(x[i]),mean(x[-i]),
    mean(x[i])/mean(x[-i]),
    mean(x[-i])/mean(x[i]),
    wilcox.test(x[i],x[-i])$p.value
  )
},i=1:ncol(mat001t.ff))))
colnames(stat012) = c('case','ctrl','fc1','fc2','p')
#
statsigup = rownames(stat001)[which(stat001$ratio-stat002$ratio>0.05)]
statsigdn = rownames(stat001)[which(stat002$ratio-stat001$ratio>0.05)]
pool001o = apply(mat001t.ff[statsigup,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=rd001t.f$res.2)
pool002o = apply(mat002t.ff[statsigup,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=rd002t.f$res.2)
pool001y = apply(mat001t.ff[statsigdn,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=rd001t.f$res.2)
pool002y = apply(mat002t.ff[statsigdn,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=rd002t.f$res.2)
pool001ostat = apply(pool001o,1,function(x){length(which(x>0.1))/length(x)})
pool002ostat = apply(pool002o,1,function(x){length(which(x>0.1))/length(x)})
pool001ystat = apply(pool001y,1,function(x){length(which(x>0.1))/length(x)})
pool002ystat = apply(pool002y,1,function(x){length(which(x>0.1))/length(x)})
pool001stat = pool001ostat/pool001ystat
pool002stat = pool002ostat/pool002ystat
lesion1old = tail(sort(pool001stat),round(length(pool001stat)/5))
lesion2old = tail(sort(pool002stat),round(length(pool002stat)/5))
lesion1yng = head(sort(pool001stat),round(length(pool001stat)/5))
lesion2yng = head(sort(pool002stat),round(length(pool002stat)/5))
pool1 = as.data.frame(t(apply(mat001t.ff,1,function(x,i,j){
  c(
    mean(x[i]),mean(x[j]),mean(x[i])/mean(x[j]),
    length(which(x[i]>0))/length(i),
    length(which(x[j]>0))/length(j),
    wilcox.test(x[i],x[j])$p.value
  )
},i=which(rd001t.f$res.2%in%names(lesion1old)),
j=which(rd001t.f$res.2%in%names(lesion1yng)))))
colnames(pool1) = c('case','ctrl','fc','r1','r2','p')
pool1$q = p.adjust(pool1$p)
pool2 = as.data.frame(t(apply(mat002t.ff,1,function(x,i,j){
  c(
    mean(x[i]),mean(x[j]),mean(x[i])/mean(x[j]),
    length(which(x[i]>0))/length(i),
    length(which(x[j]>0))/length(j),
    wilcox.test(x[i],x[j])$p.value
  )
},i=which(rd002t.f$res.2%in%names(lesion2old)),
j=which(rd002t.f$res.2%in%names(lesion2yng)))))
colnames(pool2) = c('case','ctrl','fc','r1','r2','p')
pool2$q = p.adjust(pool2$p)
#
pool3 = as.data.frame(t(apply(
  cbind(mat001t.ff[,which(rd001t.f$res.2%in%names(lesion1old))],
        mat002t.ff[,which(rd002t.f$res.2%in%names(lesion2yng))]),1,function(x,i){
          c(
            mean(x[i]),mean(x[-i]),mean(x[i])/mean(x[-i]),
            length(which(x[i]>0))/length(i),
            length(which(x[-i]>0))/(length(x)-length(i)),
            wilcox.test(x[i],x[-i])$p.value
          )
        },i=1:length(which(rd001t.f$res.2%in%names(lesion1old)))
)))
colnames(pool3) = c('case','ctrl','fc','r1','r2','p')
pool3$q = p.adjust(pool3$p)
#
# TT001, multifocal 
dat673   = readRDS('/public/workspace/lily/PS/P673_jc.RDS')
dat689   = readRDS('/public/workspace/lily/PS/P689.RDS')
dat673.t = subset(dat673,cells=which(dat673$seurat_clusters %in% c(3,5)))
dat689.t = subset(dat689,cells=which(dat689$seurat_clusters %in% c(1,3,4,6)))
mat673.t = as.matrix(dat673.t[['RNA']]@data)
mat689.t = as.matrix(dat689.t[['RNA']]@data)
dat673.tf = subset(dat673.t,cells=which(unlist(mat673.t['PTPRC',])==0))
dat689.tf = subset(dat689.t,cells=which(unlist(mat689.t['PTPRC',])==0))
dat673.tf$seurat_clusters = as.factor(as.vector(dat673.tf$seurat_clusters))
dat689.tf$seurat_clusters = as.factor(as.vector(dat689.tf$seurat_clusters))
mat673.tf = mat673.t[,which(unlist(mat673.t['PTPRC',])==0)]
mat689.tf = mat689.t[,which(unlist(mat689.t['PTPRC',])==0)]
rm(mat673.t);rm(mat689.t)
cm.gene = intersect(rownames(mat673.tf),rownames(mat689.tf))
pool6 = as.data.frame(t(apply(cbind(mat673.tf[cm.gene,],mat689.tf[cm.gene,]),1,function(x,i){
  c(
    mean(x[i]),mean(x[-i]),mean(x[i])/mean(x[-i]),
    length(which(x[i]>0))/length(i),
    length(which(x[-i]>0))/(length(x)-length(i)),
    wilcox.test(x[i],x[-i])$p.value
  )
},i=1:ncol(mat673.tf))))
colnames(pool6) = c('case','ctrl','fc','r1','r2','p')
pool6$q = p.adjust(pool6$p)
#
# TT002, multifocal 
dat912L = readRDS('/public/workspace/lily/PS/P912.L.RDS')
dat912R = readRDS('/public/workspace/lily/PS/P912.R.RDS')
dat912L.t = subset(dat912L,cells=which(dat912L$seurat_clusters %in% c(1,3,4,8,9)))
dat912R.t = subset(dat912R,cells=which(dat912R$seurat_clusters %in% c(3,4,5)))
mat912L.t = as.matrix(dat912L.t[['RNA']]@data)
mat912R.t = as.matrix(dat912R.t[['RNA']]@data)
dat912L.tf = subset(dat912L.t,cells=which(unlist(mat912L.t['PTPRC',])==0))
dat912R.tf = subset(dat912R.t,cells=which(unlist(mat912R.t['PTPRC',])==0))
dat912L.tf$seurat_clusters = as.factor(as.vector(dat912L.tf$seurat_clusters))
dat912R.tf$seurat_clusters = as.factor(as.vector(dat912R.tf$seurat_clusters))
mat912L.tf = mat912L.t[,which(unlist(mat912L.t['PTPRC',])==0)]
mat912R.tf = mat912R.t[,which(unlist(mat912R.t['PTPRC',])==0)]
rm(mat912L.t);rm(mat912R.t)
cm.gene = intersect(rownames(mat912L.tf),rownames(mat912R.tf))
pool7 = as.data.frame(t(apply(cbind(mat912L.tf[cm.gene,],mat912R.tf[cm.gene,]),1,function(x,i){
  c(
    mean(x[i]),mean(x[-i]),mean(x[i])/mean(x[-i]),
    length(which(x[i]>0))/length(i),
    length(which(x[-i]>0))/(length(x)-length(i)),
    wilcox.test(x[i],x[-i])$p.value
  )
},i=1:ncol(mat912L.tf))))
colnames(pool7) = c('case','ctrl','fc','r1','r2','p')
pool7$q = p.adjust(pool7$p)
#
# NJ02, multifocal
library(Seurat) 
h3 = readRDS('10XT1H3.RDS')
h4 = readRDS('10XT1H4.RDS')
mat3 = as.matrix(h3[['RNA']]@data)
mat4 = as.matrix(h4[['RNA']]@data)
h3$celltype = NA
h3$celltype[which(h3$seurat_clusters %in% c(0,1,2,3,4,5,7,8,9,11,12,14,15,20))] = 'malignant'
h3$celltype[which(h3$seurat_clusters %in% c(10,13,16,19))] = 'myeloid'
h3$celltype[which(h3$seurat_clusters %in% c(6,17))] = 'oligodendrocyte'
h3$celltype[which(h3$seurat_clusters %in% c(18))] = 'fibroblast/vascular'
h3$celltype[which(h3$seurat_clusters %in% c(21))] = 'T cells'
h4$celltype = NA
h4$celltype[which(h4$seurat_clusters %in% c(0,1,2,3,7,8,10,13,16,20,21))] = 'malignant'
h4$celltype[which(h4$seurat_clusters %in% c(4,6,9,11,12,14,18,19))] = 'myeloid'
h4$celltype[which(h4$seurat_clusters %in% c(5,17))] = 'oligodendrocyte'
h4$celltype[which(h4$seurat_clusters %in% c(15))] = 'fibroblast/vascular'
h4$celltype[which(h4$seurat_clusters %in% c(22))] = 'T cells'
h4$celltype[which(h4$seurat_clusters %in% c(23))] = 'Endothelial'
h4$celltype[which(h4$seurat_clusters %in% c(24))] = 'unknown'
h3$orig.ident = 'h3'
h4$orig.ident = 'h4'
h34.anchor = FindIntegrationAnchors(object.list = c(h3,h4), dims = 1:20)
h34 <- IntegrateData(anchorset=h34.anchor,dims=1:20)
h3t = subset(h34,cells=which(h34$celltype=='malignant' & h34$orig.ident=='h3'))
h4t = subset(h34,cells=which(h34$celltype=='malignant' & h34$orig.ident=='h4'))
mat3t = as.matrix(h3t[['RNA']]@data)
mat4t = as.matrix(h4t[['RNA']]@data)
h3t.f = subset(h3t,cells=which(unlist(mat3t['PTPRC',])==0))
h4t.f = subset(h4t,cells=which(unlist(mat4t['PTPRC',])==0))
mat3t.f = as.matrix(h3t.f[['RNA']]@data)
mat4t.f = as.matrix(h4t.f[['RNA']]@data)
exclude.idx = unique(c(
  which(apply(mat3t.f,1,sum)==0),
  which(apply(mat4t.f,1,sum)==0)
))
mat3t.ff = mat3t.f[-exclude.idx,]
mat4t.ff = mat4t.f[-exclude.idx,]
stat3 = as.data.frame(t(apply(mat3t.ff,1,function(x){
  c(
    length(which(x>0))/length(x),
    sd(x)/mean(x)
  )
})))
colnames(stat3) = c('ratio','cv')
stat4 = as.data.frame(t(apply(mat4t.ff,1,function(x){
  c(
    length(which(x>0))/length(x),
    sd(x)/mean(x)
  )
})))
colnames(stat4) = c('ratio','cv')
stat34 = as.data.frame(t(apply(cbind(mat3t.ff,mat4t.ff),1,function(x,i){
  c(
    mean(x[i]),mean(x[-i]),
    mean(x[i])/mean(x[-i]),
    mean(x[-i])/mean(x[i]),
    wilcox.test(x[i],x[-i])$p.value
  )
},i=1:ncol(mat3t.ff))))
colnames(stat34) = c('case','ctrl','fc1','fc2','p')
#
statsigup = rownames(stat3)[which(stat4$ratio-stat3$ratio>0.05)]
statsigdn = rownames(stat3)[which(stat3$ratio-stat4$ratio>0.05)]
pool4o = apply(mat4t.ff[statsigup,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=h4t.f$seurat_clusters)
pool3o = apply(mat3t.ff[statsigup,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=h3t.f$seurat_clusters)
pool4y = apply(mat4t.ff[statsigdn,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=h4t.f$seurat_clusters)
pool3y = apply(mat3t.ff[statsigdn,],1,function(x,y){
  unlist(tapply(x,y,function(x){length(which(x>0))/length(x)}))
},y=h3t.f$seurat_clusters)
pool4ostat = apply(pool4o,1,function(x){length(which(x>0.1))/length(x)})
pool3ostat = apply(pool3o,1,function(x){length(which(x>0.1))/length(x)})
pool4ystat = apply(pool4y,1,function(x){length(which(x>0.1))/length(x)})
pool3ystat = apply(pool3y,1,function(x){length(which(x>0.1))/length(x)})
pool4stat = pool4ostat/pool4ystat
pool3stat = pool3ostat/pool3ystat
lesion4old = tail(sort(pool4stat),round(length(pool4stat)/5))
lesion3old = tail(sort(pool3stat),round(length(pool3stat)/5))
lesion4yng = head(sort(pool4stat),round(length(pool4stat)/5))
lesion3yng = head(sort(pool3stat),round(length(pool3stat)/5))
#
pool04 = as.data.frame(t(apply(mat4t.ff,1,function(x,i,j){
  c(
    mean(x[i]),mean(x[j]),mean(x[i])/mean(x[j]),
    length(which(x[i]>0))/length(i),
    length(which(x[j]>0))/length(j),
    wilcox.test(x[i],x[j])$p.value
  )
},i=which(h4t.f$seurat_clusters%in%names(lesion4old)),
j=which(h4t.f$seurat_clusters%in%names(lesion4yng)))))
colnames(pool04) = c('case','ctrl','fc','r1','r2','p')
pool04$q = p.adjust(pool04$p)
pool03 = as.data.frame(t(apply(mat3t.ff,1,function(x,i,j){
  c(
    mean(x[i]),mean(x[j]),mean(x[i])/mean(x[j]),
    length(which(x[i]>0))/length(i),
    length(which(x[j]>0))/length(j),
    wilcox.test(x[i],x[j])$p.value
  )
},i=which(h3t.f$seurat_clusters%in%names(lesion3old)),
j=which(h3t.f$seurat_clusters%in%names(lesion3yng)))))
colnames(pool03) = c('case','ctrl','fc','r1','r2','p')
pool03$q = p.adjust(pool03$p)
#
pool34 = as.data.frame(t(apply(
  cbind(mat4t.ff[,which(h4t.f$seurat_clusters%in%names(lesion4old))],
        mat3t.ff[,which(h3t.f$seurat_clusters%in%names(lesion3yng))]),1,function(x,i){
          c(
            mean(x[i]),mean(x[-i]),mean(x[i])/mean(x[-i]),
            length(which(x[i]>0))/length(i),
            length(which(x[-i]>0))/(length(x)-length(i)),
            wilcox.test(x[i],x[-i])$p.value
          )
        },i=1:length(which(h4t.f$seurat_clusters%in%names(lesion4old)))
)))
colnames(pool34) = c('case','ctrl','fc','r1','r2','p')
pool34$q = p.adjust(pool34$p)
#
final.sig = c()
cutoff.fc=1.3 
cutoff.rd=0.05
sig34 = pool34[which(pool34$fc>cutoff.fc & pool34$r1>cutoff.rd & pool34$q<0.001),]
sigs3 = pool3[which(pool3$fc>cutoff.fc & pool3$r1>cutoff.rd & pool3$q<0.001),]
sigs6 = pool6[which(pool6$fc>cutoff.fc & pool6$r1>cutoff.rd & pool6$q<0.001),]
sigs7 = pool7[which(pool7$fc>cutoff.fc & pool7$r1>cutoff.rd & pool7$q<0.001),]
final.sig = intersect(rownames(sig34),rownames(sigs3))
final.sig = intersect(final.sig,rownames(sigs6))
final.sig = intersect(final.sig,rownames(sigs7)) 
length(final.sig)
intrinsicGene = read.table('intrinsicGliomaGene.txt',header=T,sep='\t') # Wang et al., Cancer cell, 2017
final.sig = intersect(final.sig,intrinsicGene$Gene_Symbol)
mod.generate(final.sig,'NES',out='NES.mod')






