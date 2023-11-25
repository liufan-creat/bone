
library(dplyr)
library(tidyr)
library(Seurat)
library(ggplot2)
library(Rcpp)
library(harmony)
library(clusterProfiler) ## rvcheck = 1.7 不能升级
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(readxl)
library(hdf5r)

library(GSEABase)
library(GSVA)


exp <- read_excel('~/GEO/bone/Conclusion.xlsx',sheet = 1)
exp <- exp[,c('gene_id',"log2FC...15","P_value...17")]
colnames(exp ) <- c('id','log2FC','P_value')

dir='~/GEO/bone/GSE199866_RAW/' 
samples=list.files( dir ) 
samples
# data <-data.table::fread(file.path(dir,samples[1]  ),sep = "\t",header = T,data.table = F)

dir='~/GEO/bone/GSE211407_RAW/' 
samples=list.files( dir ) 
samples
 
group <- c('C1','C2','C3','C4','C5','C6','C7','C8')


sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro) 
  sce <-  read.table(file.path(dir,pro ),sep = '\t',header = T)
  rownames(sce) <- paste0(sce[,2],"_",sce[,1])
  sce <- sce[,-c(1:9)]
  sce =CreateSeuratObject(counts = t( sce ) ,
                          project =  gsub('^GSM[0-9]*_','',pro)  ,
                          min.cells = 5,
                          min.features = 500 )
  return(sce)
})


sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)   
  sce =CreateSeuratObject(counts =  Read10X_h5(file.path(dir,pro )),
                          project =  group[which(samples== pro)],
                          min.cells = 5,
                          min.features = 500 )
  return(sce)
})
names(sceList) 

sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids =  group    )


sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-") ## human改^MT-
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sce.all <- subset(sce.all, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sce.all <- NormalizeData(sce.all,verbose =F)
sce.all <- FindVariableFeatures(sce.all, 
                                selection.method = "vst", nfeatures = 2000) 
sce.all <- ScaleData(sce.all,verbose = F) 
sce.all <- RunPCA(object = sce.all, npcs = 50, verbose = FALSE) 


seuratObj <- RunHarmony(sce.all,group.by.vars = 'orig.ident') ### ???
#
# <- RunTSNE(sce.all,dims = 1:7)

seuratObj <- RunUMAP(seuratObj,reduction = "harmony", dims = 1:20,
                     min.dist = 0.7, n.neighbors = 4L)

seuratObj <- FindNeighbors(seuratObj, reduction = "harmony", dims = 1:20)
seuratObj <- FindClusters(seuratObj, resolution = 0.6)

# 不矫正
# seuratObj <- RunTSNE(pbmc,  dims = 1:10)
# seuratObj <- FindNeighbors(seuratObj,  dims = 1:7)
# seuratObj <- FindClusters(seuratObj, resolution = 0.3)

cols <-  c('#F16300','#FF007F','#00A2DA','#68B195','#ECD12E','#00667A',
           '#C0693D','#314F56','#bcd7ec','#31A567','grey','red','green')
mycol <- c('#F16300','#ECD12E','#FF007F','#00A2DA','#72ec89','#68B195',
           '#C0693D','#0053DA','blue')

DimPlot(seuratObj,reduction = "umap",label=T, 
         split.by= "orig.ident",ncol=2,  
        label.size = 8,    
        pt.size = 1) + NoLegend() 

meta <- seuratObj@meta.data




features <- c('Comp','Acan','Fmod','Sox9','Myb','Cdh5','Ptprc','Blnk','Foxo1',
              'Cd80','Akt1s1','Strn','Pik3ca','Erbb2','Fgf9','Fgf10','Pten','Tsc2',
              'Pik3r2','Nfkb1','Raf1','Gab1')

features <- c('Scrg1','Plp1','Ptn','Gpm6b','Sox4',
              'S100b', 'Ptprz1','Prcp',
              'Cd80','Akt1s1','Strn','Pik3ca','Erbb2','Fgf9','Fgf10','Pten','Tsc2',
              'Pik3r2','Nfkb1','Raf1','Gab1','Fmod','Sox9','Myb','Cdh5','Ptprc','Blnk','Foxo1')

features <- c('Myb','Cdh5','Ptprc','Blnk','Plp1','Ptn','Gpm6b','Gd2','Tie2',
              'Uts2r','Hoxa9','Hoxd9','Etv5','Foxp2','Mecom','Nfatc4','Nr2f2','Sox4',
              'Frzb', 'Dkk','Col1a1','Prg4','Krt19','Tbxt','Pdgfra',
              'Prrx','Igf1','Tbxt','Krt8')


features <- c('Ell3','Zfp260','Hdhd5','Zscan26','Epc2','Fis1','Ltbp3',
              'Ece2','Wdfy1','Itfg1','Cyp2j10',
              'Rnase4','Rbm5','Neo1','Arpp19') 

features <- c( 'CCNB1','CCNB2','CCNA2','CDKN1C','CDKN1B','CDKN1A',
               'CD34A，PROM1','THY1','NANOG','OCT4','SOX2','SOX9')

feature <- c('CCND3', 'PDK1','SMARCA2','FOXO3', 'EZH1', 'PRDM5', 'PTOV1', 'ZFP30', 
             'ZBTB20', 'PHF1', 'CTDSP1', 'THRA', 'TEF','DICER1','A930001N09Rik', 'BCAS3',
             'DDX3Y', 'GABARAPL1', 'GLTSCR2', 'ITM2A', 'IL18', 'ZYX', 'EPHX1', 'CLSTN1',
             'GSTK1', '5730403B10Rik', 'DDT', 'IVD', 'FHL1', 'NDRG2', 'GRINA', 
             'PIK3R1', 'FYN', 'CHKB', 'PINK1', 'ULK2', 'DNAJB9', 'PFDN5', 'CTSF', 
             'CRIM1', 'SEPP1', 'GABBR1', 'GRB10', 'BBS2', 'RPS14', 'IGF2R',
             'SELENBP1', 'RNF167', 'MAP1LC3A')

down <- c('ANLN', 'BIRC5', 'CCNA2', 'CCNB1', 'CCNE2', 'SGOL1','MCM4', 'PCNA', 'RRM2',
          'TOP2A','CYCS', 'MTCH2', 'SLC25A5','H2AFZ', 'HAT1','DDX39','2810417H13Rik',
          'CAPZA1', 'HADHB', 'IDH3A', 'KPNA2', 'PGK1')

VlnPlot(seuratObj,stack = T,   
        features = s.genes)+ NoLegend() 

mark <- read_xlsx('Conclusion.xlsx',sheet = 1)
mark <- mark[rev(order(mark$Fold_change) ),]
features <- mark$gene_symbol[1:200]

s.genes <- CaseMatch(search = down,match = rownames(seuratObj))

match = as.data.frame( colnames(seuratObj) )


all <- FindAllMarkers(seuratObj)
all <- all[which(all$p_val_adj < 0.05),]

sce.all <- RunPCA(object = seuratObj , features = features,
                  npcs = 3, verbose = TRUE) 

DimPlot(sce.all,reduction = 'pca')

ps <-DimPlot(seuratObj,reduction = "tsne",label=T,  
             label.size = 6, cols= c(png,png3) ,
             # split.by= "orig.ident",ncol=4,  
             pt.size = 1) + NoLegend() +
  theme(axis.text.x = element_text(size=27 ), ### 刻度大小
        axis.title.y = element_text(size=27 ),
        axis.text.y = element_text(size=27), ### 刻度大小
        axis.title.x = element_text(size=27),
        legend.title =element_text(size=27,hjust = 0.5),
        legend.text  = element_text(size=27 ))

pl <- VlnPlot(seuratObj,stack = T,  cols= c(png,png3) , 
              features = features  )+ NoLegend()  + 
  theme(axis.ticks.x  = element_blank(), 
        axis.text.x = element_text(size=30,  angle = 270), 
        axis.title.x =  element_blank(),  
        axis.text.y = element_text(size=30,  angle = 0), 
        axis.title.y = element_blank()) 

pl$theme$strip.text  <- pl$theme$axis.text.y 

pl$theme$text$hjust <- 1
pl$theme$text$size <- 26
pl



### cycle
s.genes <- cc.genes$s.genes
s.genes <- CaseMatch(search = s.genes,match = rownames(seuratObj))

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- CaseMatch(search = g2m.genes,match = rownames(seuratObj))

seuratObj <- CellCycleScoring(seuratObj, s.features = s.genes,
                              g2m.features = g2m.genes  ,set.ident = TRUE)

meta <- meta[which(meta$Phase =='G1'),]
a <- seuratObj@assays$RNA@counts

a <- a[,rownames(meta)]

G1 <- CreateSeuratObject(counts = a )

G1 <- NormalizeData(G1,verbose =F)
G1 <- FindVariableFeatures(G1,selection.method = "vst" ) 
G1 <- ScaleData(G1,verbose = F) 
G1 <- RunPCA(object = G1, npcs = 50, verbose = FALSE) 


seuratObj <- RunHarmony(G1,group.by.vars = 'orig.ident') ### ???
#seuratObj <- RunTSNE(seuratObj, reduction = "harmony",dims = 1:7)

seuratObj <- RunUMAP(seuratObj,reduction = "harmony", dims = 1:10,
                     min.dist = 0.7, n.neighbors = 10L)

seuratObj <- FindNeighbors(seuratObj,  dims = 1:3)
seuratObj <- FindClusters(seuratObj, resolution = 0.1)

DimPlot(seuratObj  )



s.genes <- cc.genes$s.genes
s.genes <- CaseMatch(search = s.genes,match = rownames(seuratObj))

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- CaseMatch(search = g2m.genes,match = rownames(seuratObj))

seuratObj <- CellCycleScoring(seuratObj, s.features = s.genes,
                              g2m.features = g2m.genes  ,set.ident = TRUE)

seuratObj <- RunPCA(seuratObj,features = c(s.genes,g2m.genes))

DimPlot(seuratObj,reduction = 'umap' )


DimPlot(seuratObj,reduction = "umap",label=T, 
        split.by= "orig.ident",ncol=2,  
        label.size = 8,    
        pt.size = 1) + NoLegend() 


### scrna######

hs.counts <- read.csv('GSE64016_H1andFUCCI_normalized_EC.csv.gz', header=TRUE, row.names=1)

hs.G1 <- grepl("G1", colnames(hs.counts))
hs.S <- grepl("S", colnames(hs.counts))
hs.G2 <- grepl("G2", colnames(hs.counts))

counts <- cbind(hs.counts[,hs.G1],hs.counts[,hs.S],hs.counts[,hs.G2])


pbmc <- CreateSeuratObject(counts =counts,
                           project = "heart", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 4000)

pbmc <- ScaleData(pbmc )
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimPlot(pbmc)

mark <- FindAllMarkers(pbmc)
mark <- mark[which(mark$p_val_adj < 0.05 &abs(  mark$avg_log2FC) >0.5 ),]


s.genes <- mark$gene[which(mark$cluster =='S')]
s.genes <- CaseMatch(search = s.genes,match = rownames(seuratObj))

g2m.genes <- mark$gene[which(mark$cluster =='G2')]
g2m.genes <- CaseMatch(search = g2m.genes,match = rownames(seuratObj))


g1m.genes <- mark$gene[which(mark$cluster =='G1')]
g1m.genes <- CaseMatch(search = g1m.genes,match = rownames(seuratObj))

seuratObj <- CellCycleScoring(seuratObj, s.features = s.genes,
                              g2m.features = g2m.genes  ,set.ident = TRUE)

seuratObj <- RunPCA(seuratObj,features = c(s.genes, g2m.genes,g1m.genes))

DimPlot(seuratObj,reduction = 'pca')
RidgePlot(seuratObj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

## 
seuratObj$CC.Difference <- seuratObj$S.Score - seuratObj$G2M.Score 
seuratObj <- ScaleData(seuratObj, vars.to.regress = "CC.Difference")


seuratObj <- RunPCA(object = seuratObj,  features = c(s.genes, g2m.genes),
                    npcs = 40, verbose = TRUE) 

DimPlot(seuratObj,reduction = 'pca',ncol=2 )



RidgePlot(sce.all, features = c("Ccnd1", "Ccne2", "Ccna1", "Ccnb1"), ncol = 2)

####  all cc gene ##### 

m_df <- read.gmt('TFEB_TARGET_GENES.v2023.1.Hs.gmt')

b <-  as.data.frame(  unique(m_df$term) )
b$value <- 0

for (i in 1:length(b[,1])) {
  a <- strsplit(as.character(b[i,1]),split = '_')[[1]]
  if('CELL' %in% a & 'CYCLE' %in% a){b$value[i] <- 1}
}
b<- b[which(b$value == '1'),]

m_df <- m_df[which(m_df$term %in% c('GOBP_G0_TO_G1_TRANSITION',
                                    'GOBP_G1_TO_G0_TRANSITION',
                                    'REACTOME_G0_AND_EARLY_G1')) ,]

g0.genes <- CaseMatch(search = m_df$gene,match = rownames(sce.all))


cycle.genes <- CaseMatch(search = m_df$gene[which(m_df$term %in% b$`unique(m_df$term)`)],
                         match = rownames(sce.all))

a <- a[which(rownames(a) %in% cycle.genes),]

pbmc <- CreateSeuratObject(counts = a )

pbmc <- NormalizeData(pbmc,verbose =F)
pbmc <- FindVariableFeatures(pbmc,selection.method = "vst" ) 
pbmc <- ScaleData(pbmc,verbose = F) 
pbmc <- RunPCA(object = pbmc, npcs = 50, verbose = FALSE) 

seuratObj <- RunHarmony(pbmc,group.by.vars = 'orig.ident') ### ???
#seuratObj <- RunTSNE(seuratObj, reduction = "harmony",dims = 1:7)

seuratObj <- RunUMAP(seuratObj,reduction = "harmony", dims = 1:10,
                     min.dist = 0.7, n.neighbors = 4L)

seuratObj <- FindNeighbors(seuratObj,  dims = 1:10)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)

DimPlot(seuratObj,reduction = 'tsne',
        cols = c(rep('grey',6),'red',rep('grey',11) ))

cluster6 <- all[which(all$cluster == '6'),]
gene <- list(cluster6$gene,m_df$gene)


#veen
venn(gene,box = F,
     opacity = 0.7, plotsize = 10, ilcs = 1.6, sncs = 1.5 )

con_gene <- Reduce(intersect,gene)
cluster6 <- cluster6[which(cluster6$gene %in% con_gene),]


seuratObj <- CellCycleScoring(seuratObj, s.features = s.genes,replace = T,
                              g2m.features = g2m.genes  ,set.ident = TRUE)
DimPlot(seuratObj  )


meta <- seuratObj@meta.data
meta <- meta[which(rownames(meta) %in% colnames(a)),]

for (i in 1:length(meta[,1])) {
  id <- rownames(meta)[i]
  meta$tfeb[i] <- Es[1,id]
}

## 
sce.all <- ScaleData(sce.all,verbose = T) 
sce.all <- RunPCA(object = sce.all , features = c(s.genes, g2m.genes),
                  npcs = 40, verbose = TRUE) 

DimPlot(seuratObj  )

### GSVA #### 
s.genes <- as.data.frame(s.genes)
colnames(s.genes) <- 'gene'
s.genes$term <- 'S'
s.genes <- s.genes[,c('term','gene')]

g1m.genes <-  as.data.frame(g1m.genes)
colnames(g1m.genes) <- 'gene'
g1m.genes$term <- 'G1'
g1m.genes <- g1m.genes[,c('term','gene')]


g2m.genes <- as.data.frame(g2m.genes)

colnames(g2m.genes) <- 'gene'
g2m.genes$term <- 'G2'
g2m.genes <- g2m.genes[,c('term','gene')]

choose <- rbind(g1m.genes,s.genes,g2m.genes)

kegg2symbol_list <- tapply(choose[,2],as.factor(choose[,1]),function(x) x)

write.gmt <- function(geneSet=kegg2symbol_list,gmt_file='kegg2symbol.gmt'){
  sink( gmt_file )
  for (i in 1:length(geneSet)){
    cat(names(geneSet)[i])
    cat('\tNA\t')
    cat(paste(geneSet[[i]],collapse = '\t'))
    cat('\n')
  }
  sink()
} 

write.gmt(kegg2symbol_list,'choose.gmt')

keggSet <- GSEABase::getGmt('TFEB_TARGET_GENES.v2023.1.Hs.gmt')
Es <- gsva(expr=as.matrix(a), gset.idx.list=keggSet, 
           method='ssgsea',  ## gsva/ssgsea    zscore（log2(FPKM+1)）/plage (rnaseq 不适用？)
           ##kcdf ： Poisson--read count数据;log 后的CPM, RPKM, TPM等_Gaussian
           kcdf="Poisson", parallel.sz=4)
Es <- as.data.frame(Es)
rownames(Es) <- c('glu','fat')
Es <- as.data.frame(t(Es))
Es$id <- rownames(Es)


## classifier####

mark <- read_xlsx('Conclusion.xlsx',sheet = 1)

mark <- mark[,4:10]

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

for (i in 2:7) {
  mark[,i] <-fpkmToTpm( mark[,i] )
} 

mark[,2:7] <- mark[,2:7]/100
sce.all <- seuratObj@assays$RNA@counts

gene_len <- read.csv('/Users/liufan/rna/reference/rn6_len.csv',stringsAsFactors = F)
gene_len <- gene_len[,c(  "Associated.Gene.Name","width")]
gene_len <- aggregate(gene_len[2],by = list(gene_len$Associated.Gene.Name),
                      FUN = mean,na.rm =T)

gene <- list(gene_len$Group.1,rownames(sce.all))

con <- Reduce(intersect,gene)

sce.all <- sce.all[which(rownames(sce.all) %in% con),]
gene_len <- gene_len[which(gene_len$Group.1 %in% con) ,]


gene_len$Group.1 <- factor(gene_len$Group.1,
                                        levels = rownames(sce.all))
gene_len <- gene_len[order(gene_len$Group.1),]

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

for (i in 1:length(colnames(sce.all))) {
  sce.all[,i] <- countToTpm(sce.all[,i],gene_len$width)
  print(i)
}

train <- mark

con <- Reduce(intersect,list(train$gene_symbol,gene_len$Group.1))
train <- train[which(train$gene_symbol %in% con),]
gene_len <- gene_len[which(gene_len$Group.1 %in% con),]

train$gene_symbol <- factor(train$gene_symbol,
                            levels =gene_len$Group.1 )
train  <- train[order(train$gene_symbol),]


l <- 10367 
input <- layer_input(shape = l)

class <- input%>%
  layer_dense(units = 2048 ,input_shape = l) %>% 
  layer_activation_relu( )%>% 
  layer_dropout(rate=0.2)%>% 
  layer_dense(units = 512 ,input_shape = l) %>% 
  layer_activation_relu( )%>% 
  layer_dropout(rate=0.2)%>% 
  layer_dense(units = 10 ,input_shape = l) %>% 
  layer_activation_relu( )%>% 
  layer_dropout(rate=0.2)%>% 
  layer_dense(units = 1,activation = 'sigmoid') 

classifer <- keras_model(input,class)
classifer 


classifer  %>% compile(optimizer = optimizer_adamax(learning_rate = 0.0001),
                       loss = 'binary_crossentropy' )

batch <- 3
closs <- NULL

## 1 Qs 2 con
for (i in 1:500) {
  real <- train[,2:7]
  real <- as.matrix( t(real) ) 
  
  real_res <- rbind( array(runif(batch, 0.9 ,  1),dim =c(batch,1)),
                     array(runif(batch, 0,  0.1),dim =c(batch,1)))
  
  
  closs[i]  <- classifer%>%train_on_batch( real,real_res)
  print(i)
  print(closs[i])
  
}

real <- train[,2:7]
real <- as.matrix( t(real) ) 



fake <- predict_on_batch(classifer,real)

sceList <-  sce.all[which(rownames(sce.all) %in% con),]

meta <- seuratObj@meta.data
cell <- c()
cluster <-as.data.frame(  unique(meta$RNA_snn_res.0.6) )

for (i in 1:19) {
  id <- sample(rownames(meta)[which(meta$RNA_snn_res.0.6 == cluster$`unique(meta$RNA_snn_res.0.6)`[i])],50,replace = F)
  cell <- c(cell,id)
}

choose  <-  sceList[,which(colnames(sceList) %in% cell)]


for (i in 1:length(colnames(choose))) {
  choose[,i] <- countToTpm(choose[,i],gene_len$width)
  print(i)
}
choose <-choose[which(rownames(choose) %in% con),]

fake <- predict_on_batch(classifer,t(choose)/100)
fake <- as.data.frame(fake)
fake$cell <- colnames(choose)
meta$cell <- rownames(meta)
merge_data <- merge(fake,meta,all=F)
merge_data$V1 <- round(merge_data$V1,2)
merge_data <- merge_data[,c("cell"      ,      "V1", 'seurat_clusters')]


ggplot(merge_data,aes( x=seurat_clusters,y=V1,colour=seurat_clusters)) +
  geom_boxplot( size= 1,colour='#555555')+
  
  geom_point(aes(colour=seurat_clusters),size=2,shape=20)+ 
  theme_classic()  +    geom_hline(yintercept = 0.5,lty=1  ,lwd=1)+ 
  NoLegend() +
  theme(plot.title = element_text(size=27,hjust = 0.5),
        axis.text.x = element_text(size=27 ), ### 刻度大小
        axis.title.y = element_text(size=27 ),
        axis.text.y = element_text(size=27), ### 刻度大小
        axis.title.x = element_text(size=27),
        legend.title =element_text(size=27,hjust = 0.5),
        legend.text  = element_text(size=27 ))


save_model_weights_hdf5(classifer,"classifer.h5")


l <- 14932 
input <- layer_input(shape = l)

class <- input%>%
  layer_dense(units = 1024 ,input_shape = l) %>% 
  layer_activation_relu( )%>% 
  layer_dropout(rate=0.2)%>% 
  layer_dense(units = 1,activation = 'sigmoid') 

classifer <- keras_model(input,class)
classifer 
load_model_weights_hdf5(classifer,"classifer.h5")


weigt <- get_weights(classifer)

w2 <- weigt[3]
max(w2[[1]])
which(w2[[1]] == max(w2[[1]]))

w1<- as.data.frame( weigt[[1]] )
w1 <- as.data.frame( w1[,231])
w1$gene <- colnames(real)