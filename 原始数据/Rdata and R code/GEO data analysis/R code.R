
  
  library("WGCNA")                     
                               
  data=read.table("GSE24807.txt",sep="\t",header=T,check.names=F,row.names=1)     
  data<-log2(data+1)
  
  
  
  data<-data[order(apply(data,1,mad), decreasing = T)[1:5000],]                 
  
  
  datExpr0=t(data)
  
  
  gsg = goodSamplesGenes(datExpr0, verbose = 3)
  if (!gsg$allOK)
  {
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  

  sampleTree = hclust(dist(datExpr0), method = "average")
  pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
 
  dev.off()
  


  traitData=read.table("group.txt",sep="\t",header=T,check.names=F,row.names=1)     
  fpkmSamples = rownames(datExpr0)
  traitSamples =rownames(traitData)
  sameSample=intersect(fpkmSamples,traitSamples)
  datExpr0=datExpr0[sameSample,]
  datTraits=traitData[sameSample,]
  

  sampleTree2 = hclust(dist(datExpr0), method = "average")
  traitColors = numbers2colors(datTraits, signed = FALSE)
  pdf(file="2_sample_heatmap.pdf",width=12,height=12)
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  
 
  enableWGCNAThreads()   
  powers = c(1:20)       
  sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
  pdf(file="3_scale_independence.pdf",width=9,height=5)
  par(mfrow = c(1,2))
  cex1 = 0.90

  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red") 

  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  
  
  sft 
  softPower =sft$powerEstimate 
  kkk<-function(datExpr,power){
    k<-softConnectivity(datE=datExpr,
                        power= power,
                       
                        verbose = 5) 
    sizeGrWindow(10, 5)
    par(mfrow=c(1,2))
    hist(k)
    scaleFreePlot(k,main="Check Scale free topology\n")
  }
  plot(kkk(datExpr0,softPower))
  
  
  
  adjacency = adjacency(datExpr0, power = softPower)
  
 
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  
  geneTree = hclust(as.dist(dissTOM), method = "average");
  pdf(file="4_gene_clustering.pdf",width=12,height=9)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  dev.off()
  
  
 
  minModuleSize = 80     
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  pdf(file="5_Dynamic_Tree.pdf",width=8,height=6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  
  
 
  MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs);
  METree = hclust(as.dist(MEDiss), method = "average")
  pdf(file="6_Clustering_module.pdf",width=7,height=6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  MEDissThres = 0.25 
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  
  merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors
  mergedMEs = merge$newMEs
  pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  moduleColors = mergedColors
  table(moduleColors)
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  
  
  nGenes = ncol(datExpr0)
  nSamples = nrow(datExpr0)
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  pdf(file="8_Module_trait.pdf",width=8,height=8)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(5, 10, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  

  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  traitNames=names(datTraits)
  geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
  names(GSPvalue) = paste("p.GS.", traitNames, sep="")
  
  
 
  probes = colnames(datExpr0)
  geneInfo0 = data.frame(probes= probes,
                         moduleColor = moduleColors)
  for (Tra in 1:ncol(geneTraitSignificance))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                           GSPvalue[, Tra])
    names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                         names(GSPvalue)[Tra])
  }
  
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                           MMPvalue[, mod])
    names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                         names(MMPvalue)[mod])
  }
  geneOrder =order(geneInfo0$moduleColor)
  geneInfo = geneInfo0[geneOrder, ]
  write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)
  
 

  
  
  
  
  library(tidyverse)
  library(glmnet)
  source('msvmRFE.R')   
  library(VennDiagram)
  library(sigFeature)
  library(e1071)
  library(caret)
  library(randomForest)

  train<-read.table("ARGexp.txt",row.names = 1,as.is = F,header = T)
  train[1:4,1:4]
  

  
  x <- as.matrix(train[,-1])  
  (y <- ifelse(train$group == "NR", 0,1)) 
  x[1:4,1:4]
  

  set.seed(64)
  fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
  
 
  plot(fit, xvar = "dev", label = TRUE)
  
  cvfit = cv.glmnet(x, y, 
                    nfold=10, 
                    family = "binomial", type.measure = "class")
  plot(cvfit)
  
  cvfit$lambda.min 
  

  myCoefs <- coef(cvfit, s="lambda.min");
  lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
  (lasso_fea <- lasso_fea[-1])
  

  write.csv(lasso_fea,"feature_lasso.csv")
  
  
  

  predict <- predict(cvfit, newx = x[1:nrow(x),], s = "lambda.min", type = "class")
  table(predict,y)
  

  input <- train
  
  set.seed(213)
 
  svmRFE(input, k = 5, halve.above = 100) 
  
  nfold = 5
  nrows = nrow(input)
  folds = rep(1:nfold, len=nrows)[sample(nrows)]
  folds = lapply(1:nfold, function(x) which(folds == x))
  results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) 
  
  top.features = WriteFeatures(results, input, save=F) 
  head(top.features)
  
 
  write.csv(top.features,"feature_svm.csv")
  
 
  featsweep = lapply(1:10, FeatSweep.wrap, results, input) 
  featsweep
  

  no.info = min(prop.table(table(input[,1])))
  errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
  
 
  PlotErrors(errors, no.info=no.info) 
  
 
  Plotaccuracy(1-errors,no.info=no.info) 
  

  which.min(errors)
  
  (myoverlap <- intersect(lasso_fea, top.features[1:which.min(errors), "FeatureName"])) 
  
  summary(lasso_fea%in%top.features[1:which.min(errors), "FeatureName"])
  
  pdf("C_lasso_SVM_venn.pdf", width = 5, height = 3)
  grid.newpage()
  venn.plot <- venn.diagram(list(LASSO = lasso_fea, 
                                 SVM_RFE = as.character(top.features[1:which.min(errors),"FeatureName"])), NULL, 
                            fill = c("#E31A1C","#E7B800"), 
                            alpha = c(0.5,0.5), cex = 4, cat.fontface=3, 
                            category.names = c("LASSO", "SVM_RFE"), 
                            main = "Overlap")
  grid.draw(venn.plot)
  dev.off()
  
  
  
  
  
  