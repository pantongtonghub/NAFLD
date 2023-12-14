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
  
  
  
  
  
  