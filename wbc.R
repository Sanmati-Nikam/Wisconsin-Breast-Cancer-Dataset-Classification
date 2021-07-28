library(ggplot2)
library(ggcorrplot)
library(caret)
library(devtools)
library(ROCR)
library(pROC)
library(gridExtra)
library(dplyr)
library(kernlab)
library(fastDummies)    
library(rpart.plot)
library(rpart)
library(easyGgplot2)
library(factoextra)
library(nnet)
library(GGally)

#reading dataset
data<- read.csv("wbc_csv.csv")
head(data)

str(data)
summary(data)

#checking Na values
anyNA(data)

#converting to factor
data$diagnosis <- as.factor(data$diagnosis)
data<-data[-1]
#visualization
mean <- ggpairs(data[,c(2:11,1)], aes(color=diagnosis, alpha=0.75), lower=list(continuous="smooth"))+ theme_bw()+
  labs(title="Cancer_Mean")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))
mean

se <- ggpairs(data[,c(12:21,1)], aes(color=diagnosis, alpha=0.75), lower=list(continuous="smooth"))+ theme_bw()+
  labs(title="Cancer_se")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))
se

worst <- ggpairs(data[,c(22:31,1)], aes(color=diagnosis, alpha=0.75), lower=list(continuous="smooth"))+ theme_bw()+
  labs(title="Cancer_worst")+
  theme(plot.title=element_text(face='bold',color='black',hjust=0.5,size=12))
worst

a<- ggplot(data, aes(x = diagnosis, fill = diagnosis)) +
  geom_bar(stat = "count", position = "stack", show.legend = FALSE) +
  theme_minimal(base_size = 16) +
  geom_label(stat = "count", aes(label = ..count..), position = position_stack(vjust = 0.5),
             size = 5, show.legend = FALSE)
print(a)

#boxlpot
box.plot <- function(x){
  if(is.numeric(data[, x])){
    ggplot(data, aes_string('diagnosis', x)) +
      geom_boxplot( outlier.color='green',outlier.size=4,outlier.shape=5,fill=c('#FF69B4','#33FFFF')) +
      ggtitle(paste('Diagonstic plot of', x))
  }
}
col.names = names(data)
sapply(col.names, box.plot)

#HISTOGRAM PLOT
hist.plot <- function(x){
  if(is.numeric(data[, x])){
    o<-ggplot2.histogram(data=data, xName=x,
                         groupName='diagnosis', legendPosition="top",
                         alpha=0.5, addDensity=TRUE,
                         addMeanLine=TRUE, meanLineColor="white", meanLineSize=1.5,bins=30)
    print(o)
  }
}
col.names = names(data)
sapply(col.names, hist.plot)


#SCATTER PLOT
p1<-ggplot(data, aes(x = area_mean, y = perimeter_mean)) +
  geom_point(aes(color =diagnosis))
p2<-ggplot(data, aes(x = radius_mean, y = perimeter_mean)) +
  geom_point(aes(color = diagnosis))
p3<-ggplot(data, aes(x = radius_worst, y = perimeter_worst)) +
  geom_point(aes(color = diagnosis))
p4<-ggplot(data, aes(x = radius_worst, y = area_worst)) +
  geom_point(aes(color = diagnosis))
grid.arrange(p1,p2,p3,p4,nrow=2)

p1<-ggplot(data, aes(x = fractal_dimension_se, y = perimeter_worst)) +
  geom_point(aes(color =diagnosis))
p2<-ggplot(data, aes(x = fractal_dimension_mean, y = radius_se)) +
  geom_point(aes(color = diagnosis))
p3<-ggplot(data, aes(x = symmetry_se, y = texture_mean)) +
  geom_point(aes(color = diagnosis))
p4<-ggplot(data, aes(x = symmetry_se, y = smoothness_worst)) +
  geom_point(aes(color = diagnosis))
grid.arrange(p1,p2,p3,p4,nrow=2)

#Correlation
data <- dummy_cols(data, select_columns = c('diagnosis'))

corr <- round(cor(data[,2:33], use="complete.obs"), 2)
options(repr.plot.width=16, repr.plot.height=16)
ggcorrplot(corr, lab = TRUE, colors = c("mediumspringgreen", "white", "lawngreen"),
           show.legend = T, outline.color = "gray", type = "full", hc.order = T,  
           tl.cex = 10, lab_size = 2, sig.level = .3) +
  labs(fill = "Correlation")

#PCA
data <- select(data,-c(diagnosis_B,diagnosis_M))

pca = prcomp(data[-1],scale. = TRUE, center = TRUE)
print(pca)
summary(pca)
fviz_pca(pca, col.ind = data$diagnosis, col="red",
         palette = "jco", geom = "point", repel=TRUE,
         legend.title="Diagnosis", addEllipses = TRUE)

fviz_eig(pca, addlabels=TRUE, ylim=c(0,60), geom = c("bar", "line"), barfill = "mediumseagreen", barcolor="grey",linecolor = "red", ncp=10,ggtheme = theme_minimal(base_size = 12))+
  labs(title = "Cancer All Variances - PCA",
       x = "Principal Components", y = "% of explained variances")

gof = (pca$sdev)^2/sum((pca$sdev)^2)
sum(gof[1:6])

#creating dataset using pca
newdata = pca$scores[,1:6]
newdata <- as_tibble(pca$x)[1:6]
newdata = cbind(data$diagnosis,newdata)
colnames(newdata) = c("diagnosis","p1","p2","p3","p4","p5","p6")
newdata=as.data.frame(newdata)
head(newdata)

newdata$p1 = as.numeric(newdata$p1)
newdata$p2 = as.numeric(newdata$p2)
newdata$p3 = as.numeric(newdata$p3)
newdata$p4 = as.numeric(newdata$p4)
newdata$p5 = as.numeric(newdata$p5)
newdata$p6 = as.numeric(newdata$p6)
newdata$diagnosis <- as.factor(newdata$diagnosis)


#spliting the dataset
set.seed(123)
inTrain <- createDataPartition(newdata$diagnosis , p=0.75, list=FALSE)
training <- newdata[inTrain,]
testing <- newdata[-inTrain,]


#Logistic
set.seed(123)
#logistic regression
m<-glm(diagnosis~.,data=training,family =binomial)
summary(m)

m<-glm(diagnosis~p1+p2+p3+p4+p5,data=training,family =binomial)
summary(m)


#prediction for training
prob <- predict(m, training, type="response")
pred_train <- factor(prob > .7, levels=c(FALSE, TRUE), labels = c("B","M"))
cat("\ntraining data")
conf_Mat_train <- table(training$diagnosis,pred_train, dnn = c("Actual", "Predicted") )
print(conf_Mat_train)
cat("\n1.Accuracy of the model is :",sum(diag(conf_Mat_train))/sum(conf_Mat_train))
cat("\n2.Sensitivity of the model is:",conf_Mat_train[2,2]/sum(conf_Mat_train[2,]))
cat("\n3.Specificity of the model is:",conf_Mat_train[1,1]/sum(conf_Mat_train[1,]))
cat("\n4.Precision of the model is:",conf_Mat_train[2,2]/sum(conf_Mat_train[,2]))


#Prediction for testing
prob <- predict(m, testing, type="response")
pred <- factor(prob > .7, levels=c(FALSE, TRUE), labels = c("B","M"))
cat("\ntesting data")
conf_Mat <- table(testing$diagnosis,pred, dnn = c("Actual", "Predicted") )
print(conf_Mat)
acc_lr = sum(diag(conf_Mat))/sum(conf_Mat)
sen_lr = conf_Mat[2,2]/sum(conf_Mat[2,])
spec_lr = conf_Mat[1,1]/sum(conf_Mat[1,])
prec_lr = conf_Mat[2,2]/sum(conf_Mat[,2])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))


#Resampling
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats=3,
                           classProbs = T,
                           summaryFunction = twoClassSummary)


#SVM
set.seed(1234)
svm<-train(diagnosis ~.,
           data = training,
           method = "svmPoly",
           metric="ROC",
           trControl = fitControl,
           preProc=c("center","scale"))
svm
print(svm)
a<-plot(svm)
print(a)

pred_svm<-predict(svm,training)
conf_Mat<-table(actual=training$diagnosis,predicted=pred_svm)
cat("\nConfusion Matrix for testing\n")
print(conf_Mat)
acc_svm = sum(diag(conf_Mat))/sum(conf_Mat)
sen_svm = conf_Mat[2,2]/sum(conf_Mat[2,])
spec_svm = conf_Mat[1,1]/sum(conf_Mat[1,])
prec_svm = conf_Mat[2,2]/sum(conf_Mat[,2])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))




pred_svm<-predict(svm,testing)
conf_Mat<-table(actual=testing$diagnosis,predicted=pred_svm)
cat("\nConfusion Matrix for testing\n")
print(conf_Mat)
acc_svm = sum(diag(conf_Mat))/sum(conf_Mat)
sen_svm = conf_Mat[2,2]/sum(conf_Mat[2,])
spec_svm = conf_Mat[1,1]/sum(conf_Mat[1,])
prec_svm = conf_Mat[2,2]/sum(conf_Mat[,2])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))

#SVM-linear ,training=98.1,testing=97.1- overfitting
#SVM-radial ,training=97.6,testing=97.8- ok
#SVM-ploy ,training=97.1,testing=97.8- overfitting
#SVM-sigmoid ,training=100,testing=97.1- overfitting



#KNN
cat("\nKNN\n")
set.seed(123)
knnFit <- train(diagnosis ~.,
                data = training,
                method = "knn",
                metric = "ROC",
                trControl = fitControl,
                preProc=c("center","scale"),
                tuneGrid=expand.grid(k=1:40))
knnFit
print(knnFit)
a<-plot(knnFit)
print(a)
pred_knn<-predict(knnFit,testing)
conf_Mat<-table(actual=testing$diagnosis,predicted=pred_knn)
cat("\nConfusion Matrix for testing\n")
print(conf_Mat)
acc_knn = sum(diag(conf_Mat))/sum(conf_Mat)
sen_knn = conf_Mat[2,2]/sum(conf_Mat[2,])
spec_knn = conf_Mat[1,1]/sum(conf_Mat[1,])
prec_knn = conf_Mat[2,2]/sum(conf_Mat[,2])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))
varImp(knnFit)


set.seed(123)
#NEURAL NETWORKS
nnetFit <- train(diagnosis ~.,
                 data = training,
                 method = "nnet",
                 metric = "ROC",
                 trControl = fitControl,
                 tuneGrid = expand.grid(.decay = c(2,1,0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7), .size = c(3, 5, 7,10)),
                 verbose = TRUE,
                 maxit=10)

nnetFit
plot(nnetFit)
pred_nnet<-predict(nnetFit,testing)
conf_Mat<-table(actual=testing$diagnosis,predicted=pred_nnet)
cat("\nConfusion Matrix for testing\n")
print(conf_Mat)
acc_nnet = sum(diag(conf_Mat))/sum(conf_Mat)
sen_nnet = conf_Mat[2,2]/sum(conf_Mat[2,])
spec_nnet = conf_Mat[1,1]/sum(conf_Mat[1,])
prec_nnet = conf_Mat[2,2]/sum(conf_Mat[,2])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))


set.seed(123)
#DeCISION TREE
dtree <- rpart(diagnosis ~ .,
               data = training, method = "class",
               control=rpart.control(minsplit=20,
                                     maxdepth=3,
                                     cp=0.01),
)
print(dtree)
prp(dtree,
    type = 4,
    extra = 101,
    nn = TRUE,
    tweak = 1,
    space = 0.1,
    fallen.leaves = FALSE,
    roundint = FALSE,box.palette="RdBu", shadow.col="gray")

pred_dt<-predict(dtree,training,type="class")
conf_Mat<-table(actual=training$diagnosis,predicted=pred_dt)
cat("\nConfusion Matrix for testing\n")
print(conf_Mat)
acc_dt = sum(diag(conf_Mat))/sum(conf_Mat)
sen_dt = conf_Mat[2,2]/sum(conf_Mat[2,])
spec_dt = conf_Mat[1,1]/sum(conf_Mat[1,])
prec_dt = conf_Mat[2,2]/sum(conf_Mat[,2])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))


pred_dt<-predict(dtree,testing,type="class")
conf_Mat<-table(actual=testing$diagnosis,predicted=pred_dt)
cat("\nConfusion Matrix for testing\n")
print(conf_Mat)
acc_dt = sum(diag(conf_Mat))/sum(conf_Mat)
sen_dt = conf_Mat[2,2]/sum(conf_Mat[2,])
cat("\n1.Accuracy of the model is:",sum(diag(conf_Mat))/sum(conf_Mat))
cat("\n2.Sensitivity of the model is:",conf_Mat[2,2]/sum(conf_Mat[2,]))
cat("\n3.Specificity of the model is:",conf_Mat[1,1]/sum(conf_Mat[1,]))
cat("\n4.Precision of the model is:",conf_Mat[2,2]/sum(conf_Mat[,2]))

#pefomance comparision
df <- as.data.frame(df)
df <- rbind("Logistic","SVM","KNN","NNET","DTREE")
df2 <- rbind(round(acc_lr*100,4),round(acc_svm*100,4),round(acc_knn*100,4),round(acc_nnet*100,4),round(acc_dt*100))
df3 <- rbind(round(sen_lr*100,4),round(sen_svm*100,4),round(sen_knn*100,4),round(sen_nnet*100,4),round(sen_dt*100))
df4 <- rbind(round(spec_lr*100,4),round(spec_svm*100,4),round(spec_knn*100,4),round(spec_nnet*100,4),round(spec_dt*100))
df5 <- rbind(round(prec_lr*100,4),round(prec_svm*100,4),round(prec_knn*100,4),round(prec_nnet*100,4),round(prec_dt*100))
df <- cbind(df,df2,df3,df4,df5)
df <- as.data.frame(df)
names(df)[1] <- "Algorithm"
names(df)[2] <- "Accuracy"
names(df)[3] <- "Sensitivity"
names(df)[4] <- "Specificity"
names(df)[5] <- "Precision"


df$Accuracy  = as.numeric(df$Accuracy)
df$Sensitivity  = as.numeric(df$Sensitivity)
df$Specificity  = as.numeric(df$Specificity)
df$Precision = as.numeric(df$Precision)


p1 <-ggplot(df, aes(x = Algorithm,
                    y = Accuracy, fill=Algorithm)) +
  geom_bar(stat = "identity",
           position = "dodge",
           alpha=0.9) +geom_text(aes(label=Accuracy) , 
                                 position = position_dodge(width = 0.6), vjust = -0.2)

p1

p2 <-ggplot(df, aes(x = Algorithm,
                    y = Precision, fill=Algorithm)) +
  geom_bar(stat = "identity",
           position = "dodge",
           alpha=0.9) +geom_text(aes(label=Precision) , 
                                 position = position_dodge(width = 0.9), vjust = -0.2)

p2
p3 <-ggplot(df, aes(x = Algorithm,
                    y = Specificity, fill=Algorithm)) +
  geom_bar(stat = "identity",
           position = "dodge",
           alpha=0.9) +geom_text(aes(label=Specificity) , 
                                 position = position_dodge(width = 0.9), vjust = -0.2)


p4 <-ggplot(df, aes(x = Algorithm,
                    y = Sensitivity, fill=Algorithm)) +
  geom_bar(stat = "identity",
           position = "dodge",
           alpha=0.9) +geom_text(aes(label=Sensitivity) , 
                                 position = position_dodge(width = 0.9), vjust = -0.2)

grid.arrange(p2,p1,nrow=1)
grid.arrange(p3,p4,nrow=1)

# 
# lr_pred = round(predict(m, testing, type="response"),0) # Create a vector of predicitons made from the test/validation data set for the linear model.
# svm_pred = predict(svm, testing) #Prediction vector for the SVM.
# nnet_pred = predict(nnetFit,testing, type="raw") #Prediction vector for the neural network.
# knn_pred = predict(svm, testing)

pred_svm <- as.numeric(pred_svm)
pred_knn <- as.numeric(pred_knn)
pred_nnet <- as.numeric(pred_nnet)
pred_dt <- as.numeric(pred_dt)



rocobj1 <- roc(testing$diagnosis, prob)
auc <- round(auc(testing$diagnosis, prob),4)

rocobj2 <- roc(testing$diagnosis, pred_svm)
auc <- round(auc(testing$diagnosis, pred_svm),4)

rocobj3 <- roc(testing$diagnosis, pred_nnet)
auc <- round(auc(testing$diagnosis, pred_nnet),4)

rocobj4 <- roc(testing$diagnosis, pred_knn)
auc <- round(auc(testing$diagnosis, pred_knn),4)

rocobj5 <- roc(testing$diagnosis, pred_dt)
auc <- round(auc(testing$diagnosis, pred_dt),4)

#create ROC plot
g1 <- ggroc(rocobj1, colour = 'steelblue', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()


g2 <- ggroc(rocobj2, colour = 'springgreen', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()

g3 <- ggroc(rocobj3, colour = 'goldenrod1', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()

g4 <- ggroc(rocobj4, colour = 'darksalmon', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()

g5 <- ggroc(rocobj5, colour = 'purple', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()


ggroc(list(Logistic = rocobj1, SVM = rocobj2, NNET = rocobj3, KNN = rocobj4, DTREE= rocobj5))



