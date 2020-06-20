#require(doMC)
require(pROC) # for AUC-CIs
require(plotROC)
require(ROCR)
require(verification)

require(caret)
require(caTools)
require(caretEnsemble)
require(e1071)

require(pcaMethods)

require(tidyverse)
require(corrr)
require(igraph)
require(ggraph)
require(ggthemes)
require(ggpubr)

#registerDoMC(cores=2)


opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

ncoefs = function(mod_fit_obj){
  cc<-data.frame(coef.name = dimnames(coef(mod_fit_obj$finalModel,s=mod_fit_obj$bestTune$lambda))[[1]], coef.value = matrix(coef(mod_fit_obj$finalModel,s=mod_fit_obj$bestTune$lambda)))
  cc1<-cc[cc$coef.value!=0,'coef.value']
  cc2<-cc[cc$coef.value!=0,]
  # extract names of non zero elastic net coefficients
  nms<-cc[cc$coef.value!=0,1] 
  I_index<-match(nms[seq(2,length(nms))],names(dat3))
  return(nms)
}

dat_lesion_txt<-read.csv('KDI/20170620_data_list_for_Gioli_clean.csv',header=TRUE)

dat_lesion_img<-read.csv('AllSamples_Substraction_SUB_20170215_5_avg.csv',header=FALSE)

names(dat_lesion_img)[2]<-"IncognitoPatientID"

I_=seq(15,dim(dat_lesion_img)[2])

# exclude subject 23000IQQEA
dat_lesion_txt<-dat_lesion_txt[-which(dat_lesion_txt$IncognitoPatientID=='23000IQQEA'),]

dat<-merge(dat_lesion_txt,dat_lesion_img[,c(2,I_)],by="IncognitoPatientID")

dat$Mutationstatus<-gsub("BRCA1","case",dat$Mutationstatus)
dat$Mutationstatus<-gsub("BRCA2","case",dat$Mutationstatus)
dat$Mutationstatus<-as.factor(dat$Mutationstatus)
levels(dat$Mutationstatus)<-c("CASES","CONTROLS")

# rename some variables
names(dat)[names(dat)=="AlterbeiderDiagnose"]="Age"
names(dat)[names(dat)=="KontralatMaca"]="BBC"
names(dat)[names(dat)=="MaCaFstline"]="FDR-BC"
names(dat)[names(dat)=="OcVaFstline"]="FDR-OC"
names(dat)[names(dat)=="MaCaSecline"]="SDR-BC"
names(dat)[names(dat)=="OcVaSecstline"]="SDR-OC"
names(dat)[names(dat)=="MaCaThrline"]="TDR-BC"
names(dat)[names(dat)=="OcVaThrline"]="TDR-OC"
names(dat)[names(dat)=="MaCaMli_nnlich"]="MBC"


dat$Imaging_mean<-rowMeans(dat[,grep("^V",names(dat))])

dat2<-data.frame(dat[,c(1,2,seq(4,12),13,14)],scale(dat[,-c(1,2,seq(4,12),13,14)]))

imaPCA2 <- pca(dat2[,grep("^V",names(dat2))], scale = "uv", center = T, nPcs = 41, method = "svd")

imaPCA_caret=predict(preProcess(dat2[,grep("^V",names(dat2))], method=c("pca")),dat2[,grep("^V",names(dat2))])

dat3<-data.frame(dat2,imaPCA_caret)

I_pca<-grep("PC",names(dat3))

I_m <- grep("Imaging_mean",names(dat3))
I_imaging <- grep("^V",names(dat3))
I_clinical<-c(12,14)
#I_FamilyHistory_1<-c(seq(3,4))
#I_FamilyHistory_23<-c(seq(5,9))
#I_FamilyHistory_BC<-c(3,5,7)
#I_FamilyHistory_OC<-c(4,6,8)

I_FamilyHistory_1<-c(3)
I_FamilyHistory_23<-c(5,7)


I_TNBC<-c(13)

model_weights <- ifelse(dat3$Mutationstatus=="CASES",(1/table(dat3$Mutationstatus)[1])*.5,(1/table(dat3$Mutationstatus)[2])*.5)
ctrl<-trainControl(method = "repeatedcv", number=5, repeats=20, returnResamp = "final", classProbs = TRUE,verboseIter=F,summaryFunction=twoClassSummary)

# randomly permute labels
set.seed(50)
I_<-sample(41)

# TNBC
set.seed(849)
mod_fit_tnbc <- train(Mutationstatus ~ ., data=dat3[,c(2,I_TNBC)],
                                 method="glm",
                                 metric="ROC",
                                 family='binomial',
                                # weights = model_weights,
                                 trControl=ctrl)

modcvpred_tnbc <- predict(object=mod_fit_tnbc, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_tnbc  > t , 1,0)
pred <- prediction(modcvpred_tnbc , dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_TNBC <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
#correct way confusionMatrix(predict(mod_fit_tnbc,dat3),dat3$Mutationstatus)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_TNBC)
ncoefs_tnbc<-1
auc_ci_tnbc<-ci(dat3$Mutationstatus,modcvpred_tnbc$CONTROLS)

# Clinical
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=lambda.seq)
set.seed(849)

mod_fit_clinical <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_clinical)]),
                      method="glm",
                      metric="ROC",
                      family='binomial',
                     # tuneGrid=elgrid,
                      #weights = model_weights,
                      trControl=ctrl)

modcvpred_clinical <- predict(object=mod_fit_clinical, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_clinical  > t , 1,0)
pred <- prediction(modcvpred_clinical , dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_CLINICAL <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_CLINICAL)
ncoefs_clinical<-2
auc_ci_clinical<-ci(dat3$Mutationstatus,modcvpred_clinical$CONTROLS)

# Family history
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=lambda.seq)
set.seed(849)
mod_fit_FH <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_FamilyHistory_1,I_FamilyHistory_23)]),
                          method="glmnet",
                          metric="ROC",
                          family='binomial',
                          tuneGrid=elgrid,
                          weights = model_weights,
                          trControl=ctrl)

modcvpred_FH <- predict(object=mod_fit_FH, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_FH  > t , 1,0)
pred <- prediction(modcvpred_FH, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_FH <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_FH)
ncoefs_fh<-length(ncoefs(mod_fit_FH))-1
auc_ci_fh<-ci(dat3$Mutationstatus,modcvpred_FH$CONTROLS)


#Imaging only (with pca and elastic net for model selection
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=lambda.seq)
set.seed(849)

mod_fit_elastic_imaging <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_pca)]),
                                    method="glmnet",
                                    metric="ROC",
                                    family='binomial',
                                    weights = model_weights,
                                    tuneGrid=elgrid,
                                    trControl=ctrl)

modcvpred_elastic_imaging <- predict(object=mod_fit_elastic_imaging, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_elastic_imaging  > t , 1,0)
pred <- prediction(modcvpred_elastic_imaging , dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_IMAGING <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_IMAGING)
plot(varImp(mod_fit_elastic_imaging,scale = FALSE),5)
auc_ci_imaging<-ci(dat3$Mutationstatus,modcvpred_elastic_imaging$CONTROLS)

ncoefs_imaging<-length(ncoefs(mod_fit_elastic_imaging))-1
coefs_imaging<-ncoefs(mod_fit_elastic_imaging)

I_index<-match(coefs_imaging[seq(2,length(coefs_imaging))],names(dat3))

# Output for Question 1

df_tnbc<-data.frame(ID=dat3$IncognitoPatientID,D = c(1,2)[dat3$Mutationstatus]-1, D.str = dat3$Mutationstatus,TNBC = logit(modcvpred_tnbc$CONTROLS), stringsAsFactors = FALSE)

df_clinical<-data.frame(Clinical = modcvpred_clinical$CONTROLS, stringsAsFactors = FALSE)

df_fh<-data.frame(Family.history = modcvpred_FH$CONTROLS, stringsAsFactors = FALSE)

df_imaging<-data.frame(Imaging = modcvpred_elastic_imaging$CONTROLS, stringsAsFactors = FALSE)

df=cbind(df_tnbc,df_clinical,df_fh,df_imaging)

df_all <- melt_roc(df, "D", c("TNBC","Clinical","Family.history","Imaging"))

names(df_all)[3]<-"Data_Type"
df_all$Data_Type<-factor(df_all$Data_Type)

vi_tnbc<-data.frame(varImp(mod_fit_tnbc)[[1]])
vi_clinical<-data.frame(varImp(mod_fit_clinical)[[1]])
vi_fh<-data.frame(varImp(mod_fit_FH)[[1]])
vi_imaging<-data.frame(varImp(mod_fit_elastic_imaging)[[1]])

names(vi_tnbc)<-"Importance"
names(vi_clinical)<-"Importance"
names(vi_fh)<-"Importance"
names(vi_imaging)<-"Importance"

vi_tnbc$Variable<-as.factor(rownames(vi_tnbc))
vi_tnbc$Data_Type<-"TNBC"
vi_tnbc$Variable <- factor(vi_tnbc$Variable)
vi_tnbc$Importance <- 100

vi_clinical$Variable<-as.factor(rownames(vi_clinical))
vi_clinical$Data_Type<-"Clinical"
vi_clinical$Variable <- factor(vi_clinical$Variable, levels = vi_clinical$Variable[order(vi_clinical$Importance)])

vi_fh$Variable<-as.factor(rownames(vi_fh))
vi_fh$Data_Type<-"Family.history"
vi_fh$Variable <- factor(vi_fh$Variable, levels = vi_fh$Variable[order(vi_fh$Importance)])

vi_imaging$Variable<-as.factor(rownames(vi_imaging))
vi_imaging$Data_Type<-"Imaging"
vi_imaging$Variable <- factor(vi_imaging$Variable, levels = vi_imaging$Variable[order(vi_imaging$Importance,decreasing=T)])

vi_Fig2<-rbind(vi_tnbc,vi_clinical,vi_fh,vi_imaging[order(vi_imaging$Importance,decreasing=T)[1:3],])

str1=paste("AUC (TNBC) =", sprintf('%.3f',AUC_EN_TNBC),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_tnbc)[1]),"-",sprintf('%.3f',as.vector(auc_ci_tnbc)[3]))
str2=paste("AUC (Clinical) =", sprintf('%.3f',AUC_EN_CLINICAL),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_clinical)[1]),"-",sprintf('%.3f',as.vector(auc_ci_clinical)[3]))
str3=paste("AUC (Family.history) =", sprintf('%.3f',AUC_EN_FH),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_fh)[1]),"-",sprintf('%.3f',as.vector(auc_ci_fh)[3]))
str4=paste("AUC (Imaging) =", sprintf('%.3f',AUC_EN_IMAGING),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_imaging)[1]),"-",sprintf('%.3f',as.vector(auc_ci_imaging)[3]))

lp1<-ggplot(df_all, aes(d = D, m = M, group=Data_Type,color=Data_Type)) + geom_roc(n.cuts=0,size=1.5)+style_roc(guide=F,xlab = "1 - Specificity",ylab="Sensitivity")+annotate("text", size=3, x = .25, y = .12,hjust=0,label=paste(str1,str2,str3,str4,sep="\n"),color="black")
lp1<-lp1+scale_colour_Publication()+theme(legend.title = element_text(size=18, face="bold"),legend.text = element_text(size=18))+theme_Publication()

lp2<-ggplot(vi_Fig2, aes(reorder_within(Data_Type, Importance, Variable),y=Importance,fill=Data_Type,label=Variable))+coord_flip()+
  geom_col(position = position_dodge2(width = .5, preserve = "single"))+geom_text(position = position_dodge2(width = .5), angle = 0,hjust=-.1)+
  scale_x_reordered()+facet_grid(Data_Type~., scales = "free_y",space="free_y")+scale_fill_Publication()+theme_Publication()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(5,5,5,17),"mm"))+
  theme(strip.text.y = element_blank())+expand_limits(y = c(-1, 110))+guides(fill=guide_legend(nrow=1,byrow=TRUE))

lp_2_tmp<-ggarrange(lp1+rremove("legend"),lp2+rremove("y.text")+rremove("xylab")+theme(legend.position = "top",legend.title = element_blank()),labels = c("A","B"),nrow=2,ncol=1, align="h")

lp_Fig2<-annotate_figure(lp_2_tmp, bottom = text_grob("Variable Importance (%)", color = "black", size = 12))

ggsave(file=paste("C:/Users/z003tmpd.AD005/Desktop/My Documents/KDI/pics_June2019/Q1/Fig2",".jpg",sep=""),plot=lp_Fig2,width=210, height=297, units="mm",dpi=600)

#Question 2
#Imaging + clinical  (with pca and elastic net for model selection
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_imaging_tnbc <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_TNBC,I_index)]),
                                 method="glmnet",
                                 metric="ROC",
                                 family='binomial',
                                 weights = model_weights,
                                 tuneGrid=elgrid,
                                 trControl=ctrl)

modcvpred_imaging_tnbc <- predict(object=mod_fit_imaging_tnbc, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_imaging_tnbc  > t , 1,0)
pred <- prediction(modcvpred_imaging_tnbc, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_imaging_tnbc <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_imaging_tnbc)
ncoefs_imaging_tnbc<-length(ncoefs(mod_fit_imaging_tnbc))-1

auc_imaging_tnbc<-ci(dat3$Mutationstatus,modcvpred_imaging_tnbc$CONTROLS)


#Imaging + clinical  (with pca and elastic net for model selection
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_imaging_clinical <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_clinical,I_index)]),
                                      method="glmnet",
                                      metric="ROC",
                                      family='binomial',
                                      weights = model_weights,
                                      tuneGrid=elgrid,
                                      trControl=ctrl)

modcvpred_imaging_clinical <- predict(object=mod_fit_imaging_clinical, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_imaging_clinical  > t , 1,0)
pred <- prediction(modcvpred_imaging_clinical, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_imaging_clinical <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_imaging_clinical)
ncoefs_imaging_clinical<-length(ncoefs(mod_fit_imaging_clinical))-1

auc_imaging_clinical<-ci(dat3$Mutationstatus,modcvpred_imaging_clinical$CONTROLS)

#Imaging + family history  (with pca and elastic net for model selection
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_imaging_FH <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_FamilyHistory_1,I_FamilyHistory_23,I_index)]),
                                          method="glmnet",
                                          metric="ROC",
                                          family='binomial',
                                          weights = model_weights,
                                          tuneGrid=elgrid,
                                          trControl=ctrl)

modcvpred_imaging_FH <- predict(object=mod_fit_imaging_FH, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_imaging_FH  > t , 1,0)
pred <- prediction(modcvpred_imaging_FH, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_imaging_FH <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_imaging_FH)
ncoefs_imaging_FH<-length(ncoefs(mod_fit_imaging_FH))-1
auc_imaging_fh<-ci(dat3$Mutationstatus,modcvpred_imaging_FH$CONTROLS)

# Output for Question 2

df_imaging_tnbc<-data.frame(ID=dat3$IncognitoPatientID,D = c(1,2)[dat3$Mutationstatus]-1, D.str = dat3$Mutationstatus,Imaging.TNBC = modcvpred_imaging_tnbc$CONTROLS, stringsAsFactors = FALSE,check.names=F)

df_imaging_clinical<-data.frame(Imaging.Clinical = modcvpred_imaging_clinical$CONTROLS, stringsAsFactors = FALSE)

df_imaging_FH<-data.frame(Imaging.Family.history = modcvpred_imaging_FH$CONTROLS, stringsAsFactors = FALSE)

df=cbind(df_imaging_tnbc,df_imaging_clinical,df_imaging_FH)

df_all <- melt_roc(df, "D", c("Imaging.TNBC","Imaging.Clinical","Imaging.Family.history"))

names(df_all)[3]<-"Data_Type"

vi_imaging_tnbc<-data.frame(varImp(mod_fit_imaging_tnbc)[[1]])
vi_imaging_clinical<-data.frame(varImp(mod_fit_imaging_clinical)[[1]])
vi_imaging_fh<-data.frame(varImp(mod_fit_imaging_FH)[[1]])

names(vi_imaging_tnbc)<-"Importance"
names(vi_imaging_clinical)<-"Importance"
names(vi_imaging_fh)<-"Importance"


vi_imaging_tnbc$Variable<-as.factor(rownames(vi_imaging_tnbc))
vi_imaging_tnbc$Data_Type<-"Imaging.TNBC"
vi_imaging_tnbc$Variable <- factor(vi_imaging_tnbc$Variable, levels = vi_imaging_tnbc$Variable[order(vi_imaging_tnbc$Importance,decreasing=T)])

vi_imaging_clinical$Variable<-as.factor(rownames(vi_imaging_clinical))
vi_imaging_clinical$Data_Type<-"Imaging.Clinical"
vi_imaging_clinical$Variable <- factor(vi_imaging_clinical$Variable, levels = vi_imaging_clinical$Variable[order(vi_imaging_clinical$Importance,decreasing=T)])

vi_imaging_fh$Variable<-as.factor(rownames(vi_imaging_fh))
vi_imaging_fh$Data_Type<-"Imaging.Family.history"
vi_imaging_fh$Variable <- factor(vi_imaging_fh$Variable, levels = vi_imaging_fh$Variable[order(vi_imaging_fh$Importance,decreasing=T)])

#vi_tmp<-merge(vi_imaging_tnbc,vi_imaging_clinical,all=T)
#vi_Fig3<-merge(vi_imaging_fh,vi_tmp,all=T)
vi_Fig3<-rbind(vi_imaging_tnbc,vi_imaging_clinical,vi_imaging_fh)


str1=paste("AUC (Imaging.Clinical) =", sprintf('%.3f',AUC_EN_imaging_clinical),"; 95% CI:", sprintf('%.3f',as.vector(auc_imaging_clinical)[1]),"-",sprintf('%.3f',as.vector(auc_imaging_clinical)[3]))
str2=paste("AUC (Imaging.Family.history) =", sprintf('%.3f',AUC_EN_imaging_FH),"; 95% CI:", sprintf('%.3f',as.vector(auc_imaging_fh)[1]),"-",sprintf('%.3f',as.vector(auc_imaging_fh)[3]))
str3=paste("AUC (Imaging.TNBC) =", sprintf('%.3f',AUC_EN_imaging_tnbc),"; 95% CI:", sprintf('%.3f',as.vector(auc_imaging_tnbc)[1]),"-",sprintf('%.3f',as.vector(auc_imaging_tnbc)[3]))

lp1<-ggplot(df_all, aes(d = D, m = M, group=Data_Type,color=Data_Type)) + geom_roc(n.cuts=0,size=1.5)+style_roc(guide=F,xlab = "1 - Specificity",ylab="Sensitivity")+annotate("text", size=3, x = .25, y = .12,hjust=0,label=paste(str1,str2,str3,sep="\n"),color="black")
lp1<-lp1+scale_colour_Publication()+theme(legend.title = element_text(size=18, face="bold"),legend.text = element_text(size=18))+theme_Publication()

lp2<-ggplot(vi_Fig3, aes(reorder_within(Data_Type, Importance, Variable),y=Importance,fill=Data_Type,label=Variable))+coord_flip()+
  geom_col(position = position_dodge2(width = .5, preserve = "single"))+geom_text(position = position_dodge2(width = .5), angle = 0,hjust=-.1)+
  scale_x_reordered()+facet_grid(Data_Type~., scales = "free_y",space="free_y")+scale_fill_Publication()+theme_Publication()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(5,5,5,17),"mm"))+
  theme(strip.text.y = element_blank())+expand_limits(y = c(-1, 110))+guides(fill=guide_legend(nrow=1,byrow=TRUE))

lp_3_tmp<-ggarrange(lp1+rremove("legend"),lp2+rremove("y.text")+rremove("xylab")+theme(legend.position = "top",legend.title = element_blank()),labels = c("A","B"),nrow=2,ncol=1, align="h")

lp_3<-annotate_figure(lp_3_tmp, bottom = text_grob("Variable Importance (%)", color = "black", size = 12))

ggsave(file=paste("KDI/pics_June2019/Q2/Fig3",".jpg",sep=""),plot=lp_3,width=210, height=297, units="mm",dpi=600)


#Question 3

# all except imaging
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_noimaging <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_TNBC,I_clinical,I_FamilyHistory_1,I_FamilyHistory_23)]),
                                    method="glmnet",
                                    metric="ROC",
                                    family='binomial',
                                    weights = model_weights,
                                    tuneGrid=elgrid,
                                    trControl=ctrl)

modcvpred_noimaging <- predict(object=mod_fit_noimaging, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_noimaging  > t , 1,0)
pred <- prediction(modcvpred_noimaging, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_noimaging <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_noimaging)
ncoefs_noimaging<-length(ncoefs(mod_fit_noimaging))-1
auc_ci_noimaging<-ci(dat3$Mutationstatus,modcvpred_noimaging$CONTROLS)

#Imaging + no family history  (with pca and elastic net for model selection
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_imaging_noFH <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_TNBC,I_clinical,I_index)]),
                                    method="glmnet",
                                    metric="ROC",
                                    family='binomial',
                                    weights = model_weights,
                                    tuneGrid=elgrid,
                                    trControl=ctrl)

modcvpred_imaging_noFH <- predict(object=mod_fit_imaging_noFH, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_imaging_noFH  > t , 1,0)
pred <- prediction(modcvpred_imaging_noFH, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_imaging_noFH <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_imaging_noFH)
ncoefs_imaging_nofh<-length(ncoefs(mod_fit_imaging_noFH))-1
auc_ci_imaging_nofh<-ci(dat3$Mutationstatus,modcvpred_imaging_noFH$CONTROLS)

#ALL
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_all <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_TNBC,I_clinical,I_FamilyHistory_1,I_FamilyHistory_23,I_index)]),
                                      method="glmnet",
                                      metric="ROC",
                                      family='binomial',
                                      weights = model_weights,
                                      tuneGrid=elgrid,
                                      trControl=ctrl)

modcvpred_all <- predict(object=mod_fit_all, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_all  > t , 1,0)
pred <- prediction(modcvpred_all, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_all <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_all)
ncoefs_all<-length(ncoefs(mod_fit_all))-1
auc_ci_all<-ci(dat3$Mutationstatus,modcvpred_all$CONTROLS)

# Output for Question 3

df_noimaging<-data.frame(ID=dat3$IncognitoPatientID,D = c(1,2)[dat3$Mutationstatus]-1, D.str = dat3$Mutationstatus,TNBC.Clinical.Family.history = modcvpred_noimaging$CONTROLS, stringsAsFactors = FALSE)

df_imaging_noFH<-data.frame(Imaging.TNBC.Clinical = modcvpred_imaging_noFH$CONTROLS, stringsAsFactors = FALSE)

df_imaging_all<-data.frame(Imaging.TNBC.Clinical.Family.history = modcvpred_all$CONTROLS, stringsAsFactors = FALSE)

df=cbind(df_noimaging,df_imaging_noFH,df_imaging_all)

df_all <- melt_roc(df, "D", c("TNBC.Clinical.Family.history","Imaging.TNBC.Clinical","Imaging.TNBC.Clinical.Family.history"))

names(df_all)[3]<-"Data_Type"

vi_tnbc_clinical_fh<-data.frame(varImp(mod_fit_noimaging)[[1]])
vi_imaging_tnbc_clinical<-data.frame(varImp(mod_fit_imaging_noFH)[[1]])
vi_imaging_tnbc_clinical_fh<-data.frame(varImp(mod_fit_all)[[1]])

names(vi_tnbc_clinical_fh)<-"Importance"
names(vi_imaging_tnbc_clinical)<-"Importance"
names(vi_imaging_tnbc_clinical_fh)<-"Importance"

vi_tnbc_clinical_fh$Variable<-as.factor(rownames(vi_tnbc_clinical_fh))
vi_tnbc_clinical_fh$Data_Type<-"TNBC.Clinical.Family.history"
vi_tnbc_clinical_fh$Variable <- factor(vi_tnbc_clinical_fh$Variable, levels = vi_tnbc_clinical_fh$Variable[order(vi_tnbc_clinical_fh$Importance,decreasing=T)])

vi_imaging_tnbc_clinical$Variable<-as.factor(rownames(vi_imaging_tnbc_clinical))
vi_imaging_tnbc_clinical$Data_Type<-"Imaging.TNBC.Clinical"
vi_imaging_tnbc_clinical$Variable <- factor(vi_imaging_tnbc_clinical$Variable, levels = vi_imaging_tnbc_clinical$Variable[order(vi_imaging_tnbc_clinical$Importance,decreasing=T)])

vi_imaging_tnbc_clinical_fh$Variable<-as.factor(rownames(vi_imaging_tnbc_clinical_fh))
vi_imaging_tnbc_clinical_fh$Data_Type<-"Imaging.TNBC.Clinical.Family.history"
vi_imaging_tnbc_clinical_fh$Variable <- factor(vi_imaging_tnbc_clinical_fh$Variable, levels = vi_imaging_tnbc_clinical_fh$Variable[order(vi_imaging_tnbc_clinical_fh$Importance,decreasing=T)])

#vi_tmp<-merge(vi_imaging_tnbc,vi_imaging_clinical,all=T)
#vi_Fig3<-merge(vi_imaging_fh,vi_tmp,all=T)
vi_Fig4<-rbind(vi_tnbc_clinical_fh,vi_imaging_tnbc_clinical,vi_imaging_tnbc_clinical_fh)

str3=paste("AUC (TNBC.Clinical.Family.history) =", sprintf('%.3f',AUC_EN_noimaging),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_noimaging)[1]),"-",sprintf('%.3f',as.vector(auc_ci_noimaging)[3]))
str1=paste("AUC (Imaging.TNBC.Clinical) =", sprintf('%.3f',AUC_EN_imaging_noFH),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_imaging_nofh)[1]),"-",sprintf('%.3f',as.vector(auc_ci_imaging_nofh)[3]))
str2=paste("AUC (Imaging.TNBC.Clinical.Family.history) =", sprintf('%.3f',AUC_EN_all),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_all)[1]),"-",sprintf('%.3f',as.vector(auc_ci_all)[3]))

lp1<-ggplot(df_all, aes(d = D, m = M, group=Data_Type,color=Data_Type,linetype=Data_Type)) + geom_roc(n.cuts=0,size=1.5)+style_roc(guide=F,xlab = "1 - Specificity",ylab="Sensitivity")+annotate("text", size=3, x = .25, y = .12,hjust=0,label=paste(str1,str2,str3,sep="\n"),color="black")
lp1<-lp1+scale_linetype_manual(values=c("solid", "twodash", "solid"))+scale_colour_Publication()+theme(legend.title = element_text(size=18, face="bold"),legend.text = element_text(size=18))+theme_Publication()

lp2<-ggplot(vi_Fig4, aes(reorder_within(Data_Type, Importance, Variable),y=Importance,fill=Data_Type,label=Variable))+coord_flip()+
  geom_col(position = position_dodge2(width = .5, preserve = "single"))+geom_text(position = position_dodge2(width = .5), angle = 0,hjust=-.1)+
  scale_x_reordered()+facet_grid(Data_Type~., scales = "free_y",space="free_y")+scale_fill_Publication()+theme_Publication()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(5,5,5,17),"mm"))+
  theme(strip.text.y = element_blank())+expand_limits(y = c(-1, 110))+guides(fill=guide_legend(nrow=1,byrow=TRUE))

lp_3_tmp<-ggarrange(lp1+rremove("legend"),lp2+rremove("y.text")+rremove("xylab")+theme(legend.position = "top",legend.title = element_blank()),labels = c("A","B"),nrow=2,ncol=1, align="h")

lp_3<-annotate_figure(lp_3_tmp, bottom = text_grob("Variable Importance (%)", color = "black", size = 12))

ggsave(file=paste("/KDI/pics_June2019/Q3/Fig4",".jpg",sep=""),plot=lp_3,width=210, height=297, units="mm",dpi=600)

# Question 4

# tnbc,clinical,fh1,imaging (no fh23)
# no imaging
lambda.seq<-10^seq(.1,-3,length=20)
alpha.seq<-seq(.1,1, length=20)
elgrid <- expand.grid(alpha=1,lambda=.1)
set.seed(849)
mod_fit_noFH23 <- train(Mutationstatus ~ ., data=subset(dat3[,c(2,I_TNBC,I_clinical,I_FamilyHistory_1,I_index)]),
                             method="glmnet",
                             metric="ROC",
                             family='binomial',
                             weights = model_weights,
                             tuneGrid=elgrid,
                             trControl=ctrl)

modcvpred_noFH23 <- predict(object=mod_fit_noFH23, dat3, type="prob")[2]
modcvpredt <- function(t) ifelse(modcvpred_noFH23  > t , 1,0)
pred <- prediction(modcvpred_noFH23, dat3$Mutationstatus)
perf <- performance(pred,"tpr","fpr")
plot(perf)
abline(a=0, b=1)
AUC_EN_noFH23 <- performance(pred,"auc")@y.values[[1]]
ss <- opt.cut(perf, pred)
print(confusionMatrix(as.factor(modcvpredt(ss[3])),as.factor(as.numeric(dat3$Mutationstatus)-1)))
print(AUC_EN_noFH23)
ncoefs_nofh23<-length(ncoefs(mod_fit_noFH23))-1
auc_ci_nofh23<-ci(dat3$Mutationstatus,modcvpred_noFH23$CONTROLS)
# All except imaging (done -mod_fit_noimaging)

# Output for Question 4

df_noFH23<-data.frame(ID=dat3$IncognitoPatientID,D = c(1,2)[dat3$Mutationstatus]-1, D.str = dat3$Mutationstatus,Imaging.TNBC.Clinical.Family.history.1 = modcvpred_noFH23$CONTROLS, stringsAsFactors = FALSE)

df_noimaging<-data.frame(TNBC.Clinical.Family.history.1.2.3 = modcvpred_noimaging$CONTROLS, stringsAsFactors = FALSE)

df=cbind(df_noFH23,df_noimaging)

df_all <- melt_roc(df, "D", c("Imaging.TNBC.Clinical.Family.history.1","TNBC.Clinical.Family.history.1.2.3"))

names(df_all)[3]<-"Data_Type"

vi_tnbc_clinical_fh<-data.frame(varImp(mod_fit_noimaging)[[1]])
vi_imaging_tnbc_clinical_fh1<-data.frame(varImp(mod_fit_noFH23)[[1]])

names(vi_tnbc_clinical_fh)<-"Importance"
names(vi_imaging_tnbc_clinical_fh1)<-"Importance"

vi_tnbc_clinical_fh$Variable<-as.factor(rownames(vi_tnbc_clinical_fh))
vi_tnbc_clinical_fh$Data_Type<-"TNBC.Clinical.Family.history.1.2.3"
vi_tnbc_clinical_fh$Variable <- factor(vi_tnbc_clinical_fh$Variable, levels = vi_tnbc_clinical_fh$Variable[order(vi_tnbc_clinical_fh$Importance,decreasing=T)])

vi_imaging_tnbc_clinical_fh1$Variable<-as.factor(rownames(vi_imaging_tnbc_clinical_fh1))
vi_imaging_tnbc_clinical_fh1$Data_Type<-"Imaging.TNBC.Clinical.Family.history.1"
vi_imaging_tnbc_clinical_fh1$Variable <- factor(vi_imaging_tnbc_clinical_fh1$Variable, levels = vi_imaging_tnbc_clinical_fh1$Variable[order(vi_imaging_tnbc_clinical_fh1$Importance,decreasing=T)])

#vi_tmp<-merge(vi_imaging_tnbc,vi_imaging_clinical,all=T)
#vi_Fig3<-merge(vi_imaging_fh,vi_tmp,all=T)
vi_Fig5<-rbind(vi_tnbc_clinical_fh,vi_imaging_tnbc_clinical_fh1)

#I_omitfh1<-df_all$Data_Type!="Imaging.TNBC.Clinical.FH_1"

str1=paste("AUC (Imaging.TNBC.Clinical.Family.history.1) =", sprintf('%.3f',AUC_EN_noFH23),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_nofh23)[1]),"-",sprintf('%.3f',as.vector(auc_ci_nofh23)[3]))
str2=paste("AUC (TNBC.Clinical..Family.history.1.2.3) =", sprintf('%.3f',AUC_EN_noimaging),"; 95% CI:", sprintf('%.3f',as.vector(auc_ci_noimaging)[1]),"-",sprintf('%.3f',as.vector(auc_ci_noimaging)[3]))

lp1<-ggplot(df_all, aes(d = D, m = M, group=Data_Type,color=Data_Type,linetype=Data_Type)) + geom_roc(n.cuts=0,size=1.5)+style_roc(guide=F,xlab = "1 - Specificity",ylab="Sensitivity")+annotate("text", size=3, x = .25, y = .12,hjust=0,label=paste(str1,str2,sep="\n"),color="black")
lp1<-lp1+scale_linetype_manual(values=c("solid", "solid"))+scale_colour_Publication()+theme(legend.title = element_text(size=18, face="bold"),legend.text = element_text(size=18))+theme_Publication()

lp2<-ggplot(vi_Fig5, aes(reorder_within(Data_Type, Importance, Variable),y=Importance,fill=Data_Type,label=Variable))+coord_flip()+
  geom_col(position = position_dodge2(width = .5, preserve = "single"))+geom_text(position = position_dodge2(width = .5), angle = 0,hjust=-.1)+
  scale_x_reordered()+facet_grid(Data_Type~., scales = "free_y",space="free_y")+scale_fill_Publication()+theme_Publication()+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(5,5,5,17),"mm"))+
  theme(strip.text.y = element_blank())+expand_limits(y = c(-1, 110))+guides(fill=guide_legend(nrow=1,byrow=TRUE))

lp_3_tmp<-ggarrange(lp1+rremove("legend"),lp2+rremove("y.text")+rremove("xylab")+theme(legend.position = "top",legend.title = element_blank()),labels = c("A","B"),nrow=2,ncol=1)

lp_3<-annotate_figure(lp_3_tmp, bottom = text_grob("Variable Importance (%)", color = "black", size = 12))

ggsave(file=paste("KDI/pics_June2019/Q4/Fig5",".jpg",sep=""),plot=lp_3,width=210, height=297, units="mm",dpi=600)

###Correlation Plot

dat<-subset(dat3[,c(2,I_TNBC,I_clinical,I_FamilyHistory_1,I_FamilyHistory_23,I_index)])

dat$MBC<-NULL

dat$Mutationstatus<-c(1,2)[dat$Mutationstatus]-1

names(dat)[names(dat)=="Mutationstatus"]="MS"

res1 <- cor.mtest(dat,0.99)

corr.A<-cor(dat)

#corrplot.mixed(corr.A, number.cex = 1, p.mat = res1$p, insig = "blank", order = "AOE",tl.pos = "l")
#corrplot(corr.A, order = "AOE",tl.pos = "l", col = brewer.pal(n = 8, name = "RdBu"))

jpeg(file=paste("KDI/pics_June2019/CorrPlot",".jpg",sep=""),width=6,height=6,units="in",res=600)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(corr.A, method = "color", col = col(200),
         type = "upper", order = "AOE", number.cex = .7,
         addCoef.col = "black",# Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

dev.off()




