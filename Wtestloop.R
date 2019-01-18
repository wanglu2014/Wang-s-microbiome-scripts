library(reshape)
library(ggplot2)
library(gtools)
#template Groups:string cont_case:control,case
Wtest<-function(abund_table){
  name<-unique(abund_table$Groups)
  groups<-as.factor(abund_table$cont_case)
  data =as.data.frame(Filter(is.numeric, abund_table))
  W.alpha=0.01
  W.table <- data.frame()
  for (i in 1:dim(data)[2]) {
    
    W.test <- wilcox.test(data[,i]~groups)
    W.fold <- foldchange(mean(data[,i][groups=='case'],na.rm = F),mean(data[,i][groups=='control'],na.rm = F))
    W.table <- rbind(W.table,data.frame(id=names(data)[i],p.value=W.test$p.value,fold=W.fold))
  }
  
  
  W.table$adjvalue <- p.adjust(W.table $p.value,method = "fdr")
  W.table$FWER <- pbinom(q=0, p=W.table$p.value,size=dim(W.table)[1], lower.tail=FALSE)
  W.table <- W.table[order(W.table$p.value,decreasing=FALSE), ]
  W.table$q.value.factor <- dim(W.table)[1] / 1:dim(W.table)[1]
  W.table$q.value <- W.table$p.value * W.table$q.value.factor
  write.csv(W.table,paste0(name,'foldW.csv'))
}

abund_table<-read.csv("G:/CIcopy - Copy/CI/basic source/gene/0712/wuSSN9_sam01prop.csv",row.names=1,check.names=FALSE)

grouping_info<-read.csv("G:/CIcopy - Copy/CI/basic source/gene/0712/wuSSN9contro_design.csv",row.names=1,check.names=FALSE)

caseabu<-abund_table%>%
  rownames_to_column(var='Samples')%>%
  inner_join(grouping_info%>%rownames_to_column(var='Samples'))%>%
  column_to_rownames(var='Samples')%>%
  mutate(Groups=paste0(Organ,'-',`Age (group)`))%>%#craft Group label
  group_by(Groups)%>%
  do(Wtest(.))