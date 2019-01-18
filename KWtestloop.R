library(reshape)
library(ggplot2)
library(gtools)
#template Groups:string cont_case:control,case
KWtest<-function(abund_table){
	name<-unique(abund_table$Groups)
	groups<-as.factor(abund_table$cont_case)
	data =as.data.frame(Filter(is.numeric, abund_table))
	kruskal.wallis.alpha=0.01
	kruskal.wallis.table <- data.frame()
	for (i in 1:dim(data)[2]) {
		ks.test <- kruskal.test(data[,i], g=groups)
		ks.fold <- foldchange(mean(data[,i][groups=='case'],na.rm = F),mean(data[,i][groups=='control'],na.rm = F))
		kruskal.wallis.table <- rbind(kruskal.wallis.table,
																	data.frame(id=names(data)[i],
																						 p.value=ks.test$p.value,
																						 fold=ks.fold))
	}
	
	
	kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
	kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
																			size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
	kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
																										 decreasing=FALSE), ]
	kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
	kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
	write.csv(kruskal.wallis.table,paste0(name,'foldKW.csv'))
}

abund_table<-read.csv("G:/CIcopy - Copy/CI/basic source/gene/0712/wuSSN9_sam01prop.csv",row.names=1,check.names=FALSE)

grouping_info<-read.csv("G:/CIcopy - Copy/CI/basic source/gene/0712/wuSSN9contro_design.csv",row.names=1,check.names=FALSE)

caseabu<-abund_table%>%
	rownames_to_column(var='Samples')%>%
	inner_join(grouping_info%>%rownames_to_column(var='Samples'))%>%
	column_to_rownames(var='Samples')%>%
	mutate(Groups=paste0(Organ,'-',`Age (group)`))%>%#craft Group label
	group_by(Groups)%>%
	do(KWtest(.))