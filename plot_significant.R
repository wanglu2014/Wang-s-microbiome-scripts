#install.packages("ggsignif")
#How to add significant star
library(ggsignif)
attr_meta<- read.csv("~/netattrmeta.csv")
compared_list<-list(c("Chongming Nonsmoker","Shanghai Smoker"))
ggplot(data = attr_meta, mapping = aes(x=Subgroup,y =mean.degree)) +
  geom_violin(alpha = 0) +
  stat_summary(fun.y=mean, colour="darkred", geom="point", hape=18, size=3,show_guide = FALSE)+
  geom_jitter(alpha = 0.8, color = "tomato")+
  geom_signif(comparisons = compared_list,test=t.test,map_signif_level = TRUE, textsize=3, step_increase = 0.2)
