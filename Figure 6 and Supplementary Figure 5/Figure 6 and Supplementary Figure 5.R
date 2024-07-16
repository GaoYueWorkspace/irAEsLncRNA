## Figure 6.A
if(T){
  ### 参数设置
  max=ceiling(max(da$pvalue))
  min=floor(min(da$pvalue))
  ### 
  ### 画图
  ggplot(data = da,
         aes(x=id, y=pvalue, group=type, fill=type, color=type)) +
    geom_hline(data = NULL, yintercept = c(0,10,15), color="grey85") +
    geom_hline(data = NULL, yintercept = 20, color="grey15") +
    geom_segment(data = tibble(x=1:length(colnames(pdata2)), y=0, yend=(20+1)),
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey85") +
    geom_segment(data = tibble(x=1:length(colnames(pdata2)), y=20, yend=(20+1)),
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey15") +
    geom_segment(data = NULL,
                 aes(x=1, xend=1, y=0, yend=20),
                 inherit.aes = FALSE, color="grey15") +
    geom_textpath(data = da %>%
                    select(cancer,id) %>% 
                    distinct(),
                  aes(x=id, y=(20+3), label=cancer),
                  inherit.aes = FALSE, vjust=1) +
    geom_point(data = da,
               aes(x=id, y=count.y, fill=type, size=Count),
               shape=21, alpha = 0.75)+
    ggalt::stat_xspline(geom = "line", spline_shape = 0.25, linewidth=0.75) +
    ggalt::stat_xspline(geom = "area", alpha=0.25, spline_shape = 0.25, outline.type="upper") +
    scale_size(range = c(0, 10), breaks = c(10, 25, 40), 
               guide = guide_legend(title = "Count")) +
    scale_fill_manual(values = c("#17A6B1","#DEB317")) +
    scale_color_manual(values = c("#17A6B1","#DEB317")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-10, 25), expand = c(0, 0)) +
    coord_polar() +
    annotate("text", x=1, y=c(10, 15), label="-", hjust=1, size=4) +
    annotate("text", x=1.1, y=c(10, 15), label=c(10, 15), hjust=0, size=4) +
    annotate("text", x=1, y=-5, label = str_c("-log10Pvalue"), size=4.5) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}
## Figure 6.B
if(T){
  multicox<-coxph(Surv(time = OS.time,event = OS) ~ .,data =forest.da)
  ggforest(model = multicox,
           data = forest.da,
           fontsize=0.7,
           noDigits=2,
           main = i)
}
## Figure 6.C 
if(T){
  survival.km <- survfit(Surv(time = OS.time,event = OS) ~ sur.value,data = cut.da)
  ggsurvplot(survival.km, data = cut.da,
               surv.median.line = "hv",  
               pval = TRUE,
               conf.int=TRUE,
               pval.size=5,
               ggtheme=theme_light(),
               xlab="Days",
               linetype = 1,
               palette=c("#36A2E0", "#E7B800"),
               legend.labs=c("High Risk", "Low Risk"))
}
## Figure 6.D
if(T){
  survival.km <- survfit(Surv(time = PFI_time,event = PFI_status) ~ sur.value,data = cut.da)
  pvalue.km=surv_pvalue(survival.km, method = "survdiff")$pval
  ggsurvplot(survival.km, data = cut.da,
               surv.median.line = "hv", 
               pval = TRUE,
               conf.int=TRUE,
               pval.size=5,
               ggtheme=theme_light(),
               xlab="Days",
               linetype = 1,
               palette=c("#36A2E0", "#E7B800"),
               legend.labs=c("High Risk", "Low Risk"))
}
## Figure 6.E
if(T){
  cox_train1 <- cph(Surv(OS.time,OS)~Gender+Age+IRaesLnc2Gene,
                    data=cox.df,surv=T,x=T, y=T,time.inc = 365*year)
  cal_train1<-tryCatch(calibrate(cox_train1,u=365*year,cmethod='KM',m=m,B=2000),
                       error=function(e){"error"})
  if(class(cal_train1)!="calibrate") next
  x=paste0('Nomogram-Predicted Probability of ',year,'-year')
  y=paste0(year,'-year survival probability')
  plot(cal_train1,lwd=2,lty=2, 
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
       xlab=x,#x轴名称
       ylab=y,#y轴名称
       col=c(rgb(192,98,83,maxColorValue = 255)))
  abline(a = 0, b = 1, col = "grey80", lwd = 2)
  lines(cal_train1[,c('mean.predicted',"KM")],lwd = 2, col = c("#eb4d4b"), pch = 16)
}
## Figure 6.F
if(T){
  set.seed(1)
  pvalue=wilcox.test(plot.da$value[plot.da$group=="Normal"],
                     plot.da$value[plot.da$group=="Tumor"])
  ggplot(data=plot.da,
         aes(x=V1,y=value,fill=group))+
    geom_boxplot(na.rm = T)+
    geom_jitter(aes(fill=group),
                position = position_jitterdodge(0.5),
                shape = 21,size=1) +
    geom_signif(annotations = pvalue$p.value, 
                # y_position = 1.5,
                y_position = 5,
                xmin = 0.785,
                xmax = 1.225,
                tip_length = 0.02, 
                color="black",
                vjust = 0) +
    theme_bw()+
    labs(x="",y="log2(Expression)",fill="TCGA-LUAD",title = "ALAL-1") +
    theme(axis.text.x = element_text(colour = "#000000",size = 11),
          axis.text.y = element_text(colour = "#000000",size = 11))
}
## Figure 6.G
if(T){
  x.pa=ceiling((max(df$x)-min(df$x))/2)
  y.pa=max(df$y)-2
  x.min=min(df$x)
  y.min=min(df$y)
  
  p1<-ggplot(df,aes(x,y))+
    geom_point(color="grey40",alpha=0.5)+
    theme_bw()+
    geom_smooth(method = "lm",formula = "y~x",
                se=F,color="blue",linewidth=2) +
    scale_y_continuous(expand = c(0.01,0.01))+
    scale_x_continuous(expand = c(0.01,0.01))+
    labs(x="ALAL-1 Expression",y="Erlotinib IC50 value")+
    stat_cor(method = "spearman")
  
  p2<-ggplot(df,aes(x))+
    geom_density(fill="#FED439",alpha=1,lwd=NA)+
    scale_y_continuous(expand = c(0,0))+
    labs(title = "")+
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  p3<-ggplot(df,aes(y))+
    geom_density(fill="#D5E4A2",alpha=1,lwd=NA)+
    scale_y_continuous(expand = c(0,0))+
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())+
    coord_flip()
  p = p1 %>%
    aplot::insert_top(p2,height = 0.2)%>%
    aplot::insert_right(p3, 0.2) 
}
## Figure 6.H
if(T){
  survival.km <- survfit(Surv(time = PFI_time,event = PFI_status) ~ drug,data = cut.da)
  ggsurvplot(survival.km, data = cut.da,
               surv.median.line = "hv",
               pval = TRUE,
               conf.int=TRUE,
               pval.size=5,
               ggtheme=theme_light(),
               xlab="PFI Days",
               linetype = 1,
               palette=c("#36A2E0", "#E7B800"),
               legend.labs=c("High Risk", "Low Risk"))
}
## Figure 6.I
if(T){
  ggplot()+
    geom_bar(data = plot.da,
             aes(x=gene,y=AUC,fill=method),
             stat = 'identity',
             position="dodge",
             width = 0.9) +
    scale_y_continuous(expand = c(0,0),limits = c(0,1))+
    scale_fill_manual(values = c("#B195BD","#DD8C30","#B95C61","#6BA8C2"))+
    geom_text(data = plot.da, 
              aes(x=gene,y=AUC,fill=method,label = round(AUC,2)),
              position = position_dodge(0.9),
              vjust = -0.5, hjust = 0.5,color="#000000") +
    theme_bw() +
    labs(x="Different machine learning methods",y="Area Under Curve (AUC)",fill="")+
    theme(axis.text.x = element_text(colour = "#000000",size = 11),
          axis.text.y = element_text(color = "#000000",size = 11),
          legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}
## Figure 6.J
if(T){
  survival.km <- survfit(Surv(time = PFS,event = PFS_event) ~ benefit,data = survival.da)
  ggsurvplot(survival.km, data = survival.da,
               surv.median.line = "hv",
               pval = TRUE,
               conf.int=TRUE,
               pval.size=5,
               ggtheme=theme_light(),
               xlab="FPS Days",
               linetype = 1,
               palette=c("#36A2E0", "#E7B800"),
               legend.labs=c("NR", "R"))
}
## Supplementary Figure 5.A
if(T){
  max=ceiling(max(da$pvalue))
  min=floor(min(da$pvalue))
  
  ### 画图
  ggplot(data = da,
         aes(x=id, y=pvalue, group=type, fill=type, color=type)) +
    geom_hline(data = NULL, yintercept = c(0,10,15), color="grey85") +
    geom_hline(data = NULL, yintercept = 18, color="grey15") +
    geom_segment(data = tibble(x=1:length(colnames(pdata2)), y=0, yend=(18+1)),
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey85") +
    geom_segment(data = tibble(x=1:length(colnames(pdata2)), y=18, yend=(18+1)),
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey15") +
    geom_segment(data = NULL,
                 aes(x=1, xend=1, y=0, yend=18),
                 inherit.aes = FALSE, color="grey15") +
    geom_textpath(data = da %>%
                    select(cancer,id) %>% 
                    distinct(),
                  aes(x=id, y=(18+3), label=cancer),
                  inherit.aes = FALSE, vjust=1,size=3.5) +
    geom_point(data = da,
               aes(x=id, y=count.y, fill=type, size=Count),
               shape=21, alpha = 0.75) +
    ggalt::stat_xspline(geom = "line", spline_shape = 0.25, linewidth=0.75) +
    ggalt::stat_xspline(geom = "area", alpha=0.25, spline_shape = 0.25, outline.type="upper") +
    scale_size(range = c(0, 10), breaks = c(10, 25, 40), 
               guide = guide_legend(title = "Count")) +
    scale_fill_manual(values = c("#17A6B1","#DEB317")) +
    scale_color_manual(values = c("#17A6B1","#DEB317")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-10, 25), expand = c(0, 0)) +
    coord_polar() +
    annotate("text", x=1, y=c(7,10,max), label="-", hjust=1, size=4) +
    annotate("text", x=1.1, y=c(7,10,max), label=c(7,10,max), hjust=0, size=4) +
    annotate("text", x=1, y=-5, label = str_c("-log10Pvalue"), size=4.5) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}
## Supplementary Figure 5.B
if(T){
  dd<-datadist(cox.df)
  options(datadist="dd")
  cox <- cph(Surv(OS.time,OS)~Gender+Age+IRaesLnc2Gene,
             surv=T,x=T, y=T,data=cox.df)
  surv <- Survival(cox)
  sur_1_year<-function(x)surv(1*365*1,lp=x)
  sur_3_year<-function(x)surv(1*365*3,lp=x)
  sur_5_year<-function(x)surv(1*365*5,lp=x)
  condition=if_else(max(cox.df$OS.time)/(365*5)>=1,3,
                    if_else(max(cox.df$OS.time)/(365*3)>=1,2,
                            if_else(max(cox.df$OS.time)/(365*1)>=1,1,NA)))
  canshu=nomogram.input(condition)
  nom_sur <- nomogram(cox,fun=canshu$fun.list,
                      lp= F,
                      funlabel=canshu$label.c,
                      maxscale=100,
                      fun.at=c('0.9','0.7','0.5','0.3','0.1'))
  plot(nom_sur,xfrac=0.25)
}
## Supplementary Figure 5.C
if(T){
  dd<-datadist(cox.df)
  options(datadist="dd")
  cox <- cph(Surv(OS.time,OS)~Gender+Age+IRaesLnc2Gene,
             surv=T,x=T, y=T,data=cox.df)
  surv <- Survival(cox)
  sur_1_year<-function(x)surv(1*365*1,lp=x)
  sur_3_year<-function(x)surv(1*365*3,lp=x)
  sur_5_year<-function(x)surv(1*365*5,lp=x)
  condition=if_else(max(cox.df$OS.time)/(365*5)>=1,3,
                    if_else(max(cox.df$OS.time)/(365*3)>=1,2,
                            if_else(max(cox.df$OS.time)/(365*1)>=1,1,NA)))
  canshu=nomogram.input(condition)
  nom_sur <- nomogram(cox,fun=canshu$fun.list,
                      lp= F,
                      funlabel=canshu$label.c,
                      maxscale=100,
                      fun.at=c('0.9','0.7','0.5','0.3','0.1'))
  plot(nom_sur,xfrac=0.25)
}
## Supplementary Figure 5.D
if(T){
  cox_train1 <- cph(Surv(OS.time,OS)~Gender+Age+IRaesLnc2Gene,
                    data=cox.df,surv=T,x=T, y=T,time.inc = 365*year)
  cal_train1<-tryCatch(calibrate(cox_train1,u=365*year,cmethod='KM',m=m,B=2000),
                       error=function(e){"error"})
  if(class(cal_train1)!="calibrate") next
  x=paste0('Nomogram-Predicted Probability of ',year,'-year')
  y=paste0(year,'-year survival probability')
  plot(cal_train1,lwd=2,lty=2, 
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
       xlab=x,#x轴名称
       ylab=y,#y轴名称
       col=c(rgb(192,98,83,maxColorValue = 255)))
  abline(a = 0, b = 1, col = "grey80", lwd = 2)
  lines(cal_train1[,c('mean.predicted',"KM")],lwd = 2, col = c("#eb4d4b"), pch = 16)
}
## Supplementary Figure 5.E
if(T){
  ggplot(result1, aes(gene_name1, gene_name2, fill = R)) +
    geom_tile(color = 'white', size = 2) +
    geom_text(aes(label = p.sig), color="black",size = 5) +
    scale_fill_gradient2(low ="blue" , high = "red", mid = "grey90", midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "#000000",size = 12),
          axis.text.y = element_text(colour = "#000000",size = 12)) +
    labs(title = "",x = "",y = "") -> p
}
## Supplementary Figure 5.F
if(T){
  load("Supplementary Figure 5.F.Rdata")
  p=plot_dim(obj)
  ggsave("EGFR.scRNA-seq.pdf",plot=p,width = 6,height = 5)
}
## Supplementary Figure 5.G
if(T){
  set.seed(1)
  pvalue=wilcox.test(plot.da$value[plot.da$group=="Normal"],
                     plot.da$value[plot.da$group=="Tumor"])
  ggplot(data=plot.da,
         aes(x=V1,y=value,fill=group))+
    geom_boxplot(na.rm = T)+
    geom_jitter(aes(fill=group),
                position = position_jitterdodge(0.5),
                shape = 21,size=1) +
    geom_signif(annotations = pvalue$p.value, 
                # y_position = 1.5,
                y_position = 5,
                xmin = 0.785,
                xmax = 1.225,
                tip_length = 0.02, 
                color="black",
                vjust = 0) +
    theme_bw()+
    labs(x="",y="log2(Expression)",fill="TCGA-LUAD",title = "SEC61G-DT") +
    theme(axis.text.x = element_text(colour = "#000000",size = 11),
          axis.text.y = element_text(colour = "#000000",size = 11))
}
## Supplementary Figure 5.H
if(T){
  x.pa=ceiling((max(df$x)-min(df$x))/2)
  y.pa=max(df$y)-2
  x.min=min(df$x)
  y.min=min(df$y)
  
  p1<-ggplot(df,aes(x,y))+
    geom_point(color="grey40",alpha=0.5)+
    theme_bw()+
    geom_smooth(method = "lm",formula = "y~x",
                se=F,color="blue",linewidth=2) +
    scale_y_continuous(expand = c(0.01,0.01))+
    scale_x_continuous(expand = c(0.01,0.01))+
    labs(x="SEC61G-DT Expression",y="Erlotinib IC50 value")+
    stat_cor(method = "spearman")
  
  p2<-ggplot(df,aes(x))+
    geom_density(fill="#FED439",alpha=1,lwd=NA)+
    scale_y_continuous(expand = c(0,0))+
    labs(title = "")+
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  p3<-ggplot(df,aes(y))+
    geom_density(fill="#D5E4A2",alpha=1,lwd=NA)+
    scale_y_continuous(expand = c(0,0))+
    theme_minimal()+
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())+
    coord_flip()
  p = p1 %>%
    aplot::insert_top(p2,height = 0.2)%>%
    aplot::insert_right(p3, 0.2)
}
## Supplementary Figure 5.I
if(T){
  ggplot(data = plot.da,
         aes(gene, AUC, color = method, group = method))+
    geom_line(linewidth =1)+ ## 线条粗细
    geom_point(size=3)+
    labs(title = "",x = "",y = "")+
    scale_color_manual(values = c("#d75b44","#ff7500","#549688","#1685a9"))+
    theme_bw()+
    theme(axis.text.x = element_text(color = "#000000"),
          axis.text.y = element_text(color = "#000000"))
}