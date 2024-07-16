## Figure 3.A
if(T){
  export(df,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Figure 3.A.xlsx")
  
  ggplot(data=df,aes(x=Var1,y=Freq,fill=Type)) + 
    geom_bar(stat="identity",position="fill",width = 0.7) +
    scale_fill_manual(values=c("#CC0000","#D4E6AF","#25ACAE"))+
    scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                       labels = scales::percent_format())+
    labs(x="",y="Proportion",fill=" ",title="") +
    theme_bw()+
    theme(axis.title.y=element_text(size=14))+
    theme(legend.text=element_text())+
    theme(axis.text.x = element_text(size = 12, color = "black"))+
    theme(axis.text.y = element_text(size = 12, color = "black"))+
    theme(axis.ticks.length=unit(0.3,"cm"))+
    theme(axis.text.x=element_text(angle=90,size = 11))
  ggsave("D:\\DLL\\TCGA\\Part3\\0.3.1.hub and specific IRaeslnc\\hub and specific AND nohubnospecific.pdf"
         ,width = 6.5,height = 3.5)
}
## Figure 3.B
if(T){
  p=ggplot(dff,aes(x=cancer,y=median,color=type))+
    geom_errorbar(data = dff,aes(ymin = qu2, ymax=qu4,group=type), 
                  width=0.6, 
                  position=position_dodge(1), 
                  alpha = 0.7,
                  size=1)+
    geom_point(data = dff,aes(x=cancer, y=median),
               pch=19,
               position=position_dodge(1),
               size=3.4,
               alpha=1) +
    geom_signif(annotations = p.sig,
                y_position = 7,
                xmin = seq(0.75,15.75,1),
                xmax = seq(1.25,16.25,1),
                tip_length = rep(0.01,16),
                color="black") +
    geom_vline(xintercept=c(seq(1.5,15.5,by=1)), linetype="dotted")+
    labs(x="Cancers",y="Log2(Expression+0.01)")+
    theme_bw()+
    scale_color_manual(values=c("#CC0000","#53aa73","#25ACAE"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,color = "#000000"),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(color = "#000000"),
          legend.position="none",
          panel.grid.minor = element_blank()) 
  
  ggsave("D:\\DLL\\TCGA\\Part3\\0.3.2.difference of hub and specific\\标准化误差线图.pdf",
         plot = p,width = 6.5,height = 4)
}
## Figure 3.C
if(T){
  p=ggplot()+geom_segment(data =df,aes(x = zuo.x,
                                       y = zuo.y,
                                       xend = zhong.x, 
                                       yend = zhong.y),
                          size=df$segment.size,
                          color=df$segment)
  export(list(df=df,
              zuo.df=zuo.df,
              zhong.df=zhong.df,
              you.df=you.df),"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Figure 3.C.xlsx")
  
  ######  合并图层
  pp = p+geom_segment(data = df, aes(x = zhong.x, 
                                     y = zhong.y, 
                                     xend = you.x, 
                                     yend = you.y),
                      color=df$segment,
                      size=df$segment.size) +
    ## 左边点和标签
    geom_text(data = zuo.df,
              mapping = aes(x=zuo.x-0.5,y=zuo.y,label=cancer),
              hjust=0)+
    geom_point(data = zuo.df,
               mapping = aes(x=zuo.x,y=zuo.y),
               shape=16,
               size=zuo$zuo.size,
               color=zuo.color) +
    ## 中间点和标签
    geom_text(data = zhong.df,
              mapping = aes(x=zhong.x-0.1,y=zhong.y-0.5,label=da.Type),
              hjust=0)+
    geom_point(data = zhong.df,
               mapping = aes(x=zhong.x,y=zhong.y),
               shape=16,
               size=setNames(zhong$zhong.size,zhong$id),
               color="#FFFFFF") +
    geom_point(data = zhong.df,
               mapping = aes(x=zhong.x,y=zhong.y),
               shape=1,
               size=setNames((zhong$zhong.size),zhong$id),
               color=zhong.color,
               stroke =2) +
    ## 右边点和标签
    geom_text(data = you.df,
              mapping = aes(x=you.x+0.2,y=you.y,label=gene_name),
              hjust=0) +
    geom_point(data = you.df,
               mapping = aes(x=you.x,y=you.y),
               shape=16,
               size=setNames(you$you.size,you$id),
               color=you.color) +
    scale_x_continuous(limits = c(-0.6,5)) +
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.position="none")
}
## Figure 3.D
if(T){
  export(skcm.iraes.lnc.score,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Figure 3.D.xlsx")
  skcm.iraes.lnc.score %>%
    tidyr::pivot_longer(cols = -my.project,
                        names_to = "fun",
                        values_to = "value") %>%
    filter(my.project!="other")-> plot.da
  plot.da$my.project=factor(plot.da$my.project,levels = c("DERM","ENDO","GI","MUSCULO","NEURO","No_iraes","Control"))
  ggplot(data=plot.da %>%
           filter(fun==i),
         aes(x=my.project,y=value,fill=my.project))+
    scale_fill_manual(values = setNames(c("#CC789F","#A4C7DF","#E0A8A4","#CCBFD4","#3A9EC9","#f78945","#0B775E"),
                                        c("DERM","ENDO","GI","MUSCULO","NEURO","No_iraes","Control")))+
    stat_boxplot(geom = "errorbar",width=0.1,size=0.2)+theme_bw()+
    geom_boxplot(width=0.65,notch =F,outlier.fill=NA,outlier.color=NA,size=0.2)+
    geom_jitter(width =0.2,shape = 21,size=5,alpha=0.8,stroke=0.1)+
    labs(x="",y="Score of Colitis-related irAEsLnc")+
    ylim(-0.45,0.55)+
    theme(panel.grid=element_blank())+
    theme(legend.text = element_text(colour = "black",size =12))+
    theme(legend.position = "right")+
    theme(axis.text.x.top = element_blank())+
    theme(axis.ticks.y.right=element_blank())+
    theme(axis.text.y.right = element_blank())+
    theme(axis.ticks.x.top=element_blank())+
    theme(axis.title = element_text(size =14,colour = "black"))+
    theme(axis.text = element_text(size = 12,colour = "black"))+
    theme(axis.ticks = element_line(colour = "black"))+
    theme(axis.line = element_line(colour="black",size=0.4))+
    theme(axis.line.x.top = element_line(color = "black",size=0.4))+
    theme(axis.line.y.right=element_line(color = "black",size=0.4))+
    guides(x.sec="axis",y.sec = "axis")+
    theme(legend.position = "none")
}
## Figure 3.E
if(T){
  df %>%
    tidyr::pivot_wider(id_cols = -Type,
                       names_from = "cancer",
                       values_from = "lncrna.mean") %>%
    export("Figure 3.E.xlsx")
  lie=distinct(select(df,Type,iraes))
  dff=dff[c(lie$iraes[lie$Type=="PULM"],
            lie$iraes[lie$Type=="ENDO"],
            lie$iraes[lie$Type=="MUSCULO"],
            lie$iraes[lie$Type=="DERM"],
            lie$iraes[lie$Type=="NEURO"],
            lie$iraes[lie$Type=="GI"],
            lie$iraes[lie$Type=="RENAL"]),]
  ann_row=rbind(lie[lie$Type=="PULM",],
                lie[lie$Type=="ENDO",],
                lie[lie$Type=="MUSCULO",],
                lie[lie$Type=="DERM",],
                lie[lie$Type=="NEURO",],
                lie[lie$Type=="GI",],
                lie[lie$Type=="RENAL",])%>%
    as.data.frame() %>%
    tibble::column_to_rownames("iraes")
  ann_cols=list(Type=c("PULM"=adjustcolor("#AA709A",0.7),
                       "ENDO"=adjustcolor("#A4C7DF",0.7),
                       "MUSCULO"=adjustcolor("#CCBFD4",0.7),
                       "DERM"=adjustcolor("#CC789F",0.7),
                       "NEURO"=adjustcolor("#3A9EC9",0.7),
                       "GI"=adjustcolor("#E0A8A4",0.7),
                       "RENAL"=adjustcolor("#80699F",0.7)))
  pheatmap::pheatmap(as.matrix(dff),
                       cluster_cols = F,
                       cluster_rows = F,
                       na_col = "grey85",
                       border_color="#ffffff",
                       cellwidth=12,
                       cellheight = 12,
                       legend = T,
                       gaps_row = c(1,5,7,9,15,18),
                       annotation_row = ann_row,
                       annotation_colors = ann_cols,
                       color = colorRampPalette((brewer.pal(n = 9, name ="YlGnBu")))(100)[20:60])
}
## Figure 3.F
if(T){
  col1 <- grDevices::colorRampPalette(c("#440154FF","#007896FF","#3AC96DFF","#FDE725FF"))
  #####  画图
  pdf("D:\\DLL\\TCGA\\Part3\\0.3.3.share lncrna\\共享LNCRNA.pdf")
  corrplot(as.matrix(df),
           is.corr=FALSE,
           addCoef.col="grey30",
           col.lim = c(0, 240),
           type = "lower",
           method = "square",
           tl.col = "black",
           tl.cex = 1,
           number.cex = 0.8,
           cl.pos = "r",
           # outline = "#FFFFFF",
           col = rev(col1(200)))
  dev.off()
}