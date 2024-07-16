## Figure 4.A
if(T){
  name=c("YlOrRd","YlGnBu","YlGn","Reds","RdPu","Purples","PuBu","OrRd","Oranges","Greys","Greens","GnBu","BuPu","Blues")
  cols=sapply(name,function(x){
    colorRampPalette(brewer.pal(3, x))(30)
  })
  
  col2 <- colorRampPalette(c("#F4A582","#FDDBC7", "#D1E5F0",
                             "#92C5DE","#4393C3", "#2166AC", "#053061"))
  pdf("0.4.1mianyijinrun/cibersort_xcell_timer.pdf",height = 4)
  corrplot(as.matrix(t(df)),
           is.corr=FALSE,
           method="circle",
           col=cols[,2],
           col.lim = c(0,1),
           bg=brewer.pal(9,"Pastel1")[9],
           outline=F,
           cl.pos='r',
           na.label = "X",
           na.label.col = "black",
           cl.length =5,
           cl.cex = 1,
           cl.ratio = 0.5 ,
           cl.align = "r",
           tl.cex=1,
           tl.col='black')
  dev.off()
}
## Figure 4.B
if(T){
  #####  文本信息
  cancer=data.frame(cancer=unique(da$cancer),
                    cancer_y=seq(2,(2*length(unique(da$cancer))),by=2),check.names = F)
  type=data.frame(type=unique(da$cell),
                  type_x=seq(2,(2*length(unique(da$cell))),by=2),check.names = F)
  da=da %>%
    inner_join(cancer,by=c("cancer"="cancer")) %>%
    inner_join(type,by=c("cell"="type")) %>%
    select(cancer,cancer_y,cell,type_x,count,everything())
  
  ######  画图信息
  ggplot() +
    scatterpie::geom_scatterpie(data = da,aes(x=type_x, y=cancer_y,r=count/4),  
                                cols=colnames(da)[-(1:5)],color=NA) +
    scale_fill_manual(values = setNames(c("#5D9D52","#FF7F00","#77468E","#46ACC8",
                                          "#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF"),
                                        colnames(da)[-c(1:5)])) +
    scale_y_continuous(breaks = as.numeric(unique(da$cancer_y)),
                       labels = unique(da$cancer)) +
    scale_x_continuous(breaks = as.numeric(unique(da$type_x)),
                       labels = unique(da$cell)) +
    theme_bw()+
    theme(axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,color = "#000000"),
          axis.text.y = element_text(color = "#000000"),
          legend.title = element_text(),
          panel.background = element_blank(),
          panel.grid=element_blank()) 
}
## Figure 4.C && Supplementary Figure 3.B
if(T){
  
  
  p=ggplot(da,aes(x=cancer,y=value,fill=type))+
    geom_boxplot(width=0.7,
                 alpha=0.8,
                 outlier.shape=NA,
                 na.rm = T)
  plot.data=ggplot_build(p)$data[[1]]
  signif.y.position.max=max(plot.data$ymax)+0.3
  signif.y.position.min=min(plot.data$ymin)-0.3
  
  p=p + stat_compare_means(aes(group=type),
                           label = "p.signif", 
                           method = "t.test",
                           label.y=signif.y.position.max,
                           color="black") + 
    geom_signif(annotations = rep("",18),
                y_position = rep(signif.y.position.max-0.2,18),
                xmin = seq(0.75,17.75,1),
                xmax = seq(1.25,18.25,1),
                tip_length = rep(0.0055,18),
                color="black")+
    geom_vline(xintercept=c(seq(1.5,17.5,by=1)), linetype="dotted") +
    scale_y_continuous(limits = c(round(signif.y.position.min),
                                  signif.y.position.max+1),
                       expand = c(0,0))+
    scale_fill_manual(values = c(adjustcolor("#Fe4f44",alpha.f = 0.5),
                                 adjustcolor("#F6b411",alpha.f = 0.5)))+
    labs(x="",y=paste0(legend," Log2(Expression+0.01)"),fill=legend)+
    guides(colour=guide_legend(title = legend))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,colour = "#000000"),
          axis.text.y = element_text(color = "#000000"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}
## Figure 4.D
if(T){
  
  name=c("YlOrRd","YlGnBu","YlGn","Reds","RdPu","Purples","PuBu",
         "OrRd","Oranges","Greys","Greens","GnBu","BuPu","Blues")
  cols=sapply(name,function(x){ colorRampPalette(brewer.pal(3, x))(30)})
  
  type=c("CYT","MHC","TCR","IFNG","TNF","TMB")      ###  添加新的
  mm=lapply(1:length(type), function(m){
    # m=1
    print(m)
    col=cols[,m]
    m.da=-log10(xx %>% filter(type==type[m]) %>% select(-type) %>% 
                  tibble::column_to_rownames("cancer")) %>% 
      as.data.frame()
    i.da=select(m.da,c("HUB","SPECIFIC","NOHUB_SPECIFIC","ORGAN"))
    if(m==1) axis=TRUE else axis=FALSE
    
    ## 列注释
    colann.da=colann[[type[m]]] %>% select(n.p,all_n.p) %>% as.data.frame()
    rownames(colann.da)=c("HUB","SPECIFIC","NOHUB_SPECIFIC","ORGAN")
    ha <- HeatmapAnnotation(
      Percent = anno_barplot(as.matrix(colann.da),
                             gp = gpar(fill = c(adjustcolor("#F7ED11",alpha.f = 0.7),adjustcolor("#5B90B8",alpha.f = 0.7)),
                                       col=c(adjustcolor("#F7ED11",alpha.f = 0.7),adjustcolor("#5B90B8",alpha.f = 0.7))),
                             axis = axis,
                             border = TRUE,
                             bar_width = 0.8),
      show_annotation_name = FALSE,
      height=unit(1, "cm"))
    ## 拼图
    p1 = Heatmap(as.matrix(i.da),
                 cluster_columns = F,
                 cluster_rows = F,
                 na_col = "grey90",
                 col=col,
                 heatmap_legend_param = list(
                   title = type[m]
                 ),
                 rect_gp = gpar(col = "white", lty = 1, lwd = 2),
                 column_names_side = "bottom",
                 row_names_side = "left",
                 show_row_names = F,
                 show_heatmap_legend = TRUE,
                 top_annotation = ha,
                 heatmap_width = unit(2.5, "cm"), 
                 heatmap_height = unit(17, "cm"))
    return(p1)
  })
  names(mm)=type
  p=mm$CYT+mm$MHC+mm$TCR+mm$IFNG+mm$TNF+mm$TMB        ###  添加新的
  ha1=rowAnnotation(foo = anno_text(sort(rownames(rowann))))
  ha2=rowAnnotation(foo = anno_boxplot(-log2(as.matrix(rowann)),
                                       height = unit(3, "cm"),
                                       outline = F,
                                       gp = gpar(fill = rep("#00CED1",18))))
  pdf("D:\\DLL\\TCGA\\Part4\\0.4.4.ICI biomarker/ICI biomarker.pdf",width = 10,height =15)
  draw(ha1+p+ha2,ht_gap = unit(0.1, "cm"))
  dev.off()
}
## Figure 4.E
if(T){
  p=pheatmap::pheatmap(as.matrix(dff),
                       cluster_cols = F,
                       cluster_rows = F,
                       na_col = "grey90",
                       border_color="#ffffff",
                       cellwidth=12,
                       cellheight = 12,
                       legend = TRUE,
                       #annotation_row = ann_row,
                       color = colorRampPalette((brewer.pal(n = 9, name ="BuGn")))(70)[30:70])
}
## Figure 4.F && Supplementary Figure 3.E
if(T){
  
  p=ggplot(xx,aes(x=celline,y=value,fill=type))+
    geom_boxplot(width=0.6,
                 alpha=0.8,
                 outlier.shape=NA)
  
  plot.data=ggplot_build(p)$data[[1]]
  signif.y.position=max(plot.data$ymax)+0.4
  p = p + stat_compare_means(aes(group=type),
                             label = "p.signif", 
                             method = "t.test",
                             label.y=signif.y.position,
                             color="black") + 
    geom_signif(annotations = rep("",19),
                y_position = rep(signif.y.position-0.2,19),
                xmin = seq(0.75,18.75,1),
                xmax = seq(1.25,19.25,1),
                tip_length = rep(0.0055,19),
                color="black")+
    geom_vline(xintercept=c(seq(1.5,18.5,by=1)), linetype="dotted")+
    scale_y_continuous(limits = c(0,signif.y.position+1.5),expand = c(0,0))+
    scale_fill_manual(values = c(adjustcolor("#Fe4f44",alpha.f = 0.5),
                                 adjustcolor("#F6b411",alpha.f = 0.5)))+
    labs(x="",y=i,fill="Type")+
    guides(colour=guide_legend(title = "")) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,colour = "#000000"),
          axis.text.y = element_text(colour = "#000000"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="top") 
}
## Figure 4.G
if(T){
 
  set.seed(2)
  max=50
  ggplot() +
    geom_hline(data = NULL, yintercept = c(0,max), color="grey85") +
    geom_hline(data = NULL, yintercept = max, color="grey15") +
    geom_segment(data = tibble(x=seq(0.5,72.5,4)), y=max, yend=(max+0.25),
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey15") +
    geom_segment(data = tibble(x=seq(0.5,72.5,4)), y=0, yend=(max),
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey85")+
    geom_segment(data = NULL,
                 aes(x=0, xend=0, y=0, yend=max),
                 inherit.aes = FALSE, color="grey15")+
    geom_textpath(data = da2 %>% select(Var1,id) %>%  distinct(),
                  aes(x=id, y=(max+3), label=Var1),
                  inherit.aes = FALSE, vjust=1) +
    geom_col(data=da2,
             aes(x = id,y = Freq,fill = type),
             width = 0.8,alpha = 1,
             stat="fill",position = position_stack()) +
    geom_col(data=da2,
             aes(x = id,y = -1,fill = Var2),
             width = 1) +
    geom_point(data = da1,
               aes(x=id, y=y.id, fill=Var2, size=Freq),
               shape=21, alpha = 1) +
    scale_size(range = c(0, 10), breaks = c(0, 100, 200,300,400,500), 
               guide = guide_legend(title = "IRaesLnc_Count")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.5, max+4), expand = c(0, 0)) +
    scale_fill_brewer(palette = 'Paired') +
    theme_bw()+ xlab("")+ ylab("Log2Count")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid=element_blank(),
          legend.position="top")
}
## Figure 4.H
if(T){
  cancer1=cancer[1]
  ggplot()+
    geom_segment(data = tibble(x=c(-0.3,0,0.3)), y=0, yend=10,
                 aes(x=x, xend=x, y=y, yend=yend),
                 inherit.aes = FALSE, color="grey85")+
    geom_density_ridges(data=da %>% filter(type=="MYLNC") %>% filter(file==cancer1),
                        aes(x=Rvalue,y=Ptype,fill=Ptype),
                        alpha=0.75,color=NA)+
    geom_density_ridges(data=da %>% filter(type=="OTHER") %>% filter(file==cancer1),
                        aes(x=Rvalue,y=Ptype),fill="grey50",
                        alpha=0.3,color=adjustcolor("grey50",alpha.f = 0.3)) +
    scale_x_continuous(expand = c(0,0))+
    theme(legend.position = "none")+
    labs(title = cancer1)+
    theme_bw() +
    labs(y="",x="")+
    scale_fill_manual(values = setNames(c("#ffc93c","#ff9a3c","#ff6f3c","#155263"),
                                        rev(c("PDL1","PD1","LAG3","CTLA4"))))+
    theme(axis.ticks.length.x = unit(0, "mm"),
          panel.grid=element_blank(),
          # plot.margin= margin(),
          legend.position = "None",
          plot.title = element_text(hjust = 0.5)) +
    scale_y_discrete(expand = c(0.1, 0)) -> p.cancer1
  
  for(h in cancer[2:18]){
    # h=cancer[2]
    ggplot()+
      geom_segment(data = tibble(x=c(-0.3,0,0.3)), y=0, yend=10,
                   aes(x=x, xend=x, y=y, yend=yend),
                   inherit.aes = FALSE, color="grey85")+
      geom_density_ridges(data=da %>% filter(type=="MYLNC") %>% filter(file==h),
                          aes(x=Rvalue,y=Ptype,fill=Ptype),
                          alpha=0.75,color=NA)+
      geom_density_ridges(data=da %>% filter(type=="OTHER") %>% filter(file==h),
                          aes(x=Rvalue,y=Ptype),fill="grey50",
                          alpha=0.3,color=adjustcolor("grey50",alpha.f = 0.3)) +
      scale_x_continuous(expand = c(0,0))+
      theme(legend.position = "none")+
      labs(title = h)+
      theme_bw() +
      labs(y="",x="")+
      scale_fill_manual(values = setNames(c("#ffc93c","#ff9a3c","#ff6f3c","#155263"),
                                          rev(c("PDL1","PD1","LAG3","CTLA4"))))+
      theme(axis.ticks.length.x = unit(0, "mm"),
            panel.grid=element_blank(),
            # plot.margin= margin(),
            legend.position = "None",
            plot.title = element_text(hjust = 0.5)) +
      scale_y_discrete(expand = c(0.1, 0)) -> p.cancer
    p.cancer1=p.cancer1+p.cancer
  }
}
## Supplementary Figure 3.A
if(T){
  
  df1<-data.frame(x = seq(1.5,length(unique(df$cancer))-0.5,1),
                  xend = seq(1.5,length(unique(df$cancer))-0.5,1),
                  y = -Inf,
                  yend = Inf)
  df2<-data.frame(y = seq(1.5,length(unique(df$name))-0.5,1),
                  yend = seq(1.5,length(unique(df$name))-0.5,1),
                  x = -Inf,
                  xend = Inf)
  df3=mutate(df,count=if_else(is.na(value),NA,count+10))
  ggplot() +
    geom_point(data=df,
               aes(x=cancer,y=name,size = count,color=value),
               shape = 16,alpha=1)+
    geom_point(data=df3,
               aes(x=cancer,y=name,size = count),
               shape = 1,alpha=1,color="black")+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color="black"),
          axis.ticks = element_blank(),
          axis.text.x=element_text(angle = 90,vjust = 0.5,hjust = 1))+
    geom_segment(data=df1,aes(x=x,xend=xend,y=y,yend=yend),
                 color="black")+
    geom_segment(data=df2,aes(x=x,xend=xend,y=y,yend=yend),
                 color="black")+
    scale_colour_gradient(low = "#F08080" ,high = "#FFE4B5") +
    labs(col="Pvalue",size="Count",x="",y="")+
    scale_size(range = c(4,8),
               limits = c(0,160),
               labels = c(40,80,120,160),
               breaks = c(40,80,120,160))+
    theme(legend.position="top",
          axis.text.x = element_text(color="#000000",size=12),
          axis.text.y = element_text(color="#000000",size=12))
}
## Supplementary Figure 3.C
if(T){
  ggplot(data=df,aes(x=cancer,y=value,fill=Type)) + 
    geom_bar(stat="identity",position="fill",width = 0.8) +
    scale_fill_manual(values=c("Mylnc"=adjustcolor("#F7ED11",alpha.f = 0.7),
                               "Otherlnc"=adjustcolor("#5B90B8",alpha.f = 0.7)))+
    scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                       labels = scales::percent_format())+
    labs(x="",y="Proportion",fill=" ",title="") +
    coord_flip()+
    theme_bw()+
    theme(legend.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}
## Supplementary Figure 3.D
if(T){
  #####  文本信息
  cancer=data.frame(cancer=unique(da$cancer),
                    cancer_y=seq(2,(2*length(unique(da$cancer))),by=2),check.names = F)
  type=data.frame(type=unique(da$cell),
                  type_x=seq(2,(2*length(unique(da$cell))),by=2),check.names = F)
  da=da %>%
    inner_join(cancer,by=c("cancer"="cancer")) %>%
    inner_join(type,by=c("cell"="type")) %>%
    select(cancer,cancer_y,cell,type_x,count,everything())
  
  
  ######  画图信息
  ggplot() +
    geom_scatterpie(data = da,aes(x=type_x, y=cancer_y,r=count/4),  
                    cols=colnames(da)[-(1:5)],color=NA) +
    scale_fill_manual(values = setNames(c("#FFFFFF","#5D9D52","#FFFFFF","#77468E",
                                          "#FFFFFF","#FF7F00","#FFFFFF","#46ACC8"),
                                        colnames(da)[-c(1:5)])) +
    scale_y_continuous(breaks = as.numeric(unique(da$cancer_y)),
                       labels = unique(da$cancer)) +
    scale_x_continuous(breaks = as.numeric(unique(da$type_x)),
                       labels = unique(da$cell)) +
    theme_bw()+
    theme(axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(color = "#000000"),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5,colour = "#000000"),
          legend.title = element_text(),
          panel.background = element_blank(),
          panel.grid=element_blank()) 
}
## Supplementary Figure 3.F
if(T){
  
  ggplot()+
    geom_bar(data = df3,aes(x=`_file`,y=count,fill=type),
             stat = 'identity',
             position="dodge") +
    scale_y_continuous(expand = c(0,0),limits = c(0,700))+
    scale_fill_manual(values = rev(c("#db7b8a","#7bbadb","#dbcc7b")))+
    theme_bw() +
    labs(x="",y="The number of IRaes LncRNAs",fill="Group")+
    theme(axis.text.x = element_text(angle = 45, hjust =0.5 ,vjust = 0.5,colour = "#000000"),
          axis.text.y = element_text(color = "#000000"))
}
## Supplementary Figure 3.G
if(T){
  ggplot() +
    geom_point(data=df,
               aes(x=cancer,y=name,size = count,color=value),
               shape = 15,alpha=1)+
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color="black"),
          axis.ticks = element_blank(),
          axis.text.y=element_text(colour = "#000000"),
          axis.text.x=element_text(angle = 90,vjust = 0.5,hjust = 1,color="#000000"))+
    geom_segment(data=df1,aes(x=x,xend=xend,y=y,yend=yend),
                 color="black")+
    geom_segment(data=df2,aes(x=x,xend=xend,y=y,yend=yend),
                 color="black")+
    scale_colour_gradient(low = "#FFE4B5" ,high = "#F08080") +
    labs(col="-Log10(Pvalue)",size="Count")+
    theme(legend.position="top")
}