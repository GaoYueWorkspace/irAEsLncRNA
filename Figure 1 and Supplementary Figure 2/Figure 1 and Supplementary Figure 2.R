## Figure 1.E
if(T){
  
  ###  画图
  p=ggplot(df2) +
    geom_bar(data=distinct(select(df2,y,group1,percentage1)),
             aes(x=y,y=percentage1,fill = group1),
             stat = 'identity', 
             width = 1.3,color = NA,alpha = 0.5) +
    geom_text(data=text1,
              aes(x=y,y=percentage1,label = label1), 
              position = position_stack(vjust = 0.5),size=3)+
    geom_bar(aes(x, percentage22,fill=group22),
             stat = 'identity', width = .8, color = 'white',alpha = 0.3) +
    geom_text(data=text2,
              aes(x=x, y=percentage22,label = label2),
              size = 2.5, color = 'black') +
    scale_fill_manual(values = yanse) +
    scale_y_continuous(labels = scales::percent) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = 'none') 
  ggsave(file.path(output,"1.B.筛选后的LNCrna个数展示.pdf"),
         plot = p,width = 8,height = 8)
  num=import("D:\\DLL\\TCGA\\Part1\\之前错的\\2.Lncrna-ICP\\Gene_id.lncRNA.csv") %>%
    mutate(LNC=str_replace(gene_id,"\\.\\d.*",""))
  da=da %>% mutate(LNC=str_replace(LNC,"\\.\\d.*",""))
  su=length(unique(da$LNC))/length(unique(num$gene_id))
  df=data.frame(group=c("P≥0.05","P<0.05"),
                percentage=c(1-su,su),check.names = F) %>%
    mutate(label=paste0(group,"\n",round(percentage*100,2),"%"))
  percent=rev(df$percentage)
  
  ### 画图
  p=ggplot(df,aes(x="",y=percentage, fill=group)) +
    geom_bar(stat = "identity",color="white",width = 0.2) +
    geom_text(aes(y= cumsum(percent)-percent/2, x= 1),label=df$label,size=4)+
    theme_classic() +
    xlab("")+
    theme(legend.position = 'none') 
  ggsave(file.path(output,"1.B.筛选前后的LNCrna个数展示.pdf"),
         plot = p,width = 8,height = 8)
}
## Figure 1.F
if(T){
  
  ggplot(da,aes(x=IRaes_LNC_R,y=type1,fill=type1,height = after_stat(density)))+
    geom_density_ridges(alpha=0.7,stat = "density")+
    theme(legend.position = "none")+
    scale_fill_manual(values = yanse) +
    theme_classic() +
    labs(x="",y="",fill="Organ Group") +
    theme(axis.text.x = element_text(colour = "#000000",size=12),
          axis.text.y = element_text(colour = "#000000",size=12)) -> p2
  ggsave(file.path(output,"1.D.山脊图.pdf"),
         plot = p2,width = 10,height = 7)
}
## Supplementary Figure 2.A
if(T){
  
  p=ggplot(data = data,aes(x=sigma,y=value,color=cancer,label = cancer))+
    geom_line(linewidth=1.5,alpha=0.5) +
    facet_zoom(xlim = c(9.1, 10), ylim = c(0.000355, 0.00048), split = T,zoom.size = 2)+
    labs(x="Sigma Value",y="Root Mean Squared Error")+theme_bw() +
    scale_color_manual(name = "Cancers",values=cols)+
    theme(axis.text.y =  element_text(color="black"),
          axis.title.y=element_text(size=14),
          axis.text.x =  element_text(color="black"),
          axis.title.x=element_text(size=14),
          legend.title = element_text(size=14),
          legend.background = element_blank(),
          legend.key = element_blank(),
          # panel.grid=element_blank(),
          legend.box.background = element_rect(fill=NA,color = "black",linetype = 1))
  p
  ggsave("D:\\DLL\\TCGA\\Part1\\0.0.sigma.plot\\sigma.vis.pdf",plot = p)
}
## Supplementary Figure 2.B
if(T){
  ggplot()+
    geom_bar(data = result_lnc,
             aes(x=b,y=n),
             stat = 'identity',
             fill="#256491",
             width = 0.9) +
    geom_text(data = result_lnc, 
              aes(x=b,y=n,label = n),
              position = position_dodge(0.9),
              vjust = -0.5, hjust = 0.5,color="#000000")+
    theme_bw() +
    labs(x="",y="The number of irAEsLnc")+
    theme(axis.text.x = element_text(angle = 45, hjust =0.5 ,vjust = 0.5,colour = "#000000"),
          axis.text.y = element_text(color = "#000000")) -> p
  ggsave("补充图1.D.pdf",width = 6,height = 4,plot=p)
}
## Supplementary Figure 2.C
if(T){
  set.seed(1)
  vol.plot <- ggplot()+
    geom_col(background.dat,mapping=aes(x.local,y.localup),
             fill="grey50",alpha=0.2,width=0.9,just = 0.5)+
    geom_col(background.dat,mapping=aes(x.local,y.localdown),
             fill="grey50",alpha=0.2,width=0.9,just = 0.5)+
    geom_jitter(dat.plot,mapping=aes(x.local,log2FC,
                                     color=significance,
                                     fill=significance),
                size=1.5,width = 0.4,alpha= 0.4)+
    scale_color_manual(values = c("#5390b5","#eaebea","#d56e5e")) +
    geom_tile(dat.infor,mapping=aes(x.local,y.infor,fill=compared.group,
                                    color = compared.group),
              height=log2Foldchang*1.5,
              color = color.pals[1:length(unique(dat.plot$compared.group))],
              fill = color.pals[1:length(unique(dat.plot$compared.group))],
              alpha = 0.6,
              width=0.9)+guides(size=guide_legend(title="Count"))+ 
    labs(x=NULL,y="log2 Fold change")+
    geom_text(dat.infor,mapping=aes(x.local,y.infor,label=compared.group),size=4)+
    ggrepel::geom_label_repel(dat.marked.up,mapping=aes(x.local,log2FC,label=proteins,color=significance),
                              force = 2,size=3, 
                              max.overlaps = max_overlaps,
                              seed = 233,
                              min.segment.length = 0,
                              force_pull = 2,
                              box.padding = 0.1,
                              segment.linetype = 3, 
                              segment.color = 'black', 
                              segment.alpha = 0.5, 
                              direction = "x", 
                              hjust = 0.5)+
    ggrepel::geom_label_repel(dat.marked.down,mapping=aes(x.local,log2FC,label=proteins,color=significance),
                              force = 2,size=3, 
                              max.overlaps = max_overlaps,
                              seed = 233,
                              min.segment.length = 0,
                              force_pull = 2,
                              box.padding = 0.1,
                              segment.linetype = 3, 
                              segment.color = 'black', 
                              segment.alpha = 0.5, 
                              direction = "x", 
                              hjust = 0.5)+
    annotate("text", x=8, y=max(background.dat$y.localup)+2, 
             label=paste0("|log2FC|>=",log2Foldchang," & FDR<0.05"))+
    theme_classic()+
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.spacing.x=unit(0.1,'cm'),
      legend.key.width=unit(0.5,'cm'),
      legend.key.height=unit(0.5,'cm')
    )
  ggsave('D:\\DLL\\TCGA\\1.Plot\\补充图/各癌症差异基因.pdf', width = 10, height = 5,plot = vol.plot)
}
## Supplementary Figure 2.D
if(T){
  ggplot()+
    geom_bar(data = plot.da,aes(x=Var1,y=Freq,fill=Var2),
             stat = 'identity',
             position="dodge") +
    scale_y_continuous(expand = c(0,0),limits = c(0,700))+
    scale_fill_manual(values = rev(c("#db7b8a","#7bbadb","#dbcc7b")))+
    theme_bw() +
    labs(x="",y="The number of IRaes LncRNAs",fill="Times")+
    theme(axis.text.x = element_text(angle = 45, hjust =0.5 ,vjust = 0.5,colour = "#000000"),
          axis.text.y = element_text(color = "#000000")) -> p
  geom_text(data = plot.da, 
            aes(x=method,y=AUC,fill=gene,label = round(AUC,2)),
            position = position_dodge(0.9),
            vjust = -0.5, hjust = 0.5,color="#000000")ggsave(file.path(output,"补充图1.c.pdf"),width = 6,height = 4,plot = p)
}
## Supplementary Figure 2.F
if(T){
  
  cancers=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","GBM",
            "LIHC","LUAD","LUSC","OV","PAAD","PRAD",
            "READ","SARC","SKCM","STAD","THCA")
  sigma_da=import("D:\\DLL\\TCGA\\Part1\\0.0.sigma.plot\\TCGA-RMSE.sigma.tsv")
  lapply(cancers, function(cancer){
    print(cancer)
    # cancer="READ"
    result.path=grep(cancer,dir("D:\\DLL\\TCGA\\Part1\\0.result.P0.05",full.names = T,pattern = "^Res"),value = T)
    if(length(result.path)!=0){
      result=import(result.path) %>%
        mutate(id=paste0(LNC,"_",file))
      df=import(grep(cancer,dir("D:\\DLL\\TCGA\\Part1\\0.Score",full.names = T),value = T)) %>%
        mutate(id=paste0(LNC,"_",file)) %>%
        filter(id %in% result$id,abs(IRaes_LNC_R)>0.6) %>%
        filter(sigma == as.character(sigma_da$state[grep(cancer,sigma_da$cancer)])) %>%
        select(IRaes_LNC_R,pre) %>%
        setNames(c("x","y"))
      # return(df)
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
        labs(x="Observed value",y="Predictive value")+
        stat_cor(method = "pearson")

      p2<-ggplot(df,aes(x))+
        geom_density(fill="#FED439",alpha=1,lwd=NA)+
        scale_y_continuous(expand = c(0,0))+
        labs(title = cancer)+
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
        insert_top(p2,height = 0.2)%>%
        insert_right(p3, 0.2)
      ggsave(file.path("D:\\DLL\\TCGA\\1.Plot\\补充图\\每种癌症预测和真实相关性",paste0(cancer,".pdf")),
             plot=p,width = 5,height = 5)
      cor=cor.test(df$x,df$y,method = "pearson")
      return(data.frame(cancer=cancer,cor=cor$estimate,pvalue=cor$p.value))
    }
  }) -> plots
  names(plots)=cancers
}