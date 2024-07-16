## Figure 1.E
if(T){
  rm(list = ls())
  library(rio)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggridges)
  
  # iraes.lnc="D:\\DLL\\TCGA\\Part1\\0.shangceng.LNC/"
  # filename=dir("D:\\DLL\\TCGA\\Part1\\DESeq2",full.names = T) %>%
  #   map(.,~unlist(str_split(.,"\\/"))[2] %>%
  #         str_replace(.,".tsv","")) %>% unlist()
  # iraes.organ=import("D:\\DLL\\TCGA\\Part1\\之前错的\\1.irAES.data\\分类\\irAEs_guildeline_2.csv")
  # 
  # output="D:\\DLL\\TCGA\\1.Plot\\Result 1"
  # shaixuan.lnc=import_list(
  #   dir("D:\\DLL\\TCGA\\Part1\\0.result.P0.05",
  #       full.names = T,
  #       pattern = "Result"),rbind = T
  # )
  # 
  # da=import_list(dir(iraes.lnc,full.names = T),rbind = T) %>%
  #   select(LNC,file) %>%
  #   tidyr::separate(col = "file",into = c("type","type2"),sep = "_") %>%
  #   left_join(iraes.organ,by=c("type"="IRaes")) %>%
  #   select(LNC,Type,type2) %>%
  #   setNames(c("LNC","type1","type2")) %>%
  #   filter(LNC %in% shaixuan.lnc$LNC)
  # 
  # df=df1=as.data.frame(table(da$type1,da$type2)) %>%
  #   arrange(Var1) %>%
  #   group_by(Var1) %>%
  #   mutate(percentage2 = Freq/sum(Freq)) %>%
  #   inner_join(data.frame(table(da$type1)) %>%
  #                mutate(percentage1=Freq/sum(Freq))
  #              ,by=c("Var1"="Var1")) %>%
  #   select(-Freq.x,-Freq.y) %>%
  #   setNames(c("group1","group2","percentage2","percentage1")) %>%
  #   mutate(group22=paste0(group1,"_",group2)) %>% ungroup() %>%
  #   group_by(group1) %>%
  #   mutate(percentage22=percentage2*percentage1) %>% ungroup()
  # df2=df1;df2$x="b";df2$y="a";df1=df2
  # df1$label1 = paste0(df1$group1,'\n',round(df1$percentage1*100,2),"%");df2=df1
  # df2$label2 = paste0(df2$group2,'\n',round(df1$percentage2*100,2),"%")
  # 
  # yanse1=c("CARDIO"="#CFE1BA","RENAL"="#80699F","GI"="#E0A8A4","NEURO"="#3A9EC9",
  #          "DERM"="#CC789F","MUSCULO"="#CCBFD4","ENDO"="#A4C7DF","PULM"="#AA709A",
  #          "CTLA4"="#008B8B","PD1"="#BDB76B","PDL1"="#6495ED")
  # yanse2=data.frame(x=df2$group22);yanse2$y=0
  # yanse2$y[grep("CTLA-4",yanse2$x)]="#008B8B"
  # yanse2$y[grep("PD-1",yanse2$x)]="#BDB76B"
  # yanse2$y[grep("PD-L1",yanse2$x)]="#6495ED"
  # yanse2=setNames(yanse2$y,yanse2$x)
  # yanse=c(yanse1,yanse2)
  # text1=distinct(select(df2,y,percentage1,label1))
  # text1$percentage1=rev(text1$percentage1)
  # text1$label1=rev(text1$label1)
  # pp=ggplot(df2) +
  #   geom_bar(aes(x, percentage22,fill=group22),stat = 'identity', width = .8, color = 'white')
  # text2=distinct(select(df2,x,group22,percentage22,label2))
  # text2$percentage22=((ggplot_build(pp)$data[[1]]$ymax-
  #                        ggplot_build(pp)$data[[1]]$ymin)/2)+ggplot_build(pp)$data[[1]]$ymin
  
  export(list(df2=df2,text1=text1,text2=text2),
         "D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Figure 1.E.xlsx")
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
  rm(list = ls())
  # shaixuan.lnc=import_list(
  #   dir("D:\\DLL\\TCGA\\Part1\\0.result.P0.05",
  #       full.names = T,
  #       pattern = "Result"),rbind = T
  # )
  # da=import_list(dir(iraes.lnc,full.names = T),rbind = T) %>%
  #   select(LNC,IRaes_LNC_R,file) %>%
  #   tidyr::separate(col = "file",into = c("type","type2"),sep = "_") %>%
  #   left_join(iraes.organ,by=c("type"="IRaes")) %>%
  #   select(LNC,IRaes_LNC_R,Type,type2) %>%
  #   setNames(c("LNC","IRaes_LNC_R","type1","type2")) %>%
  #   filter(LNC %in% shaixuan.lnc$LNC) %>%
  #   distinct()
  
  export(da,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Figure 1.F.xlsx")
  
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
  library(tidyverse)
  library(ggforce)
  library(nycflights13)
  
  # qujian=import("D:\\DLL\\TCGA\\Part1\\0.0.sigma.plot\\TCGA-RMSE.sigma.tsv")
  # data=import("D:\\DLL\\TCGA\\Part1\\0.0.sigma.plot\\Sigma_RMSE.txt")
  # data=split(data,data$cancer)
  # lapply(names(data),function(x){
  #   data[[x]] %>% filter(sigma<=qujian$end[qujian$cancer==x]) %>% return()
  # }) %>% bind_rows() -> data
  # 
  # cols=c(c("#B6CCD7","#A6CEE3FF","deepskyblue2","#1F78B4FF","#20ACBD",'blue'),
  #        c("#9AB294","#F4A2A3", "#FD6467","#CCBFD4","#D1817E"),
  #        c("#FFFF99FF","#E1BD6D","darkgoldenrod2","#E58601", "#B15928FF"),
  #        c("#FDBF6FFF","#FF7F00FF","#9C964A","#959897" ,"#5B1A18"),
  #        c("#CAB2D6FF","#AF98B5"  ,"#7294D4","#6A3D9AFF", "#35274A"),
  #        c("#B8D9A9FF","#cadbba"  ,"seagreen"),
  #        c("#cdb262","#17692CFF","#0B775E", "#AA709A"))
  
  export(data,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Supplementary Figure 2.A.xlsx")
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
  # cancer=c("ACC","BLCA","BRCA","CESC","CHOL","COAD",
  #          "GBM","LIHC","LUAD","LUSC","OV","PAAD",
  #          "PRAD","READ","SARC","SKCM","STAD","THCA")
  # id=import("D:\\DLL\\data\\TCGA-GENEID.tsv") %>%
  #   slice(-c(1:4)) %>%
  #   select(1:2) %>%
  #   mutate(gene_id=str_replace(gene_id,"\\..*","")) 
  # iraes_lnc=import_list(grep("Result|Deseq",
  #                            dir("D:\\DLL\\TCGA\\Part1\\0.result.P0.05",full.names = T),
  #                            value=T,invert=T),
  #                       rbind = T) %>%
  #   filter(fdr<0.05) %>%
  #   select(LNC) %>%
  #   distinct()
  # hall_lnc=import("D:\\DLL\\TCGA\\Part1\\0.hallmark\\hallmarker.padj0.05.tsv") %>%
  #   select(gene_name) %>%
  #   distinct() %>%
  #   inner_join(id,by=c("gene_name"="gene_name"))
  # diff_lnc=import_list(dir("D:\\DLL\\TCGA\\Part1\\DESeq2",full.names = T),rbind=T) %>%
  #   filter(padj<0.01 & abs(log2FoldChange)>=2)
  # diff_lnc$`_file`=purrr::map(diff_lnc$`_file`,~unlist(str_split(.,"-"))[2] %>% str_replace(.,".tsv","")) %>% unlist()
  # diff_lnc=filter(diff_lnc,`_file` %in% cancer) %>%
  #   select(Genes) %>%
  #   distinct()
  # 
  # result_lnc=import_list(dir("D:\\DLL\\TCGA\\Part1\\0.result.P0.05",pattern = "Result",full.names = T),rbind = T) %>%
  #   tidyr::separate(`_file`,sep = "\\/",into = c("a","b")) %>%
  #   mutate(b=str_replace_all(b,"Result_|.tsv","")) %>%
  #   select(LNC,b) %>%
  #   distinct() %>%
  #   group_by(b) %>% 
  #   mutate(n=n()) %>%
  #   select(b,n) %>%
  #   distinct()
  
  export(result_lnc,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Supplementary Figure 2.B.xlsx")
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
  library(tidyverse)
  # lncid=import("D:\\DLL\\data\\TCGA-GENEID.tsv") %>%
  #   slice(-c(1:4)) %>%
  #   select(gene_id,gene_name) %>%
  #   mutate(gene_id=str_replace(gene_id,"\\..*",""))
  # cancers=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","GBM",
  #           "LIHC","LUAD","LUSC","MESO","OV","PAAD","PRAD",
  #           "READ","SARC","SKCM","STAD","THCA")
  # dat<-import_list(dir("D:\\DLL\\TCGA\\Part1\\DESeq2",full.names = T),rbind = T) %>%
  #   inner_join(lncid,by=c("Genes"="gene_id")) 
  # dat$`_file`=purrr::map(dat$`_file`,~unlist(str_split(.,"\\/|\\.|-"))[3]) %>% unlist()
  # dat=select(dat,log2FoldChange,padj,`_file`,gene_name) %>%
  #   filter(`_file` %in% cancers)
  # colnames(dat)=c("log2FC","fdr","compared.group","proteins")
  # 
  # log2Foldchang=2
  # adjp=0.05
  # top_marker=1
  # max_overlaps=5
  # dat.plot <- dat %>% mutate(
  #   "significance"=case_when(fdr < adjp & log2FC>= log2Foldchang  ~ 'up',
  #                            fdr < adjp &log2FC<= -log2Foldchang  ~ 'down',
  #                            TRUE ~ 'insig'))
  # dat.plot$compared.group <- factor(dat.plot$compared.group,
  #                                   levels = unique(dat.plot$compared.group))
  # background.dat <- data.frame(
  #   dat.plot %>% group_by(compared.group) %>% filter(log2FC>2) %>% 
  #     summarise("y.localup"=max(log2FC)),
  #   dat.plot %>% group_by(compared.group) %>% filter(log2FC<2) %>% 
  #     summarise("y.localdown"=min(log2FC)),
  #   x.local=seq(1:length(unique(dat.plot$compared.group)))
  # ) %>% select(-compared.group.1)
  # names(background.dat)
  # 
  # x.number <- background.dat %>% select(compared.group,x.local) 
  # dat.plot <- dat.plot%>% left_join(x.number,by = "compared.group")
  # names(dat.plot)
  # dat.marked.up <- dat.plot %>% filter(significance=="up") %>% 
  #   group_by(compared.group) %>% arrange(-log2FC) %>% 
  #   top_n(top_marker,abs(log2FC))
  # dat.marked.down <- dat.plot %>% filter(significance=="down") %>% 
  #   group_by(compared.group) %>% arrange(log2FC) %>% 
  #   top_n(top_marker,abs(log2FC))
  # dat.marked <- dat.marked.up %>% bind_rows(dat.marked.down)
  # dat.infor <- background.dat %>% 
  #   mutate("y.infor"=rep(0,length(compared.group)))
  # names(dat.infor)
  # color.pals=c(c("#B6CCD7","#A6CEE3FF","deepskyblue2","#1F78B4FF","#20ACBD",'blue'),
  #              c("#9AB294","#F4A2A3", "#FD6467","#CCBFD4","#D1817E"),
  #              c("#FFFF99FF","#E1BD6D","darkgoldenrod2","#E58601", "#B15928FF"),
  #              c("#FDBF6FFF","#FF7F00FF","#9C964A","#959897" ,"#5B1A18"),
  #              c("#CAB2D6FF","#AF98B5"  ,"#7294D4","#6A3D9AFF", "#35274A"),
  #              c("#B8D9A9FF","#cadbba"  ,"seagreen"),
  #              c("#cdb262","#17692CFF","#0B775E", "#AA709A"))
  
  export(list(background.dat=background.dat,
              dat.plot=dat.plot,
              dat.infor=dat.infor,
              dat.marked.up=dat.marked.up,
              dat.marked.down=dat.marked.down
              ),"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Supplementary Figure 2.C.xlsx")
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
  # da=import("D:\\DLL\\TCGA\\Part1\\0.hallmark\\hallmarker.padj0.05.tsv")
  # table(da$`_file`,da$sum) %>%
  #   as.data.frame()-> plot.da
  
  export(plot.da,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Supplementary Figure 2.D.xlsx")
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
  library(dplyr)
  library(rio)
  library(stringr)
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(aplot)
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
  export(plots,"D:\\DLL\\TCGA\\0.文章\\终终终\\数据和代码\\代码\\图代码\\Supplementary Figure 2.F.xlsx")
  bind_rows(plots) %>% 
    export(.,"D:\\DLL\\TCGA\\1.Plot\\补充图\\每种癌症预测和真实相关性\\pearson.tsv")
}