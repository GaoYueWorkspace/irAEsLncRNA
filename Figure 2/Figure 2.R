## Figure 2.A
if(T){
  
  ######  画图信息
  ggplot() +
    geom_scatterpie(data = data,aes(x=type_x, y=cancer_y,r=count),  
                    cols=colnames(data)[-(1:5)],color=NA) +
    scale_fill_manual(values = setNames(distinct(xxx[,c(2,4)])$color,
                                        distinct(xxx[,c(2,4)])$iraes)) +
    scale_y_continuous(breaks = as.numeric(unique(data$cancer_y)),
                       labels = unique(data$cancer))+
    scale_x_continuous(breaks = as.numeric(unique(data$type_x)),
                       labels = unique(data$Type)) +
    geom_text(data = text,aes(x=type_x, y=cancer_y,label=count),size=3)+
    theme_bw()+
    theme(axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          legend.title = element_text(),
          legend.position = "none") +
    ((legends[[1]]+legends[[3]])/(legends[[2]]+legends[[4]])/
       (legends[[5]]+(legends[[6]]/legends[[7]])))
  ggsave("D:\\DLL\\TCGA\\Part2\\0.lncrna_iraes\\1.每种癌症不同器官不同iraes对应的lncrna个数展示.pdf",
         width = 7,height = 9)
}
## Figure 2.B
if(T){
  ############ 和多个iraes相关的lncrna   #############
  name="LUAD"
  ggplot(da,
         aes(x=factor(x,level=c("ORGAN","IRAES","LNCRNA")),
             y=y,
             stratum=factor(stratum,levels = unique(stratum)),
             alluvium=Cohort,
             fill=stratum,
             label=stratum)) +
    #geom_flow()+ ##  连线的宽度
    geom_flow(width = 0.5) +
    geom_stratum(width=0.5,
                 linetype=1,
                 size=0.5,
                 alpha=1,
                 color="white") +
    geom_text(stat="stratum",
              size=2) +
    scale_x_discrete(limits=c()) +
    theme_bw() +
    theme(legend.position="none",
          axis.title=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())  +
    scale_fill_manual(values = setNames(da$color.y,da$stratum))
  
  ggsave(file.path("D:\\DLL\\TCGA\\Part2\\0.lncrna_iraes\\output",paste0("2.",name,".pdf")),
         width = 5,height = (length(unique(data$gene_name)))*7/33,limitsize = FALSE)
  
  #########  只和一个iraes相关的lncrna   #############
  
  vertices_name<-unique(c(as.character(edges$from), as.character(edges$to)))
  vertices<-data.frame(name = vertices_name, value =0.5)
  rownames(vertices)<-vertices_name
  
  d2<-d2%>% 
    left_join(vertices ,by = c("to" = "name")) %>% 
    arrange(order1, desc(group_sum),order2,desc(value))
  edges=rbind(d0,d1[,1:2], d2[,1:2])
  list_unique<-unique(c(as.character(edges$from), as.character(edges$to)))
  vertices = data.frame(
    name = list_unique, 
    value = vertices[list_unique,'value']) 
  vertices$group<-edges$from[ match( vertices$name, edges$to ) ]
  vertices$id<-NA
  myleaves<-which(is.na( match(vertices$name, edges$from) ))
  nleaves<-length(myleaves)
  vertices$id<-90
  vertices$id[myleaves]<-0
  vertices$angle<-90 - vertices$id 
  mygraph <- graph_from_data_frame( edges, vertices = vertices, directed = TRUE )
  
  ####  颜色数据
  col2=data.frame(name=c("NEURO","MUSCULO","DERM","ENDO","GI","RENAL","CARDIO","PULM"),
                  color=c("#3A9EC9","#CCBFD4","#CC789F","#A4C7DF","#E0A8A4","#80699F","#CFE1BA","#AA709A"))
  col=left_join(vertices,color,by=c("name"="iraes")) %>%
    left_join(color,by=c("group"="iraes"))
  col=col[-1,]
  col$color.x[1:length(unique(to2$Type))]=col2$color[which(col2$name %in% col$name)]
  col$color.x[is.na(col$color.x)]=col$color.y[!is.na(col$color.y)]
  col=col %>% 
    mutate(group0=seq(1.001,1+(0.001*(dim(vertices)[1]-1)),0.001)) %>%
    select(name,group,group0,color.x)
  
  #####  画图
  p=ggraph(mygraph, layout = 'dendrogram', circular = FALSE) +
    geom_edge_diagonal(aes(colour=..index..)) +
    scale_edge_colour_gradient2(low = "#9DD1C8",
                                mid = "#9DD1C8",
                                high = "#9DD1C8") +
    geom_node_text(aes(x = x, y=y,  angle = 90,  label=name,color=group),
                   hjust=1.25,size=2.3, alpha=1) +
    geom_node_point(aes( x = x, y=y, fill=group, size=value),
                    shape=21,stroke=0.00,color='black',alpha=1) +
    scale_colour_manual(values= setNames(rep("black",length(c(name,col$name))),c(name,col$name))) +
    scale_fill_manual(values= setNames(c(col$color.x),c(col$name))) +
    scale_size_continuous( range = c(0.1,7) ) +
    expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))+
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0),"cm"))
  
  ggsave(file.path("D:\\DLL\\TCGA\\Part2\\0.lncrna_iraes\\output",paste0("1.",name,".pdf")),plot = p,
         width = dim(d2)[1]/4,height = 5)
}
## Figure 2.C
if(T){
  
  cols=c("#F4B4B4","#C67DA3","#715D70","#46A265","#6F72B6","#FAEACA")
  
  col2=data.frame(id=c(unique(data$Sub_Class)),
                  col=cols[1:length(unique(data$Sub_Class))])
  inner_join(data,col2,by=c("Sub_Class"="id")) -> data
  
  ## 画图
  data %>% 
    ggplot()+
    geom_segment(aes(x=x,xend=xend,y=y.id,yend=y.id),size=1,color="pink")+
    geom_bar(data=select(data,label,value) %>% distinct(),
             aes(x=-2.5,y=value+0.25,fill=label),
             width = 4,stat = "identity",position = "stack") +
    geom_text(data = data.frame(x=rep(-2.5,3),
                                y=c(2.5,14,26.5),
                                text=c("Specific","Nohubnospecific","Hub")),
              aes(x=x,y=y,label=text),color="#FFFFFF",angle=270,size=5) +
    geom_point(aes(x=xend,y=y.id,size=GeneRatio,fill=Sub_Class),shape=21)+
    scale_size(range = c(3,7),
               labels = c(0.1,0.3,0.5),
               breaks = c(0.1,0.3,0.5)) +
    scale_y_continuous(name = "",
                       breaks = c(1:length(data$id)),
                       labels = data$ID,
                       expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0),
                       limits = c(-5,35)) +
    scale_fill_manual(values = setNames(c(col2$col,"#CC0000","#D4E6AF","#25ACAE"),
                                        c(col2$id,unique(data$label)))) +
    theme_classic() +
    theme(axis.ticks.y  = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.text.y = element_text(colour = data$col,size=10),
          axis.text.x = element_text(colour = "#000000",size=10))+
    geom_vline(xintercept = 0,color="black")
  ggsave("D:\\DLL\\TCGA\\Part3\\0.3.1.hub and specific IRaeslnc\\enricher/可视富集分析结果.pdf",
         width = 8,height = 6)
}
## Figure 2.D
if(T){
  
  p=ggplot(data=df,aes(x=Type,y=value,fill=ID)) + 
    geom_bar(stat="identity",position="fill",width = 0.8) +
    scale_fill_manual(values=c("CGGA-GBM"="#b2dec4","Intersect"="#f6b3b7","TCGA-GBM"="#62C7CC",
                               "ICGC-LIHC"="#40E0D0","TCGA-LIHC"="#5F9EA0",
                               "ICGC-LUAD"="#87bb99","TCGA-LUAD"="#53aa73"))+
    scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                       labels = scales::percent_format())+
    labs(x="",y="Proportion",fill=" ",title="") +
    annotate("text", x = 1, y = 0.5, label = paste0("P<","0.05"))+
    annotate("text", x = 2, y = 0.5, label = paste0("P<","0.05"))+  
    annotate("text", x = 3, y = 0.5, label = paste0("P<","0.05"))+  
    theme_bw()+
    theme(axis.title.y=element_text(size=14))+
    theme(legend.text=element_text())+
    theme(axis.text.x = element_text(size = 12, color = "black"))+
    theme(axis.text.y = element_text(size = 12, color = "black"))+
    theme(axis.ticks.length=unit(0.3,"cm"))+
    theme(axis.text.x=element_text(angle=0,size = 11))
  ggsave(file.path("D:\\DLL\\TCGA\\Part2",paste0("4.独立数据集柱状图.pdf")),
         plot=p,width = 5,height = 4)
}
## Figure 2.E
if(T){
  ggplot(data=da,aes(x=Type,y=value,fill=ID)) + 
    geom_bar(stat="identity",position="fill",width = 0.8) +
    scale_fill_manual(values=c("#b2dec4","#f6b3b7","#62C7CC"))+
    scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                       labels = scales::percent_format())+
    annotate("text", x = 1, y = 0.5, label = "****")+
    annotate("text", x = 2, y = 0.5, label = "****")+
    annotate("text", x = 3, y = 0.5, label = "****")+
    annotate("text", x = 4, y = 0.5, label = "NS")+
    labs(x="",y="Proportion",
         fill=" ",title="") +
    theme_bw()+
    theme(axis.title.y=element_text(size=14))+
    theme(legend.text=element_text(size=10))+
    theme(axis.text.x = element_text(size = 12, color = "black"))+
    theme(axis.text.y = element_text(size = 12, color = "black"))+
    theme(axis.ticks.length=unit(0.3,"cm"))+
    theme(axis.text.x=element_text(angle=90,size = 11))
  ggsave(file.path("D:\\DLL\\TCGA\\Part2\\0.lncrna_iraes","4.毒性癌症疾病柱状图展示.pdf"),
         width = 6,height = 4)
}
