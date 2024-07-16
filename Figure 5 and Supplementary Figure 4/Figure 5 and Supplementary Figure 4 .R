## Figure 5.A
if(T){
  ggpubr::ggboxplot(plot.data,
                    x="sample", y="IRaeslnc", width = 0.6, 
                    color = "black",
                    fill="celltype",
                    xlab = "",
                    bxp.errorbar=T,
                    bxp.errorbar.width=0.5, 
                    size=0.5, 
                    outlier.shape=NA,
                    legend = "right")  -> p1
  p1+
    scale_fill_manual(values = color_plot) +
    labs(y="",fill="") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top")
}
## Figure 5.B
if(T){
  ggboxplot(plot.da, x="celltype", y="INFLAMMATORY_RESPONSE", width = 0.6, 
            color = "black",
            fill="celltype",
            palette = col,
            xlab = "", 
            bxp.errorbar=T,
            bxp.errorbar.width=0.5,
            size=0.5,
            outlier.shape=NA,
            legend = "right")  +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top")
}
## Figure 5.C 
if(T){
  ggplot(result1, aes(irAEsLnc, Mac, fill = R)) +
    geom_tile(color = 'white', size = 2) +
    geom_text(aes(label = p.sig), color="black",size = 5) +
    scale_fill_gradient2(low ="blue" , high = "red", mid = "grey90", midpoint = 0) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "#000000",size = 12),
          axis.text.y = element_text(colour = "#000000",size = 12)) +
    labs(title = "",x = "",y = "")
}
## Figure 5.D
if(T){
  ggplot(plot.da,aes(x=subcelltype,y=IRaeslnc_IRgene,fill=sample))+
    geom_violin(position = "dodge",
                width=0.7,trim = F,color=NA)+
    geom_boxplot(aes(x=subcelltype,y=IRaeslnc_IRgene,fill=sample),
                 width=0.2,alpha=1,outlier.shape=NA,color="#FFFFFF",
                 position = position_dodge(width = 0.7),
                 na.rm = T)+
    geom_signif(annotations = wilcox.sig$p.sig[wilcox.sig$type=="w1"],
                y_position = c(rep(0.4,11)),
                xmin = seq(0.785,10.785,1),
                xmax = seq(1.225,11.225,1),
                tip_length = c(rep(c(0.01, 0.01),11)),
                color="black",
                vjust = 0.3)+
    geom_signif(annotations = wilcox.sig$p.sig[wilcox.sig$type=="w2"],
                y_position = c(rep(0.38,11)),
                xmin = seq(1,11,1),
                xmax = seq(1.225,11.225,1),
                tip_length = c(rep(c(0.01, 0.01),11)),
                color="black",
                vjust = 0.3)+
    geom_signif(annotations = wilcox.sig$p.sig[wilcox.sig$type=="w3"],
                y_position = c(rep(0.36,11)),
                xmin = seq(0.785,10.785,1),
                xmax = seq(1,11,1),
                tip_length = c(rep(c(0.01, 0.01),11)),
                color="black",
                vjust = 0.3)+
    scale_fill_manual(values =rev(c("#34ADB4","#46A265","#F27073"))) +
    labs(x="",y="IRaesLnc-IRgene AUCell",fill="Sample")+
    theme_bw()+
    theme(axis.text.y = element_text(size = 12,colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size = 12,colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}
## Figure 5.E
if(T){
  top_annotation <- HeatmapAnnotation(
    Type = anno_block(
      gp = gpar(
        fill = c("#0f444c","#334435"),
        fontsize = 9), 
      labels = c("Colitis","IRaes"),
    ),
    Pathway=c(rep("INF",3),rep("INTa",3),rep("INTg",3),
              rep("INF",3),rep("INTa",3),rep("INTg",3)),
    Group=rep(c("+CPI Colitis","+CPI no Colitis","Control"),6),
    col=list(
      Group=c("Control"=adjustcolor("#34ADB4",alpha.f = 0.8),
              "+CPI no Colitis"=adjustcolor("#46A265",alpha.f = 0.8),
              "+CPI Colitis"=adjustcolor("#F27073",alpha.f = 0.8)),
      Pathway=c("INF"=adjustcolor("#40779d",alpha.f = 0.8),
                "INTa"=adjustcolor("#32666c",alpha.f = 0.6),
                "INTg"=adjustcolor("#87bb99",alpha.f = 0.8),
                "INF"=adjustcolor("#40779d",alpha.f = 0.8),
                "INTa"=adjustcolor("#32666c",alpha.f = 0.6),
                "INTg"=adjustcolor("#87bb99",alpha.f = 0.8)),
      Type=c("Colitis"="#0f444c","IRaes"="#334435")
    )
  )
  col_fun = colorRamp2(c(0, 0.3,0.5,0.8), c("#FFFFFF","#e0faff","#7dcaf1", "#0088cc"))
  pdf("IRaes.Colitis_ssGSEA_VS_AUcell.Cor.pdf",width = 12,height = 6)
  Heatmap(
    as.matrix(plot.da), 
    name = "Cor",
    col = col_fun,
    na_col = "grey85",
    cluster_rows = F,
    cluster_columns = F,
    width = ncol(plot.da)*unit(10, "mm"), 
    height = nrow(plot.da)*unit(10, "mm"),
    show_column_names = F,
    row_names_side = "left",
    column_split = rep(c("IRaes","Colitis"),each=9),
    rect_gp = gpar(col = "white", lwd = 2),
    top_annotation = top_annotation,
    cell_fun = function(j, i, x, y, width, height, fill){
      if(!is.na(plot.da[i, j])){
        grid.text(sprintf("%.2f", plot.da[i, j]), 
                  x,y,gp = gpar(fontsize = 10))
      }
    })
  dev.off()
}
## Figure 5.F
if(T){
  plot %>%
    ggplot(aes(sample, mean, linetype = irae, color = inf.patheay, group = pathway))+
    geom_line(linewidth =1)+ ## 线条粗细
    geom_point(aes(sample, mean,color=sample),size=3)+
    labs(title = "",x = "",y = "")+
    scale_color_manual(values = setNames(c(adjustcolor("#40779d",alpha.f = 1),
                                           adjustcolor("#32666c",alpha.f = 1),
                                           adjustcolor("#87bb99",alpha.f = 1),
                                           "#34ADB4","#46A265","#F27073"),
                                         c("inf","inta","intg","Control",
                                           "+CPI no Colitis","+CPI Colitis")))+
    theme_bw()+
    guides(color = guide_legend(reverse = TRUE,title = "ccc"),
           linetype = guide_legend(reverse = TRUE,title = "aaa")) +
    theme(axis.text.x = element_text(color = "#000000"),
          axis.text.y = element_text(color = "#000000"))
}
## Figure 5.G
if(T){
  col_fun = colorRamp2(c(-0.2, 0,0.2,0.4), c("#b8e6fa","#FFFFFF","#fab8b8","#ed1a1a"))
  pdf("D:\\DLL\\TCGA\\1.单细胞代谢\\代谢通路相关性.pdf",width = 12,height = 5)
  Heatmap(
    as.matrix(plot.da), 
    name = "Cor",
    col = col_fun,
    na_col = "grey85",
    cluster_rows = F,
    cluster_columns = F,
    width = ncol(plot.da)*unit(10, "mm"), 
    height = nrow(plot.da)*unit(10, "mm"),
    show_column_names = T,
    row_names_side = "left",
    rect_gp = gpar(col = "white", lwd = 2),
    cell_fun = function(j, i, x, y, width, height, fill){
      if(!is.na(plot.da[i, j])){
        grid.text(sprintf("%.2f", plot.da[i, j]), 
                  x,y,gp = gpar(fontsize = 10))
      }
    })
  dev.off()
}
## Supplementary Figure 4.E
if(T){
  p=ggboxplot(da@meta.data, x="celltype", y=i, width = 0.6, 
              color = "black",
              fill="sample",
              palette = rev(c("#34ADB4","#46A265","#F27073")),
              xlab = "", 
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=0.5, 
              outlier.shape=NA, 
              legend = "right")  + 
    stat_compare_means(aes(group = sample), label = "p.signif",method = "kruskal.test") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "top")
}
## Supplementary Figure 4.F
if(T){
  p=ggboxplot(da@meta.data, x="celltype", y=i, width = 0.6, 
              color = "black",
              fill="celltype",
              palette = cols,
              xlab = "", 
              bxp.errorbar=T,
              bxp.errorbar.width=0.5, 
              size=1, 
              outlier.shape=NA, 
              legend = "right")  + 
    theme_bw()
}
## Supplementary Figure 4.G
if(T){
  result=do.tissueDist(cellInfo.tb = cellinfo,
                       meta.cluster = cellinfo$celltype,
                       colname.patient = "sample",
                       loc = cellInfo.tb$sample,
                       out.prefix=sprintf("%sColitis.OR.dist",out.prefix),
                       pdf.width=3,
                       pdf.height=5,
                       verbose=0)
}
## Supplementary Figure 4.H
if(T){
  result=do.tissueDist(cellInfo.tb = cellinfo,
                       meta.cluster = cellinfo$subcelltype,
                       colname.patient = "sample",
                       loc = cellInfo.tb$sample,
                       out.prefix=sprintf("%sColitis_Tcell.OR.dist",out.prefix),
                       pdf.width=3,
                       pdf.height=5,
                       verbose=0,z.hi = 8)
}
## Supplementary Figure 4.I
if(T){
  p=DotPlot.metabolism(obj = countexp.Seurat, 
                       pathway = input.pathway, 
                       phenotype = "subcelltype", # 修拉对象元数据中包含的特征之一。
                       norm = "y")  + # 按行或列缩放数值。"x" "y" "na"
    labs(x="") +
    theme(axis.text.x = element_text(color = "#000000",size = 10,angle = 90),
          axis.text.y = element_text(color = "#000000",size = 10,angle = 0))

  p=DotPlot.metabolism(obj = test, 
                       pathway = input.pathway, 
                       phenotype = "test", # 修拉对象元数据中包含的特征之一。
                       norm = "y")  + # 按行或列缩放数值。"x" "y" "na"
    labs(x="") +
    theme(axis.text.x = element_text(color = "#000000",size = 10,angle = 90),
          axis.text.y = element_text(color = "#000000",size = 10,angle = 0))
}