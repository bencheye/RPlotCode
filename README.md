# RPlotCode Note

> 关于R语言画图的一些脚本

---

## enrichmentAnalysisCode

> 使用KEGG、GO进行富集分析，并重建自定义画图函数，组合4个癌症的结果。

- **重点**

  1. 主题中设置了保留**绘图边框**

  2. 获取对象的值 通过@

  3. ggplot2 自定义画气泡图

     ```R
     p <- ggplot(data = res[1:10, ], 
              mapping = aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) + 
         geom_point() + mytheme + labs(y='') + 
         scale_color_gradient(low="#1f77b4", high="#cde64c")
     ```

  4. 自定义映射颜色

  5. 取消Y轴标签

  6. ggpubr的ggarrange函数将多个图组合到一个绘图区

---

## networkFeature

> 硕士论文网络拓扑中心性等频分箱表征癌症驱动基因更区域分布特征
>
> 最后使用**分面**的方式画柱状图来可视化，并添加显著富集参考线

- 重点

  1. 画不同分组条件下，**各模块分布条形图**，数据整理方法

     `geom_bar(stat = 'identity', position = 'dodge') `

  2. 如何在画图中改变字符串的默认排序（字母顺序）

     `factor(tmpDat$Level, levels = colnames(dat)[2:dim(dat)[2]])`  levels中是指定的字符串顺序

  3. 手动修改各个模块映射的颜色

     `scale_fill_manual(breaks = factor(dat$Portion), values=rev(mycol)) `

  4. 添加辅助线

     `geom_hline(aes(yintercept = numPvalue), colour = 'red', linetype='dashed')`

---

## plotnetCode

> 硕士论文原始PPI网络中比较**Driver**与**TF**、**DEG**、**Essential**五种中心性的差异
>
> 加入指定不同组别的显著性检验*P*值

- 重点

  1. 画不同分组中各模块**差异比较箱线图**；

     ```R
     ggboxplot(cyto_topo, x = "Group", y = centralityName,
                 color="Group", palette = "jco") 
     # palette = "jco" 参数不用修改
     ```

  2. 增加显著性检验*P*值，可以用*表示，可以更改检验方法，可以加入**多组间比较**

     ```R
     my_comparisons <- list(c("DriverGene", "Essential"), c("DriverGene","TF"), 
                              c("DriverGene", "DEG"))
     # my_comparisons 指定需要检验的分组
     stat_compare_means(comparisons = my_comparisons)
     ```
  3. 给组合图增加标题等注释
    ```R
    annotate_figure(pp1, bottom=text_grob("Degree",color="black", size = 15))
    ```
---

## topology_distribution

>   硕士论文PPI网络拓扑分布点图
>
>   ggplot2如何自定义使用两组数据不同节点画图

-   重点

    1.  自定义标注特殊节点方法

        **使用另外数据集再画一个图层**

        ```R
        pp <- ggplot(plotDat, aes_string('Order', centralityType))+labs(x='') +
              geom_point(size = 2.0, shape = 16, color = 'grey') + mytheme +
              geom_point(data = xTmp, aes_string(x= 'Order', y = centralityType), color = "#FC4E07") 
        ```

---

