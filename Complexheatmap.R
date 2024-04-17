library(devtools)
library(circlize)
library(ComplexHeatmap)

mat <- read.table("C:/Users/86183/Desktop/临时文件/ARGs.txt", sep="\t",quote="",
                 header = T, row.names = 1)

print(mat)
mat <- log10(mat)
as.matrix(mat)
mat[is.na(mat)]<-0
print(mat)
mat <- t(mat)
col_fun = colorRamp2(c(0, 50, 100, 150), c("grey","yellow","blue", "red"))
ha = rowAnnotation(boxplot = anno_boxplot(mat, height = unit(0.5, "cm"), 
                                          box_width = 0.5, outline = FALSE,
                                          gp = gpar(fill = 1:2)))
p_heatmap <- Heatmap(mat,
        width = unit(9, "cm"), height = unit(6, "cm"), 
        name="Abundance",
        col=col_fun,
        #column_title = "I am a column title",
        #row_title = "I am a row title",
        cluster_rows = F,
        cluster_columns = T,
        column_names_rot = 45,
        #row_km = 2,
        row_names_side = "left",
        #column_km = 2, column_title_gp = gpar(fill = c("red", "skyblue"), font = 1:4),
        #top_annotation = ha,
        right_annotation = ha,
        show_parent_dend_line = T,
        #row_order = order(as.numeric(gsub("row", "", rownames(mat)))),
        #column_order = order(as.numeric(gsub("column", "", colnames(mat)))),
        rect_gp = gpar(col = "white", lwd = 2)
        )
p_heatmap
ggsave("C:/Users/86183/Desktop/临时文件/p_heatmap.pdf", p_heatmap, width = 10, height = 6)
p2
