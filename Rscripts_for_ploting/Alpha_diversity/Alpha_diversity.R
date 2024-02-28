library(vegan)
library(picante)      
library(ggplot2)
library(ggpubr)

alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #注意，这里是Gini-Simpson 指数
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  #保留四位小数
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  
  
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
 
  result
}
alpha <- read.table("C:/Users/86183/Desktop/临时文件/all_ARGs_abundance.txt", sep="\t",quote="",
                  header = T, row.names = 1)
alpha_g <- read.table("C:/Users/86183/Desktop/临时文件/ARGs_group.txt", sep="\t",quote="",
                    header = T, row.names = 1)
alpha_g
alpha <- t(alpha)
alpha_plasmid <- alpha_diversity (alpha)
alpha_plasmid$group <- alpha_g$group
alpha_plasmid

ggplot(alpha_plasmid, aes(group, x = factor(group,
                                            levels = c("IF",
                                                       "AS")), y = Chao1, fill = group)) +
  #geom_boxplot(outlier.shape = NA)+
  theme_classic() +
  geom_boxplot(width = 0.45, outlier.color = "red", outlier.shape = 5, outlier.size = 5) + 
  geom_jitter(shape = 16, size = 5, position = position_jitter(0.2)) +
  labs(x= " ", y = "Chao1")+
  #stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.")+
  stat_compare_means(comparisons = list(c("IF", "AS")
  ),
  label = "p.signif", method = "wilcox.test") + # 使用wilcox.test方法进行非参数检验
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 20, vjust = 0.5, hjust = 1, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        legend.key.size = unit(c(1, 1), "cm"),
        legend.text = element_text(size = 20),
        panel.grid.major = element_line(colour = NA))
alpha_plasmid