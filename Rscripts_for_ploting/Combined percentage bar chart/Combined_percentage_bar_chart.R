library(ggplot2)
library(cols4all)
library(patchwork)
library(readxl)
library(tidyverse)

# 读取数据
dat01 <- read_excel("C:/Users/86183/Downloads/Figure1b.xlsx", sheet = "long_size", na = "NA")
dat02 <- read_excel("C:/Users/86183/Downloads/Figure1b.xlsx", sheet = "hybrid_size", na = "NA")

# 创建新的 x 轴变量
dat01 <- dat01 %>% 
  mutate(new_x = rep(paste0('X', formatC(1:25, width = 2, flag = 0)), each=5))

dat02 <- dat02 %>% 
  mutate(new_x = rep(paste0('X', formatC(1:25, width = 2, flag = 0)), each=5))

# 计算每个 new_x 的总数并添加回原始数据框
dat01 <- dat01 %>%
  group_by(new_x) %>%
  mutate(total_n = sum(n)) %>%
  ungroup()

dat02 <- dat02 %>%
  group_by(new_x) %>%
  mutate(total_n = sum(n)) %>%
  ungroup()



# 定义 fill 映射的顺序，确保 >50 kbp 在底部
fill_levels <- c("1-5 kbp","5-10 kbp","10-25 kbp", "25-50 kbp", "> 50 kbp")

# 初始图形 p1
p2 <- ggplot(data = dat02, aes(x = n, y = new_x, fill = factor(rlCodes, levels = fill_levels))) +
  geom_bar(stat = "identity", position = "fill", alpha = 1) +
  theme_classic() +
  theme(axis.title.y = element_blank(), # 去掉 y 轴标题
        axis.text.y = element_blank(),  # 去掉 y 轴标签
        axis.title.x = element_blank()) + 
  scale_fill_manual(values = c("1-5 kbp" = "#98d09d", "5-10 kbp" = "#BBCF85",
                               "10-25 kbp" = "#61D8D6", "25-50 kbp" = "#D4BBFC",
                               "> 50 kbp" = "#f7a895"),
                    limits = fill_levels,
                    name = "")

# 调整图形 p2
p1 <- ggplot(data = dat01, aes(x = n, y = new_x, fill = factor(rlCodes, levels = fill_levels))) +
  geom_bar(stat = "identity", position = "fill", alpha = 1) +
  theme_classic() +
  theme(axis.title.y = element_blank(), # 去掉 y 轴标题
        axis.text.y = element_blank(),  # 去掉 y 轴标签
        axis.title.x = element_blank()) + 
  scale_fill_manual(values = c("1-5 kbp" = "#98d09d", "5-10 kbp" = "#BBCF85",
                               "10-25 kbp" = "#61D8D6", "25-50 kbp" = "#D4BBFC",
                               "> 50 kbp" = "#f7a895"),
                    limits = fill_levels,
                    name = "") +
  scale_y_discrete(position = "right") + # y 轴逆转
  scale_x_reverse(expand = c(0,0)) +     # 数值逆转
  theme(legend.position = "top",         # 图例放顶部
        legend.title = element_blank())  # 去掉图例标题

# 打印最终图形
p1 + p2
