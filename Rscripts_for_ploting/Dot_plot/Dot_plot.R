library(ggplot2)
plot_data <- read.table("C:/Users/86183/Desktop/��ʱ�ļ�/plot.txt", sep="\t",quote="",
                        header = T, row.names = 1)
ggplot(data = plot_data,aes(x=GC_content,y=coverage,colour=group,shape=Species))+
  geom_point(size=4) +
  theme(axis.line = element_line(colour = "black"))+
  labs(x= " ", y = "Accurracy (%)")+
  theme(
    text = element_text(size=20),  #���������СΪ20
    axis.text.x = element_text(angle=60, vjust = 1, colour = "black",hjust =1),
    axis.text.y = element_text(colour = "black"),#X������90�Ȱڷ�,ˮƽ΢��0.5
    # axis.title.x=element_blank(),  #ɾ��X�����
    # axis.title.y=element_blank(),  #ɾ��Y�����
    panel.background = element_blank(), #ɾ����ɫ����
    #scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8))
  )+
  scale_shape_manual(values = c(15, 16, 17, 18, 19, 7, 9, 8))+
  scale_color_manual(values = c("red", "green"))