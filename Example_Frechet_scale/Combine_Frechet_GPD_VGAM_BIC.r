rm(list = ls())

library(ggplot2)
library(reshape2)
library(patchwork)

#Combinacion Frechet GPD Pixel Plot -----

load("../Example_Frechet_scale//MC_study_p5/MC_study_VGAM_BIC_n200/Matrix_n200.Rdata")
matrix_result_n200 = matrix_result
load("../Example_Frechet_scale/MC_study_p5/MC_study_VGAM_BIC_n400/Matrix_n400.Rdata")
matrix_result_n400 = matrix_result
load("../Example_Frechet_scale/MC_study_p5/MC_study_VGAM_BIC_n800/Matrix_n800.Rdata")
matrix_result_n800 = matrix_result


# Suponiendo que matrix_result_n200 y matrix_result_n400 son matrices ya cargadas
# Convertir ambas matrices a formato largo
df1 <- melt(matrix_result_n200)
df2 <- melt(matrix_result_n400)
df3 <- melt(matrix_result_n800)

# Nombres para los ejes, como factores para mantener el orden
df1$Var1 <- factor(df1$Var1, levels = rev(rownames(matrix_result_n200)))
df1$Var2 <- factor(df1$Var2, levels = colnames(matrix_result_n200))

df2$Var1 <- factor(df2$Var1, levels = rev(rownames(matrix_result_n400)))
df2$Var2 <- factor(df2$Var2, levels = colnames(matrix_result_n400))

df3$Var1 <- factor(df3$Var1, levels = rev(rownames(matrix_result_n800)))
df3$Var2 <- factor(df3$Var2, levels = colnames(matrix_result_n800))


# Crear ambos plots
p1 <- ggplot(df1, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.01)) +
  labs(title = "n = 200", x = NULL, y = "Frechet \n True Model") +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 15)
  )


p2 <- ggplot(df2, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.01)) +
  labs(title = "n = 400", x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 15)
  )

p3 <- ggplot(df3, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.01)) +
  labs(title = "n = 800", x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 15)
  )


load("../Example_Pareto_Generalized_scale/MC_study_p5/MC_study_VGAM_BIC_n200/Matrix_n200.Rdata")
matrix_result_n200 = matrix_result
load("../Example_Pareto_Generalized_scale/MC_study_p5/MC_study_VGAM_BIC_n400/Matrix_n400.Rdata")
matrix_result_n400 = matrix_result
load("../Example_Pareto_Generalized_scale/MC_study_p5/MC_study_VGAM_BIC_n800/Matrix_n800.Rdata")
matrix_result_n800 = matrix_result



# Suponiendo que matrix_result_n200 y matrix_result_n400 son matrices ya cargadas
# Convertir ambas matrices a formato largo
df4 <- melt(matrix_result_n200)
df5 <- melt(matrix_result_n400)
df6 <- melt(matrix_result_n800)

# Nombres para los ejes, como factores para mantener el orden
df4$Var1 <- factor(df1$Var1, levels = rev(rownames(matrix_result_n200)))
df4$Var2 <- factor(df1$Var2, levels = colnames(matrix_result_n200))

df5$Var1 <- factor(df2$Var1, levels = rev(rownames(matrix_result_n400)))
df5$Var2 <- factor(df2$Var2, levels = colnames(matrix_result_n400))

df6$Var1 <- factor(df3$Var1, levels = rev(rownames(matrix_result_n800)))
df6$Var2 <- factor(df3$Var2, levels = colnames(matrix_result_n800))


# Crear ambos plots
p4 <- ggplot(df4, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.01)) +
  labs(title = NULL, x = "Model Estimate", y = "GP \n True Model") +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 15)
  )


p5 <- ggplot(df5, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.01)) +
  labs(title = NULL, x = "Model Estimate", y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 15)
  )

p6 <- ggplot(df6, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "gray90") +
  scale_fill_gradient(name = "Prob", low = "white", high = "deepskyblue", limits = c(0, 1.01)) +
  labs(title = NULL, x = "Model Estimate", y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 15)
  )


# Combinar horizontalmente
final_plot = p1 + p2 + p3 +  p4 + p5 + p6 + plot_layout(guides = "collect") & theme(legend.position = "right")
final_plot

path_save = "../Example_Frechet_scale/MC_study_p5/Img_Study/"
ggsave(filename = paste(path_save,"Combine_Pixel_plot_VGAM_BIC", ".png", sep= ""), plot = final_plot,
       width = 15, height = 10, dpi = 500)