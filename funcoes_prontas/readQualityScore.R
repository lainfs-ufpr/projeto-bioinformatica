# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("ggplot2")
bioc_pkgs <- c("ShortRead")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)
library(ggplot2)

load("cores_lainfs.RData")

# -------------------------------------------
# Função para gerar gráfico da qualidade
# Param -  qa_output: output da funcao qa()
# Return - qualidade_p: ggplot
# -------------------------------------------
readQualityScore <- function(qa_output) {
  
  # Pega tabela com scores de qualidade
  df <- qa_output[["readQualityScore"]]
  
  # Filtrar apenas tipo "read"
  df_filtrada <- subset(df, type == "read")
  
  colnames(df_filtrada) <- c("Qualidade", "Densidade", "Arquivo", "Tipo")
  
  # Gerar gráfico com ggplot2
  p_qualidade <- ggplot(df_filtrada,
                        aes(x = Qualidade, y = Densidade, fill = Arquivo, color = Arquivo)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(linewidth = 1) +
    scale_fill_manual(values = rep_len(cores_lainfs,
                                       length.out = length(unique(df_filtrada$Arquivo)))) +
    scale_color_manual(values = rep_len(cores_lainfs,
                                        length.out = length(unique(df_filtrada$Arquivo)))) +
    labs(
      title = "Distribuição da Qualidade Média por Read",
      x = "Score médio (Phred)",
      y = "Densidade relativa",
      fill = "Arquivo",
      color = "Arquivo"
    ) +
    theme(
      legend.position = "right",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey"))
  
  return(p_qualidade)
}