# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("dplyr", "ggplot2")
bioc_pkgs <- c("ShortRead")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(dplyr)
library(ggplot2)
library(stats)

load("cores_lainfs.RData")

# -------------------------------------------
# Função para verificar a frequência das bases por read
# Param -  caminho_fastq: caminho para o arquivo .fastq
# Return - p_contagens: ggplot com boxplot para contagem de bases
# -------------------------------------------
plotOcurrences <- function(qa_output) {
  df <- qa_output[["sequenceDistribution"]]
  df_ocorrencias <- df %>% filter(type == "read")
  
  # Ordenar por n Occurrences e calcular a frequência cumulativa
  df_ocorrencias <- df_ocorrencias %>%
    arrange(nOccurrences) %>%
    mutate(
      freq = nReads / sum(nReads),
      cumFreq = cumsum(freq))
  
  colnames(df_ocorrencias) <- c("nOccurrences", "nReads", "type",
                                "Arquivo", "freq", "cumFreq")
  
  p_ocorrencias <- ggplot(df_ocorrencias, aes(x = log10(nOccurrences),
                                              y = cumFreq, group=Arquivo)) +
    geom_line(aes(color = Arquivo), linewidth = 1.2, alpha = 0.7) +
    geom_point(aes(color = Arquivo), size = 2, alpha = 0.7) +
    scale_color_manual(values = cores_lainfs) +
    labs(
      title = "Distribuição Cumulativa de Leituras por Frequência de Ocorrência",
      x = "Número de ocorrências por leitura (log10)",
      y = "Proporção cumulativa das leituras") +
    theme(
      legend.position = "right",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey"))
  
  return(p_ocorrencias)
}