# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("dplyr", "ggplot2", "plotly")
bioc_pkgs <- c("ShortRead")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(dplyr)
library(ggplot2)
library(stats)
library(plotly)

# -------------------------------------------
# Função para verificar a frequência das bases por read
# Param -  caminho_fastq: caminho para o arquivo .fastq
# Param -  paleta_cores: string com nome da paleta escolhida pelo usuario
# Param -  nomes_arquivos: nomes originais dos arquivos de entrada
# Return - lista_plots: lista nomeada com plot estático e plot interativo
# -------------------------------------------
plotOcurrences <- function(qa_output, paleta_cores="viridis", nomes_arquivos) {
  df <- qa_output[["sequenceDistribution"]]
  df_ocorrencias <- df %>% filter(type == "read")
  
  # Ordenar por n Occurrences e calcular a frequência cumulativa
  df_ocorrencias <- df_ocorrencias %>%
    arrange(nOccurrences) %>%
    mutate(
      freq = nReads / sum(nReads),
      cumFreq = cumsum(freq))
  
  colnames(df_ocorrencias) <- c("nOccurrences", "nReads", "Tipo",
                                "Arquivo", "freq", "cumFreq")
  # Ajusta nomes dos arquivos
  df_ocorrencias$Arquivo <- nomes_arquivos[
    match(df_ocorrencias$Arquivo,
          unique(df_ocorrencias$Arquivo))
  ]
  
  # Gerar cores 
  n_cores <- length(unique((df_ocorrencias$Arquivo)))
  cores <- viridis::viridis(n_cores, option = paleta_cores)
  
  # Montar plot
  p_estatico <- ggplot(df_ocorrencias, aes(x = log10(nOccurrences),
                                              y = cumFreq, group=Arquivo)) +
    geom_line(aes(color = Arquivo), linewidth = 2, alpha = 0.5) +
    geom_point(aes(color = Arquivo), size = 3, alpha = 0.5) +
    scale_color_manual(values = cores) +
    labs(
      title = "Distribuição Cumulativa por Número de Ocorrências",
      x = "Número de ocorrências por leitura (log10)",
      y = "Proporção cumulativa das leituras") +
    theme(
      legend.position = "right",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey"))
  
  # Cria gráfico interativo
  p_interativo <- ggplotly(p_estatico)
  p_interativo <- p_interativo %>%
    layout(
      title = list(
        font = list(
          size = 16,
          color = "black"
        )
      ),
      xaxis = list(
        titlefont = list(
          size = 14,
          color = "black"
        ),
        showgrid = FALSE,      
        zeroline = FALSE,   
        showline = TRUE,
        linecolor = "gray",      
        linewidth = 1         
      ),
      yaxis = list(
        titlefont = list(
          size = 14,
          color = "black"
        ),
        showgrid = FALSE,      
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "gray",    
        linewidth = 1
      )
    )
  
  # Retorna lista com os dois gráficos (estático e interativo)
  lista_plots <- list(p_interativo = p_interativo, p_estatico = p_estatico)
  
  return(lista_plots)
}