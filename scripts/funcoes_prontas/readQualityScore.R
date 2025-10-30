# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("scripts/instalar_pacotes.R")
cran_pkgs <- c("ggplot2", "plotly")
bioc_pkgs <- c("ShortRead")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)
library(ggplot2)
library(plotly)

# -------------------------------------------
# Função para gerar gráfico da qualidade
# Param -  qa_output: output da funcao qa()
# Param -  paleta_cores: string com nome da paleta escolhida pelo usuario
# Param -  nomes_arquivos: nomes originais dos arquivos de entrada
# Return - lista_plots: lista nomeada com plot estático e plot interativo
# -------------------------------------------
readQualityScore <- function(qa_output, paleta_cores="viridis", nomes_arquivos) {
  
  # Pega tabela com scores de qualidade
  df <- qa_output[["readQualityScore"]]
  
  # Filtrar apenas tipo "read"
  df_filtrada <- subset(df, type == "read")
  colnames(df_filtrada) <- c("Qualidade", "Densidade", "Arquivo", "Tipo")
  
  # Ajusta nomes dos arquivos
  df_filtrada$Arquivo <- nomes_arquivos[
    match(df_filtrada$Arquivo,
          unique(df_filtrada$Arquivo))
  ]
  
  # Gerar cores 
  n_cores <- length(unique(df_filtrada$Arquivo))
  cores <- viridis::viridis(n_cores, option = paleta_cores)
  
  # Expandir se a paleta for menor que n_cores
  if (length(cores) < n_cores) {
    cores <- rep(cores, length.out = n_cores)
  }

  # Gerar gráfico com ggplot2
  p_qualidade <- ggplot(df_filtrada,
                        aes(x = Qualidade, y = Densidade, fill = Arquivo, color = Arquivo)) +
    geom_area(alpha = 0.3, position = "identity") +
    geom_line(linewidth = 2, show.legend = FALSE, alpha = 0.4) +
    scale_fill_manual(values = rep_len(cores,
                                       length.out=length(unique(df_filtrada$Arquivo)))) +
    scale_color_manual(values = rep_len(cores,
                                       length.out=length(unique(df_filtrada$Arquivo)))) +
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
  
  # Cria gráfico interativo
  p_interativo <- ggplotly(p_qualidade)
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
  return(list(
    p_estatico = p_qualidade,
    p_interativo = p_interativo))
}