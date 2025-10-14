# -------------------------------------------
# Carregamento das Bibliotecas
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("ggplot2", "plotly", "htmlwidgets", "dplyr")
bioc_pkgs <- c("ShortRead")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)
library(dplyr)
library(plotly)
library(htmlwidgets)
library(ggplot2)

# -------------------------------------------
# Função para gerar gráfico interativo e estático
# da qualidade do sequenciamento por ciclo 
# Param -  qa_output: output da funcao qa()
# Param -  paleta_cores: cor para linhas principais do plot
# Param -  nomes_arquivos: nomes originais dos arquivos de entrada
# Return - lista_plots: lista nomeada com plot estático e plot interativo
# -------------------------------------------
plotCycleQuality <- function(qa_output, paleta_cores="viridis", nomes_arquivos) {
  
  # Organiza dataframe
  df <- qa_output[["perCycle"]][["quality"]]
  
  df_mediana <- df %>%
    group_by(lane, Cycle) %>%
    summarise(Median = quantile(rep(Score, Count), 0.5),
              .groups = "drop")
  
  colnames(df_mediana) <- c("Arquivo", "Ciclo", "Mediana")
  x_max <- max(df_mediana$Ciclo)
  
  # Ajusta nomes dos arquivos
  df_mediana$Arquivo <- nomes_arquivos[
    match(df_mediana$Arquivo,
          unique(df_mediana$Arquivo))
  ]
  
  # Gerar cores 
  n_cores <- length(unique((df_mediana$Arquivo)))
  cores_paleta <- viridis::viridis(n_cores, option = paleta_cores)
  
  # Separar cores em arquivos
  arquivo_niveis <- unique(df_mediana$Arquivo)
  cores <- setNames(cores_paleta, arquivo_niveis)
  
  # Monta plot interativo
  p_interativo <- plot_ly(
    data = df_mediana,
    x = ~Ciclo,
    y = ~Mediana,
    split = ~Arquivo,
    type = 'scatter',
    mode = 'lines',
    color = ~Arquivo,  
    colors = cores,  
    hoverinfo = 'text',
    text = ~paste("Amostra:", Arquivo,
                  "<br>Ciclo:", Ciclo,
                  "<br>Score Mediano:", round(Mediana, 2))
  ) %>% plotly::layout(
    shapes = list(
      list(type = "rect", x0 = 0, x1 = x_max, y0 = 0, y1 = 20,
           fillcolor = "rgba(255, 0, 0, 0.05)", line = list(width = 0), layer = "below"),
      list(type = "rect", x0 = 0, x1 = x_max, y0 = 20, y1 = 28,
           fillcolor = "rgba(255, 255, 0, 0.05)", line = list(width = 0), layer = "below"),
      list(type = "rect", x0 = 0, x1 = x_max, y0 = 28, y1 = 40,
           fillcolor = "rgba(0, 255, 0, 0.05)", line = list(width = 0), layer = "below")
    ),
    plot_bgcolor = "white",  
    paper_bgcolor = "white",  
    title = list(
      text = "Qualidade das Sequências por Ciclo",
      font = list(
        size = 16,
        color = "black"
      ),
      xanchor = "left",
      x = 0.07,
      y = 0.98
    ),
    xaxis = list(
      title = "Ciclo",
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
      title = "Score de Qualidade (Phred)",
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
    legend = list(
      title = list(text = "Arquivo"),
      font = list(
        color = "black"  
      ),
      x = 1,        
      y = 0.5,     
      xanchor = "left", 
      yanchor = "middle" 
    )
  )
  
  # Monta plot estático
  p_estatico <- ggplot(df_mediana, aes(x = Ciclo, y = Mediana, group = Arquivo)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 28, ymax = 42),
              fill = "#d0f0d0", alpha = 0.015) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 28),
              fill = "#f0f0d0", alpha = 0.015) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 20),
              fill = "#f0d0d0", alpha = 0.015) +
    geom_line(aes(color = Arquivo), linewidth = 0.7, show.legend = TRUE) +
    labs(
      title = "Qualidade das Sequências por Ciclo",
      x = "Ciclo",
      y = "Score de Qualidade (Phred)"
    ) +
    scale_y_continuous(limits = c(0, 42), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = cores) +
    theme(legend.position = "right",
          legend.text = element_text(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "grey"))
  
  # Retorna lista com os dois gráficos (estático e interativo)
  lista_plots <- list(p_interativo = p_interativo, p_estatico = p_estatico)
  
  return(lista_plots)
}
