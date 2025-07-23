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
#          cor_linhas: cor para linhas principais do plot
# Return - lista com plot interativo e plot estático
# -------------------------------------------

plotCycleQuality <- function(qa_output, cor_linhas="royalblue") {
  
  # Pegar qualidade por ciclo sem processamento
  df <- qa_output[["perCycle"]][["quality"]]
  
  # Calcular mediana por lane e ciclo
  df_median <- df %>%
    group_by(lane, Cycle) %>%
    summarise(
      Median = quantile(rep(Score, Count), 0.5),
      .groups = "drop"
    )
  
  # Montar plot interativo
  x_max <- max(df_median$Cycle)
  p_interativo <- plot_ly(
    data = df_median,
    x = ~Cycle,
    y = ~Median,
    split = ~lane, # Dividir linhas por 'lane'
    type = 'scatter',
    mode = 'lines',
    color = I(cor_linhas),
    hoverinfo = 'text', # Informações que aparecem ao passar o mouse
    text = ~paste("Amostra:", lane, # Texto exibido ao passar o mouse
                  "<br>Ciclo:", Cycle,
                  "<br>Score Mediano:", round(Median, 2))
  ) %>% layout(
    # Adicionar retângulos de fundo para indicar qualidade
    shapes = list(
      # Faixa vermelha clara: qualidade ruim (0–20)
      list(type = "rect", x0 = 0, x1 = x_max, y0 = 0, y1 = 20,
           fillcolor = "rgba(255, 0, 0, 0.1)",
           line = list(width = 0), layer = "below"),
      # Faixa amarela clara: qualidade intermediária (20–28)
      list(type = "rect", x0 = 0, x1 = x_max, y0 = 20, y1 = 28,
           fillcolor = "rgba(255, 255, 0, 0.1)",
           line = list(width = 0), layer = "below"),
      # Faixa verde clara: boa qualidade (28–40)
      list(type = "rect", x0 = 0, x1 = x_max, y0 = 28, y1 = 40,
           fillcolor = "rgba(0, 255, 0, 0.1)",
           line = list(width = 0), layer = "below")
    ),
    title = "Qualidade das Sequências por Ciclo",
    xaxis = list(title = "Ciclo da Sequência"),
    yaxis = list(title = "Score de Qualidade (Phred)", range = c(0, 42)),
    showlegend = FALSE
  )
  
  # Montar plot estático
  p_estatico <- ggplot(df_median, aes(x = Cycle, y = Median, group = lane)) +
    
    # Camada 1: Retângulos de fundo para indicar a qualidade
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 28, ymax = 42),
              fill = "#d0f0d0", alpha = 0.5) + # Verde
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 20, ymax = 28),
              fill = "#f0f0d0", alpha = 0.5) + # Amarelo
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = 20),
              fill = "#f0d0d0", alpha = 0.5) +   # Vermelho
    
    # Camada 2: Linhas de qualidade para cada amostra
    # O 'group = lane' garante que cada amostra tenha sua própria linha.
    geom_line(color = cor_linhas, linewidth = 0.7) +
    
    # Camada 3: Títulos, legendas e tema visual
    labs(
      title = "Qualidade das Sequências por Ciclo",
      x = "Ciclo da Sequência",
      y = "Score de Qualidade (Phred)"
    ) +
    scale_y_continuous(limits = c(0, 42), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_bw()
  
  return(list(p_interativo=p_interativo, 
              p_estatico=p_estatico))
}
