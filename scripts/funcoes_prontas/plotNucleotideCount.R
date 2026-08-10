# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("scripts/instalar_pacotes.R")
cran_pkgs <- c("tidyverse", "dplyr", "ggplot2", "tidyr", "scales", "plotly")
instalar_pacotes(cran_pkgs, install.packages)

library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(plotly)

# -------------------------------------------
# Função para verificar a frequência das bases por read
#
# Reaproveita os dados já calculados pelo ShortRead::qa() (componente
# perCycle$baseCall), evitando reabrir/reler os arquivos FASTQ do disco.
#
# Param -  qa_output: objeto retornado por ShortRead::qa()
# Param -  paleta_cores: string com nome da paleta escolhida pelo usuario
# Param -  nomes_arquivos: nomes originais dos arquivos de entrada
# Return - lista_plots: lista nomeada com plot estático e plot interativo
# -------------------------------------------
plotNucleotideCount <- function(qa_output, paleta_cores = "viridis", nomes_arquivos) {
  
  # perCycle$baseCall: Cycle, Base, Count, lane (lane = nome do arquivo)
  df_base <- qa_output[["perCycle"]][["baseCall"]]
  
  df_base <- df_base %>%
    dplyr::filter(Base %in% c("A", "T", "C", "G", "N")) %>%
    dplyr::mutate(Base = factor(Base, levels = c("A", "T", "C", "G", "N")),
                  Arquivo = as.character(lane))
  
  # Frequência de cada base dentro de cada ciclo (posição) e arquivo
  df_frequencias <- df_base %>%
    dplyr::group_by(Arquivo, Cycle) %>%
    dplyr::mutate(Frequency = Count / sum(Count)) %>%
    dplyr::ungroup()
  
  # Ajusta nomes dos arquivos
  df_frequencias$Arquivo <- nomes_arquivos[
    match(df_frequencias$Arquivo,
          unique(df_frequencias$Arquivo))
  ]
  
  n_arquivos <- length(unique(df_frequencias$Arquivo))
  
  # Um arquivo: distribuição da frequência de cada base ao longo dos ciclos
  if (n_arquivos == 1) {
    
    n_cores <- length(levels(df_frequencias$Base))
    cores <- viridis::viridis(n_cores, option = paleta_cores)
    if (length(cores) < n_cores) {
      cores <- rep(cores, length.out = n_cores)
    }
    
    p_estatico <- ggplot2::ggplot(df_frequencias, ggplot2::aes(x = Base,
                                                               y = Frequency,
                                                               fill = Base)) +
      ggplot2::geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
      ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
      ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
      ggplot2::scale_fill_manual(values = cores) +
      ggplot2::labs(
        title = "Distribuição da Frequência de Bases por Ciclo",
        x = "Base",
        y = "Frequência (pseudo-log)"
      ) +
      ggplot2::theme(
        legend.position = "none",
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "grey")
      )
    
    # Mais de um arquivo: frequência média de cada base por arquivo
  } else {
    
    df_todos <- df_frequencias %>%
      dplyr::group_by(Arquivo, Base) %>%
      dplyr::summarise(MeanFrequency = mean(Frequency), .groups = "drop")
    
    cores <- viridis::viridis(n_arquivos, option = paleta_cores)
    
    p_estatico <- ggplot(df_todos, aes(x = Base, y = MeanFrequency,
                                       color = Arquivo, group = Arquivo)) +
      geom_point(size = 2, alpha = 0.5, position = position_dodge(width = 0.4)) +
      scale_color_manual(values = cores) +
      labs(
        title = "Média da Frequência de Bases",
        x = "Base",
        y = "Frequência Média"
      ) +
      theme(plot.title = element_text(size = 20, color = "black"),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16), legend.position = "right",
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "grey"))
  }
  
  # Cria gráfico interativo
  p_interativo <- ggplotly(p_estatico)
  p_interativo <- p_interativo %>%
    layout(
      title = list(
        font = list(
          size = 20,
          color = "black"
        )
      ),
      xaxis = list(
        titlefont = list(
          size = 16,
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
          size = 16,
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