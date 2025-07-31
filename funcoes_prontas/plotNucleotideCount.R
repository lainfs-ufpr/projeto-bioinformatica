# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("tidyverse", "dplyr", "ggplot2", "tidyr", "scales")
bioc_pkgs <- c("ShortRead")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(stats)

load("cores_lainfs.RData")

# -------------------------------------------
# Função para verificar a frequência das bases por read
# Param -  caminho_fastq: caminho para o arquivo .fastq
# Return - p_contagens: ggplot com boxplot para contagem de bases
# -------------------------------------------
plotNucleotideCount <- function(caminho_fastq) {
  
  # É um único arquivo, um plot só
  if (is.character(caminho_fastq) && length(caminho_fastq) == 1) {
    # Ler os reads
    fq <- readFastq(caminho_fastq)
    reads <- as.character(sread(fq))
    
    # Contar bases por read
    contagens <- lapply(reads, function(seq) {
      table(factor(strsplit(seq, "")[[1]], levels = c("A", "T", "C", "G", "N")))
    })
    
    # Converter em data.frame
    df_contagens <- do.call(rbind, contagens)
    rownames(df_contagens) <- paste0("Read", seq_along(reads))
    df_contagens <- as.data.frame(df_contagens)
    
    # Transformar para formato longo e calcular frequência por read
    df_frequencias <- df_contagens %>%
      rownames_to_column(var = "Sample") %>%
      pivot_longer(cols = A:N, names_to = "Base", values_to = "Count") %>%
      mutate(Base = factor(Base, levels = c("A", "T", "C", "G", "N"))) %>% 
      group_by(Sample) %>%
      mutate(Frequency = Count / sum(Count)) %>%
      ungroup()
    
    df_frequencias$Arquivo <- rep(basename(caminho_fastq), nrow(df_frequencias))
    
    p_contagens <- ggplot(df_frequencias, aes(x = Base, y = Frequency, fill = Base)) +
      geom_boxplot(alpha = 0.9, outlier.shape = NA,
                   color = "black") + 
      geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1)) +
      scale_fill_manual(values = cores_lainfs[1:length(unique(df_frequencias$Base))]) +
      labs(
        title = "Distribuição da Frequência de Bases por Read",
        x = "Base",
        y = "Frequência (pseudo-log)"
      ) +
      theme(legend.position = "none",
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "grey"))
  } else {
    # vetor de arquivos
    lista_dfs <- lapply(caminho_fastq, function(arquivo) {
      fq <- readFastq(arquivo)
      reads <- as.character(sread(fq))
      
      contagens <- lapply(reads, function(seq) {
        table(factor(strsplit(seq, "")[[1]], levels = c("A", "T", "C", "G", "N")))
      })
      
      df_contagens <- do.call(rbind, contagens)
      df_contagens <- as.data.frame(df_contagens)
      
      # Frequência por read
      df_frequencias <- df_contagens %>%
        pivot_longer(cols = A:N, names_to = "Base", values_to = "Count") %>%
        group_by(Base) %>%
        summarise(MeanFrequency = mean(Count / rowSums(df_contagens))) %>%
        mutate(Arquivo = basename(arquivo))
      
      return(df_frequencias)
    })
    
    df_todos <- do.call(rbind, lista_dfs)
    df_todos$Base <- factor(df_todos$Base, levels = c("A", "T", "C", "G", "N"))
    
    p_contagens <- ggplot(df_todos, aes(x = Base, y = MeanFrequency,
                                        color = Arquivo, group = Arquivo)) +
      geom_point(size = 3, position = position_dodge(width = 0.4)) +
      scale_color_manual(values = rep_len(cores_lainfs, length(unique(df_todos$Arquivo)))) +
      labs(
        title = "Média da Frequência de Bases",
        x = "Base",
        y = "Frequência Média"
      ) +
      theme(legend.position = "right",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "grey"))
  }
  
  return(p_contagens)
}