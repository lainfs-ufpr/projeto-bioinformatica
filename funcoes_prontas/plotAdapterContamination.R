# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("reshape2", "dplyr", "ggplot2", "scales")
bioc_pkgs <- c("ShortRead", "Biostrings")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(Biostrings)    
library(reshape2)  
library(dplyr)         
library(ggplot2)      
library(scales)        

load("cores_lainfs.RData")

# -------------------------------------------
# Função para verificar a existência de adaptadores
# Param -  caminho_fastq: caminho para o arquivo .fastq
# Return - p_adapter: plot com proporção cumulativa de adaptadores
# -------------------------------------------
plotAdapterContamination <- function(caminho_fastq){
  
  constroi_df <- function(caminho_fastq){
    # Ler sequências do fastq
    fq <- readFastq(caminho_fastq)
    reads <- sread(fq)
    
    # Para cada adaptador, um vetor com 0 ou 1 por read
    adapter_matriz <- sapply(adapter_fasta, function(adapter_seq) {
      as.integer(vcountPattern(adapter_seq, reads) > 0)
    })
    
    # Adicionar index dos reads
    adapter_df <- data.frame(Read = 1:length(reads), adapter_matriz)
    
    # Tamanho das reads
    tamanho_reads <- nchar(as.character(sread(fq)))  
    
    # Final de cada read (posição cumulativa)
    final_reads <- cumsum(tamanho_reads)
    
    # Início de cada read
    inicio_reads <- c(1, head(final_reads + 1, -1))
    
    # Passar para df
    adapter_df$Inicio <- inicio_reads
    adapter_df$Fim <- final_reads
    

    
    # Reorganizar o df de formato wide para long: uma linha por read por adaptador
    adapter_df_long <- melt(adapter_df, id.vars = c("Inicio", "Fim", "Read"),
                            variable.name = "Adaptador", value.name = "Hit")
    
    # Soma cumulativa para proporção dos adaptadores
    adapter_df_long <- adapter_df_long %>%
      dplyr::group_by(Adaptador) %>%
      dplyr::arrange(Inicio) %>%
      dplyr::mutate(
        CumHits = cumsum(Hit),
        Prop = CumHits / max(Read)
      )
    
    return(adapter_df_long)
  }
  
  # Ler arquivo de adaptadores
  adapter_fasta <- readDNAStringSet("adapters.fasta")

  # É um único arquivo, um plot só
  if (is.character(caminho_fastq) && length(caminho_fastq) == 1) {
    
    adapter_df_long <- constroi_df(caminho_fastq)
    
    # Separar entre adaptadores com hit e sem hit
    adapter_hit_info <- adapter_df_long %>%
      dplyr::group_by(Adaptador) %>%
      dplyr::summarise(has_hits = any(Hit > 0), .groups = "drop")
    
    # Com hit
    hit_adapters <- adapter_hit_info$Adaptador[adapter_hit_info$has_hits]
    adapter_com_hits <- adapter_df_long %>%
      dplyr::filter(Adaptador %in% hit_adapters)
    adapter_com_hits$Adaptador <- factor(adapter_com_hits$Adaptador)
    
    # Sem hit
    adapter_sem_hits <- adapter_df_long %>%
      dplyr::filter(!Adaptador %in% hit_adapters)
    
    # Checa se houve hits
    houve_hits <- any(adapter_com_hits$Prop > 0)
    
    if (houve_hits) {
      # Plot com os adaptadores
      max_y <- max(adapter_com_hits$Prop)
      top_break <- ceiling(max_y * 20) / 20
      breaks_y <- seq(0, top_break, by = 0.05)
      
      # Fazer plot
      p_adapter <- ggplot() +
        geom_line(data = adapter_com_hits,
                  aes(x = Inicio, y = Prop, group = Adaptador, color = Adaptador),
                  alpha = 0.8, linewidth = 1, show.legend = TRUE) +
        scale_y_continuous(
          limits = c(0, top_break),
          breaks = breaks_y,
          labels = scales::percent_format(accuracy = 1)
        )
      
    } else {
      # Plot sem hits
      p_adapter <- ggplot() +
        geom_line(data = adapter_sem_hits,
                  aes(x = Inicio, y = 0, group = Adaptador),
                  color = "grey70", alpha = 0.5, linewidth = 0.5, show.legend = FALSE) +
        scale_y_continuous(
          limits = c(0, NA),
          breaks = 0,       
          labels = scales::percent_format(accuracy = 1)
        )
    }
    
    # Adiciona os elementos comuns ao plot
    p_adapter <- p_adapter +
      labs(
        title = if (houve_hits) "Contaminação Cumulativa de Adaptadores" else 
          "Contaminação Cumulativa de Adaptadores (sem hits detectados)",
        x = "Posição (par de bases)",
        y = "Proporção de reads contaminados"
      ) +
      scale_color_manual(values = cores_lainfs) +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "grey")
      )
  # Mais de um arquivo, um plot com cada
  } else if (is.character(caminho_fastq) && length(caminho_fastq) > 1) {
    
    p_adapter <- ggplot()
    
    for (fastq in caminho_fastq) {
      adapter_df_long <- constroi_df(fastq)
      
      # Separar entre adaptadores com hit e sem hit
      adapter_hit_info <- adapter_df_long %>%
        dplyr::group_by(Adaptador) %>%
        dplyr::summarise(has_hits = any(Hit > 0), .groups = "drop")
      
      # Com hit
      hit_adapters <- adapter_hit_info$Adaptador[adapter_hit_info$has_hits]
      adapter_com_hits <- adapter_df_long %>%
        dplyr::filter(Adaptador %in% hit_adapters)
      adapter_com_hits$Adaptador <- factor(adapter_com_hits$Adaptador)
      
      # Sem hit
      adapter_sem_hits <- adapter_df_long %>%
        dplyr::filter(!Adaptador %in% hit_adapters)
      
      # Checa se houve hits
      houve_hits <- any(adapter_com_hits$Prop > 0)
      
      if (houve_hits) {
        # Fazer plot
        p_adapter <- p_adapter +
          geom_line(data = adapter_com_hits,
                    aes(x = Inicio, y = Prop,
                        group = Adaptador, color = Adaptador),
                    alpha = 0.3, linewidth = 2, show.legend = TRUE) +
          geom_point(data = adapter_com_hits,
                    aes(x = Inicio, y = Prop,
                        color = Adaptador,shape = Adaptador),
                    size = 2, alpha = 0.5)
      } 
    }
    # Adiciona os elementos comuns ao plot
    p_adapter <- p_adapter +
      labs(
        title = "Contaminação Cumulativa de Adaptadores",
        x = "Posição (par de bases)",
        y = "Proporção de reads contaminados"
      ) +
      scale_color_manual(values = cores_lainfs) +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "grey")
      )
  }
  return(p_adapter)
}
