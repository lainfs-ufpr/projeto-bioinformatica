# -------------------------------------------
# Carregamento das Bibliotecas
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("ggplot2", "scales", "dplyr", "reshape2")
bioc_pkgs <- c("ShortRead", "Biostrings")

instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(Biostrings)    
library(reshape2)  
library(dplyr)         
library(ggplot2)      
library(scales)

# -------------------------------------------
# Função para verificar a presença de adaptadores e gerar plot
# Param -  arquivo_fastq: caminho do arquivo fastq 
# Return - adapter_p: ggplot dos adapters encontrados
# -------------------------------------------
adapterContamination <- function(arquivo_fastq) {
  
  # Ler arquivo de adapters
  adapter_fasta <- readDNAStringSet("adapters.fasta")
  
  # Ler sequencias do fastq
  arquivo_fastq_seqs <- sread(arquivo_fastq)
  
  # Para cada adapter, um vetor com 0 ou 1 por read
  adapter_matriz <- sapply(adapter_fasta, function(adapter_seq) {
    as.integer(vcountPattern(adapter_seq, arquivo_fastq_seqs) > 0)
  })
  
  # Adicionar index dos reads
  adapter_df <- data.frame(Read = 1:length(arquivo_fastq_seqs), adapter_matriz)
  
  # Tamanho das reads
  tamanho_reads <- nchar(as.character(sread(arquivo_fastq)))  
  
  # Posicao inicial e final de cada read  
  final_reads <- cumsum(tamanho_reads)
  inicio_reads <- c(1, head(final_reads + 1, -1))
  
  # Passar para df
  adapter_df$Inicio <- inicio_reads
  adapter_df$Fim <- final_reads
  
  # Reorganizar o df de formato wide para long: uma linha por read por adapter
  adapter_df_long <- melt(adapter_df, id.vars = c("Inicio", "Fim", "Read"),
                          variable.name = "Adapter", value.name = "Hit")
  
  # Soma cumulativa para proporca dos adapters
  adapter_df_long <- adapter_df_long %>%
    group_by(Adapter) %>%
    arrange(Inicio) %>%
    mutate(
      CumHits = cumsum(Hit),
      Prop = CumHits / max(Read)
    )
  
  # Separar entre adapters com hit e sem hit
  adapter_hit_info <- adapter_df_long %>%
    group_by(Adapter) %>%
    summarise(has_hits = any(Hit > 0), .groups = "drop")
  
  # Com hit
  hit_adapters <- adapter_hit_info$Adapter[adapter_hit_info$has_hits]
  adapter_com_hits <- adapter_df_long %>%
    filter(Adapter %in% hit_adapters)
  adapter_com_hits$Adapter <- factor(adapter_com_hits$Adapter)
  
  # Sem hit
  adapter_sem_hits <- adapter_df_long %>%
    filter(!Adapter %in% hit_adapters)
  
  # Checa se houve hits
  # Se nao, usa uma escala menor para o grafico
  houve_hits <- any(adapter_com_hits$Prop > 0)
  if (houve_hits) {
    max_y <- max(adapter_com_hits$Prop)
    top_break <- ceiling(max_y * 20) / 20  # arredondar
    breaks_y <- seq(0, top_break, by = 0.05)
    
    y_scale <- scale_y_continuous(
      limits = c(0, top_break),
      breaks = breaks_y,
      labels = scales::percent_format(accuracy = 1)
    )
  } else {
    y_scale <- scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0),  
      labels = scales::percent_format(accuracy = 1)
    )
  }
  
  # Plot
  adapter_p <- ggplot() +
    geom_line(data = adapter_sem_hits,
              aes(x = Inicio, y = Prop, group = Adapter),
              color = "grey70", alpha = 0.5, linewidth = 0.5, show.legend = FALSE) +
    geom_line(data = adapter_com_hits,
              aes(x = Inicio, y = Prop, group = Adapter, color = Adapter),
              alpha = 0.8, linewidth = 1, show.legend = TRUE) +
    # Apply to y-axis
    y_scale +
    labs(title = "Contaminação por adaptador cumulativa nos reads",
         x = "Posição do par de bases",
         y = "Proporção de reads contaminados") +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_blank()
    )
  
  return(adapter_p)
}