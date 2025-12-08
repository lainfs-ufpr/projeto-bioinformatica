# -------------------------------------------
# Carregamento das Bibliotecas
# -------------------------------------------
source("scripts/instalar_pacotes.R")
bioc_pkgs <- c("ShortRead", "Biostrings")

instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(Biostrings)
library(stats)

# -------------------------------------------
# Função para verificar as sequências mais frequentes e adaptadores
# Param -  qa_output: output da funcao qa()
#          read: tipo da leitura, para buscar em df$type
#          n_frequentes: número de sequências mais frequentes para retornar
# Return - freqSequences_table: tabela com sequências mais frequentes
# -------------------------------------------
freqSequences <- function(qa_output, n_frequentes = 20) {
  
  # Ler arquivo de adapters
  adapter_fasta <- readDNAStringSet("data/adapters.fasta")
  
  # Contagens de todos os arquivos
  all_counts <- qa_output[["readCounts"]]

  # Sequências frequentes (de todos os arquivos) com frequencia > 1
  df <- qa_output[["frequentSequences"]] %>%
          filter(count > 1)

  resultados <- list()
  # Itera sobre os arquivos (ou lanes)
  for (lane_name in rownames(all_counts)) {
    
    # Filtra apenas as sequências da lane atual
    df_filtrada <- df[df$type == "read" & df$lane == lane_name, ]
    
    # Pega total de reads dessa lane
    total_reads <- all_counts[lane_name, "read"]
    
    # Calcula a proporção
    df_filtrada[["ppn"]] <- df_filtrada[["count"]] / total_reads
    
    if (nrow(df_filtrada) == 0) {
      # Adiciona entrada informando ausência de sequências frequentes
      resultados[[lane_name]] <- data.frame(
        Sequência = NA,
        Frequência = NA,
        Adaptador = "Sem seqs com freq > 1",
        Arquivo = lane_name,
        stringsAsFactors = FALSE
      )
      next
    }
    # Seleciona as mais frequentes
    top_df <- head(df_filtrada[order(df_filtrada$count, decreasing = TRUE),
                               c("sequence", "count")],
                   n_frequentes)
    
    # Verifica se são adaptadores conhecidos
    top_df$possible_hit <- sapply(top_df$sequence, function(seq) {
      matches <- vcountPattern(DNAString(seq), adapter_fasta, fixed = FALSE)
      if (any(matches > 0)) {
        names(adapter_fasta)[which.max(matches)]
      } else {
        "Sem hit"
      }
    })
    
    # Adiciona coluna do nome do arquivo (ou lane)
    top_df$Arquivo <- lane_name
    colnames(top_df) <- c("Sequência", "Frequência", "Adaptador", "Arquivo")
    
    # Adiciona aos resultados
    resultados[[lane_name]] <- top_df
    
  }

  # Junta tudo em uma única tabela
  df_final <- do.call(rbind, resultados)
  
  # Organiza nomes
  rownames(df_final) <- NULL

  # Reordena colunas para deixar "Arquivo" primeiro
  df_final <- df_final[, c("Arquivo", "Sequência", "Frequência", "Adaptador")]
  df_final$Arquivo <- as.character(df_final$Arquivo)

  return(df_final)
}
