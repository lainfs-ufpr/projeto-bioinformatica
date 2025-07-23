# -------------------------------------------
# Carregamento das Bibliotecas
# -------------------------------------------
source("instalar_pacotes.R")
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
freqSequences <- function(qa_output, read="read", n_frequentes=20) {
  
  # Ler arquivo de adapters
  adapter_fasta <- readDNAStringSet("adapters.fasta")
  
  # Contagens
  cnt <- qa_output[["readCounts"]]
  
  # Frequencia das sequencias
  df <- qa_output[["frequentSequences"]]
  
  # Filtra as sequencias do tipo especificado (por padrao: "read")
  df_filtrada <- df[df$type==read,]
  
  # Calcula a proporcao da sequencia em relacao ao total de leituras 
  df_filtrada[["ppn"]] <- df_filtrada[["count"]] / cnt[df_filtrada[["lane"]], read]
  
  # Seleciona as sequencias mais frequentes (por padrao: 20), com colunas especificas
  freqSequences_df <- head(df_filtrada[order(df_filtrada$count,
                                             decreasing=TRUE),
                                             c("sequence", "count")],
                           n_frequentes)
  
  # Remove nome das linhas pra simplificar
  rownames(freqSequences_df) <- NULL
  
  # Verifica se sao adapters conhecidos
  freqSequences_df$possible_hit <- sapply(freqSequences_df$sequence, function(seq) {
    matches <- vcountPattern(DNAString(seq), adapter_fasta, fixed=FALSE)
    if (any(matches > 0)) {
      # Retorna o nome do adapter que tem mais hits
      names(adapter_fasta)[which.max(matches)] 
    } else {
      "no hit"
    }
  })
  
  freqSequences_tabela <- table(freqSequences_df)
  
  return(freqSequences_tabela)
}