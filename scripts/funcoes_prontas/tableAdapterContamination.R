# -------------------------------------------
# Carregamento das bibliotecas, funções e objetos necessários
# -------------------------------------------
source("scripts/instalar_pacotes.R")
cran_pkgs <- c("reshape2", "dplyr")
bioc_pkgs <- c("ShortRead", "Biostrings")
instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(Biostrings)    
library(reshape2)  
library(dplyr)          

# -------------------------------------------
# Função para verificar a existência de adaptadores
#
# Como o qa() do ShortRead não guarda as sequências das reads (só
# estatísticas agregadas), essa checagem precisa reler os arquivos.
# Para arquivos grandes, isso pode ser bem mais rápido usando uma
# amostra de reads (como o FastQC faz) em vez do arquivo inteiro.
#
# Param -  caminho_fastq: caminho (ou vetor de caminhos) para arquivo(s) .fastq
# Param -  amostragem: se TRUE (padrão), analisa uma amostra aleatória de reads
#                       via FastqSampler; se FALSE, lê o arquivo inteiro
# Param -  n_amostra: número de reads da amostra quando amostragem = TRUE
#                      (se o arquivo tiver menos reads que isso, lê todas)
# Return - df_adapters_final: data.frame com Arquivo e Adaptador(es) encontrados
# -------------------------------------------
tableAdapterContamination <- function(caminho_fastq, amostragem = TRUE, n_amostra = 100000) {
  
  # Ler arquivo de adaptadores (uma única vez, fora do loop de arquivos)
  adapter_fasta <- readDNAStringSet("data/adapters/adapters.fasta")
  
  constroi_df <- function(caminho_fastq) {
    
    if (amostragem) {
      # Amostra aleatória de reads, sem precisar ler o arquivo inteiro.
      # Se o arquivo tiver menos reads que n_amostra, retorna todas.
      sampler <- ShortRead::FastqSampler(caminho_fastq, n = n_amostra)
      fq <- ShortRead::yield(sampler)
      close(sampler)
      reads <- ShortRead::sread(fq)
      
    } else {
      # Lê o arquivo inteiro, em blocos, até acabar (loop corrigido:
      # antes só se lia o primeiro bloco do FastqStreamer).
      stream <- ShortRead::FastqStreamer(caminho_fastq)
      reads_list <- list()
      repeat {
        chunk <- ShortRead::yield(stream)
        if (length(chunk) == 0) break
        reads_list[[length(reads_list) + 1]] <- ShortRead::sread(chunk)
      }
      close(stream)
      reads <- do.call(c, reads_list)
    }
    
    # Para cada adaptador, verifica se ocorre em pelo menos uma read da amostra/arquivo
    presenca_adapters <- sapply(adapter_fasta, function(adapter_seq) {
      any(Biostrings::vcountPattern(adapter_seq, reads) > 0)
    })
    
    if (!any(presenca_adapters)) {
      adapters_presentes <- "Sem hits"
    } else {
      adapters_presentes <- names(adapter_fasta)[presenca_adapters]
    }
    
    # Montar a linha da tabela
    data.frame(
      Arquivo = basename(caminho_fastq),
      Adaptador = paste(adapters_presentes, collapse = ", ")
    )
  }
  
  # Constrói a tabela para um ou vários arquivos (evita crescer data.frame em loop)
  lista_dfs <- lapply(caminho_fastq, constroi_df)
  df_adapters_final <- do.call(rbind, lista_dfs)
  
  return(df_adapters_final)
}