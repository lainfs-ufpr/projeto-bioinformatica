# Função para gerar sequências de RNA e qualidades aleatórias em formato FASTQ
generate_random_fastq <- function(num_sequences = 10, seq_length = 20, file_name = "random_output.fastq") {
  
  # Gerar sequências aleatórias de RNA
  generate_rna_sequence <- function(length) {
    nucleotides <- c("A", "U", "G", "C")
    return(paste0(sample(nucleotides, length, replace = TRUE), collapse = ""))
  }
  
  # Gerar sequências aleatórias de qualidades no formato Phred (ASCII 33-73)
  generate_phred_quality <- function(length) {
    quality_scores <- sample(33:73, length, replace = TRUE)
    return(intToUtf8(quality_scores))  # Converter para caracteres ASCII
  }
  
  # Abrir o arquivo para escrita
  file_conn <- file(file_name, "w")
  
  # Gerar e escrever cada sequência e qualidade no arquivo FASTQ
  for (i in 1:num_sequences) {
    # Identificador da sequência
    identifier <- paste0("@seq", i)
    
    # Gerar uma sequência de RNA aleatória
    sequence <- generate_rna_sequence(seq_length)
    
    # Gerar uma sequência de qualidade aleatória
    quality <- generate_phred_quality(seq_length)
    
    # Escrever no arquivo FASTQ
    writeLines(c(identifier, sequence, "+", quality), file_conn)
  }
  
  # Fechar o arquivo
  close(file_conn)
  
  cat("Arquivo FASTQ aleatório gerado com sucesso:", file_name, "\n")
}

# Exemplo de uso:


for (i in 1:50){
  generate_random_fastq(num_sequences = 100, seq_length = 200, 
                        file_name = paste0("fastq_gerados/sequencia_gerada", i, ".fastq"))
}

