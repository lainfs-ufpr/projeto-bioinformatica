# -------------------------------------------
# Carregamento das Bibliotecas
# -------------------------------------------
source("instalar_pacotes.R")
bioc_pkgs <- c("ShortRead", "tidyverse")

instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)    
library(tidyverse)
library(stats)

# -------------------------------------------
# Função para verificar as sequências mais frequentes e adaptadores
# Param -  qa_output: output da funcao qa()
#          read: tipo da leitura, para buscar em df$type
#          n_frequentes: número de sequências mais frequentes para retornar
# Return - freqSequences_table: tabela com sequências mais frequentes
# -------------------------------------------
plotNucleotideCount <- function(qa_output) {
  
  # Pegar a contagem de cada base
  contagem_bases <- qa_output[["baseCalls"]]
  
  df <- contagem_bases %>%
    as.data.frame() %>%         
    
    # Transforma nomes das linhas em uma coluna 'Sample'
    rownames_to_column(var = "Sample") %>%  
    
    # Transforma as colunas em formato longo e renomeia colunas
    pivot_longer(cols = A:N,                    
                 names_to = "Base",             
                 values_to = "Count") %>%      
    group_by(Sample) %>%   
    
    # Calcula a frequência relativa de cada base
    mutate(Frequency = Count / sum(Count)) %>%   
    ungroup()   
  
  p_contagens <- ggplot(df, aes(x = Base, y = Frequency, fill = Base)) +
                  geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
                  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  
                  scale_y_continuous(trans = pseudo_log_trans(sigma = 0.5)) +
                  scale_fill_brewer(palette = "Set1") +  
                  labs(
                    title = "Distribuição da Frequência de Bases por Amostra (Escala Pseudo-Log)",
                    x = "Base",
                    y = "Frequência (pseudo-log)"
                  ) +
                  theme_minimal(base_size = 12) +
                  theme(legend.position = "none")
  
  return(p_contagens)
}