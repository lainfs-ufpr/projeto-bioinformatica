# -------------------------------------------
# Carregamento das Bibliotecas
# -------------------------------------------
source("instalar_pacotes.R")
cran_pkgs <- c("ggplot2")
bioc_pkgs <- c("ShortRead")

instalar_pacotes(cran_pkgs, install.packages)
instalar_pacotes(bioc_pkgs, BiocManager::install)

library(ShortRead)
library(ggplot2)

# -------------------------------------------
# Função para gerar gráfico da qualidade
# Param -  qa_output: output da funcao qa()
# Return - qualidade_p: ggplot
# -------------------------------------------
readQualityScore <- function(qa_output) {
  
  # Pega scores de qualidade
  df <- qa_output[["readQualityScore"]]
  
  # Filtrar apenas o tipo "read" (qualidade média por leitura)
  df_filtrada <- subset(df, type == "read")
  
  # Gerar gráfico com ggplot2
  p_qualidade <- ggplot(df_filtrada,
                      aes(x = quality, y = density, fill = lane, color = lane)) +
                geom_area(alpha = 0.3, position = "identity") +
                geom_line(linewidth = 1) +
                theme_minimal() +
                labs(
                  title = "Distribuição da Qualidade Média por Leitura - por Amostra",
                  x = "Score médio (Phred)",
                  y = "Densidade relativa",
                  fill = "Amostra (lane)",
                  color = "Amostra (lane)"
                )
  
  return(p_qualidade)
}
