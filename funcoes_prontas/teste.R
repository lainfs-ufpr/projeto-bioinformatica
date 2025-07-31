source("funcoes_prontas/plotAdapterContamination.R")
source("funcoes_prontas/tableAdapterContamination.R")
source("funcoes_prontas/freqSequences.R")
source("funcoes_prontas/plotCycleQuality.R")
source("funcoes_prontas/plotNucleotideCount.R")
source("funcoes_prontas/readQualityScore.R")
source("funcoes_prontas/plotOcurrences.R")

# Carrega arquivos fastq

caminho_fastq <- "teste_funcoes/adapterContamination/mock.fastq"
caminho_fastq2 <- "teste_funcoes/adapterContamination/mock2.fastq"
caminho_fastq3 <- "teste_funcoes/adapterContamination/mock3.fastq"
caminho_fastq4 <- "fastq_gerados/output1.fastq"
caminho_fastq5 <- "fastq_gerados/output51.fastq"

resultado_qa <- qa(caminho_fastq5, type = "fastq")
resultado_qa2 <- qa(c(caminho_fastq, caminho_fastq2), type = "fastq")
resultado_qa3 <- qa(c(caminho_fastq, caminho_fastq2, caminho_fastq3), type = "fastq")

# Testa funcoes

plot_adapter <- plotAdapterContamination(c(caminho_fastq, caminho_fastq2, caminho_fastq4))
tabela_adapter <- tableAdapterContamination(c(caminho_fastq, caminho_fastq2, caminho_fastq4))

# ---

tabela_frequencias <- freqSequences(resultado_qa)
tabela_frequencias <- freqSequences(resultado_qa2)

# ---

plots_qualidade_ciclo <- plotCycleQuality(resultado_qa2)
qualidade_ciclo_interativo <- plots_qualidade_ciclo$p_interativo
qualidade_ciclo_estatico <- plots_qualidade_ciclo$p_estatico

# ---

plot_contagens <- plotNucleotideCount(caminho_fastq)

plot_contagens <- plotNucleotideCount(c(caminho_fastq, caminho_fastq2, caminho_fastq3))

# ---
namess <- c("oi", "tchau")
plot_qualidade_media <- readQualityScore(resultado_qa2)
plot_qualidade_media$data$Arquivo <- namess[
  match(plot_qualidade_media$data$Arquivo,
        unique(plot_qualidade_media$data$Arquivo))
]

# ---

plot_ocorrencias <- plotOcurrences(resultado_qa2)
plot_ocorrencias <- plotOcurrences(resultado_qa3)

