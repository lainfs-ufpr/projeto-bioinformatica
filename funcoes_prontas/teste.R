source("funcoes_prontas/plotAdapterContamination.R")
source("funcoes_prontas/tableAdapterContamination.R")
source("funcoes_prontas/freqSequences.R")
source("funcoes_prontas/plotCycleQuality.R")
source("funcoes_prontas/plotNucleotideCount.R")
source("funcoes_prontas/readQualityScore.R")
source("funcoes_prontas/plotOcurrences.R")

# Carrega arquivos fastq
caminho_fastq <- "fastq_exemplos/output1.fastq"
caminho_fastq2 <- "fastq_exemplos/output2.fastq"
caminho_fastq3 <- "fastq_exemplos/adapters_exemplo1.fastq"
caminho_fastq4 <- "fastq_exemplos/adapters_exemplo2.fastq"

resultado_qa <- qa(caminho_fastq, type = "fastq")
resultado_qa2 <- qa(c(caminho_fastq, caminho_fastq2, caminho_fastq3), type = "fastq")

# Testa funcoes
plot_adapter <- plotAdapterContamination(c(caminho_fastq3))
p_estatico <- plot_adapter$p_estatico
p_interativo <- plot_adapter$p_interativo
tabela_adapter <- tableAdapterContamination(caminho_fastq5)
# --
tabela_frequencias <- freqSequences(resultado_qa)
tabela_frequencias <- freqSequences(resultado_qa2)
# ---
plots_qualidade_ciclo <- plotCycleQuality(resultado_qa2)
qualidade_ciclo_interativo <- plots_qualidade_ciclo$p_interativo
qualidade_ciclo_estatico <- plots_qualidade_ciclo$p_estatico
# ---
plot_contagens <- plotNucleotideCount(c(caminho_fastq, caminho_fastq2, caminho_fastq3))
plot_contagens_est <- plot_contagens$p_estatico
plot_contagens_int <- plot_contagens$p_interativo
# ---
plot_qualidade_media <- readQualityScore(resultado_qa2)
plot_qualidade_media_est <- plot_qualidade_media$p_estatico
plot_qualidade_media_int <- plot_qualidade_media$p_interativo
# ---
plot_ocorrencias <- plotOcurrences(resultado_qa2)
plot_ocorrencias_est <- plot_ocorrencias$p_estatico
plot_ocorrencias_int <- plot_ocorrencias$p_interativo

