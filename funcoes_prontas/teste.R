source("funcoes_prontas/adapterContamination.R")
source("funcoes_prontas/freqSequences.R")
source("funcoes_prontas/plotCycleQuality.R")
source("funcoes_prontas/plotNucleotideCount.R")
source("funcoes_prontas/readQualityScore.R")

caminho_fastq <- "fastq_gerados/output2.fastq"
caminho_fastq <- "teste_funcoes/adapter_contamination/mock2.fastq"

resultado_qa <- qa(caminho_fastq, type = "fastq")

adapter_plot <- adapterContamination(caminho_fastq)

frequencias_tabela <- freqSequences(resultado_qa)

plots_qualidade_ciclo <- plotCycleQuality(resultado_qa)
qualidade_ciclo_interativo <- plots_qualidade_ciclo$p_interativo
qualidade_ciclo_estatico <- plots_qualidade_ciclo$p_estatico

plot_contagens <- plotNucleotideCount(resultado_qa)

plot_qualidade_media <- readQualityScore(resultado_qa)
