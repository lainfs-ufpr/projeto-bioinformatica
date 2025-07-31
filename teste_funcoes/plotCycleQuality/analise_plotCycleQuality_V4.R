library(ShortRead)
library(dplyr)
library(plotly)

# Carrega os dados de qualidade por ciclo
qa_result <- qa("fastq_gerados/", pattern = ".fastq", type = "fastq")
df <- qa_result[["perCycle"]][["quality"]]

# Calcula a mediana ponderada por ciclo e amostra
df_median <- df %>%
  group_by(lane, Cycle) %>%
  summarise(
    Median = quantile(rep(Score, Count), 0.5),
    .groups = "drop"
  )

# Cria o plotly básico (linhas por amostra)
p <- plot_ly()

# Adiciona uma linha para cada amostra
for (sample_name in unique(df_median$lane)) {
  df_sample <- df_median %>% filter(lane == sample_name)
  p <- add_trace(p,
                 data = df_sample,
                 x = ~Cycle,
                 y = ~Median,
                 type = 'scatter',
                 mode = 'lines',
                 name = sample_name,
                 hoverinfo = 'text',
                 text = ~paste("Sample:", sample_name,
                               "<br>Cycle:", Cycle,
                               "<br>Median Score:", Median))
}

# Adiciona os retângulos de fundo (estilo FastQC)
p <- layout(p,
            shapes = list(
              list(type = "rect", x0 = 0, x1 = max(df_median$Cycle), y0 = 0, y1 = 20,
                   fillcolor = "rgba(255,0,0,0.1)", line = list(width = 0), layer = "below"),
              list(type = "rect", x0 = 0, x1 = max(df_median$Cycle), y0 = 20, y1 = 28,
                   fillcolor = "rgba(255,255,0,0.1)", line = list(width = 0), layer = "below"),
              list(type = "rect", x0 = 0, x1 = max(df_median$Cycle), y0 = 28, y1 = 40,
                   fillcolor = "rgba(0,255,0,0.1)", line = list(width = 0), layer = "below")
            ),
            title = "Per base sequence quality (FastQC-style)",
            xaxis = list(title = "Cycle"),
            yaxis = list(title = "Quality Score", range = c(0, 42))
)

# Exibe
p

# Salvar
library(htmlwidgets)

saveWidget(p, "teste_funcoes/plot_cycle_quality/grafico_interativo_all_samples.html", selfcontained = TRUE)
