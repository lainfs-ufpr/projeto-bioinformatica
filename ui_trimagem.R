# Função UI do Módulo de Trimagem
trimagemUI <- function(id) {
  ns <- NS(id) 
  
  sidebarLayout(
    sidebarPanel(
      class = "sidebar-trim",
      fileInput(ns("fastq_trim"), "Selecione arquivo FASTQ", accept = c(".fastq", ".fq")),
      textInput(ns("adapters"), "Arquivo de adaptadores",
                value = "Trimmomatic-0.39/adapters/TruSeq3-SE.fa"),
      numericInput(ns("sliding_window"), "Tamanho da janela (SLIDINGWINDOW)", 4, min = 1),
      numericInput(ns("quality_cutoff"), "Qualidade mínima", 20, min = 1),
      numericInput(ns("minlen"), "Comprimento mínimo (MINLEN)", 50, min = 1),
      textInput(ns("output_trim"), "Nome do arquivo de saída", "reads_trimmed.fastq"),
      actionButton(ns("run_trim"), "Rodar Trimmomatic", class = "btnQA-custom")
    ),
    mainPanel(
      div(class = "trim-main-panel",
          verbatimTextOutput(ns("trim_log")),
          div(style = "width: fit-content; margin-left: auto; margin-right: auto;",
              tableOutput(ns("trim_stats"))
          )
      )
    )
  )
}