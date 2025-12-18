trimagemUI <- function(id) {
  ns <- NS(id) 
  
  sidebarLayout(
    sidebarPanel(
      class = "sidebar-trim",
      
      radioButtons(ns("mode"), "Tipo de leitura:",
                   choices = c("Single-end" = "SE", "Paired-end" = "PE"),
                   selected = "SE"),
      
      conditionalPanel(
        condition = "input.mode == 'SE'",
        ns = ns,
        div(class = "botao-arquivo",
            fileInput(ns("r1"), 
                      "FASTQ R1 (gzip aceito .fastq(.gz))", 
                      accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz")
            )
        )
      ),
      
      conditionalPanel(
        condition = "input.mode == 'PE'",
        ns = ns,
        div(class = "lado-lado",
            div(class = "botao-arquivo",
                fileInput(ns("r1"), 
                          "FASTQ R1 (gzip aceito .fastq(.gz))", 
                          accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz")
                )
            ),
            div(class = "botao-arquivo",
                fileInput(ns("r2"), 
                          "FASTQ R2 (gzip aceito .fastq(.gz))", 
                          accept = c(".fastq", ".fq", ".fastq.gz", ".fq.gz")
                ) 
            )
        )
      ),
      
      textInput(ns("trimmomatic_path"), "Caminho do Trimmomatic .jar", 
                value = "Trimmomatic-0.39/trimmomatic-0.39.jar"),
      
      selectInput(ns("adapter_choice"), "Escolher adaptador:",
                  choices = c("TruSeq3-PE" = "TruSeq3-PE.fa",
                              "TruSeq3-SE" = "TruSeq3-SE.fa",
                              "TruSeq2-PE" = "TruSeq2-PE.fa",
                              "TruSeq2-SE" = "TruSeq2-SE.fa",
                              "NexteraPE"  = "NexteraPE-PE.fa",
                              "Custom (inserir manualmente)" = "custom"),
                  selected = "TruSeq3-PE.fa"),
      
      conditionalPanel(
        condition = "input.adapter_choice == 'custom'", 
        ns = ns, 
        textInput(ns("custom_adapter"), "Caminho completo do arquivo de adaptadores:",
                  value = "")
      ),
      
      div(class = "coluna",
          div(class = "lado-lado",
              numericInput(ns("window_size"), "SLIDINGWINDOW - tamanho da janela", value = 4, min = 1),
              numericInput(ns("qual_cut"), "SLIDINGWINDOW - cutoff de qualidade", value = 20, min = 2),
              ),
          numericInput(ns("minlen"), "MINLEN - comprimento mínimo (bp)", value = 25, min = 1),
          div(class = "lado-lado",
              numericInput(ns("leading"), "LEADING - corta na esquerda (Q)", value = 3, min = 0),
              numericInput(ns("trailing"), "TRAILING - corta na direita (Q)", value = 3, min = 0),
              ),
      ),
      
      actionButton(ns("run"), "Rodar Trimmomatic", class = "btnQA-custom"), 
      hr(),
      
      conditionalPanel(
        condition = "input.mode == 'PE'",
        ns = ns,
        downloadButton(ns("download_r1_paired"), 
                        "Download R1_paired",
                        class = "btn-trim-download"),
        downloadButton(ns("download_r2_paired"), 
                       "Download R2_paired",
                       class = "btn-trim-download") 
      ),
      
      conditionalPanel(
        condition = "input.mode == 'SE'",
        ns = ns,
        downloadButton(ns("download_single_trimmed"), 
                        "Download SE_trimmed", 
                        class = "btn-trim-download")
        )
    ),
    
    mainPanel(
      div(class = "trim-main-panel",
          verbatimTextOutput(ns("log")),
          h4("Estatísticas (antes x depois)"),
          tableOutput(ns("stats_table")), 
          h4("Plot: qualidade média por ciclo (R1)"),
          plotOutput(ns("qual_plot"), height = "360px") 
      )
    )
  )
}