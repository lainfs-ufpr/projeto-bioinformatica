library(shiny)
library(shinyFiles)
library(ShortRead)
library(ggplot2) 
library(dplyr)  
library(tidyr)   
library(scales)  
library(viridis) 
library(RColorBrewer) 

options(shiny.maxRequestSize = 100 * 1024^4)

ui <- fluidPage(
  titlePanel("Seleção de Pasta, Upload e Controle de Qualidade"),
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("diretorio", "Escolher Pasta", "Selecionar"),
      verbatimTextOutput("caminhoPasta"),
      fileInput("arquivos", "Escolha os arquivos",
                multiple = TRUE,
                accept = c(".fasta", ".fa", ".fastq", ".fq", "text/plain")),
      actionButton("run_analysis", "Rodar Controle de Qualidade")
    ),
    mainPanel(
      tabsetPanel( 
        tabPanel("Conteúdo Arquivos",
                 tableOutput("conteudoArquivos")),
        tabPanel("QA Output",
                 verbatimTextOutput("qa_output")),
        tabPanel("Contagem de Bases por Ciclo (Log)",
                 plotOutput("plot_count_log")),
        tabPanel("Heatmap de Frequência de Bases", 
                 plotOutput("plot_heatmap")),
        tabPanel("Contagem de Bases (Barras)", 
                 plotOutput("plot_bar_count")),
        tabPanel("Distribuição de Contagens por Base", 
                 plotOutput("plot_violin_dist")),
        tabPanel("Proporção de Bases por Ciclo (Área)",
                 plotOutput("plot_area_prop")),
        tabPanel("Count por Ciclo - Base e Lane (Log)", 
                 plotOutput("plot_facet_base_lane")),
        tabPanel("Proporção de N's por Ciclo e Lane",
                 plotOutput("plot_proportion_N"))
      )
    )
  )
)

server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "C:" = "C:/", "D:" = "D:/")
  shinyDirChoose(input, "diretorio", roots = volumes, session = session)
  
  dir_path <- reactiveVal(NULL)
  qa_obj_reactive <- reactiveVal(NULL) 
  df_log_reactive <- reactiveVal(NULL) 
  df_prop_reactive <- reactiveVal(NULL) 
  df_proportion_N_reactive <- reactiveVal(NULL) 
  
  observeEvent(input$diretorio, {
    path <- parseDirPath(volumes, input$diretorio)
    dir_path(path)
  })
  
  output$caminhoPasta <- renderText({
    req(dir_path())
    paste("Pasta selecionada:", dir_path())
  })
  
  output$conteudoArquivos <- renderTable({
    req(input$arquivos)
    data.frame(Nome_Arquivo = input$arquivos$name, Caminho = input$arquivos$datapath)
  })
  
  observeEvent(input$run_analysis, {
    showNotification("Iniciando análise de controle de qualidade...", type = "message") 
    
    
    fls <- NULL
    if (!is.null(input$arquivos)) {
      fls <- input$arquivos$datapath
    } else if (!is.null(dir_path())) {
      
      fls <- dir_path()
    } else {
      showNotification("Nenhum arquivo ou pasta fornecido para análise!", type = "error")
      return()
    }
    
    if (is.null(fls) || length(fls) == 0) {
      showNotification("Nenhum arquivo FASTQ válido encontrado para análise!", type = "error")
      return()
    }
    
    tryCatch({
      #executa o controle de qualidade
      qa_result <- qa(fls, type = "fastq")
      qa_obj_reactive(qa_result) 
      
      qa_cycle <- qa_result[["perCycle"]]
      df_cycle <- qa_cycle$baseCall
      
      df_log_temp <- df_cycle %>%
        filter(Base %in% c("A", "C", "G", "T", "N")) %>% 
        mutate(
          Base = factor(Base, levels = c("A", "C", "G", "T", "N")),
          logCount = log10(Count + 1) #evita log10(0)
        )
      df_log_reactive(df_log_temp) 
      
  
      df_prop_temp <- df_cycle %>%
        group_by(Cycle, lane) %>%
        mutate(Proportion = Count / sum(Count)) %>%
        ungroup()
      df_prop_reactive(df_prop_temp) 
      
  
      df_proportion_N_temp <- df_log_temp %>% 
        group_by(Cycle, lane) %>%
        mutate(Total_Count_Per_Cycle_Lane = sum(Count)) %>%
        ungroup() %>%
        filter(Base == "N") %>%
        mutate(Proportion_N = Count / Total_Count_Per_Cycle_Lane) %>%
        select(Cycle, lane, Proportion_N)
      df_proportion_N_reactive(df_proportion_N_temp)
      

      showNotification("Análise de controle de qualidade concluída e dados de plot preparados!", type = "message") 
      

      
    }, error = function(e) {
      showNotification(paste("Erro na análise:", e$message), type = "error", duration = 10)
      qa_obj_reactive(NULL)
      df_log_reactive(NULL)
      df_prop_reactive(NULL)
      df_proportion_N_reactive(NULL)
    })
  })
  

  
  output$plot_count_log <- renderPlot({
    req(df_log_reactive())
    df_plot <- df_log_reactive() %>% filter(Base %in% c("A", "C", "G", "T")) 
    ggplot(df_plot, aes(x = Cycle, y = logCount, color = Base)) +
      geom_line(stat = "identity", alpha = 0.7) +
      facet_wrap(~ lane) +
      scale_color_manual(values = c("A" = "#1B9E77", "C" = "#D95F02",
                                    "G" = "#7570B3", "T" = "#E7298A")) +
      labs(
        title = "Contagem de Bases por Ciclo (Log)",
        x = "Ciclo",
        y = "log10(Count)",
        color = "Base"
      ) +
      theme_minimal() +
      theme(strip.text = element_text(face = "bold"))
  })
  
  output$plot_heatmap <- renderPlot({
    req(df_log_reactive())
    df_plot <- df_log_reactive() %>% filter(Base %in% c("A", "C", "G", "T"))
    ggplot(df_plot, aes(x = Cycle, y = Base, fill = logCount)) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(title = "Heatmap: Frequência log10 de Bases por Ciclo") +
      theme_minimal()
  })
  
  output$plot_bar_count <- renderPlot({
    req(df_log_reactive())
    df_plot <- df_log_reactive() %>% filter(Base %in% c("A", "C", "G", "T"))
    ggplot(df_plot, aes(x = Cycle, y = Count, fill = Base)) +
      geom_bar(stat = "identity", position = "dodge") +
      labs(title = "Contagem de Bases por Ciclo") +
      theme_minimal()
  })
  
  output$plot_violin_dist <- renderPlot({
    req(df_log_reactive())
    df_plot <- df_log_reactive() %>% filter(Base %in% c("A", "C", "G", "T"))
    ggplot(df_plot, aes(x = Base, y = logCount, fill = Base)) +
      geom_violin() +
      labs(title = "Distribuição das Contagens por Base") +
      theme_minimal()
  })
  
  output$plot_area_prop <- renderPlot({
    req(df_prop_reactive())
    ggplot(df_prop_reactive(), aes(x = Cycle, y = Proportion, fill = Base)) +
      geom_area() +
      labs(
        title = "Proporção de Bases por Ciclo (Área Empilhada)",
        x = "Ciclo",
        y = "Proporção"
      ) +
      theme_minimal() +
      scale_fill_brewer(palette = "Set1")
  })
  
  output$plot_facet_base_lane <- renderPlot({
    req(df_log_reactive())
    df_plot <- df_log_reactive() %>% filter(Base %in% c("A", "C", "G", "T"))
    ggplot(df_plot, aes(x = Cycle, y = logCount)) +
      geom_line(aes(color = Base)) +
      facet_grid(Base ~ lane) +
      labs(
        title = "Count por Ciclo - Base e Lane",
        x = "Ciclo",
        y = "log10(Count)"
      ) +
      theme_minimal()
  })
  
  output$plot_proportion_N <- renderPlot({
    req(df_proportion_N_reactive()) 
    ggplot(df_proportion_N_reactive(), aes(x = Cycle, y = Proportion_N)) +
      geom_line(aes(color = lane)) +
      facet_wrap(~ lane, scales = "free_y") +
      labs(
        title = "Proporção de Bases 'N' por Ciclo e Lane",
        x = "Ciclo",
        y = "Proporção de 'N's",
        color = "Lane"
      ) +
      theme_minimal() +
      scale_y_continuous(labels = scales::percent)
  })
  
}

shinyApp(ui = ui, server = server)