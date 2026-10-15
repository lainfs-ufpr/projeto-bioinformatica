# Módulo de interface para gerenciamento de projetos do INTRA
# Requer: shiny, bslib
# Carregue com: source("projetos.R")

projetosUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("tela"))
  )
}

projetosServer <- function(id, diretorio_projetos = "dados/projetos") {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

    if (!dir.exists(diretorio_projetos)) {
      dir.create(diretorio_projetos, recursive = TRUE, showWarnings = FALSE)
    }

    tela_atual <- shiny::reactiveVal("lista")
    projeto_atual <- shiny::reactiveVal(NULL)
    atualizar_lista <- shiny::reactiveVal(0)
    erro <- shiny::reactiveVal(NULL)
    projeto_para_excluir <- shiny::reactiveVal(NULL)

    caminho_projeto <- function(id_projeto) {
      file.path(diretorio_projetos, id_projeto)
    }

    ler_metadados <- function(pasta) {
      arq <- file.path(pasta, "projeto.rds")
      if (!file.exists(arq)) return(NULL)
      tryCatch(readRDS(arq), error = function(e) NULL)
    }

    salvar_metadados <- function(projeto) {
      projeto$atualizado_em <- Sys.time()
      saveRDS(projeto, file.path(caminho_projeto(projeto$id), "projeto.rds"))
      projeto
    }

    listar_projetos <- function() {
      pastas <- list.dirs(diretorio_projetos, recursive = FALSE, full.names = TRUE)
      itens <- lapply(pastas, ler_metadados)
      itens <- Filter(Negate(is.null), itens)
      if (length(itens) == 0) return(list())
      itens[order(vapply(itens, function(x) as.numeric(x$atualizado_em), numeric(1)),
                  decreasing = TRUE)]
    }

    contar_arquivos <- function(id_projeto) {
      pasta <- file.path(caminho_projeto(id_projeto), "arquivos_originais")
      if (!dir.exists(pasta)) return(0L)
      length(list.files(pasta, recursive = TRUE, all.files = FALSE))
    }

    arquivos_projeto <- function(id_projeto) {
      pasta <- file.path(caminho_projeto(id_projeto), "arquivos_originais")
      if (!dir.exists(pasta)) return(data.frame(
        Arquivo = character(), Tamanho = character(), Modificado = character(),
        stringsAsFactors = FALSE
      ))
      caminhos <- list.files(pasta, recursive = TRUE, full.names = TRUE)
      caminhos <- caminhos[!dir.exists(caminhos)]
      if (!length(caminhos)) return(data.frame(
        Arquivo = character(), Tamanho = character(), Modificado = character(),
        stringsAsFactors = FALSE
      ))
      info <- file.info(caminhos)
      data.frame(
        Arquivo = substring(caminhos, nchar(pasta) + 2),
        Tamanho = vapply(info$size, function(x) {
          if (is.na(x)) return("—")
          if (x < 1024^2) paste0(round(x / 1024, 1), " KB")
          else paste0(round(x / 1024^2, 2), " MB")
        }, character(1)),
        Modificado = format(info$mtime, "%d/%m/%Y %H:%M"),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }

    formatar_data <- function(x) {
      if (is.null(x) || length(x) == 0 || is.na(x)) return("—")
      format(x, "%d/%m/%Y %H:%M")
    }

    output$tela <- shiny::renderUI({
      atualizar_lista()
      projeto <- projeto_atual()

      if (identical(tela_atual(), "novo")) {
        return(shiny::tagList(
          shiny::div(class = "container-fluid py-4",
            shiny::actionButton(ns("voltar_lista_novo"), "← Meus projetos",
                                class = "btn btn-outline-secondary mb-3"),
            shiny::h2("Criar novo projeto"),
            shiny::p("Organize os arquivos e os resultados de uma análise em um único espaço."),
            shiny::wellPanel(
              shiny::textInput(ns("nome_projeto"), "Nome do projeto *",
                               placeholder = "Ex.: Amostras do experimento 01"),
              shiny::textAreaInput(ns("descricao_projeto"), "Descrição (opcional)",
                                   placeholder = "Objetivo, amostras ou observações...",
                                   rows = 3, width = "100%"),
              shiny::fileInput(ns("arquivos_iniciais"), "Adicionar arquivos (opcional)",
                               multiple = TRUE,
                               accept = c(".fastq", ".fq", ".fasta", ".fa", ".fastq.gz", ".fq.gz")),
              shiny::p(class = "text-muted small",
                       "Você também pode criar um projeto vazio e adicionar arquivos depois."),
              shiny::div(class = "d-flex gap-2",
                shiny::actionButton(ns("cancelar_novo"), "Cancelar",
                                    class = "btn btn-outline-secondary"),
                shiny::actionButton(ns("criar_projeto"), "Criar projeto",
                                    icon = shiny::icon("folder-plus"),
                                    class = "btn btn-primary")
              )
            )
          )
        ))
      }

      if (identical(tela_atual(), "workspace") && !is.null(projeto)) {
        arquivos <- arquivos_projeto(projeto$id)
        return(shiny::tagList(
          shiny::div(class = "container-fluid py-4",
            shiny::actionButton(ns("voltar_lista_workspace"), "← Meus projetos",
                                class = "btn btn-outline-secondary mb-3"),
            shiny::div(class = "d-flex justify-content-between align-items-start flex-wrap gap-3",
              shiny::div(
                shiny::h2(projeto$nome),
                shiny::p(if (nzchar(projeto$descricao)) projeto$descricao else
                           "Este projeto ainda não tem descrição.",
                         class = "text-muted"),
              ),
              shiny::actionButton(ns("abrir_edicao"), "Editar descrição",
                                  icon = shiny::icon("pen"),
                                  class = "btn btn-outline-primary")
            ),
            bslib::navset_card_tab(
              nav_panel("Arquivos",
                shiny::h4("Arquivos do projeto"),
                shiny::p("Os arquivos enviados são copiados para a pasta deste projeto."),
                shiny::fileInput(ns("arquivos_workspace"), "Adicionar arquivos",
                                 multiple = TRUE,
                                 accept = c(".fastq", ".fq", ".fasta", ".fa",
                                            ".fastq.gz", ".fq.gz")),
                shiny::actionButton(ns("adicionar_arquivos"), "Adicionar ao projeto",
                                    icon = shiny::icon("upload"),
                                    class = "btn btn-primary"),
                shiny::br(), shiny::br(),
                if (nrow(arquivos) == 0) {
                  shiny::div(class = "alert alert-light",
                             "Nenhum arquivo foi adicionado a este projeto.")
                } else {
                  shiny::tableOutput(ns("tabela_arquivos"))
                }
              ),
              nav_panel("Qualidade",
                        qualidadeUI("qualidade_id") 
              ),
              nav_panel("Trimagem",
                        trimagemUI("trimagem_id")
              ),
              nav_panel("Histórico",
                shiny::div(class = "p-3",
                  shiny::h4("Histórico do projeto"),
                  shiny::p("Criado em: ", formatar_data(projeto$criado_em)),
                  shiny::p("Última atualização: ", formatar_data(projeto$atualizado_em)),
                  shiny::p(class = "text-muted",
                           "O histórico detalhado de análises será implementado junto à integração dos módulos.")
                )
              ),
              shiny::div(
                shiny::p(class = "small text-muted",
                         paste("Criado em", formatar_data(projeto$criado_em),
                               "· Atualizado em", formatar_data(projeto$atualizado_em)))
              )
            )
          )
        ))
      }

      projetos <- listar_projetos()
      filtro <- input$busca_projetos %||% ""
      if (nzchar(trimws(filtro)) && length(projetos)) {
        projetos <- Filter(function(p) {
          grepl(filtro, p$nome, ignore.case = TRUE) ||
            grepl(filtro, p$descricao, ignore.case = TRUE)
        }, projetos)
      }

      shiny::tagList(
        shiny::div(class = "container-fluid py-4",
          shiny::div(class = "d-flex justify-content-between align-items-center flex-wrap gap-3 mb-4",
            shiny::div(
              shiny::h2("Meus projetos"),
              shiny::p("Gerencie seus arquivos e acompanhe as análises em um só lugar.",
                       class = "text-muted")
            ),
            shiny::actionButton(ns("novo_projeto"), "Novo projeto",
                                icon = shiny::icon("plus"),
                                class = "btn btn-primary btn-lg")
          ),
          shiny::textInput(ns("busca_projetos"), NULL,
                           placeholder = "Buscar projeto por nome ou descrição",
                           width = "100%"),
          shiny::uiOutput(ns("mensagem_erro")),
          shiny::uiOutput(ns("cards_projetos"))
        )
      )
    })

    # Renderização separada dos cards: permite atualizar a busca sem depender
    # de reconstruir toda a tela.
    output$cards_projetos <- shiny::renderUI({
      atualizar_lista()
      shiny::req(identical(tela_atual(), "lista"))
      projetos <- listar_projetos()
      filtro <- input$busca_projetos %||% ""
      if (nzchar(trimws(filtro))) {
        projetos <- Filter(function(p) {
          grepl(filtro, p$nome, ignore.case = TRUE) ||
            grepl(filtro, p$descricao, ignore.case = TRUE)
        }, projetos)
      }

      if (!length(projetos)) {
        return(shiny::div(class = "text-center p-5 border rounded",
          shiny::icon("folder-open", class = "fa-3x text-muted"),
          shiny::h4("Nenhum projeto encontrado"),
          shiny::p("Crie seu primeiro projeto para começar a organizar suas análises.")
        ))
      }

      shiny::tagList(lapply(projetos, function(p) {
        shiny::div(class = "card mb-3 shadow-sm",
          shiny::div(class = "card-body",
            shiny::div(class = "d-flex justify-content-between align-items-start flex-wrap gap-3",
              shiny::div(
                shiny::h4(p$nome, class = "card-title"),
                shiny::p(if (nzchar(p$descricao)) p$descricao else "Sem descrição.",
                         class = "card-text text-muted"),
                shiny::p(class = "small text-muted",
                  paste0("Atualizado em ", formatar_data(p$atualizado_em),
                         " · ", contar_arquivos(p$id), " arquivo(s)"))
              ),
              shiny::div(class = "d-flex gap-2",
                shiny::tags$button(
                  type = "button", class = "btn btn-primary",
                  onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'});",
                                    ns("selecionar_projeto"), p$id),
                  shiny::icon("folder-open"), " Abrir"
                ),
                shiny::tags$button(
                  type = "button", class = "btn btn-outline-danger",
                  onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'});",
                                    ns("solicitar_exclusao"), p$id),
                  shiny::icon("trash"), " Excluir"
                )
              )
            )
          )
        )
      }))
    })

    output$mensagem_erro <- shiny::renderUI({
      if (is.null(erro())) return(NULL)
      shiny::div(class = "alert alert-danger", erro())
    })

    output$tabela_arquivos <- shiny::renderTable({
      shiny::req(identical(tela_atual(), "workspace"), !is.null(projeto_atual()))
      arquivos_projeto(projeto_atual()$id)
    }, striped = TRUE, bordered = FALSE, hover = TRUE, spacing = "s")

    shiny::observeEvent(input$novo_projeto, {
      erro(NULL)
      tela_atual("novo")
    })

    shiny::observeEvent(input$voltar_lista_novo, tela_atual("lista"))
    shiny::observeEvent(input$cancelar_novo, tela_atual("lista"))
    shiny::observeEvent(input$voltar_lista_workspace, {
      projeto_atual(NULL)
      tela_atual("lista")
      atualizar_lista(isolate(atualizar_lista()) + 1)
    })

    shiny::observeEvent(input$criar_projeto, {
      nome <- trimws(input$nome_projeto %||% "")
      descricao <- trimws(input$descricao_projeto %||% "")
      if (!nzchar(nome)) {
        erro("Informe um nome para o projeto.")
        return()
      }

      id <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_",
                   paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = ""))
      pasta <- caminho_projeto(id)
      dir.create(file.path(pasta, "arquivos_originais"), recursive = TRUE,
                 showWarnings = FALSE)
      dir.create(file.path(pasta, "qualidade"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(pasta, "trimagem"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(pasta, "historico"), recursive = TRUE, showWarnings = FALSE)

      projeto <- list(
        id = id,
        nome = nome,
        descricao = descricao,
        criado_em = Sys.time(),
        atualizado_em = Sys.time(),
        status = "Em andamento"
      )

      uploads <- input$arquivos_iniciais
      if (!is.null(uploads) && nrow(uploads) > 0) {
        destino <- file.path(pasta, "arquivos_originais")
        ok <- file.copy(uploads$datapath, file.path(destino, basename(uploads$name)),
                        overwrite = FALSE)
        if (any(!ok)) {
          erro("O projeto foi criado, mas alguns arquivos não puderam ser copiados.")
        }
      }

      saveRDS(projeto, file.path(pasta, "projeto.rds"))
      projeto_atual(projeto)
      tela_atual("workspace")
      atualizar_lista(isolate(atualizar_lista()) + 1)
      shiny::showNotification("Projeto criado com sucesso.", type = "message")
    })

    shiny::observeEvent(input$selecionar_projeto, {
      id_projeto <- input$selecionar_projeto
      projeto <- ler_metadados(caminho_projeto(id_projeto))
      if (is.null(projeto)) {
        erro("Não foi possível abrir os dados deste projeto.")
        return()
      }
      projeto_atual(projeto)
      erro(NULL)
      tela_atual("workspace")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$solicitar_exclusao, {
      id_projeto <- input$solicitar_exclusao
      projeto <- ler_metadados(caminho_projeto(id_projeto))
      if (is.null(projeto)) {
        shiny::showNotification("Não foi possível localizar este projeto.",
                                type = "error")
        return()
      }
      projeto_para_excluir(id_projeto)
      shiny::showModal(shiny::modalDialog(
        title = "Excluir projeto?",
        shiny::p("Esta ação excluirá o projeto e os arquivos armazenados nele."),
        shiny::p(shiny::strong(projeto$nome)),
        footer = shiny::tagList(
          shiny::modalButton("Cancelar"),
          shiny::actionButton(ns("confirmar_exclusao"),
                              "Excluir definitivamente",
                              class = "btn btn-danger")
        ),
        easyClose = TRUE
      ))
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$confirmar_exclusao, {
      id_projeto <- projeto_para_excluir()
      shiny::req(id_projeto)
      unlink(caminho_projeto(id_projeto), recursive = TRUE, force = TRUE)
      projeto_para_excluir(NULL)
      shiny::removeModal()
      atualizar_lista(isolate(atualizar_lista()) + 1)
      shiny::showNotification("Projeto excluído.", type = "message")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$adicionar_arquivos, {
      projeto <- projeto_atual()
      shiny::req(projeto)
      uploads <- input$arquivos_workspace
      if (is.null(uploads) || nrow(uploads) == 0) {
        shiny::showNotification("Selecione pelo menos um arquivo.", type = "warning")
        return()
      }
      destino <- file.path(caminho_projeto(projeto$id), "arquivos_originais")
      dir.create(destino, recursive = TRUE, showWarnings = FALSE)
      nomes <- basename(uploads$name)
      resultados <- file.copy(uploads$datapath, file.path(destino, nomes),
                              overwrite = FALSE)
      if (any(!resultados)) {
        shiny::showNotification(
          "Alguns arquivos não foram copiados. Verifique se já existem arquivos com o mesmo nome.",
          type = "warning", duration = 8
        )
      } else {
        shiny::showNotification("Arquivos adicionados ao projeto.", type = "message")
      }
      projeto <- salvar_metadados(projeto)
      projeto_atual(projeto)
      atualizar_lista(isolate(atualizar_lista()) + 1)
    })

    shiny::observeEvent(input$abrir_edicao, {
      projeto <- projeto_atual()
      shiny::req(projeto)
      shiny::showModal(shiny::modalDialog(
        title = "Editar projeto",
        shiny::textInput(ns("editar_nome"), "Nome do projeto", value = projeto$nome),
        shiny::textAreaInput(ns("editar_descricao"), "Descrição",
                             value = projeto$descricao, rows = 3, width = "100%"),
        footer = shiny::tagList(
          shiny::modalButton("Cancelar"),
          shiny::actionButton(ns("salvar_edicao"), "Salvar alterações",
                              class = "btn btn-primary")
        ),
        easyClose = TRUE
      ))
    })

    shiny::observeEvent(input$salvar_edicao, {
      projeto <- projeto_atual()
      shiny::req(projeto)
      nome <- trimws(input$editar_nome %||% "")
      if (!nzchar(nome)) {
        shiny::showNotification("O nome do projeto não pode ficar vazio.", type = "error")
        return()
      }
      projeto$nome <- nome
      projeto$descricao <- trimws(input$editar_descricao %||% "")
      projeto <- salvar_metadados(projeto)
      projeto_atual(projeto)
      shiny::removeModal()
      atualizar_lista(isolate(atualizar_lista()) + 1)
      shiny::showNotification("Projeto atualizado.", type = "message")
    })
  })
}
