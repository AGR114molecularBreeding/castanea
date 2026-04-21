# ===============================================================
# SERVER
# ===============================================================
server <- function(input, output, session) {
  
  ## To download the proteins from the icon
  fasta_modal <- reactiveVal(NULL)
  protein_modal <- reactiveVal(NULL)
  
  dir.create("data/sessions", recursive = TRUE, showWarnings = FALSE)
  dir.create("data/hidden_sessions", recursive = TRUE, showWarnings = FALSE)
  
  stored_data <- reactiveVal(NULL)
  
  refresh_sessions <- function() {
    gsub("\\.rds$", "",
         list.files("data/sessions", pattern = "\\.rds$"))
  }
  
  observe({
    updateSelectInput(session, "load_session",
                      choices = refresh_sessions())
  })
  
  # -------------------------------------------------------------
  # Read CSV or stored data
  # -------------------------------------------------------------
  reactive_dataframe <- reactive({
    
    req(input$data_source)
    
    # -------------------------
    # USE STORED SESSION
    # -------------------------
    if (input$data_source == "session") {
      if (is.null(stored_data())) return(NULL)
      return(stored_data())
    }
    
    # -------------------------
    # USE UPLOADED CSV
    # -------------------------
    if (is.null(input$upload_csv)) return(NULL)
    
    df <- read_delim(
      file = input$upload_csv$datapath,
      delim = NULL,
      locale = locale(encoding = "UTF-8"),
      col_types = cols(.default = "c"),
      trim_ws = TRUE
    )
    df <- as.data.frame(df, check.names = FALSE)
    
    validate(
      need(all(c("Protein", "Sequence") %in% colnames(df)),
           "CSV must contain Protein and Sequence columns")
    )
    
    df
  })
  
  
  
  # -------------------------------------------------------------
  # Load session
  # -------------------------------------------------------------
  observeEvent(input$load_session, {
    req(input$load_session)
    stored_data(readRDS(
      file.path("data/sessions",
                paste0(input$load_session, ".rds"))
    ))
  })
  
  # -------------------------------------------------------------
  # Save session
  # -------------------------------------------------------------
  observeEvent(input$save_session, {
    req(reactive_dataframe(), input$session_name)
    
    saveRDS(
      reactive_dataframe(),
      file.path("data/sessions",
                paste0(input$session_name, ".rds"))
    )
    
    updateSelectInput(session, "load_session",
                      choices = refresh_sessions())
    
    showNotification("Session saved", type = "message")
  })
  
  # -------------------------------------------------------------
  # Merge session
  # -------------------------------------------------------------
  observeEvent(input$merge_session, {
    req(input$upload_csv, input$load_session)
    
    # -------------------------
    # Read uploaded CSV
    # -------------------------
    csv_df <- read.csv(
      input$upload_csv$datapath,
      sep = ";",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      encoding = "UTF-8"
    )
    
    # -------------------------
    # Read selected session
    # -------------------------
    session_df <- readRDS(
      file.path("data/sessions",
                paste0(input$load_session, ".rds"))
    )
    
    # -------------------------
    # VALIDATION: same structure
    # -------------------------
    if (!identical(colnames(csv_df), colnames(session_df))) {
      showAlert(
        title = "Merge error",
        text = "CSV and selected session must have the same columns (same names and order).",
        type = "error"
      )
      return(NULL)
    }
    
    # -------------------------
    # Merge CSV + session
    # -------------------------
    merged <- rbind(session_df, csv_df)
    
    # -------------------------
    # Save OVER the selected session
    # -------------------------
    saveRDS(
      merged,
      file.path("data/sessions",
                paste0(input$load_session, ".rds"))
    )
    
    # -------------------------
    # Update app state
    # -------------------------
    stored_data(merged)
    updateSelectInput(session, "data_source", selected = "session")
    
    showNotification("CSV successfully merged with selected session",
                     type = "message")
  })
  
  
  # -------------------------------------------------------------
  # Delete session
  # -------------------------------------------------------------
  observeEvent(input$delete_session, {
    req(input$load_session)
    
    file.remove(
      file.path("data/sessions",
                paste0(input$load_session, ".rds"))
    )
    
    stored_data(NULL)
    updateSelectInput(session, "load_session",
                      choices = refresh_sessions())
  })
  
  # -------------------------------------------------------------
  # Update checkbox columns
  # -------------------------------------------------------------
  observe({
    df <- reactive_dataframe()
    req(df)
    
    updateCheckboxGroupInput(
      session,
      "show_vars_protein",
      choices = colnames(df),
      selected = colnames(df)
    )
  })
  
  # -------------------------------------------------------------
  # Prepare data
  # -------------------------------------------------------------
  final_dataframe_renamed <- reactive({
    df <- reactive_dataframe()
    req(df)
    
    df$Sequence_real <- df$Sequence
    
    df$Sequence <- paste0(
      '<i class="fas fa-dna" style="cursor:pointer;color:#28a745;font-size:20px;" ',
      'data-row="', seq_len(nrow(df)), '"></i>'
    )
    
    df
  })
  
  
  # -------------------------------------------------------------
  # Render table
  # -------------------------------------------------------------
  output$proteinTable <- renderDT({
    df_all <- final_dataframe_renamed()
    req(input$show_vars_protein)
    
    valid_cols <- intersect(input$show_vars_protein, colnames(df_all))
    df <- df_all[, valid_cols, drop = FALSE]
    
    format_go <- function(x) {
      sapply(x, function(v) {
        if (is.na(v) || v == "") return("")
        paste0("• ", trimws(unlist(strsplit(v, ";"))),
               collapse = "<br>")
      })
    }
    
    if ("GO Terms" %in% colnames(df)) df$`GO Terms` <- format_go(df$`GO Terms`)
    if ("GO IDs" %in% colnames(df)) df$`GO IDs` <- format_go(df$`GO IDs`)
    
    dt <- datatable(
      df,
      options = list(
        pageLength = 5,
        search = list(regex = FALSE),
        columnDefs = list(
          list(className = "dt-center", targets = "_all")
        )
      ),
      escape = FALSE,
      class = "compact cell-border",
      selection = "none"
    )
    
    if ("Specie" %in% colnames(df))
      dt <- dt %>% formatStyle("Specie", fontStyle = "italic")
    if ("Genus" %in% colnames(df))
      dt <- dt %>% formatStyle("Genus", fontStyle = "italic")
    if ("Gene" %in% colnames(df))
      dt <- dt %>% formatStyle("Gene", fontStyle = "italic")
    
    dt
  })
  
  # -------------------------------------------------------------
  # Download
  # -------------------------------------------------------------
  output$downloadDATABASE <- downloadHandler(
    filename = function() "protein_data_filtered.csv",
    
    content = function(file) {
      
      selected_columns <- input$show_vars_protein
      req(selected_columns)
      
      # Real Dataframe
      df <- reactive_dataframe()
      
      # Exact rows that the user sees (after filtering/searching)
      rows <- input$proteinTable_rows_all
      
      if (!is.null(rows) && length(rows) > 0) {
        df <- df[rows, , drop = FALSE]
      }
      
      # Clean columns with HTML bullets
      cols_to_clean <- c("Conserved domains", "GO IDs", "GO Terms")
      for (col in intersect(cols_to_clean, colnames(df))) {
        df[[col]] <- sapply(df[[col]], function(x) {
          if (is.na(x) || x == "") return("")
          terms <- unlist(strsplit(x, "[•â€¢]", perl = TRUE))
          terms <- trimws(terms[terms != ""])
          paste(terms, collapse = ", ")
        })
      }
      
      write.csv(
        df[, selected_columns, drop = FALSE],
        file,
        row.names = FALSE,
        fileEncoding = "UTF-8"
      )
    }
  )
  
  # -------------------------------------------------------------
  # FASTA modal
  # -------------------------------------------------------------
  observeEvent(input$selected_protein_row, {
    
    df <- final_dataframe_renamed()
    i <- input$selected_protein_row
    
    protein <- df$Protein[i]
    seq <- df$Sequence_real[i]
    
    fasta <- paste0(">", protein, "\n", seq)
    
    fasta_modal(fasta)
    protein_modal(protein)
    
    showModal(
      modalDialog(
        title = paste("Protein sequence:", protein),
        pre(
          style = "white-space: pre-wrap; word-break: break-all;",
          fasta
        ),
        footer = tagList(
          modalButton("Close"),
          downloadButton("download_modal_fasta", "Download FASTA")
        ),
        easyClose = TRUE
      )
    )
    output$download_modal_fasta <- downloadHandler(
      filename = function() paste0(protein, ".fasta"),
      content = function(file) writeLines(fasta, file)
    )
    
  })
  
  ############################################################################################################################################################
  # Protein Sequences module (INDEPENDENT)
  #####################################################################
  
  ############################################################
  # MODULE 2 — COMPLETELY INDEPENDENT STATE
  ############################################################
  
  stored_data_seq <- reactiveVal(NULL)
  
  load_seq_session <- function(session_name) {
    stored_data_seq(
      readRDS(
        file.path("data/sessions", paste0(session_name, ".rds"))
      )
    )
  }
  # -------------------------------------------------------------
  # Sessions list
  # -------------------------------------------------------------
  refresh_sessions_seq <- function() {
    gsub("\\.rds$", "", list.files("data/sessions", pattern = "\\.rds$"))
  }
  
  observe({
    updateSelectInput(
      session,
      "load_session_seq",
      choices = refresh_sessions_seq()
    )
  })
  
  # -------------------------------------------------------------
  # Load selected session (MODULE 2 ONLY)
  # -------------------------------------------------------------
  observeEvent(input$load_session_seq, {
    req(input$data_source_seq == "session", input$load_session_seq)
    load_seq_session(input$load_session_seq)
  })
  
  observeEvent(input$load_session_seq, {
    
    # Reset selections and domains
    updateCheckboxGroupInput(
      session,
      "seq_protein_ids",
      choices = character(0),
      selected = character(0)
    )
    
    shinyjs::reset("upload_domains")
    
  }, ignoreInit = TRUE)
  
  
  observeEvent(input$data_source_seq, {
    
    if (input$data_source_seq == "session") {
      req(input$load_session_seq)
      load_seq_session(input$load_session_seq)
    }
    
    if (input$data_source_seq == "current") {
      stored_data_seq(reactive_dataframe())
    }
    
    # -------------------------
    # RESET STATE (VERY IMPORTANT)
    # -------------------------
    updateCheckboxGroupInput(
      session,
      "seq_protein_ids",
      choices = character(0),
      selected = character(0)
    )
    
    shinyjs::reset("upload_domains")
    
  }, ignoreInit = TRUE)
  
  # -------------------------------------------------------------
  # Use currently loaded dataset (from app global state)
  # -------------------------------------------------------------
  observe({
    if (input$data_source_seq == "current") {
      req(exists("reactive_dataframe"))
      stored_data_seq(reactive_dataframe())
    }
  })
  
  # -------------------------------------------------------------
  # Load & clean CD-search file (Batch format supported)
  # -------------------------------------------------------------
  domains_df <- reactive({
    req(input$upload_domains)
    
    lines <- readLines(input$upload_domains$datapath)
    lines <- lines[!grepl("^#", lines)]
    lines <- lines[nzchar(lines)]
    
    header_line <- grep("^Query", lines)
    req(length(header_line) == 1)
    
    header <- strsplit(lines[header_line], "\t")[[1]]
    data_lines <- lines[(header_line + 1):length(lines)]
    
    df <- do.call(rbind, lapply(data_lines, function(x) strsplit(x, "\t")[[1]]))
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    colnames(df) <- header
    
    df$Protein <- trimws(sub("^.*>", "", df$Query))
    df$From <- as.numeric(df$From)
    df$To   <- as.numeric(df$To)
    
    df
  })
  
  # -------------------------------------------------------------
  # Filter domains by type
  # -------------------------------------------------------------
  filtered_domains_df <- reactive({
    df <- domains_df()
    req(df, input$domain_types)
    df[df$`Hit type` %in% input$domain_types, , drop = FALSE]
  })
  
  # -------------------------------------------------------------
  # Proteins AVAILABLE = intersection dataset × domains
  # -------------------------------------------------------------
  available_proteins <- reactive({
    req(stored_data_seq())
    
    # No domains uploaded → no proteins
    if (is.null(input$upload_domains)) {
      return(character(0))
    }
    
    dom <- filtered_domains_df()
    req(nrow(dom) > 0)
    
    sort(
      intersect(
        stored_data_seq()$Protein,
        unique(dom$Protein)
      )
    )
  })
  
  # -------------------------------------------------------------
  # Update protein selector (MULTI)
  # -------------------------------------------------------------
  observe({
    prots <- available_proteins()
    
    updateCheckboxGroupInput(
      session,
      "seq_protein_ids",
      choices = prots,
      selected = character(0)
    )
  })
  observeEvent(input$select_all_proteins, {
    updateCheckboxGroupInput(
      session,
      "seq_protein_ids",
      selected = available_proteins()
    )
  })
  
  observeEvent(input$deselect_all_proteins, {
    updateCheckboxGroupInput(
      session,
      "seq_protein_ids",
      selected = character(0)
    )
  })
  # -------------------------------------------------------------
  # Selected proteins
  # -------------------------------------------------------------
  selected_proteins_df <- reactive({
    df <- stored_data_seq()
    req(df, input$seq_protein_ids)
    df[df$Protein %in% input$seq_protein_ids, , drop = FALSE]
  })
  
  # -------------------------------------------------------------
  # Selected domains
  # -------------------------------------------------------------
  selected_domains_df <- reactive({
    df <- filtered_domains_df()
    req(df, input$seq_protein_ids)
    df[df$Protein %in% input$seq_protein_ids, , drop = FALSE]
  })
  
  # -------------------------------------------------------------
  # Protein summary
  # -------------------------------------------------------------
  output$protein_fasta_view <- renderUI({
    df <- selected_proteins_df()
    req(nrow(df) > 0)
    
    tagList(
      lapply(seq_len(nrow(df)), function(i) {
        
        header <- paste0(">", df$Protein[i])
        
        
        tags$pre(
          style = "
          white-space: pre-wrap;
          word-break: break-all;
          font-family: monospace;
          background-color: #f8f9fa;
          padding: 12px;
          border-radius: 6px;
          margin-bottom: 20px;
          border: 1px solid #dee2e6;
        ",
          paste(header, df$Sequence[i], sep = "\n")
        )
      })
    )
  })
  
  
  # -------------------------------------------------------------
  # Domains plot
  # -------------------------------------------------------------
  domains_plot <- reactive({
    dom <- selected_domains_df()
    prot <- selected_proteins_df()
    req(nrow(dom) > 0, nrow(prot) > 0)
    
    prot_len <- data.frame(
      Protein = prot$Protein,
      start = 0,
      end = as.numeric(prot$`N. AA`)
    )
    
    # 1. CALCULATION OF THE EXACT LIMIT (+100 AA)
    max_total <- max(prot_len$end)
    extension_eje <- max_total + 100 
    
    ggplot() +
      # Protein
      geom_rect(
        data = prot_len,
        aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
        fill = "grey90"
      ) +
      # Protein domains
      geom_rect(
        data = dom,
        aes(xmin = From, xmax = To, ymin = 0, ymax = 1, fill = `Short name`),
        color = "black"
      ) +
      # Final numbering of each sequence
      geom_text(
        data = prot_len,
        aes(x = end, y = 0.5, label = end),
        hjust = -0.3, 
        size = 5,
        fontface = "bold"
      ) +
      facet_grid(Protein ~ .) + 
      scale_x_continuous(
        limits = c(0, extension_eje),
        breaks = seq(0, extension_eje, by = 100),
        expand = c(0, 0)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        # Sequences names
        strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", size = 14),
        
        # SCALE SEPARATION
        axis.text.x = element_text(size = 13, face = "bold", margin = margin(t = 30)), 
        axis.title.x = element_text(size = 15, face = "bold", margin = margin(t = 30)),
        
        # Leyend
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        
        # Spacing and aesthetics
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.8),
        axis.ticks.x = element_line(color = "black"),
        plot.margin = margin(10, 50, 10, 10)
      ) +
      labs(
        title = "Conserved Domain Architecture",
        x = "Distance (Amino Acids)", 
        fill = "Conserved Domains")
  }) 
  output$protein_domains_plot <- renderPlot({
    domains_plot()
  })
  
  output$download_domains_plot <- downloadHandler(
    filename = function() {
      suffix <- if (isTRUE(input$card_dominios_full_screen)) "_expanded" else ""
      paste0("conserved_domains_", Sys.Date(), suffix, ".png")
    },
    content = function(file) {
      # 1. Detect if it is expanded
      is_full <- isTRUE(input$card_dominios_full_screen)
      
      # 2. We define dynamic dimensions
      # If it's expanded, we increase its width and height so it doesn't look pixelated.
      ancho <- if (is_full) 16 else 10
      alto  <- if (is_full) 10 else 6
      
      ggsave(
        filename = file,
        plot = domains_plot(), # We use the reagent that we already have
        width = ancho,
        height = alto,
        units = "in",
        dpi = 300,
        bg = "white"
      )
    }
  )
  
  
  
  # -------------------------------------------------------------
  # FASTA download — MODULE 2 (ONLY ONE)
  # -------------------------------------------------------------
  output$download_all_visible_fasta <- downloadHandler(
    filename = function() {
      paste0("selected_proteins_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      df <- selected_proteins_df()
      req(nrow(df) > 0)
      
      fasta <- unlist(
        lapply(seq_len(nrow(df)), function(i) {
          header <- paste0(">", df$Protein[i])
          paste(header, df$Sequence[i], sep = "\n")
        })
      )
      
      writeLines(fasta, file)
    }
  )
  
  # ===============================================================
  # DYNAMIC CONTROL OF DOWNLOAD BUTTONS
  # ===============================================================
  observe({
    # --- BUTTON 1: Database Module (downloadDATABASE) ---
    # It is enabled only if the current dataframe has rows
    df_db <- tryCatch(reactive_dataframe(), error = function(e) NULL)
    shinyjs::toggleState("downloadDATABASE", condition = !is.null(df_db) && nrow(df_db) > 0)
    
    # --- BUTTON 2: Sequences Module - FASTA (download_all_visible_fasta) ---
    # It is only enabled if protein IDs are selected in the checkboxes
    shinyjs::toggleState("download_all_visible_fasta", 
                         condition = !is.null(input$seq_protein_ids) && length(input$seq_protein_ids) > 0)
    
    # --- BUTTON 3: Sequences Module - Graph (download_domains_plot) ---
    # It is enabled only if domains are loaded and proteins are selected.
    # We use tryCatch to avoid errors if the reactant does not yet exist
    doms <- tryCatch(selected_domains_df(), error = function(e) NULL)
    shinyjs::toggleState("download_domains_plot", condition = !is.null(doms) && nrow(doms) > 0)
  })
  ######################################################################### phylogenia
  
  # Create session folder
  dir.create("data/Tree_session", recursive = TRUE, showWarnings = FALSE)
  
  # Tree ID activo
  itol_tree_manual  <- reactiveVal(NULL)
  itol_tree_session <- reactiveVal(NULL)
  itol_tree_id      <- reactiveVal(NULL)
  
  observeEvent(input$tree_source, {
    
    if (input$tree_source == "manual") {
      itol_tree_id(NULL)
      updateTextInput(session, "itol_id_input", value = "")
    }
    
  })
  
  # List saved sessions
  refresh_tree_sessions <- function() {
    gsub(
      "\\.rds$", "",
      list.files("data/Tree_session", pattern = "\\.rds$")
    )
  }
  
  # Update session selector
  observe({
    updateSelectInput(
      session,
      "load_tree_session",
      choices = refresh_tree_sessions()
    )
  })
  
  # Load manual tree
  observeEvent(input$load_itol_tree, {
    req(input$tree_source == "manual", input$itol_id_input)
    
    itol_tree_manual(trimws(input$itol_id_input))
    itol_tree_id(itol_tree_manual())
  })
  
  # Save session
  observeEvent(input$save_tree_session, {
    req(itol_tree_id(), input$tree_session_name)
    
    saveRDS(
      itol_tree_id(),
      file.path(
        "data/Tree_session",
        paste0(input$tree_session_name, ".rds")
      )
    )
    
    updateSelectInput(
      session,
      "load_tree_session",
      choices = refresh_tree_sessions(),
      selected = input$tree_session_name
    )
    
    showNotification("Tree session saved", type = "message")
  })
  
  # Load session
  observe({
    req(input$tree_source == "session", input$load_tree_session)
    
    itol_tree_session(
      readRDS(
        file.path(
          "data/Tree_session",
          paste0(input$load_tree_session, ".rds")
        )
      )
    )
    
    itol_tree_id(itol_tree_session())
  })
  
  # Switch between manual / session trees WITHOUT losing state
  # -------------------------------------------------------------
  observeEvent(input$tree_source, {
    
    if (input$tree_source == "manual") {
      itol_tree_id(itol_tree_manual())
    }
    
    if (input$tree_source == "session") {
      itol_tree_id(itol_tree_session())
    }
    
  })
  
  # Delete session
  observeEvent(input$delete_tree_session, {
    req(input$load_tree_session)
    
    file.remove(
      file.path(
        "data/Tree_session",
        paste0(input$load_tree_session, ".rds")
      )
    )
    
    itol_tree_id(NULL)
    
    updateSelectInput(
      session,
      "load_tree_session",
      choices = refresh_tree_sessions()
    )
    
    showNotification("Tree session deleted", type = "warning")
  })
  
  # Show tree
  output$itol_frame <- renderUI({
    req(itol_tree_id())
    
    tags$div(
      style = "height: calc(100vh - 120px); padding: 10px;",
      tags$iframe(
        src = paste0("https://itol.embl.de/tree/", itol_tree_id()),
        width = "100%",
        height = "100%",
        frameborder = "0"
      )
    )
  })
  
  ##################################################### SERVER 7. BLASTP  MENU  ##########################################################
  
  # ================================
  # Load proteome files dynamically
  # ================================
  observe({
    
    proteome_files <- list.files(
      path = "./proteomes",
      pattern = "\\.fasta$",
      full.names = FALSE
    )
    
    display_names <- gsub("\\.fasta$", "", proteome_files)
    choices <- setNames(proteome_files, display_names)
    
    updateSelectInput(
      session,
      "species_selection",
      choices = choices
    )
    
  })
  
  
  output$Tableprotein <- renderDT({
    req(input$run_blast)
    datatable(run_blast()$results, options = list(pageLength = 5), selection = 'single', editable = FALSE)
  })
  
  output$blast_instructions <- renderText({
    instructions <- "1. Select the algorithm (BLASTp or DIAMOND).
2. Select the target species to compare against your query protein.
3. Enter your amino acid sequence in FASTA format (starting with the \">\" symbol).
4. Click \"Run Analysis\" to start the search.
5. To visualize alignments, select a specific entry from the results table.
6. Click \"Download data\" to save the results as a CSV file. 
7. Click  \"Download Subject Sequence\" to export the matching sequence."
    return(instructions)
  })
  
  run_blast <- eventReactive(input$run_blast, {
    req(input$fasta_file)  
    
    session$sendCustomMessage(type = "resetTableSelection", message = NULL)
    
    output$alignments <- renderUI({
      HTML("<pre>No alignment selected yet.</pre>")
    })
    
    fasta_file_path <- input$fasta_file$datapath
    fasta_sequences <- readLines(fasta_file_path, warn = FALSE)
    
    if (length(fasta_sequences) > 0 && !grepl("\\n$", fasta_sequences[length(fasta_sequences)])) {
      fasta_sequences <- c(fasta_sequences, "")  
    }
    
    if (length(fasta_sequences) == 0 || all(trimws(fasta_sequences) == "")) {
      output$sequence_error_message <- renderText("Error: The uploaded file is empty or invalid.")
      return(NULL)
    }
    
    sequence_lines <- fasta_sequences[!grepl("^>", fasta_sequences)]
    combined_sequence <- paste(sequence_lines, collapse = "")
    
    if (grepl("^[ATGCU]+$", combined_sequence, ignore.case = TRUE)) {
      output$sequence_error_message <- renderText("Error: Nucleotide sequence detected. Please upload an amino acid sequence.")
      return(NULL)
    }
    
    if (!grepl("^[ARNDCQEGHILKMFPSTWYVX\\*]+$", combined_sequence, ignore.case = TRUE)) {
      output$sequence_error_message <- renderText("Error: Invalid amino acid sequence detected.")
      return(NULL)
    }
    
    output$sequence_error_message <- renderText("")
    
    query_fasta_file <- tempfile(fileext = ".fasta")
    writeLines(fasta_sequences, query_fasta_file)
    
    subject_fasta_file <- file.path("./proteomes", input$species_selection)
    subject_fasta_file <- normalizePath(subject_fasta_file)
    
    if (!file.exists(subject_fasta_file)) {
      showNotification("The subject file does not exist.", type = "error")
      return(NULL)
    }
    
    ### Conditional logic for Blastp vs Diamond
    if (input$program_selection == "BLASTp") {
      
      blastp_executable <- file.path("./data/bin2/Blastp/blastp.exe")
      
      blast_result <- system(paste(
        shQuote(blastp_executable),
        "-query", shQuote(query_fasta_file),
        "-subject", shQuote(subject_fasta_file),
        "-outfmt 5" 
      ), intern = TRUE)
      
      blast_xml <- read_xml(paste(blast_result, collapse = "\n"))
      
    } else if (input$program_selection == "DIAMOND") {
      
      diamond_executable <- file.path("./data/bin2/Diamond/diamond.exe")
      
      # 1. Create a temporary Diamond database from FASTA
      temp_db <- tempfile()
      system(paste(
        shQuote(diamond_executable), "makedb",
        "--in", shQuote(subject_fasta_file),
        "-d", shQuote(temp_db)
      ), ignore.stdout = TRUE, ignore.stderr = TRUE)
      
      # 2.Run the Diamond analysis
      temp_xml_out <- tempfile(fileext = ".xml")
      system(paste(
        shQuote(diamond_executable), "blastp",
        "-q", shQuote(query_fasta_file),
        "-d", shQuote(temp_db),
        "--outfmt", "5",       #Diamond supports Outfm 5 (XML format, same as BLASTp)
        "-o", shQuote(temp_xml_out)
      ), ignore.stdout = TRUE, ignore.stderr = TRUE)
      
      # Read the XML file generated by Diamond
      blast_xml <- read_xml(temp_xml_out)
    }
    
    # Parsing the BLAST/Diamond results in XML format
    extract_info <- function(node) {
      query_ID <- xml_text(xml_find_all(node, ".//Iteration_query-def"))
      subject_ID <- xml_text(xml_find_all(node, ".//Hit_def"))
      bitscore <- xml_text(xml_find_all(node, ".//Hsp_bit-score"))
      evalue <- xml_text(xml_find_all(node, ".//Hsp_evalue"))
      Hsp_identity <- xml_text(xml_find_all(node, ".//Hsp_identity"))
      Hsp_align <- xml_text(xml_find_all(node, ".//Hsp_align-len"))
      perc_id <- signif((as.numeric(Hsp_identity) / as.numeric(Hsp_align)) * 100, digits = 4)
      
      max_length <- max(length(query_ID), length(subject_ID), length(bitscore), length(evalue), length(Hsp_identity), length(Hsp_align), na.rm = TRUE)
      
      query_ID <- rep(query_ID, length.out = max_length)
      subject_ID <- rep(subject_ID, length.out = max_length)
      bitscore <- rep(bitscore, length.out = max_length)
      evalue <- rep(evalue, length.out = max_length)
      perc_id <- rep(perc_id, length.out = max_length)
      
      query_length <- as.numeric(xml_text(xml_find_all(node, ".//Iteration_query-len")))
      query_start <- as.numeric(xml_text(xml_find_all(node, ".//Hsp_query-from")))
      query_end <- as.numeric(xml_text(xml_find_all(node, ".//Hsp_query-to")))
      
      query_cover <- signif(((query_end - query_start + 1) / query_length) * 100, digits = 4)
      species_selected <- rep(input$species_selection, length.out = max_length)
      
      data.frame(Query_ID = query_ID, 
                 Subject_ID = subject_ID,
                 Specie_Selected = species_selected,
                 Score = bitscore,
                 Query_Cover = query_cover,
                 E_Value = evalue,
                 Identity = perc_id)
    }
    
    results <- xml_find_all(blast_xml, "//Iteration")
    result_list <- lapply(results, extract_info)
    result_df <- do.call(rbind, result_list)
    
    alignments <- xml_find_all(blast_xml, "//Hit_hsps/Hsp")
    formatted_alignments <- sapply(alignments, function(x) format_alignment(x), USE.NAMES = FALSE)
    
    shinyalert("Success!", "Analysis complete.", type = "success", animation = TRUE)
    
    list(results = result_df, alignments = formatted_alignments)
  })
  
  output$downloadButtonUI <- renderUI({
    req(run_blast())  
    if (nrow(run_blast()$results) > 0) {
      downloadButton('downloadData', 'Download Data')
    } else {
      return(NULL)  
    }
  })
  
  output$downloadSubjectButtonUI <- renderUI({
    req(input$Tableprotein_rows_selected)  
    if (!is.null(input$Tableprotein_rows_selected)) {
      downloadButton('downloadSequence', 'Download Subject Sequence')
    } else {
      return(NULL)  
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "results.csv"
    },
    content = function(file) {
      blast_result <- run_blast()$results
      write.csv(blast_result, file, row.names = FALSE)
    }
  )
  
  output$downloadSequence <- downloadHandler(
    filename = function() {
      selected_subject <- run_blast()$results$Subject_ID[input$Tableprotein_rows_selected]
      paste(selected_subject, ".fasta", sep = "")
    },
    content = function(file) {
      selected_subject <- run_blast()$results$Subject_ID[input$Tableprotein_rows_selected]
      fasta_file <- file.path("./proteomes", run_blast()$results$Specie_Selected[input$Tableprotein_rows_selected])
      sequences <- readAAStringSet(fasta_file)
      selected_sequence <- toString(sequences[selected_subject])
      
      writeLines(c(paste(">", gsub(" ", "", selected_subject), sep = ""), selected_sequence), file)
    }
  )
  
  observeEvent(input$Tableprotein_rows_selected, {
    output$alignments <- renderUI({
      req(input$run_blast)
      
      if (is.null(input$Tableprotein_rows_selected) || length(input$Tableprotein_rows_selected) == 0) {
        return(HTML("<pre>No alignment selected.</pre>"))
      }
      
      selected_row <- input$Tableprotein_rows_selected
      alignments <- run_blast()$alignments
      
      if (length(alignments) >= selected_row) {
        alignment_to_show <- alignments[selected_row]
        
        if (nchar(alignment_to_show) > 0) {
          HTML(paste("<pre>", alignment_to_show, "</pre>"))
        } else {
          HTML("<pre>No alignment found for the selected row.</pre>")
        }
      } else {
        HTML("<pre>No alignment data available.</pre>")
      }
    })
  })
  
  ##################################################### SERVER 8. ALIGN SEQUENCES  ##########################################################  
  
  # On app load, hide spinner and alignment container
  observe({
    hide("loading_spinner")
    hide("alignment_output_container")
  })
  
  # Reactive variable to store the alignment
  aligned_seqs_reactive <- reactiveVal(NULL)
  
  observeEvent(input$align_btn, {
    req(input$align_fasta_file)  # Make sure the file has been uploaded
    
    # Clear the message before a new alignment
    output$message <- renderText("")
    
    aligned_seqs_reactive(NULL) # Clear the previous result
    
    # Show spinner and hide previous result
    show("loading_spinner")
    hide("alignment_output_container")
    
    # Read fasta file according to the selected stream type
    align_fasta_file <- input$align_fasta_file$datapath
    seq_data <- NULL
    detected_type <- NULL
    warning_detected <- FALSE  # Variable to control specific warnings
    
    tryCatch({
      options(warn = 2)
      
      if (input$seq_type == "DNA") {
        seq_data <- readDNAStringSet(align_fasta_file)
        detected_type <- "DNA"
        colorscheme <- "nucleotide" 
      } else if (input$seq_type == "AA") {
        seq_data <- readAAStringSet(align_fasta_file)
        detected_type <- "Aminoácidos"
        colorscheme <- "clustal"
      }
      
      options(warn = 0)
    }, warning = function(w) {
      warning_message <- conditionMessage(w)
      if (grepl("invalid one-letter sequence codes", warning_message)) {
        output$message <- renderText("Warning: The sequence contains invalid characters. Please check your file.")
        warning_detected <<- TRUE
      } else {
        output$message <- renderText(paste("Advertencia:", warning_message))
      }
    }, error = function(e) {
      output$message <- renderText("Error reading the Fasta file. Please check the format and try again.")
    }, finally = {
      options(warn = 0)
    })
    
    # If a specific warning was detected, stop execution
    if (warning_detected) {
      hide("loading_spinner")
      return(NULL)
    }
    
    # Validate that there are at least two sequences
    if (is.null(seq_data) || length(seq_data) < 2) {
      output$message <- renderText("The file must contain at least two sequences to perform the alignment.")
      hide("loading_spinner")
      return(NULL)
    }
    
    # Perform the alignment
    alignment <- msa(seq_data, method = input$align_method)
    aligned_seqs <- as(alignment, "XStringSet")
    aligned_seqs_reactive(aligned_seqs)
    
    # Display success message with detected sequence type
    output$message <- renderText(paste("Alignment completed successfully."))
    
    # Show alignment in the interface
    output$alignment_output_container <- renderUI({
      msaR(
        aligned_seqs_reactive(),
        colorscheme = colorscheme,
        labelNameLength = 150,
        rowheight = 15,
        alignmentHeight = length(aligned_seqs_reactive()) * 15
      )
    })
    
    # Hide the spinner and show the result
    hide("loading_spinner")
    show("alignment_output_container")
  })
}

# ===============================================================
# Run app
# ===============================================================
shinyApp(ui, server)