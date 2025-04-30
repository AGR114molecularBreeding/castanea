# ===================================================================================================================================================
######################################################### SERVER OF THE APP #########################################################
# ===================================================================================================================================================

server <- function(input, output, session) {
  
  ########################################################### SERVER 1. RATIONALE ####################################################################  
  
  output$rationale_text <- renderUI({
    tags$p(
      "The increasing impact of climate change on forest ecosystems REQUIRES efficient tools to study key stress-response mechanisms. Among these, drought stress is a major threat to Fagaceae species, affecting their survival, distribution, and ecological roles. Despite the recognized importance of HSP90 proteins in mitigating stress responses, research on these genes in Fagaceae has been hindered by scattered genomic information across multiple databases.",
      tags$br(), tags$br(),
      "HSP90-Fagaceae was developed to address this gap by providing a centralized and user-friendly platform for accessing, visualizing, and analyzing HSP90-related data in major Fagaceae species. By integrating genomic resources in one place, this application facilitates comparative studies, evolutionary analyses, and research on drought tolerance, ultimately supporting conservation efforts and breeding programs."
    )
  })
  
  output$funding_text <- renderUI({
    tags$p(
      "This work is funded by the project ProyExcel_00351 (Proyectos de investigación de Excelencia, Junta de Andalucía)."
    )
  })
  
  
  ########################################################### SERVER 2. TAXONOMIC INFO ####################################################################  
  
  # Set default selection on selectinput
  observe({
    updateSelectInput(session, "var_species", selected = final_dataframe$SpecieII[1])
  })
  
  
  # Filter the information according to the selected species
  especie_seleccionada_info <- reactive({
    filtered_data <- final_dataframe[final_dataframe$SpecieII == input$var_species, ]
    return(filtered_data$Description[1])
  }) 
  
  # Show description of selected species
  output$info_box <- renderText({     
    descripcion <- especie_seleccionada_info()
    if (!is.null(descripcion)) {
      paste("<div style='padding: 10px; border: 1px solid #ccc; background-color: #f9f9f9; border-radius: 5px; word-wrap: break-word; margin: 0 20px;'>",
            "<div style='font-size: 17px;'>", descripcion, "</div>",
            "</div>")
    } else {
      "No information available for this species."
    }
  })  
  
  # Display CSS styles for text
  output$css <- renderUI({
    tags$style(HTML(".box-body {word-wrap: break-word;}"), scoped = TRUE)
  })
  
  # Store the index of the current image
  current_image_index <- reactiveVal(1)
  
  # Display the image of the selected species
  output$species_image <- renderImage({
    list(src = image_data[[input$var_species]]$paths[current_image_index()],
         alt = "Species Image",
         width = 500,
         height = 460)
  }, deleteFile = FALSE)
  
  # Change image on button click
  observeEvent(input$change_image_btn, {
    current_index <- current_image_index()
    current_index <- ifelse(current_index == length(image_data[[input$var_species]]$paths), 1, current_index + 1)
    current_image_index(current_index)
  })
  
  # Reset image index when species changes
  observeEvent(input$var_species, {
    current_image_index(1)
  })
  
  output$image_metadata <- renderUI({
    species <- input$var_species
    authors <- image_data[[species]]$authors[current_image_index()]
    source <- image_data[[species]]$source
    license <- image_data[[species]]$license
    
    # Adjust for cases where 'license' is a vector
    if (is.vector(license)) {
      license <- license[current_image_index()]
    }
    
    # Generate correctly formatted HTML
    HTML(
      paste0(
        "<span style='white-space: nowrap;'>",
        "Photograph by <strong>", authors, "</strong>, available on <strong>", source, "</strong>, ",
        "under a <a href='", license, "' target='_blank'>CC BY-NC 4.0 license</a>.",
        "</span>"
      )
    )
  })
  
  # Display map with Leaflet
  output$species_map <- renderLeaflet({
    req(input$var_species)
    species_shape <- shapefiles[[input$var_species]]
    leaflet(data = species_shape) %>%
      addTiles() %>%
      addPolygons(
        color = ~ifelse(id == 1, "green", "red"),
        weight = 1,
        opacity = 0.8,
        fillOpacity = 0.5,
        popup = ~paste("id:", id)
      ) %>%
      addLegend("bottomright", colors = c("green", "red"), labels = c("Natural", "Introduced"), title = "Area")
  })
  
  
  ##################################################### SERVER 3. HSP90 FAGACEAE DATABASE  ##########################################################
  
  # Reactive function that renames columns consistently
  final_dataframe_renamed <- reactive({
    df <- final_dataframe
    colnames(df) <- gsub("Conserved.domains", "Conserved domains", colnames(df))
    colnames(df) <- gsub("N..AA", "N. AA", colnames(df))
    return(df)
  })
  
  # Render the interactive table with renamed columns
  output$proteinTable <- renderDT({
    selected_columns <- input$show_vars_protein
    
    # If no columns are selected, display a message instead of the table
    if (length(selected_columns) == 0) {
      return(datatable(
        data.frame(Message = "No columns selected. Please select at least one column."),
        options = list(dom = 't'), # Hide table controls
        rownames = FALSE
      ))
    }
    
    # Create the unformatted table
    datatable_object <- datatable(
      final_dataframe_renamed()[, selected_columns, drop = FALSE],
      options = list(
        pageLength = 5,
        search = list(regex = FALSE, smart = FALSE, caseInsensitive = TRUE),
        columnDefs = list(
          list(className = "dt-center", targets = "_all")
        )
      ),
      class = "compact cell-border",
      extensions = 'Buttons',
      selection = "none"
    )
    
    # Apply formatting only if 'Specie' column is selected
    if ("Specie" %in% selected_columns) {
      datatable_object <- datatable_object %>%
        formatStyle('Specie', fontStyle = 'italic')
    }
    
    return(datatable_object)
  })
  
  # Handling table unloading with renamed columns
  output$downloadDATABASE <- downloadHandler(
    filename = function() {
      paste("protein_data", ".csv", sep = "")
    },
    content = function(file) {
      selected_columns <- input$show_vars_protein
      
      # Check if columns are selected
      if (length(selected_columns) == 0) {
        write.csv(data.frame(Message = "No columns selected."), file, row.names = FALSE)
        return()
      }
      
      # Check if the search box is empty
      if (input$proteinTable_search == "") {
        # If the search box is empty, download all rows
        data_to_download <- final_dataframe_renamed()
      } else {
        # If there is text in the search box, download only the displayed rows
        rows_to_download <- input$proteinTable_rows_all
        data_to_download <- final_dataframe_renamed()[rows_to_download, ]
      }
      
      # Download selected data
      write.csv(data_to_download[, selected_columns, drop = FALSE], file, row.names = FALSE)
    }
  )
  
  ########################################### SERVER 4.  HSP90 PROTEINS SEQUENCES ############################################################
  output$protein_selection <- renderUI({
    species <- input$Species
    proteins <- final_dataframe$Protein[final_dataframe$Specie == species]
    selectInput("Protein", 
                label = "Choose a protein", 
                choices = proteins)
  })
  
  output$output_sequence <- renderUI({
    selected_protein <- input$Protein
    selected_sequence <- final_dataframe$Sequence[final_dataframe$Protein == selected_protein]
    
    if (!is.null(selected_sequence)) {
      result <- paste(">", selected_protein, "\n", selected_sequence, sep = "")
      pre(
        style = "font-size: 17px; white-space: pre-wrap;",  # Add white-space for line break
        result
      )
    } else {
      HTML("No information found for the selected protein")
    }
  })
  
  # Function to generate the content of the FASTA file
  generate_fasta_content <- function() {
    selected_protein <- input$Protein
    selected_sequence <- final_dataframe$Sequence[final_dataframe$Protein == selected_protein]
    if (!is.null(selected_sequence)) {
      return(paste(">", selected_protein, "\n", selected_sequence, sep = ""))
    } else {
      return("No information found for the selected protein")
    }
  }
  
  #Function to generate the content of the FASTA file for all the sequences of the selected species.
  generate_all_fasta_content <- function() {
    selected_species <- input$Species
    species_data <- final_dataframe[final_dataframe$Specie == selected_species, ]
    fasta_content <- apply(species_data, 1, function(row) {
      paste(">", row["Protein"], "\n", row["Sequence"], sep = "")
    })
    return(paste(fasta_content, collapse = "\n"))
  }
  
  # Action when pressing the button to download the FASTA file of a sequence
  output$downloadamino <- downloadHandler(
    filename = function() {
      paste("sequence_", input$Protein, ".fasta", sep = "")
    },
    content = function(file) {
      fasta_content <- generate_fasta_content()
      writeLines(fasta_content, file)
    }
  )
  
  # Dynamically generate download button for all sequences of the selected species
  output$downloadall_button <- renderUI({
    species <- input$Species
    if (!is.null(species) && species != "") {
      downloadButton('downloadall', label = HTML(paste('Download all sequences of <i>', species, '</i>')), class = "btn btn-default")
    }
  })
  
  # Action when the button is pressed to download the FASTA file of all sequences of the selected species.
  output$downloadall <- downloadHandler(
    filename = function() {
      paste("all_sequences_", input$Species, ".fasta", sep = "")
    },
    content = function(file) {
      fasta_content <- generate_all_fasta_content()
      writeLines(fasta_content, file)
    }
  )  
  
  
  
  ##################################################### SERVER 5. PHYLOGENETIC MENU  ##########################################################
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$gg_image_output <- renderPlot({
    gg_image(filo, ranges)
  })
  
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  observeEvent(input$plot2_dblclick, {
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges_castanea$x <- c(brush$xmin, brush$xmax)
      ranges_castanea$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges_castanea$x <- NULL
      ranges_castanea$y <- NULL
    }
  })
  
  output$image_caption <- renderUI({
    tags$div(style = "text-align: center; width: 700px; margin: 0 auto;",
             tags$i("Ccr: Castanea crenata. Cd: Castanea dentata. Cmnh7: Castanea mollisima H7. Cmn11: Castanea mollisima N-11. Cmv: Castanea mollisima vanuxem. CSH1: Castanea sativa. Ch: Castanopisis hystrix. Qd: Quercus dentata. Q.ilex: Quercus ilex. Ql: Quercus lobata. Qr: Quercus robur. Qs: Quercus suber. Qv: Quercus variabilis. Ci: Carya illinoinensis.")
    )
  })
  
  output$inforr_box <- renderUI({
    tags$div(
      style = "text-align: left; padding: 10px; border: 1px solid #ccc; background-color: #f9f9f9; border-radius: 5px; word-wrap: break-word; margin: 20px 0;",
      tags$div(
        style = "font-size: 17px;",
        HTML("This phylogeny tree has been constructed using the Maximum Likelihood Tree method by generating 1,000 replicates with the Bootstrap method in Mega X. Additionally, the least reliable nodes were removed from the tree by setting a threshold for bootstrap support values at 50%. The tree visualization was further enhanced using <a href='https://itol.embl.de/' target='_blank'>iTOL</a> (<a href='https://doi.org/10.1093/nar/gkae268' target='_blank'>Letunic and Bork 2024</a>).")
      )
    )
  })
  
  
  # Download button UI
  output$download_button <- renderUI({
    downloadButton(outputId = "download_tree", label = "Download Tree.Newick", class = "btn-primary")
  })
  
  # Download handler
  output$download_tree <- downloadHandler(
    filename = function() {
      "Tree.Newick"
    },
    content = function(file) {
      file.copy("./Data/Tree.Newick", file)  # Ensure Tree.Newick is in the working directory
    }
  )
  
  ##################################################### SERVER 6. BLASTP  MENU  ##########################################################
  
  # Output for BLASTP results and protein table
  output$Tableprotein <- renderDT({
    req(input$run_blast)
    datatable(run_blast()$results, options = list(pageLength = 5), selection = 'single', editable = FALSE)
  })
  
  # Output for BLASTP usage instructions
  output$blast_instructions <- renderText({
    instructions <- "1. Select the species with which you want to compare your protein.
2. Enter your amino acid sequence with the description sequence using the \">\" symbol.
3. To run the Blastp, press the button \"Run Blast\".
4. To visualize the alignments, select the option of your interest in the results table.
5. To download the Blastp results in CSV, press the button \"Download data\". 
6. To download the Subject Sequence, press the button \"Download Subject Sequence\"."
    return(instructions)
  })
  
  # Reactive function to execute the BLAST
  run_blast <- eventReactive(input$run_blast, {
    req(input$fasta_file)  # Check that the file is uploaded
    
    # Reset table selection
    session$sendCustomMessage(type = "resetTableSelection", message = NULL)
    
    # Clearing the contents of the alignments
    output$alignments <- renderUI({
      HTML("<pre>No alignment selected yet.</pre>")
    })
    
    # Read sequences from loaded FASTA file
    fasta_file_path <- input$fasta_file$datapath
    fasta_sequences <- readLines(fasta_file_path, warn = FALSE)
    
    # Make sure there is a line break at the end of the file
    if (length(fasta_sequences) > 0 && !grepl("\\n$", fasta_sequences[length(fasta_sequences)])) {
      fasta_sequences <- c(fasta_sequences, "")  # Add an empty line if missing
    }
    
    if (length(fasta_sequences) == 0 || all(trimws(fasta_sequences) == "")) {
      output$sequence_error_message <- renderText("Error: The uploaded file is empty or invalid.")
      return(NULL)
    }
    
    # Extract sequences (ignore lines starting with '>')
    sequence_lines <- fasta_sequences[!grepl("^>", fasta_sequences)]
    
    # Concatenate the sequence lines to analyze their content
    combined_sequence <- paste(sequence_lines, collapse = "")
    
    # Check if it is a nucleotide sequence (A, T, G, C, U)
    if (grepl("^[ATGCU]+$", combined_sequence, ignore.case = TRUE)) {
      output$sequence_error_message <- renderText("Error: Nucleotide sequence detected. Please upload an amino acid sequence.")
      return(NULL)
    }
    
    # Check for invalid characters for amino acids
    if (!grepl("^[ARNDCQEGHILKMFPSTWYVX\\*]+$", combined_sequence, ignore.case = TRUE)) {
      output$sequence_error_message <- renderText("Error: Invalid amino acid sequence detected.")
      return(NULL)
    }
    
    # If everything is ok, we clear the error messages
    output$sequence_error_message <- renderText("")
    
    # Writes the sequences to a temporary file to use as a query
    query_fasta_file <- tempfile(fileext = ".fasta")
    writeLines(fasta_sequences, query_fasta_file)
    
    # Build the path to the reference file
    subject_fasta_file <- paste0(input$species_selection, ".fasta")
    subject_fasta_file <- normalizePath(subject_fasta_file)
    
    if (!file.exists(subject_fasta_file)) {
      showNotification("The subject file does not exist.", type = "error")
      return(NULL)
    }
    
    # Path to blastp executable
    blastp_executable <- file.path("./blastp")
    
    # Execute the blastp command in R
    blast_result <- system(paste(
      shQuote(blastp_executable),
      "-query", shQuote(query_fasta_file),
      "-subject", shQuote(subject_fasta_file),
      "-outfmt 5"  # XML output format
    ), intern = TRUE)
    
    # Parsing the BLAST results in XML format
    blast_xml <- read_xml(paste(blast_result, collapse = "\n"))
    
    # Function to extract information from the results
    extract_info <- function(node) {
      query_ID <- xml_text(xml_find_all(node, ".//Iteration_query-def"))
      subject_ID <- xml_text(xml_find_all(node, ".//Hit_def"))
      bitscore <- xml_text(xml_find_all(node, ".//Hsp_bit-score"))
      evalue <- xml_text(xml_find_all(node, ".//Hsp_evalue"))
      Hsp_identity <- xml_text(xml_find_all(node, ".//Hsp_identity"))
      Hsp_align <- xml_text(xml_find_all(node, ".//Hsp_align-len"))
      perc_id <- signif((as.numeric(Hsp_identity) / as.numeric(Hsp_align)) * 100, digits = 4)
      
      # Get the maximum length of the vectors
      max_length <- max(length(query_ID), length(subject_ID), length(bitscore), length(evalue), length(Hsp_identity), length(Hsp_align), na.rm = TRUE)
      
      # Fill the vectors with NA if necessary
      query_ID <- rep(query_ID, length.out = max_length)
      subject_ID <- rep(subject_ID, length.out = max_length)
      bitscore <- rep(bitscore, length.out = max_length)
      evalue <- rep(evalue, length.out = max_length)
      perc_id <- rep(perc_id, length.out = max_length)
      
      # Get the length of the query sequence
      query_length <- as.numeric(xml_text(xml_find_all(node, ".//Iteration_query-len")))
      
      # Get the values ​​of Query_start and Query_end
      query_start <- as.numeric(xml_text(xml_find_all(node, ".//Hsp_query-from")))
      query_end <- as.numeric(xml_text(xml_find_all(node, ".//Hsp_query-to")))
      
      # Calculate query coverage
      query_cover <- signif(((query_end - query_start + 1) / query_length) * 100, digits = 4)
      
      # Create a column with the name of the selected species
      species_selected <- rep(input$species_selection, length.out = max_length)
      
      data.frame(Query_ID = query_ID, 
                 Subject_ID = subject_ID,
                 Specie_Selected = species_selected,
                 Score = bitscore,
                 Query_Cover = query_cover,
                 E_Value = evalue,
                 Identity = perc_id)
    }
    
    # Extract information from each result
    results <- xml_find_all(blast_xml, "//Iteration")
    result_list <- lapply(results, extract_info)
    result_df <- do.call(rbind, result_list)
    
    # Get the list of all alignments
    alignments <- xml_find_all(blast_xml, "//Hit_hsps/Hsp")
    formatted_alignments <- sapply(alignments, function(x) format_alignment(x), USE.NAMES = FALSE)
    
    # Show notification that the search has ended
    shinyalert("Success!", "BLAST complete.", type = "success", animation = TRUE)
    
    list(results = result_df, alignments = formatted_alignments)
  })
  
  # Show the "Download Data" button only if the table has results
  output$downloadButtonUI <- renderUI({
    req(run_blast())  # Make sure the run_blast() function has been executed
    if (nrow(run_blast()$results) > 0) {
      downloadButton('downloadData', 'Download Data')
    } else {
      return(NULL)  # If there are no results, the button is not displayed
    }
  })
  
  # Display the "Download Subject Sequence" button only if a row is selected in the table
  output$downloadSubjectButtonUI <- renderUI({
    req(input$Tableprotein_rows_selected)  # Requires a row to be selected
    if (!is.null(input$Tableprotein_rows_selected)) {
      downloadButton('downloadSequence', 'Download Subject Sequence')
    } else {
      return(NULL)  # If no row is selected, the button is not displayed.
    }
  })
  
  # Download BLASTP results as CSV file
  output$downloadData <- downloadHandler(
    filename = function() {
      "blast_results.csv"
    },
    content = function(file) {
      blast_result <- run_blast()$results
      write.csv(blast_result, file, row.names = FALSE)
    }
  )
  
  # Download the sequence of the selected subject
  output$downloadSequence <- downloadHandler(
    filename = function() {
      selected_subject <- run_blast()$results$Subject_ID[input$Tableprotein_rows_selected]
      paste(selected_subject, ".fasta", sep = "")
    },
    content = function(file) {
      selected_subject <- run_blast()$results$Subject_ID[input$Tableprotein_rows_selected]
      fasta_file <- file.path(paste0(run_blast()$results$Specie_Selected[input$Tableprotein_rows_selected], ".fasta"))
      sequences <- readAAStringSet(fasta_file)
      selected_sequence <- toString(sequences[selected_subject])
      
      # Write the file making sure there are no extra line breaks
      writeLines(c(paste(">", gsub(" ", "", selected_subject), sep = ""), selected_sequence), file)
    }
  )
  
  # Observe Event to display the selected alignments in the table
  observeEvent(input$Tableprotein_rows_selected, {
    output$alignments <- renderUI({
      req(input$run_blast)
      
      # If no row is selected, display a message
      if (is.null(input$Tableprotein_rows_selected) || length(input$Tableprotein_rows_selected) == 0) {
        return(HTML("<pre>No alignment selected.</pre>"))
      }
      
      # Get alignment based on selected row
      selected_row <- input$Tableprotein_rows_selected
      alignments <- run_blast()$alignments
      
      # Check that the alignments are not empty
      if (length(alignments) >= selected_row) {
        alignment_to_show <- alignments[selected_row]
        
        # If the alignment has content, display it
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
  
  
  
  ##################################################### SERVER 6. ALIGN SEQUENCES  ##########################################################  
  
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


shinyApp(ui, server)

