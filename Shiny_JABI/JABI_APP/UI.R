# ===============================================================
# Load required packages 
# ===============================================================
library(Biostrings)
library(DT)
library(grid)
library(ggplot2)
library(leaflet)
library(msa)
library(msaR)
library(openxlsx)
library(png)
library(sf)
library(shiny)
library(shinyalert)
library(shinycssloaders)
library(shinyjs)
library(xml2)
library(bslib)
library(readr)
library(shinyWidgets)

# ===============================================================
# App options
# ===============================================================
options(shiny.host = "0.0.0.0")
options(shiny.port = 80)
options(shiny.maxRequestSize = 2*1024^2)


# =============================================================================================
# MENU BLASTP
# =============================================================================================

###############################  Formatting the alignment  ############################################
format_alignment <- function(alignment) {
  query_seq <- xml_text(xml_find_all(alignment, "./Hsp_qseq")[1])
  midline <- xml_text(xml_find_all(alignment, "./Hsp_midline")[1])
  subject_seq <- xml_text(xml_find_all(alignment, "./Hsp_hseq")[1])
  
  max_width <- 60
  
  split_sequence <- function(sequence) {
    n <- nchar(sequence)
    num_lines <- ceiling(n / max_width)
    split_seq <- sapply(1:num_lines, function(i) {
      start <- (i - 1) * max_width + 1
      end <- min(i * max_width, n)
      substr(sequence, start, end)
    })
    split_seq
  }
  
  query_lines <- split_sequence(query_seq)
  midline_lines <- split_sequence(midline)
  subject_lines <- split_sequence(subject_seq)
  
  query_length <- nchar(query_seq)
  
  formatted_text <- "<pre>"
  for (i in 1:length(query_lines)) {
    query_line <- query_lines[i]
    midline_line <- midline_lines[i]
    subject_line <- subject_lines[i]
    
    start_pos <- (i - 1) * max_width + 1
    end_pos <- min(i * max_width, query_length)
    
    formatted_text <- paste0(formatted_text, 
                             sprintf("Query  %d  %s  %d\nMidln    %s%s\nSbjct  %d  %s  %d\n\n",
                                     start_pos, 
                                     paste(mapply(function(x, y) {
                                       if (x == y) {
                                         x
                                       } else {
                                         paste('<span style="background-color:red">', x, '</span>', sep='')
                                       }
                                     }, strsplit(query_line, '')[[1]], strsplit(subject_line, '')[[1]]), collapse=""),
                                     end_pos,
                                     strrep(" ", nchar(as.character(start_pos))), 
                                     midline_line,
                                     start_pos, 
                                     paste(mapply(function(x, y) {
                                       if (x == y) {
                                         y
                                       } else {
                                         paste('<span style="background-color:red">', y, '</span>', sep='')
                                       }
                                     }, strsplit(query_line, '')[[1]], strsplit(subject_line, '')[[1]]), collapse=""), 
                                     end_pos))
  }
  formatted_text <- paste0(formatted_text, "</pre>")
  
  formatted_text
}


# ===============================================================
# UI
# ===============================================================
ui <- tagList(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .navbar-nav .nav-item .nav-link {
        margin-left: 20px;
        margin-right: 20px;
      }
    ")),
    tags$style(HTML("
  .btn.disabled, .btn[disabled] {
    cursor: not-allowed !important;
    opacity: 0.5;
  }
")),
    tags$script(HTML("
$(document).on('click', '.fa-dna', function() {
  Shiny.setInputValue(
    'selected_protein_row',
    $(this).data('row'),
    {priority: 'event'}
  );
});
    "))
  ),
  
  
  ############################################################ UI 1. Database ############################################################ 
  page_navbar(
    title = tags$span("Protein Database", style = "color: black;"),
    theme = bs_theme(
      bootswatch = "flatly",
      base_font = font_google("Space Mono"),
      heading_font = font_google("Space Mono"),
      bg = "#FFFBEB",
      fg = "#000000"
    ),
    navbarMenu(
      title = tags$samp("Database", style = "color:black;"),
      icon = icon("info-circle", style = "color: black;"),
      
      nav_panel(
        title = tags$span("Protein database", style = "color: black;"),
        icon = icon("table"),
        
        page_sidebar(
          sidebar = sidebar(
            
            fileInput("upload_csv", "Upload your CSV file:", accept = ".csv"),
            selectInput(
              "data_source",
              "Select data source:",
              choices = c(
                "Uploaded CSV" = "csv",
                "Loaded session" = "session"
              ),
              selected = "csv"
            ),
            hr(),
            
            textInput("session_name", "Session name"),
            
            actionButton("save_session", "Save session", icon = icon("save"),
                         class = "btn-success"),
            
            selectInput("load_session", "Load session", choices = NULL),
            
            actionButton("merge_session", "Merge uploaded CSV",
                         icon = icon("plus"), class = "btn-warning"),
            
            actionButton("delete_session", "Delete session permanently",
                         icon = icon("trash"), class = "btn-danger"),
            
            hr(),
            
            checkboxGroupInput(
              "show_vars_protein",
              "Columns in Protein Database to show:",
              choices = NULL
            )
          ),
          
          fluidPage(
            
            # ---------------- INSTRUCTIONS ----------------
            card(
              card_header("Instructions"),
              tags$p(
                "You can upload a CSV file with any number of columns. ",
                "However, the following columns are recommended:"
              ),
              tags$p(
                strong("Recommended columns: "),
                "Specie, Genus, Gene, Annotation, Location, Coordinates, Strand, Protein, ",
                "N. AA, Subcellular location, Conserved domains, GO IDs, GO Terms, Sequence."
              ),
              tags$ul(
                tags$li("Protein: Protein ID (mandatory)"),
                tags$li("Sequence: Protein sequence (mandatory)"),
                tags$li("GO Terms: separate multiple terms using ;"),
                tags$li("GO IDs: separate multiple IDs using ;")
              )
            ),
            
            br(),
            
            # ---------------- TABLE ----------------
            card(
              card_header(
                tags$div(
                  style = "text-align: center; font-weight: bold; font-size: 20px;",
                  "Database"
                )
              ),
              div(style = "overflow-x: auto;", DTOutput("proteinTable")),
              
              fluidRow(
                column(4, downloadButton("downloadDATABASE", "Download Data"))
              )
            )
          )
        )
      ),
      ############################################################ UI 5. Proteins sequences ############################################################
      ############################################################
      # UI 5. Protein sequences (INDEPENDENT)
      ############################################################
      nav_panel(
        title = tags$span("HSP90 Protein sequences", style = "color: black;"),
        icon = icon("dna"),
        
        page_sidebar(
          sidebar = sidebar(
            
            selectInput(
              "data_source_seq",
              "Select data source:",
              choices = c(
                "Currently loaded dataset" = "current",
                "Saved session" = "session"
              ),
              selected = "current"
            ),
            
            conditionalPanel(
              condition = "input.data_source_seq == 'session'",
              selectInput(
                "load_session_seq",
                "Load session:",
                choices = NULL
              )
            ),
            
            hr(),
            
            fileInput(
              "upload_domains",
              "Upload CD-search results (hitsStandardResults)",
              accept = c(".txt", ".tsv")
            ),
            
            checkboxGroupInput(
              "domain_types",
              "Domain types to display:",
              choices = c(
                "specific",
                "superfamily",
                "non-specific"
              ),
              selected = c("specific", "superfamily")
            ),
            
            hr(),
            
            tags$div(
              style = "
    max-height: 220px;     
    overflow-y: auto;
    border: 1px solid #ced4da;
    border-radius: 5px;
    padding: 8px;
    background-color: #ffffff;
  ",
              
              checkboxGroupInput(
                "seq_protein_ids",
                "Proteins with conserved domains:",
                choices = NULL
              )
            ),
            
            fluidRow(
              column(
                6,
                actionButton(
                  "select_all_proteins",
                  "Select all",
                  icon = icon("check"),
                  class = "btn-success btn-sm"
                )
              ),
              column(
                6,
                actionButton(
                  "deselect_all_proteins",
                  "Deselect all",
                  icon = icon("times"),
                  class = "btn-danger btn-sm"
                )
              )
            ),
            
            hr(),),
          
          fluidPage(
            
            accordion(
              id = "acc_results",
              open = "Conserved domains", 
              
              accordion_panel(
                title = "Selected proteins & sequences (FASTA format)",
                icon = icon("dna"),
                uiOutput("protein_fasta_view"),
                br(),
                downloadButton(
                  "download_all_visible_fasta",
                  "Download all visible proteins (FASTA)",
                  icon = icon("download")
                )
              )),
            
            
            br(),
            
            card(
              id = "card_dominios",
              collapsible = TRUE,
              full_screen = TRUE,
              card_header(
                "Conserved domains",
                collapsible = TRUE
              ),
              plotOutput("protein_domains_plot", height = "400px"),
              br(),
              downloadButton(
                "download_domains_plot",
                "Download conserved domains plot",
                icon = icon("download")
              )
            )
          )
        )
      )
    ),
    
    
    ############################################# Phylogenetic tree ########################################################################    
    nav_panel(
      title = tags$samp("iTOL phylogenetic tree", style = "color:black;"),
      icon = icon("sitemap", lib = "font-awesome", style = "color:black;"),
      
      page_sidebar(
        sidebar = sidebar(
          
          selectInput(
            "tree_source",
            "Tree source:",
            choices = c(
              "Manual iTOL ID" = "manual",
              "Saved tree session" = "session"
            ),
            selected = "manual"
          ),
          
          hr(),
          
          conditionalPanel(
            condition = "input.tree_source == 'manual'",
            
            textInput(
              "itol_id_input",
              "iTOL Tree ID",
              placeholder = "e.g. R9KJpM"
            ),
            
            actionButton(
              "load_itol_tree",
              "Load phylogenetic tree",
              icon = icon("play"),
              class = "btn-success"
            )
          ),
          
          hr(),
          
          textInput(
            "tree_session_name",
            "Tree session name"
          ),
          
          actionButton(
            "save_tree_session",
            "Save tree session",
            icon = icon("save"),
            class = "btn-primary"
          ),
          
          hr(),
          
          conditionalPanel(
            condition = "input.tree_source == 'session'",
            
            selectInput(
              "load_tree_session",
              "Load saved tree session",
              choices = NULL
            )
          ),
          
          actionButton(
            "delete_tree_session",
            "Delete tree session",
            icon = icon("trash"),
            class = "btn-danger"
          )
        ),
        
        uiOutput("itol_frame")
      )
    ),
    ############################################################ UI. BLASTP ################################################################
    
    nav_panel(
      title = tags$samp("BLASTp / Diamond", style = "color:black"),
      icon = icon("search", style = "color: black;"),
      fluidPage(
        layout_columns(
          card(
            style = "height: auto;",  
            
            ### NUEVO ###: Selector para el programa
            selectInput(
              "program_selection",
              label = "Select Search Algorithm:",
              choices = c("BLASTp", "DIAMOND")
            ),
            
            selectInput(
              "species_selection",
              label = "Select Species",
              choices = NULL
            ),
            
            fileInput("fasta_file",
                      label = "Upload FASTA file",
                      accept = c(".fasta", ".txt")),
            actionButton("run_blast", "Run Analysis"),  
            textOutput("sequence_error_message")  
            
          ),
          card(
            card_header(
              tags$div(
                style = "text-align: center; font-weight: bold; font-size: 20px;",
                "Instruction for use"
              )
            ),
            verbatimTextOutput("blast_instructions")
          ),
          card(
            card_header(
              tags$div(
                style = "text-align: center; font-weight: bold; font-size: 20px;",
                "Alignments"
              )
            ),
            tags$div(
              DTOutput("Tableprotein") %>% withSpinner(),  
              uiOutput("downloadButtonUI"),  
              htmlOutput("alignments"),
              uiOutput("downloadSubjectButtonUI")
            )) ,
          col_widths = c(3,9,12)
        ) )          
    ),
    nav_panel(
      title = tags$samp("Sequence Aligner", style = "color:black"),
      icon = shiny::icon("dna", style = "color:black"),
      fluidPage(
        useShinyjs(),
        layout_columns(
          card(
            wellPanel(
              radioButtons("seq_type", "Sequence type:",
                           choices = c("Amino-acid" = "AA", "Nucleotid" = "DNA")),
              selectInput("align_method", "Choose the alignment method",
                          choices = c("ClustalW" = "ClustalW",
                                      "ClustalOmega" = "ClustalOmega",
                                      "Muscle" = "Muscle")),
              fileInput("align_fasta_file", "Upload your FASTA file",
                        accept = c(".fasta", ".fa")),
              actionButton("align_btn", "Align"),
              verbatimTextOutput("message")
            )
          ),
          
          card(
            fluidPage(
              card_header(
                tags$div(
                  style = "text-align: center; font-weight: bold; font-size: 20px;",
                  "Aligned sequences"
                )
              ),
              hidden(div(
                id = "loading_spinner",
                tags$p("Processing... Please wait.", style = "color: #28a745; font-weight: bold;"),
                tags$div(class = "spinner-border text-success", role = "status")
              )),
              hidden(uiOutput("alignment_output_container"))
            )),
          col_widths = c(12, 12)
        )
      )
    )
    
  )
)