# Load required packages
library(shiny)
library(msa)
library(msaR)
library(Biostrings)
library(DT) ##For interative datatable 
library(openxlsx)
library(shinyjs)
library(sf) #Neccesary for the interative map
library(leaflet) #Package to add maps 
library(shinycssloaders)
library(shinyalert) ## Package for blastp end
library(ggplot2) ##Graphics display package
library(grid)
library(png) ##Package that allows to display image
library(xml2)        ## Package for work with blast outmf 5
library(bslib)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 80)
options(shiny.maxRequestSize = 2*1024^2)  # This limit the size of uploaded file 
final_dataframe <- read.csv("./final_dataframe.csv")  #Dataframe wwith all information
filo <- file.path("./atodos.png")    #Phylogenetic PNG



# =============================================================================================
# PHYLOGENETIC TREE MENU
# =============================================================================================
#Create the image with ggplot2 and adjust the size
gg_image <- function(file_path, ranges) {
  img <- readPNG(file_path)
  g <- rasterGrob(img, interpolate = TRUE)
  img_width <- ncol(img)
  img_height <- nrow(img)
  
  ggplot() +
    annotation_custom(g, xmin = 0, xmax = img_width, ymin = 0, ymax = img_height) +
    coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, img_width)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, img_height)) +
    theme_void()  
}

# =============================================================================================
# LIST OF TREE IMAGES AND LICENSE
# =============================================================================================
image_data <- list(
  `Castanea crenata` = list(
    paths = paste0("./Images/Ccrenata_", 1:4, ".jpg"),
    authors = c("Yves SPM", "Trap Hers", "Emmanuel Bouchard", "Andrzej Konstantynowicz"),
    source = "Pl@ntNet",
    license = c("https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Castanea dentata` = list(
    paths = paste0("./Images/Cdentata", 1:3, ".jpg"),
    authors = c("M Amy", "Roby Lupi", "Willow Coville"),
    source = "Pl@ntNet",
    license = c("https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Castanea mollissima` = list(
    paths = paste0("./Images/Cmollissima", 1:3, ".jpg"),
    authors = c("Dieter Albrecht","Brett Bissell","Wenjing Guo"),
    source = "Pl@ntNet" ,
    license = c("https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Castanea sativa` = list(
    paths = paste0("./Images/Csativa", 1:4, ".jpg"),
    authors = c("Alain Bigou","Dominique Ponton","Alain Bigou","Alain Bigou"),
    source = "Pl@ntNet" ,
    license = c("https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/",
                "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Castanopsis hystrix` = list(
    paths = paste0("./Images/Chystrix", 1:3, ".jpg"),
    authors = c("yoavdanielbarness", "hhn126", "pasangk"),
    source = "Natusfera",
    license = c("https://creativecommons.org/licenses/by/4.0/", 
                "https://creativecommons.org/licenses/by-sa/4.0/", 
                "https://creativecommons.org/licenses/by-sa/4.0/")
  ),
  `Quercus dentata` = list(
    paths = paste0("./Images/Qdentata", 1:3, ".jpg"),
    authors = c("Kai Best", "Kai Best", "Зобенко Андрей"),
    source = "Pl@ntNet",
    license =  c("https://creativecommons.org/licenses/by-nc/4.0/",
                 "https://creativecommons.org/licenses/by-nc/4.0/",
                 "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Quercus ilex` = list(
    paths = paste0("./Images/Qilex", 1:3, ".jpg"),
    authors = c("Alain Bigou", "Dieter Albrecht", "Dieter Albrecht"),
    source = "Pl@ntNet",
    license =  c("https://creativecommons.org/licenses/by-nc/4.0/",
                 "https://creativecommons.org/licenses/by-nc/4.0/",
                 "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Quercus lobata`= list(
    paths= paste0("./Images/Qlobata", 1:3, ".jpg"),
    authors = c("Dieter Albrecht", "Angela Zhou", "Loidi Javier"),
    source = "Pl@ntNet",
    license =  c("https://creativecommons.org/licenses/by-nc/4.0/",
                 "https://creativecommons.org/licenses/by-nc/4.0/",
                 "https://creativecommons.org/licenses/by-nc/4.0/")
  ),
  `Quercus robur` = list(
    paths = paste0("./Images/Qrobur", 1:3, ".jpg"),
    authors = c("Alain Bigou", "Andrzej Konstantynowicz", "bergvelt"),
    source = "Pl@ntNet",
    license = c("https://creativecommons.org/licenses/by-sa/4.0/",
                "https://creativecommons.org/licenses/by-sa/4.0/",
                "https://creativecommons.org/licenses/by-sa/4.0/")
  ),
  `Quercus suber` = list( 
    paths = paste0("./Images/Qsuber", 1:3, ".jpg"),
    authors = c("Anna LL", "Viktor Nettoyeur", "Almudena Jiménez"),
    source = "Pl@ntNet",
    license = c("https://creativecommons.org/licenses/by-sa/4.0/",
                "https://creativecommons.org/licenses/by-sa/4.0/",
                "https://creativecommons.org/licenses/by-sa/4.0/")
  ),
  
  `Quercus variabilis` = list(
    paths= paste0("./Images/Qvariabilis", 1:3, ".jpg"),
    authors = c("harum.koh", "Cheng Te Hsu", "harum.koh"),
    source = "Pl@ntNet",
    license = c("https://creativecommons.org/licenses/by-sa/4.0/",
                "https://creativecommons.org/licenses/by-sa/4.0/",
                "https://creativecommons.org/licenses/by-sa/4.0/")
  )
)

# =============================================================================================
# DEFINE THE POLYGONS FROM SHP FILE
# =============================================================================================
# load the shapefiles with sf
shapefiles <- list(
  "Castanea crenata" = st_read("./Maps/Ccrenata.shp"),
  "Castanea dentata" = st_read("./Maps/Cdent.shp"),
  "Castanea mollissima" = st_read("./Maps/Cmolli.shp"),
  "Castanea sativa" = st_read("./Maps/SATIVA.shp"),
  "Castanopsis hystrix" = st_read("./Maps/Chystrix.shp"),
  "Quercus dentata" = st_read("./Maps/Qdent.shp"),
  "Quercus ilex" = st_read("./Maps/QilexARTICULO.shp"),
  "Quercus lobata" = st_read("./Maps/Quercus_lobata_disuelto.shp"),
  "Quercus robur" = st_read("./Maps/Querc_ROB_ART.shp"),
  "Quercus suber" = st_read("./Maps/Qerc_SUBER__ART.shp"),
  "Quercus variabilis" = st_read("./Maps/Qvari.shp")
)
# =============================================================================================
# MENU BLASTP
# =============================================================================================

###############################  Formatting the alignment  ############################################
format_alignment <- function(alignment) {
  # Change getNodeSet and xmlValue to xml_find_all and xml_text
  query_seq <- xml_text(xml_find_all(alignment, "./Hsp_qseq")[1])
  midline <- xml_text(xml_find_all(alignment, "./Hsp_midline")[1])
  subject_seq <- xml_text(xml_find_all(alignment, "./Hsp_hseq")[1])
  
  # Define maximum width for displaying characters per line
  max_width <- 60
  
  # Function for splitting a sequence into lines of maximum length
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
  
  # Get the total length of the query sequence
  query_length <- nchar(query_seq)
  
  # Create text formatted in HTML format
  formatted_text <- "<pre>"
  for (i in 1:length(query_lines)) {
    query_line <- query_lines[i]
    midline_line <- midline_lines[i]
    subject_line <- subject_lines[i]
    
    #Calculate start and end item numbers
    start_pos <- (i - 1) * max_width + 1
    end_pos <- min(i * max_width, query_length)
    
    # Add formatted line to HTML text
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
                                     strrep(" ", nchar(as.character(start_pos))),  #Add spaces based on the number of digits in start_pos
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
  
  # Return formatted HTML text
  formatted_text
}

# ===================================================================================================================================================
######################################################### UI OF THE APP #########################################################
# ===================================================================================================================================================

ui <- page_navbar(
  title = tags$span("Protein Database of Fagaceae :", style = "color: black;"),  
  theme = bs_theme(
    bootswatch = "flatly",
    base_font = font_google("Roboto"),
    heading_font = font_google("Montserrat"),
    primary = "#28a745",
    bg = "#e7e2da",
    fg = "#000000"  
  ),
  
  # Add more space betwween the titles of the nav_panels which are in the navbarMenu
  header = tags$style(HTML("
    .navbar-nav .nav-item .nav-link {
      margin-left: 20px;  /* Adds space to the left of the titles */
      margin-right: 20px; /* Adds space to the right of the titles */
    }
  ")),
  
  
  ########################################################### UI 1. RATIONALE ####################################################################
  
  nav_panel(
    title = tags$span("Rationale", style = "color: black;"),
    icon = icon("lightbulb", style = "color: black;"),
    
    mainPanel(
      style = "width: 100%;",
      
      layout_columns(
        card(
          card_header(
            tags$div(
              style = "text-align: center; font-weight: bold; font-size: 20px;",
              "Rationale"
            )
          ),
          div(
            style = "background-color: white; padding: 15px; border-radius: 10px;", # White background and rounded edges
            uiOutput("rationale_text")
          )
        ),
        
        card(
          
          div(
            style = "background-color: white; padding: 15px; border-radius: 10px;", # White background and rounded edges
            uiOutput("funding_text")
          )
        ),
        
        col_widths = c(12) 
      )
    )
  ),
  
  ############################################################ UI 2. TAXONOMIC INFO ####################################################################
  navbarMenu(
    title = tags$samp("Fagaceae information", style = "color:black;"),
    icon = icon("info-circle", style = "color: black;"),
    
    #submenu for Taxonomic information
    nav_panel(
      title = tags$span("General information", style = "color: black;"),  # Add black color to tab title
      icon = icon("info-circle", style = "color: black;"),
      page_sidebar(
        sidebar = sidebar(
          div(
            style = "padding-left: 19px;",
            tags$style(HTML("#var_species + .selectize-control .selectize-input, #var_species + .selectize-control .selectize-dropdown-content, #var_species + .selectize-control .selectize-dropdown-content .option:hover {font-style: italic;}")),
            style = "font-size: 18px;", 
            selectInput("var_species", label = "Choose a species", choices = unique(final_dataframe$SpecieII), width = "auto")
          )
        ),
        mainPanel(
          style = "width: 100%;",  # Ensures that the container occupies the entire available width
          
          layout_columns(
            card(
              card_header(
                tags$div(
                  style = "text-align: center; font-weight: bold; font-size: 20px;",
                  "Taxonomic Info"
                )
              ),
              uiOutput("info_box"),
              uiOutput("css")
            ),
            
            card(
              card_header(
                tags$div(
                  style = "text-align: center; font-weight: bold; font-size: 20px;",
                  "Images"
                )
              ),
              div(
                style = "text-align: center;",  # Center the image and the text below it
                imageOutput("species_image", width = "100%", height = "auto"),  # Dynamic size adjustment
                actionButton(
                  "change_image_btn", NULL,
                  icon = icon("arrow-right"),
                  style = "position: absolute; top: 50%; right: 10px; transform: translateY(-50%);" 
                  # This allow the button to stay in the mid-right part of the image and maade it immobile
                ),
                htmlOutput("image_metadata", style = "margin-top: 10px;")  # Top margin between image and metadata
              )
              
            ),
            
            card(
              card_header(
                tags$div(
                  style = "text-align: center; font-weight: bold; font-size: 20px;",
                  "Habitat map"
                )
              ),
              leafletOutput("species_map", width = "100%", height = "400px")
            ),
            
            col_widths = c(12, 6, 6)
          )
        )
      )
    ),
    
    ############################################################ UI 3. HSP90 Fagaceae database ############################################################
    nav_panel(
      title = tags$span("HSP90 Protein database", style = "color: black;"),
      icon = icon("table"),
      page_sidebar(
        sidebar = sidebar(
          div(
            checkboxGroupInput("show_vars_protein", "Columns in Protein Database to show:",
                               choices = c("Specie", "Protein", "Annotation", "Gene", "Chr",
                                           "Genus", "Conserved domains", "N. AA"),
                               selected = c("Specie", "Protein", "Annotation", "Gene", "Chr", "Genus", "Conserved domains", "N. AA"))
          )
        ),
        fluidPage(
          layout_columns(
            card(
              card_header(
                tags$div(
                  style = "text-align: center; font-weight: bold; font-size: 20px;",
                  "HSP90 database"
                )
              ),
              div(style = "overflow-x: auto;", DTOutput("proteinTable")),
              conditionalPanel(
                condition = "input.show_vars_protein.length > 0",  # If there are selected columns
                downloadButton('downloadDATABASE', 'Download Data')
              )
            ),
            col_widths = c(12)
          )
        )
      )
    ),
    
    ############################################################ UI 4. HSP90 proteins sequences ############################################################
    
    nav_panel(title = tags$span("HSP90 Protein sequences", style = "color: black;"),  # Add black color to tab title
              icon = icon("font", class = "fa-solid"),  # Letters icon (ABC)
              page_sidebar(
                sidebar = sidebar(
                  div(
                    tags$style(HTML("#Species + .selectize-control .selectize-input, #Species + .selectize-control .selectize-dropdown-content, #Species + .selectize-control .selectize-dropdown-content .option:hover {font-style: italic;}")),
                    selectInput("Species", 
                                label = "Choose a species", 
                                choices = unique(final_dataframe$Specie)),
                    uiOutput("protein_selection")
                  )
                ),
                layout_columns(
                  card(
                    card_header(
                      tags$div(
                        style = "text-align: center; font-weight: bold; font-size: 20px;",
                        "Amino acid sequence"
                      )
                    ),
                    uiOutput("output_sequence"),
                    downloadButton('downloadamino', 'Download sequence', class = "btn btn-default"),
                    uiOutput("downloadall_button")),
                  
                  col_widths = c(12))
              )  
    )
  ),
  ############################################################ UI 5. PHYLOGENETIC TREE ############################################################
  
  nav_panel(
    title = tags$samp("Phylogenetic tree", style = "color:black;"),
    icon = icon("tree", lib = "font-awesome", style = "color:black;"),
    fluidPage(
      layout_columns(
        card(
          card_header(
            tags$div(
              style = "text-align: center; font-weight: bold; font-size: 20px;",
              (p("Phylogenetic tree of HSP90 proteins from the main species of the ", tags$em("Fagaceae"), " family"))
            )
          ),
          div(style = "margin-top: 20px; text-align: center;",
              div(style = "display: inline-block;",   #  Centered container
                  plotOutput("gg_image_output", height = "1000px", width = "1400px",  # Adjust height and width 
                             dblclick = "plot1_dblclick",
                             brush = brushOpts(
                               id = "plot1_brush",
                               resetOnNew = TRUE
                             ))
              )
          ),
          uiOutput("image_caption"),   # Image caption
          uiOutput("inforr_box"),       # Additional information
          div(style = "text-align: center; margin-top: 20px;",  # Center the button
              uiOutput("download_button"))
        )))
  ),
  
  ############################################################ UI 6. BLASTP ################################################################
  
  nav_panel(
    title = tags$samp("BLASTp", style = "color:black"),
    icon = icon("search", style = "color: black;"),
    fluidPage(
      layout_columns(
        card(
          selectInput("species_selection", label = "Select Species", choices = c("Castanea crenata", "Castanea dentata", 
                                                                                 "Castanea mollissima H7", "Castanea mollissima N11",
                                                                                 "Castanea mollissima vanuxem", "Castanea sativa",
                                                                                 "Castanopsis hystrix", "Quercus dentata", 
                                                                                 "Quercus lobata", "Quercus robur", 
                                                                                 "Quercus suber", "Quercus variabilis")),
          
          fileInput("fasta_file",
                    label = "Upload FASTA file",
                    accept = c(".fasta", ".txt")),
          actionButton("run_blast", "Run Blast"),  # Button to start the Blastp
          textOutput("sequence_error_message")  # Messages will be displayed here
          
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
            DTOutput("Tableprotein") %>% withSpinner(),  # Spinner for the table
            uiOutput("downloadButtonUI"),  # Replace button with uiOutput
            htmlOutput("alignments"),
            uiOutput("downloadSubjectButtonUI")
          )) ,
        col_widths = c(2,10,12)
      ) )          
  ),
  
  
  ############################################################ UI 7. ALIGNMENT ################################################################
  nav_panel(
    title = tags$samp("Sequence Aligner", style = "color:black"),
    icon = shiny::icon("dna", style = "color:black"),
    fluidPage(
      useShinyjs(),  # Activate shinyjs for the dinamic control
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
            hidden(uiOutput("alignment_output_container"))  # Container hidden by default
          )),
        col_widths = c(12, 12)
      )
    )
  )
)
