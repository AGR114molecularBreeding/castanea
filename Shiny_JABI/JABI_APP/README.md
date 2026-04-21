## HSP90-_Fagaceae_: An Interactive R Shiny App for Centralizing Genomic Data of Heat Shock Protein in Fagaceae Species

This project is an interactive application built with R Shiny designed to centralize genomic protein data at the family and species level for user-defined studies. The tool was developed to address the fragmentation of genomic information across multiple servers. Its main objective is to provide a centralized, intuitive, and user-friendly interface that allows users to efficiently retrieve, visualize, and analyze genomic data. By integrating different data sources into a single platform, this application facilitates smoother and more efficient access to protein-related genomic information.


## 🌐 Availability
## 💻 Run Locally in R

## How to use locally?
### **Prerequisites**
- **R** (>= 4.0.0): [Download](https://cran.r-project.org/)
- **RStudio** (recommended): [Download](https://www.rstudio.com/products/rstudio/download/)

### **Steps**

<br>

1. **Download the folder [JABI_APP](https://github.com/AGR114molecularBreeding/castanea/tree/main/Shiny_JABI/JABI_APP)** , which contains:
   - **Core Scripts**:
     - [`ui.R`](https://github.com/AGR114molecularBreeding/castanea/blob/main/Shiny_JABI/JABI_APP/UI.R)
     - [`server.R`](https://github.com/AGR114molecularBreeding/castanea/blob/main/Shiny_JABI/JABI_APP/Server.R)
<br>

2. **Folder Structure**  
   To ensure proper functionality, maintain this exact structure:
   ```
   your_shiny_app/ # Main folder (any name)
   ├── app.R    # Main R script (combined UI + Server)
   ├── Data/    # Data folder
   │   ├── bin2/ # Folder contain BLASTp and DIAMOND algortihm
   │       └── Blastp/ # Folder
   │           └── blastp.exe # Blastp file
   │ 
   │       └── Diamond/ # Folder
   │           └── Diamond.exe # DIAMOND file
   │ 
   │   └── sessions/ # Folder that will contain the database session
   │       └──... (all .rds files)
   │
   │   └── Tree_sessions/ # Folder that will contain the phylogeny tree session
   │
   └────── proteomes/ # Images folder
           ├── 1.fasta # Genome coding sequence files
           ├── 2.fasta # Genome coding sequence files
           └──... (all FASTA files)


 3. **BLAST+ and DIAMOND Configuration** 🔧 (**Critical!**)
The executables for BLASTp and DIAMOND depend on your operating system (Windows, macOS, or Linux). Make sure to place the correct version inside:

Data/bin2/Blastp/
Data/bin2/Diamond/

#### 💻 **Windows Users**  
- Use `blastp.exe` from the downloaded files.
- It is also necessary to have the support libraries ncbi-vdb-md.dll and nghttp2.dll in the main folder for blastp.exe to run properly. These files are located in the Blastp folder.
- **Edit Line 926 in `server.R`**:  
  ```r  
  blastp_executable <- file.path("./blastp.exe")  # Write this    
  ```  
- **Edit Line 939 in `server.R`**:
  ```r  
  diamond_executable <- file.path("./data/bin2/Diamond/diamond.exe")  # Write this    
  ```
  
#### 🐧 **Linux Users**  
- Use `blastp` from the downloaded files.  
- Grant permissions:  
  ```bash  
  chmod +x blastp  
  ```  
- **Edit Line 926 in `server.R`**:  
  ```r  
  # blastp_executable <- file.path("./blastp.exe")  # Write this   
  ```
  
  - **Edit Line 939 in `server.R`**:
  ```r  
  diamond_executable <- file.path("./data/bin2/Diamond/diamond.exe")  # Write this    
  ```
  <br>

4. **Install R Packages**:
To run this tool it is necessary to install all the libraries that are indicated in the UI.R file.
   ```r
   # 1. Install CRAN packages
   
   cran_packages <- c(
   "bslib", "DT", "ggplot2", "grid", "leaflet", "openxlsx",
   "png", "sf", "shiny", "shinyalert",
   "shinycssloaders", "shinyjs", "xml2",
   "shinyWidgets", "readr"
   )

   for (pkg in cran_packages) {
   if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    }
   }



   # 2. Install Bioconductor packages
   
   if (!requireNamespace("BiocManager", quietly = TRUE)) {
   install.packages("BiocManager")
   }

   bioc_packages <- c("Biostrings","msa")

   for (pkg in bioc_packages) {
   if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
    }
   }

   # 3. Install msaR from GitHub (not available on CRAN or Bioconductor)

   
   if (!require("msaR", character.only = TRUE)) {
   if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
   }
   devtools::install_github("lcolladotor/msaR")
   }
   ```
<br>
  
5. **Run the App**
## To start the application, click the **Run** button at the top right:
![Figure 1](https://github.com/AGR114molecularBreeding/castanea/blob/main/HSP90/Wiki%20figures/7.PNG)

## Once clicked, the app will launch as shown below: 
![Figure 2](https://github.com/AGR114molecularBreeding/castanea/blob/main/HSP90/Wiki%20figures/8.PNG)
