# HSP90-Fagaceae
## HSP90-Fagaceae: An Interactive R Shiny App for Centralizing Genomic Data of Heat Shock Protein in Fagaceae Species

Welcome to HSP90-Fagaceae: An Interactive R Shiny App for Centralizing Genomic Data of Heat Shock Protein in Fagaceae Species. This tool was designed to visualize and analyze genomic information on the main Fagaceae species currently available in fragmented form on servers. The purpose is to offer a user-friendly and intuitive interface that allows users to obtain data centrally, thus enabling more fluid and efficient access to information. 
---
## 🌐 Availability
This app is available in two ways:
### **Option 1: Online**
[![Live Demo](https://img.shields.io/badge/HSP90_Fagaceae-Available-green)](https://hsp90.ext.uco.es/)  

### **Option 2: Locally in R** 
Run it locally for full customization or offline use:  

### **Prerequisites**
- **R** (>= 4.0.0): [Download](https://cran.r-project.org/)
- **RStudio** (recommended): [Download](https://www.rstudio.com/products/rstudio/download/)
### **Steps**

1. **Download Files**:
   - **Core Scripts**:
     - [`ui.R`](https://github.com/your-username/your-repo/raw/main/ui.R)
     - [`server.R`](https://github.com/your-username/your-repo/raw/main/server.R)
   - **Essential Folders**:
     - 🧬 [Coding Sequences](https://github.com/AGR114molecularBreeding/castanea/tree/main/HSP90/Proteomes)
     - 🌿 [Images](https://github.com/AGR114molecularBreeding/castanea/tree/main/HSP90/Images)
     - 🗺️ [Maps](https://github.com/AGR114molecularBreeding/castanea/tree/main/HSP90/Maps)
     - 📊 [HSP90 Data](https://github.com/your-username/your-repo/tree/main/data)
     - ⚙️ [Blastp](https://github.com/AGR114molecularBreeding/castanea/tree/main/HSP90/Blastp)

2. **Folder Structure**  
   To ensure proper functionality, maintain this exact structure:
   ```
   your_shiny_app/              # Main folder (any name)
   ├── app.R                    # Main R script (combined UI + Server)
   ├── blastp/blastp.exe        # Blastp file for Windows or Linux
   ├── Castanea_crenata.fasta
   ├── Castanea_dentata.fasta
   ├── ... (all FASTA files)
   │
   ├── images/
   │   ├── Ccrenata_1.jpg
   │   ├── Ccrenata_2.jpg
   │   ├──... (all images files)
   │   └── atodos.png
   │
   ├── data/
   │   ├── final_dataframe.csv
   │   └── Tree.Newick
   │
   └── maps/
       ├── Ccrenata.cpg
       ├── Ccrenata.dbf
       ├── ... (all map files)
   ```
3. **BLAST+ Configuration** 🔧 (**Critical!**)  
#### 💻 **Windows Users**  
- Use `blastp.exe` from the downloaded files.  
- **Edit Line 411 in `server.R`**:  
  ```r  
  blastp_executable <- file.path("./blastp.exe")  # Write this    
  ```  

#### 🐧 **Linux Users**  
- Use `blastp` from the downloaded files.  
- Grant permissions:  
  ```bash  
  chmod +x blastp  
  ```  
- **Edit Line 411 in `server.R`**:  
  ```r  
  # blastp_executable <- file.path("./blastp.exe")  # Write this   
  ```  
4. **Install R Packages**:
To run this tool it is necessary to install all the libraries that are indicated in the UI.R file
   ```r
   install.packages(c("shiny", "BiocManager", "ggplot2", "seqinr"))
   BiocManager::install("Biostrings")  # For sequence analysis
   ```

5. **Run the App**

