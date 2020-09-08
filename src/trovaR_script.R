"_                       _____
| |                     |  __ \
| |_ _ __ _____   ____ _| |__) |
| __| '__/ _ \ \ / / _` |  _  /
| |_| | | (_) \ V / (_| | | \ \
\___|_|  \___/ \_/ \__,_|_|  \_\
"

# ---------------
# Copyright © 2020 Daniel León Periñán (danielleon.tk) All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#   * Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above copyright notice, this
#     list of conditions and the following disclaimer in the documentation and/or
#     other materials provided with the distribution.
#   * Neither the name of Daniel León Periñán nor the names of its contributors may be used to
#     endorse or promote products derived from this software without specific prior
#     written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# ------------------------------
### PROGRAM PARAMETERS
# Directory containing txt SGA files
DI.ROOT <- "~/Escritorio/SGA_Alfonso/data/allPictures/analysis/"
# List of phenotypes
PH.LIST <- c("WT", "1020_d")
# List of conditions
# (if a name begins with number, must be preceded by X (ex: 30 -> "X30"))
CD.LIST <- c("FOR")


### Significance values
# p-value for all analyses
P.VALUE <- c(1, 0.001, 0.005, 0.05)
# log2FC frontier
ABS.LFC <- c(0, 2, 1, 0.3)
# Identifier for each of the values above
VAL.NAM <- c("none", "high", "med", "low")


### Plotting settings
# Show heatmap
PL.HEAT <- TRUE
# Heatmap size (default c(7,12))
PL.SIZE <- c(7, 12)
# Show scatterplots
PL.SCTT <- TRUE
# Folder to store plots (default: DI.ROOT)
PL.SAVE <- DI.ROOT


### Storage of value tables
# Create a table to save results
FS.CREA <- TRUE
# Folder to store tables (default: DI.ROOT)
FS.SAVE <- DI.ROOT

# -------------------------------
# MAIN PROGRAM

### Required libraries for loading
packages <-
  c("dplyr",
    "tidyr",
    "textshape",
    "gplots",
    "ggplot2",
    "mixtools",
    "stringr")
# fct for libraries loading and/or installing
libraries <- function(packages) {
  for (package in packages) {
    #checks if package is installed
    if (!require(package, character.only = TRUE)) {
      #If package does not exist, then it will install
      install.packages(package, dependencies = TRUE)
      #Loads package
      library(package, character.only = TRUE)
    }
  }
}

# call function
libraries(packages)

# Create the main dataframe
SGA_data <- data.frame()

# Set directory to given root path
setwd(DI.ROOT)

# List directories in given root path
directories <-
  list.dirs(path = ".",
            full.names = TRUE,
            recursive = TRUE)

# Collapse the list of conditions to perform search
conditions = CD.LIST

# Iterate over conditions
for (cond in conditions) {
  # Iterate over directories
  for (dir in directories) {
    if (str_count(dir, "/") < 3) {
      next
    }
    
    strain_pattern = paste(PH.LIST, collapse="|")
    if (!grepl(strain_pattern, dir)) {
      next
    }
    
    setwd(dir)
    for (f in list.files(pattern = "*_non-filtered.txt")) {
      dir <- gsub('^\\.|\\.$', '', gsub("/", ".", dir))
      colname_current <- paste(cond, dir, sep = "")
      single_data <- read.csv(f, header = TRUE, sep = '\t')
      single_data <-
        single_data %>% filter(gene != "#N/A") %>% filter(gene != "")
      single_data <- single_data %>%
        select(gene, ORF, plate, Ring, starts_with(cond))
      single_data <- single_data %>%
        mutate(Ring = ifelse(Ring > 6, 6, Ring))
      single_data <- single_data %>%
        group_by(plate, Ring) %>%
        mutate(!!cond := !!as.name(cond) / median(!!as.name(cond))) %>%
        ungroup()
      single_data <-
        single_data %>% dplyr::rename(!!colname_current := !!as.name(cond))
      current <- single_data %>% select(!!colname_current)
    }
    if (dim(SGA_data)[1] == 0) {
      SGA_data <- single_data
    } else {
      SGA_data <- cbind(SGA_data, current)
    }
    setwd(DI.ROOT)
  }
}

# Calculate pvalues and Fold Changes for each gene
for (cond in CD.LIST) {
  SGA_cond <- SGA_data %>%
    select(matches(cond))
  wt_cols <- SGA_cond %>%
    select(matches(PH.LIST[1]))
  mutant_cols <- SGA_cond %>%
    select(matches(PH.LIST[2]))

  pval_name <- paste("pvalue", cond, sep = ".")
  l2fc_name <- paste("log2FC", cond, sep = ".")
  wmnn_name <- paste("wmean", cond, sep = ".")
  mmnn_name <- paste("mmean", cond, sep = ".")
  for (i in 1:nrow(SGA_data))
  {
    data_WT = as.numeric(wt_cols[i,])
    data_mut = as.numeric(mutant_cols[i,])
    t.test(data_WT, data_mut)$p.value -> pval
    SGA_data[i, pval_name] <- pval
    SGA_data[i, l2fc_name] <-
      log2((mean(data_mut) + 0.001) / (mean(data_WT) + 0.001))
    SGA_data[i, wmnn_name] <- mean(data_WT)
    SGA_data[i, mmnn_name] <- mean(data_mut)
  }
}

# Set the frontiers for each phenotype (H,S,L classification)
for (cond in CD.LIST) {
  wmnn_name <- paste("wmean", cond, sep = ".")
  mmnn_name <- paste("mmean", cond, sep = ".")
  wnph_name <- paste("wt.pheno", cond, sep = ".")
  mtph_name <- paste("mt.pheno", cond, sep = ".")
  lb_mt <-
    quantile(SGA_data[, mmnn_name], 0.25) - 1.5 * IQR(SGA_data[, mmnn_name])
  lb_wt <-
    quantile(SGA_data[, wmnn_name], 0.25) - 1.5 * IQR(SGA_data[, wmnn_name])
  SGA_data <- SGA_data %>%
    mutate(
      !!wnph_name := ifelse(
        !!as.name(wmnn_name) > lb_wt,
        "H",
        ifelse(!!as.name(wmnn_name) < 0.3, "L", "S")
      ),
      !!mtph_name := ifelse(
        !!as.name(mmnn_name) > lb_mt,
        "H",
        ifelse(!!as.name(mmnn_name) < 0.3, "L", "S")
      )
    )
}

# Write non-filtered results to text file
if (FS.CREA) {
  setwd(FS.SAVE)
  write.csv(SGA_data, "results_non_filtered.csv")
  write.csv(
    SGA_data %>% select(-starts_with(CD.LIST)) %>% select(!plate) %>% select(!Ring),
    "results_non_filtered_minified.csv"
  )
}

# Plotting of heatmap
heat_plotting <- function(route, name, df) {
  if (PL.HEAT) {
    # Set gene as rowname
    rownames(df) <- df$gene
    x  <- as.matrix(df %>% select(!gene))

    # Select directory
    setwd(route)

    # Generate the heatmap to htmp
    pdf(paste0("heatmap_", name),
        width = PL.SIZE[1],
        height = PL.SIZE[2])
    heatmap.2(
      x,
      col = bluered(75),
      trace = "none",
      Colv = FALSE,
      dendrogram = "row",
      distfun = dist,
      margins = c(10, 10),
      scale = "none",
      breaks = seq(0, 1.3, length.out = 76)
    ) -> htmp
    dev.off()

  }
}

# Write filtered results to text file
save_value_table <- function(route, name, df) {
  if (FS.CREA) {
    setwd(FS.SAVE)
    write.csv(df, paste0("results_",name))
  }
}

for (i in 1:length(P.VALUE)) {
  # Deselect uninteresting columns and filter by pvalue and log2FC
  SGA_analysis <- SGA_data %>%
    select(-plate, -Ring, -contains("sign")) %>%
    filter(select(., contains("pvalue")) <= P.VALUE[i] &
             abs(select(., contains("log2FC"))) >= ABS.LFC[i])
  SGA_analysis_text <-
    SGA_analysis %>% select(-starts_with(CD.LIST))
  SGA_analysis <- SGA_analysis %>%
    select(
      -contains("mean"),-contains("sign"),-contains("pheno"),-ORF,-contains("pvalue"),-contains("log2FC")
    )
  heat_plotting(PL.SAVE, VAL.NAM[i], SGA_analysis)
  save_value_table(FS.SAVE, VAL.NAM[i], SGA_analysis_text)
}

