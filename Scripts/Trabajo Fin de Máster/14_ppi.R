#### ================================
#### 1. Configuraciones previas
#### ================================

# Directorios de entrada y salida que debe elegir el usuario:

# Directorio con los resultados del análisis de expresión diferenciall
input <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/DE/"
# Directorio de salida
output <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/ppi/"

# Inicialización del temporizador de tiempo de ejecución
t0 <- proc.time()

# Crear o sobreescribir carpeta de salida
if (dir.exists(output)) {
  unlink(output, recursive = TRUE, force = TRUE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Importar librerías
library(STRINGdb)
library(dplyr)

#### ================================
#### 2. Carga de datos
#### ================================

# Coger las tablas con los DEGs de cada comparación
comparison_list_DESeq2 <- list.files(path = paste0(input, "DESeq2/"), pattern = "_deg_all", recursive = TRUE, full.names  = TRUE)
comparison_list_edgeR <- list.files(path = paste0(input, "edgeR/"), pattern = "_deg_all", recursive = TRUE, full.names  = TRUE)
comparison_list <- c(comparison_list_DESeq2, comparison_list_edgeR)

#### ================================
#### 3. Análisis de interacciones proteína-proteína
#### ================================

# Seleccionar especie homo sapiens
species_taxid <- 9606

# Crear objeto STRING
string_db <- STRINGdb$new(
  version      = "12",
  species      = species_taxid,
  score_threshold = 400,
  input_directory = ""
)

# Ejecutar análisis de cada comparación
for (comparison_file in comparison_list) {
  
  # Nombre de la comparación
  comparison <- gsub("_deg_all.csv", "", basename(comparison_file))
  groups <- strsplit(comparison, "_vs_")[[1]]
  
  # Mensaje de progreso
  cat("\n", "Realizando análisis PPI de:", comparison, "\n")
  
  # Cargar tabla de DEGs
  table <- read.csv2(comparison_file, header = TRUE, check.names = FALSE)
  
  # Filtrar todos los genes cuyo pvalor o log2fc sea Na y sustituir pvalue=0 por el minimo
  table <- table[!is.na(table$pvalue) & !is.na(table$log2FoldChange), ]
  table$pvalue[table$pvalue == 0] <- .Machine$double.xmin
  
  # Crear tabla de genes ordenados por el valor absoluto de score = −log10(p value) x (log2FC)
  p_values <- table$pvalue
  log2fc <- table$log2FoldChange
  score <- -log10(p_values)*log2fc
  table$score <- score
  table <- table[order(abs(table$score), decreasing = TRUE), ]
  
  # Quedarse con los 2000 de mayor score
  table <- table[1:2000, ]
  
  # Eliminar duplicados
  table <- table[!duplicated(table$gene_id), ]
  
  # Coger lista de DEGs totales, de up genes y de down genes
  gene_set_list <- list()
  deg_genes <- table
  up_genes <- table[table$log2FoldChange > 1, ]
  down_genes <- table[table$log2FoldChange < -1, ]
  
  gene_set_list[["all"]] <- deg_genes
  gene_set_list[["up"]] <- up_genes
  gene_set_list[["down"]] <- down_genes
  
  for (gene_set_name in names(gene_set_list)) {
    
    # Escoger la lista de genes
    degs <- as.data.frame(gene_set_list[[gene_set_name]])
    
    # Mapear genes a IDs de STRING
    degs_mapped <- string_db$map(
      degs,
      "gene_id",
      removeUnmappedRows = TRUE
    )
    
    # Seleccionar el conjunto de proteínas para la red
    hits <- degs_mapped$STRING_id
    
    # Obtener interacciones PPI desde STRING
    ppi <- string_db$get_interactions(hits)
    ppi <- ppi %>%
      distinct(from, to, combined_score, .keep_all = TRUE)
    
    # Coger correpondencias entre IDs
    map_ids <- degs_mapped[, c("STRING_id", "gene_id")]
    colnames(map_ids) <- c("STRING_id", "gene")
    
    # Crear tabla final
    ppi_table <- ppi %>%
      
      # Renombrar columnas que contienen los IDs de las proteínas que interaccionan
      rename(node1_protein = from,
             node2_protein = to,
             score         = combined_score) %>%
      
      # Añadir nombre de gen para node1
      left_join(map_ids,
                by = c("node1_protein" = "STRING_id")) %>%
      rename(node1_gene = gene) %>%
      
      # Añadir nombre de gen para node2
      left_join(map_ids,
                by = c("node2_protein" = "STRING_id")) %>%
      rename(node2_gene = gene) %>%
      
      # Reordenar columnas
      dplyr::select(node1_protein, node2_protein,
             node1_gene,    node2_gene,
             score)
    
    # Ordenar tabla en función del score
    ppi_table <- ppi_table[order(ppi_table$score, decreasing = TRUE), ]
    
    # Eliminar el prefijo de los IDs de proteínas
    ppi_table$node1_protein <- gsub("9606.", "", ppi_table$node1_protein)
    ppi_table$node2_protein <- gsub("9606.", "", ppi_table$node2_protein)
    
    # Exportar resultados
    out_ppi <- paste0(output, comparison, "/")
    dir.create(out_ppi, recursive = TRUE)
    write.csv2(ppi_table, paste0(out_ppi, comparison, ".", gene_set_name, "_ppi.csv"), row.names = FALSE)
  
  }
}

# Mostrar tiempo de ejecución
t_total <- (proc.time() - t0)/60
cat("Script ejecutado correctamente.\n")
cat("Tiempo de ejecución:", round(t_total["elapsed"], 2), "minutos\n")