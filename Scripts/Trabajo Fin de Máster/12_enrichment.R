#### ================================
#### 1. Configuraciones previas
#### ================================

# Directorios de entrada y salida que debe elegir el usuario:

# Directorio con los resultados del análisis de expresión diferencial
input <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/DE/"
# Directorio de salida
output <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/enrichment/"
# Directorio con la tabla de correspondencias entre Ensembl IDs y Entrez IDs
corresp_table_dir <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"

# Inicialización del temporizador de tiempo de ejecución
t0 <- proc.time()

# Crear o sobreescribir carpeta de salida
if (dir.exists(output)) {
  unlink(output, recursive = TRUE, force = TRUE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Importar librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

#### ================================
#### 2. Carga de datos
#### ================================

# Coger todos los genes de estudio (universo)
all_comparison_list <- list.files(path = paste0(input, "DESeq2/"), pattern = "_all.csv", recursive = TRUE)
all_table <- read.csv2(paste0(input, "DESeq2/", all_comparison_list[1]), header = TRUE, check.names = FALSE)

# Convertir IDs de Ensembl a IDs de ENTREZ para los genes del universo
entrez2ensembl <- read.table(file = paste0(corresp_table_dir,"entrez2ensembl.txt"), sep = "\t", header = TRUE)
all_table_entrez <- merge(all_table, entrez2ensembl, by.x="gene_id", by.y="Ensembl", all.x=TRUE)
all_gene_entrez <- as.character(all_table_entrez$Entrez)

# Coger tablas con los DEGs de cada comparación
deg_comparison_list_edgeR <- list.files(path = paste0(input, "edgeR/"), pattern = "_deg.csv", recursive = TRUE)
deg_comparison_list_deseq <- list.files(path = paste0(input, "DESeq2/"), pattern = "_deg.csv", recursive = TRUE)
deg_comparison_list <- c(deg_comparison_list_deseq, deg_comparison_list_edgeR)

#### ================================
#### 3. Creación de función necesaria posteriormente y establecer umbrales
#### ================================

# Función para obtener IDs originales y nombres originales separados por "/"
map_entrez_to_names <- function(entrez_string, all_table_entrez) {
  if (is.na(entrez_string) || entrez_string == "") return(NA_character_)
  ids <- strsplit(entrez_string, "/")[[1]]
  sub_table <- all_table_entrez[match(ids, all_table_entrez$Entrez), ]
  original_id  <- paste(na.omit(sub_table$gene_id), collapse = "/")
  gene_name    <- paste(na.omit(sub_table$gene_name),   collapse = "/")
  list(original_id = original_id, gene_name = gene_name)
}

# Umbrales de filtrado
padj_thresh <- 0.05
lfc_thresh <- 1

#### ================================
#### 4. Análisis de enriquecimiento GO y KEGG para cada comparación
#### ================================

for (comparison_file in deg_comparison_list) {
  
  comparison <- gsub("_deg_all.csv", "", basename(comparison_file))
  
  # Mensaje de progreso
  cat("Realizando análisis de enriquecimiento de:", comparison, "\n")
  
  # Cargar CSV de DEGs
  if (grepl("SRR", comparison) == TRUE) {
    table <- read.csv2(paste0(input, "edgeR/", comparison_file), header = TRUE, check.names = FALSE)
  } else{
    table <- read.csv2(paste0(input, "DESeq2/", comparison_file), header = TRUE, check.names = FALSE)
  }
  
  # Traducir para la tabla de DEGs los Ensembl IDs a Entrez IDs
  table_entrez <- merge(table, entrez2ensembl, by.x="gene_id", by.y="Ensembl", all.x=TRUE)
  
  # Separar en DEGs totales, up genes y down genes
  gene_set_list <- list()
  deg_genes <- as.character(table_entrez$Entrez)
  up_genes <- as.character(table_entrez$Entrez[table_entrez$log2FoldChange > lfc_thresh])
  down_genes <- as.character(table_entrez$Entrez[table_entrez$log2FoldChange < -lfc_thresh])
  
  gene_set_list[["all"]] <- deg_genes
  gene_set_list[["up"]] <- up_genes
  gene_set_list[["down"]] <- down_genes
  
  # Realizar un análisis GO y un análisis KEGG para cada conjunto de genes
  for (gene_set_name in names(gene_set_list)) {
    
    # Escoger la lista de genes
    gene_entrez <- gene_set_list[[gene_set_name]]
    
    #### ================================
    #### 4.1. Análisis de enriquecimiento GO
    #### ================================
    
    # Ejecutar enrichGO
    go <- enrichGO(gene          = gene_entrez,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   universe      = all_gene_entrez,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   readable      = FALSE)
    
    if (is.null(go) || nrow(as.data.frame(go)) == 0) {
      next
    }
    
    # Elaborar la tabla
    go_df <- as.data.frame(go)
  
    # Añadir columna de símbolos de genes de cada pathway separados por '/'  
    go_df <- go_df %>%
      rowwise() %>%
      mutate(
        gene_name = map_entrez_to_names(geneID, all_table_entrez)$gene_name
      ) %>%
      ungroup()
    
    # Seleccionar y renombrar columnas
    go_table <- go_df %>%
      transmute(
        Category    = ONTOLOGY,
        GO_ID       = ID,
        Description = Description,
        Gene_Ratio  = GeneRatio,
        Bg_Ratio    = BgRatio,
        pvalue      = pvalue,
        padj        = p.adjust,
        gene_ID     = geneID,
        gene_name   = gene_name,
        count       = Count
      )
    
    # Filtrar solo los pathways significativos
    filt_go_table <- go_table[go_table$padj <= padj_thresh, ]
    
    # Crear ruta de salida si no existe
    out_dir <- paste0(output, "GO/", comparison, "/", gene_set_name, "/")
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    # Guardar tabla con todos los pathways y tabla con los pathways significativos
    write.csv2(go_table, paste0(out_dir, "GO_table_", comparison, ".csv"), row.names = FALSE)
    write.csv2(filt_go_table, paste0(out_dir, "signif_GO_table_", comparison, ".csv"), row.names = FALSE)
    
    # Guardar gráfico dotplot
    png(paste0(out_dir, "dotplot_GO - ", comparison, ".png"), width = 550, height = 700)
    print(dotplot(go, showCategory = 20))
    dev.off()
  
    #### ================================
    #### 4.2. Análisis de enriquecimiento KEGG
    #### ================================
    
    # Enriquecimiento KEGG
    kegg <- enrichKEGG(gene          = gene_entrez,
                       organism      = "hsa",
                       universe      = all_gene_entrez,
                       pvalueCutoff  = 1,
                       qvalueCutoff  = 1)
    
    if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0) {
      next
    }
    
    # Elaborar la tabla
    kegg_df <- as.data.frame(kegg)
    
    # Añadir columnas de IDs y símbolos de los genes de cada pathway separados por '/'
    kegg_df <- kegg_df %>%
      rowwise() %>%
      mutate(
        original_id = map_entrez_to_names(geneID, all_table_entrez)$original_id,
        gene_name   = map_entrez_to_names(geneID, all_table_entrez)$gene_name,
        kegg_id     = paste(paste0("hsa:", strsplit(geneID, "/")[[1]]), collapse = "/")
      ) %>%
      ungroup()
    
    # Seleccionar y renombrar columnas
    kegg_table <- kegg_df %>%
      transmute(
        KEGG_ID     = ID,
        Description = Description,
        Gene_Ratio  = GeneRatio,
        BgRatio     = BgRatio,
        pvalue      = pvalue,
        padj        = p.adjust,
        gene_id     = original_id,
        gene_name   = gene_name,
        kegg_id     = kegg_id,
        count       = Count
      )
    
    # Filtrar solo los pathways significativos
    filt_kegg_table <- kegg_table[kegg_table$padj <= padj_thresh, ]
    
    # Crear ruta de salida si no existe
    out_dir <- paste0(output, "KEGG/", comparison, "/", gene_set_name, "/")
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    # Guardar tabla con todos los pathways y tabla con los pathways significativos
    write.csv2(kegg_table, paste0(out_dir, "KEGG_table_", comparison, ".csv"), row.names = FALSE)
    write.csv2(filt_kegg_table, paste0(out_dir, "signif_KEGG_table_", comparison, ".csv"), row.names = FALSE)
    
    # Guardar gráfico dotplot
    png(paste0(out_dir, "dotplot_KEGG - ", comparison, ".png"), width = 550, height = 700)
    print(dotplot(kegg, showCategory = 20))
    dev.off()
  }
}

# Mostrar tiempo de ejecución
t_total <- (proc.time() - t0)/60
cat("Script ejecutado correctamente.\n")
cat("Tiempo de ejecución:", round(t_total["elapsed"], 2), "minutos\n")