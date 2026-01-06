#### ================================
#### 1. Configuraciones previas
#### ================================

# Directorios de entrada y salida que debe elegir el usuario:

# Directorio con los resultados del análisis de expresión diferencial
input <- "C:/Users/Alexinander/Desktop/TFM/DE/"
gsea_home <- "C:/Program Files/GSEA_4.4.0/"
# Directorio de salida
output <- "C:/enrichment/"
# Directorio con la tabla de correspondencias entre Ensembl IDs y Entrez IDs
corresp_table_dir <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"

# Inicialización del temporizador de tiempo de ejecución
t0 <- proc.time()

# Importar librerías
library(msigdbr)
library(KEGGREST)
library(dplyr)
library(org.Hs.eg.db)
library(tools)

# Crear directorio para los archivos GMT (pathways) y RNK (listas de genes)
gmt_rnk_dir <- paste0(output, "gmt_rnk/")
if (!dir.exists(gmt_rnk_dir)) {dir.create(gmt_rnk_dir, recursive = TRUE)}

#### ================================
#### 2. Carga de datos
#### ================================

# Coger las tablas de cada comparación con los resultados de expresión diferencial para todo el universo de genes
comparison_list_DESeq2 <- list.files(path = paste0(input, "DESeq2/"), pattern = "_all.csv", recursive = TRUE, full.names  = TRUE)
comparison_list_edgeR <- list.files(path = paste0(input, "edgeR/"), pattern = "_all.csv", recursive = TRUE, full.names  = TRUE)
comparison_list <- c(comparison_list_DESeq2, comparison_list_edgeR)

# Usar gene2ensembl para hacer la correspondecia entre Entrez IDs y Ensembl IDs
entrez2ensembl <- read.table(file = paste0(corresp_table_dir,"entrez2ensembl.txt"), sep = "\t", header = TRUE)

#### ================================
#### 3. Obtención de pathways GO que se emplearán en GSEA
#### ================================

# Obtener las vías de GO desde MSigDB
go_df <- msigdbr(species = 'Homo sapiens', collection = 'C5')

# Excluir vias HPO para dejar solo las de GO
go_df <- go_df[grepl("^GO:", go_df$gs_subcollection), ]

# Crear tabla con el código y la descripción de cada pathway
go_df_codes <- go_df[!duplicated(go_df$gs_name), ]
go_df_codes <- go_df_codes[, c("gs_name", "gs_exact_source", "gs_description")]

# Construir tabla agrupada para GMT
gmt_go <- go_df %>%
  group_by(gs_name, gs_description) %>%
  summarise(genes = list(unique(ncbi_gene)), .groups = "drop")

# Exportar tabla a GMT
con <- file(paste0(gmt_rnk_dir, "MSigDB_GO_all.gmt"), "w")
apply(gmt_go, 1, function(row) {
  line <- c(
    row[["gs_name"]],
    row[["gs_description"]],
    unlist(row[["genes"]])
  )
  writeLines(paste(line, collapse = "\t"), con)
})
close(con)

#### ================================
#### 4. Obtención de pathways KEGG que se emplearán en GSEA
#### ================================

# Descargar KEGG_LEGACY
m_kegg_legacy <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")

# Descargar KEGG_MEDICUS
m_kegg_medicus <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_MEDICUS")

# Unir ambas colecciones
kegg_df <- bind_rows(m_kegg_legacy, m_kegg_medicus)

# Crear tabla con el código y la descripción de cada pathway
kegg_df_codes <- kegg_df[!duplicated(kegg_df$gs_name), ]
kegg_df_codes <- kegg_df_codes[, c("gs_name", "gs_exact_source", "gs_description")]

# Construir tabla agrupada para GMT
gmt_kegg_all <- kegg_df %>%
  group_by(gs_name, gs_description) %>%
  summarise(genes = list(unique(ncbi_gene)), .groups = "drop")

# Exportar tabla a GMT
con <- file(paste0(gmt_rnk_dir, "MSigDB_KEGG_all.gmt"), "w")
apply(gmt_kegg_all, 1, function(row) {
  line <- c(
    row[["gs_name"]],
    row[["gs_description"]],
    unlist(row[["genes"]])
  )
  writeLines(paste(line, collapse = "\t"), con)
})
close(con)

#### ================================
#### 4. Ejecutar GSEA para cada comparación
#### ================================

# Cambiar directorio de trabajo al que contiene el programa de GSEA
oldwd <- getwd()
setwd(gsea_home)

# Iterar sobre cada comparación
for (comparison_file in comparison_list) {
  
  # Nombre de la comparación
  comparison <- gsub("_all.csv", "", basename(comparison_file))
  groups <- strsplit(comparison, "_vs_")[[1]]
  
  # Mensaje de progreso
  cat("Realizando GSEA de:", comparison, "\n")
  
  # Cargar tabla de esa comparación con todos los genes del universo
  table <- read.csv2(comparison_file, header = TRUE, check.names = FALSE)

  # Filtrar todos los genes cuyo pvalor o log2fc sea Na y sustituir pvalue=0 por el mínimo
  table <- table[!is.na(table$pvalue) & !is.na(table$log2FoldChange), ]
  table$pvalue[table$pvalue == 0] <- .Machine$double.xmin
  
  # Usar gene2ensembl para traducir Ensembl IDs a Entrez IDs
  table_entrez <- merge(table, entrez2ensembl, by.x="gene_id", by.y="Ensembl", all.x=TRUE)
  entrez_ids <- as.character(table_entrez$Entrez)
  
  # Crear tabla RNK de genes ordenados por score = −log10(p value) x (log2FC)
  p_values <- table_entrez$pvalue
  log2fc <- table_entrez$log2FoldChange
  gene_list <- -log10(p_values)*log2fc
  names(gene_list) <- entrez_ids
  gene_list <- sort(gene_list, decreasing = TRUE)
  rnk <- cbind(names(gene_list), gene_list)
  colnames(rnk) <- NULL
  write.table(rnk, paste0(gmt_rnk_dir, comparison, "_ranked_list.rnk"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Crear directorios de salida
  go_out_dir <- paste0(output, "GO/", comparison, "/", "GSEA/")
  if (dir.exists(go_out_dir)) {
    unlink(go_out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(go_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  kegg_out_dir <- paste0(output, "KEGG/", comparison, "/", "GSEA/")
  if (dir.exists(kegg_out_dir)) {
    unlink(kegg_out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(kegg_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Comando para ejecutar GSEA Broad con pathways GO
  comando_gsea_go <- paste(
    "java",
    "--module-path=modules",
    "-Xmx4g",
    "-Djava.awt.headless=true",
    "-Djava.util.logging.config.file=logging.properties",
    "@gsea.args",
    "--module=org.gsea_msigdb.gsea/xapps.gsea.CLI",
    "GSEAPreranked",
    "-rpt_label", shQuote(comparison),
    "-rnk",       shQuote(paste0(gmt_rnk_dir, comparison, "_ranked_list.rnk")),
    "-gmx",       shQuote(paste0(gmt_rnk_dir, "MSigDB_GO_all.gmt")),
    "-out",       shQuote(go_out_dir)
  )
  
  # Comando para ejecutar GSEA Broad con pathways KEGG
  comando_gsea_kegg <- paste(
    "java",
    "--module-path=modules",
    "-Xmx4g",
    "-Djava.awt.headless=true",
    "-Djava.util.logging.config.file=logging.properties",
    "@gsea.args",
    "--module=org.gsea_msigdb.gsea/xapps.gsea.CLI",
    "GSEAPreranked",
    "-rpt_label", shQuote(comparison),
    "-rnk",       shQuote(paste0(gmt_rnk_dir, comparison, "_ranked_list.rnk")),
    "-gmx",       shQuote(paste0(gmt_rnk_dir, "MSigDB_KEGG_all.gmt")),
    "-out",       shQuote(kegg_out_dir)
  )
  
  # Ejecutar los comandos
  system(comando_gsea_go)
  system(comando_gsea_kegg)
  
  #### ================================
  #### 5. Creación de tablas de resultados para análisis GO
  #### ================================
  
  # Localizar carpeta de salida de GSEA de análisis GO
  all_dirs <- list.dirs(go_out_dir, recursive = TRUE, full.names = TRUE)
  good_dirs <- vapply(
    all_dirs,
    function(d) {
      any(grepl("^gsea_report_for_.*\\.tsv$", list.files(d)))
    },
    logical(1)
  )
  gsea_broad_dir <- all_dirs[good_dirs][1]
  
  # Listar todos los archivos de esa carpeta
  all_files <- list.files(gsea_broad_dir, recursive = TRUE, full.names = TRUE)
  
  #### ================================
  #### 5.1. Tabla ranked_gene_list de GO (igual en KEGG)
  #### ================================
  
  # Extraer los IDs de genes y su score e incorporar a la tabla
  df_rank <- data.frame(
    GENE_ID    = names(gene_list),
    SCORE      = as.numeric(gene_list),
    stringsAsFactors = FALSE
  )
  
  # Extraer mediante org.Hs.eg.db los símbolos y los títulos de cada gen
  ann <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = df_rank$GENE_ID,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  )
  
  # Incoporar símbolo y título a la tabla
  df_rank <- merge(df_rank, ann, by.x = "GENE_ID", by.y = "ENTREZID", all.x = TRUE)
  df_rank <- merge(df_rank, table_entrez, by.x = "GENE_ID", by.y = "Entrez", all.x = TRUE)
  df_rank <- df_rank[!duplicated(df_rank$GENE_ID), ]
  
  # Crear de la tabla final ranked_gene_list tanto para GSEA GO como para GSEA KEGG
  ranked_gene_list <- data.frame(
    NAME        = df_rank$gene_id,
    GENE_SYMBOL = df_rank$SYMBOL,
    GENE_TITLE  = df_rank$GENENAME,
    SCORE       = df_rank$SCORE
  )
  
  # Ordenar por score
  ranked_gene_list <- ranked_gene_list[order(ranked_gene_list$SCORE, decreasing = TRUE), ]
  
  # Exportar tabla
  write.csv2(ranked_gene_list, paste0(go_out_dir, "ranked_gene_list_", comparison, ".csv"), row.names = FALSE)
  write.csv2(ranked_gene_list, paste0(kegg_out_dir, "ranked_gene_list_", comparison, ".csv"), row.names = FALSE)
  
  #### ================================
  #### 5.2. Tabla gene_set_sizes de GO
  #### ================================
  
  # Cargar tabla y eliminar columna no válida
  gene_set_sizes <- read.table(paste0(gsea_broad_dir, "/gene_set_sizes.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gene_set_sizes <- gene_set_sizes[, !(names(gene_set_sizes) %in% "X")]
  
  # Fusionar con la tabla que tiene el código y la descripción de cada ruta y mantener el orden
  gene_set_sizes$.__ord__ <- seq_len(nrow(gene_set_sizes))
  pathway_codes <- merge(gene_set_sizes, go_df_codes, by.x = "NAME", by.y = "gs_name", all.x = TRUE)
  pathway_codes <- pathway_codes[order(pathway_codes$.__ord__), ]
  gene_set_sizes$.__ord__ <- NULL
  
  # Añadir el código entre paréntesis al lado del nombre
  gene_set_sizes$NAME <- gsub("GO.*?_", "", paste0(gene_set_sizes$NAME, " (", pathway_codes$gs_exact_source, ")"))
  colnames(gene_set_sizes) <- gsub("\\.", " ", colnames(gene_set_sizes))
  write.csv2(gene_set_sizes, paste0(go_out_dir, "gene_set_sizes.csv"), row.names = FALSE)
  
  #### ================================
  #### 5.3. Tablas gsea_report de GO
  #### ================================
  
  # Localizar tablas
  gsea_report_find <- all_files[
    grepl("gsea_report_for", basename(all_files), ignore.case = TRUE) &
      grepl("\\.tsv$", all_files, ignore.case = TRUE)
  ]
  
  # Cargar tablas
  gsea_report_g1 <- read.table(gsea_report_find[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gsea_report_g2 <- read.table(gsea_report_find[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Eliminar columnas no válidas
  gsea_report_g1 <- gsea_report_g1[, !(names(gsea_report_g1) %in% "X")]
  gsea_report_g2 <- gsea_report_g2[, !(names(gsea_report_g2) %in% "X")]
  
  # Mantener el orden para fusionar con la tabla que contiene código y descripción
  gsea_report_g1$.__ord__ <- seq_len(nrow(gsea_report_g1))
  gsea_report_g2$.__ord__ <- seq_len(nrow(gsea_report_g2))
  
  # Añadir código y decsripción a la tabla del grupo 1
  pathway_codes <- merge(gsea_report_g1, go_df_codes, by.x = "NAME", by.y = "gs_name", all.x = TRUE)
  pathway_codes <- pathway_codes[order(pathway_codes$.__ord__), ]
  gsea_report_g1$.__ord__ <- NULL
  gsea_report_g1$NAME <- gsub("GO.*?_", "", paste0(gsea_report_g1$NAME, " (", pathway_codes$gs_exact_source, ")"))
  gsea_report_g1$GS.DETAILS <- pathway_codes$gs_description
  
  # Añadir código y decsripción a la tabla del grupo 2
  pathway_codes <- merge(gsea_report_g2, go_df_codes, by.x = "NAME", by.y = "gs_name", all.x = TRUE)
  pathway_codes <- pathway_codes[order(pathway_codes$.__ord__), ]
  gsea_report_g2$.__ord__ <- NULL
  gsea_report_g2$NAME <- gsub("GO.*?_", "", paste0(gsea_report_g2$NAME, " (", pathway_codes$gs_exact_source, ")"))
  gsea_report_g2$GS.DETAILS <- pathway_codes$gs_description
  
  # Añadir link de cada pathway a las tablas
  gsea_report_g1$GS.br..follow.link.to.MSigDB = paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/", gsea_report_g1$GS.br..follow.link.to.MSigDB, ".html")
  gsea_report_g2$GS.br..follow.link.to.MSigDB = paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/", gsea_report_g2$GS.br..follow.link.to.MSigDB, ".html")
  
  # Crear directorios correspondientes y exportar
  m1_dir <- paste0(go_out_dir, groups[1], "/")
  m2_dir <- paste0(go_out_dir, groups[2], "/")
  m1_pathway_dir <- paste0(m1_dir, "pathways/")
  m2_pathway_dir <- paste0(m2_dir, "pathways/")
  if (!dir.exists(m1_dir)) {dir.create(m1_dir, recursive = TRUE)}
  if (!dir.exists(m2_dir)) {dir.create(m2_dir, recursive = TRUE)}
  if (!dir.exists(m1_pathway_dir)) {dir.create(m1_pathway_dir, recursive = TRUE)}
  if (!dir.exists(m2_pathway_dir)) {dir.create(m2_pathway_dir, recursive = TRUE)}
  
  colnames(gsea_report_g1) <- gsub("\\.", " ", colnames(gsea_report_g1))
  colnames(gsea_report_g2) <- gsub("\\.", " ", colnames(gsea_report_g2))
  
  write.csv2(gsea_report_g1, paste0(m1_dir, "gsea_report_", groups[1], ".csv"), row.names = FALSE)
  write.csv2(gsea_report_g2, paste0(m2_dir, "gsea_report_", groups[2], ".csv"), row.names = FALSE)
  
  #### ================================
  #### 5.4. Tablas de los 20 mejores pathways GO de cada grupo
  #### ================================
  
  # Localizar tablas
  pathway_find <- all_files[
    grepl("GO", basename(all_files), ignore.case = TRUE) &
      grepl("\\.tsv$", all_files, ignore.case = TRUE)
  ]
  
  # Cargar tablas
  for (pathway in pathway_find){
    pathway_table <- read.table(paste0("\\\\?\\", normalizePath(pathway, winslash = "\\", mustWork = FALSE)), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Añadir código entre paréntesis a cada nombre de pathway
    pathway_name <- file_path_sans_ext(basename(pathway))
    pathway_code <- go_df_codes$gs_exact_source[go_df_codes$gs_name == pathway_name]
    pathway_true_name <- gsub("GO.*?_", "", paste0(pathway_name, " (", pathway_code, ")"))
    
    # Crear directorios de salida y nombres de los archivos
    pathway_save_name_m1 <- paste0(m1_pathway_dir, gsub(":", "_", pathway_true_name), ".csv")
    pathway_save_name_m2 <- paste0(m2_pathway_dir, gsub(":", "_", pathway_true_name), ".csv")
    
    # Recortar nombres si los directorios fueran demasiado largos
    if (nchar(pathway_save_name_m1) > 252){
      cut_len <- nchar(pathway_save_name_m1) - 252
      cut_pathway_true_name <- gsub("GO.*?_", "", paste0(substr(pathway_name, 1, nchar(pathway_name)-cut_len), "(etc) (", pathway_code, ")"))
      pathway_save_name_m1 <- paste0(m1_pathway_dir, gsub(":", "_", cut_pathway_true_name), ".csv")
    }
    if (nchar(pathway_save_name_m2) > 252){
      cut_len <- nchar(pathway_save_name_m2) - 252
      cut_pathway_true_name <- gsub("GO.*?_", "", paste0(substr(pathway_name, 1, nchar(pathway_name)-cut_len), "(etc) (", pathway_code, ")"))
      pathway_save_name_m2 <- paste0(m2_pathway_dir, gsub(":", "_", cut_pathway_true_name), ".csv")
    }
    
    # Fusionar con otras tablas para añadir información
    pathway_table_temp <- merge(pathway_table, go_df, by.x = "SYMBOL", by.y = "ncbi_gene", all.x = TRUE)
    pathway_table_temp <- merge(pathway_table_temp, ann[, c("ENTREZID", "GENENAME")], by.x = "SYMBOL", by.y = "ENTREZID", all.x = TRUE)
    pathway_table_temp <- pathway_table_temp[!duplicated(pathway_table_temp$SYMBOL), ]
    
    # Crear cada tabla final
    pathway_table <- data.frame(
      NAME                     = pathway_table_temp$NAME,
      PROBE                    = pathway_table_temp$ensembl_gene,
      GENE_SYMBOL              = pathway_table_temp$gene_symbol,
      GENE_TITLE               = pathway_table_temp$GENENAME,
      RANK_IN_GENE_LIST        = pathway_table_temp$RANK.IN.GENE.LIST,
      RANK_METRIC_SCORE        = pathway_table_temp$RANK.METRIC.SCORE,
      CORE_ENRICHMENT          = pathway_table_temp$CORE.ENRICHMENT
    )
    colnames(pathway_table) <- gsub("_", " ", colnames(pathway_table))
    
    # Ordenar tal y como las ordenó GSEA
    pathway_table$.__ord__ <- as.numeric(gsub("row_", "", pathway_table$NAME))
    pathway_table <- pathway_table[order(pathway_table$.__ord__), ]
    pathway_table$.__ord__ <- NULL
    
    # Exportar
    if (pathway_true_name %in% gsea_report_g1$NAME) {
      write.csv2(pathway_table, pathway_save_name_m1, row.names = FALSE)
    } else if (pathway_true_name %in% gsea_report_g2$NAME) {
      write.csv2(pathway_table, pathway_save_name_m2, row.names = FALSE)
    }
  }
  
  #### ================================
  #### 6. Creación de tablas de resultados para análisis KEGG
  #### ================================
  
  # Localizar carpeta de salida de GSEA de análisis KEGG
  all_dirs <- list.dirs(kegg_out_dir, recursive = TRUE, full.names = TRUE)
  good_dirs <- vapply(
    all_dirs,
    function(d) {
      any(grepl("^gsea_report_for_.*\\.tsv$", list.files(d)))
    },
    logical(1)
  )
  gsea_broad_dir <- all_dirs[good_dirs][1]
  
  # Listar todos los archivos de esa carpeta
  all_files <- list.files(gsea_broad_dir, recursive = TRUE, full.names = TRUE)
  
  #### ================================
  #### 6.1. Tabla gene_set_sizes de KEGG
  #### ================================
  
  # Cargar tabla y eliminar columna no válida
  gene_set_sizes <- read.table(paste0(gsea_broad_dir, "/gene_set_sizes.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gene_set_sizes <- gene_set_sizes[, !(names(gene_set_sizes) %in% "X")]
  
  # Fusionar con la tabla que tiene el código y la descripción de cada ruta y mantener el orden
  gene_set_sizes$.__ord__ <- seq_len(nrow(gene_set_sizes))
  pathway_codes<- merge(gene_set_sizes, kegg_df_codes, by.x = "NAME", by.y = "gs_name", all.x = TRUE)
  pathway_codes <- pathway_codes[order(pathway_codes$.__ord__), ]
  gene_set_sizes$.__ord__ <- NULL
  
  # Añadir el código entre paréntesis al lado del nombre
  gene_set_sizes$NAME <- gsub("KEGG.*?_", "", paste0(gene_set_sizes$NAME, " (", pathway_codes$gs_exact_source, ")"))
  colnames(gene_set_sizes) <- gsub("\\.", " ", colnames(gene_set_sizes))
  write.csv2(gene_set_sizes, paste0(kegg_out_dir, "gene_set_sizes.csv"), row.names = FALSE)
  
  #### ================================
  #### 6.2. Tablas gsea_report de KEGG
  #### ================================
  
  # Localizar tablas
  gsea_report_find <- all_files[
    grepl("gsea_report_for", basename(all_files), ignore.case = TRUE) &
      grepl("\\.tsv$", all_files, ignore.case = TRUE)
  ]
  
  # Cargar tablas
  gsea_report_g1 <- read.table(gsea_report_find[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gsea_report_g2 <- read.table(gsea_report_find[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Eliminar columnas no válidas
  gsea_report_g1 <- gsea_report_g1[, !(names(gsea_report_g1) %in% "X")]
  gsea_report_g2 <- gsea_report_g2[, !(names(gsea_report_g2) %in% "X")]

  # Mantener el orden para fusionar con la tabla que contiene código y descripción
  gsea_report_g1$.__ord__ <- seq_len(nrow(gsea_report_g1))
  gsea_report_g2$.__ord__ <- seq_len(nrow(gsea_report_g2))
  
  # Añadir código y descripción a la tabla del grupo 1
  pathway_codes <- merge(gsea_report_g1, kegg_df_codes, by.x = "NAME", by.y = "gs_name", all.x = TRUE)
  pathway_codes <- pathway_codes[order(pathway_codes$.__ord__), ]
  gsea_report_g1$.__ord__ <- NULL
  gsea_report_g1$NAME <- gsub("KEGG.*?_", "", paste0(gsea_report_g1$NAME, " (", pathway_codes$gs_exact_source, ")"))
  gsea_report_g1$GS.DETAILS <- pathway_codes$gs_description
  
  # Añadir código y descripción a la tabla del grupo 2
  pathway_codes <- merge(gsea_report_g2, kegg_df_codes, by.x = "NAME", by.y = "gs_name", all.x = TRUE)
  pathway_codes <- pathway_codes[order(pathway_codes$.__ord__), ]
  gsea_report_g2$.__ord__ <- NULL
  gsea_report_g2$NAME <- gsub("KEGG.*?_", "", paste0(gsea_report_g2$NAME, " (", pathway_codes$gs_exact_source, ")"))
  gsea_report_g2$GS.DETAILS <- pathway_codes$gs_description
  
  # Añadir link de cada pathway a las tablas
  gsea_report_g1$GS.br..follow.link.to.MSigDB = paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/", gsea_report_g1$GS.br..follow.link.to.MSigDB, ".html")
  gsea_report_g2$GS.br..follow.link.to.MSigDB = paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/", gsea_report_g2$GS.br..follow.link.to.MSigDB, ".html")
  
  # Crear directorios correspondientes y exportar
  m1_dir <- paste0(kegg_out_dir, groups[1], "/")
  m2_dir <- paste0(kegg_out_dir, groups[2], "/")
  m1_pathway_dir <- paste0(m1_dir, "pathways/")
  m2_pathway_dir <- paste0(m2_dir, "pathways/")
  if (!dir.exists(m1_dir)) {dir.create(m1_dir, recursive = TRUE)}
  if (!dir.exists(m2_dir)) {dir.create(m2_dir, recursive = TRUE)}
  if (!dir.exists(m1_pathway_dir)) {dir.create(m1_pathway_dir, recursive = TRUE)}
  if (!dir.exists(m2_pathway_dir)) {dir.create(m2_pathway_dir, recursive = TRUE)}
  
  colnames(gsea_report_g1) <- gsub("\\.", " ", colnames(gsea_report_g1))
  colnames(gsea_report_g2) <- gsub("\\.", " ", colnames(gsea_report_g2))
  
  write.csv2(gsea_report_g1, paste0(m1_dir, "gsea_report_", groups[1], ".csv"), row.names = FALSE)
  write.csv2(gsea_report_g2, paste0(m2_dir, "gsea_report_", groups[2], ".csv"), row.names = FALSE)
  
  #### ================================
  #### 6.3. Tablas de los 20 mejores pathways KEGG de cada grupo
  #### ================================
  
  # Localizar tablas
  pathway_find <- all_files[
    grepl("KEGG", basename(all_files), ignore.case = TRUE) &
      grepl("\\.tsv$", all_files, ignore.case = TRUE)
  ]
  
  # Cargar tablas
  for (pathway in pathway_find){
    pathway_table <- read.table(paste0("\\\\?\\", normalizePath(pathway, winslash = "\\", mustWork = FALSE)), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Añadir código entre paréntesis a cada nombre de pathway
    pathway_name <- file_path_sans_ext(basename(pathway))
    pathway_code <- kegg_df_codes$gs_exact_source[kegg_df_codes$gs_name == pathway_name]
    pathway_true_name <- gsub("KEGG.*?_", "", paste0(pathway_name, " (", pathway_code, ")"))
    
    # Crear directorios de salida y nombres de los archivos
    pathway_save_name_m1 <- paste0(m1_pathway_dir, gsub(":", "_", pathway_true_name), ".csv")
    pathway_save_name_m2 <- paste0(m2_pathway_dir, gsub(":", "_", pathway_true_name), ".csv")
    
    # Recortar nombres si los directorios fueran demasiado largos
    if (nchar(pathway_save_name_m1) > 252){
      cut_len <- nchar(pathway_save_name_m1) - 252
      cut_pathway_true_name <- gsub("KEGG.*?_", "", paste0(substr(pathway_name, 1, nchar(pathway_name)-cut_len), "(etc) (", pathway_code, ")"))
      pathway_save_name_m1 <- paste0(m1_pathway_dir, gsub(":", "_", cut_pathway_true_name), ".csv")
    }
    if (nchar(pathway_save_name_m2) > 252){
      cut_len <- nchar(pathway_save_name_m2) - 252
      cut_pathway_true_name <- gsub("KEGG.*?_", "", paste0(substr(pathway_name, 1, nchar(pathway_name)-cut_len), "(etc) (", pathway_code, ")"))
      pathway_save_name_m2 <- paste0(m2_pathway_dir, gsub(":", "_", cut_pathway_true_name), ".csv")
    }
    
    # Fusionar con otras tablas para añadir información
    pathway_table_temp <- merge(pathway_table, kegg_df, by.x = "SYMBOL", by.y = "ncbi_gene", all.x = TRUE)
    pathway_table_temp <- merge(pathway_table_temp, ann[, c("ENTREZID", "GENENAME")], by.x = "SYMBOL", by.y = "ENTREZID", all.x = TRUE)
    pathway_table_temp <- pathway_table_temp[!duplicated(pathway_table_temp$SYMBOL), ]
    
    # Crear cada tabla final
    pathway_table <- data.frame(
      NAME                     = pathway_table_temp$NAME,
      PROBE                    = pathway_table_temp$ensembl_gene,
      GENE_SYMBOL              = pathway_table_temp$gene_symbol,
      GENE_TITLE               = pathway_table_temp$GENENAME,
      RANK_IN_GENE_LIST        = pathway_table_temp$RANK.IN.GENE.LIST,
      RANK_METRIC_SCORE        = pathway_table_temp$RANK.METRIC.SCORE,
      CORE_ENRICHMENT          = pathway_table_temp$CORE.ENRICHMENT
    )
    colnames(pathway_table) <- gsub("_", " ", colnames(pathway_table))

    # Ordenar tal y como las ordenó GSEA
    pathway_table$.__ord__ <- as.numeric(gsub("row_", "", pathway_table$NAME))
    pathway_table <- pathway_table[order(pathway_table$.__ord__), ]
    pathway_table$.__ord__ <- NULL
    
    # Exportar
    if (pathway_true_name %in% gsea_report_g1$NAME) {
      write.csv2(pathway_table, pathway_save_name_m1, row.names = FALSE)
    } else if (pathway_true_name %in% gsea_report_g2$NAME) {
      write.csv2(pathway_table, pathway_save_name_m2, row.names = FALSE)
    }
  }
}

#### ================================
#### 7. Ajustes finales
#### ================================

# Cambiar directorio de trabajo al original
setwd(oldwd)

# Borrar archivos y directorios no deseados
dirs <- list.dirs(output, recursive = TRUE, full.names = TRUE)
dirs_delete <- dirs[grepl("GseaPreranked", basename(dirs), ignore.case = TRUE)]
unlink(dirs_delete, recursive = TRUE, force = TRUE)

# Mostrar tiempo de ejecución
t_total <- (proc.time() - t0)/60
cat("Script ejecutado correctamente.\n")
cat("Tiempo de ejecución:", round(t_total["elapsed"], 2), "minutos\n")