#### ================================
#### 1. Configuraciones previas
#### ================================

# Directorios de entrada y salida que debe elegir el usuario:

# Donde se encuentra gene_counts_final.csv
input <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"
# Directorio de salida
output <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/DE/"

# Inicialización del temporizador de tiempo de ejecución
t0 <- proc.time()

# Crear o sobreescribir carpeta de salida
if (dir.exists(output)) {
  unlink(output, recursive = TRUE, force = TRUE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Importar librerías
library(DESeq2)
library(edgeR)

#### ================================
#### 2. Carga de datos
#### ================================

# Cargar tabla de conteos de la cuantificación y guardar las columnas de conteos
counts <- read.csv2(paste0(input, "gene_counts_final.csv"), header = TRUE, check.names = FALSE)
counts_sample <- counts[, grep("SRR", names(counts))]
colnames(counts_sample) <- gsub(" ", "_", colnames(counts_sample))
sample_names <- colnames(counts_sample)

# Cargar tabla de metadatos (muestras y grupos a los que pertenecen)
metadata <- read.csv2(paste0(input, "metadata.csv"), header = TRUE, check.names = FALSE)
groups <- unique(metadata$groups)

# Cargar tablas de fpkm y tpm
fpkm <- read.csv2(paste0(input, "fpkm_table.csv"), header = TRUE, check.names = FALSE)
fpkm_sample <- fpkm[, grep("SRR", names(fpkm))]
tpm <- read.csv2(paste0(input, "tpm_table.csv"), header = TRUE, check.names = FALSE)
tpm_sample <- tpm[, grep("SRR", names(tpm))]

# Cargar counts promedio
counts_group <- read.csv2(paste0(input, "counts_group_table.csv"), header = TRUE, check.names = FALSE)
counts_group_sample <- counts_group[colnames(counts_group) %in% groups]

#### ================================
#### 3. Separar en 2 análisis (DESeq2 y edgeR), inicializar tablas de resultados y umbrales
#### ================================

# Separar las muestras que solo tienen una réplica de las que tienen varias
unique_samples <- metadata$samples[isUnique(metadata$groups)]
final_groups <- unique(metadata$groups[!metadata$samples %in% unique_samples])

# Comparar muestras con muestras y grupos con grupos
unique_samples_comb <- combn(unique_samples, 2, simplify = FALSE)
final_groups_comb <- combn(final_groups, 2, simplify = FALSE)
combinations <- c(final_groups_comb, unique_samples_comb)

# Preparar tabla de conteos por muestra
gene_ids <- counts$gene_id
counts_mat <- as.matrix(counts_sample)
rownames(counts_mat) <- gene_ids

# Inicializar tabla de resultados finales
results_table <- cbind(
  counts[,1, drop = FALSE],
  counts_sample,
  fpkm_sample,
  tpm_sample
)

colnames(results_table) <- c(
  colnames(counts[,1, drop = FALSE]),
  paste0(colnames(counts_sample), "_counts"),
  paste0(colnames(fpkm_sample), "_fpkm"),
  paste0(colnames(tpm_sample), "_tpm")
)

# Umbrales de filtrado
FDR_thresh <- 0.05
padj_thresh <- 0.05
lfc_thresh <- 1
lcpm_thresh <- 1

# Inicializar listas de resultados
all_results <- list()
filt_res_list <- list()
up_res_list <- list()
down_res_list <- list()
all_results_deseq <- list()
filt_res_list_deseq <- list()
up_res_list_deseq <- list()
down_res_list_deseq <- list()

#### ================================
#### 4. Ejecutar edgeR para cada comparación entre muestras y DESeq2 para cada comparación entre grupos
#### ================================

for (pair in combinations) {
  
  # Mensaje de progreso
  comparison <- paste(pair, collapse="_vs_")
  cat("Comparando:", comparison, "\n")
  
  ## Realizar comparaciones entre muestras con edgeR
  if (pair[1] %in% metadata$samples) {
  
    # Seleccionar los counts correspondientes al par de muestras a comparar
    counts_subset <- counts_mat[, c(pair[1], pair[2])]
    
    # Crear el objeto DGEList
    dge <- DGEList(counts=counts_subset, group=pair)
    
    # Calcular los factores de normalización
    dge <- calcNormFactors(dge)
    
    # Estimación manual de la dispersión típica cuando no hay réplicas
    dge$common.dispersion <- 0.4   # Valor aproximado para humanos
    
    # Prueba exacta tipo Fisher
    et <- exactTest(dge, dispersion=dge$common.dispersion)
    
    # Resultados
    res <- topTags(et, n = Inf)
    
    # Reordenar
    res <- res[match(rownames(counts_subset), rownames(res)), ]
    
    # Guardar la tabla de resultados ordenados
    res_table <- res$table
    
    # Renombrar los encabezados para la tabla final
    colnames(res_table) <- c("log2FoldChange", "log2CPM", "pvalue", "FDR")
    res_table_final <- res_table
    colnames(res_table_final) <- paste0(comparison, "_", colnames(res_table_final))
    
    # Añadir a la tabla global
    results_table <- cbind(results_table, res_table_final)
    
    # Añadir anotaciones génicas
    res_table <- cbind(counts[,1, drop = FALSE], data.frame(res_table), counts[,(length(sample_names)+2):ncol(counts)])
    
    # Añadir a la lista de resultados sin filtrar, filtrados, up y down
    all_results[[comparison]] <- res_table
    filt_res_list[[comparison]] <- res_table[!is.na(res_table$FDR) & res_table$FDR <= FDR_thresh & abs(res_table$log2FoldChange) >= lfc_thresh & abs(res_table$log2CPM) >= lcpm_thresh, ]
    up_res_list[[comparison]] <- res_table[!is.na(res_table$FDR) & res_table$FDR <= FDR_thresh & res_table$log2FoldChange >= lfc_thresh & abs(res_table$log2CPM) >= lcpm_thresh, ]
    down_res_list[[comparison]] <- res_table[!is.na(res_table$FDR) & res_table$FDR <= FDR_thresh & res_table$log2FoldChange <= -lfc_thresh & abs(res_table$log2CPM) >= lcpm_thresh, ]
    
    ## Realizar comparaciones entre grupos con DESeq2
  } else {
    
    # Seleccionar los conteos correspondientes al par de muestras a comparar
    samples_subset <- metadata$samples[metadata$groups %in% pair]
    counts_subset <- counts_mat[, c(samples_subset)]
    metadata_subset <- metadata[metadata$groups %in% pair,]
    
    
    # Ejecutar DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts_subset,
                                  colData   = metadata_subset,
                                  design    = ~ groups)
    dds <- DESeq(dds)
    
    res <- results(dds)
    res <- res[, -c(1, 3, 4)]
    
    # Reordenar
    res_table <- res[match(rownames(counts_subset), rownames(res)), ]
    
    # Renombrar los encabezados para la tabla final
    colnames(res_table) <- c("log2FoldChange", "pvalue", "padj")
    res_table_final <- res_table
    colnames(res_table_final) <- paste0(comparison, "_", colnames(res_table_final))
    
    # Añadir a la tabla global
    results_table <- cbind(results_table, res_table_final)
    
    # Añadir anotaciones génicas
    res_table <- cbind(counts[,1, drop = FALSE], data.frame(res_table), counts[,(length(sample_names)+2):ncol(counts)])
    
    # Añadir a la lista de resultados sin filtrar, filtrados, up y down
    all_results_deseq[[comparison]] <- res_table
    filt_res_list_deseq[[comparison]] <- res_table[!is.na(res_table$padj) & res_table$padj <= padj_thresh & abs(res_table$log2FoldChange) >= lfc_thresh, ]
    up_res_list_deseq[[comparison]] <- res_table[!is.na(res_table$padj) & res_table$padj <= padj_thresh & res_table$log2FoldChange >= lfc_thresh, ]
    down_res_list_deseq[[comparison]] <- res_table[!is.na(res_table$padj) & res_table$padj <= padj_thresh & res_table$log2FoldChange <= -lfc_thresh, ]
  }
  
}

# Añadir anotaciones génicas a la tabla global
results_table <- cbind(results_table, counts[,(length(sample_names)+2):ncol(counts)])

# Generar tabla resúmen de genes up y down de edegR
summary_table <- data.frame()

for (name in names(filt_res_list)) {
  res <- filt_res_list[[name]]
  
  total_genes <- nrow(res)
  up_genes    <- sum(res$log2FoldChange >  lfc_thresh)
  down_genes  <- sum(res$log2FoldChange < -lfc_thresh)
  
  summary_table <- rbind(summary_table, data.frame(
    compare   = name,
    all       = total_genes,
    up        = up_genes,
    down      = down_genes,
    threshold = paste0("FDR<", FDR_thresh, " & |log2FC|>", lfc_thresh)
  ))
}

# Generar tabla resúmen de genes up y down de DESeq2 y unirla a la anterior
for (name in names(filt_res_list_deseq)) {
  res <- filt_res_list_deseq[[name]]
  
  total_genes <- nrow(res)
  up_genes    <- sum(res$log2FoldChange >  lfc_thresh)
  down_genes  <- sum(res$log2FoldChange < -lfc_thresh)
  
  summary_table <- rbind(summary_table, data.frame(
    compare   = name,
    all       = total_genes,
    up        = up_genes,
    down      = down_genes,
    threshold = paste0("padj<", padj_thresh, " & |log2FC|>", lfc_thresh)
  ))
}

#### ================================
#### 5. Exportar resultados
#### ================================

# Crear subcarpetas de salida
out_dir_edgeR <- paste0(output, "edgeR/")
if (!dir.exists(out_dir_edgeR)) {
  dir.create(out_dir_edgeR, recursive = TRUE)
}
out_dir_deseq <- paste0(output, "DESeq2/")
if (!dir.exists(out_dir_deseq)) {
  dir.create(out_dir_deseq, recursive = TRUE)
}

# Guardar tabla global
write.csv2(results_table, paste0(output, "Global_table_DE.csv"), row.names = FALSE)

# Guardar tabla resumen 
write.csv2(summary_table, paste0(output, "summary_table_DE.csv"), row.names = FALSE)

# Guardar cada comparación con todos los genes
i=1
for (comparison in names(all_results)){
  dir.create(paste0(out_dir_edgeR, comparison), recursive = TRUE)
  write.csv2(all_results[[i]], paste0(out_dir_edgeR, comparison, "/", comparison, "_all.csv"), row.names = FALSE)
  i<- i + 1
}
i=1
for (comparison in names(all_results_deseq)){
  dir.create(paste0(out_dir_deseq, comparison), recursive = TRUE)
  write.csv2(all_results_deseq[[i]], paste0(out_dir_deseq, comparison, "/", comparison, "_all.csv"), row.names = FALSE)
  i<- i + 1
}

# Guardar cada comparación filtrada (DEGs)
i=1
for (comparison in names(filt_res_list)){
  write.csv2(filt_res_list[[i]], paste0(out_dir_edgeR, comparison, "/", comparison, "_deg_all.csv"), row.names = FALSE)
  i<- i + 1
}
i=1
for (comparison in names(filt_res_list_deseq)){
  write.csv2(filt_res_list_deseq[[i]], paste0(out_dir_deseq, comparison, "/", comparison, "_deg_all.csv"), row.names = FALSE)
  i<- i + 1
}

# Guardar DEGs up
i=1
for (comparison in names(up_res_list)){
  write.csv2(up_res_list[[i]], paste0(out_dir_edgeR, comparison, "/", comparison, "_deg_up.csv"), row.names = FALSE)
  i<- i + 1
}
i=1
for (comparison in names(up_res_list_deseq)){
  write.csv2(up_res_list_deseq[[i]], paste0(out_dir_deseq, comparison, "/", comparison, "_deg_up.csv"), row.names = FALSE)
  i<- i + 1
}

# Guardar DEGs down
i=1
for (comparison in names(down_res_list)){
  write.csv2(down_res_list[[i]], paste0(out_dir_edgeR, comparison, "/", comparison, "_deg_down.csv"), row.names = FALSE)
  i<- i + 1
}
i=1
for (comparison in names(down_res_list_deseq)){
  write.csv2(down_res_list_deseq[[i]], paste0(out_dir_deseq, comparison, "/", comparison, "_deg_down.csv"), row.names = FALSE)
  i<- i + 1
}

# Mostrar tiempo de ejecución
t_total <- (proc.time() - t0)/60
cat("Script ejecutado correctamente.\n")
cat("Tiempo de ejecución:", round(t_total["elapsed"], 2), "minutos\n")