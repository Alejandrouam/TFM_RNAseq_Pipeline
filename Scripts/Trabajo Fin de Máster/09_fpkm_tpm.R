#### ================================
#### 1. Configuraciones previas
#### ================================

# Directorios de entrada y salida que debe elegir el usuario:

# Donde se encuentra gene_counts_final.csv
input <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"
# Directorio de salida
output <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"

# Inicialización del temporizador de tiempo de ejecución
t0 <- proc.time()

# Importar librerías
library(dplyr)
library(tidyr)
library(tibble)

#### ================================
#### 2. Carga de datos
#### ================================

# Cargar tabla de conteos de la cuantificación
counts <- read.csv2(paste0(input, "gene_counts_final.csv"), header = TRUE, check.names = FALSE)

# Coger columnas correspondientes a los conteos por muestra
counts_sample <- counts[, grep("SRR", names(counts))]
colnames(counts_sample) <- gsub(" ", "_", colnames(counts_sample))
sample_names <- colnames(counts_sample)

#### ================================
#### 3. Cálculo de FPKM
#### ================================

# Normalización por muestra
rpkm_sc_factor <- colSums(counts_sample) / 1e6
rpm <- sweep(counts_sample, 2, rpkm_sc_factor, "/")

# Normalización por longitud del gen
gene_lengths <- as.numeric(counts$gene_length) / 1000
fpkm <- sweep(rpm, 1, gene_lengths, "/")

#### ================================
#### 4. Cálculo de TPM
#### ================================

# Normalización por longitud del gen
rpk <- sweep(counts_sample, 1, gene_lengths, "/")

# Normalización por muestra normalizada
tpm_sc_factor <- colSums(rpk) / 1e6
tpm <- sweep(rpk, 2, tpm_sc_factor, "/")

#### ================================
#### 5. Exportar resultados de FPKM y TPM por muestra
#### ================================

# Exportar a csv
fpkm_table <- cbind(counts[,1, drop = FALSE], fpkm, counts[,(length(sample_names)+2):ncol(counts)])
write.csv2(fpkm_table, paste0(output, "fpkm_table.csv"), row.names = FALSE)

tpm_table <- cbind(counts[,1, drop = FALSE], tpm, counts[,(length(sample_names)+2):ncol(counts)])
write.csv2(tpm_table, paste0(output, "tpm_table.csv"), row.names = FALSE)

#### ================================
#### 6. Cálculo de FPKM y TPM por grupo y conteos por grupo
#### ================================

# Asignar un grupo a cada muestra
groups <- c()

for (sample_name in sample_names) {
  group <- gsub("_.*", "", sample_name)
  groups <- c(groups, group)
}

sample_groups <- data.frame(samples = sample_names, groups = groups)
sample_groups$samples <- gsub(" ", "_", sample_groups$samples)
sample_groups$groups <- gsub(" ", "_", sample_groups$groups)
group_names <- unique(sample_groups$groups)

# Guardar tabla de muestras y grupos
write.csv2(sample_groups, paste0(output, "metadata.csv"), row.names = FALSE)

# Cálculo de FPKM por grupo
fpkm_group <- fpkm
for (group in group_names) {
  cols <- sample_groups$samples[sample_groups$groups == group]
  if (length(cols) > 1) {
    fpkm_group[[group]] <- rowMeans(fpkm_group[, all_of(cols)], na.rm = TRUE)
  } else {
    fpkm_group[[group]] <- fpkm_group[, cols]
  }
}
fpkm_group <- fpkm_group[, (length(sample_names) + 1):ncol(fpkm_group)]

# Cálculo de TPM por grupo
tpm_group <- tpm
for (group in group_names) {
  cols <- sample_groups$samples[sample_groups$groups == group]
  if (length(cols) > 1) {
    tpm_group[[group]] <- rowMeans(tpm_group[, all_of(cols)], na.rm = TRUE)
  } else {
    tpm_group[[group]] <- tpm_group[, cols]
  }
}
tpm_group <- tpm_group[, (length(sample_names) + 1):ncol(tpm_group)]

# Cálculo de los conteos promedio de cada grupo
counts_group <- counts_sample
for (group in group_names) {
  cols <- sample_groups$samples[sample_groups$groups == group]
  if (length(cols) > 1) {
    counts_group[[group]] <- rowMeans(counts_group[, all_of(cols)], na.rm = TRUE)
  } else {
    counts_group[[group]] <- counts_group[[cols]]
  }
}
counts_group <- counts_group[, (length(sample_names) + 1):ncol(counts_group)]

#### ================================
#### 7. Añadir anotaciones génicas y exportar resultados
#### ================================

fpkm_group_table <- cbind(counts[,1, drop = FALSE], fpkm_group, counts[,(length(sample_names)+2):ncol(counts)])
write.csv2(fpkm_group_table, paste0(output, "fpkm_group_table.csv"), row.names = FALSE)

tpm_group_table <- cbind(counts[,1, drop = FALSE], tpm_group, counts[,(length(sample_names)+2):ncol(counts)])
write.csv2(tpm_group_table, paste0(output, "tpm_group_table.csv"), row.names = FALSE)

counts_group_table <- cbind(counts[,1, drop = FALSE], counts_group, counts[,(length(sample_names)+2):ncol(counts)])
write.csv2(counts_group_table, paste0(output, "counts_group_table.csv"), row.names = FALSE)

# Mostrar tiempo de ejecución
t_total <- (proc.time() - t0)/60
cat("Script ejecutado correctamente.\n")
cat("Tiempo de ejecución:", round(t_total["elapsed"], 2), "minutos\n")
