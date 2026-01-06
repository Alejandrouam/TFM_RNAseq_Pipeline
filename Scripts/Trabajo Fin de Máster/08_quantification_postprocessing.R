#### ================================
#### 1. Configuraciones previas
#### ================================

# Directorios de entrada y salida que debe elegir el usuario:

# Donde se encuentra gene_counts.csv generado durante la cuantificación
input <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"
# Directorio de salida
output <- "//wsl.localhost/Ubuntu/home/alex/rnaseq/featurecounts_out/"

# Inicialización del temporizador de tiempo de ejecución
t0 <- proc.time()

# Importar librerías
library(biomaRt)

#### ================================
#### 2. Descarga y formateo de archivo de correspondencias Entrez-Ensembl
#### ================================

# Descargar y descomprimir el archivo gene2ensembl desde NCBI FTP
gene2ensembl_file <- "gene2ensembl"
if (!file.exists(gene2ensembl_file)) {
  url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz"
  destfile <- "gene2ensembl.gz"
  download.file(url, destfile)
  R.utils::gunzip(destfile, remove=FALSE, overwrite=TRUE)
}

# Guardar únicamente las columnas que contienen los Entrez IDs y los Ensembl IDs
ens_corresp <- read.table("gene2ensembl", header=TRUE, sep="\t", quote="", comment.char="", fill=TRUE)
ens_subset <- ens_corresp[, c("GeneID", "Ensembl_gene_identifier")]
colnames(ens_subset) <- c("Entrez", "Ensembl")

# Eliminar genes de especies no humanas
ens_subset$keep_ENSG <- grepl("^ENSG[0-9]+$", ens_subset$Ensembl)
ens_subset <- ens_subset[ens_subset$keep_ENSG == TRUE, ]

# Quitar duplicados de Ensembl 
ens_by_ensembl <- ens_subset[order(ens_subset$Ensembl,
                                   ens_subset$Entrez), ]

ens_unique_ensembl <- ens_by_ensembl[!duplicated(ens_by_ensembl$Ensembl),
                                     c("Entrez", "Ensembl")]

# Quitar duplicados de Entrez priorizando por orden alfabético
ens_by_gene <- ens_unique_ensembl[order(ens_unique_ensembl$Entrez,
                                        ens_unique_ensembl$Ensembl), ]

ens_unique_gene <- ens_by_gene[!duplicated(ens_by_gene$Entrez),
                               c("Entrez", "Ensembl")]

# Guardar tabla de correpondencias entre Ensembl y Entrez
write.table(ens_unique_gene, file = paste0(output,"entrez2ensembl.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#### ================================
#### 3. Añadir columna de Ensembl ID a tabla de conteos de la cuantificación
#### ================================

# Cargar tabla de conteos de genes generada durante la cuantificación
counts <- read.csv2(paste0(input, "gene_counts.csv"), header = TRUE, check.names = FALSE)

# Descartar duplicados de Entrez ID priorizando los que se encuentran en locus principales
counts$keep_NC <- grepl("^NC", counts$gene_chr)

counts_by_chr <- counts[order(!counts$keep_NC), ]

counts_unique <- counts_by_chr[!duplicated(counts_by_chr$gene_id), ]

# Añadir una columna con el ID Ensembl de cada gen
counts_trueID <- merge(counts_unique, ens_unique_gene, by.x="gene_id", by.y="Entrez", all.x=TRUE)

#### ================================
#### 4. Añadir columna con una descripción para cada gen
#### ================================

# Usar useEnsembl para extraer una descripción de cada gen a través del Ensembl ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c("ensembl_gene_id", "description"),
  filters = "ensembl_gene_id",
  values = counts_trueID$Ensembl,
  mart = ensembl
)

# Añadir una columna con las descripciones obtenidas
merged <- merge(counts_trueID, annotations,
                by.x = "Ensembl", by.y = "ensembl_gene_id",
                all.x = TRUE)

merged$keep_NC <- NULL
names(merged)[names(merged) == "description"] <- "gene_description"

#### ================================
#### 5. Conservar solo los cromosomas principales y traducir los identificadores de cromosoma a números
#### ================================

# Extraer de la tabla todos los tipos de identificadores de cromosoma principal
chrs <- as.data.frame(merged$gene_chr)
colnames(chrs) <- "name"
chrs$NC <- grepl("^NC", chrs$name)
chrs_NC <- chrs[chrs$NC == TRUE, ]
chrs_NC$name <- gsub("\\,.*", "", chrs_NC$name)
chrs_NC <- chrs_NC[order(chrs_NC$name), ]
chrs_f <- unique(chrs_NC$name)

# Diccionario de identificadores a nombres de cromosomas
dicc <- c(seq(22),"X","Y","MT")
names(dicc) <- chrs_f

# Traducir identificadores de cromosomas a los valores del diccionario
merged <- merged[startsWith(as.character(merged$gene_chr), "NC"), ]
merged$gene_chr <- gsub("\\,.*", "",  merged$gene_chr)
merged$gene_chr <- dicc[merged$gene_chr]

# Quedarse solo con los cromosomas principales
merged$gene_start <- gsub(",.*", "", merged$gene_start)
merged$gene_end <- gsub(",.*", "", merged$gene_end)
merged$gene_strand <- gsub(",.*", "", merged$gene_strand)

#### ================================
#### 6. Exportar resultados descartando genes que no tengan Ensembl ID
#### ================================

# Crear tabla que servirá para añadir las caracteristicas de los genes en stringtie
merged_entrez <- merged
names(merged_entrez)[1:2] <- c("ensembl_id", "entrez_id")

# Eliminar Entrez ID y eliminar filas que no tengan Ensembl ID para la tabla final
merged_final <- merged[, -2]
names(merged_final)[1] <- "gene_id"
merged_final <- merged_final[!is.na(merged_final$gene_id), ]

# Exportar
write.csv2(merged_entrez, paste0(output, "gene_counts_entrez.csv"), row.names = FALSE)
write.csv2(merged_final, paste0(output, "gene_counts_final.csv"), row.names = FALSE)

# Mostrar tiempo de ejecución
t_total <- (proc.time() - t0)/60
cat("Script ejecutado correctamente.\n")
cat("Tiempo de ejecución:", round(t_total["elapsed"], 2), "minutos\n")
