# TFM_RNAseq_Pipeline
El presente repositorio recoge el conjunto de scripts que componen el pipeline diseñado durante las Prácticas de Emptresa y el Trabajo Fin de Máster. Además, también contiene algunos resultados generados durante su ejecución para 14 muestras del proyecto SRP017465.

## Distribución de contenidos
- Contenido de la asignatura Prácticas de Empresa: Etapas A, B y D
- Contenido del Trabajo Fin de Máster: Etapas C, E, F, G y H

## Orden del flujo de trabajo
- A. Control de calidad
  - 1. `01_qc`
  - 2. `02_export_qc`
  - 3. `03_check_adapters`
- B. Alineamiento al genoma de referencia  
  - 4. `04_mapping`  
  - 5. `05_export_mapping`  
- C. Cuantificación de la expresión  
  - 6. `06_fix_gtf`  
  - 7. `07_quantification`  
  - 8. `08_quantification_postprocessing.R`  
  - 9. `09_fpkm_tpm.R`  
- D. Ensamblaje de transcritos  
  - 10. `10_assembly`  
- E. Análisis de expresión diferencial  
  - 11. `11_de.R`  
- F. Análisis de enriquecimiento GO/KEGG  
  - 12. `12_enrichment.R`  
- G. Análisis GSEA  
  - 13. `13_gsea.R`  
- H. Análisis de interacciones proteína-proteína  
  - 14. `14_ppi.R`  
