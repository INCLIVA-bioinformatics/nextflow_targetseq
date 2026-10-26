<p align="center"><img src="/mnt/zonahpc/home/bioinformatica/comun/software/targetseq_pipeline/workflow_targetseq.nf" width="600"></p>

# INCLIVA Target-seq pipeline (Nextflow)

Pipeline de rutina para análisis de paneles de secuenciación dirigida (target-seq), adaptado a Nextflow (DSL2) a partir del script bash de rutina (`inc_pipe_for_panels_only_tumors.sh`). Está basado en una versión modificada de nf-core/Sarek.

**Versión:** 1.2.0

## Descripción

El pipeline procesa FASTQ de paneles de genes en tres subworkflows:

1. **Preprocesado** (`subworkflows/preprocessing.nf`)
   - FastQC sobre los FASTQ crudos.
   - Cutadapt: eliminación de adaptadores y descarte de lecturas más cortas que `min_read_length`.
   - Cutadapt: filtrado por primers (solo si `filter_by_primer = true`).
   - Prinseq: trimado por calidad.
   - FastQC sobre las lecturas preprocesadas.
2. **Mapeo y postprocesado** (`subworkflows/mapping_postprocessing.nf`)
   - BWA-MEM2 + filtro por calidad de mapeo (`samtools view -q min_mapping_quality`) y ordenación.
   - Según el tipo de panel:
     - **Amplicones** (`amplicon = true`): soft-clipping de primers con BAMClipper (no se eliminan duplicados).
     - **Captura** (`amplicon = false`): eliminación de duplicados con Picard MarkDuplicates.
   - Picard CleanSam.
   - GATK BaseRecalibrator + ApplyBQSR (sobre el BAM completo).
   - LoFreq indelqual.
   - Profundidad en las regiones target (samtools depth) y generación de los BED de regiones cubiertas y no cubiertas a la cobertura mínima `cov`.
3. **Llamada de variantes** (`subworkflows/variant_calling.nf`)
   - Mutect2 + FilterMutectCalls y LoFreq.
   - Normalización, selección de variantes, anotación con información del BAM, ordenación e indexado.
   - Unión de los VCF de ambos callers en un VCF final por muestra.
   - Unión de los VCF de todas las muestras en un VCF multi-muestra del run (`bcftools merge`).
   - **Anotación** (`subworkflows/annotation.nf`, equivalente a `annotate_variants` de la rutina; solo si se indica `vep_cache`):
     - VEP 110 con SpliceAI (`spliceai_cutoff`).
     - Cancer hotspots, COSMIC v91 (con paso del ID de COSMIC a la anotación), mixed_db, base de datos interna (`inc_db`), conservación y OncoKB (SnpSift).
     - Eliminación de registros duplicados (vcfuniq) y conversión a CSV (aminoácidos en código de una letra; tránscritos de `list_of_transcripts` si el panel la tiene).
     - Matriz de variantes compartidas entre muestras (`sample_contamination.csv`) y CSV de variantes por muestra.

Para reproducir los resultados de la rutina bash se fijan las mismas versiones de GATK (4.4.0.0) y LoFreq (2.1.3.1), el mismo read group (`ID:AgilentXTHS`) y el tamaño de lote de BWA-MEM2 (`-K`), de modo que el resultado no depende del número de CPUs.

### Scripts de rutina modificados

Los scripts de `bin/` son los de la rutina, salvo los de estadísticas por amplicón, que se han modificado sin cambiar su salida (`STATS_BY_AMPLICON` en `modules/incliva/statistics/main.nf`):

| Script del pipeline | Script de rutina | Se usa en |
|---|---|---|
| `bin/inc_amplicons_BRCA_universal_v1.1.py` | `inc_amplicons_BRCA_universal_v1.0.py` | Paneles con `brca = true` (BRC2) |
| `bin/inc_amplicons_TP53_universal_v1.1.py` | `inc_amplicons_TP53_universal_v1.0.py` | Panel TP53 |

- **Optimización (ambos scripts).** Para cada muestra, la versión de rutina recorría el SAM completo una vez por cada amplicón del panel (en BRC2, alrededor de un millón de lecturas por muestra, leídas 219 veces). La v1.1 lee el SAM una sola vez por muestra, agrupa las lecturas por cromosoma con su posición y su longitud (el máximo de la CIGAR) y, para cada amplicón, busca solo las lecturas que caen en su rango. Los criterios de selección no cambian (margen de ±40 pb en BRCA y ±50 pb en TP53; en BRCA, longitud de al menos el 60 % del amplicón), ni el cálculo posterior de profundidad y estadísticas. En BRC2, `STATS_BY_AMPLICON` pasa de 28,5 min a 1 min, con un TSV idéntico byte a byte al de la v1.0.

## Requisitos

- Linux con Slurm (executor definido en `config/base.config`).
- [Nextflow](https://www.nextflow.io/) (DSL2).
- Conda/Mamba y Singularity.
- Entorno conda: `sarek_nf`.
- Referencias del genoma (GRCh38) y ficheros del panel, definidos en el fichero de configuración del panel.

## Estructura del repositorio

```text
.
├── workflow_targetseq.nf      # Workflow principal
├── incliva_modules.nf         # Procesos propios de INCLIVA
├── config/
│   └── base.config            # Recursos, argumentos y publishDir de cada proceso
├── modules/
│   ├── incliva/               # Módulos propios o adaptados
│   ├── local/                 # Módulos locales heredados de Sarek
│   └── nf-core/               # Módulos de nf-core
├── subworkflows/              # preprocessing, mapping_postprocessing, variant_calling
├── panels/
│   ├── panel_configs/         # Configuración por panel (conf_<PANEL>.config)
│   └── panel_designs/         # Diseño de cada panel (BED, primers, amplicones)
└── test_data/                 # FASTQ y samplesheet de test (panel TOX3)
```

## Paneles disponibles

| Panel | Tipo       | Filtrado por primers | Eliminación duplicados / clipping | Config                                  |
|-------|------------|----------------------|-----------------------------------|-----------------------------------------|
| TOX3  | Amplicones | Sí (error 0.3)       | BAMClipper (19/19)                | `panels/panel_configs/conf_TOX3.config` |
| BRC2  | Amplicones | Sí (error 0.2)       | BAMClipper (18/18)                | `panels/panel_configs/conf_BRC2.config` |
| ONC4  | Captura    | No                   | Picard MarkDuplicates             | `panels/panel_configs/conf_ONC4.config` |
| TP53  | Amplicones | Sí (error 0.3)       | BAMClipper (18/18)                | `panels/panel_configs/conf_TP53.config` |

TP53 usa `af = 0.01` (Proyecto HEMATO). Para muestras de Oncocarta, lánzalo con `--af 0.02`.

Cada diseño en `panels/panel_designs/<PANEL>/` contiene, entre otros, `capture.bed`, `target.bed` y, para paneles de amplicones, `R1_primers.fasta`, `R2_primers.fasta` y `amplicons.pe.bed` (usado por BAMClipper).

## Lanzar el pipeline en el clúster

Crea un script bash, por ejemplo `run_pipeline.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=targetseq_pipeline
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=targetseq_%j.out
#SBATCH --error=targetseq_%j.err

source /mnt/zonahpc/home/bioinformatica/comun/software/miniconda3/etc/profile.d/conda.sh
conda activate sarek_nf

nextflow run /ruta/al/repo/workflow_targetseq.nf \
    -c /ruta/al/repo/panels/panel_configs/conf_<PANEL>.config \
    -profile singularity
```

Y lánzalo con:

```bash
sbatch run_pipeline.sh
```

> El job de Nextflow solo coordina: cada proceso se envía a Slurm como un job independiente según el executor y los recursos definidos en `config/base.config`. Por eso el job principal necesita pocos recursos.

Para retomar una ejecución interrumpida, añade `-resume` al comando.

## Test con el panel TOX3

El repositorio incluye FASTQ de test del panel TOX3 (`test_data/rawdata_test/TOX3/`) y su samplesheet (`test_data/samplesheet_TOX3.csv`), que ya está referenciado en `conf_TOX3.config`. Para lanzarlo:

```bash
conda activate sarek_nf

nextflow run workflow_targetseq.nf \
    -c panels/panel_configs/conf_TOX3.config \
    -profile singularity
```

Los resultados se escriben en `TOX3_outdir/` dentro del directorio de lanzamiento.

## Analizar otras muestras

### 1. Samplesheet

Crea un `samplesheet.csv` con las siguientes columnas:

| Columna   | Descripción                                        |
|-----------|----------------------------------------------------|
| `patient` | ID del paciente                                    |
| `status`  | Estado de la muestra (`0` = normal, `1` = tumor)   |
| `sample`  | ID de la muestra                                   |
| `lane`    | Lane de secuenciación                              |
| `fastq_1` | Ruta absoluta al FASTQ R1                          |
| `fastq_2` | Ruta absoluta al FASTQ R2                          |

Ejemplo:

```csv
patient,status,sample,lane,fastq_1,fastq_2
TOX3-DDPD-1119,1,TOX3-DDPD-1119_S3,1,/ruta/TOX3-DDPD-1119_S3_L001_R1_001.fastq.gz,/ruta/TOX3-DDPD-1119_S3_L001_R2_001.fastq.gz
TOX3-DDPD-1120,1,TOX3-DDPD-1120_S4,1,/ruta/TOX3-DDPD-1120_S4_L001_R1_001.fastq.gz,/ruta/TOX3-DDPD-1120_S4_L001_R2_001.fastq.gz
```

### 2. Fichero de configuración del panel

Copia la configuración del panel que quieres analizar desde `panels/panel_configs/` y edita al menos:

- `samplesheet`: ruta al samplesheet.
- `study_name`: nombre del estudio/run.
- `outdir`: directorio de salida.

Parámetros de procesado principales (ya ajustados para cada panel):

| Parámetro             | Descripción                                                      |
|-----------------------|------------------------------------------------------------------|
| `amplicon`            | `true` para paneles de amplicones (BAMClipper), `false` para captura (MarkDuplicates) |
| `filter_by_primer`    | Filtrar lecturas por primers con Cutadapt                        |
| `cut_primer_error` / `cut_primer_overlap` | Tolerancia y solapamiento mínimo del filtrado por primers |
| `min_read_length`     | Longitud mínima de lectura tras el trimado                       |
| `min_trim_qual`       | Calidad mínima para el trimado con Prinseq                       |
| `min_mapping_quality` | MAPQ mínima tras el mapeo                                        |
| `up_bamclip` / `down_bamclip` | Ventana upstream/downstream de BAMClipper                |
| `cov` / `af`          | Cobertura mínima y frecuencia alélica mínima                     |
| `has_umi`             | Uso de UMIs (todavía no implementado)                            |

## Salida

Dentro de `outdir`:

| Directorio          | Contenido                                                       |
|---------------------|-----------------------------------------------------------------|
| `QC/rawdata/`       | Informes FastQC (html y zip) por muestra                         |
| `QC/mapping/`       | BED de regiones cubiertas y no cubiertas a `cov` x por muestra (`<sample>_covered_regions_<cov>x.bed`, `<sample>_non_covered_<cov>x.bed`); con anotación, `<study_name>_sample_contamination.csv` |
| `mapping/`          | BAM y BAI finales por muestra y `versions.yml`                   |
| `variant_calling/`  | VCF final por muestra y VCF multi-muestra del run; con anotación, además `<study_name>_annotated.vcf.gz`, `<study_name>_annotated.csv` y un CSV por muestra |
| `pipeline_info/`    | Informes de ejecución de Nextflow (report, timeline, trace, DAG) |

Por defecto los ficheros se publican como enlaces simbólicos (`publish_dir_mode = "symlink"`) al directorio `work/`. Si se va a borrar `work/`, usa `publish_dir_mode = "copy"`.

## Configuración del entorno

Los paths de caché de Singularity y Conda dependen de cada clúster. Las configuraciones de panel incluyen los del clúster de INCLIVA; si trabajas en otro entorno, ajústalos en la config del panel o en tu `nextflow.config` local (por ejemplo `~/.nextflow/config`):

```groovy
singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = "/ruta/a/tu/singularity_cacheDir"
}

conda {
    enabled  = true
    useMamba = false
    cacheDir = "/ruta/a/tu/.nextflow/cache"
}
```

Asegúrate de que los directorios de caché existen y son accesibles desde los nodos de cómputo.

## Bases de datos y referencias

Las bases de datos están en el clúster (no se suben a Git):

- `/mnt/zonahpc/home/bioinformatica/comun/databases`
  - Genoma de referencia e índice BWA-MEM2: `bwa_mem2/Homo_sapiens.GRCh38.fa` (+ `.fai`, `.dict` e índices).
  - Mills and 1000G gold standard indels (hg38).
  - `inc_db.vcf` (base de datos interna).
  - dbSNP 153: `dbSNP/153/All_20200219.vcf.gz`
  - Anotación (todos con su `.tbi`):
    - Caché de VEP 110 (con `Plugins/`): `v110_last/`
    - SpliceAI: `SpliceAI/spliceai_scores.raw.{snv,indel}.hg38.vcf.gz`
    - Cancer hotspots: `cancer_hotspots/20180723/cancerhospots.bed.gz`
    - COSMIC v91: `cosmicv91/AllMutations_COSMIC_v91.vcf.gz`
    - Mixed DB: `DB/mixed_db.vcf.gz`
    - OncoKB: `dbSNP155.cosmic.oncokb.noIDs.vcf.gz`

## Room for improvement (próximas versiones)

Este pipeline es una adaptación inicial a Nextflow (DSL2) de la rutina bash. Las siguientes mejoras están previstas para futuras versiones:

- **Más paneles.** Actualmente se incluyen TOX3, BRC2, ONC4 y TP53. Se añadirán el resto de paneles usados en la rutina.
- **Soporte completo de UMIs.** Añadir el manejo de UMIs (umi_tools o fgbio) en preprocesado y deduplicación, con parámetros `has_umi` y `umi_nt` equivalentes a la rutina.
- **MultiQC**