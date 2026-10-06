# INCLIVA Target-seq pipeline (Nextflow)

Pipeline de rutina para análisis de paneles de secuenciación dirigida (target-seq), adaptado a Nextflow (DSL2) a partir del script bash de rutina. Está basado en una versión modificada de nf-core/Sarek.

**Versión:** 1.0.0

## Descripción

El pipeline procesa FASTQ de paneles de genes y genera:

- Control de calidad y preprocesado de lecturas (FastQC, Cutadapt, Prinseq).
- Mapeo y postprocesado del BAM (BWA-MEM2, Picard, GATK BQSR, LoFreq indelqual).
- Llamada de variantes (Mutect2 y LoFreq), normalización y filtrado.
- VCF final por muestra.

<!-- TODO: completar con anotación (VEP, etc.) si aplica en esta versión -->

## Requisitos

- Linux con Slurm (o el executor que se use en el clúster).
- [Nextflow](https://www.nextflow.io/) (DSL2).
- Conda/Mamba, y Singularity (o Docker).
- Entorno conda: `sarek_nf`.
- Referencias del genoma (hg38) y ficheros del panel, definidos en el fichero de configuración del panel.

## Estructura del repositorio

```text
.
├── workflow_targetseq.nf      # Workflow principal
├── conf/                      # Configuración de módulos/procesos
├── modules/                   # Módulos
├── subworkflows/              # Subworkflows
├── panels/
│   └── panel_configs/         # Plantillas de configuración por panel
│   └── panel_designs/         # Configuración por panel
└── test_data/                 # FASTQ y samplesheet de test (panel TOX3)
```

<!-- TODO: ajustar a la estructura real -->

## Lanzar el pipeline en el clúster

Crea un script bash, por ejemplo `run_pipeline.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=targetseq_pipeline
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=199G
#SBATCH --time=48:00:00
#SBATCH --output=targetseq_%j.out
#SBATCH --error=targetseq_%j.err

source /mnt/zonahpc/home/bioinformatica/comun/software/miniconda3/etc/profile.d/conda.sh
conda activate sarek_nf

nextflow run /ruta/al/repo/workflow_targetseq.nf \
    -c /ruta/al/repo/panels/panel_configs/<panel>.config \
    -profile singularity 
    ```

Y lánzalo con:

```bash
sbatch run_pipeline.sh
```

> El job de Nextflow solo coordina. Los procesos se envían a Slurm según el executor definido en la configuración.

## Test con el panel TOX3

El repositorio incluye FASTQ de test del panel TOX3. Para lanzarlo:

```bash
conda activate sarek_nf

nextflow run workflow_targetseq.nf \
    -c conf_tox3.config \
    -profile singularity
```

<!-- TODO: confirmar nombres de ficheros y parámetros -->

## Analizar otras muestras

### 1. Samplesheet

Crea un `samplesheet.csv` con las siguientes columnas:

| Columna   | Descripción                              |
|-----------|------------------------------------------|
| `patient` | ID del paciente                          |
| `status`  | Estado de la muestra (`<TODO valores>`)  |
| `sample`  | ID de la muestra                         |
| `lane`    | Lane de secuenciación                    |
| `fastq_1` | Ruta absoluta al FASTQ R1                |
| `fastq_2` | Ruta absoluta al FASTQ R2                |

Ejemplo:

```csv
patient,status,sample,lane,fastq_1,fastq_2
P1,0,TOX3-DDPD-1119,L001,/ruta/muestra1_R1.fastq.gz,/ruta/muestra1_R2.fastq.gz
P2,0,TOX3-DDPD-1121,L001,/ruta/muestra2_R1.fastq.gz,/ruta/muestra2_R2.fastq.gz
```

### 2. Fichero de configuración del panel

Selecciona el fichero de configuración que corresponda al panel que quieres analizar. Las plantillas están en:

```text
panels/panel_configs/
```

Copia la plantilla y edita: samplesheet, study_name y out_dir.

## Salida

Dentro de `--outdir`:

| Directorio      | Contenido                                  |
|-----------------|--------------------------------------------|
| `mapping/`      | BAM y BAI por muestra  y `versions.yml`    |
| `QC/rawdata`    | fastqc html y zip por muestra              |
| `variant_calling/` | VCF por muestra y `versions.yml`        |
| `pipeline_info/`| Informes de ejecución de Nextflow          |

## Configuración del entorno

Los paths de caché de Singularity y Conda dependen de cada clúster. No están incluidos en el repositorio.

Cada usuario debe configurarlos en su `nextflow.config` local (por ejemplo `~/.nextflow/config` o un `nextflow.local.config` que no se sube a Git):

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

## Paths databases
Puedes encontrar los paths para las bases de datos en el clúster: /mnt/zonahpc/home/bioinformatica/comun/software/targetseq_pipeline/databases

## Room for improvement (próximas versiones)

Este pipeline es una adaptación inicial a Nextflow (DSL2) de la rutina bash. Las siguientes mejoras están previstas para futuras versiones:

- **Configuraciones de más paneles.** Actualmente solo se incluye la configuración de TOX3. Se añadirán las configs de ONC4, BRC2 y otros paneles usados en la rutina.
- **Módulo de estadísticas.** Implementar un proceso/subworkflow que genere las mismas estadísticas que la rutina (cobertura por región, lecturas mapeadas, duplicados, inserción, etc.) y los plots asociados.
- **Soporte completo de UMIs.** Añadir el manejo de UMIs (umitools o fgbio) en preprocesado y deduplicación, con parámetros `has_umi` y `umi_nt` equivalentes a la rutina.
- **VCFs anotados.** Actualmente el pipeline obtiene los VCF sin anotar; añadir los VCF anotados y los CSV por muestra.
