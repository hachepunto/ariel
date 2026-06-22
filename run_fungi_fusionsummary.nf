#!/usr/bin/env nextflow
/*
 * Re-ejecuta Fungi y FusionSummary usando los archivos ya publicados en resultsDir.
 * No agregar al repositorio.
 *
 * Uso:
 *   nextflow run ~/run_fungi_fusionsummary.nf \
 *     --resultsDir      /ruta/a/resultados \
 *     --runSampleSheet  /ruta/a/sample_sheet.tsv \
 *     --referenceDir    /ruta/a/GRCh38_no_alt \
 *     --arielDir        /ruta/a/ariel          # directorio del repo (default: ~/ARIEL)
 *
 * Parametros opcionales (con sus defaults):
 *   --threadsFungi          15
 *   --fungiFilterEnsembl    "same_gene,homologs"
 *   --fungiFilterDb         "banned,paralog"
 *   --fungiFilterMinCount   0
 */

nextflow.enable.dsl=2

params.arielDir             = "$HOME/pipelines/ariel"
params.threadsFungi         = 50
params.fungiFilterEnsembl   = "same_gene,homologs"
params.fungiFilterDb        = "banned,paralog"
params.fungiFilterMinCount  = 0

include { Fungi         } from '/scratch/home/acerverat/pipelines/ariel/modules/modules.nf'
include { FusionSummary } from '/scratch/home/acerverat/pipelines/ariel/modules/modules.nf'

workflow {

  cicero_ch = Channel
    .fromPath("${params.resultsDir}/fusions/cicero/*_final_fusions.txt")
    .map { f -> [f.simpleName.replace('_final_fusions', ''), 'cicero', f] }

  arriba_ch = Channel
    .fromPath("${params.resultsDir}/fusions/arriba/*_fusions.tsv")
    .map { f -> [f.simpleName.replace('_fusions', ''), 'arriba', f] }

  fusioncatcher_ch = Channel
    .fromPath("${params.resultsDir}/fusions/fusioncatcher/*_final-list_candidate-fusion-genes.txt")
    .map { f -> [f.name.replace('_final-list_candidate-fusion-genes.txt', ''), 'fusioncatcher', f] }

  fusions_ch = cicero_ch
    .concat(arriba_ch, fusioncatcher_ch)
    .collect()

  Fungi(fusions_ch)

  rascall_ch = Channel
    .fromPath("${params.resultsDir}/rascall/**/*_final_variants.csv")
    .collect()

  FusionSummary(
    file(params.runSampleSheet),
    Fungi.out.bp,
    rascall_ch
  )
}
