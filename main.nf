#!/usr/bin/env nextflow
/*
========================================================================================
                         mpozuelo/exome_mosdepth
========================================================================================
mpozuelo/exome_mosdepth Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/mpozuelo/exome_mosdepth
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info mpozueloHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run mpozuelo/exome_mosdepth --input '*.csv' -profile docker

    Mandatory arguments:
      --input [file]                Samplesheet samples information (only samples that want to show in the same MultiQC)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity.

    Other options
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


 // Has the run name been specified by the user?
 //  this has the bonus effect of catching both -name and --name
 custom_runName = params.name
 if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
   custom_runName = workflow.runName
 }
 else{
   workflow.runName = params.user + " " + params.timestamp
   custom_runName = workflow.runName
 }


if (!params.outdir) {
  params.outdir = params.run
}

// Mandatory arguments for the publishDir option (cluster_path has a default but project is mandatory as input)
cluster_path = params.cluster_path
project = params.project


// Validate inputs

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }
ch_genome = file("${cluster_path}/References/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", checkIfExists: true)
ch_genome_index = file("${cluster_path}/References/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai", checkIfExists: true)
//ch_genome = file("${cluster_path}/References/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa", checkIfExists: true)
//ch_genome_index = file("${cluster_path}/References/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa.fai", checkIfExists: true)

// Stage multiqc config files
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_image_docs            = file("$baseDir/assets/figures/Logo_IdisNA_CIMA.png", checkIfExists: true)



// Header log info
log.info mpozueloHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Project'] = params.project
summary['Input'] = params.input
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['User'] = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'mpozuelo-exome_mosdepth-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'mpozuelo/exome_mosdepth Workflow Summary'
    section_href: 'https://github.com/mpozuelo/exome_mosdepth'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}




 /*
 * Parse software version numbers
 */
 process get_software_versions {
   publishDir "${params.outdir}/pipeline_info", mode: 'copy',
   saveAs: { filename ->
     if (filename.indexOf(".csv") > 0) filename
     else null
   }

   output:
   file 'software_versions_mqc.yaml' into software_versions_yaml
   file "software_versions.csv"

   script:
   """
   echo $workflow.manifest.version &> v_ngi_QC.txt
   echo $workflow.nextflow.version &> v_nextflow.txt
   mosdepth --version > v_mosdepth.txt
   scrape_software_versions.py &> software_versions_mqc.yaml
   """
 }


 /*
 * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing

 /*
 * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing
 It contains the following columns:
 1- Name given to the sample
 2- Path to file
 3- Experiment
 */


 process modify_samplesheet {
   publishDir "${cluster_path}/data/05_QC/${project}/samplesheet/", mode: params.publish_dir_mode

   input:
   path samplesheet from ch_input

   output:
   path "samplesheet_bed.csv" into ch_samplesheet, ch_samplesheet_percentage

   script:
   out = "samplesheet_bed.csv"

   """
   modify_samplesheet.py $samplesheet $out
   """
 }


def count_bam(LinkedHashMap sample) {
  def sample_id = sample.sampleID
  def bam = sample.bam

  def array = []
  array = [ sample_id, file(bam, checkIfExists: true) ]
}


ch_samplesheet
  .splitCsv(header:true, sep:',')
  .map{ count_bam(it) }
  .set { ch_bam_count }



process count_total {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/data/05_QC/${project}/numberReadsBAM/", mode: params.publish_dir_mode

  input:
  set val(sample), file(bam) from ch_bam_count

  output:
  path ("*.tsv") into ch_total_reads

  script:
  """
  samtools view -c $bam > "${sample}.tsv"
  printf "%s\t%s\n" "$sample" "\$(echo \$(cat ${sample}.tsv))" > "${sample}_total_reads.tsv"
  """
}


process percentages {
  publishDir "${cluster_path}/data/05_QC/${project}/", mode: params.publish_dir_mode,
  saveAs: { filename ->
      filename.endsWith(".tsv") ? "numberReadsBAM/$filename" : "samplesheet/$filename"
    }

  input:
  path samplesheet from ch_samplesheet_percentage
  path ("totalReads/*") from ch_total_reads.collect().ifEmpty([])

  output:
  path "samplesheet_percentage.csv" into ch_samplesheet_count
  path "total_counts.tsv"

  script:
  out = "samplesheet_percentage.csv"

  """
  printf "sampleID\tcounts\n" > total_counts.tsv
  cat totalReads/*.tsv >> total_counts.tsv
  percentage_samplesheet.py $samplesheet $out -c total_counts.tsv
  """
}




 def validate_input(LinkedHashMap sample) {
     def sample_id = sample.sampleID
     def bam = sample.bam
     def protocol = sample.protocol
     def bed = sample.bed
     def interval = sample.interval
     def percentage = sample.percentage

     def array = []
     array = [ sample_id, protocol, file(bam, checkIfExists: true), file(bed, checkIfExists: true), file(interval, checkIfExists: true), percentage ]

     return array
 }

 /*
  * Create channels for input fastq files
  */
 ch_samplesheet_count
     .splitCsv(header:true, sep:',')
     .map { validate_input(it) }
     .set { ch_samtools }



/*ch_samplesheet
 .splitCsv(header:true, sep:',')
 .map { it = ["${it[0]}", "${it[1]}", "${it[2]}", "${it[3]}"]}
 .set { ch_bam_index }
 */



process samtools {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/data/05_QC/${project}/samtools/${sample}", mode: params.publish_dir_mode


  input:
  set val(sample), val(experiment), path(bam), path(bed), path(interval), val(percentage) from ch_samtools

  output:
  set val(sample), path(bam), path("${bam.baseName}.bam.bai"), path("${bam.baseName}_${percentage}.bam"), path("${bam.baseName}_${percentage}.bam.bai"), val(experiment), path(bed), path(interval) into ch_mosdepth
  set val(sample), path(bam), path("${bam.baseName}.bam.bai"), path("${bam.baseName}_${percentage}.bam"), path("${bam.baseName}_${percentage}.bam.bai"), path(interval) into ch_picard_hsmetrics,
                                                                                                                                                                             ch_picard_alignmentmetrics

  script:
  subset = "${bam.baseName}_${percentage}.bam"
  """
  samtools index $bam
  samtools view -s $percentage -b $bam > $subset
  samtools index $subset
  """
}


 /*
 * STEP 1 - Mosdepth to get the coverage
 */

 process mosdepth {
   tag "$sample"
   label 'process_low'
   publishDir "${cluster_path}/data/05_QC/${project}/mosdepth/${sample}", mode: params.publish_dir_mode


   input:
   set val(sample), path(bam), path(bai), path(bam_subset), path(bai_subset), val(experiment), path(bed), path(interval) from ch_mosdepth

   output:
   set val(sample), path("${prefix}.thresholds.bed"), path("${prefix_subset}.thresholds.bed") into ch_ontarget_coverage
   //path "*.mosdepth.global.dist.txt" into ch_mosdepth_mqc
   path "*.{txt,gz,csi}"

   script:

   prefix = "${bam.baseName}"
   prefix_subset = "${bam_subset.baseName}"
   threshold = params.threshold

   """
   mosdepth \\
   --by $bed \\
   --fast-mode \\
   --thresholds 0,1,5,10,25,50,100,200,300,400,500 \\
   -Q 20 \\
   -t $task.cpus \\
   ${prefix} \\
   $bam

   pigz -dk ${prefix}.thresholds.bed.gz

   mosdepth \\
   --by $bed \\
   --fast-mode \\
   --thresholds 0,1,5,10,25,50,100,200,300,400,500 \\
   -Q 20 \\
   -t $task.cpus \\
   ${prefix_subset} \\
   $bam_subset

   pigz -dk ${prefix_subset}.thresholds.bed.gz
   """
 }


process picard_hsmetrics {
  tag "$sample"
  label 'process_medium'
  publishDir "${cluster_path}/data/05_QC/${project}/HSmetrics/${sample}", mode: params.publish_dir_mode

  input:
  set val(sample), path(bam), path(bai), path(bam_subset), path(bai_subset), path(interval) from ch_picard_hsmetrics
  file(genome) from ch_genome
  file(index) from ch_genome_index

  output:
  path("${bam.baseName}.hybrid_selection_metrics.txt") into ch_merge_original_HSmetrics
  path("${bam_subset.baseName}.hybrid_selection_metrics.txt") into ch_merge_subset_HSmetrics

  script:
  outfile = "${bam.baseName}.hybrid_selection_metrics.txt"
  outfile_subset = "${bam_subset.baseName}.hybrid_selection_metrics.txt"
  java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""

  """
  picard ${java_options} CollectHsMetrics \
  INPUT=$bam \
  OUTPUT=$outfile \
  TARGET_INTERVALS=$interval \
  BAIT_INTERVALS=$interval \
  REFERENCE_SEQUENCE=$genome \
  TMP_DIR=tmp

  picard ${java_options} CollectHsMetrics \
  INPUT=$bam_subset \
  OUTPUT=$outfile_subset \
  TARGET_INTERVALS=$interval \
  BAIT_INTERVALS=$interval \
  REFERENCE_SEQUENCE=$genome \
  TMP_DIR=tmp
  """
}

process picard_alignmentmetrics {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/data/05_QC/${project}/HSmetrics/${sample}", mode: params.publish_dir_mode

  input:
  set val(sample), path(bam), path(bai), path(bam_subset), path(bai_subset), path(interval) from ch_picard_alignmentmetrics
  file(genome) from ch_genome
  file(index) from ch_genome_index

  output:
  path("${bam.baseName}.alignment_metrics.txt") into ch_merge_original_Alignmentmetrics
  path("${bam_subset.baseName}.alignment_metrics.txt") into ch_merge_subset_Alignmentmetrics

  script:
  outfile = "${bam.baseName}.alignment_metrics.txt"
  outfile_subset = "${bam_subset.baseName}.alignment_metrics.txt"
  java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2 )+"g "+ "-Xmx" + (task.memory.toGiga() - 1)+ "g\""

  """
  picard ${java_options} CollectAlignmentSummaryMetrics \
  INPUT=$bam \
  OUTPUT=$outfile \
  R=$genome

  picard ${java_options} CollectAlignmentSummaryMetrics \
  INPUT=$bam_subset \
  OUTPUT=$outfile_subset \
  R=$genome
  """
}


process merge_metrics {
  label 'process_low'
  publishDir "${cluster_path}/data/05_QC/${project}/mergedmetrics/", mode: params.publish_dir_mode

  input:
  path("originalHSMetrics/*") from ch_merge_original_HSmetrics.collect()
  path("subsetHSMetrics/*") from ch_merge_subset_HSmetrics.collect()
  path("originalAlignmentMetrics/*") from ch_merge_original_Alignmentmetrics.collect()
  path("subsetAlignmentMetrics/*") from ch_merge_subset_Alignmentmetrics.collect()

  output:
  path("*_metrics.txt")

  script:
  """
  merge_Metrics.sh "originalHSMetrics" "Original.hybrid_selection_metrics_tmp.txt"
  merge_Metrics.sh "subsetHSMetrics" "Subset.hybrid_selection_metrics_tmp.txt"
  merge_Metrics.sh "originalAlignmentMetrics" "Original.Alignment_metrics_tmp.txt"
  merge_Metrics.sh "subsetAlignmentMetrics" "Subset.Alignment_metrics_tmp.txt"
  """
}


process ontarget_coverage {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/data/05_QC/${project}/coverageTables/", mode: params.publish_dir_mode

  input:
  set val(sample), path(bed), path(bed_subset) from ch_ontarget_coverage

  output:
  set path("${sample}_summary.tsv"), path("${sample}_subset_summary.tsv") into ch_collect_tables

  script:
  out = "${sample}_summary.tsv"
  out_subset = "${sample}_subset_summary.tsv"
  sample_subset = "${sample}_subset"

  """
  mosdepth_table.py $bed $out -s $sample
  mosdepth_table.py $bed_subset $out_subset -s $sample_subset
  """
}


process cat_summary {
  label 'process_low'
  publishDir "${cluster_path}/data/05_QC/${project}/coverageMergedTables/", mode: params.publish_dir_mode

  input:
  set path("tables/*"), path("subset_tables/") from ch_collect_tables.collect().ifEmpty([])

  output:
  path("*.tsv")

//  printf "%s\s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "0" "1" "5" "10" "20" "40" "80" "100" "200" "300" "400" "500" "600" > coverages.tsv
//"0,1,5,10,25,50,100,200,300,400,500"
  script:
  """
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "" "1X" "5X" "10X" "25X" "50X" "100X" "200X" "300X" "400X" "500X" > coverages.tsv

  for f in \$(find tables -name "*.tsv"); do \\
  cat \$f >> coverages.tsv;
  done


  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "" "1X" "5X" "10X" "25X" "50X" "100X" "200X" "300X" "400X" "500X" > coverages_subset.tsv

  for f in \$(find subset_tables -name "*.tsv"); do \\
  cat \$f >> coverages_subset.tsv;
  done
  """
}


// MultiQC
/*
process multiqc {
    label 'process_low'
    publishDir "${cluster_path}/data/05_QC/${project}/multiqc", mode: params.publish_dir_mode,

    when:
    !params.skip_multiqc

    input:
    path (multiqc_config) from ch_multiqc_config
    path workflow_summary from create_workflow_summary(summary)
    path ('mosdepth/*') from ch_mosdepth_mqc.collect().ifEmpty([])
    path "*" from ch_image_docs

    output:
    path "*multiqc_report.html" into ch_multiqc_report
    path "*_data"
    path "*.tsv"

    script:
    rtitle = custom_runName ? "--title \"$project\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : "$project"
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[mpozuelo/exome_mosdepth] Successful: $workflow.runName"

    if (!workflow.success) {
      subject = "[mpozuelo/exome_mosdepth] FAILED: $workflow.runName"
    }



    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";


    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[mpozuelo/exome_mosdepth]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[mpozuelo/exome_mosdepth]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def mpozueloHeader() {
  // Log colors ANSI codes
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";


  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  __  __  __   __  ___         ${c_reset}
  ${c_blue}  | \\/ | |__| |  |  /  |  |     ${c_reset}
  ${c_blue}  |    | |    |__| /__ |__|         ${c_reset}
  ${c_white}  mpozuelo/exome_mosdepth v${workflow.manifest.version}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}


def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
          "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
          "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
          "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
          "============================================================"
        }
      }
    }
  }
}
