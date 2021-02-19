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


// Validate inputs

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }
ch_genome = file("${cluster}/References/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa", checkIfExists: true)

if (!params.outdir) {
  params.outdir = params.run
}

// Mandatory arguments for the publishDir option (cluster_path has a default but project is mandatory as input)

cluster_path = params.cluster_path
project = params.project


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
   publishDir "${cluster_path}/05_QC/${project}/samplesheet/", mode: params.publish_dir_mode

   input:
   path samplesheet from ch_input

   output:
   path "samplesheet_bed.csv" into ch_samplesheet

   script:
   out = "samplesheet_bed.csv"

   """
   modify_samplesheet.py $samplesheet $out
   """
 }



 def validate_input(LinkedHashMap sample) {
     def sample_id = sample.sampleID
     def bam = sample.bam
     def protocol = sample.protocol
     def bed = sample.bed
     def interval = sample.interval

     def array = []
     array = [ sample_id, protocol, file(bam, checkIfExists: true), file(bed, checkIfExists: true), file(interval, checkIfExists: true) ]

     return array
 }

 /*
  * Create channels for input fastq files
  */
 ch_samplesheet
     .splitCsv(header:true, sep:',')
     .map { validate_input(it) }
     .set { ch_bam_index }



/*ch_samplesheet
 .splitCsv(header:true, sep:',')
 .map { it = ["${it[0]}", "${it[1]}", "${it[2]}", "${it[3]}"]}
 .set { ch_bam_index }
 */



process samtools {
  tag "$sample"
  label 'process_low'

  input:
  set val(sample), val(experiment), path(bam), path(bed), path(interval) from ch_bam_index

  output:
  set val(sample), path(bam), path("*.bam.bai"), val(experiment), path(bed), path(interval) into ch_mosdepth

  script:
  """
  samtools index $bam
  """
}


 /*
 * STEP 1 - Mosdepth to get the coverage
 */

 process mosdepth {
   tag "$sample"
   label 'process_low'
   publishDir "${cluster_path}/05_QC/${project}/mosdepth/", mode: params.publish_dir_mode


   input:
   tuple val(sample), path(bam), path(bai), val(experiment), path(bed), path(interval) from ch_mosdepth

   output:
   tuple val(sample), path(bam), path(bai), path(interval) into ch_picard_hsmetrics
   tuple val(sample), path("*.thresholds.bed") into ch_ontarget_coverage
   path "*.mosdepth.global.dist.txt" into ch_plot_distances, ch_mosdepth_mqc
   path "*.{txt,gz,csi}"

   script:

   prefix = "${sample}"
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
   """
 }


process picard_hsmetrics {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/05_QC/${project}/HSmetrics/", mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bam), path(bai), path(interval) from ch_picard_hsmetrics
  file(genome) from ch_genome

  output:
  path("*.hybrid_selection_metrics.txt") into ch_merge_metrics

  script:
  outfile = sample + ".hybrid_selection_metrics.txt"

  """
  picard -Xmx${task.memory.toGiga()-mem_adjust}G CollectHsMetrics \
  INPUT=$bam \
  OUTPUT=$outfile \
  TARGET_INTERVALS=$interval \
  BAIT_INTERVALS=$interval \
  REFERENCE_SEQUENCE=$genome
  """
}


process merge_metrics {
label 'process_low'
publishDir "${cluster_path}/05_QC/${project}/mergedHSmetrics/", mode: params.publish_dir_mode

input:
path("metrics/*") from ch_merge_metrics.collect().ifEmpty[()]

output:
path("hybrid_selection_metrics.txt")

script:
"""
metrics=(\$(echo -e `awk 'FNR>1 || NR==1' metrics/*.txt`))
echo \$metrics > hybrid_selection_metrics.txt
"""
}


process ontarget_coverage {
  tag "$sample"
  label 'process_low'
  publishDir "${cluster_path}/05_QC/${project}/coverageTables/", mode: params.publish_dir_mode

  input:
  tuple val(sample), path(bed) from ch_ontarget_coverage

  output:
  path("*summary.tsv") into ch_collect_tables

  script:
  out = "${sample}_summary.tsv"

  """
  mosdepth_table.py $bed $out -s $sample
  """
}


process cat_summary {
  label 'process_low'
  publishDir "${cluster_path}/05_QC/${project}/coverageMergedTables/", mode: params.publish_dir_mode

  input:
  path("tables/*") from ch_collect_tables.collect().ifEmpty([])

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
  """
}


// MultiQC

process multiqc {
    label 'process_low'
    publishDir "${cluster_path}/05_QC/${project}/multiqc", mode: params.publish_dir_mode,

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
    title = custom_runName ? "--title \"$project\"" : ''
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
