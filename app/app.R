library(shiny)
library(bslib)

app_dir <- normalizePath('.')
source(file.path(app_dir,'R','engine.R'),local=TRUE)
source(file.path(app_dir,'R','examples.R'),local=TRUE)
source(file.path(app_dir,'R','report.R'),local=TRUE)
if(!all(c('barbac_palettes','barbac_palette') %in% getNamespaceExports('barbac')))
  stop('Studio needs the current barbac release. Launch with Rscript app/run.R.')
options(shiny.maxRequestSize=as.numeric(Sys.getenv('BARBAC_STUDIO_UPLOAD_MB','512'))*1024^2)
# Explicit background processes even on hosts that report a single CPU.
future::plan(future::multisession,workers=I(as.integer(Sys.getenv('BARBAC_STUDIO_WORKERS','2'))))
options(future.globals.maxSize=1024^3)

mark <- function() HTML('<svg viewBox="0 0 40 40" fill="none" aria-hidden="true"><path d="M8 7v26M16 12v16M24 7v26M32 12v16" stroke="currentColor" stroke-width="4" stroke-linecap="round"/><path d="M8 14h16M16 26h16" stroke="#d8bf87" stroke-width="3"/></svg>')
arrow <- function() span('\u2197',class='arrow',`aria-hidden`='true')
eyebrow <- function(text) div(class='eyebrow',text)
panel_heading <- function(number,title,subtitle=NULL) div(class='panel-heading',span(number,class='step-number'),div(h3(title),if(!is.null(subtitle))p(subtitle)))
empty <- function(title,text) div(class='empty-state',div(mark(),class='empty-mark'),h3(title),p(text))
hero_art <- HTML('<svg class="hero-art" viewBox="0 0 470 190" role="img" aria-label="Stylized barcode lineages flowing through time"><defs><clipPath id="flow-clip"><rect x="0" y="0" width="470" height="190" rx="18"/></clipPath></defs><g clip-path="url(#flow-clip)"><path d="M0 0H470V38C350 34 240 80 0 28Z" fill="#76978a"/><path d="M0 28C240 80 350 34 470 38V76C310 75 250 103 0 60Z" fill="#a6b8a3"/><path d="M0 60C250 103 310 75 470 76V111C350 119 200 122 0 85Z" fill="#d8bd82"/><path d="M0 85C200 122 350 119 470 111V138C320 157 200 137 0 123Z" fill="#bc8f69"/><path d="M0 123C200 137 320 157 470 138V167C330 184 210 163 0 154Z" fill="#d5a494"/><path d="M0 154C210 163 330 184 470 167V190H0Z" fill="#4a7770"/></g><g stroke="#fff" stroke-opacity=".2" stroke-dasharray="3 5"><path d="M70 0v190M180 0v190M290 0v190M400 0v190"/></g></svg>')

ui <- page_fluid(
  title='barbac Studio · From sequences to lineages',
  theme=bs_theme(version=5,bg='#f6f5f0',fg='#193f3b',primary='#23685f',
                 base_font=font_collection('Avenir Next','Avenir','Segoe UI','sans-serif')),
  tags$head(tags$link(rel='stylesheet',href='studio.css'),tags$script(src='studio.js')),
  div(class='studio-shell',
    tags$aside(class='sidebar',
      div(class='brand',div(mark(),class='brand-mark'),div(span('barbac',class='wordmark'),span('STUDIO',class='brand-sub'))),
      div(class='sidebar-caption','A workspace for barcode discovery'),
      div(class='nav-label','YOUR ANALYSIS'),
      radioButtons('page',NULL,choices=c('01   Data workspace'='data','02   Cluster sequences'='cluster',
        '03   Explore results'='results','04   Downloads'='downloads'),selected='data'),
      div(class='sidebar-bottom',div(class='engine-status',span(class='status-dot'),'Native barbac engine'),
        p('From extracted sequences to the lineages that matter.'),
        tags$a(href='https://loukesio.github.io/barbac/',target='_blank',rel='noopener','Package documentation ',arrow()),
        div(class='sidebar-version',paste('v',packageVersion('barbac'),' · Local workspace',sep='')))),
    tags$main(class='main-content',
      div(class='topbar',div(span('WORKSPACE',class='topbar-label'),span('/',class='topbar-slash'),textOutput('page_title',inline=TRUE)),
        div(class='topbar-right',uiOutput('dataset_badge',inline=TRUE),actionButton('reset','New analysis',class='button-quiet'))),
      uiOutput('job_banner'),
      conditionalPanel("input.page === 'data'",
        div(class='hero',div(class='hero-copy',eyebrow('BARCODE LINEAGE ANALYSIS'),
          h1('Every barcode.',tags$br(),span('A clearer story.',class='hero-italic')),
          p('Bring your sequences. Discover your lineages. Move from barcode counts to an analysis you can explore and share.')),
          div(class='hero-visual',hero_art,div(class='hero-caption',span('SEQUENCES'),span('LINEAGES'),span('DISCOVERY')))),
        div(class='section-top',div(h2('Start with your data'),p('Use extracted barcode counts, or begin with raw sequencing reads.')),
          div(class='section-tag','01 / INPUT')),
        radioButtons('source_type',NULL,c('Extracted barcodes'='counts','Raw FASTQ reads'='fastq'),inline=TRUE),
        conditionalPanel("input.source_type === 'counts'",
          div(class='input-grid',
            div(class='studio-card upload-card',panel_heading('A','Upload barcode counts','CSV or TSV · one or several sample files'),
              fileInput('counts_files',NULL,multiple=TRUE,accept=c('.csv','.tsv','.csv.gz','.tsv.gz'),buttonLabel='Choose files',placeholder='Drop your barcode-count files here'),
              div(class='schema-hint',span('Required'),tags$code('barcode'),tags$code('counts'),span('Optional'),tags$code('sample'),tags$code('time'),tags$code('population')),
              tags$details(tags$summary('Add sample metadata'),p('A CSV with sample, time and population. Population labels define independent libraries or replicate series.'),
                fileInput('metadata',NULL,accept='.csv',buttonLabel='Choose metadata',placeholder='Optional sample metadata.csv')),
              uiOutput('import_message')),
            div(class='studio-card example-card',eyebrow('TAKE A LOOK AROUND'),h3('A small dataset.',tags$br(),'A complete story.'),
              p('Explore a synthetic population over eight timepoints, with known barcode variants.'),
              actionButton('load_demo',tagList('Load example',arrow()),class='button-secondary'),
              downloadLink('demo_csv','Download example CSV',class='text-link'),
              downloadLink('metadata_csv','Download metadata template',class='text-link')))),
        conditionalPanel("input.source_type === 'fastq'",
          div(class='studio-card raw-card',panel_heading('A','Extract from sequencing reads','One barcode locus per run · single-end or overlapping paired-end reads'),
            div(class='raw-grid',fileInput('raw_r1','R1 or single-end FASTQ',multiple=TRUE,accept=c('.fastq','.fq','.fastq.gz','.fq.gz')),
              fileInput('raw_r2','R2 FASTQ · optional',multiple=TRUE,accept=c('.fastq','.fq','.fastq.gz','.fq.gz')),
              fileInput('reference','Barcode-cassette reference',accept=c('.fasta','.fa','.fna'))),
            p(class='field-note','For multiple paired samples, select R1 and R2 files in matching order. The sample name comes from the R1 filename. Paired reads must overlap; unmerged pairs are excluded.'),
            div(class='raw-grid',selectInput('extract_mode','Extraction method',c('Exact flanks · preserve indels'='flanks','Fixed reference coordinates'='fixed')),
              numericInput('start','Barcode start · 1-based',171,min=2,step=1),numericInput('end','Barcode end · inclusive',196,min=3,step=1)),
            conditionalPanel("input.extract_mode === 'flanks'",div(class='raw-grid',textInput('left_flank','Left flank · exact DNA'),textInput('right_flank','Right flank · exact DNA'),
              div(class='two-fields',numericInput('min_length','Min. length',24,min=1,max=250),numericInput('max_length','Max. length',28,min=1,max=250)))),
            conditionalPanel("input.extract_mode === 'fixed'",p(class='field-note','Fixed-coordinate extraction follows reference positions and does not preserve observed insertion/deletion lengths. Use exact flanks when length variation matters.')),
            fileInput('raw_metadata','Sample metadata · optional',accept='.csv',placeholder='sample, time, population'),
            div(class='raw-actions',actionButton('extract','Extract barcodes',class='button-primary job-action'),
              downloadLink('demo_fastq','Download synthetic FASTQ example',class='text-link'),
              actionButton('demo_flanks','Use example coordinates & flanks',class='button-quiet')),
            p(class='field-note','Supply the reference and settings for your own construct. This route applies mapping and barcode extraction; it does not deduplicate UMIs or reproduce study-specific quality filters.'))),
        div(class='section-top compact',div(h2('Your dataset'),uiOutput('input_caption')),uiOutput('input_actions')),
        uiOutput('input_metrics'),
        div(class='studio-card table-card',DT::DTOutput('input_table'))),
      conditionalPanel("input.page === 'cluster'",
        div(class='page-heading',eyebrow('02 / CLUSTER SEQUENCES'),h1('Find the shared signal.'),p('Choose how barcode variants become lineages. Your original counts stay intact.')),
        div(class='settings-grid',div(class='studio-card',panel_heading('01','Distance model'),
          radioButtons('method',NULL,c('Levenshtein · substitutions and indels'='lv','Hamming · fixed-length substitutions'='hamming')),
          p(class='field-note','LV handles insertions, deletions and positional shifts. Use Hamming only for fixed-length A/C/G/T barcodes up to 32 bases when indels are excluded.'),
          numericInput('distance','Maximum edit distance',3,min=1,max=5,step=1),
          tags$details(tags$summary('Advanced clustering settings'),
            numericInput('merge_ratio','Count-ratio guard',20,min=1,max=1000),
            numericInput('error_rate','Per-base error rate',.005,min=.000001,max=.249,step=.001),
            selectInput('tie_break','Equal-count ordering',c('Sequence · package default'='sequence','Neighbour support'='support')),
            conditionalPanel("input.method === 'lv'",checkboxInput('poisson','Experimental Poisson indel option',FALSE),
              p(class='field-note','Opt-in merge-guard exception for repeated-base single indels. May merge genuine length variants; validate against independent controls.')))),
          div(class='studio-card run-card',eyebrow('READY WHEN YOU ARE'),h3('One shared identity.',tags$br(),'Across every timepoint.'),
            p('Within each population, we pool counts once, cluster the sequences, and map every sample back to those memberships. Independent populations stay separate.'),
            uiOutput('run_summary'),actionButton('run','Run clustering',class='button-primary job-action'),
            p(class='field-note','A background worker keeps this workspace responsive. Large, diverse libraries can take minutes.')))),
      conditionalPanel("input.page === 'results'",
        div(class='page-heading results-heading',div(eyebrow('03 / EXPLORE RESULTS'),h1('Your lineages, revealed.'),uiOutput('result_caption')),
          actionButton('go_downloads',tagList('Export analysis',arrow()),class='button-secondary')),
        uiOutput('result_metrics'),
        conditionalPanel('output.has_results',
          div(class='studio-card plot-card',div(class='plot-heading',div(h3('Lineages through time'),p('Every inferred lineage. All barcode counts.')),
            div(class='plot-controls',selectInput('population',NULL,choices=NULL),selectInput('palette',NULL,choices=names(barbac::barbac_palettes()),selected='alger'))),
            uiOutput('trajectory'),div(class='plot-footnote','Frequencies use all barcode counts per sample. Missing counts are zero. No lineages are grouped into a remainder.')),
          div(class='results-lower',div(class='studio-card',h3('Abundance profile'),p(class='field-note','All clusters in the selected population, ranked by barcode count.'),plotOutput('rank_plot',height='255px')),
            div(class='studio-card',h3('Clustering at a glance'),uiOutput('stats_detail'),p(class='field-note','These are descriptive summaries of inferred clusters, not accuracy estimates.'))),
          div(class='studio-card table-card',div(class='table-heading',div(h3('Explore the clusters'),p('Search a centroid sequence or cluster identifier.')),
              radioButtons('table_view',NULL,c('Centroids'='centroids','Memberships'='memberships','Sample counts'='time_series'),inline=TRUE)),DT::DTOutput('result_table'))),
        conditionalPanel('!output.has_results',empty('Your results will appear here','Load barcode counts and run clustering to explore your analysis.'))),
      conditionalPanel("input.page === 'downloads'",
        div(class='page-heading',eyebrow('04 / DOWNLOADS'),h1('Take the whole story.'),p('Keep your data, inspect every assignment, and share a reproducible analysis.')),
        div(class='download-grid',
          div(class='studio-card download-card',span('CSV',class='file-type'),h3('Extracted barcode counts'),p('The complete input table, ready to cluster again in Studio or R.'),downloadButton('download_input','Download counts',class='button-secondary')),
          div(class='studio-card download-card',span('ZIP',class='file-type'),h3('Complete analysis'),p('Centroids, memberships, sample counts, statistics, settings and the saved R result.'),downloadButton('download_all','Download analysis',class='button-primary')),
          div(class='studio-card download-card',span('HTML',class='file-type'),h3('An analysis to share'),p('A self-contained Quarto report with settings, cluster summaries and lineage figures.'),actionButton('build_report','Prepare report',class='button-secondary'),uiOutput('report_download_ui')),
          div(class='studio-card download-card',span('PDF / CSV',class='file-type'),h3('Figures & individual tables'),p('Export the current population with your chosen palette, or download complete tables.'),
            downloadLink('download_plot','Lineage figure · PDF',class='text-link'),downloadLink('download_centroids','Centroids · CSV',class='text-link'),
            downloadLink('download_memberships','Memberships · CSV',class='text-link'),downloadLink('download_series','Sample counts · CSV',class='text-link'))),
        div(class='studio-card provenance-card',h3('A record of your analysis'),uiOutput('provenance'),
          p(class='field-note','Your workspace is session-specific. Download results before closing this browser tab; temporary uploads and outputs are removed when the session ends.'))),
      tags$footer(class='footer',span('barbac Studio'),span('Built for careful science. Designed for discovery.')))))

server <- function(input,output,session) {
  session_dir <- tempfile('barbac-studio-');dir.create(session_dir)
  data <- reactiveVal(studio_validate(studio_demo()))
  label <- reactiveVal('Synthetic example · 8 timepoints')
  synthetic <- reactiveVal(TRUE)
  result <- reactiveVal(NULL)
  extraction <- reactiveVal(NULL)
  import_note <- reactiveVal(NULL)
  report_path <- reactiveVal(NULL)
  progress_dir <- reactiveVal(NULL)
  worker <- ExtendedTask$new(function(job,directory) {
    engine <- file.path(app_dir,'R','engine.R');report_code <- file.path(app_dir,'R','report.R')
    template <- file.path(app_dir,'report.qmd')
    promises::future_promise(local({
      on.exit(if(file.exists(file.path(dirname(directory),'session-closed')))
        unlink(dirname(directory),recursive=TRUE),add=TRUE)
      source(engine,local=TRUE)
      if(job$kind=='cluster') list(kind='cluster',value=studio_cluster(job$data,job$settings,directory,job$provenance))
      else if(job$kind=='extract') list(kind='extract',value=studio_extract(job,directory))
      else {
        source(report_code,local=TRUE)
        target <- file.path(directory,'report.html')
        studio_report(job$result,target,template,job$palette)
        list(kind='report',value=target)
      }
    }),globals=list(job=job,directory=directory,engine=engine,report_code=report_code,template=template),seed=TRUE)
  })
  # Every operation shares this worker; the UI prevents competing submissions.
  busy <- reactive(identical(worker$status(),'running'))
  can_plot <- reactive({
    r<-result();if(is.null(r)||is.null(input$population))return(FALSE)
    ts<-r$time_series[r$time_series$population==input$population,]
    nrow(ts)>0 && !anyNA(ts$time) && length(unique(ts$time))>=2 &&
      length(unique(ts$cluster_id))*length(unique(ts$time))<=2000000
  })
  session$onSessionEnded(function(){
    # The worker removes its own files on exit if the browser has disconnected.
    writeLines('closed',file.path(session_dir,'session-closed'))
    if(!isTRUE(isolate(busy()))) unlink(session_dir,recursive=TRUE)
  })
  set_page <- function(page) updateRadioButtons(session,'page',selected=page)
  reset_result <- function() {result(NULL);report_path(NULL)}
  notify_error <- function(e) showNotification(conditionMessage(e),type='error',duration=12)
  begin <- function(job) {
    if(busy()) {showNotification('An analysis is already running.',type='message');return(invisible(NULL))}
    directory <- tempfile('run-',tmpdir=session_dir);dir.create(directory)
    studio_progress(directory,'Starting background worker')
    progress_dir(directory);worker$invoke(job,directory)
  }
  observe({session$sendCustomMessage('studio-busy',busy())})
  observe({session$sendCustomMessage('studio-download-state',list(input=!is.null(data()),results=!is.null(result()),plot=can_plot()))})
  observeEvent(worker$status(),{
    if(worker$status()=='success') {
      done <- worker$result()
      if(done$kind=='cluster') {
        result(done$value);report_path(NULL)
        updateSelectInput(session,'population',choices=unique(done$value$centroids$population))
        set_page('results')
      } else if(done$kind=='extract') {
        data(done$value$input);extraction(done$value);synthetic(FALSE)
        label('Extracted from your FASTQ reads');reset_result();set_page('cluster')
      } else report_path(done$value)
      showNotification(if(done$kind=='report')'Your report is ready to download.' else 'Analysis complete.',type='message',duration=4)
    }
    if(worker$status()=='error') tryCatch(worker$result(),error=notify_error)
  },ignoreInit=TRUE)
  output$page_title <- renderText(switch(input$page,data='Data workspace',cluster='Cluster sequences',results='Explore results',downloads='Downloads'))
  output$dataset_badge <- renderUI(span(class=if(synthetic())'badge-example' else 'badge-live',if(synthetic())'EXAMPLE DATA' else 'YOUR DATA'))
  output$job_banner <- renderUI({
    if(!busy())return(NULL)
    invalidateLater(800)
    path <- file.path(progress_dir(),'progress.txt')
    msg <- if(file.exists(path))paste(readLines(path,warn=FALSE),collapse=' ') else 'Working'
    div(class='job-banner',span(class='spinner'),strong(msg),span('You can keep exploring this workspace.'))
  })
  output$import_message <- renderUI(if(!is.null(import_note()))div(class='inline-message',import_note()))
  observeEvent(list(input$counts_files,input$metadata),{
    if(is.null(input$counts_files)||busy())return()
    tryCatch({
      x <- studio_import(input$counts_files,if(!is.null(input$metadata))input$metadata$datapath else NULL)
      data(x);synthetic(FALSE);label(paste(nrow(input$counts_files),'uploaded file(s)'));extraction(NULL);reset_result()
      import_note('Files checked. Exact duplicate barcodes are summed within each sample.')
    },error=function(e){data(NULL);synthetic(FALSE);label('Input needs attention');reset_result();import_note(conditionMessage(e));notify_error(e)})
  },ignoreInit=TRUE)
  observeEvent(input$load_demo,{if(busy())return();data(studio_validate(studio_demo()));synthetic(TRUE);label('Synthetic example · 8 timepoints');extraction(NULL);reset_result();import_note(NULL)})
  observeEvent(input$reset,{if(busy())return();data(NULL);synthetic(FALSE);label('No data loaded');extraction(NULL);reset_result();import_note(NULL);set_page('data');session$sendCustomMessage('studio-clear-files',TRUE)})
  observeEvent(input$continue,set_page('cluster'))
  observeEvent(input$go_downloads,set_page('downloads'))
  output$input_caption <- renderUI(p(class='field-note',label()))
  output$input_actions <- renderUI(if(!is.null(data()))actionButton('continue',tagList('Continue to clustering',arrow()),class='button-primary'))
  metric <- function(title,value,note) div(class='metric-card',span(title,class='metric-label'),strong(value,class='metric-value'),span(note,class='metric-note'))
  fmt <- function(x)format(x,big.mark=',',scientific=FALSE,trim=TRUE)
  output$input_metrics <- renderUI({
    x <- data();if(is.null(x))return(empty('A fresh workspace','Upload barcode counts or load the example to get started.'))
    div(class='metric-grid',metric('BARCODE READS',fmt(sum(x$counts)),'Total supplied counts'),
      metric('UNIQUE SEQUENCES',fmt(length(unique(x$barcode))),'Before error correction'),
      metric('SAMPLES',fmt(length(unique(x$sample))),'Original identities retained'),
      metric('BARCODE LENGTH',paste(range(nchar(x$barcode)),collapse='–'),'Observed bases'))
  })
  table_options <- list(pageLength=8,lengthChange=FALSE,scrollX=TRUE,dom='ftp',autoWidth=TRUE,
    language=list(search='',searchPlaceholder='Search sequences or samples…',emptyTable='No data loaded'))
  output$input_table <- DT::renderDT({req(data());DT::datatable(data(),rownames=FALSE,options=table_options,escape=TRUE)},server=TRUE)
  output$run_summary <- renderUI({x<-data();if(is.null(x))return(p('Load data first.'));div(class='run-summary',
    div(strong(fmt(length(unique(x$barcode)))),span('unique sequences')),div(strong(fmt(length(unique(x$population)))),span('independent populations')),
    div(strong(fmt(sum(x$counts))),span('barcode reads')))})
  observeEvent(input$run,{
    if(busy())return()
    tryCatch({
      if(is.null(data()))studio_error('Load a barcode-count table before clustering.')
      settings <- studio_settings(input$method,input$distance,input$merge_ratio,input$error_rate,input$tie_break,
        if(input$method=='lv'&&isTRUE(input$poisson))'poisson' else 'none')
      begin(list(kind='cluster',data=data(),settings=settings,provenance=list(synthetic=synthetic(),
        extraction=if(!is.null(extraction()))extraction()$provenance else NULL,
        extraction_stats=if(!is.null(extraction()))extraction()$extraction_stats else NULL)))
    },error=notify_error)
  })
  observeEvent(input$extract,{
    if(busy())return()
    tryCatch({
      if(is.null(input$raw_r1)||is.null(input$reference))studio_error('Choose R1/single-end reads and a reference FASTA.')
      files <- c(input$raw_r1$datapath,input$raw_r2$datapath,input$reference$datapath)
      if(sum(file.info(files)$size)>as.numeric(Sys.getenv('BARBAC_STUDIO_UPLOAD_MB','512'))*1024^2)
        studio_error('Combined raw-read uploads exceed the configured size limit.')
      # Copy uploads before starting: a later browser upload must not replace worker inputs.
      staging <- tempfile('uploads-',tmpdir=session_dir);dir.create(staging)
      stage_files <- function(z,prefix) {if(is.null(z))return(NULL);for(i in seq_len(nrow(z))){p<-file.path(staging,paste0(prefix,i));file.copy(z$datapath[i],p);z$datapath[i]<-p};z}
      r1<-stage_files(input$raw_r1,'r1-');r2<-stage_files(input$raw_r2,'r2-')
      ref<-file.path(staging,'reference.fasta');file.copy(input$reference$datapath,ref)
      metadata<-NULL;if(!is.null(input$raw_metadata)){metadata<-file.path(staging,'metadata.csv');file.copy(input$raw_metadata$datapath,metadata)}
      begin(list(kind='extract',r1=r1,r2=r2,reference=ref,metadata=metadata,
        settings=list(mode=input$extract_mode,start=input$start,end=input$end,left=input$left_flank,right=input$right_flank,min_length=input$min_length,max_length=input$max_length)))
    },error=notify_error)
  })
  output$has_results <- reactive(!is.null(result()))
  outputOptions(output,'has_results',suspendWhenHidden=FALSE)
  output$result_caption <- renderUI({r<-result();if(is.null(r))return(p('Run an analysis to begin.'));p(
    if(isTRUE(r$provenance$synthetic))'Synthetic example · ' else '',
    if(r$settings$method=='lv')'Levenshtein' else 'Hamming',' · distance ',r$settings$distance,' · ',sprintf('%.2f',r$seconds),' s clustering')})
  output$result_metrics <- renderUI({r<-result();if(is.null(r))return(NULL);div(class='metric-grid',
    metric('INFERRED CLUSTERS',fmt(nrow(r$centroids)),'Across all populations'),
    metric('BARCODE READS',fmt(sum(r$input$counts)),'All counts conserved'),
    metric('SEQUENCES ABSORBED',fmt(nrow(r$memberships)-nrow(r$centroids)),'Non-centroid sequences'),
    metric('CLUSTERING TIME',paste0(sprintf('%.2f',r$seconds),' s'),'Excludes upload & extraction'))})
  output$trajectory <- renderUI({
    req(result(),input$population)
    ts<-result()$time_series;ts<-ts[ts$population==input$population,]
    if(anyNA(ts$time)||length(unique(ts$time))<2)return(empty('A timeline needs timepoints','Provide a time column, or sample metadata with at least two timepoints in each population.'))
    if(length(unique(ts$cluster_id))*length(unique(ts$time))>2000000)return(empty('A larger story','Download the complete count table to plot this large lineage grid in R.'))
    if(length(unique(ts$cluster_id))>1500)plotOutput('area_static',height='365px') else ggiraph::girafeOutput('area_interactive',height='365px')
  })
  output$area_interactive <- ggiraph::renderGirafe({
    req(can_plot());width<-session$clientData$output_area_interactive_width
    studio_area(result(),input$population,input$palette,'ggiraph',compact=!is.null(width)&&width<600)
  })
  output$area_static <- renderPlot({req(can_plot());studio_area(result(),input$population,input$palette)},res=110)
  output$rank_plot <- renderPlot({
    req(result(),input$population)
    x<-result()$centroids;x<-x[x$population==input$population,];x<-x[order(-x$sum_counts),];x$rank<-seq_len(nrow(x))
    ggplot2::ggplot(x,ggplot2::aes(rank,sum_counts))+ggplot2::geom_line(colour='#2c756a',linewidth=.8)+
      ggplot2::geom_point(colour='#b49a64',size=1.7)+ggplot2::scale_y_log10(labels=scales::label_number())+
      ggplot2::labs(x='Cluster rank',y='Barcode reads · log scale')+ggplot2::theme_minimal(base_size=11)+
      ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),plot.background=ggplot2::element_rect(fill='white',colour=NA))
  },res=110)
  output$stats_detail <- renderUI({req(result(),input$population);s<-result()$stats;s<-s[s$population==input$population,];
    div(class='stats-list',div(span('Clusters'),strong(fmt(s$n_clusters))),div(span('Singleton clusters'),strong(fmt(s$n_singletons))),
      div(span('Largest cluster share'),strong(sprintf('%.2f%%',100*s$largest_cluster_frac))),div(span('Top 10 cluster share'),strong(sprintf('%.2f%%',100*s$top10_frac))),
      div(span('Total barcode reads'),strong(fmt(s$total_reads))))})
  output$result_table <- DT::renderDT({req(result());x<-result()[[input$table_view]];if(!is.null(input$population))x<-x[x$population==input$population,,drop=FALSE];
    DT::datatable(x,rownames=FALSE,options=table_options,escape=TRUE)},server=TRUE)
  output$provenance <- renderUI({r<-result();if(is.null(r))return(p('Complete a clustering run to record its settings and engine version.'));
    div(class='provenance-grid',div(span('ENGINE'),tags$code(r$provenance$build_id)),div(span('COMPLETED'),strong(r$provenance$created_utc)),
      div(span('METHOD'),strong(paste(r$settings$method,'· distance',r$settings$distance,'·',r$settings$tie_break,'ordering'))),
      div(span('INPUT FINGERPRINT'),tags$code(r$provenance$input_sha256)))})
  output$demo_csv <- downloadHandler('barbac-example.csv',function(file)readr::write_csv(studio_demo(),file))
  output$metadata_csv <- downloadHandler('sample-metadata.csv',function(file)readr::write_csv(unique(studio_demo()[c('sample','time','population')]),file))
  output$demo_fastq <- downloadHandler('barbac-fastq-example.zip',function(file){d<-tempfile('fastq-example-');dir.create(d);on.exit(unlink(d,recursive=TRUE));studio_fastq_example(d);zip::zipr(file,list.files(d,full.names=TRUE),root=d)})
  observeEvent(input$demo_flanks,{d<-tempfile();dir.create(d);cfg<-studio_fastq_example(d);unlink(d,recursive=TRUE);
    updateNumericInput(session,'start',value=cfg$start);updateNumericInput(session,'end',value=cfg$end);updateNumericInput(session,'min_length',value=cfg$min_length);updateNumericInput(session,'max_length',value=cfg$max_length);
    updateTextInput(session,'left_flank',value=cfg$left);updateTextInput(session,'right_flank',value=cfg$right);updateSelectInput(session,'extract_mode',selected='flanks')})
  output$download_input <- downloadHandler('extracted_barcodes.csv',function(file){req(data());readr::write_csv(data(),file)})
  output$download_all <- downloadHandler('barbac-analysis.zip',function(file){req(result());r<-result();files<-list.files(r$directory,pattern='[.](csv|json|rds)$',full.names=TRUE);if(!is.null(report_path()))files<-c(files,report_path());zip::zipr(file,files,mode='cherry-pick')})
  output$download_centroids <- downloadHandler('centroids.csv',function(file){req(result());readr::write_csv(result()$centroids,file)})
  output$download_memberships <- downloadHandler('memberships.csv',function(file){req(result());readr::write_csv(result()$memberships,file)})
  output$download_series <- downloadHandler('time_series.csv',function(file){req(result());readr::write_csv(result()$time_series,file)})
  output$download_plot <- downloadHandler('lineage_trajectories.pdf',function(file){req(can_plot());p<-studio_area(result(),input$population,input$palette);ggplot2::ggsave(file,p,width=10,height=5,device=grDevices::pdf)})
  observeEvent(input$build_report,{if(busy())return();if(is.null(result())){showNotification('Run clustering before preparing a report.',type='message');return()};begin(list(kind='report',result=result(),palette=input$palette))})
  output$report_download_ui <- renderUI(if(!is.null(report_path()))downloadButton('download_report','Download HTML report',class='button-primary'))
  output$download_report <- downloadHandler('barbac-report.html',function(file){req(report_path());file.copy(report_path(),file)})
}

shinyApp(ui,server)
