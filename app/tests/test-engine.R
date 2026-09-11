library(testthat)
app <- if(file.exists('app/R/engine.R')) normalizePath('app') else normalizePath('..')
.libPaths(c(file.path(app,'.runtime','library'),.libPaths()))
source(file.path(app,'R','engine.R'))
source(file.path(app,'R','examples.R'))
source(file.path(app,'R','report.R'))

test_that('imports reject malformed counts and ambiguous sample identities',{
  x <- data.frame(barcode=c('acgt','ACGT'),counts=c(4,6))
  expect_equal(studio_validate(x)$counts,10)
  expect_equal(studio_validate(x)$barcode,'ACGT')
  for(bad in list(NA,-1,0,Inf,1.2,'text',2147483648)) {
    x$counts<-c(4,bad)
    expect_error(studio_validate(x))
  }
  expect_error(studio_validate(data.frame(barcode='AC;GT',counts=1)),'A, C, G')
  x<-studio_demo();x$time[1]<-999
  expect_error(studio_validate(x),'exactly one')
  x<-studio_demo();x$time[x$time==8]<-0
  expect_error(studio_validate(x),'share a timepoint')
  x<-studio_demo();x$time[x$time==8]<-NA
  expect_error(studio_validate(x),'every sample')
})

test_that('CSV and TSV imports preserve counts and require complete metadata',{
  folder<-tempfile();dir.create(folder);on.exit(unlink(folder,recursive=TRUE))
  x<-data.frame(barcode=c('ACGT','TGCA'),counts=c(4,6))
  readr::write_csv(x,file.path(folder,'a.csv'));readr::write_tsv(x,file.path(folder,'b.tsv'))
  files<-data.frame(name=c('a.csv','b.tsv'),datapath=file.path(folder,c('a.csv','b.tsv')))
  got<-studio_import(files)
  expect_equal(sum(got$counts),20)
  expect_setequal(got$sample,c('a','b'))
  metadata<-file.path(folder,'metadata.csv')
  readr::write_csv(data.frame(sample=c('a','b'),time=c(0,1),population='p'),metadata)
  expect_equal(sort(unique(studio_import(files,metadata)$time)),c(0,1))
  readr::write_csv(data.frame(sample='a',time=0,population='p'),metadata)
  expect_error(studio_import(files,metadata),'missing samples')
  writeLines(c('barcode,counts','ACGT,1,2'),file.path(folder,'a.csv'))
  expect_error(suppressWarnings(studio_import(files)),'wrong number')
})

test_that('the app exactly preserves native memberships and every sample count',{
  x<-studio_validate(studio_demo());folder<-tempfile()
  on.exit(unlink(folder,recursive=TRUE))
  settings<-studio_settings()
  r<-studio_cluster(x,settings,folder,list(synthetic=TRUE))
  pooled<-aggregate(counts~barcode,x,sum)
  native<-do.call(barbac::super_cluster2,c(list(input_path=pooled,verbose=FALSE),settings))
  expect_identical(r$centroids$central_barcode,native$central_barcode)
  expect_equal(r$centroids$sum_counts,native$sum_counts)
  expect_identical(r$memberships$barcode,unlist(native$all_barcodes,use.names=FALSE))
  expect_equal(nrow(r$centroids),18)
  expect_equal(sum(r$centroids$sum_counts),192000)
  expect_equal(tapply(r$time_series$counts,r$time_series$sample,sum),tapply(x$counts,x$sample,sum))
  expect_equal(as.numeric(tapply(r$time_series$frequency,r$time_series$sample,sum)),rep(1,8))
  expect_true(isTRUE(r$provenance$synthetic))
  expect_true(all(file.exists(file.path(folder,c('centroids.csv','memberships.csv','time_series.csv','analysis.json','analysis.rds')))))
  expect_equal(sum(read.csv(file.path(folder,'time_series.csv'))$counts),sum(x$counts))
  expect_s3_class(studio_area(r,'Example population'),'ggplot')
  expect_s3_class(studio_area(r,'Example population',interactive='ggiraph'),'girafe')
  native_plot<-barbac::barbac_ts_area(r$time_series,id_col='cluster_id',
    min_total_count=0,fill_missing='zero',palette='alger')
  app_plot<-studio_area(r,'Example population')
  cols<-c('x','y','ymin','ymax','group','fill')
  expect_equal(ggplot2::ggplot_build(app_plot)$data[[1]][cols],
               ggplot2::ggplot_build(native_plot)$data[[1]][cols])
  expect_equal(app_plot$scales$get_scales('x')$breaks,seq(0,56,8))
})

test_that('independent populations stay separate and row order does not change clusters',{
  x<-studio_demo();y<-x;y$sample<-paste0('rep2-',y$sample);y$population<-'Replicate 2'
  x<-studio_validate(rbind(x,y));folder<-tempfile();other<-tempfile()
  on.exit(unlink(c(folder,other),recursive=TRUE))
  a<-studio_cluster(x,studio_settings(tie_break='support'),folder)
  b<-studio_cluster(x[nrow(x):1,],studio_settings(tie_break='support'),other)
  expect_equal(nrow(a$centroids),36)
  expect_equal(sort(a$centroids$sum_counts),sort(b$centroids$sum_counts))
  expect_equal(length(unique(a$centroids$cluster_id)),36)
  expect_setequal(a$centroids$population,c('Example population','Replicate 2'))
  expect_equal(sum(a$time_series$counts),384000)
})

test_that('Hamming is constrained to compatible libraries and options are validated',{
  x<-studio_demo();x$barcode[1]<-paste0(x$barcode[1],'A')
  expect_error(studio_cluster(x,studio_settings(method='hamming'),tempfile()),'Use LV')
  expect_error(studio_settings(distance=2.1),'whole number')
  expect_error(studio_settings(error_rate=0),'Error rate')
  expect_error(studio_settings(method='hamming',indel_model='poisson'),'requires LV')
})

test_that('single-timepoint data cluster without inventing a timeline',{
  x<-studio_demo();x<-x[x$time==0,c('barcode','counts')]
  folder<-tempfile();on.exit(unlink(folder,recursive=TRUE))
  r<-studio_cluster(x,studio_settings(),folder)
  expect_equal(nrow(r$centroids),18)
  expect_error(studio_area(r,'Library 1'),'at least two')
})

test_that('real CLI extraction recovers known single-end and paired-end barcodes',{
  fixture<-tempfile();cfg<-studio_fastq_example(fixture)
  folder<-tempfile();dir.create(folder)
  on.exit(unlink(c(fixture,folder),recursive=TRUE))
  upload<-function(name)data.frame(name=name,datapath=file.path(fixture,name))
  base<-list(kind='extract',reference=file.path(fixture,'reference.fasta'),settings=cfg,metadata=NULL)
  truth<-read.csv(file.path(fixture,'expected_barcodes.csv'))
  for(paired in c(FALSE,TRUE)) {
    job<-base;job$r1<-upload(if(paired)'example_R1.fastq.gz' else 'example_single.fastq.gz')
    job$r2<-if(paired)upload('example_R2.fastq.gz') else NULL
    r<-studio_extract(job,file.path(folder,if(paired)'paired' else 'single'))
    expect_equal(r$input$counts[match(truth$barcode,r$input$barcode)],truth$counts)
    expect_equal(r$extraction_stats$input_reads_or_pairs,140)
    expect_equal(r$extraction_stats$mapped_primary,140)
    expect_equal(r$extraction_stats$extracted_counts,140)
  }
  job$settings$mode<-'fixed'
  r<-studio_extract(job,file.path(folder,'fixed'))
  expect_equal(r$input$counts[match(truth$barcode,r$input$barcode)],truth$counts)
  bad<-file.path(fixture,'bad.fastq');writeLines(c('@bad','ACGT','+','III'),bad)
  expect_error(studio_check_fastq(bad),'Invalid FASTQ')
})

test_that('the Quarto export is self-contained and includes actual analysis settings',{
  folder<-tempfile();dir.create(folder);on.exit(unlink(folder,recursive=TRUE))
  r<-studio_cluster(studio_demo(),studio_settings(),file.path(folder,'analysis'),list(synthetic=TRUE))
  target<-file.path(folder,'report.html')
  studio_report(r,target,file.path(app,'report.qmd'))
  html<-paste(readLines(target,warn=FALSE),collapse='\n')
  expect_match(html,'192,000',fixed=TRUE)
  expect_match(html,'barbac-2026-09-10-eligible-parent-v14',fixed=TRUE)
  expect_match(html,'data:image/png;base64',fixed=TRUE)
  expect_match(html,'synthetic demonstration',fixed=TRUE)
})
