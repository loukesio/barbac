#!/usr/bin/env Rscript
# R drives Chrome for UI checks; JavaScript inspects rendered elements only.
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
here <- dirname(normalizePath(script));dir.create(file.path(here,'qa'),showWarnings=FALSE)
source(file.path(here,'../provenance_helpers.R'))
if('--wait' %in% commandArgs(TRUE)) {
  ready <- function() {
    tryCatch({
      v <- jsonlite::read_json(file.path(here,'validation.json'))
      assert_hash_map(v$hashes,here)
    },error=function(e)FALSE)
  }
  deadline <- Sys.time()+7200
  while(!ready()) {if(Sys.time()>deadline)stop('Report build deadline exceeded');Sys.sleep(10)}
}
assert_hash_map(jsonlite::read_json(file.path(here,'validation.json'))$hashes,here)
endpoints <- jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory('http://127.0.0.1:9223/json/list')$content))
ws <- websocket::WebSocket$new(endpoints$webSocketDebuggerUrl[endpoints$type=='page'][1],maxMessageSize=64*1024^2)
responses <- new.env(parent=emptyenv());messages <- list();id <- 0L
ws$onMessage(function(event) {
  value <- jsonlite::fromJSON(event$data,simplifyVector=FALSE)
  if(!is.null(value$id)) responses[[as.character(value$id)]] <- value else messages[[length(messages)+1L]] <<- value
})
wait_for <- function(predicate,timeout=45) {
  start <- Sys.time()
  while(!predicate()) {
    if(as.numeric(difftime(Sys.time(),start,units='secs'))>timeout) stop('Browser QA timed out')
    later::run_now(.05)
  }
}
wait_for(function()ws$readyState()==1)
rpc <- function(method,params=list()) {
  id <<- id+1L;key <- as.character(id);command <- list(id=id,method=method)
  if(length(params)) command$params <- params
  ws$send(jsonlite::toJSON(command,auto_unbox=TRUE,null='null'))
  wait_for(function()exists(key,responses,inherits=FALSE))
  result <- responses[[key]]
  if(!is.null(result$error))stop(jsonlite::toJSON(result$error,auto_unbox=TRUE))
  result$result
}
js <- function(expression) {
  value <- rpc('Runtime.evaluate',list(expression=expression,returnByValue=TRUE,awaitPromise=TRUE))
  if(!is.null(value$exceptionDetails))stop(jsonlite::toJSON(value$exceptionDetails,auto_unbox=TRUE))
  value$result$value
}
screenshot <- function(name) {
  js('new Promise(resolve=>setTimeout(resolve,600))')
  result <- rpc('Page.captureScreenshot',list(format='png',captureBeyondViewport=FALSE))
  writeBin(base64enc::base64decode(result$data),file.path(here,'qa',paste0(name,'.png')))
}
rpc('Page.enable');rpc('Page.navigate',list(url='about:blank'))
rpc('Runtime.enable');rpc('Network.enable')
rpc('Network.setBlockedURLs',list(urls=c('http://*','https://*')))
rpc('Emulation.setDeviceMetricsOverride',list(width=1440L,height=1100L,deviceScaleFactor=1,mobile=FALSE))
messages <- list()
rpc('Page.navigate',list(url=paste0('file://',file.path(here,'report.html'))))
wait_for(function()isTRUE(js("document.readyState==='complete' && typeof jQuery==='function' && !!jQuery.fn.dataTable && jQuery.fn.dataTable.tables().length===1")))
labels <- c(paste('Chloramphenicol · replicate',1:3),paste('No antibiotic · replicate',1:3))
for(label in labels) {
  js(sprintf('Array.from(document.querySelectorAll(".nav-tabs a")).find(e=>e.textContent.trim()===%s).click()',jsonlite::toJSON(label,auto_unbox=TRUE)))
  later::run_now(.2)
}
js('document.querySelector(".nav-tabs a").click();window.scrollTo(0,0)')
screenshot('desktop_summary')
inventory <- js("({title:document.title,gt:document.querySelectorAll('.gt_table').length,tableRows:jQuery(jQuery.fn.dataTable.tables()[0]).DataTable().rows().count(),girafe:document.querySelectorAll('.girafe').length,girafeWithPlots:Array.from(document.querySelectorAll('.girafe')).filter(e=>e.querySelector('svg [data-id]')).length,plotly:document.querySelectorAll('.js-plotly-plot').length,images:Array.from(document.querySelectorAll('img')).map(x=>({alt:x.alt,loaded:x.complete&&x.naturalWidth>0})),overflow:document.documentElement.scrollWidth>window.innerWidth+5})")
jsonlite::write_json(inventory,file.path(here,'qa/browser_inventory.json'),pretty=TRUE,auto_unbox=TRUE)
stopifnot(inventory$gt==19,inventory$tableRows==1620,inventory$girafe==6,
  inventory$girafeWithPlots==6,inventory$plotly==5,!inventory$overflow,
  all(vapply(inventory$images,`[[`,logical(1),'loaded')))
stopifnot(isTRUE(js("document.querySelectorAll('section.level2').length===10 && Array.from(document.querySelectorAll('section.level2')).every(e=>e.closest('main'))")))
js("Array.from(document.querySelectorAll('.callout')).filter(e=>e.querySelector('.diagnostic-columns'))[1].querySelector('[data-bs-toggle=collapse]').click();Array.from(document.querySelectorAll('.callout')).filter(e=>e.querySelector('.diagnostic-columns'))[0].scrollIntoView()")
wait_for(function()isTRUE(js("!!Array.from(document.querySelectorAll('.callout')).filter(e=>e.querySelector('.diagnostic-columns'))[1].querySelector('.collapse.show')")))
screenshot('desktop_diagnostics')
js("document.querySelector('#following-the-mixture-over-time').scrollIntoView()")
screenshot('desktop_composition')
point <- js("(()=>{const e=document.querySelector('.girafe [data-id=\"All remaining barcodes\"]');const b=e.getBoundingClientRect();return {x:b.x+b.width/2,y:b.y+b.height/2};})()")
rpc('Input.dispatchMouseEvent',list(type='mouseMoved',x=point$x,y=point$y))
wait_for(function()isTRUE(js("Array.from(document.querySelectorAll('[class*=tooltip]')).some(e=>e.getBoundingClientRect().width>0&&e.textContent.includes('All remaining barcodes'))")))
screenshot('desktop_composition_hover')
js("Array.from(document.querySelectorAll('.nav-tabs a')).find(e=>e.textContent.trim()==='No antibiotic · replicate 3').click()")
stopifnot(isTRUE(js("Array.from(document.querySelectorAll('.nav-tabs a.active')).some(e=>e.textContent.trim()==='No antibiotic · replicate 3')")))
screenshot('desktop_control_replicate3')
js("document.querySelector('#what-agrees-with-the-published-measurements').scrollIntoView()")
screenshot('desktop_agreement')
js("document.querySelectorAll('.js-plotly-plot')[2].scrollIntoView()")
screenshot('desktop_diversity')
js("document.querySelectorAll('.js-plotly-plot')[3].scrollIntoView()")
screenshot('desktop_richness_comparison')
js("document.querySelectorAll('.js-plotly-plot')[4].scrollIntoView()")
screenshot('desktop_frequency_comparison')
barcode <- readRDS(file.path(here,'../generated/report_data.rds'))$explorer$Barcode[1]
rows_expected <- sum(readRDS(file.path(here,'../generated/report_data.rds'))$explorer$Barcode==barcode)
js(sprintf("jQuery(jQuery.fn.dataTable.tables()[0]).DataTable().search('%s').draw();true",barcode))
stopifnot(js("jQuery(jQuery.fn.dataTable.tables()[0]).DataTable().rows({search:'applied'}).count()") == rows_expected)
js("document.querySelector('#inspect-the-observations-behind-the-plots').scrollIntoView()")
screenshot('desktop_barcode_search')
download <- tempfile('barbac-report-downloads-');dir.create(download)
rpc('Browser.setDownloadBehavior',list(behavior='allow',downloadPath=download))
js('jQuery(jQuery.fn.dataTable.tables()[0]).DataTable().button(0).trigger();true')
wait_for(function()length(list.files(download,pattern='[.]csv$'))==1)
csv <- read.csv(list.files(download,pattern='[.]csv$',full.names=TRUE)[1],check.names=FALSE)
stopifnot(nrow(csv)==rows_expected,all(csv$Barcode==barcode),
  all(csv[['Barbac / input reads (%)']]>=0 & csv[['Barbac / input reads (%)']]<=100))
original <- readRDS(file.path(here,'../generated/report_data.rds'))$explorer
original <- original[original$Barcode==barcode,]
key <- function(pop,passage)paste(pop,passage)
idx <- match(key(csv$Population,csv$Passage),key(original$Population,original$Passage))
stopifnot(!anyNA(idx),
  isTRUE(all.equal(csv[['Barbac / input reads (%)']],original$Barbac_frequency_input_percent[idx])),
  isTRUE(all.equal(csv[['Barbac / extracted reads (%)']],original$Barbac_frequency_extracted_percent[idx])))
rpc('Emulation.setDeviceMetricsOverride',list(width=390L,height=844L,deviceScaleFactor=1,mobile=TRUE))
js('window.scrollTo(0,0)');screenshot('mobile_summary')
mobile <- js('({width:window.innerWidth,scrollWidth:document.documentElement.scrollWidth,overflow:document.documentElement.scrollWidth>window.innerWidth+5})')
stopifnot(!mobile$overflow)
exceptions <- Filter(function(x)identical(x$method,'Runtime.exceptionThrown'),messages)
stopifnot(length(exceptions)==0)
jsonlite::write_json(list(status='passed',offline_network_blocked=TRUE,desktop=inventory,mobile=mobile,
  sample_expansion=TRUE,composition_tab_switch=TRUE,composition_hover_tooltip=TRUE,
  barcode_search=TRUE,CSV_download_verified=TRUE,javascript_exceptions=length(exceptions),
  report_sha256=digest::digest(file=file.path(here,'report.html'),algo='sha256'),
  screenshots=list.files(file.path(here,'qa'),pattern='[.]png$')),
  file.path(here,'qa/browser_validation.json'),pretty=TRUE,auto_unbox=TRUE)
validation <- jsonlite::read_json(file.path(here,'validation.json'))
assert_hash_map(validation$hashes,here)
validation$status <- 'passed';validation$browser_validation <- 'passed'
validation$browser_checks_sha256 <- digest::digest(file=file.path(here,'browser_checks.R'),algo='sha256')
validation$browser_receipt_sha256 <- digest::digest(file=file.path(here,'qa/browser_validation.json'),algo='sha256')
jsonlite::write_json(validation,file.path(here,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
ws$close();cat('Offline rendering, tables, panels, tabs, hover, search, CSV and mobile checks passed.\n')
