#!/usr/bin/env Rscript
# Drive an existing Chrome debugging session. No test-only server endpoints.
repo<-normalizePath('.');app<-file.path(repo,'app');qa<-file.path(app,'.qa')
dir.create(qa,showWarnings=FALSE)
source(file.path(app,'R','examples.R'))
endpoints<-jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory('http://127.0.0.1:9223/json/list')$content))
ws<-websocket::WebSocket$new(endpoints$webSocketDebuggerUrl[endpoints$type=='page'][1],maxMessageSize=128*1024^2)
responses<-new.env(parent=emptyenv());events<-list();id<-0L
ws$onMessage(function(event){v<-jsonlite::fromJSON(event$data,simplifyVector=FALSE);if(!is.null(v$id))responses[[as.character(v$id)]]<-v else events[[length(events)+1L]]<<-v})
wait_for<-function(predicate,timeout=90){start<-Sys.time();while(!isTRUE(predicate())){if(as.numeric(difftime(Sys.time(),start,units='secs'))>timeout)stop('Browser check timed out');later::run_now(.08)}}
wait_for(function()ws$readyState()==1)
rpc<-function(method,params=list()){
  id<<-id+1L;key<-as.character(id);msg<-list(id=id,method=method)
  if(length(params))msg$params<-params
  ws$send(jsonlite::toJSON(msg,auto_unbox=TRUE,null='null'));wait_for(function()exists(key,responses,inherits=FALSE))
  value<-responses[[key]];if(!is.null(value$error))stop(jsonlite::toJSON(value$error,auto_unbox=TRUE));value$result
}
js<-function(code){v<-rpc('Runtime.evaluate',list(expression=code,returnByValue=TRUE,awaitPromise=TRUE));if(!is.null(v$exceptionDetails))stop(jsonlite::toJSON(v$exceptionDetails,auto_unbox=TRUE));v$result$value}
nap<-function()js('new Promise(resolve=>setTimeout(resolve,500))')
shot<-function(name){nap();js('document.querySelectorAll(".shiny-notification-close").forEach(e=>e.click());true');r<-rpc('Page.captureScreenshot',list(format='png',captureBeyondViewport=FALSE));writeBin(base64enc::base64decode(r$data),file.path(qa,paste0(name,'.png')))}
input<-function(name,value)js(sprintf('Shiny.setInputValue(%s,%s,{priority:"event"});true',jsonlite::toJSON(name,auto_unbox=TRUE),jsonlite::toJSON(value,auto_unbox=TRUE)))
click<-function(id)js(sprintf('document.getElementById(%s).click();true',jsonlite::toJSON(id,auto_unbox=TRUE)))
page<-function(value){input('page',value);nap()}
upload<-function(id,files){doc<-rpc('DOM.getDocument');node<-rpc('DOM.querySelector',list(nodeId=doc$root$nodeId,selector=paste0('#',id)));rpc('DOM.setFileInputFiles',list(nodeId=node$nodeId,files=unname(as.list(normalizePath(files)))));nap()}
no_errors<-function(){errors<-js('Array.from(document.querySelectorAll(".shiny-output-error")).filter(e=>e.offsetParent!==null).map(e=>e.innerText)');if(length(errors))stop(paste(errors,collapse='; '))}
if('--site' %in% commandArgs(TRUE)) {
  invisible(rpc('Page.enable'));invisible(rpc('Runtime.enable'))
  invisible(rpc('Emulation.setDeviceMetricsOverride',list(width=1440L,height=1080L,deviceScaleFactor=1,mobile=FALSE)))
  invisible(rpc('Page.navigate',list(url=paste0('file://',file.path(repo,'docs','index.html')))))
  wait_for(function()isTRUE(js('document.readyState==="complete" && document.title.includes("barbac")')))
  stopifnot(isTRUE(js('[...document.images].every(i=>i.complete && i.naturalWidth>0)')))
  shot('publication-readme-desktop')
  stopifnot(isTRUE(js('document.body.innerText.includes("R1-only") && document.body.innerText.includes("barbac Studio")')))
  invisible(rpc('Emulation.setDeviceMetricsOverride',list(width=390L,height=844L,deviceScaleFactor=1,mobile=TRUE)))
  shot('publication-readme-mobile')
  stopifnot(isTRUE(js('document.documentElement.scrollWidth<=window.innerWidth+2')))
  ws$close();cat('Publication README images and desktop/mobile layouts passed.\n');quit(status=0)
}
if('--media' %in% commandArgs(TRUE)) {
  invisible(rpc('Page.enable'));invisible(rpc('Runtime.enable'));invisible(rpc('Network.enable'))
  invisible(rpc('Network.setBlockedURLs',list(urls=list('http://*','https://*'))))
  invisible(rpc('Emulation.setDeviceMetricsOverride',list(width=1440L,height=1080L,deviceScaleFactor=1,mobile=FALSE)))
  invisible(rpc('Page.navigate',list(url=paste0('file://',file.path(app,'media','watch.html')))))
  wait_for(function()isTRUE(js('document.querySelectorAll("video").length===2 && [...document.querySelectorAll("video")].every(v=>v.readyState>=1)')))
  stopifnot(isTRUE(js('[...document.querySelectorAll("video")].every(v=>!v.error && v.videoWidth===1440 && v.videoHeight===960)')))
  for(index in 0:1) {
    stopifnot(isTRUE(js(sprintf('(async()=>{const v=document.querySelectorAll("video")[%d];v.muted=true;await v.play();await new Promise(r=>setTimeout(r,700));v.pause();return v.currentTime>0;})()',index))))
  }
  shot('video-viewer-desktop')
  invisible(rpc('Emulation.setDeviceMetricsOverride',list(width=390L,height=844L,deviceScaleFactor=1,mobile=TRUE)))
  shot('video-viewer-mobile')
  stopifnot(isTRUE(js('document.documentElement.scrollWidth<=window.innerWidth+2')))
  remote_requests<-Filter(function(e)identical(e$method,'Network.requestWillBeSent')&&grepl('^https?://',e$params$request$url),events)
  stopifnot(length(remote_requests)==0)
  jsonlite::write_json(list(status='passed',offline=TRUE,videos=2,playback=TRUE,mobile=TRUE,network_requests=0),
    file.path(qa,'media-validation.json'),pretty=TRUE,auto_unbox=TRUE)
  invisible(rpc('Network.setBlockedURLs',list(urls=list('https://*'))))
  ws$close();cat('Offline videos play successfully; desktop and mobile checks passed.\n');quit(status=0)
}
if('--inspect' %in% commandArgs(TRUE)) {
  print(js('({url:location.href,ready:document.readyState,shiny:typeof Shiny,connected:typeof Shiny!=="undefined"&&!!Shiny.shinyapp?.isConnected(),values:Shiny.shinyapp.$values,inputs:Shiny.shinyapp.$inputValues,text:document.body.innerText.slice(0,9000),tables:document.querySelectorAll("table").length})'))
  shot('inspect');ws$close();quit(status=0)
}
if('--diagnose' %in% commandArgs(TRUE)) {
  rpc('Network.enable');rpc('Runtime.enable');rpc('Page.navigate',list(url='http://127.0.0.1:3838'))
  js('new Promise(resolve=>setTimeout(resolve,5000))')
  for(e in events)if(grepl('webSocketFrame|exceptionThrown',e$method))print(e)
  ws$close();quit(status=0)
}
rpc('Page.enable');rpc('Runtime.enable');rpc('Network.enable')
rpc('Network.setCacheDisabled',list(cacheDisabled=TRUE))
rpc('Network.setBlockedURLs',list(urls=list('https://*')))
rpc('Emulation.setDeviceMetricsOverride',list(width=1440L,height=1080L,deviceScaleFactor=1,mobile=FALSE))
rpc('Page.navigate',list(url='http://127.0.0.1:3838'))
wait_for(function()isTRUE(js('typeof Shiny!=="undefined" && Shiny.shinyapp && Shiny.shinyapp.isConnected() && !!document.querySelector("#input_table table")')))
shot('desktop-data');no_errors()
stopifnot(isTRUE(js('document.body.innerText.includes("192,000")')))
page('cluster');shot('desktop-cluster');click('run')
wait_for(function()isTRUE(js('!!document.querySelector("#result_metrics .metric-value")')))
wait_for(function()isTRUE(js('!!document.querySelector("#area_interactive svg [data-id]")')))
shot('desktop-results');no_errors()
stopifnot(isTRUE(js('document.querySelector("#result_metrics .metric-value").textContent==="18"')))
input('palette','dora');nap()
stopifnot(isTRUE(js('!!document.querySelector("#area_interactive svg [data-id]")')))
download<-tempfile('downloads-',tmpdir=qa);dir.create(download,showWarnings=FALSE)
rpc('Browser.setDownloadBehavior',list(behavior='allow',downloadPath=download))
page('downloads');shot('desktop-downloads');click('download_all')
wait_for(function()file.exists(file.path(download,'barbac-analysis.zip')))
unpack<-tempfile('download-inspection-',tmpdir=qa);dir.create(unpack,showWarnings=FALSE)
utils::unzip(file.path(download,'barbac-analysis.zip'),exdir=unpack)
stopifnot(sum(read.csv(file.path(unpack,'time_series.csv'))$counts)==192000,nrow(read.csv(file.path(unpack,'centroids.csv')))==18)
click('build_report');wait_for(function()isTRUE(js('document.querySelector("#download_report")?.getAttribute("href")?.includes("/download/")')),timeout=150)
click('download_report');wait_for(function()file.exists(file.path(download,'barbac-report.html')))
stopifnot(file.info(file.path(download,'barbac-report.html'))$size>10000)
page('results');rpc('Emulation.setDeviceMetricsOverride',list(width=390L,height=844L,deviceScaleFactor=1,mobile=TRUE));nap();shot('mobile-results')
stopifnot(isTRUE(js('document.documentElement.scrollWidth<=window.innerWidth+2')))
js('document.querySelector(".plot-card").scrollIntoView({block:"start"});true');shot('mobile-plot')
page('data');shot('mobile-data');no_errors()
stopifnot(isTRUE(js('document.documentElement.scrollWidth<=window.innerWidth+2')))
rpc('Emulation.setDeviceMetricsOverride',list(width=1440L,height=1080L,deviceScaleFactor=1,mobile=FALSE))
bad<-file.path(qa,'invalid-counts.csv');writeLines(c('barcode,counts','ACGT,-2'),bad)
upload('counts_files',bad)
wait_for(function()isTRUE(js('document.querySelector("#import_message").textContent.includes("positive whole")')))
good<-file.path(qa,'uploaded-example.csv');readr::write_csv(studio_demo(),good)
upload('counts_files',good)
wait_for(function()isTRUE(js('document.querySelector("#input_metrics").textContent.includes("192,000")')))
stopifnot(isTRUE(js('document.querySelector("#dataset_badge").textContent.includes("YOUR DATA")')))
# Exercise real paired FASTQ upload and extraction, not just server input mocks.
fixture<-file.path(qa,'raw-fixture');cfg<-studio_fastq_example(fixture)
input('source_type','fastq');nap()
upload('raw_r1',file.path(fixture,'example_R1.fastq.gz'))
upload('raw_r2',file.path(fixture,'example_R2.fastq.gz'))
upload('reference',file.path(fixture,'reference.fasta'))
click('demo_flanks');nap();shot('desktop-raw');click('extract')
wait_for(function()isTRUE(js('document.querySelector("#input_caption").textContent.includes("Extracted from your FASTQ")')),timeout=150)
stopifnot(isTRUE(js('document.querySelector("#input_metrics").textContent.includes("140")')))
click('run');wait_for(function()isTRUE(js('document.querySelector("#result_metrics .metric-value")?.textContent==="3"')))
no_errors()
exceptions<-Filter(function(x)identical(x$method,'Runtime.exceptionThrown'),events)
stopifnot(length(exceptions)==0)
jsonlite::write_json(list(status='passed',desktop=TRUE,mobile=TRUE,upload_validation=TRUE,
  synthetic_clusters=18,synthetic_counts=192000,paired_extraction_counts=140,paired_clusters=3,
  csv_zip_download=TRUE,quarto_download=TRUE,javascript_exceptions=length(exceptions)),
  file.path(qa,'browser-validation.json'),pretty=TRUE,auto_unbox=TRUE)
cat('Browser checks passed. Screenshots:',qa,'\n')
ws$close()
