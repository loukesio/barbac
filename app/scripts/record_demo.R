#!/usr/bin/env Rscript
# Record actual local app interactions; output contains synthetic examples only.
repo<-normalizePath('.');app<-file.path(repo,'app');media<-file.path(app,'media')
dir.create(media,showWarnings=FALSE);work<-tempfile('video-',tmpdir=file.path(app,'.qa'));dir.create(work)
source(file.path(app,'R','examples.R'))
pages<-jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory('http://127.0.0.1:9223/json/list')$content))
ws<-websocket::WebSocket$new(pages$webSocketDebuggerUrl[pages$type=='page'][1],maxMessageSize=64*1024^2)
responses<-new.env(parent=emptyenv());id<-0L
ws$onMessage(function(event){v<-jsonlite::fromJSON(event$data,simplifyVector=FALSE);if(!is.null(v$id))responses[[as.character(v$id)]]<-v})
wait_for<-function(predicate,timeout=120){start<-Sys.time();while(!isTRUE(predicate())){if(difftime(Sys.time(),start,units='secs')>timeout)stop('Recording step timed out');later::run_now(.05)}}
wait_for(function()ws$readyState()==1)
rpc<-function(method,params=list()){
  id<<-id+1L;key<-as.character(id);msg<-list(id=id,method=method);if(length(params))msg$params<-params
  ws$send(jsonlite::toJSON(msg,auto_unbox=TRUE,null='null'));wait_for(function()exists(key,responses,inherits=FALSE))
  x<-responses[[key]];if(!is.null(x$error))stop(jsonlite::toJSON(x$error));x$result
}
js<-function(code){x<-rpc('Runtime.evaluate',list(expression=code,returnByValue=TRUE,awaitPromise=TRUE));if(!is.null(x$exceptionDetails))stop(jsonlite::toJSON(x$exceptionDetails));x$result$value}
set_input<-function(id,value)invisible(js(sprintf('Shiny.setInputValue(%s,%s,{priority:"event"});true',jsonlite::toJSON(id,auto_unbox=TRUE),jsonlite::toJSON(value,auto_unbox=TRUE))))
click<-function(id)invisible(js(sprintf('document.getElementById(%s).click();true',jsonlite::toJSON(id,auto_unbox=TRUE))))
nav<-function(value)invisible(js(sprintf('document.querySelector("#page input[value=%s]").click();window.scrollTo(0,0);true',value)))
upload<-function(id,path){doc<-rpc('DOM.getDocument');node<-rpc('DOM.querySelector',list(nodeId=doc$root$nodeId,selector=paste0('#',id)));invisible(rpc('DOM.setFileInputFiles',list(nodeId=node$nodeId,files=list(normalizePath(path)))))}
caption<-function(title,subtitle){
  invisible(js(sprintf('(()=>{let e=document.getElementById("recording-caption");if(!e){e=document.createElement("div");e.id="recording-caption";e.style.cssText="position:fixed;z-index:99999;left:280px;right:36px;bottom:22px;background:rgba(20,60,54,.97);border:1px solid #799887;border-radius:13px;padding:17px 25px;color:#fff;box-shadow:0 8px 30px #163f3a26;pointer-events:none;font-family:Avenir,Arial,sans-serif";document.body.appendChild(e);}e.replaceChildren();const h=document.createElement("div");h.style.cssText="font-size:18px;font-weight:600;letter-spacing:-.3px";h.textContent=%s;const p=document.createElement("div");p.style.cssText="font-size:12px;color:#bfd0be;margin-top:5px";p.textContent=%s;e.append(h,p);return true;})()',jsonlite::toJSON(title,auto_unbox=TRUE),jsonlite::toJSON(subtitle,auto_unbox=TRUE))))
}
frames<-character();times<-numeric();folder<-NULL
frame<-function(){
  path<-file.path(folder,sprintf('frame-%05d.jpg',length(frames)+1))
  x<-rpc('Page.captureScreenshot',list(format='jpeg',quality=88L,captureBeyondViewport=FALSE))
  writeBin(base64enc::base64decode(x$data),path);frames<<-c(frames,path);times<<-c(times,as.numeric(Sys.time()))
}
hold<-function(seconds){deadline<-Sys.time()+seconds;while(Sys.time()<deadline){frame();later::run_now(.13)};invisible(NULL)}
await_record<-function(predicate,timeout=120){start<-Sys.time();while(!isTRUE(js(predicate))){if(difftime(Sys.time(),start,units='secs')>timeout)stop('Recorded job timed out');frame();later::run_now(.15)}}
start_clip<-function(name){folder<<-file.path(work,name);dir.create(folder);frames<<-character();times<<-numeric()}
finish_clip<-function(name){
  durations<-c(diff(times),.3)
  entries<-unlist(Map(function(path,d)c(paste0("file '",path,"'"),paste('duration',sprintf('%.6f',d))),frames,durations))
  manifest<-file.path(folder,'frames.txt');writeLines(c(entries,paste0("file '",tail(frames,1),"'")),manifest)
  output<-file.path(media,paste0(name,'.mp4'))
  status<-system2(Sys.which('ffmpeg'),c('-y','-hide_banner','-loglevel','error','-f','concat','-safe','0','-i',shQuote(manifest),
    '-vf',shQuote('fps=24,format=yuv420p'),'-c:v','libx264','-preset','slow','-crf','23','-movflags','+faststart',shQuote(output)))
  if(status!=0)stop('Video encoding failed')
  cat('Recorded',output,'\n')
  list(file=basename(output),duration_seconds=sum(durations),frames=length(frames),synthetic_data=TRUE,continuous_capture=TRUE)
}
invisible(rpc('Page.enable'));invisible(rpc('Runtime.enable'));invisible(rpc('Network.enable'))
invisible(rpc('Network.setCacheDisabled',list(cacheDisabled=TRUE)))
invisible(rpc('Emulation.setDeviceMetricsOverride',list(width=1440L,height=960L,deviceScaleFactor=1,mobile=FALSE)))
invisible(rpc('Page.navigate',list(url='http://127.0.0.1:3838')))
wait_for(function()isTRUE(js('!!document.querySelector("#input_table table")')))
download<-file.path(work,'downloads');dir.create(download)
invisible(rpc('Browser.setDownloadBehavior',list(behavior='allow',downloadPath=download)))
example<-file.path(work,'example.csv');readr::write_csv(studio_demo(),example)
start_clip('studio-walkthrough')
caption('barbac Studio','A complete barcode analysis, from sequences to a shareable result.');hold(3)
caption('01 · Bring your extracted barcode counts','CSV or TSV, with optional sample and timepoint metadata.');upload('counts_files',example);hold(3)
invisible(js('document.querySelector("#input_metrics").scrollIntoView({block:"center",behavior:"smooth"});true'));hold(3)
nav('cluster');caption('02 · Cluster with the native barbac engine','Choose LV or Hamming. Keep advanced options explicit.');hold(4)
click('run');caption('Clustering runs in the background','The original sample counts and independent populations are preserved.')
await_record('!!document.querySelector("#result_metrics .metric-value")');hold(1)
await_record('!!document.querySelector("#area_interactive svg [data-id]")')
caption('03 · Explore every inferred lineage','18 clusters. All 192,000 synthetic barcode reads conserved.');hold(4)
set_input('palette','dora');caption('Choose your palette','Native LTC colours. The data and clustering stay the same.');hold(3)
invisible(js('document.querySelector(".results-lower").scrollIntoView({block:"start",behavior:"smooth"});true'));hold(3)
invisible(js('document.querySelector("#result_table").scrollIntoView({block:"center",behavior:"smooth"});true'));set_input('table_view','memberships')
caption('Inspect the assignments','Search centroids, every member sequence, and counts in every sample.');hold(4)
nav('downloads');caption('04 · Take the complete analysis','Download counts, memberships, statistics, settings and publication figures.');hold(3)
click('download_all');hold(1)
click('build_report');caption('Create a report you can share','Quarto packages your settings, results and lineage figures into one HTML file.')
await_record('document.querySelector("#download_report")?.getAttribute("href")?.includes("/download/")');hold(2)
click('download_report');caption('Ready to share with your team','Complete results, inspectable settings, and the same native barbac engine.');hold(3)
one<-finish_clip('studio-walkthrough')

fixture<-file.path(work,'raw');cfg<-studio_fastq_example(fixture)
nav('data');set_input('source_type','fastq');invisible(js('document.querySelector(".raw-card").scrollIntoView({block:"start"});true'))
start_clip('studio-fastq')
caption('Start from raw sequencing reads','Single-end or overlapping paired-end FASTQ.gz files.');hold(3)
upload('raw_r1',file.path(fixture,'example_R1.fastq.gz'));upload('raw_r2',file.path(fixture,'example_R2.fastq.gz'))
upload('reference',file.path(fixture,'reference.fasta'));hold(3)
click('demo_flanks');caption('Define your barcode cassette','Supply a reference, barcode coordinates, and exact flanks that preserve indels.');hold(4)
click('extract');caption('Merge · map · extract','PEAR, minimap2 and samtools prepare the reads for native barbac extraction.')
await_record('document.querySelector("#input_caption").textContent.includes("Extracted from your FASTQ")')
nav('data');set_input('source_type','counts');invisible(js('document.querySelector("#input_metrics").scrollIntoView({block:"center"});true'))
caption('Extracted counts are ready','This synthetic fixture recovers 140 molecules across three known barcodes.');hold(4)
nav('downloads');click('download_input');caption('Download now, or continue to clustering','Save extracted barcode counts without having to run clustering.');hold(3)
nav('cluster');click('run');await_record('document.querySelector("#result_metrics .metric-value")?.textContent==="3"')
caption('The same workflow, from either starting point','Extracted sequences or raw reads. All counts remain inspectable.');hold(4)
two<-finish_clip('studio-fastq')
invisible(js('document.getElementById("recording-caption")?.remove();true'))
jsonlite::write_json(list(recorded_utc=format(Sys.time(),tz='UTC',usetz=TRUE),source='Real local Shiny app',
  clips=list(one,two),note='Synthetic fixtures only. No timing speedup, invented screen or fabricated result.'),
  file.path(media,'recording.json'),pretty=TRUE,auto_unbox=TRUE)
ws$close()
