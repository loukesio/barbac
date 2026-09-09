# Optional R-driven browser QA. Start a temporary headless Chrome on port 9223.
# JavaScript expressions below inspect the UI only; all analysis is done in R.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
here <- dirname(normalizePath(script))
dir.create(file.path(here, "qa"), showWarnings = FALSE)
endpoints <- jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory("http://127.0.0.1:9223/json/list")$content))
ws <- websocket::WebSocket$new(endpoints$webSocketDebuggerUrl[endpoints$type == "page"][1], maxMessageSize = 64*1024^2)
responses <- new.env(parent = emptyenv()); messages <- list(); id <- 0L
ws$onMessage(function(event) {
  value <- jsonlite::fromJSON(event$data, simplifyVector = FALSE)
  if (!is.null(value$id)) responses[[as.character(value$id)]] <- value else messages[[length(messages)+1L]] <<- value
})
wait_for <- function(predicate, timeout = 45) {
  start <- Sys.time()
  while (!predicate()) {
    if (as.numeric(difftime(Sys.time(), start, units = "secs")) > timeout) stop("Browser QA timed out")
    later::run_now(.05)
  }
}
wait_for(function() ws$readyState() == 1)
rpc <- function(method, params = list()) {
  id <<- id + 1L; key <- as.character(id)
  command <- list(id = id, method = method)
  if (length(params)) command$params <- params
  ws$send(jsonlite::toJSON(command, auto_unbox = TRUE, null = "null"))
  wait_for(function() exists(key, responses, inherits = FALSE))
  result <- responses[[key]]
  if (!is.null(result$error)) stop(jsonlite::toJSON(result$error, auto_unbox = TRUE))
  result$result
}
js <- function(expression) {
  value <- rpc("Runtime.evaluate", list(expression = expression, returnByValue = TRUE, awaitPromise = TRUE))
  if (!is.null(value$exceptionDetails)) stop(jsonlite::toJSON(value$exceptionDetails, auto_unbox = TRUE))
  value$result$value
}
rpc("Page.enable")
rpc("Page.navigate", list(url = "about:blank"))
rpc("Runtime.enable"); rpc("Network.enable")
rpc("Network.setBlockedURLs", list(urls = c("http://*", "https://*")))
rpc("Emulation.setDeviceMetricsOverride", list(width = 1440L, height = 1100L, deviceScaleFactor = 1, mobile = FALSE))
messages <- list()
rpc("Page.navigate", list(url = paste0("file://", file.path(here, "report.html"))))
wait_for(function() isTRUE(js("document.readyState === 'complete'")))
jsonlite::write_json(list(state = js("({jquery:typeof jQuery,DT:typeof jQuery==='function'&&!!jQuery.fn.dataTable,widgets:document.querySelectorAll('.html-widget').length,body:document.body.innerText.slice(0,300)})"),
  exceptions = Filter(function(x) identical(x$method, "Runtime.exceptionThrown"), messages)),
  file.path(here, "qa", "initial_browser_state.json"), pretty = TRUE, auto_unbox = TRUE)
wait_for(function() isTRUE(js("document.readyState === 'complete' && typeof jQuery === 'function' && !!jQuery.fn.dataTable && jQuery.fn.dataTable.tables().length >= 5")))
# Quarto initializes widgets in hidden tabs when those tabs are first opened.
for (label in c("LV · R1", "LV · R2", "Hamming · R1", "Hamming · R2", "Plot label lookup", "LV counts", "Hamming counts")) {
  js(sprintf("Array.from(document.querySelectorAll('.nav-tabs a')).find(e=>e.textContent.trim()===%s).click()", jsonlite::toJSON(label, auto_unbox = TRUE)))
  later::run_now(.25)
}
wait_for(function() js("jQuery.fn.dataTable.tables().length") >= 8)
js("Array.from(document.querySelectorAll('.nav-tabs a')).find(e=>e.textContent.trim()==='LV · R1').click(); Array.from(document.querySelectorAll('.nav-tabs a')).find(e=>e.textContent.trim()==='Plot label lookup').click(); window.scrollTo(0,0)")
later::run_now(1)
screenshot <- function(name) {
  js("new Promise(resolve=>setTimeout(resolve,600))")
  result <- rpc("Page.captureScreenshot", list(format = "png", captureBeyondViewport = FALSE))
  writeBin(base64enc::base64decode(result$data), file.path(here, "qa", paste0(name, ".png")))
}
screenshot("desktop_summary")
inventory <- js("({title:document.title, tables:Array.from(jQuery.fn.dataTable.tables()).map(t=>({id:t.id, caption:t.closest('.dataTables_wrapper').querySelector('caption')?.textContent||'', rows:jQuery(t).DataTable().rows().count()})), girafe:document.querySelectorAll('.girafe').length, girafeWithPlots:Array.from(document.querySelectorAll('.girafe')).filter(e=>e.querySelector('svg [data-id]')).length, plotly:document.querySelectorAll('.js-plotly-plot').length, images:Array.from(document.querySelectorAll('img')).map(x=>({alt:x.alt,loaded:x.complete&&x.naturalWidth>0})), horizontalOverflow:document.documentElement.scrollWidth>window.innerWidth+5, overflowElements:Array.from(document.querySelectorAll('body *')).filter(e=>e.getBoundingClientRect().width>0&&e.getBoundingClientRect().right>window.innerWidth+5).slice(0,12).map(e=>({tag:e.tagName,id:e.id,class:e.className,right:e.getBoundingClientRect().right}))})")
jsonlite::write_json(inventory, file.path(here, "qa", "browser_inventory.json"), pretty = TRUE, auto_unbox = TRUE)
stopifnot(length(inventory$tables) >= 8, inventory$girafe == 4, inventory$girafeWithPlots == 4,
  inventory$plotly == 3, all(vapply(inventory$images, function(x) x$loaded, logical(1))),
  !inventory$horizontalOverflow)
js("Array.from(document.querySelectorAll('h2')).find(e=>e.textContent.startsWith('Extraction retains')).scrollIntoView()")
screenshot("desktop_extraction")
js("Array.from(document.querySelectorAll('details'))[1].open=true")
stopifnot(isTRUE(js("Array.from(document.querySelectorAll('details'))[1].open")))
js("Array.from(document.querySelectorAll('h2')).find(e=>e.textContent.startsWith('Major barcode')).scrollIntoView()")
screenshot("desktop_composition")
point <- js("(()=>{const e=document.querySelector('.girafe [data-id=\"Other\"]');const b=e.getBoundingClientRect();return {x:b.x+b.width/2,y:b.y+b.height/2};})()")
rpc("Input.dispatchMouseEvent", list(type = "mouseMoved", x = point$x, y = point$y))
wait_for(function() isTRUE(js("Array.from(document.querySelectorAll('[class*=tooltip]')).some(e=>e.getBoundingClientRect().width>0&&e.textContent.includes('Other'))")))
screenshot("desktop_composition_hover")
js("Array.from(document.querySelectorAll('.nav-tabs a')).find(e=>e.textContent.trim()==='Hamming · R2').click()")
stopifnot(isTRUE(js("Array.from(document.querySelectorAll('.nav-tabs li.active a,.nav-tabs a.active')).some(e=>e.textContent.trim()==='Hamming · R2')")))
screenshot("desktop_hamming_r2")
js("Array.from(document.querySelectorAll('.nav-tabs a')).find(e=>e.textContent.trim()==='LV counts').click()")
lv_id <- Filter(function(x) grepl("Complete LV barcode", x$caption), inventory$tables)[[1]]$id
ham_id <- Filter(function(x) grepl("Complete Hamming barcode", x$caption), inventory$tables)[[1]]$id
stopifnot(js(sprintf("jQuery('#%s').DataTable().rows().count()", lv_id)) == 13687,
          js(sprintf("jQuery('#%s').DataTable().rows().count()", ham_id)) == 29526)
barcode <- "ATAAAAAAGCACAAGCCTTTTGACGA_CACCTAAATTAGTTATCCATTCGGCT"
js(sprintf("jQuery('#%s').DataTable().search('%s').draw(); true", lv_id, barcode))
stopifnot(js(sprintf("jQuery('#%s').DataTable().rows({search:'applied'}).count()", lv_id)) == 1)
js("Array.from(document.querySelectorAll('h2')).find(e=>e.textContent.startsWith('Every paired')).scrollIntoView()")
screenshot("desktop_barcode_search")
download <- file.path(tempdir(), "barbac-report-downloads"); dir.create(download, showWarnings = FALSE)
rpc("Browser.setDownloadBehavior", list(behavior = "allow", downloadPath = download))
js(sprintf("jQuery('#%s').DataTable().button(0).trigger(); true", lv_id))
wait_for(function() length(list.files(download, pattern = "\\.csv$")) > 0)
csv <- suppressWarnings(read.csv(list.files(download, pattern = "\\.csv$", full.names = TRUE)[1], check.names = FALSE))
stopifnot(nrow(csv) == 1, csv$Barcode == barcode, csv$pooled_molecules == 1531711)
rpc("Emulation.setDeviceMetricsOverride", list(width = 390L, height = 844L, deviceScaleFactor = 1, mobile = TRUE))
js("window.scrollTo(0,0)"); later::run_now(.5); screenshot("mobile_summary")
mobile <- js("({width:window.innerWidth, scrollWidth:document.documentElement.scrollWidth, overflow:document.documentElement.scrollWidth>window.innerWidth+5})")
stopifnot(!mobile$overflow)
exceptions <- Filter(function(x) identical(x$method, "Runtime.exceptionThrown"), messages)
stopifnot(length(exceptions) == 0)
result <- list(status = "passed", offline_network_blocked = TRUE, desktop = inventory, mobile = mobile,
  sample_expansion = TRUE, composition_tab_switch = TRUE, complete_table_row_counts = TRUE,
  composition_hover_tooltip = TRUE, barcode_search = TRUE, CSV_download_verified = TRUE, javascript_exceptions = length(exceptions),
  screenshots = list.files(file.path(here, "qa"), pattern = "\\.png$"))
jsonlite::write_json(result, file.path(here, "qa", "browser_validation.json"), pretty = TRUE, auto_unbox = TRUE)
ws$close()
cat("Browser rendering, offline use, sample/tab controls, full tables, search and CSV export passed.\n")
