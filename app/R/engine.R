# Data and analysis adapter. All clustering uses the installed barbac engine.
studio_error <- function(...) stop(..., call. = FALSE)
studio_sum <- function(x, groups) stats::aggregate(x, groups, sum)

studio_read_table <- function(path, name = basename(path)) {
  delim <- if (grepl('[.](tsv|txt)([.]gz)?$', name, ignore.case = TRUE)) '\t' else ','
  x <- readr::read_delim(path, delim = delim, col_types = readr::cols(.default = 'c'),
    name_repair = 'minimal', progress = FALSE, show_col_types = FALSE)
  if (nrow(readr::problems(x))) studio_error('Some rows have the wrong number of columns in ', name, '.')
  if (anyDuplicated(names(x))) studio_error('Column names must be unique in ', name, '.')
  names(x) <- tolower(trimws(names(x)))
  if (anyDuplicated(names(x))) studio_error('Column names differ only by case in ', name, '.')
  as.data.frame(x, stringsAsFactors = FALSE)
}

studio_validate <- function(x) {
  if (!is.data.frame(x) || !nrow(x)) studio_error('The barcode table is empty.')
  if (!all(c('barcode', 'counts') %in% names(x)))
    studio_error('Your table needs columns named barcode and counts. Optional: sample, time, population.')
  if (nrow(x) > as.numeric(Sys.getenv('BARBAC_STUDIO_MAX_ROWS', '2000000')))
    studio_error('This table exceeds the configured row limit. Use the local R workflow for larger inputs.')
  x$barcode <- toupper(trimws(as.character(x$barcode)))
  if (anyNA(x$barcode) || any(!grepl('^[ACGTN]+$', x$barcode)))
    studio_error('Barcodes must contain A, C, G, T or N only, with no empty sequences.')
  if (any(nchar(x$barcode) > 250)) studio_error('Barcodes longer than 250 bases are outside this app’s input range.')
  x$counts <- suppressWarnings(as.numeric(x$counts))
  if (anyNA(x$counts) || any(!is.finite(x$counts) | x$counts <= 0 | x$counts != floor(x$counts)))
    studio_error('Counts must be positive whole numbers. Remove empty or zero-count rows.')
  if (!'sample' %in% names(x)) x$sample <- 'Sample 1'
  if (!'population' %in% names(x)) x$population <- 'Library 1'
  for (field in c('sample', 'population')) {
    x[[field]] <- trimws(as.character(x[[field]]))
    if (anyNA(x[[field]]) || any(!nzchar(x[[field]])) || any(nchar(x[[field]]) > 100))
      studio_error('Every row needs a nonempty ', field, ' label of at most 100 characters.')
  }
  if (!'time' %in% names(x)) x$time <- NA_real_
  supplied <- !is.na(x$time) & nzchar(trimws(as.character(x$time)))
  times <- suppressWarnings(as.numeric(x$time))
  if (any(supplied & !is.finite(times))) studio_error('Timepoints must be numeric (for example 0, 8, 16).')
  x$time <- times
  metadata <- unique(x[c('sample', 'population', 'time')])
  if (anyDuplicated(metadata$sample))
    studio_error('Each sample must map to exactly one population and one timepoint. Use unique sample names.')
  for (pop in unique(x$population)) {
    part <- metadata[metadata$population == pop, ]
    if (anyNA(part$time) && !all(is.na(part$time)))
      studio_error('Supply timepoints for every sample in ', pop, ', or leave all its timepoints blank.')
    if (anyDuplicated(part$time[!is.na(part$time)]))
      studio_error('Two samples share a timepoint in ', pop, '. Give independent replicate series different population labels.')
    if (sum(x$counts[x$population == pop]) > .Machine$integer.max)
      studio_error('A population exceeds the current engine’s 32-bit count capacity. Split independent populations before clustering.')
  }
  # Preserve sample identities while collapsing exact duplicate input rows.
  agg <- studio_sum(x$counts, x[c('sample', 'barcode')])
  names(agg)[3] <- 'counts'
  idx <- match(agg$sample, metadata$sample)
  agg$population <- metadata$population[idx]
  agg$time <- metadata$time[idx]
  agg[c('sample', 'population', 'time', 'barcode', 'counts')]
}

studio_import <- function(files, metadata_file = NULL) {
  if (is.null(files) || !nrow(files)) studio_error('Choose at least one barcode-count table.')
  parts <- lapply(seq_len(nrow(files)), function(i) {
    x <- studio_read_table(files$datapath[i], files$name[i])
    if (!'sample' %in% names(x)) x$sample <- sub('[.](csv|tsv|txt)([.]gz)?$', '', files$name[i], ignore.case = TRUE)
    if (!'population' %in% names(x)) x$population <- 'Library 1'
    if (!'time' %in% names(x)) x$time <- NA_real_
    if (!all(c('barcode','counts') %in% names(x))) studio_error('Missing barcode or counts in ', files$name[i], '.')
    x[c('sample','population','time','barcode','counts')]
  })
  x <- do.call(rbind, parts)
  if (!is.null(metadata_file)) x <- studio_metadata(x, metadata_file)
  studio_validate(x)
}

studio_metadata <- function(x, path) {
  m <- studio_read_table(path)
  if (!all(c('sample','time','population') %in% names(m)))
    studio_error('The metadata CSV needs sample, time and population columns.')
  m$sample <- trimws(m$sample)
  if (anyDuplicated(m$sample)) studio_error('Metadata contains duplicate sample names.')
  idx <- match(x$sample, m$sample)
  if (anyNA(idx)) studio_error('Metadata is missing samples: ', paste(unique(x$sample[is.na(idx)]), collapse = ', '))
  x$time <- m$time[idx]; x$population <- m$population[idx]
  x
}

studio_settings <- function(method = 'lv', distance = 3, merge_ratio = 20,
                            error_rate = .005, tie_break = 'sequence', indel_model = 'none') {
  if (!method %in% c('lv','hamming') || !tie_break %in% c('sequence','support') ||
      !indel_model %in% c('none','poisson')) studio_error('Unknown clustering option.')
  if (length(distance) != 1 || !is.finite(distance) || distance < 1 || distance > 5 || distance != floor(distance))
    studio_error('Maximum distance must be a whole number from 1 to 5.')
  if (length(merge_ratio) != 1 || !is.finite(merge_ratio) || merge_ratio < 1 || merge_ratio > 1000)
    studio_error('Merge ratio must be between 1 and 1,000.')
  if (length(error_rate) != 1 || !is.finite(error_rate) || error_rate <= 0 || error_rate >= .25)
    studio_error('Error rate must be greater than 0 and less than 0.25.')
  if (method == 'hamming' && indel_model != 'none') studio_error('The Poisson option requires LV.')
  list(method = method, distance = distance, merge_ratio = merge_ratio,
       error_rate = error_rate, tie_break = tie_break, indel_model = indel_model)
}

studio_progress <- function(directory, message) {
  target <- file.path(directory, 'progress.txt')
  tmp <- paste0(target, '.tmp')
  writeLines(message, tmp); file.rename(tmp, target)
  invisible(NULL)
}

studio_cluster <- function(data, settings, directory, provenance = list()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  data <- studio_validate(data)
  settings <- do.call(studio_settings, settings)
  if (settings$method == 'hamming' &&
      (any(grepl('N', data$barcode)) || any(nchar(data$barcode) > 32) ||
       any(vapply(split(nchar(data$barcode), data$population), function(z) length(unique(z)) != 1L, logical(1)))))
    studio_error('Use LV for variable-length barcodes, N-containing sequences, or sequences longer than 32 bases. Hamming in Studio is restricted to fixed-length A/C/G/T libraries.')
  clock <- proc.time()[['elapsed']]
  populations <- unique(data$population)
  outputs <- memberships <- trajectories <- stats <- list()
  for (i in seq_along(populations)) {
    pop <- populations[i]
    studio_progress(directory, paste('Clustering population', i, 'of', length(populations)))
    part <- data[data$population == pop, ]
    pooled <- studio_sum(part$counts, list(barcode = part$barcode)); names(pooled)[2] <- 'counts'
    result <- do.call(barbac::super_cluster2, c(list(input_path = pooled, verbose = FALSE), settings))
    ids <- sprintf('P%02d-C%06d', i, seq_len(nrow(result)))
    members <- data.frame(population = pop, cluster_id = rep(ids, lengths(result$all_barcodes)),
      central_barcode = rep(result$central_barcode, lengths(result$all_barcodes)),
      barcode = unlist(result$all_barcodes, use.names = FALSE),
      pooled_counts = as.numeric(unlist(result$all_counts, use.names = FALSE)))
    if (anyDuplicated(members$barcode) || !setequal(members$barcode, pooled$barcode))
      studio_error('Clustering membership failed completeness checks.')
    if (sum(result$sum_counts) != sum(part$counts)) studio_error('Clustering did not conserve counts.')
    rows <- match(part$barcode, members$barcode)
    mapped <- data.frame(sample = part$sample, cluster_id = members$cluster_id[rows], counts = part$counts)
    ts <- studio_sum(mapped$counts, mapped[c('sample', 'cluster_id')]); names(ts)[3] <- 'counts'
    meta <- unique(part[c('sample','time')])
    ts$time <- meta$time[match(ts$sample, meta$sample)]
    ts$population <- pop
    ts$central_barcode <- members$central_barcode[match(ts$cluster_id, members$cluster_id)]
    input_totals <- tapply(part$counts, part$sample, sum)
    output_totals <- tapply(ts$counts, ts$sample, sum)
    if (!identical(input_totals, output_totals)) studio_error('Sample-level counts did not reconcile.')
    ts$frequency <- ts$counts / as.numeric(input_totals[ts$sample])
    outputs[[i]] <- data.frame(population = pop, cluster_id = ids,
      central_barcode = result$central_barcode, member_sequences = lengths(result$all_barcodes),
      sum_counts = as.numeric(result$sum_counts))
    memberships[[i]] <- members
    trajectories[[i]] <- ts
    stats[[i]] <- cbind(population = pop, as.data.frame(barbac::cluster_stats(result, verbose = FALSE)))
  }
  elapsed <- proc.time()[['elapsed']] - clock
  out <- list(input = data, centroids = do.call(rbind, outputs),
    memberships = do.call(rbind, memberships), time_series = do.call(rbind, trajectories),
    stats = do.call(rbind, stats), settings = settings, seconds = elapsed,
    provenance = c(list(created_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
      package_version = as.character(utils::packageVersion('barbac')),
      build_id = barbac:::barbac_build_id(),
      source_sha256 = Sys.getenv('BARBAC_STUDIO_SOURCE_SHA256', unset = NA_character_),
      input_sha256 = digest::digest(data, algo = 'sha256'),
      synthetic = isTRUE(provenance$synthetic),
      pooling = 'All supplied timepoints pooled within each population; independent populations clustered separately.',
      normalization = 'All supplied barcode counts per sample. No abundance cutoff; missing lineage counts are zero.',
      clustering_seconds = elapsed), provenance[setdiff(names(provenance), 'synthetic')]))
  studio_progress(directory, 'Preparing downloads')
  for (name in c('centroids', 'memberships', 'time_series', 'stats'))
    readr::write_csv(out[[name]], file.path(directory, paste0(name, '.csv')))
  readr::write_csv(data, file.path(directory, 'extracted_barcodes.csv'))
  jsonlite::write_json(list(settings = settings, provenance = out$provenance),
                      file.path(directory, 'analysis.json'), auto_unbox = TRUE, pretty = TRUE, na = 'null')
  saveRDS(out, file.path(directory, 'analysis.rds'))
  studio_progress(directory, 'Complete')
  out$directory <- directory
  out
}

studio_area <- function(result, population, palette = 'alger', interactive = FALSE, compact = FALSE) {
  ts <- result$time_series[result$time_series$population == population, ]
  if (anyNA(ts$time) || length(unique(ts$time)) < 2)
    studio_error('Add at least two numeric timepoints to this population to show lineage trajectories.')
  n <- length(unique(ts$cluster_id))
  if (n * length(unique(ts$time)) > 2000000)
    studio_error('This lineage grid is too large for the browser view. Download all counts for plotting in R.')
  theme <- ggplot2::theme_minimal(base_size = 12, base_family = 'sans') +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = 'white', colour = NA),
      axis.title = ggplot2::element_text(colour = '#476361'))
  p <- barbac::barbac_ts_area(ts, id_col = 'cluster_id', min_total_count = 0,
    fill_missing = 'zero', include_late = TRUE, palette = palette,
    x_breaks = if(length(unique(ts$time))<=12)sort(unique(ts$time)) else pretty(range(ts$time),n=8),
    theme = theme, interactive = FALSE, x_lab = 'Timepoint', y_lab = 'Fraction of barcode reads')
  if(isFALSE(interactive))return(p)
  # Add hover metadata to the native plot's already prepared data. Keep its
  # stacking, normalized frequencies, scales and palette; widen only the device.
  p$data$.tooltip <- paste0(as.character(p$data$cluster_id),'\nTime ',p$data$time,
                           ': ',sprintf('%.2f%%',100*p$data$.freq))
  p$layers[[1]] <- ggiraph::geom_area_interactive(
    ggplot2::aes(tooltip=.data$.tooltip,data_id=.data$cluster_id),position='stack',colour=NA)
  widget <- ggiraph::girafe(ggobj=p,width_svg=if(compact)6 else 12,height_svg=if(compact)4.8 else 4.6,
    options=list(ggiraph::opts_hover(css='stroke:#163f3a;stroke-width:.8px;'),
                 ggiraph::opts_hover_inv(css='opacity:.4;')))
  ggiraph::girafe_options(widget,ggiraph::opts_sizing(rescale=TRUE,width=1))
}

studio_tool <- function(name) {
  p <- barbac:::.barbac_tool(name)
  if (!nzchar(p)) studio_error('Raw-read extraction needs ', name,
    '. Install the barbac command-line environment or make the tool available on PATH.')
  p
}

studio_command <- function(exe, args, log, stdout = log) {
  status <- system2(exe, vapply(as.character(args), shQuote, character(1)), stdout = stdout, stderr = log)
  if (status != 0L) studio_error(basename(exe), ' failed. Check that your reads and reference match. ',
    paste(tail(readLines(log, warn = FALSE), 5), collapse = ' '))
}

studio_check_fastq <- function(path, expected_ids = NULL) {
  con <- gzfile(path, 'rt'); on.exit(close(con))
  count <- 0L
  repeat {
    z <- readLines(con, n = 40000L, warn = FALSE)
    if (!length(z)) break
    if (length(z) %% 4L) studio_error('FASTQ records must contain four lines. The uploaded file may be truncated.')
    k <- seq.int(1L, length(z), 4L)
    if (any(!startsWith(z[k], '@')) || any(!startsWith(z[k+2L], '+')) ||
        any(nchar(z[k+1L]) != nchar(z[k+3L])) || any(!grepl('^[ACGTNacgtn]+$', z[k+1L])))
      studio_error('Invalid FASTQ record: check identifiers, DNA sequences and quality lengths.')
    count <- count + length(k)
  }
  if (!count) studio_error('The uploaded FASTQ is empty.')
  count
}

studio_extract <- function(job, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  cfg <- job$settings
  minimap <- studio_tool('minimap2'); samtools <- studio_tool('samtools')
  paired <- !is.null(job$r2) && nrow(job$r2) > 0
  if (paired && nrow(job$r1) != nrow(job$r2)) studio_error('Supply the same number of R1 and R2 files, in matching order.')
  sample_labels <- sub('([_-]R?1)?[.](fastq|fq)([.]gz)?$', '', job$r1$name, ignore.case=TRUE)
  if(anyDuplicated(sample_labels))studio_error('FASTQ filenames produce duplicate sample names. Rename the files to identify each sample uniquely.')
  pear <- if (paired) studio_tool('pear') else NULL
  refs <- Biostrings::readDNAStringSet(job$reference)
  if (length(refs) != 1L) studio_error('Supply one barcode-cassette sequence in the reference FASTA.')
  ref_name <- strsplit(names(refs)[1], '[[:space:]]')[[1]][1]
  if (!is.finite(cfg$start) || !is.finite(cfg$end) || cfg$start < 2 || cfg$end <= cfg$start ||
      cfg$end >= Biostrings::width(refs)[1] || cfg$start != floor(cfg$start) || cfg$end != floor(cfg$end))
    studio_error('Use one-based barcode start/end positions inside the reference, with a base on either side.')
  pattern <- NULL
  if (cfg$mode == 'flanks') {
    left <- toupper(trimws(cfg$left)); right <- toupper(trimws(cfg$right))
    if (!grepl('^[ACGT]{4,60}$', left) || !grepl('^[ACGT]{4,60}$', right))
      studio_error('Each exact flank must contain 4–60 A/C/G/T bases.')
    if (!is.finite(cfg$min_length) || !is.finite(cfg$max_length) || cfg$min_length < 1 ||
        cfg$max_length > 250 || cfg$min_length > cfg$max_length ||
        cfg$min_length != floor(cfg$min_length) || cfg$max_length != floor(cfg$max_length))
      studio_error('Check minimum and maximum observed barcode lengths (whole numbers, 1–250).')
    pattern <- paste0(left, '([ACGTN]{', cfg$min_length, ',', cfg$max_length, '})', right)
  }
  output <- receipts <- list()
  for (i in seq_len(nrow(job$r1))) {
    sample_name <- sub('([_-]R?1)?[.](fastq|fq)([.]gz)?$', '', job$r1$name[i], ignore.case = TRUE)
    if (!nzchar(sample_name)) sample_name <- paste0('Sample ', i)
    studio_progress(directory, paste('Checking reads:', sample_name))
    reads <- studio_check_fastq(job$r1$datapath[i])
    sample_dir <- file.path(directory, sprintf('sample-%03d', i)); dir.create(sample_dir)
    fq <- job$r1$datapath[i]
    if (paired) {
      mate_reads <- studio_check_fastq(job$r2$datapath[i])
      if (mate_reads != reads) studio_error('R1 and R2 have different read counts for ', sample_name, '.')
      # Verify mate ordering without loading the whole FASTQ into memory.
      a <- gzfile(fq, 'rt'); b <- gzfile(job$r2$datapath[i], 'rt')
      tryCatch(repeat {
        z <- readLines(a, n = 40000L, warn = FALSE); w <- readLines(b, n = 40000L, warn = FALSE)
        if (!length(z)) break
        k <- seq.int(1L, length(z), 4L)
        clean <- function(v) sub('/[12]$', '', sub('[[:space:]].*$', '', v))
        if (!identical(clean(z[k]), clean(w[k]))) studio_error('Paired FASTQ identifiers are not in matching order.')
      }, finally = {close(a); close(b)})
      studio_progress(directory, paste('Merging paired reads:', sample_name))
      prefix <- file.path(sample_dir, 'merged')
      studio_command(pear, c('-f', fq, '-r', job$r2$datapath[i], '-o', prefix, '-j', 2), file.path(sample_dir,'pear.log'))
      fq <- paste0(prefix, '.assembled.fastq')
      if (!file.exists(fq) || file.info(fq)$size == 0) studio_error('No overlapping pairs merged for ', sample_name, '.')
    }
    processed <- if (paired) studio_check_fastq(fq) else reads
    studio_progress(directory, paste('Mapping reads:', sample_name))
    sam <- file.path(sample_dir,'mapped.sam'); bam <- file.path(sample_dir,'mapped.bam')
    studio_command(minimap, c('-a','-x','sr','--secondary=no','-t','2', job$reference, fq), file.path(sample_dir,'mapping.log'), stdout = sam)
    primary <- file.path(sample_dir,'primary.bam')
    studio_command(samtools,c('view','-b','-F','2304','-o',primary,sam),file.path(sample_dir,'primary.log'))
    studio_command(samtools, c('sort','-@','2','-o',bam,primary), file.path(sample_dir,'sort.log'))
    studio_command(samtools, c('index',bam), file.path(sample_dir,'index.log'))
    mapped <- as.numeric(system2(samtools, c('view','-c','-F','2308',shQuote(bam)), stdout = TRUE))
    if (length(mapped) != 1L || !is.finite(mapped)) studio_error('Could not count mapped primary alignments.')
    studio_progress(directory, paste('Extracting barcodes:', sample_name))
    target <- file.path(sample_dir,'barcodes.csv')
    barbac::barbac_xtr(bam, ref_name = ref_name, start_pos = cfg$start, end_pos = cfg$end,
      output_file = target, flank_pattern = pattern, min_count = 1, verbose = FALSE)
    x <- read.csv(target)
    if (!nrow(x)) studio_error('No barcodes were extracted for ', sample_name, '. Check the reference, coordinates and flanks.')
    x$sample <- sample_name
    output[[i]] <- x[c('sample','barcode','counts')]
    receipts[[i]] <- data.frame(sample = sample_name, input_reads_or_pairs = reads,
      reads_mapped_input = processed, mapped_primary = mapped, extracted_counts = sum(x$counts),
      extraction_mode = cfg$mode, paired = paired)
  }
  x <- do.call(rbind, output)
  if (!is.null(job$metadata)) x <- studio_metadata(x, job$metadata)
  x <- studio_validate(x)
  receipt <- do.call(rbind, receipts)
  if (any(receipt$extracted_counts > receipt$mapped_primary)) studio_error('Extraction counts exceeded primary mapped reads.')
  readr::write_csv(x, file.path(directory,'extracted_barcodes.csv'))
  readr::write_csv(receipt, file.path(directory,'extraction_stats.csv'))
  provenance <- list(extraction = cfg, reference_sha256 = digest::digest(file=job$reference,algo='sha256'),
    files = lapply(seq_len(nrow(job$r1)), function(i) list(name=job$r1$name[i],
      sha256=digest::digest(file=job$r1$datapath[i],algo='sha256'),
      mate_name=if(paired)job$r2$name[i] else NULL,
      mate_sha256=if(paired)digest::digest(file=job$r2$datapath[i],algo='sha256') else NULL)),
    note = 'Single-end reads map directly. Paired-end counts use overlapping PEAR-merged molecules only; unmerged pairs are excluded. No UMI deduplication or study-specific quality filtering is applied.')
  jsonlite::write_json(provenance,file.path(directory,'extraction.json'),auto_unbox=TRUE,pretty=TRUE)
  studio_progress(directory,'Complete')
  list(input = x, extraction_stats = receipt, provenance = provenance, directory = directory)
}
