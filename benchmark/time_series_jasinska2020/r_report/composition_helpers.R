# Full lineage composition: one polygon per barcode, with no display grouping.
# barbac_ts_area supplies the frequencies, fixed factor order, palette and theme.
# Drawing the completed stack as one vectorised grid grob avoids ggplot's
# per-group interpolation overhead for ~200,000 lineages per population.
composition_stack <- function(plot) {
  d <- plot$data
  ids <- levels(d$barcode)
  times <- sort(unique(d$time))
  stopifnot(nrow(d)==length(ids)*length(times),!anyDuplicated(d[c('barcode','time')]))
  freq <- matrix(0,length(ids),length(times),dimnames=list(ids,as.character(times)))
  freq[cbind(as.integer(d$barcode),match(d$time,times))] <- d$.freq
  stopifnot(all(is.finite(freq)),all(freq>=0),all(abs(colSums(freq)-1)<1e-10))
  # ggplot's default stack places the last factor level at the bottom.
  upper <- apply(freq[nrow(freq):1,,drop=FALSE],2,cumsum)[nrow(freq):1,,drop=FALSE]
  lower <- upper-freq
  list(ids=ids,times=times,freq=freq,lower=lower,upper=upper)
}

draw_composition_stack <- function(plot,stack) {
  n <- length(stack$ids); nt <- length(stack$times)
  x <- (stack$times-min(stack$times))/diff(range(stack$times))
  colours <- plot$scales$get_scales('fill')$palette(n)
  stopifnot(identical(names(colours),stack$ids))
  polygons <- grid::polygonGrob(
    x=rep(c(x,rev(x)),each=n),
    y=as.vector(cbind(stack$upper,stack$lower[,nt:1,drop=FALSE])),
    id=rep(seq_len(n),2*nt),default.units='npc',
    gp=grid::gpar(fill=unname(colours),col=NA))
  # Reuse native barbac scales/labels/theme; replace only polygon rendering.
  plot$data <- data.frame(time=range(stack$times),.freq=c(0,1))
  plot$mapping <- ggplot2::aes(x=time,y=.freq)
  plot$layers <- list()
  plot+ggplot2::geom_blank()+ggplot2::annotation_custom(polygons)
}

make_composition <- function(series,colours,here) {
  m <- series$count_matrix
  relative <- file.path('figures',paste0(series$manifest$well[1],'_composition.png'))
  cache <- file.path(here,'../generated/report_diagnostics',
    paste0(series$manifest$well[1],'_composition.rds'))
  signature <- digest::digest(list(counts=m,passages=series$passages,colours=colours,
    source=digest::digest(file=file.path(here,'composition_helpers.R'),algo='sha256'),
    barbac_source=digest::digest(file=file.path(here,'../../../R/12_barbac_ts_area.R'),algo='sha256')),
    algo='sha256')
  if(file.exists(cache) && file.exists(file.path(here,relative))) {
    old <- readRDS(cache)
    if(identical(signature,old$signature) && identical(old$receipt$sha256,
      digest::digest(file=file.path(here,relative),algo='sha256'))) return(old$receipt)
  }
  input <- data.frame(barcode=rep(rownames(m),ncol(m)),
    time=rep(series$passages,each=nrow(m)),counts=as.vector(m))
  ids <- sort(rownames(m))
  title <- sub('No drug r','No antibiotic · replicate ',
    sub('Low CMP r','Chloramphenicol · replicate ',series$population,fixed=TRUE),fixed=TRUE)
  message('Preparing every barcode: ',series$population,' (',nrow(m),' lineages)')
  p <- barbac::barbac_ts_area(input,min_total_count=0,include_late=TRUE,
    fill_missing='zero',time_zero_shift=FALSE,palette=unname(colours[ids]),
    show_legend=FALSE,x_breaks=c(0,6,12,18,24,30),x_lab='Passage (culture transfers)',
    y_lab='Barcode frequency among extracted reads',title=title,
    theme=barbac::theme_barbac(base_size=12,family='Arial')+
      ggplot2::theme(plot.margin=ggplot2::margin(10,25,12,14)))+
    ggplot2::labs(subtitle=paste(format(nrow(m),big.mark=',',trim=TRUE),
      'individual barcode bands · plasma palette · all extracted barcode reads'))
  stack <- composition_stack(p)
  expected <- sweep(m[stack$ids,match(stack$times,series$passages),drop=FALSE],2,
    colSums(m)[match(stack$times,series$passages)],'/')
  stopifnot(identical(stack$ids,ids),isTRUE(all.equal(unname(stack$freq),unname(expected),tolerance=1e-12)),
    all(abs((stack$upper-stack$lower)-expected)<1e-12),
    identical(unname(p$scales$get_scales('fill')$palette(length(ids))),unname(colours[ids])))
  ggplot2::ggsave(file.path(here,relative),draw_composition_stack(p,stack),
    width=11,height=6,dpi=220,device=ragg::agg_png)
  message('Saved ',relative)
  receipt <- list(path=relative,population=series$population,lineages=nrow(m),
    timepoints=ncol(m),input_cells=length(m),polygon_count=length(ids),
    all_barcodes_individual=TRUE,min_total_count=0,fill_missing='zero',
    denominator='All extracted barcode reads',palette='plasma',
    every_frequency_and_band_width_verified=TRUE,
    sha256=digest::digest(file=file.path(here,relative),algo='sha256'))
  saveRDS(list(signature=signature,receipt=receipt),cache)
  receipt
}
