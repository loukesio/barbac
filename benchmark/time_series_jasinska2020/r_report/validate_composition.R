#!/usr/bin/env Rscript
# The large-data drawing adapter must preserve native ggplot stack geometry.
kit <- 'benchmark/time_series_jasinska2020'
source(file.path(kit,'load_release.R'));invisible(load_release(kit))
source(file.path(kit,'r_report/composition_helpers.R'))
d <- data.frame(barcode=rep(c('A','B','C'),each=3),
  time=rep(c(0,2,5),3),counts=c(100,30,0,0,40,80,0,30,20))
p <- barbac::barbac_ts_area(d,min_total_count=0,fill_missing='zero',
  palette=viridisLite::plasma(3))
s <- composition_stack(p)
b <- ggplot2::ggplot_build(p)$data[[1]]
stopifnot(identical(s$ids,c('A','B','C')),s$freq['B','0']==0,
  isTRUE(all.equal(unname(s$freq['B',]),c(0,.4,.8))))
for(i in seq_along(s$times)) {
  z <- b[b$x==s$times[i],];z <- z[order(z$group),]
  stopifnot(isTRUE(all.equal(z$ymin,unname(s$lower[,i]),tolerance=1e-12)),
    isTRUE(all.equal(z$ymax,unname(s$upper[,i]),tolerance=1e-12)))
}
# Missing early cells and duplicate input records retain native semantics.
d <- rbind(d[!(d$barcode=='B' & d$time==0),],d[1,])
p <- barbac::barbac_ts_area(d,min_total_count=0,fill_missing='zero',
  palette=viridisLite::plasma(3))
s <- composition_stack(p)
stopifnot(s$freq['B','0']==0,s$freq['A','0']==1,
  all(abs(colSums(s$freq)-1)<1e-12))
cat('Native-stack equivalence, late arrivals, duplicate cells and frequency checks passed.\n')
