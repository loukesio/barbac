library(tidyverse)

library(reclin2)

vec <- c("apple","aple","banan","bananan")
stringdist(vec,useNames="strings", method = c("lv"))


# Put vector in data.frame as that is what reclin2 expects
dta <- data.table(str = vec)
dta
p <- pair_minsim(dta, on = "str", deduplication = TRUE,
                 default_comparator = jaro_winkler(), minsim = 0)
p

p[, selection := simsum > 0.8]
deduplicate_equivalence(p, "group", "selection")

my_dist <- function(x, y, ...) stringdist(x, y, method=c("lv"))
p <- pair_minsim(dta, on = "str", deduplication = TRUE,
                 default_comparator = my_dist, minsim = 0)

p

vec <- c("apple","aple","banan","bananan")
stringdistmatrix(vec, useNames = "strings")


library(tidyverse)
library(stringdist)

vec <- c("apple","aple","banan","bananan")
stringdistmatrix(vec, useNames = "strings")




library(lvec)

a <- sample(c("jan", "pier", "tjorres", "korneel"), 10, replace = TRUE)
a
b <- sample(c("jan", "pier", "joris", "corneel"), 10, replace = TRUE)
b

chunks <- lvec::chunk(a, chunk_size = 1E1)
i=seq(chunks[1])
j <- seq_along(b)
res <- expand.grid(i=i, j=j)
res
res$dist <- stringdist(a[res$i], b[res$j])
res$dist
res <- res[res$dist <= 2, ]
res

dist <- lapply(chunks, function(chunk, a, b, threshold, ...) {
  i <- seq(chunk[1], chunk[2])
  j <- seq_along(b)
  res <- expand.grid(i=i, j=j)
  res$dist <- stringdist(a[res$i], b[res$j])
  res <- res[res$dist <= threshold, ]
  res
}, a=a, b=b, threshold = 2)

dist <- do.call(rbind, dist)
dist








library(stringdist)
lookup <- c('Dog', 'Cat', 'Bear')
data <- c('Do g', 'Do gg', 'Caat')
d.matrix <- stringdistmatrix(a = lookup, b = data, useNames="strings",method="cosine")
d.matrix
#list of minimun cosine.s
cosines<-apply(d.matrix, 2, min)
cosines
#return list of the row number of the minimum value
minlist<-apply(d.matrix, 2, which.min) 
#return list of matching values
matchwith<-lookup[minlist]

#final answer
answer<-data.frame(data, matchwith, cosines)



