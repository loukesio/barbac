library(tidyverse)

df <- tibble(position=c(100,200,300),
             correction=c("62M89S", 
                     "8M1D55M88S",
                     "1S25M1P36M89S"))

df1 <- df %>% 
  separate(correction, into = str_c("col", 1:5), 
           sep = "(?<=\\D)(?=\\d)", fill = "left", remove = FALSE)

df1

reprex()


df1 %>% 
  mutate(across(starts_with("col")), 
                ~case_when(grepl("*M") | grepl("*S") | grepl("*D")   ~ "",
                           TRUE ~ 0))

         
         
#_______________________________

library(tidyverse)

df1 <- tibble(.x=c(334,335,395),
              .y=c(574,600,466))

df1


df2 <- tibble(id=c(334,335,395,466,574,600),
              fruits=c("apple","banana","ananas","pear","cherry","orange"))

df2

df1 %>%
  mutate(across(everything(), .names="fruits_{col}", ~df2$fruits[match(., df2$id)]))

reprex()

library(reprex)

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



