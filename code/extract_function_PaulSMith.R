guess I have solved the problem, @LDT!
  
  The idea is to use _all_ variants of the the following patterns to extract the wanted middle string:
  
  ```
"(?<=TGCA).*(?=CTGC)"
"AAAAAA"
"(?<=TGCA).*(?=CTGC)"
"AAAAAA"
"(?<=TGCA).*(?=CTGC)"
"AAAAAA"
"(?<=TGCA).*(?=CTGC)"
"AAAAAA"
"(?<=TGCA).*(?=CTGC)"
"AAAAAA"
"(?<=[A-Z]GCA).*(?=CTGC)"
"AAAAAA"
"(?<=[A-Z]GCA).*(?=[A-Z]TGC)"
...

```

The code below accomplishes that:
  
  ``` r
library(tidyverse)

s <- "AAAAAATGCAAAAAAACTGC"

delim1 <- "TGCA"
delim2 <- "CTGC"

for (i in 0:nchar(delim1))
{
  for (j in 0:nchar(delim2))
  {
    sdelim1 <- str_split(delim1, "") %>% unlist
    sdelim2 <- str_split(delim2, "") %>% unlist
    
    if (i != 0)
      sdelim1[i] <- "[A-Z]"
    
    if (j != 0)
      sdelim2[i] <- "[A-Z]"
    
    pattern <- str_c("(?<=",str_c(sdelim1, collapse = ""),").*(?=",str_c(sdelim2, collapse = ""),")")
    
    print(str_extract(s,pattern))
    
  }
}
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
#> [1] "AAAAAA"
```

I would like to collaborate with you as long as my free time allows me.

