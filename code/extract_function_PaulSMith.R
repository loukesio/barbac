library(tidyverse)



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
