library(tidyverse)
# adding everything I need to do and updating the github

df2 <- tibble(col1= c("apple","apple","banana","banana","banana"), 
              col2 = c("appl","aple","banan","bananb","banat"),
              count_col1=c(1,1,4,4,4), count_col2=c(3,4,1,1,1))


df2 %>% 
  group_by(col1) %>% 
  mutate(sum_col2=sum(count_col2))  %>% 
  mutate(Max = pmax(count_col1, count_col2)) %>% 
  filter(Max == max(Max)) %>% 
  ungroup %>%
  select(-Max) %>% 
  rowwise() %>% 
  mutate(all_sum= sum(c_across(c(sum_col2, count_col1))))


library(tidyverse)
df3 <- tibble(col1 = c("apple",rep("banana",3)),
              col2 = c("aple", "banan","bananb","banat"), 
              count_col1 = c(1,4,4,4), 
              count_col2 = c(4,1,1,1))
df3

reprex()

df3 %>% 
  group_by(col1) %>% 
  mutate(case_when(count_col2 > count_col1 ~ col1==NA,
                   count_col1 > count_col2 ~ col2==NA ))


reprex()
