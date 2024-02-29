nqnts=10000
#sdf = (nqnts)/((nqnts)/0.0001)
#sdf

library(tidyverse)
library(foreach)

o3 =
foreach(nqnts = c(10,100,1000,1000), .combine = "rbind")%do%{
  rnorm(nqnts, 0, rexp(1, rate = nqnts*10) ) %>%
    data.frame(s=.) %>% 
    mutate(effect = case_when( s > 0.0001 ~ "large",
                               s <= 0.0001 & s > 0.00001 ~ "mid",
                               s <= 0.00001 ~ "small"
    )) %>% .$effect %>% 
    table() %>% prop.table()*100 -> o
  
  data.frame(nqnts, o) -> o2
  names(o2) = c("nqnts", "eff", "Freq")
  return(o2)
}
