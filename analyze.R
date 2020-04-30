library(tidyverse)
library(here)

data <- readRDS("input/data.rds")

commondata <- data$commondata %>%
   mutate(schulbildung_min_12_j = as.factor(schulbildung_min_12_j),
          blutdruckmedis = as.numeric(blutdruckmedis)) 

#write_rds(commondata, path = paste0(here("input"),"/commondata.rds"), compress = "none")

commondata %>% count(delirium, group)

# diffx <- function(x){
#    x <- last(x) - first(x)
# }
# 
# descriptive <- commondata %>%
#    dplyr::select(matches(
#        "delirium|age|del_recall_|m3_|ospan|recall_|vt_number_rt_lb1|vt_number_trail_ie_cal|antihyper"
#     )) %>% 
#    dplyr::select(-matches("_correct_|_false_|_diff")) %>% 
#    group_split(delirium, keep = TRUE) %>% 
#    map(~ .x %>%
#           dplyr::mutate(delirium = as.numeric(delirium)-1) %>% 
#           dplyr::select_if(is.numeric) %>% 
#           skimr::skim() %>% as_tibble() %>% 
#           dplyr::select(matches("_variable|numeric.")) %>% 
#           dplyr::select(1:8) %>% 
#           mutate_if(is.numeric, ~ round(., digits = 2))  
#    ) %>% reduce(bind_rows) %>% arrange(skim_variable) %>% 
#    rename(variable = skim_variable) %>% 
#    group_by(variable) %>% 
#       summarise_if(is.numeric, list(diff = diffx)) %>% 
#    select_all(~ str_replace_all(., "numeric.","")) %>% 
#    filter(!p50_diff == 0) %>% 
#    filter(!variable == "delirium")

selected <- commondata %>%
   dplyr::select(delirium, age, matches("m3_|ospan_|calculated")) %>% 
   dplyr::select(-matches("correct|false|avg|diff|trials|switch")) %>% 
   dplyr::rename(
      pre_M3_ie = m3_pre_ie, 
      post_M3_ie = m3_post_ie,
      pre_OSPAN_ie = ospan_ie_pre,
      post_OSPAN_ie = ospan_ie_post,
      pre_TRAIL_ie = vt_number_trail_ie_calculated_pre,
      post_TRAIL_ie = vt_number_trail_ie_calculated_post
      ) %>% 
   mutate_if(is.numeric, as.integer)

selected %>% write_csv(path = here("input/selected.csv"))

upsampled <- selected %>% 
   mutate(count = ifelse(delirium == 1, 23, 1)) %>% # upscales the data, foctor 23x
   tidyr::uncount(count) %>% 
   sample_frac() # randomizes the order

upsampled %>% write_csv(path = here("input/upsampled.csv"))

   


### scaled data for dimensionality reduction   
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)



#####
library(C50)
set.seed(111)
tree_mod_pre <- C5.0(delirium ~  age + pre_M3_ie + pre_OSPAN_ie + pre_TRAIL_ie, 
                 control = C5.0Control(winnow = TRUE),
                 trials = 3,
                 data = upsampled)

tree_mod_post <- C5.0(delirium ~  age + post_M3_ie + post_OSPAN_ie + post_TRAIL_ie, 
                      #rules = TRUE,
                      data = upsampled)
 





