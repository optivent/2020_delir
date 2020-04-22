
library(tidyverse)
library(here)

data <- readRDS("input/data.rds")

commondata <- data$commondata %>%
   mutate(schulbildung_min_12_j = as.factor(schulbildung_min_12_j),
          blutdruckmedis = as.numeric(blutdruckmedis)) 

write_rds(commondata, path = paste0(here("input"),"/commondata.rds"), compress = "none")

library(readr)

commondata %>% select(
   number, group, delirium,
   ospan_ie_pre, ospan_ie_post, ospan_ie_diff,
   ospan_avg_math_rt_pre, ospan_avg_math_rt_post, ospan_avg_math_rt_diff,
   vt_number_rt_lb_last_yes_pre, vt_number_rt_lb_last_yes_post, vt_number_rt_lb_last_yes_diff,
   m3_pre_ie, m3_post_ie, m3_diff_ie,
   trail_b_trail_a_pre_ie, trail_b_trail_a_post_ie, trail_b_trail_a_diff_ie,
   recall_pre_ie, recall_post_ie, recall_diff_ie,
   del_recall_pre_ie, del_recall_post_ie, del_recall_diff_ie,
   vt_number_trail_ie_calculated_pre, vt_number_trail_ie_calculated_post, vt_number_trail_ie_calculated_diff
) %>% write_csv(path = paste0(here("input"),"/delirium_target.csv"))


delirium_targetOP <- commondata %>% 
   filter(group == "Op") %>% 
   select(
   number, group, delirium,
   ospan_ie_pre, ospan_ie_post, ospan_ie_diff,
   ospan_avg_math_rt_pre, ospan_avg_math_rt_post, ospan_avg_math_rt_diff,
   vt_number_rt_lb_last_yes_pre, vt_number_rt_lb_last_yes_post, vt_number_rt_lb_last_yes_diff,
   m3_pre_ie, m3_post_ie, m3_diff_ie,
   trail_b_trail_a_pre_ie, trail_b_trail_a_post_ie, trail_b_trail_a_diff_ie,
   recall_pre_ie, recall_post_ie, recall_diff_ie,
   del_recall_pre_ie, del_recall_post_ie, del_recall_diff_ie,
   vt_number_trail_ie_calculated_pre, vt_number_trail_ie_calculated_post, vt_number_trail_ie_calculated_diff
) 

delirium_targetOP %>% write_csv(path = paste0(here("input"),"/delirium_targetOP.csv"))

deliriumOP <- commondata %>% filter(group == "Op") %>% 
                  mutate(delirium = as.numeric(delirium)) %>% 
                  select_if(is.numeric) %>% 
                  mutate(delirium = as.factor(delirium))

full.model <- glm(delirium ~.,
                  data = deliriumOP,
                  family = binomial)


delirium_selected <- commondata %>% filter(group == "Op") %>% 
   mutate(
      ospan_false_math_diff = ospan_false_math_pre - ospan_false_math_post,
      vt_number_lb_last_diff = vt_number_lb_last_pre - vt_number_lb_last_post,
      m3_correct_diff = m3_correct_pre - m3_correct_post,
      m3_diff_ie = m3_pre_ie - m3_post_ie,
      recall_diff_ie = recall_pre_ie - recall_post_ie   
   ) %>% 
   dplyr::select(
      delirium, 
      bpi_5, m3_diff_ie, ospan_false_math_diff, vt_number_lb_last_diff,
      age, m3_correct_diff, recall_diff_ie
   ) 

delirium_selected %>% write_csv(path = paste0(here("input"),"/delirium_selectedOP.csv"))



delirium_selected_conserv <- commondata %>% 
   filter(group == "VG") %>% 
   mutate(
      ospan_false_math_diff = ospan_false_math_pre - ospan_false_math_post,
      vt_number_lb_last_diff = vt_number_lb_last_pre - vt_number_lb_last_post,
      m3_correct_diff = m3_correct_pre - m3_correct_post,
      m3_diff_ie = m3_pre_ie - m3_post_ie
   ) %>% 
   dplyr::select(
      delirium, 
      bpi_5, m3_diff_ie, ospan_false_math_diff, vt_number_lb_last_diff,
      age, m3_correct_diff
   ) 

delirium_selected_conserv %>%
   write_csv(path = paste0(here("input"),"/delirium_selected_conserv.csv"))

######################################################
commondata %>% count(delirium, group)

diffx <- function(x){
   x <- last(x) - first(x)
}

descriptive <- commondata %>%
   dplyr::select(matches(
       "delirium|age|del_recall_|m3_|ospan|recall_|vt_number_rt_lb1|vt_number_trail_ie_cal|antihyper"
    )) %>% 
   dplyr::select(-matches("_correct_|_false_|_diff")) %>% 
   group_split(delirium, keep = TRUE) %>% 
   map(~ .x %>%
          dplyr::mutate(delirium = as.numeric(delirium)-1) %>% 
          dplyr::select_if(is.numeric) %>% 
          skimr::skim() %>% as_tibble() %>% 
          dplyr::select(matches("_variable|numeric.")) %>% 
          dplyr::select(1:8) %>% 
          mutate_if(is.numeric, ~ round(., digits = 2))  
   ) %>% reduce(bind_rows) %>% arrange(skim_variable) %>% 
   rename(variable = skim_variable) %>% 
   group_by(variable) %>% 
      summarise_if(is.numeric, list(diff = diffx)) %>% 
   select_all(~ str_replace_all(., "numeric.","")) %>% 
   filter(!p50_diff == 0) %>% 
   filter(!variable == "delirium")

cols_to_select <- descriptive$variable

commondata_filtered <- commondata %>% 
   mutate_at(vars(gender), ~ ifelse(. == "m", 1, 0)) %>% 
   mutate_at(vars(group), ~ ifelse(. == "VG", 0 , 1)) %>% 
   dplyr::select(
      number,
      delirium, group, gender, schulbildung_min_12_j,
      lachs_1:lachs_3, lachs_5:lachs_15,
      parkinson:ppi,
      age,
      del_recall_post_ie, del_recall_pre_ie,
      m3_avg_rt_pre, m3_pre_ie, m3_trials_pre,
      ospan_avg_math_rt_post, ospan_avg_math_rt_pre, ospan_ie_post, ospan_ie_pre,
      recall_post_ie, summe_antihypertensiva,
      vt_number_rt_lb1_yes_post, vt_number_rt_lb1_yes_pre, vt_number_rt_lb13_yes_pre, vt_number_trail_ie_calculated_pre
      ) %>% 
   mutate_all(~ as.numeric(as.character(.)))

### scaled data for dimensionality reduction   
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

commondata_ML <- commondata_filtered %>% 
   #mutate_at(vars(age:vt_number_trail_ie_calculated_pre), scale2, na.rm = TRUE) %>% 
   mutate(count = ifelse(delirium == 1, 23, 1)) %>% 
   tidyr::uncount(count)
   
commondata_ML %>% write_csv(path = paste0(here("input"),"/commondata_ML.csv"))



############ autoxgboost

restrained.task = makeClassifTask(data = boruta_restrained, target = "delirium")
ctrl = makeMBOControl()
ctrl = setMBOControlTermination(ctrl, iters = 1000L) 
res = autoxgboost(restrained.task, control = ctrl, tune.threshold = FALSE)
res



