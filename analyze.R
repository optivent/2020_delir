library(tidyverse)
library(here)

data <- readRDS("input/data.rds")

commondata <- data$commondata %>%
   mutate(schulbildung_min_12_j = as.factor(schulbildung_min_12_j),
          blutdruckmedis = as.numeric(blutdruckmedis)) 

#write_rds(commondata, path = paste0(here("input"),"/commondata.rds"), compress = "none")

commondata %>% write_csv(path = here("input/commondata.csv"))

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

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

Boruta_search <- commondata %>%
   dplyr::rename(
      pre_M3_ie = m3_pre_ie, 
      post_M3_ie = m3_post_ie,
      pre_OSPAN_ie = ospan_ie_pre,
      post_OSPAN_ie = ospan_ie_post,
      pre_TRAIL_ie = vt_number_trail_ie_calculated_pre,
      post_TRAIL_ie = vt_number_trail_ie_calculated_post
   ) %>% 
   mutate_if(is.numeric, as.integer) %>% 
   mutate(count = ifelse(delirium == 1, 2, 1),
          delirium = as.integer(delirium)) %>% 
   uncount(weights = count, .remove = TRUE) %>% 
   select_if(is.numeric) %>% 
   mutate(delirium = as.factor(delirium))

library(Boruta)

Boruta_results <- Boruta(delirium ~ . ,
                         doTrace = 2,
                         pValue = 0.05,
                         maxRuns = 5000, 
                         data = Boruta_search) %>%
   attStats() %>% rownames_to_column() %>% 
   filter(decision == "Confirmed") %>% 
   arrange(desc(meanImp))

Boruta_results %>% filter(str_detect(rowname, 'pre'), meanImp > 3) %>% pull(rowname) %>% paste0(collapse =  " + " ) 

scaled_data <- commondata %>%
   dplyr::rename(
      pre_M3_ie = m3_pre_ie, 
      post_M3_ie = m3_post_ie,
      pre_OSPAN_ie = ospan_ie_pre,
      post_OSPAN_ie = ospan_ie_post,
      pre_TRAIL_ie = vt_number_trail_ie_calculated_pre,
      post_TRAIL_ie = vt_number_trail_ie_calculated_post
   ) %>% 
   mutate_if(is.numeric, as.integer) %>% 
   mutate(count = ifelse(delirium == 1, 6, 1), delirium = as.integer(as.character(delirium))) %>% 
   select_if(is.numeric) %>% 
   mutate(delirium = as.factor(delirium)) %>% 
   mutate_at(vars(-c(delirium, count)), scale2)

pre_data <- scaled_data %>% 
   select(count, delirium, age,
          vt_color_pre_ie , vt_number_rt_lb_last_yes_pre , vt_switch_rt_lb13_yes_pre , pre_OSPAN_ie , vt_switch_trail_ie_calculated_pre , m3_false_pre)

Delir_pre_scallogreg <- glm(delirium ~ age + 
                               vt_color_pre_ie + 
                               vt_number_rt_lb_last_yes_pre + 
                               vt_switch_rt_lb13_yes_pre +
                               pre_OSPAN_ie +
                               vt_switch_trail_ie_calculated_pre +
                               m3_false_pre,
                            weights = count,
                            family = binomial(link = "logit"),
                            data = scaled_data) 


library(MASS)
step.model <- Delir_pre_scallogreg  %>% stepAIC(trace = TRUE)






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
   mutate_if(is.numeric, as.integer) %>% 
   mutate(count = ifelse(delirium == 1, 6, 1)) %>% 
   mutate_at(vars(-c(delirium, count)), log2)


library(gt)
library(stargazer)

Delir_pre_scallogreg <- glm(delirium ~ age + pre_OSPAN_ie + pre_M3_ie + pre_TRAIL_ie,
                        weights = count,
                        family = binomial(link = "logit"),
                         data = selected) 








Delir_post_scallogreg <- glm(delirium ~ age + post_M3_ie + post_OSPAN_ie + post_TRAIL_ie,
                        weights = count,
                        family = binomial(link = "logit"),
                        data = selected)



selected %>% write_csv(path = here("input/selected.csv"))

upsampled <- selected %>% 
   mutate(count = ifelse(delirium == 1, 23, 1)) %>% # upscales the data, foctor 23x
   tidyr::uncount(count) %>% 
   sample_frac()   # randomize the rows


upsampled %>% write_csv(path = here("input/upsampled.csv"))

   


# ### scaled data for dimensionality reduction   
# scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)



#####

# library(C50)
# set.seed(111)
# tree_mod_pre <- C5.0(delirium ~  age + pre_M3_ie + pre_OSPAN_ie + pre_TRAIL_ie, 
#                  #control = C5.0Control(winnow = TRUE),
#                  #trials = 3,
#                  rules = TRUE,
#                  data = upsampled)
# 
# tree_mod_post <- C5.0(delirium ~  age + post_M3_ie + post_OSPAN_ie + post_TRAIL_ie, 
#                       #rules = TRUE,
#                       rules = TRUE,
#                       data = upsampled)

###
library(rpart)
library(rpart.plot)
library(party)

ctree_pre <- party::ctree(delirium ~  age + pre_M3_ie + pre_OSPAN_ie + pre_TRAIL_ie,
                          data = selected %>% 
                             mutate(count = ifelse(delirium == 1, 8 , 1)) %>% # upscales the data, foctor 23x
                             tidyr::uncount(count) %>% 
                             sample_frac())

plot(ctree_pre)
 

#####
# postoperative cognitive disfunction
## himca / ospan / m3 / vt_switch / vt_number / recall delay

# Boruta search for Delirium






## Boruta search
set.seed(111)
library(Boruta)
boruta_search <- Boruta(group ~ . ,
                        doTrace = 2,
                        pValue = 0.05,
                        maxRuns = 5000,
                        data = commondata %>% select(-category,-number)
                        ) 
boruta_results <- boruta_search %>% attStats() %>% rownames_to_column() %>% 
   filter(decision == "Confirmed")


plot(boruta_search, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta_search$ImpHistory),function(i)
   boruta_search$ImpHistory[is.finite(boruta_search$ImpHistory[,i]),i])
names(lz) <- colnames(boruta_search$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(boruta_search$ImpHistory), cex.axis = 0.7)

