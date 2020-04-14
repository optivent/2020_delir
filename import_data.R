#### prep steps ####
# library(usethis)
# ?use_github
# edit_r_environ()
# use_github(protocol = 'https', auth_token = Sys.getenv("GITHUB_PAT"))

#### import the data ####
library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(skimr)

import <- read_excel("input/DatensatzExcel10.04..xlsx") %>% 
   janitor::clean_names()

import_demograph <- import %>% select(number:category) %>% 
   mutate_at(vars(category), ~ str_replace_all(., "Nase/NNH","NaseNNH")) %>% 
   mutate_at(vars(gender, delirium, group, category), as.factor) 

import_lachs <- import %>% dplyr::select(number, matches("lachs")) %>% 
   mutate_at(vars(-lachs_summe, -number), as.factor) 

import_himca <- import %>%
   dplyr::select(number, matches("himca|schulbildung")) %>% 
   mutate_at(vars(schulbildung_min_12_j),
             ~ ifelse(
                is.na(.) | . %in% c("0", "1"), . , 0
                ) 
             ) %>% 
   mutate_at(vars(-number), as.integer) %>% 
   mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% 
   rowwise() %>% 
      mutate(himca_summe = sum(
             himca_vis_spa +
             himca_benennen + 
             himca_aufmerksamkeit + 
             himca_sprache + 
             himca_abstraktion +
             himca_erinnerung + 
             himca_orientierung + 
             schulbildung_min_12_j
      )
   ) %>%
   ungroup()

import_phqbpi <- import %>% select(number, phq_2a:bpi_15) %>% 
   rowwise() %>% mutate(phq_2 = sum(phq_2a, phq_2b)) %>% ungroup() %>% 
   select(number, matches("phq"), everything()) %>% 
   mutate_at(vars(bpi_1,bpi_8), as.factor)

import_VEundMeds <- import %>% select(number, parkinson:summe_antihypertensiva) %>% 
   mutate_at(vars(-number, - summe_antihypertensiva), ~ as.factor(.))


#########  Split in operative and conservative    
import_operativ <- import %>%
   select(number, narkosedauer:r_rdia_premed, ponv_premed:bosa_score) %>% 
   filter(is.na(narkosedauer) | narkosedauer %in% c("S","G")) %>% 
   mutate(narkosedauer = coalesce(narkosedauer, schnitt_naht_zeit)) %>% 
   mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% 
   mutate_at(vars(narcosis_type), ~ replace_na(.,"Kombi") %>% as.factor())
   
import_konservativ <- import %>%
   select(number, narkosedauer:r_rdia_premed, ponv_premed:bosa_score) %>% 
   filter(narkosedauer == "K") %>% 
   mutate_all(~ replace_na(., "K")) %>% 
   mutate(number = as.numeric(number))
##############################################

import_ala <- import %>% select(number, m3_pre_ie:trail_b_trail_a_post_ie_2) 

import_ala <- import_ala %>% 
   select(one_of(as_tibble(skim(import_ala)) %>% filter(n_missing < 20) %>% pull(skim_variable) )) %>% 
   mutate_if(is.numeric, ~ ifelse(is.na(.), median(., na.rm = TRUE), .)) %>% 
   select_all(~ str_replace_all(., "_1",""))

diff <- import_ala %>% dplyr::select(number, matches("_pre")) %>%
   pivot_longer(cols = -number,
                values_to = "value_preOP",
                names_to = "name_preOP") %>% 
   mutate(row = row_number()) %>% 
inner_join(
   import_ala %>% dplyr::select(number, matches("_post")) %>%
      pivot_longer(cols = -number,
                   values_to = "value_postOP",
                   names_to = "name_postOP") %>% 
   mutate(row = row_number())
) %>% 
rowwise() %>% 
   mutate(value_diff = value_postOP - value_preOP,
   name_diff = str_replace_all(name_postOP, "post", "diff")       
   ) %>% 
ungroup() %>% 
   select(number, name_diff, value_diff) %>% 
   pivot_wider(id_cols = number, names_from = name_diff, values_from = value_diff)

import_ala <- import_ala %>% full_join(diff) 
import_ala <- import_ala %>% select(sort(tidyselect::peek_vars())) %>%
   select(number, everything())
rm(diff, import)



operativ <- list(
   import_demograph,
   import_operativ,
   import_himca,
   import_lachs, 
   import_phqbpi,
   import_VEundMeds, 
   import_ala) %>% reduce(full_join) %>% 
   filter(group == "Op")
   
conservativ <- list(
   import_demograph,
   import_konservativ,
   import_himca,
   import_lachs, 
   import_phqbpi,
   import_VEundMeds, 
   import_ala) %>% reduce(full_join) %>% 
   filter(group == "VG")
   
commondata <- list(
   import_demograph,
   import_himca,
   import_lachs, 
   import_phqbpi,
   import_VEundMeds, 
   import_ala) %>%
   reduce(full_join)

alldata <- rbind(operativ, conservativ)

alldata %>% select(one_of(setdiff(names(alldata), names(commondata)))) %>% glimpse()

alldata <- alldata %>%
   mutate_at(
      vars(narkosedauer,schnitt_naht_zeit,ponv_premed, tranquillizer_preop, analgesic_awr,
           nausea_awr, antiemetics, rr_reduction_periop, atropin_periop, catecholamines_periop),
      as.factor
   ) %>% 
   mutate_if(is.character, as.numeric) %>% 
   mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .))
# some brute imputation (median) to merge the operative and non-operative

operativ %>% select(one_of(setdiff(names(operativ), names(commondata)))) %>% glimpse()

operativ <- operativ %>% 
   mutate_at(
      vars(narkosedauer,schnitt_naht_zeit,ponv_premed, tranquillizer_preop, analgesic_awr,
           nausea_awr, antiemetics, rr_reduction_periop, atropin_periop, catecholamines_periop),
      as.factor
   ) %>% 
   mutate_if(is.character, as.numeric) %>% 
   mutate_all(~ ifelse(is.na(.), median(., na.rm = TRUE), .))



rm(list=ls(pattern="import"))

delir <- list(alldata, commondata, conservativ, operativ)
delir <- delir %>% map(~ .x %>% mutate_if(is.numeric, as.double))

names(delir) <- c("alldata","commondata", "conservativ", "operativ")

write_rds(delir, path = here("input/data.rds"), "xz", compression = 9L)

rm(alldata, commondata, conservativ, operativ)

