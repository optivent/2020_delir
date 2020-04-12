clean.it <- function() {
  basic.packages <- c("package:stats","package:graphics",
                      "package:grDevices","package:utils",
                      "package:datasets","package:methods",
                      "package:base")
  package.list <- dplyr::setdiff( search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)] , basic.packages)
  if (length(package.list)>0)  for(package in package.list) detach(package, character.only=TRUE)
  if(!require(pacman))install.packages("pacman"); require(here)
  
  rm(list = dplyr::setdiff( ls(envir = globalenv()),
                            c("clean.it", "path")
  ),
  envir = globalenv())
  #gc() # or sessionInfo()
  
}
####### check data ##########
library(tidyverse)
library(mosaic)

check_struct <- function(df){
  suppressWarnings(
    df <- df %>% purrr::map_dfr(suppressWarnings(mosaic::favstats)) %>%
      cbind(names = names(df)) %>% 
      dplyr::select(names, everything()) %>%  
      mutate(NAs = 100*missing/nrow(df)) %>% 
      cbind(NA_NULL = map_dbl(df, ~ 100*sum(is.na(dplyr::na_if(.,0)))/nrow(df))) %>% 
      dplyr::select(-c(n, missing)) %>% 
      #mutate(imp = sd*(100-perc_NA)) %>% mutate(imp = imp/sum(imp)) %>% 
      mutate_if(is.numeric, ~ round(., digits = 1)) 
  )
  return(df)
}

library(missRanger)
impute_RF <- function(df) {
  require(missRanger)
  df %>% missRanger::missRanger(
    num.trees = 1000, maxiter = 100, pmm.k = 3)
}

####### random-forest-wrapper ######
library(tidyverse)
library(janitor)
library(Boruta)
library(furrr)

scale01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))} # ... allows "na.rm = TRUE"

corr_RF <- function(df, iter){
  set.seed(111)
  plan(multiprocess)
  suppressWarnings(
    
    result <- df %>% names() %>% 
      future_map_dfr(
        ~ attStats(
          Boruta::Boruta(
            formula(
              paste0('`', ., '` ~ ', paste(names(df), collapse = " + "))
            ),
            data = df,
            mcAdj = TRUE, doTrace = 0, holdHistory = TRUE, 
            pValue = min(1/(ncol(df)*nrow(df)), 0.01),
            maxRuns = iter
          )
        ) %>% 
          rownames_to_column() %>% 
          mutate(
            Score = ifelse(decision == "Rejected", 0, 
                           scale01(medianImp * normHits, na.rm = TRUE)),
            Score = na_if(Score, 1)) %>%
          dplyr::select(rowname, Score) %>% 
          pivot_wider(names_from = rowname, values_from = Score),
        .progress = TRUE
      ) %>% mutate(target = colnames(df)) %>% column_to_rownames(var = "target") %>% 
      na_if(0) 
  )
  return(result)
}


matrix_RF <- function(df, sensibility, clusters){
  full_join(enframe(colMeans(df, na.rm = TRUE)) %>% rename(feature_wise = value),
            enframe(rowMeans(df, na.rm = TRUE)) %>% transmute(predictor_wise = value, name = colnames(df))
  ) %>% mutate(product = feature_wise*predictor_wise) %>% 
    top_frac(sensibility, product) %>% pull(name) -> selection
  
  subset(df, rownames(df) %in% selection) %>% subset(select = selection) %>% 
    pheatmap::pheatmap(display_numbers = TRUE, cutree_rows = clust, cutree_cols = clust)
}



library(viridis)
dff <- df %>% 
  mutate(target = names(df)) %>% na_if(0) %>% 
  pivot_longer(-target, names_to = "variable", values_to = "RFcorr") 

ggplot(mutate(dff,
              target = factor(target, levels = rev(colnames(dff))),
              variable = factor(variable, levels = colnames(dff))),
       aes(variable, target, fill = RFcorr %>% replace_na(0))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
  geom_tile(aes(fill = RFcorr)) + 
  geom_text(aes(label = round(RFcorr, 1)), na.rm=TRUE) +
  scale_fill_viridis(name = "RF_imp") +
  coord_equal()



library(scales)
library(viridis)
ggplot(na.omit(dff),
       aes(RFcorr)) +
  geom_histogram(aes(y = stat(count) / sum(count)), bins = 20) +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("Percentage Histogram of non-zero RF_imp)") +
  ylab("Frequency") +
  scale_fill_viridis() +
  theme(plot.title = element_text(hjust = 0.5), 
        text = element_text(size=16)) 

dff %>% na.omit() %>% 
  mutate(bins = cut(RFcorr, breaks = seq(1, 0, by = -0.1))) %>% 
  count(bins) %>% 
  rownames_to_column() %>% arrange(desc(rowname)) %>% rename(nr_of_RFcorr = n) %>% 
  mutate(
    sensibility = 
      paste0(round(cumsum(100*nr_of_RFcorr/nrow(na.omit(dff))),1), " %"),
    bins = str_sub(bins, 2, -2)
  ) %>% 
  separate(bins, c("from_RF_imp", "to_RF_imp"), sep = ",") %>% 
  dplyr::select(-rowname) %>% 
  knitr::kable("html", align = "c") %>% 
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "central")) 


input_procent_of_corr <- 50 # put a slider here :D 

undirected <- rnd_forests %>% 
  mutate(
    RFcorr = replace_na(RFcorr, 0),
    connection = map2_chr( variable, target, ~ toString( sort(c(.x,.y))) )
  ) %>% 
  group_by(connection) %>%
  #summarise(Weight = mean(RFcorr, na.rm = TRUE)) %>% 
  summarise(Weight = min(RFcorr, na.rm = TRUE)) %>% 
  separate(connection, c("variable", "target")) %>% 
  modify_if(~is.numeric(.), ~round(., 3)) %>%
  dplyr::filter( Weight > (max(c(0,NA), na.rm = TRUE)) ) %>% 
  top_frac(ifelse(input_procent_of_corr > 100, 100,
                  ifelse(input_procent_of_corr < 1, NA, input_procent_of_corr)) %>%
             min(c(., 100), na.rm = TRUE)*0.01)


nodes <- intersect( colnames(df),
                    dplyr::select(undirected, variable, target) %>% flatten_chr() %>% unique()
) %>%
  enframe(name = "id", value = "label")

library(tidygraph)
library(ggraph)
library(igraph)
library(extrafont)

ggraph(
  tbl_graph(
    nodes = nodes,
    edges = bind_cols(
      from = left_join(undirected, rename(nodes, variable = label))$id,
      to = left_join(undirected, rename(nodes, target = label))$id,
      weight = undirected$Weight),
    directed = FALSE
  ),
  layout = "linear", circular = TRUE
) + 
  geom_edge_arc(aes(width = weight, colour = weight, alpha = weight)) + 
  scale_edge_width(
    range = c(min(undirected$Weight)*pi, max(undirected$Weight)*pi), guide = FALSE
  ) +
  scale_edge_color_viridis(
    alpha = 0.75,
    begin = min(undirected$Weight)/2,
    end = 1,
    discrete = FALSE, option = "D",direction = 1
  ) +
  scale_edge_alpha(
    range = c(min(undirected$Weight), max(undirected$Weight)),
    guide = FALSE
  ) + 
  geom_node_text(aes(label = label), size = 6, repel = TRUE) +
  labs(edge_colour = "pRF importance") +
  theme_graph() +
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(color = "black", size = 16, vjust = +1),
        legend.text = element_text(color = "black", size = 12)
  ) +
  ggtitle("Graph of permutated random forests") 



########## varrank ################################
library(mlbench)
data("Sonar")
library(tidyverse)

df <- Sonar
target <- "Class"
rm(Sonar)

glimpse(df)

library(VIM)
aggr(df, col=c('navyblue','red'), cex.axis=.7, gap=3,
     numbers=TRUE, sortVars=TRUE, labels=names(df),
     ylab=c("Histogram of missing data","Pattern"))
DataExplorer::plot_missing(df)

library(varrank)
varrank.df<- varrank(data.df = df, method = "estevez",
                     variable.important = df[[target]],
                     discretization.method = "sturges",
                     algorithm = "forward", scheme="mid", verbose = FALSE)

summary(varrank.PimaIndiansDiabetes)

##############
library(purrr)
library(purrrlyr)
library(broom)
library(furrr)
library(htmlTable)


df %>% select_if(is.numeric)  %>%
  future_map_dfr(~tidy(summary(.x))) %>% 
  htmlTable()

mtcars %>%
  select_if(is.numeric)  %>%
  map(~tidy(summary(.x))) %>%  # compute tidy summary of each var
  #map_if(., has_n_col, add_na_col) %>%   # add na-col if missing
  do.call(rbind, .) %>% 
  htmlTable()


library(tidyverse)
library(viridis)
ggplot(dplyr::mutate(na.omit(rnd_forests), bins = cut(RFcorr, breaks = seq(0, 1, by = 0.1))),
       aes(x = RFcorr, fill = bins)) +
  geom_histogram()+
  scale_x_continuous(limits =c(0,1)) +
  scale_fill_manual(
    drop=FALSE,
    values = c(viridis(nlevels(.$bins)/2), 
               viridis(nlevels(.$bins)/2, direction = -1)))








