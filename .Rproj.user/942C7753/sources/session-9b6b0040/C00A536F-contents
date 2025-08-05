#!/usr/bin Rscript
# LOAD LIBRARIES ----
suppressPackageStartupMessages({
   library(DT)
   library(DataExplorer)
   library(MASS)
   library(bslib)
   library(data.table)
   library(dplyr)
   library(foreach)
   library(immunotation)
   library(modeldata)
   library(openxlsx)
   library(plyr)
   library(shiny)
   library(shinyBS)
   library(shinyalert)
   library(shinycssloaders)
   library(shinydashboard)
   library(shinyjs)
   library(shinythemes)
   library(stringr)
   library(tidyverse)
   library(vroom)
})

# Set global options
options(sass.cache = FALSE)

# Increase the maximum request size limit to 50 MB
options(shiny.maxRequestSize = 50 * 1024^2) # Set to 30 megabytes (MB)

# LOAD DATASETS ----
outputDir <- "responses"
load("data/organ.match.rda")
load("data/cwd.alleles.rda")
load("data/check.list.known.rda")
load("data/ep_ABC.rda")
load("data/ep_D_all.rda")
load("data/exceptions_df.rda")

suppressMessages(
   organ.match.melt <- organ.match %>%
      pivot_longer(
         data = .,
         values_to = "Candidate_antigen",
         names_to = "Locus",
         cols = everything()
      ) %>%
      na.omit() %>%
      distinct() %>%
      mutate(Candidate_antigen = as.character(lapply(Candidate_antigen, function(x) {
         unlist(strsplit(x, "[/]"))[1]
      }))) %>%
      mutate(Candidate_antigen = stringr::str_trim(as.character(
         paste0(Locus, Candidate_antigen)
      ))) %>%
      filter(!grepl("N$", Candidate_antigen)) %>%
      mutate(Allele.1.field = as.character(lapply(strsplit(Candidate_antigen, ":"), function(x) {
         ifelse(length(x) >= 1, x[[1]], "NA")
      })))
)

colnames(ep_ABC) <- paste(sep = "_", as.character(unlist(ep_ABC[2, ])), as.character(unlist(ep_ABC[1, ])))
ep_ABC <- ep_ABC[-c(1:3), ]
ep_ABC <- ep_ABC[, -c(1)]

ep_ABC <- ep_ABC[, !duplicated(names(ep_ABC))] ## added for 2024 data
ep_ABC <- ep_ABC %>%
   # dplyr::rename("Allele" = "_") %>%## 2023
   dplyr::rename("Allele" = "NA_NA") %>% ## 2024
   subset(., select = !grepl("NA", names(.)))

## Check empty columns
colSums(is.na(ep_ABC) | ep_ABC == "")
nrow(ep_ABC) # Num of rows

empty_columns <- colSums(is.na(ep_ABC) |
   ep_ABC == "") == nrow(ep_ABC)
ep_ABC.complete <- ep_ABC[, !empty_columns]
ep_ABC.complete.long <- ep_ABC.complete %>%
   pivot_longer(.,
      names_to = "exon.position",
      values_to = "eplet",
      cols = 2:ncol(.)
   ) %>%
   na.omit() %>%
   filter(eplet != "") %>%
   distinct()


colnames(ep_D_all) <- paste(sep = "_", as.character(unlist(ep_D_all[2, ])), as.character(unlist(ep_D_all[1, ])))
ep_D_all <- ep_D_all[-c(1:3), ]
ep_D_all <- ep_D_all[, -c(1)]

ep_D_all <- ep_D_all[, !duplicated(names(ep_D_all))] ## added for 2024 data

ep_D_all <- ep_D_all %>%
   # dplyr::rename("Allele" = "_") %>%
   dplyr::rename("Allele" = "NA_NA") %>%
   subset(., select = !grepl("NA", names(.)))

## Check empty columns
colSums(is.na(ep_D_all) | ep_D_all == "")
nrow(ep_D_all) # Num of rows

empty_columns <- colSums(is.na(ep_D_all) |
   ep_D_all == "") == nrow(ep_D_all)
ep_D_all.complete <- ep_D_all[, !empty_columns]
ep_D_all.complete.long <- ep_D_all.complete %>%
   pivot_longer(.,
      names_to = "exon.position",
      values_to = "eplet",
      cols = 2:ncol(.)
   ) %>%
   na.omit() %>%
   filter(eplet != "") %>%
   distinct()

antibody_analysis.df <- rbind(ep_D_all.complete.long, ep_ABC.complete.long)

fx <- function(x) {
   paste(x, collapse = ", ")
}
fx.eplet <- function(x) {
   paste(x, collapse = "; ")
}
remove_brackets <- function(x) {
   as.character(lapply(x, function(y) {
      as.character(gsub("\\s*\\([^\\)]+\\)", "", y))
   }))
}

process_cwd_data <- function(data) {
   df <- data %>%
      as.data.frame() %>%
      mutate(key = paste(PatientID, MFI, CWD_2.0, sep = "_")) %>%
      InterMineR::simplifyResult(., index_column = "key", values_column = "Candidate_antigen") %>%
      rownames_to_column(var = "key") %>%
      set_names(c("key", "Candidate_antigen")) %>%
      separate(.,
         key,
         sep = "_",
         into = c("PatientID", "MFI", "CWD_2.0")
      ) %>%
      pivot_wider(
         names_from = "CWD_2.0",
         values_from = "Candidate_antigen",
         id_cols = c(PatientID),
         values_fn = fx,
         unused_fn = fx
      )

   if ("C" %in% colnames(df)) {
      df <- df %>%
         dplyr::rename(Common = C)
   } else {
      df$Common <- NA
   }

   if ("NA" %in% colnames(df)) {
      df <- df %>%
         dplyr::rename(Not_CWD = "NA")
   } else {
      df$Not_CWD <- NA
   }

   if ("WD" %in% colnames(df)) {
      df <- df %>%
         dplyr::rename(Well.Defined = WD)
   } else {
      df$Well.Defined <- NA
   }

   df <- as.data.frame(df) %>%
      dplyr::rename(MFI_with_Allele = MFI) %>%
      mutate_all(as.character) %>%
      as_tibble()

   return(df)
}

process_genotypes <- function(typing_data, id_column_name) {
   ids_to_fill <- toupper(c("a", "b", "c", "drb1", "dqb1", "dqa1", "dpb1", "dpa1"))

   typing <- typing_data %>%
      pivot_longer(
         cols = 2:ncol(.),
         names_to = "Locus",
         values_to = "Allele"
      ) %>%
      mutate(Locus.geno = gsub("_1$|_2$", "", Locus)) %>%
      group_by_at(vars({{ id_column_name }})) %>%
      mutate(Allele = if_else(Locus.geno %in% ids_to_fill, zoo::na.locf(Allele), Allele)) %>%
      pivot_wider(
         id_cols = {{ id_column_name }},
         names_from = Locus,
         values_from = Allele
      )

   return(typing)
}

add_space_after_comma <- function(input_string) {
   sorted_unique_elements <- sort(unique(trimws(unlist(
      strsplit(input_string, ",")
   ))))
   output_string <- paste(sorted_unique_elements, collapse = ", ")
   return(output_string)
}
