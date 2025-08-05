#!/usr/bin Rscript
# UA assessmentapp serevr logic.
# This app takes MFI and HLA
server <- function(input, output, session) {
   options(warn = -1)
   output$table_display <- renderDataTable(NULL)
   output$geno_display <- renderDataTable(NULL)
   output$cleaned_table <- renderDataTable(NULL)
   output$organ_match_tb <- renderDataTable(NULL)
   output$cwd_table <- renderDataTable(NULL)
   output$eplets <- renderDataTable(NULL)
   output$eplets_geno <- renderDataTable(NULL)
   output$shortlisted_eplets <- renderDataTable(NULL)
   output$donor_of_interest <- renderDataTable(NULL)

   shinyjs::toggle("progress_div", anim = TRUE)

   observeEvent(input$submit, {
      validate_inputs <- function(input) {
         # Scenario 1: Check compulsory inputs
         if (is.null(input$mfi_upload) || is.null(input$genotype_upload)) {
            shinyalert::shinyalert(
               title = "Error!",
               text = "Please upload both MFI and Genotype files.",
               type = "error"
            )
            return(FALSE) # Stop processing
         }

         # Scenario 2: Check optional input (Donor of Interest)
         if (input$donor_interest_checkbox) {
            if (is.null(input$donor_of_interest_genotype)) {
               shinyalert::shinyalert(
                  title = "Error!",
                  text = "Please upload Donor of Interest Genotype file or uncheck the option to proceed.",
                  type = "error"
               )
               return(FALSE) # Stop processing
            } else {
               shinyalert::shinyalert(
                  title = "Success!",
                  text = "Data submitted successfully!",
                  type = "success"
               )
               return(TRUE) # Proceed with processing
            }
         }

         # If compulsory inputs are provided and no checkbox for Donor of Interest
         shinyalert::shinyalert(
            title = "Success!",
            text = "Data submitted successfully!",
            type = "success"
         )
         return(TRUE) # Proceed with processing
      }

      # Validate inputs before proceeding
      if (!validate_inputs(input)) {
         return() # Exit event if validation fails
      }

      mfi_data <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         ext <- tools::file_ext(tolower(input$mfi_upload$name))
         switch(ext,
            xls = readxl::read_excel(
               input$mfi_upload$datapath,
               sheet = 1,
               col_types = "text"
            ),
            xlsx = readxl::read_excel(
               input$mfi_upload$datapath,
               sheet = 1,
               col_types = "text"
            ),
            csv = vroom::vroom(
               input$mfi_upload$datapath,
               delim = ",",
               show_col_types = FALSE,
               col_types = vroom::cols(.default = "c")
            ),
            tsv = vroom::vroom(
               input$mfi_upload$datapath,
               delim = "\t",
               show_col_types = FALSE,
               col_types = vroom::cols(.default = "c")
            ),
            txt = vroom::vroom(
               input$mfi_upload$datapath,
               delim = "\t",
               show_col_types = FALSE,
               col_types = vroom::cols(.default = "c")
            ),
            stop("Invalid file; Please upload .xls, .xlsx, .csv, .txt, or .tsv")
         )
      })

      p_serology_clean <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         req(input$mfi_cutoff)
         suppressMessages(
            mfi_data() %>%
               dplyr::select(c("PatientID", "Class", "Final Assignment ALLELES with MFIs")) %>%
               dplyr::rename(MFI = `Final Assignment ALLELES with MFIs`) %>%
               na.omit() %>%
               mutate(MFI = str_replace_all(MFI, " \\(", "\\(")) %>%
               arrange(PatientID) %>%
               tidyr::separate_rows(., MFI, sep = ",", convert = T) %>%
               mutate_if(is.character, stringr::str_trim) %>%
               mutate(
                  Num = str_extract(MFI, "(?<=\\()\\d+(\\.\\d+)?(?=\\))"),
                  Threshold = ifelse(as.numeric(Num) >= input$mfi_cutoff, "PASS", "FAIL") # Use input threshold
               ) %>%
               dplyr::select(-Num) %>%
               mutate(
                  Allele = as.character(lapply(strsplit(MFI, "[(]"), function(x) {
                     ifelse(length(x) >= 1, x[[1]], "NA")
                  })),
                  Allele.1.field = as.character(lapply(strsplit(MFI, ":"), function(x) {
                     ifelse(length(x) >= 1, x[[1]], "NA")
                  })),
                  Locus = as.character(lapply(strsplit(MFI, "\\*"), function(x) {
                     return(x[[1]])
                  }))
               ) %>%
               mutate(Allele = trimws(Allele))
         )
      })

      organ_match <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         suppressMessages(
            p_serology_clean() %>%
               plyr::join(., organ.match.melt) %>%
               distinct() %>%
               dplyr::filter(
                  Threshold == "PASS",
                  Allele.1.field != Candidate_antigen
               ) %>%
               dplyr::filter(stringr::str_trim(Allele) != stringr::str_trim(Candidate_antigen)) %>%
               dplyr::left_join(
                  .,
                  (exceptions_df %>%
                     dplyr::rename(Allele = Recipient_Antibody))
               )
         )
      })

      cwd_match <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         suppressMessages(
            organ_match() %>%
               dplyr::left_join(
                  .,
                  (cwd.alleles %>%
                     dplyr::select(1, 3, 4, 6) %>%
                     set_names(c("Locus", "Candidate_antigen", "CWD_1.0", "CWD_2.0")))
               ) %>%
               dplyr::filter(Candidate_antigen != Allele.1.field) %>%
               dplyr::select(
                  PatientID,
                  MFI,
                  Threshold,
                  Allele,
                  Donor_HLA_exception,
                  Locus,
                  Candidate_antigen,
                  CWD_1.0,
                  CWD_2.0
               )
         )
      })

      cwd_match_summary <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         cwd_match() %>% process_cwd_data()
      })

      antibody_analysis <- reactive({
         df1 <- rbind(ep_D_all.complete.long, ep_ABC.complete.long)
         return(df1)
      })

      # Eplets Donor v/s Bead
      eplets_match <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)

         # Start parallel back-end
         cl <- makeCluster(detectCores())
         registerDoParallel(cl)

         antibody_analysis.df <- antibody_analysis()
         cwd <- cwd_match()

         suppressMessages(
            combined.cwd.cleanOrganMatch <- cwd[, 1:7] %>%
               mutate(Organ.match.1st.ALT = ifelse(
                  grepl("/", Candidate_antigen),
                  sub("\\/.*", "", Candidate_antigen),
                  Candidate_antigen
               )) %>%
               dplyr::left_join(
                  .,
                  cwd.alleles %>%
                     dplyr::select(1, 3, 4, 6) %>%
                     set_names(c("Locus", "Candidate_antigen", "CWD_1.0", "CWD_2.0"))
               ) %>%
               distinct(.) %>%
               mutate(
                  OrganMatch.2nd.F = ifelse(
                     !endsWith(Organ.match.1st.ALT, "G"),
                     as.character(lapply(strsplit(Organ.match.1st.ALT, ":"), function(x) {
                        ifelse(length(x) >= 2, paste(x[[1]], x[[2]], sep = ":"), x)
                     })),
                     Organ.match.1st.ALT
                  )
               ) %>%
               as_tibble()
         )

         class1 <- c("A", "B", "C")
         class2 <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")

         class1.alleles <- combined.cwd.cleanOrganMatch %>%
            dplyr::filter(
               .,
               grepl(paste0("^(", paste(class1, collapse = "|"), ")"), Candidate_antigen)
            ) %>%
            dplyr::filter(
               !endsWith(Candidate_antigen, "G"),
               !grepl(":[A-Z]", Candidate_antigen)
            ) %>%
            distinct(Candidate_antigen) %>%
            pull(Candidate_antigen)

         class2.alleles <- combined.cwd.cleanOrganMatch %>%
            dplyr::filter(
               .,
               grepl(paste0("^(", paste(class2, collapse = "|"), ")"), Candidate_antigen)
            ) %>%
            dplyr::filter(
               !endsWith(Candidate_antigen, "G"),
               !grepl(":[A-Z]", Candidate_antigen)
            ) %>%
            distinct(Candidate_antigen) %>%
            pull(Candidate_antigen)

         x <- get_serotypes(class1.alleles, mhc_type = "MHC-I") %>%
            enframe() %>%
            dplyr::rename(MRO_Complex = name, MRO_Serotype = value) %>%
            na.omit(MRO_Serotype) %>%
            mutate(
               OrganMatch.2nd.F = gsub(" protein complex| serotype|HLA-", "", MRO_Complex),
               MRO_Serotype = gsub(" protein complex| serotype|HLA-", "", MRO_Serotype)
            )

         y <- get_serotypes(class2.alleles, mhc_type = "MHC-II") %>%
            enframe() %>%
            dplyr::rename(MRO_Complex = name, MRO_Serotype = value) %>%
            na.omit(MRO_Serotype) %>%
            mutate(OrganMatch.2nd.F = gsub(" protein complex| serotype|HLA-", "", MRO_Complex)) %>%
            separate_rows(OrganMatch.2nd.F, sep = "/") %>%
            mutate(
               MRO_Serotype = gsub(" protein complex| serotype|HLA-", "", MRO_Serotype)
            )

         Serotype <- bind_rows(x, y) %>%
            na.omit(MRO_Serotype) %>%
            distinct() %>%
            as.data.frame() %>%
            InterMineR::simplifyResult(.,
               index_column = "OrganMatch.2nd.F",
               values_column = "MRO_Serotype",
               returnInDataframe = T
            ) %>%
            mutate(simplified_results = sapply(
               strsplit(as.character(simplified_results), ","), function(x) {
                  paste(sort(unique(x)), collapse = ",")
               }
            )) %>%
            dplyr::select(-MRO_Serotype) %>%
            dplyr::rename(MRO_Serotype = simplified_results) %>%
            arrange(OrganMatch.2nd.F, MRO_Serotype) %>%
            group_by(OrganMatch.2nd.F) %>%
            slice(1) %>%
            mutate(MRO_Serotype = as.character(MRO_Serotype))

         combined.cwd.cleanOrganMatch <- combined.cwd.cleanOrganMatch %>%
            dplyr::left_join(., Serotype, relationship = "many-to-many")

         eplets_shared <- data.frame()

         # Global functions
         fx <- function(x) paste(x, collapse = ", ")
         remove_brackets <- function(x) {
            as.character(
               lapply(x, function(y) {
                  as.character(gsub("\\s*\\([^\\)]+\\)", "", y))
               })
            )
         }

         collapse_eplet_list <- function(x) paste(x, collapse = "; ")

         process_cwd_data <- function(mfi_data) {
            df <- mfi_data %>%
               as.data.frame() %>%
               mutate(key = paste(PatientID, MFI, CWD_2.0, sep = "_")) %>%
               InterMineR::simplifyResult(., index_column = "key", values_column = "Candidate_antigen") %>%
               rownames_to_column(var = "key") %>%
               set_names(c("key", "Candidate_antigen")) %>%
               separate(., key, sep = "_", into = c("PatientID", "MFI", "CWD_2.0")) %>%
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

         df <- foreach(i = 1:nrow(combined.cwd.cleanOrganMatch), .combine = rbind) %dopar% {
            x <- combined.cwd.cleanOrganMatch[i, ]
            a <- stringr::str_trim(x[["Allele"]], side = c("both"))
            b <- stringr::str_trim(x[["OrganMatch.2nd.F"]], side = c("both"))

            if (endsWith(b, "G") || a == b) {
               x["shared_eplets(position)"] <- "NOT_CHECKED_identical donor-recipient or G-Group alleles)"
               x["Bead_Specific_Eplet"] <- "NOT_CHECKED_identical donor-recipient or G-Group alleles)"
               x["Eplet_of_Candidate_antigen"] <- "NOT_CHECKED_identical donor-recipient or G-Group alleles)"

               eplets_shared <- rbind(eplets_shared, x)
            }

            suppressMessages(
               if (a != b && !endsWith(b, "G")) {
                  df1 <- dplyr::filter(antibody_analysis.df, Allele == a)
                  df1 <- subset(df1, select = c(2:3))
                  df2 <- dplyr::filter(antibody_analysis.df, Allele == b)
                  df2 <- subset(df2, select = c(2:3))
                  res <- dplyr::full_join(df1, df2)
                  df1.1 <- dplyr::filter(antibody_analysis.df, Allele == a)
                  df2.1 <- dplyr::filter(antibody_analysis.df, Allele == b)

                  res.uniq.p <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
                  res.uniq.p <- dplyr::filter(res.uniq.p, is.na(Allele.y))
                  res.uniq.o <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
                  res.uniq.o <- dplyr::filter(res.uniq.o, is.na(Allele.x))

                  if (nrow(res) >= 2 | nrow(res.uniq.p) >= 2 | nrow(res.uniq.o) >= 2) {
                     res$PatientID <- x[[1]]
                     res.uniq.p$PatientID <- x[[1]]
                     res.uniq.o$PatientID <- x[[1]]
                     a <- ifelse(nrow(res) >= 1,
                        collapse_eplet_list(
                           paste0(res$eplet, "(", res$exon.position, ")")
                        ),
                        "Eplet data missing"
                     )
                     b <- ifelse(nrow(res.uniq.p) >= 1,
                        collapse_eplet_list(
                           paste0(res.uniq.p$eplet, "(", res.uniq.p$exon.position, ")")
                        ),
                        "Eplet data missing"
                     )
                     c <- ifelse(nrow(res.uniq.o) >= 1,
                        collapse_eplet_list(
                           paste0(res.uniq.o$eplet, "(", res.uniq.o$exon.position, ")")
                        ),
                        "Eplet data missing"
                     )

                     # Split the strings in a, b, and c into individual elements
                     a_elements <- strsplit(a, "; ")[[1]]
                     b_elements <- strsplit(b, "; ")[[1]]
                     c_elements <- strsplit(c, "; ")[[1]]

                     # Remove elements of b and c from a
                     filtered_elements <- a_elements[!(a_elements %in% b_elements | a_elements %in% c_elements)]

                     # Join the elements back into a string
                     filtered_a <- paste(filtered_elements, collapse = "; ")
                     x["shared_eplets(position)"] <- filtered_a
                     x["Bead_Specific_Eplet"] <- b
                     x["Eplet_of_Candidate_antigen"] <- c
                     eplets_shared <- rbind(eplets_shared, x)
                  }
               }
            )
         }

         suppressMessages(
            shared_eplets <- distinct(df) %>%
               dplyr::left_join(., (exceptions_df %>% dplyr::rename(Allele = Recipient_Antibody)))
         )

         # Stop parallel backend
         stopCluster(cl)
         return(shared_eplets)
      })
      # ------------------------------END EPLETS DONOR v/s BEAD -----------------------------------

      ## Patient genotype
      patient_geno <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)

         ext <- tools::file_ext(tolower(input$genotype_upload$name))
         genotype <- switch(ext,
            xls = read_excel(
               input$genotype_upload$datapath,
               sheet = 1,
               col_types = "text"
            ),
            xlsx = read_excel(
               input$genotype_upload$datapath,
               sheet = 1,
               col_types = "text"
            ),
            csv = vroom::vroom(
               input$genotype_upload$datapath,
               delim = ",",
               show_col_types = FALSE,
               col_types = cols(.default = "c")
            ),
            tsv = vroom::vroom(
               input$genotype_upload$datapath,
               delim = "\t",
               show_col_types = FALSE,
               col_types = cols(.default = "c")
            ),
            txt = vroom::vroom(
               input$genotype_upload$datapath,
               delim = "\t",
               show_col_types = FALSE,
               col_types = cols(.default = "c")
            ),
            stop("Invalid file; Please upload .xls, .xlsx, .csv, .txt, or .tsv")
         )

         genotype <- genotype %>%
            as_tibble() %>%
            dplyr::select(
               PatientID,
               A_1 = HLA1_A,
               A_2 = HLA2_A,
               B_1 = HLA1_B,
               B_2 = HLA2_B,
               C_1 = HLA1_C,
               C_2 = HLA2_C,
               DRB1_1 = HLA1_DRB1,
               DRB1_2 = HLA2_DRB1,
               DRB3_1 = HLA1_DRB3,
               DRB3_2 = HLA2_DRB3,
               DRB4_1 = HLA1_DRB4,
               DRB4_2 = HLA2_DRB4,
               DRB5_1 = HLA1_DRB5,
               DRB5_2 = HLA2_DRB5,
               DQB1_1 = HLA1_DQB1,
               DQB1_2 = HLA2_DQB1,
               DQA1_1 = HLA1_DQA1,
               DQA1_2 = HLA2_DQA1,
               DPB1_1 = HLA1_DPB1,
               DPB1_2 = HLA2_DPB1,
               DPA1_1 = HLA1_DPA1,
               DPA1_2 = HLA2_DPA1
            ) %>%
            dplyr::filter_all(any_vars(!is.na(.))) %>%
            mutate_at(vars(2:ncol(.)), ~ gsub("\\\\\\*", ",", .)) %>%
            mutate_at(vars(2:ncol(.)), ~ gsub("\\*", "", .)) %>%
            as_tibble()

         genotype_fill <- process_genotypes(genotype, PatientID)
         return(genotype_fill)
      })

      ## Donor of interest genotype ||depends on check box state!!
      donor_geno <- reactive({
         req(input$donor_of_interest_genotype)

         ext <- tools::file_ext(tolower(input$donor_of_interest_genotype$name))
         genotyped <- switch(ext,
            xls = read_excel(
               input$donor_of_interest_genotype$datapath,
               sheet = 1,
               col_types = "text"
            ),
            xlsx = read_excel(
               input$donor_of_interest_genotype$datapath,
               sheet = 1,
               col_types = "text"
            ),
            csv = vroom::vroom(
               input$donor_of_interest_genotype$datapath,
               delim = ",",
               show_col_types = FALSE,
               col_types = cols(.default = "c")
            ),
            tsv = vroom::vroom(
               input$donor_of_interest_genotype$datapath,
               delim = "\t",
               show_col_types = FALSE,
               col_types = cols(.default = "c")
            ),
            txt = vroom::vroom(
               input$donor_of_interest_genotype$datapath,
               delim = "\t",
               show_col_types = FALSE,
               col_types = cols(.default = "c")
            ),
            stop("Invalid file; Please upload .xls, .xlsx, .csv, .txt, or .tsv")
         )

         genotyped <- genotyped %>%
            as_tibble() %>%
            dplyr::select(
               DonorID,
               A_1 = HLA1_A,
               A_2 = HLA2_A,
               B_1 = HLA1_B,
               B_2 = HLA2_B,
               C_1 = HLA1_C,
               C_2 = HLA2_C,
               DRB1_1 = HLA1_DRB1,
               DRB1_2 = HLA2_DRB1,
               DRB3_1 = HLA1_DRB3,
               DRB3_2 = HLA2_DRB3,
               DRB4_1 = HLA1_DRB4,
               DRB4_2 = HLA2_DRB4,
               DRB5_1 = HLA1_DRB5,
               DRB5_2 = HLA2_DRB5,
               DQB1_1 = HLA1_DQB1,
               DQB1_2 = HLA2_DQB1,
               DQA1_1 = HLA1_DQA1,
               DQA1_2 = HLA2_DQA1,
               DPB1_1 = HLA1_DPB1,
               DPB1_2 = HLA2_DPB1,
               DPA1_1 = HLA1_DPA1,
               DPA1_2 = HLA2_DPA1
            ) %>%
            dplyr::filter_all(any_vars(!is.na(.))) %>%
            mutate_at(vars(2:ncol(.)), ~ gsub("\\\\\\*", ",", .)) %>%
            mutate_at(vars(2:ncol(.)), ~ gsub("\\*", "", .)) %>%
            as_tibble()

         genotyped_fill <- process_genotypes(genotyped, DonorID)

         return(genotyped_fill)
      })

      # Eplets donor vs recipient
      eplets_match_geno <- reactive({
         req(input$genotype_upload)
         req(input$mfi_upload)

         # Start parallel backend
         cl <- makeCluster(detectCores())
         registerDoParallel(cl)

         suppressMessages(
            patient_geno.long <- patient_geno() %>%
               pivot_longer(
                  cols = 2:ncol(.),
                  names_to = "Locus_allele",
                  values_to = "Typing"
               ) %>%
               mutate(
                  allele = as.character(
                     lapply(Locus_allele, function(x) {
                        if (grepl("_1$", x, ignore.case = TRUE)) {
                           return("allele_1")
                        }
                        if (grepl("_2$", x, ignore.case = TRUE)) {
                           return("allele_2")
                        }
                     })
                  )
               ) %>%
               mutate(
                  Locus = gsub("_1$|_2$", "", Locus_allele),
                  Typing = ifelse(!is.na(Typing), paste0(Locus, "*", Typing), Typing),
                  Patient_genotype.2f = as.character(
                     lapply(strsplit(Typing, ":"), function(x) {
                        ifelse(length(x) >= 2, paste(x[[1]], x[[2]], sep = ":"), "NA")
                     })
                  )
               )
         )

         shared_eplets <- eplets_match()

         suppressMessages(
            combine_bead_patient_typing <- dplyr::left_join(
               shared_eplets %>% dplyr::select(PatientID, Locus, OrganMatch.2nd.F) %>% distinct(),
               patient_geno.long %>%
                  mutate(PatientID = stringr::str_trim(PatientID, side = c("both"))) %>%
                  dplyr::select(PatientID, allele, Locus, Patient_genotype.2f),
               relationship = "many-to-many"
            ) %>%
               dplyr::filter(!is.na(PatientID)) %>%
               pivot_wider(names_from = allele, values_from = Patient_genotype.2f)
         )


         fx <- function(x) paste(x, collapse = ", ")
         remove_brackets <- function(x) {
            as.character(lapply(x, function(y) {
               as.character(gsub("\\s*\\([^\\)]+\\)", "", y))
            }))
         }

         collapse_eplet_list <- function(x) paste(x, collapse = "; ")

         antibody_analysis.df <- antibody_analysis()

         eplets_shared_geno <- data.frame()

         df.geno_allele1 <- foreach(i = 1:nrow(combine_bead_patient_typing), .combine = rbind) %dopar% {
            ## Get patient and organ match allele from row
            row <- combine_bead_patient_typing[i, ]
            patient.allele1 <- stringr::str_trim(row[["allele_1"]], side = c("both"))
            organ.m.allele <- stringr::str_trim(row[["OrganMatch.2nd.F"]], side = c("both"))

            ## Add 3 new columns to row
            row[["shared_eplets_patient_donor_1(position)"]] <- NA
            row[["Recipient_Specific_Eplet_1"]] <- NA
            row[["Eplet_of_Candidate_antigen_1"]] <- NA

            ## Perform if-else conditions
            suppressMessages(
               if (is.na(patient.allele1) | is.na(organ.m.allele)) {
                  row[["shared_eplets_patient_donor_1(position)"]] <- "No typing result"
                  row[["Recipient_Specific_Eplet_1"]] <- "No typing result"
                  row[["Eplet_of_Candidate_antigen_1"]] <- "No typing result"
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               } else if (endsWith(organ.m.allele, "G")) {
                  row[["shared_eplets_patient_donor_1(position)"]] <- "NOT_CHECKED_G-Group donor alleles"
                  row[["Recipient_Specific_Eplet_1"]] <- "NOT_CHECKED_G-Group donor alleles"
                  row[["Eplet_of_Candidate_antigen_1"]] <- "NOT_CHECKED_G-Group donor alleles"
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               } else if (patient.allele1 == organ.m.allele) {
                  row[["shared_eplets_patient_donor_1(position)"]] <- "NOT_CHECKED_identical donor-recipient alleles"
                  row[["Recipient_Specific_Eplet_1"]] <- "NOT_CHECKED_identical donor-recipient alleles"
                  row[["Eplet_of_Candidate_antigen_1"]] <- "NOT_CHECKED_identical donor-recipient alleles"
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               } else if (!is.na(patient.allele1) | !is.na(organ.m.allele)) {
                  df1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele1)
                  df1 <- subset(df1, select = c(2:3))
                  df2 <- dplyr::filter(antibody_analysis.df, Allele == organ.m.allele)
                  df2 <- subset(df2, select = c(2:3))
                  res <- dplyr::full_join(df1, df2)
                  df1.1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele1)
                  df2.1 <- dplyr::filter(antibody_analysis.df, Allele == organ.m.allele)

                  res.uniq.p <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
                  res.uniq.p <- dplyr::filter(res.uniq.p, is.na(Allele.y))
                  res.uniq.o <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
                  res.uniq.o <- dplyr::filter(res.uniq.o, is.na(Allele.x))

                  if (nrow(res) >= 2 || nrow(res.uniq.p) >= 2 || nrow(res.uniq.o) >= 2) {
                     res$PatientID <- row[["PatientID"]]
                     res.uniq.p$PatientID <- row[["PatientID"]]
                     res.uniq.o$PatientID <- row[["PatientID"]]

                     a <- ifelse(nrow(res) >= 1,
                        collapse_eplet_list(paste0(res$eplet, "(", res$exon.position, ")")),
                        "Eplet data missing"
                     )
                     b <- ifelse(nrow(res.uniq.p) >= 1,
                        collapse_eplet_list(paste0(res.uniq.p$eplet, "(", res.uniq.p$exon.position, ")")),
                        "Eplet data missing"
                     )
                     c <- ifelse(nrow(res.uniq.o) >= 1,
                        collapse_eplet_list(paste0(res.uniq.o$eplet, "(", res.uniq.o$exon.position, ")")),
                        "Eplet data missing"
                     )

                     # Split the strings in a, b, and c into individual elements
                     a_elements <- strsplit(a, "; ")[[1]]
                     b_elements <- strsplit(b, "; ")[[1]]
                     c_elements <- strsplit(c, "; ")[[1]]

                     # Remove elements of b and c from a
                     filtered_elements <- a_elements[!(a_elements %in% b_elements | a_elements %in% c_elements)]

                     # Join the elements back into a string
                     filtered_a <- paste(filtered_elements, collapse = "; ")

                     row[["shared_eplets_patient_donor_1(position)"]] <- filtered_a
                     row[["Recipient_Specific_Eplet_1"]] <- b
                     row[["Eplet_of_Candidate_antigen_1"]] <- c
                     eplets_shared_geno <- rbind(eplets_shared_geno, row)
                  }
               }
            )
         }

         eplets_shared_geno <- data.frame()

         df.geno_allele2 <- foreach(i = 1:nrow(combine_bead_patient_typing), .combine = rbind) %dopar% {
            ## Get patient and organ match allele from row
            row <- combine_bead_patient_typing[i, ]
            patient.allele2 <- stringr::str_trim(row[["allele_2"]], side = c("both"))
            organ.m.allele <- stringr::str_trim(row[["OrganMatch.2nd.F"]], side = c("both"))

            ## Add 3 new columns to row
            row[["shared_eplets_patient_donor_2(position)"]] <- NA
            row[["Recipient_Specific_Eplet_2"]] <- NA
            row[["Eplet_of_Candidate_antigen_2"]] <- NA

            ## Perform if-else conditions
            if (is.na(patient.allele2) | is.na(organ.m.allele)) {
               row[["shared_eplets_patient_donor_2(position)"]] <- "No typing result"
               row[["Recipient_Specific_Eplet_2"]] <- "No typing result"
               row[["Eplet_of_Candidate_antigen_2"]] <- "No typing result"
               eplets_shared_geno <- rbind(eplets_shared_geno, row)
            } else if (endsWith(organ.m.allele, "G")) {
               row[["shared_eplets_patient_donor_2(position)"]] <- "NOT_CHECKED_G-Group donor alleles"
               row[["Recipient_Specific_Eplet_2"]] <- "NOT_CHECKED_G-Group donor alleles"
               row[["Eplet_of_Candidate_antigen_2"]] <- "NOT_CHECKED_G-Group donor alleles"
               eplets_shared_geno <- rbind(eplets_shared_geno, row)
            } else if (patient.allele2 == organ.m.allele) {
               row[["shared_eplets_patient_donor_2(position)"]] <- "NOT_CHECKED_identical donor-recipient alleles"
               row[["Recipient_Specific_Eplet_2"]] <- "NOT_CHECKED_identical donor-recipient alleles"
               row[["Eplet_of_Candidate_antigen_2"]] <- "NOT_CHECKED_identical donor-recipient alleles"
               eplets_shared_geno <- rbind(eplets_shared_geno, row)
            } else if (!is.na(patient.allele2) | !is.na(organ.m.allele)) {
               df1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele2)
               df1 <- subset(df1, select = c(2:3))
               df2 <- dplyr::filter(antibody_analysis.df, Allele == organ.m.allele)
               df2 <- subset(df2, select = c(2:3))
               res <- dplyr::full_join(df1, df2)
               df1.1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele2)
               df2.1 <- dplyr::filter(antibody_analysis.df, Allele == organ.m.allele)

               res.uniq.p <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
               res.uniq.p <- dplyr::filter(res.uniq.p, is.na(Allele.y))
               res.uniq.o <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
               res.uniq.o <- dplyr::filter(res.uniq.o, is.na(Allele.x))

               if (nrow(res) >= 2 || nrow(res.uniq.p) >= 2 || nrow(res.uniq.o) >= 2) {
                  res$PatientID <- row[["PatientID"]]
                  res.uniq.p$PatientID <- row[["PatientID"]]
                  res.uniq.o$PatientID <- row[["PatientID"]]

                  a <- ifelse(nrow(res) >= 1,
                     collapse_eplet_list(paste0(res$eplet, "(", res$exon.position, ")")),
                     "Eplet data missing"
                  )
                  b <- ifelse(nrow(res.uniq.p) >= 1,
                     collapse_eplet_list(paste0(res.uniq.p$eplet, "(", res.uniq.p$exon.position, ")")),
                     "Eplet data missing"
                  )
                  c <- ifelse(nrow(res.uniq.o) >= 1,
                     collapse_eplet_list(paste0(res.uniq.o$eplet, "(", res.uniq.o$exon.position, ")")),
                     "Eplet data missing"
                  )

                  # Split the strings in a, b, and c into individual elements
                  a_elements <- strsplit(a, "; ")[[1]]
                  b_elements <- strsplit(b, "; ")[[1]]
                  c_elements <- strsplit(c, "; ")[[1]]

                  # Remove elements of b and c from a
                  filtered_elements <- a_elements[!(a_elements %in% b_elements | a_elements %in% c_elements)]

                  # Join the elements back into a string
                  filtered_a <- paste(filtered_elements, collapse = "; ")

                  row[["shared_eplets_patient_donor_2(position)"]] <- filtered_a
                  row[["Recipient_Specific_Eplet_2"]] <- b
                  row[["Eplet_of_Candidate_antigen_2"]] <- c
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               }
            }
         }

         # Stop parallel back-end
         stopCluster(cl)

         df.geno_allele1 <- df.geno_allele1 %>%
            distinct() %>%
            arrange(PatientID, Locus, OrganMatch.2nd.F)

         df.geno_allele2 <- df.geno_allele2 %>%
            distinct() %>%
            arrange(PatientID, Locus, OrganMatch.2nd.F)

         combined.geno.eplets <- full_join(df.geno_allele1, df.geno_allele2, relationship = "many-to-many") %>%
            arrange(PatientID, Locus, OrganMatch.2nd.F) %>%
            distinct() %>%
            dplyr::filter(!is.na(PatientID)) %>%
            group_by(PatientID, Locus, OrganMatch.2nd.F, allele_1, allele_2) %>%
            mutate(across(
               !c(PatientID, Locus, OrganMatch.2nd.F, allele_1, allele_2),
               ~ ifelse(is.na(.), lag(.), .)
            ))

         shared_eplets_bead_donor_recipient <- full_join(shared_eplets, combined.geno.eplets)

         return(shared_eplets_bead_donor_recipient)
      })
      # ------------------------------END EPLETS DONOR v/s RECIPIENT -----------------------------------

      # Eplets donor vs recipient (positions removed)
      eplets_match_geno_minimal <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         req(eplets_match_geno())

         suppressMessages(
            eplets_match_geno_no_pos <- eplets_match_geno() %>%
               mutate(
                  `shared_eplets(position)` = remove_brackets(`shared_eplets(position)`),
                  Bead_Specific_Eplet = remove_brackets(Bead_Specific_Eplet),
                  Eplet_of_Candidate_antigen = remove_brackets(Eplet_of_Candidate_antigen),
                  # allele 1
                  `shared_eplets_patient_donor_1(position)` = remove_brackets(`shared_eplets_patient_donor_1(position)`),
                  Recipient_Specific_Eplet_1 = remove_brackets(Recipient_Specific_Eplet_1),
                  Eplet_of_Candidate_antigen_1 = remove_brackets(Eplet_of_Candidate_antigen_1),
                  # allele 2
                  `shared_eplets_patient_donor_2(position)` = remove_brackets(`shared_eplets_patient_donor_2(position)`),
                  Recipient_Specific_Eplet_2 = remove_brackets(Recipient_Specific_Eplet_2),
                  Eplet_of_Candidate_antigen_2 = remove_brackets(Eplet_of_Candidate_antigen_2)
               ) %>%
               dplyr::select(
                  PatientID,
                  MFI,
                  Donor_HLA_exception,
                  Candidate_antigen,
                  MRO_Complex,
                  MRO_Serotype,
                  `shared_eplets(position)`,
                  Bead_Specific_Eplet,
                  Eplet_of_Candidate_antigen,
                  allele_1,
                  `shared_eplets_patient_donor_1(position)`,
                  Recipient_Specific_Eplet_1,
                  Eplet_of_Candidate_antigen_1,
                  allele_2,
                  `shared_eplets_patient_donor_2(position)`,
                  Recipient_Specific_Eplet_2,
                  Eplet_of_Candidate_antigen_2,
                  CWD_1.0,
                  CWD_2.0
               ) %>%
               dplyr::rename(
                  shared_eplets = `shared_eplets(position)`,
                  shared_eplets_patient_donor_1 = `shared_eplets_patient_donor_1(position)`,
                  shared_eplets_patient_donor_2 = `shared_eplets_patient_donor_2(position)`
               )
         )
         eplets_match_geno_no_pos$Predicted_reactive_eplet <- apply(eplets_match_geno_no_pos, 1, function(row) {
            # Split the values in the columns by "; "
            shared_eplets <- unlist(strsplit(as.character(row["shared_eplets"]), "; "))
            shared_eplets_patient_donor1 <- unlist(strsplit(as.character(row["shared_eplets_patient_donor_1"]), "; "))
            Recipient_Specific_Eplet1 <- unlist(strsplit(as.character(row["Recipient_Specific_Eplet_1"]), "; "))
            shared_eplets_patient_donor2 <- unlist(strsplit(as.character(row["shared_eplets_patient_donor_2"]), "; "))
            Recipient_Specific_Eplet2 <- unlist(strsplit(as.character(row["Recipient_Specific_Eplet_2"]), "; "))

            # Find unique elements in shared_eplets that are not in other columns
            unique_elements <- sort(
               setdiff(
                  shared_eplets,
                  unique(c(
                     shared_eplets_patient_donor1,
                     Recipient_Specific_Eplet1,
                     shared_eplets_patient_donor2,
                     Recipient_Specific_Eplet2
                  ))
               )
            )

            # Combine the unique elements into a single string separated by "; "
            unique_elements_string <- paste(unique_elements, collapse = "; ")

            return(unique_elements_string)
         })

         eplets_match_geno_no_pos <- eplets_match_geno_no_pos %>%
            mutate_all(as.character) %>%
            dplyr::select(
               PatientID,
               MFI,
               Donor_HLA_exception,
               Candidate_antigen,
               Predicted_reactive_eplet,
               shared_eplets,
               Bead_Specific_Eplet,
               Eplet_of_Candidate_antigen,
               allele_1,
               shared_eplets_patient_donor_1,
               Recipient_Specific_Eplet_1,
               Eplet_of_Candidate_antigen_1,
               allele_2,
               shared_eplets_patient_donor_2,
               Recipient_Specific_Eplet_2,
               Eplet_of_Candidate_antigen_2,
               CWD_1.0,
               CWD_2.0,
               MRO_Complex,
               MRO_Serotype
            ) %>%
            dplyr::filter(!grepl("^NOT_CHECKED", Predicted_reactive_eplet)) %>%
            mutate(
               Candidate_antigen = as.character(
                  lapply(strsplit(Candidate_antigen, ":"), function(x) {
                     if (length(x) >= 2) paste(x[[1]], x[[2]], sep = ":") else "NA"
                  })
               )
            ) %>%
            distinct(PatientID, Candidate_antigen, .keep_all = TRUE) %>%
            arrange(PatientID, Candidate_antigen)
         return(eplets_match_geno_no_pos)
      })
      # ------------------------------- END Eplet donor recipient no position ------------------------------

      # Eplets donor vs recipient (positions removed)
      eplets_match_geno_dedup <- reactive({
         req(eplets_match_geno_minimal())
         eplets_dedup <- eplets_match_geno_minimal() %>%
            dplyr::filter(!grepl("^NOT_CHECKED", Predicted_reactive_eplet)) %>%
            mutate(Candidate_antigen = as.character(
               lapply(strsplit(Candidate_antigen, ":"), function(x) {
                  ifelse(length(x) >= 2, paste(x[[1]], x[[2]], sep = ":"), "NA")
               })
            )) %>%
            distinct(PatientID, Candidate_antigen, .keep_all = TRUE) %>%
            dplyr::select(
               PatientID,
               MFI,
               Donor_HLA_exception,
               Candidate_antigen,
               Predicted_reactive_eplet,
               shared_eplets,
               Bead_Specific_Eplet,
               Eplet_of_Candidate_antigen,
               allele_1,
               shared_eplets_patient_donor_1,
               Recipient_Specific_Eplet_1,
               Eplet_of_Candidate_antigen_1,
               allele_2,
               shared_eplets_patient_donor_2,
               Recipient_Specific_Eplet_2,
               Eplet_of_Candidate_antigen_2
            ) %>%
            arrange(PatientID, Candidate_antigen)
         return(eplets_dedup)
      })
      # ------------------------- END Eplet donor recipient no position de-duplicated ------------------------------

      eplets_match_shortlisted <- reactive({
         req(input$genotype_upload)
         req(input$mfi_upload)

         antibody_analysis.df <- antibody_analysis()

         dedup_eplets <- eplets_match_geno_dedup()
         mfi_df <- p_serology_clean()

         known <- as.vector(unique(check.list.known$Known))
         mfi <- as.vector(unique(mfi_df$Allele))
         not_in_mfi <- data.frame(known[!known %in% mfi]) %>%
            set_names("Allele") %>%
            as_tibble()

         ## create data.frame with overlaps from
         negative_eplets <- dplyr::left_join(not_in_mfi, antibody_analysis.df) %>%
            distinct()
         unique_eplets_neg <- sort(unique(negative_eplets$eplet))

         priority_eplets <- dedup_eplets %>%
            dplyr::select(1:6) %>%
            mutate(
               Shortlisted_reactive_eplets = sapply(Predicted_reactive_eplet, function(eplets) {
                  eplets_list <- str_split(eplets, ";\\s*")[[1]]
                  shortlisted <- setdiff(eplets_list, unique_eplets_neg)
                  paste(shortlisted, collapse = "; ")
               })
            ) %>%
            dplyr::select(1, 2, 4, 3, 7, everything())
         return(priority_eplets)
      })
      # -------------------------  END Eplet donor recipient shortlisted ------------------------------

      eplets_donor_of_interest <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         req(input$donor_of_interest_genotype)
         antibody_analysis.df <- antibody_analysis()

         ## Create eplet list fordonor of interest
         donor_of_interest_geno.long <- donor_geno() %>%
            pivot_longer(cols = 2:ncol(.), names_to = "Locus_allele", values_to = "Typing") %>%
            mutate(
               allele = as.character(
                  lapply(
                     Locus_allele,
                     function(x) {
                        if (grepl("_1$", x, ignore.case = TRUE)) {
                           return("allele_1")
                        }
                        if (grepl("_2$", x, ignore.case = TRUE)) {
                           return("allele_2")
                        }
                     }
                  )
               )
            ) %>%
            mutate(
               Locus = gsub("_1$|_2$", "", Locus_allele),
               Typing = ifelse(!is.na(Typing), paste0(Locus, "*", Typing), Typing),
               Allele.2f = as.character(
                  lapply(
                     strsplit(Typing, ":"), function(x) {
                        ifelse(
                           length(x) >= 2,
                           paste(x[[1]], x[[2]], sep = ":"),
                           "NA"
                        )
                     }
                  )
               )
            )

         organ.match.melt.doi <- donor_of_interest_geno.long %>%
            ungroup() %>%
            dplyr::select(DonorID, Locus, Typing) %>%
            dplyr::rename(Candidate_antigen = Typing) %>%
            mutate(
               Allele.1.field = as.character(
                  lapply(
                     strsplit(Candidate_antigen, ":"),
                     function(x) {
                        ifelse(length(x) >= 1, x[[1]], "NA")
                     }
                  )
               )
            ) %>%
            na.omit()

         organ_match_donor <- plyr::join(p_serology_clean(), organ.match.melt.doi) %>%
            distinct() %>%
            dplyr::filter(
               Threshold == "PASS",
               Allele.1.field != Candidate_antigen
            ) %>%
            mutate(
               OrganMatch.2nd.F = ifelse(
                  !endsWith(Candidate_antigen, "G"),
                  as.character(lapply(
                     strsplit(Candidate_antigen, ":"),
                     function(x) {
                        ifelse(length(x) >= 2,
                           paste(x[[1]], x[[2]], sep = ":"), x
                        )
                     }
                  )),
                  Candidate_antigen
               )
            ) %>%
            dplyr::filter(stringr::str_trim(Allele) != stringr::str_trim(Candidate_antigen)) %>%
            dplyr::left_join(., (exceptions_df %>% dplyr::rename(Allele = Recipient_Antibody))) %>%
            dplyr::left_join(
               .,
               (cwd.alleles %>%
                  dplyr::select(1, 3, 4, 6) %>%
                  set_names(c("Locus", "Candidate_antigen", "CWD_1.0", "CWD_2.0"))
               )
            ) %>%
            dplyr::filter(Candidate_antigen != Allele.1.field) %>%
            as_tibble()

         class1 <- c("A", "B", "C")
         class2 <- c("DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")

         class1.alleles.d <- organ_match_donor %>%
            filter(
               grepl(paste0("^(", paste(class1, collapse = "|"), ")"), Candidate_antigen)
            ) %>%
            filter(
               !endsWith(Candidate_antigen, "G"),
               !grepl(":[A-Z]", Candidate_antigen)
            ) %>%
            distinct(Candidate_antigen) %>%
            pull(Candidate_antigen)

         class2.alleles.d <- organ_match_donor %>%
            filter(
               grepl(paste0("^(", paste(class2, collapse = "|"), ")"), Candidate_antigen)
            ) %>%
            filter(
               !endsWith(Candidate_antigen, "G"),
               !grepl(":[A-Z]", Candidate_antigen)
            ) %>%
            distinct(Candidate_antigen) %>%
            pull(Candidate_antigen)

         a <- get_serotypes(class1.alleles.d, mhc_type = "MHC-I") %>%
            enframe() %>%
            dplyr::rename(MRO_Complex = name, MRO_Serotype = value) %>%
            na.omit(MRO_Serotype) %>%
            mutate(
               OrganMatch.2nd.F = gsub(" protein complex| serotype|HLA-", "", MRO_Complex),
               MRO_Serotype = gsub(" protein complex| serotype|HLA-", "", MRO_Serotype)
            )

         b <- get_serotypes(class2.alleles.d, mhc_type = "MHC-II") %>%
            enframe() %>%
            dplyr::rename(MRO_Complex = name, MRO_Serotype = value) %>%
            na.omit(MRO_Serotype) %>%
            mutate(OrganMatch.2nd.F = gsub(" protein complex| serotype|HLA-", "", MRO_Complex)) %>%
            separate_rows(OrganMatch.2nd.F, sep = "/") %>%
            mutate(
               MRO_Serotype = gsub(" protein complex| serotype|HLA-", "", MRO_Serotype)
            )

         Serotype.d <- bind_rows(a, b) %>%
            na.omit(MRO_Serotype) %>%
            distinct() %>%
            as.data.frame() %>%
            InterMineR::simplifyResult(.,
               index_column = "OrganMatch.2nd.F",
               values_column = "MRO_Serotype",
               returnInDataframe = T
            ) %>%
            mutate(simplified_results = sapply(
               strsplit(as.character(simplified_results), ","), function(x) {
                  paste(sort(unique(x)), collapse = ",")
               }
            )) %>%
            dplyr::select(-MRO_Serotype) %>%
            dplyr::rename(MRO_Serotype = simplified_results) %>%
            arrange(OrganMatch.2nd.F, MRO_Serotype) %>%
            group_by(OrganMatch.2nd.F) %>%
            slice(1) %>%
            mutate(MRO_Serotype = as.character(MRO_Serotype))

         organ_match_donor <- organ_match_donor %>%
            left_join(., Serotype.d, relationship = "many-to-many") %>%
            distinct() %>%
            mutate_all(as.character)

         collapse_eplet_list <- function(x) paste(x, collapse = "; ")

         cl <- makeCluster(detectCores())
         registerDoParallel(cl)

         eplets_shared_doi <- data.frame()

         df_doi <- foreach(i = 1:nrow(organ_match_donor), .combine = rbind) %dopar% {
            x <- organ_match_donor[i, ]
            a <- stringr::str_trim(x[["Allele"]], side = c("both"))
            b <- stringr::str_trim(x[["OrganMatch.2nd.F"]], side = c("both"))

            if (endsWith(b, "G") || a == b) {
               x["shared_eplets(position)"] <- "NOT_CHECKED (identical donor-recipient or G-Group alleles)"
               x["Bead_Specific_Eplet"] <- "NOT_CHECKED (identical donor-recipient or G-Group alleles)"
               x["Eplet_of_Candidate_antigen"] <- "NOT_CHECKED (identical donor-recipient or G-Group alleles)"
               eplets_shared_doi <- rbind(eplets_shared_doi, x)
            }

            if (a != b && !endsWith(b, "G")) {
               df1 <- dplyr::filter(antibody_analysis.df, Allele == a)
               df1 <- subset(df1, select = c(2:3))
               df2 <- dplyr::filter(antibody_analysis.df, Allele == b)
               df2 <- subset(df2, select = c(2:3))
               res <- dplyr::full_join(df1, df2)
               df1.1 <- dplyr::filter(antibody_analysis.df, Allele == a)
               df2.1 <- dplyr::filter(antibody_analysis.df, Allele == b)

               res.uniq.p <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
               res.uniq.p <- dplyr::filter(res.uniq.p, is.na(Allele.y))
               res.uniq.o <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
               res.uniq.o <- dplyr::filter(res.uniq.o, is.na(Allele.x))

               if (nrow(res) >= 2 | nrow(res.uniq.p) >= 2 | nrow(res.uniq.o) >= 2) {
                  res$PatientID <- x[[1]]
                  res.uniq.p$PatientID <- x[[1]]
                  res.uniq.o$PatientID <- x[[1]]

                  a <- ifelse(nrow(res) >= 1,
                     collapse_eplet_list(
                        paste0(res$eplet, "(", res$exon.position, ")")
                     ),
                     "Eplet data missing"
                  )
                  b <- ifelse(nrow(res.uniq.p) >= 1,
                     collapse_eplet_list(
                        paste0(res.uniq.p$eplet, "(", res.uniq.p$exon.position, ")")
                     ),
                     "Eplet data missing"
                  )
                  c <- ifelse(nrow(res.uniq.o) >= 1,
                     collapse_eplet_list(
                        paste0(res.uniq.o$eplet, "(", res.uniq.o$exon.position, ")")
                     ),
                     "Eplet data missing"
                  )

                  a_elements <- strsplit(a, "; ")[[1]]
                  b_elements <- strsplit(b, "; ")[[1]]
                  c_elements <- strsplit(c, "; ")[[1]]

                  filtered_elements <- a_elements[!(a_elements %in% b_elements | a_elements %in% c_elements)]

                  filtered_a <- paste(filtered_elements, collapse = "; ")

                  x["shared_eplets(position)"] <- filtered_a
                  x["Bead_Specific_Eplet"] <- b
                  x["Eplet_of_Candidate_antigen"] <- c
                  eplets_shared_doi <- rbind(eplets_shared_doi, x)
               }
            }
         }

         stopCluster(cl)

         shared_eplets_doi <- distinct(df_doi) %>%
            full_join(
               .,
               (exceptions_df %>%
                  dplyr::rename(Allele = Recipient_Antibody)
               )
            ) %>%
            dplyr::rename(
               "DonorOfInterest.2nd.F" = `OrganMatch.2nd.F`,
               Donor_Specific_Eplet = Bead_Specific_Eplet
            ) %>%
            filter(
               !is.na(MFI),
               !is.na(PatientID)
            )

         return(shared_eplets_doi)
      })

      # Eplets donor-of-interest vs recipient -------------------------------

      eplets_match_geno_doi <- reactive({
         req(input$mfi_upload)
         req(input$genotype_upload)
         req(input$donor_of_interest_genotype)
         antibody_analysis.df <- antibody_analysis()

         suppressMessages(
            patient_geno.long <- patient_geno() %>%
               pivot_longer(
                  cols = 2:ncol(.),
                  names_to = "Locus_allele",
                  values_to = "Typing"
               ) %>%
               mutate(
                  allele = as.character(
                     lapply(Locus_allele, function(x) {
                        if (grepl("_1$", x, ignore.case = TRUE)) {
                           return("allele_1")
                        }
                        if (grepl("_2$", x, ignore.case = TRUE)) {
                           return("allele_2")
                        }
                     })
                  )
               ) %>%
               mutate(
                  Locus = gsub("_1$|_2$", "", Locus_allele),
                  Typing = ifelse(!is.na(Typing), paste0(Locus, "*", Typing), Typing),
                  Patient_genotype.2f = as.character(
                     lapply(strsplit(Typing, ":"), function(x) {
                        ifelse(length(x) >= 2, paste(x[[1]], x[[2]], sep = ":"), "NA")
                     })
                  )
               )
         )

         shared_eplets_doi <- eplets_donor_of_interest()

         suppressMessages(
            combine_doi_patient_typing <- dplyr::left_join(
               shared_eplets_doi %>%
                  dplyr::select(PatientID, Locus, DonorOfInterest.2nd.F) %>%
                  distinct(),
               patient_geno.long %>%
                  mutate(PatientID = stringr::str_trim(PatientID, side = c("both"))) %>%
                  dplyr::select(PatientID, allele, Locus, Patient_genotype.2f),
               relationship = "many-to-many"
            ) %>%
               dplyr::filter(!is.na(PatientID)) %>%
               pivot_wider(names_from = allele, values_from = Patient_genotype.2f)
         )


         fx <- function(x) paste(x, collapse = ", ")
         remove_brackets <- function(x) {
            as.character(lapply(x, function(y) {
               as.character(gsub("\\s*\\([^\\)]+\\)", "", y))
            }))
         }

         collapse_eplet_list <- function(x) paste(x, collapse = "; ")

         antibody_analysis.df <- antibody_analysis()

         eplets_shared_geno <- data.frame()

         # Start parallel backend
         cl <- makeCluster(detectCores())
         registerDoParallel(cl)

         doi.df.geno_allele1 <- foreach(i = 1:nrow(combine_doi_patient_typing), .combine = rbind) %dopar% {
            ## Get patient and organ match allele from row
            row <- combine_doi_patient_typing[i, ]
            patient.allele1 <- stringr::str_trim(row[["allele_1"]], side = c("both"))
            doi_allele <- stringr::str_trim(row[["DonorOfInterest.2nd.F"]], side = c("both"))

            ## Add 3 new columns to row
            row[["shared_eplets_doi_patient_donor_1(position)"]] <- NA
            row[["Recipient_Specific_Eplet_1"]] <- NA
            row[["Eplet_of_Candidate_antigen_1"]] <- NA

            ## Perform if-else conditions
            suppressMessages(
               if (is.na(patient.allele1) | is.na(doi_allele)) {
                  row[["shared_eplets_doi_patient_donor_1(position)"]] <- "No typing result"
                  row[["Recipient_Specific_Eplet_1"]] <- "No typing result"
                  row[["Eplet_of_Candidate_antigen_1"]] <- "No typing result"
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               } else if (endsWith(doi_allele, "G")) {
                  row[["shared_eplets_doi_patient_donor_1(position)"]] <- "NOT_CHECKED_G-Group donor alleles"
                  row[["Recipient_Specific_Eplet_1"]] <- "NOT_CHECKED_G-Group donor alleles"
                  row[["Eplet_of_Candidate_antigen_1"]] <- "NOT_CHECKED_G-Group donor alleles"
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               } else if (patient.allele1 == doi_allele) {
                  row[["shared_eplets_doi_patient_donor_1(position)"]] <- "NOT_CHECKED_identical donor-recipient alleles"
                  row[["Recipient_Specific_Eplet_1"]] <- "NOT_CHECKED_identical donor-recipient alleles"
                  row[["Eplet_of_Candidate_antigen_1"]] <- "NOT_CHECKED_identical donor-recipient alleles"
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               } else if (!is.na(patient.allele1) | !is.na(doi_allele)) {
                  df1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele1)
                  df1 <- subset(df1, select = c(2:3))
                  df2 <- dplyr::filter(antibody_analysis.df, Allele == doi_allele)
                  df2 <- subset(df2, select = c(2:3))
                  res <- dplyr::full_join(df1, df2)
                  df1.1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele1)
                  df2.1 <- dplyr::filter(antibody_analysis.df, Allele == doi_allele)

                  res.uniq.p <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
                  res.uniq.p <- dplyr::filter(res.uniq.p, is.na(Allele.y))
                  res.uniq.o <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
                  res.uniq.o <- dplyr::filter(res.uniq.o, is.na(Allele.x))

                  if (nrow(res) >= 2 || nrow(res.uniq.p) >= 2 || nrow(res.uniq.o) >= 2) {
                     res$PatientID <- row[["PatientID"]]
                     res.uniq.p$PatientID <- row[["PatientID"]]
                     res.uniq.o$PatientID <- row[["PatientID"]]

                     a <- ifelse(nrow(res) >= 1,
                        collapse_eplet_list(paste0(res$eplet, "(", res$exon.position, ")")),
                        "Eplet data missing"
                     )
                     b <- ifelse(nrow(res.uniq.p) >= 1,
                        collapse_eplet_list(paste0(res.uniq.p$eplet, "(", res.uniq.p$exon.position, ")")),
                        "Eplet data missing"
                     )
                     c <- ifelse(nrow(res.uniq.o) >= 1,
                        collapse_eplet_list(paste0(res.uniq.o$eplet, "(", res.uniq.o$exon.position, ")")),
                        "Eplet data missing"
                     )

                     # Split the strings in a, b, and c into individual elements
                     a_elements <- strsplit(a, "; ")[[1]]
                     b_elements <- strsplit(b, "; ")[[1]]
                     c_elements <- strsplit(c, "; ")[[1]]

                     # Remove elements of b and c from a
                     filtered_elements <- a_elements[!(a_elements %in% b_elements | a_elements %in% c_elements)]

                     # Join the elements back into a string
                     filtered_a <- paste(filtered_elements, collapse = "; ")

                     row[["shared_eplets_doi_patient_donor_1(position)"]] <- filtered_a
                     row[["Recipient_Specific_Eplet_1"]] <- b
                     row[["Eplet_of_Candidate_antigen_1"]] <- c
                     eplets_shared_geno <- rbind(eplets_shared_geno, row)
                  }
               }
            )
         }

         eplets_shared_geno <- data.frame()

         doi.df.geno_allele2 <- foreach(i = 1:nrow(combine_doi_patient_typing), .combine = rbind) %dopar% {
            ## Get patient and organ match allele from row
            row <- combine_doi_patient_typing[i, ]
            patient.allele2 <- stringr::str_trim(row[["allele_2"]], side = c("both"))
            doi_allele <- stringr::str_trim(row[["DonorOfInterest.2nd.F"]], side = c("both"))

            ## Add 3 new columns to row
            row[["shared_eplets_doi_patient_donor_2(position)"]] <- NA
            row[["Recipient_Specific_Eplet_2"]] <- NA
            row[["Eplet_of_Candidate_antigen_2"]] <- NA

            ## Perform if-else conditions
            if (is.na(patient.allele2) | is.na(doi_allele)) {
               row[["shared_eplets_doi_patient_donor_2(position)"]] <- "No typing result"
               row[["Recipient_Specific_Eplet_2"]] <- "No typing result"
               row[["Eplet_of_Candidate_antigen_2"]] <- "No typing result"
               eplets_shared_geno <- rbind(eplets_shared_geno, row)
            } else if (endsWith(doi_allele, "G")) {
               row[["shared_eplets_doi_patient_donor_2(position)"]] <- "NOT_CHECKED_G-Group donor alleles"
               row[["Recipient_Specific_Eplet_2"]] <- "NOT_CHECKED_G-Group donor alleles"
               row[["Eplet_of_Candidate_antigen_2"]] <- "NOT_CHECKED_G-Group donor alleles"
               eplets_shared_geno <- rbind(eplets_shared_geno, row)
            } else if (patient.allele2 == doi_allele) {
               row[["shared_eplets_doi_patient_donor_2(position)"]] <- "NOT_CHECKED_identical donor-recipient alleles"
               row[["Recipient_Specific_Eplet_2"]] <- "NOT_CHECKED_identical donor-recipient alleles"
               row[["Eplet_of_Candidate_antigen_2"]] <- "NOT_CHECKED_identical donor-recipient alleles"
               eplets_shared_geno <- rbind(eplets_shared_geno, row)
            } else if (!is.na(patient.allele2) | !is.na(doi_allele)) {
               df1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele2)
               df1 <- subset(df1, select = c(2:3))
               df2 <- dplyr::filter(antibody_analysis.df, Allele == doi_allele)
               df2 <- subset(df2, select = c(2:3))
               res <- dplyr::full_join(df1, df2)
               df1.1 <- dplyr::filter(antibody_analysis.df, Allele == patient.allele2)
               df2.1 <- dplyr::filter(antibody_analysis.df, Allele == doi_allele)

               res.uniq.p <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
               res.uniq.p <- dplyr::filter(res.uniq.p, is.na(Allele.y))
               res.uniq.o <- dplyr::full_join(df1.1, df2.1, by = c("exon.position", "eplet"))
               res.uniq.o <- dplyr::filter(res.uniq.o, is.na(Allele.x))

               if (nrow(res) >= 2 || nrow(res.uniq.p) >= 2 || nrow(res.uniq.o) >= 2) {
                  res$PatientID <- row[["PatientID"]]
                  res.uniq.p$PatientID <- row[["PatientID"]]
                  res.uniq.o$PatientID <- row[["PatientID"]]

                  a <- ifelse(nrow(res) >= 1,
                     collapse_eplet_list(paste0(res$eplet, "(", res$exon.position, ")")),
                     "Eplet data missing"
                  )
                  b <- ifelse(nrow(res.uniq.p) >= 1,
                     collapse_eplet_list(paste0(res.uniq.p$eplet, "(", res.uniq.p$exon.position, ")")),
                     "Eplet data missing"
                  )
                  c <- ifelse(nrow(res.uniq.o) >= 1,
                     collapse_eplet_list(paste0(res.uniq.o$eplet, "(", res.uniq.o$exon.position, ")")),
                     "Eplet data missing"
                  )

                  # Split the strings in a, b, and c into individual elements
                  a_elements <- strsplit(a, "; ")[[1]]
                  b_elements <- strsplit(b, "; ")[[1]]
                  c_elements <- strsplit(c, "; ")[[1]]

                  # Remove elements of b and c from a
                  filtered_elements <- a_elements[!(a_elements %in% b_elements | a_elements %in% c_elements)]

                  # Join the elements back into a string
                  filtered_a <- paste(filtered_elements, collapse = "; ")

                  row[["shared_eplets_doi_patient_donor_2(position)"]] <- filtered_a
                  row[["Recipient_Specific_Eplet_2"]] <- b
                  row[["Eplet_of_Candidate_antigen_2"]] <- c
                  eplets_shared_geno <- rbind(eplets_shared_geno, row)
               }
            }
         }

         # Stop parallel back-end
         stopCluster(cl)

         doi.df.geno_allele1 <- doi.df.geno_allele1 %>%
            distinct() %>%
            arrange(PatientID, Locus, DonorOfInterest.2nd.F)

         doi.df.geno_allele2 <- doi.df.geno_allele2 %>%
            distinct() %>%
            arrange(PatientID, Locus, DonorOfInterest.2nd.F)

         doi.combined.geno.eplets <- full_join(doi.df.geno_allele1, doi.df.geno_allele2, relationship = "many-to-many") %>%
            arrange(PatientID, Locus, DonorOfInterest.2nd.F) %>%
            distinct() %>%
            dplyr::filter(!is.na(PatientID)) %>%
            group_by(PatientID, Locus, DonorOfInterest.2nd.F, allele_1, allele_2) %>%
            mutate(across(
               !c(PatientID, Locus, DonorOfInterest.2nd.F, allele_1, allele_2),
               ~ ifelse(is.na(.), lag(.), .)
            ))

         suppressMessages(
            shared_eplets_doi_recipient <- full_join(shared_eplets_doi, doi.combined.geno.eplets) %>%
               # dplyr::rename(Donor_Specific_Eplet = Bead_Specific_Eplet) %>%
               filter(
                  !is.na(MFI),
                  !is.na(PatientID)
               ) %>%
               mutate_all(as.character)
         )

         return(shared_eplets_doi_recipient)
      })


      ## format eplets_match_geno_doi() further

      eplets_match_geno_doi_clean <- reactive({
         req(input$genotype_upload)
         req(input$mfi_upload)

         eplets_match_geno_doi.df <- eplets_match_geno_doi()
         if (nrow(eplets_match_geno_doi.df > 0)) {
            eplets_match_geno_doi.df <- eplets_match_geno_doi.df %>%
               filter(
                  !is.na(MFI),
                  !is.na(PatientID)
               ) %>%
               # dplyr::rename(Donor_Specific_Eplet = Bead_Specific_Eplet) %>%
               mutate(
                  `shared_eplets(position)` = remove_brackets(`shared_eplets(position)`),
                  Donor_Specific_Eplet = remove_brackets(Donor_Specific_Eplet),
                  Eplet_of_Candidate_antigen = remove_brackets(Candidate_antigen),
                  Donor_HLA_exception = add_space_after_comma(Donor_HLA_exception),
                  `shared_eplets_doi_patient_donor_1(position)` = remove_brackets(`shared_eplets_doi_patient_donor_1(position)`),
                  Recipient_Specific_Eplet_1 = remove_brackets(`Recipient_Specific_Eplet_1`)
               ) %>%
               filter(
                  !Candidate_antigen %in% check.list.known$Known,
                  !grepl("^NOT_CHECKED", `shared_eplets(position)`)
               ) %>%
               as_tibble() %>%
               filter(!is.na(PatientID))
            eplets_match_geno_doi.df$Predicted_reactive_eplet <- apply(eplets_match_geno_doi.df, 1, function(row) {
               # Split the values in the columns by "; "
               shared_eplets <- unlist(strsplit(as.character(row["shared_eplets(position)"]), "; "))
               shared_eplets_patient_doi <- unlist(
                  strsplit(
                     as.character(row["shared_eplets_doi_patient_donor_1(position)"]), "; "
                  )
               )
               Recipient_Specific_Eplet_1 <- unlist(
                  strsplit(
                     as.character(row["Recipient_Specific_Eplet_1"]), "; "
                  )
               )

               # Find unique elements in shared_eplets that are not in other columns
               unique_elements <- sort(setdiff(
                  shared_eplets,
                  unique(
                     shared_eplets_patient_doi,
                     Recipient_Specific_Eplet_1
                  )
               ))
               # Combine the unique elements into a single string separated by "; "
               unique_elements_string <- paste(unique_elements, collapse = "; ")

               return(unique_elements_string)
            })

            eplets_match_geno_doi.df <- eplets_match_geno_doi.df %>%
               dplyr::rename(shared_eplets = `shared_eplets(position)`) %>%
               filter(!is.na(PatientID))
         }

         antibody_analysis.df <- antibody_analysis()
         mfi_df <- p_serology_clean()
         known <- as.vector(unique(check.list.known$Known))
         mfi <- as.vector(unique(mfi_df$Allele))
         not_in_mfi <- data.frame(known[!known %in% mfi]) %>%
            set_names("Allele") %>%
            as_tibble()

         ## create data.frame with overlaps from
         negative_eplets <- dplyr::left_join(not_in_mfi, antibody_analysis.df) %>%
            distinct()
         unique_eplets_neg <- sort(unique(negative_eplets$eplet))

         priority_eplets <- eplets_match_geno_doi.df %>%
            # dplyr::select(1:6) %>%
            mutate(
               Shortlisted_reactive_eplets = sapply(Predicted_reactive_eplet, function(eplets) {
                  eplets_list <- str_split(eplets, ";\\s*")[[1]]
                  shortlisted <- setdiff(eplets_list, unique_eplets_neg)
                  paste(shortlisted, collapse = "; ")
               })
            ) %>%
            dplyr::select(
               PatientID,
               MFI,
               DonorID,
               DonorOfInterest.2nd.F,
               MRO_Serotype,
               Shortlisted_reactive_eplets,
               Predicted_reactive_eplet,
               shared_eplets,
               Candidate_antigen,
               Eplet_of_Candidate_antigen,
               # Donor_Specific_Eplet,
               everything()
            )

         ## TODO
         ## IF priotity list is empty, return previous list
         return(priority_eplets)
         print(priority_eplets)
      })

      # -------------------------  END Eplet analysis for "donor of interest" ------------------------------

      ## ++++++++++++++++++ Render output tables in DT::datatable() for display in TABS++++++++++++++++++++++++++++++

      output$table_display <- renderDataTable({
         f <- mfi_data()

         suppressMessages(
            f <- f %>%
               dplyr::select(c("PatientID", "Class", "Final Assignment ALLELES with MFIs")) %>%
               dplyr::rename("MFI" = "Final Assignment ALLELES with MFIs") %>%
               na.omit() %>%
               arrange(PatientID)
         )

         DT::datatable(
            f,
            options = list(
               pageLength = 20,
               scrollY = 500,
               orderClasses = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 1)
            ),
            rownames = F
         )
      })

      output$geno_display <- renderDataTable({
         f0 <- patient_geno() %>%
            as_tibble()

         DT::datatable(
            f0,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               pageLength = 20,
               scrollY = 500,
               orderClasses = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 1)
            ),
            rownames = FALSE
         )
      })

      output$cleaned_table <- renderDataTable({
         f1 <- p_serology_clean()

         DT::datatable(
            f1,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               pageLength = 20,
               scrollY = 500,
               orderClasses = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 1)
            ),
            rownames = F
         )
      })

      output$organ_match_tb <- renderDataTable({
         f2 <- organ_match()

         DT::datatable(
            f2,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               pageLength = 20,
               scrollY = 500,
               orderClasses = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 1)
            ),
            rownames = F
         )
      })

      output$cwd_table <- renderDataTable({
         f3 <- cwd_match()

         DT::datatable(
            f3,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               deferRender = F,
               dom = "t",
               pageLength = 20,
               scrollY = 500,
               scroller = TRUE,
               scrollX = T,
               fixedColumns = list(leftColumns = 2)
            ),
            rownames = F
         )
      })

      output$eplets <- renderDataTable({
         f4 <- eplets_match()

         suppressMessages(
            f4 <- f4 %>%
               dplyr::filter(!is.na(MFI), !is.na(PatientID)) %>%
               dplyr::filter(!OrganMatch.2nd.F %in% check.list.known$Known) %>%
               dplyr::filter(!grepl("^NOT_CHECKED", `shared_eplets(position)`)) %>%
               mutate(
                  shared_eplets = remove_brackets(`shared_eplets(position)`),
                  Bead_Specific_Eplet = remove_brackets(Bead_Specific_Eplet),
                  Eplet_of_Candidate_antigen = remove_brackets(Eplet_of_Candidate_antigen)
               ) %>%
               dplyr::select(
                  PatientID,
                  MFI,
                  Donor_HLA_exception,
                  Candidate_antigen,
                  MRO_Complex,
                  MRO_Serotype,
                  `shared_eplets(position)`,
                  Bead_Specific_Eplet,
                  Eplet_of_Candidate_antigen,
                  Threshold,
                  CWD_1.0,
                  CWD_2.0,
                  OrganMatch.2nd.F,
                  Allele,
                  Locus
               )
         )

         DT::datatable(
            f4,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               deferRender = F,
               dom = "t",
               pageLength = 20,
               scrollY = 500,
               scroller = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 2)
            ),
            rownames = F
         )
      })

      output$eplets_geno <- renderDataTable({
         f5 <- eplets_match_geno_minimal()

         DT::datatable(
            f5,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               deferRender = F,
               dom = "t",
               pageLength = 20,
               scrollY = 500,
               scroller = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 2)
            ),
            rownames = F
         )
      })

      output$deduplicated <- renderDataTable({
         f6 <- eplets_match_geno_dedup()

         DT::datatable(
            f6,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               deferRender = F,
               dom = "t",
               pageLength = 20,
               scrollY = 500,
               scroller = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 2)
            ),
            rownames = F
         )
      })

      output$shortlisted_eplets <- renderDataTable({
         f6 <- eplets_match_shortlisted()
         f6 <- f6 %>% filter(!is.na(Shortlisted_reactive_eplets) | Shortlisted_reactive_eplets != "")

         DT::datatable(
            f6,
            extensions = c("FixedColumns", "Scroller"),
            options = list(
               deferRender = F,
               dom = "t",
               pageLength = 20,
               scrollY = 500,
               scroller = TRUE,
               scrollX = TRUE,
               fixedColumns = list(leftColumns = 3)
            ),
            rownames = F
         )
      })

      output$donor_of_interest <- renderDataTable({
         if (input$donor_interest_checkbox) {
            f7 <- eplets_match_geno_doi_clean()
            if (nrow(f7 > 0)) {
               rownames(f7) <- NULL
            }

            if (nrow(f7) == 0) {
               return(as.data.frame("No shortlisted reactive eplets for donor of interest!"))
            } else {
               DT::datatable(
                  f7,
                  extensions = c("FixedColumns", "Scroller"),
                  options = list(
                     deferRender = FALSE,
                     dom = "t",
                     pageLength = 20,
                     scrollY = 500,
                     scroller = TRUE,
                     scrollX = TRUE,
                     fixedColumns = list(leftColumns = 1)
                  ),
                  rownames = FALSE
               )
            }
         } else {
            return(as.data.frame("Donor of interest data not provided"))
         }
      })


      ## ++++++++++++++++++++++++++++ capture output tables into user download button +++++++++++++++++++++++++++++
      output$downloadData0 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_organmatch.tsv")
         },
         content = function(file) {
            df <- organ_match()
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData1 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_CWD_match.tsv")
         },
         content = function(file) {
            df <- cwd_match()
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData2 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_CWD_summary.tsv")
         },
         content = function(file) {
            df <- cwd_match_summary()
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData3 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_bead_donor.tsv")
         },
         content = function(file) {
            df <- eplets_match()
            df <- df %>%
               dplyr::filter(!Candidate_antigen %in% check.list.known$Known) %>%
               dplyr::filter(!grepl("^NOT_CHECKED", `shared_eplets(position)`)) %>%
               dplyr::select(
                  PatientID,
                  MFI,
                  Donor_HLA_exception,
                  Candidate_antigen,
                  MRO_Complex,
                  MRO_Serotype,
                  `shared_eplets(position)`,
                  Bead_Specific_Eplet,
                  Eplet_of_Candidate_antigen,
                  Threshold,
                  CWD_1.0,
                  CWD_2.0,
                  OrganMatch.2nd.F,
                  Allele,
                  Locus
               )
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData3xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_bead_donor.xlsx")
         },
         content = function(file) {
            df <- eplets_match()
            df <- df %>%
               dplyr::filter(!Candidate_antigen %in% check.list.known$Known) %>%
               dplyr::filter(!grepl("^NOT_CHECKED", `shared_eplets(position)`)) %>%
               dplyr::select(
                  PatientID,
                  MFI,
                  Donor_HLA_exception,
                  Candidate_antigen,
                  MRO_Complex,
                  MRO_Serotype,
                  `shared_eplets(position)`,
                  Bead_Specific_Eplet,
                  Eplet_of_Candidate_antigen,
                  Threshold,
                  CWD_1.0,
                  CWD_2.0,
                  OrganMatch.2nd.F,
                  Allele,
                  Locus
               )
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData4 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_bead_donor_no_position.tsv")
         },
         content = function(file) {
            df <- suppressMessages(
               eplets_match() %>%
                  dplyr::filter(!Candidate_antigen %in% check.list.known$Known) %>%
                  dplyr::filter(!grepl("^NOT_CHECKED", `shared_eplets(position)`)) %>%
                  dplyr::select(
                     PatientID,
                     MFI,
                     Donor_HLA_exception,
                     Candidate_antigen,
                     MRO_Complex,
                     MRO_Serotype,
                     `shared_eplets(position)`,
                     Bead_Specific_Eplet,
                     Eplet_of_Candidate_antigen,
                     Threshold,
                     CWD_1.0,
                     CWD_2.0,
                     OrganMatch.2nd.F,
                     Allele,
                     Locus
                  ) %>%
                  mutate(
                     shared_eplets = remove_brackets(`shared_eplets(position)`),
                     Bead_Specific_Eplet = remove_brackets(Bead_Specific_Eplet),
                     Eplet_of_Candidate_antigen = remove_brackets(Eplet_of_Candidate_antigen)
                  )
            )
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData4xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_bead_donor_no_position.xlsx")
         },
         content = function(file) {
            df <- suppressMessages(
               eplets_match() %>%
                  dplyr::filter(!Candidate_antigen %in% check.list.known$Known) %>%
                  dplyr::filter(!grepl("^NOT_CHECKED", `shared_eplets(position)`)) %>%
                  dplyr::select(
                     PatientID,
                     MFI,
                     Donor_HLA_exception,
                     Candidate_antigen,
                     MRO_Complex,
                     MRO_Serotype,
                     `shared_eplets(position)`,
                     Bead_Specific_Eplet,
                     Eplet_of_Candidate_antigen,
                     Threshold,
                     CWD_1.0,
                     CWD_2.0,
                     OrganMatch.2nd.F,
                     Allele,
                     Locus
                  ) %>%
                  mutate(
                     shared_eplets = remove_brackets(`shared_eplets(position)`),
                     Bead_Specific_Eplet = remove_brackets(Bead_Specific_Eplet),
                     Eplet_of_Candidate_antigen = remove_brackets(Eplet_of_Candidate_antigen)
                  )
            )
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData5 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_donor_recipient.tsv")
         },
         content = function(file) {
            df <- eplets_match_geno()
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData5xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_donor_recipient.xlsx")
         },
         content = function(file) {
            df <- eplets_match_geno()
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData6 <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_donor_recipient_no_position.tsv")
         },
         content = function(file) {
            df <- eplets_match_geno_minimal()
            write_tsv(df, file, col_names = TRUE, na = "NA", quote = "needed", num_threads = 4)
         }
      )

      output$downloadData6xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_donor_recipient_no_position_complete.xlsx")
         },
         content = function(file) {
            df <- eplets_match_geno_minimal()
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData7xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shared_eplets_donor_recipient_no_position.xlsx")
         },
         content = function(file) {
            df <- eplets_match_geno_dedup()
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData8xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shortlisted_reactive_eplets_summary.xlsx")
         },
         content = function(file) {
            df <- eplets_match_shortlisted() %>%
               filter(!is.na(Shortlisted_reactive_eplets) | Shortlisted_reactive_eplets != "")
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData9xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_shortlisted_reactive_eplets.xlsx")
         },
         content = function(file) {
            df <- eplets_match_geno_dedup()
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      output$downloadData10xls <- downloadHandler(
         filename = function() {
            outfile_name <- gsub("_MFI", "", tools::file_path_sans_ext(input$mfi_upload$name))
            paste0(outfile_name, "_donor_of_interest_complete.xlsx")
         },
         content = function(file) {
            df <- eplets_match_geno_doi_clean()
            openxlsx::write.xlsx(df, file, asTable = TRUE, overwrite = TRUE)
         }
      )

      ## ++++++++++++++++++++++ CALL THEREACTIVE FUNCTIONS ++++++++++++++++++++++++++++++++
      mfi_data()
      patient_geno()
      p_serology_clean()
      organ_match()
      cwd_match()
      cwd_match_summary()
      antibody_analysis()
      eplets_match()
      eplets_match_geno()
      eplets_match_geno_minimal()
      eplets_match_geno_dedup()
      eplets_match_shortlisted()
      donor_geno()
      eplets_donor_of_interest()
      eplets_match_geno_doi()
   })


   ## ================= DEFINE RESET FUNCTIONS FOR CLEAR BUTTON ============================
   observeEvent(input$resetButton, {
      # Clear output tables
      output$table_display <- NULL
      output$geno_display <- NULL
      output$cleaned_table <- NULL
      output$organ_match_tb <- NULL
      output$cwd_table <- NULL
      output$eplets <- NULL
      output$eplets_geno <- NULL
      output$shortlisted_eplets <- NULL
      output$donor_of_interest <- NULL

      # Reload the session and clear everything
      session$reload()

      # Clear download handlers
      output$downloadData0 <- NULL
      output$downloadData1 <- NULL
      output$downloadData2 <- NULL
      output$downloadData3 <- NULL
      output$downloadData4 <- NULL
      output$downloadData4xls <- NULL
      output$downloadData5 <- NULL
      output$downloadData5xls <- NULL
      output$downloadData6 <- NULL
      output$downloadData6xls <- NULL
      output$downloadData7xls <- NULL
      output$downloadData8xls <- NULL
      output$downloadData9xls <- NULL
      output$downloadData10xls <- NULL
   })

   # Hide progress div once data is processed
   shinyjs::hide("progress_div")

   observeEvent(input$help, {
      shinyalert::shinyalert(
         title = "Usage information",
         text = '<div style="color: black; text-align: justify;">
            This app takes the patient’s HLA antibody profile and HLA genotype files,
            and an optional donor-of-interest HLA genotype file for cross-comparison.<br>
            The inputs are as follows:
            <ul style="list-style-type: square; color: blue; text-align: justify;">
               <li>the antibody results and HLA genotype of the recipient from the
               <a href="https://pubmed.ncbi.nlm.nih.gov/36851856/" target="_blank">OLI SAB assay</a>
               and SOFT Laboratory Information System, respectively,</li>
               <li>population data of deceased donor HLA genotypes extracted from the OrganMatch portal, and</li>
               <li>HLA eplets from the HLA Matchmaker tool (<a href="http://www.epitopes.net/"
               target="_blank">Duquesnoy 2002</a>).</li>
            </ul><br>
            You can download sample input files here:
            <ul style="list-style-type: square; color: blue; text-align: justify;">
               <li><a href="Donor_of_interest.csv" download>Donor of Interest (CSV)</a></li>
               <li><a href="SamplePatient_HLA.CSV" download>Sample Patient HLA (CSV)</a></li>
               <li><a href="SamplePatient_MFI.csv" download>Sample Patient MFI (CSV)</a></li>
            </ul><br>
            The software performs the following analyses:
            <ul style="list-style-type: square; color: blue; text-align: justify;">
               <li>identifies alleles that are not included in OLI SAB but are present in our deceased donor population,</li>
               <li>determines whether these identified alleles have shared eplets or serotypes with the patient’s high-risk anti-HLA specificities (with an MFI exceeding 4000 or custom input threshold in strength),</li>
               <li>excludes any self-eplets that are shared between the recipient and candidate donor,</li>
               <li>excludes any eplets on the negative/unassigned beads, and</li>
               <li>consolidates the final list of potential eplets and candidate high-risk donor antigens.</li>
            </ul>
            </div>',
         type = "info",
         closeOnEsc = TRUE,
         closeOnClickOutside = TRUE,
         showConfirmButton = TRUE,
         html = TRUE,
         size = "m"
      )
   })

   observeEvent(input$terms_of_use_link, {
      shinyalert::shinyalert(
         title = "Terms of Use",
         text = '<div style="color: black; text-align: justify;">
            The copyright owner(s) and developer(s) assume no responsibility for any injury to person or
            damage to persons or property, real or perceived, arising out of or related to any download or use of
            HLA-PANDORA.<br>
            We offer no warranty that:
            <ul style="list-style-type: square; color: blue; text-align: justify;">
                <li>the application meets your requirements</li>
                <li>the application is free of errors or omissions</li>
                <li>access to the services will be uninterrupted, timely, secure, or error-free</li>
                <li>results obtained from the application will be accurate or reliable</li>
                <li>the application will be compatible with all hardware and software you may use</li>
                <li>the sites or services on which this application is hosted will be free of viruses or harmful code</li>
                <li>the application will be free from any claim of infringement of any third party intellectual or rights</li>
                <li>decisions or conclusions made pursuant to your use of this application will meet your expectations, and</li>
                <li>any errors on this application will be corrected.</li>
            </ul>
            This application was created with a &#x1F61C; and &#10084; by <a href="https://www.linkedin.com/in/fmobegi/"
                target="_blank">Fredrick M. Mobegi, PhD</a>.
        </div>',
         type = "warning",
         closeOnEsc = TRUE,
         closeOnClickOutside = TRUE,
         showConfirmButton = TRUE,
         html = TRUE,
         size = "m"
      )
   })
}
