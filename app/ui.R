#!/usr/bin Rscript
ui <- fluidPage(
   tags$head(
      tags$link(rel = "icon", sizes = "32x32", href = "favicon.ico"),
      includeCSS("www/appstyle.css"),
      tags$style(
         HTML(
            "
            #title-container {
                display: flex;
                flex-direction: column;
                align-items: center;
                width: 100%;
                margin-bottom: 20px;
            }
            #logo-row {
                display: flex;
                justify-content: space-between;
                align-items: center;
                width: 100%;
                padding: 10px 0;
            }
            #logo, #uwa_logo {
                height: 100px;
                width: auto;
            }
            #title-text {
                text-align: center;
                width: 100%;
                margin-top: 10px;
                font-family: 'Arial', sans-serif;
                font-weight: bold;
            }
        "
         )
      )
   ),

   ## App title ----
   titlePanel(
      title = div(
         id = "title-container",
         div(
            id = "logo-row",
            tags$a(
               href = "https://pathwest.health.wa.gov.au/Our-Services/Clinical-Services/Immunology",
               # URL for PathWest logo
               target = "_blank",
               tags$img(src = "PathWestlogo.png", id = "logo")
            ),
            tags$a(
               href = "https://www.uwa.edu.au",
               # URL for UWA logo
               target = "_blank",
               tags$img(src = "UWA_logo.png", id = "uwa_logo")
            )
         ),
         div(id = "title-text", h1(
            "HLA-EP-RESOLVER: DONOR-RECIPIENT EPLET ANALYSIS"
         ))
      ),
      windowTitle = "EpletMatch"
   ),
   theme = bslib::bs_theme(bootswatch = "minty"),

   ## App side panel ----
   sidebarLayout(
      sidebarPanel(
         width = 2,
         position = c("left", "justify"),
         fluid = TRUE,
         helpText(a(
            HTML(
               "A program for determining alloreactivity or tolerance to rare, missing and novel HLA alleles.<br><br>
               For detailed usage information, please use the help link below"
            )
         )),
         helpText(h5(HTML("Select input files (Max 50Mb)"))),
         helpText(HTML("Acceptable formats: .xls, .xlsx, .csv, .txt, and .tsv<br>")),

         ## Required inputs (MFA and Genotype for patients)
         helpText(h4(HTML("Required (recipient)"))),
         fileInput(
            "mfi_upload",
            NULL,
            buttonLabel = "MFi",
            multiple = FALSE
         ),
         fileInput(
            "genotype_upload",
            NULL,
            buttonLabel = "Genotype",
            multiple = FALSE
         ),
         tableOutput("output"),

         # Optional input (Donor of Interest)
         div(
            style = "margin: 0;",
            helpText(h4(HTML("Optional (donor)"))),
            tags$div(
               class = "checkbox",
               tags$input(type = "checkbox", id = "donor_interest_checkbox", style = "margin: 0; color:red"),
               tags$label(
                  "Donor of interest ",
                  `for` = "donor_interest_checkbox",
                  tags$i(
                     class = "glyphicon glyphicon-info-sign slim-info-icon", # Added class for styling
                     title = "If selected, provide HLA typing (genotype) for a donor of interest to be processed alongside the other inputs"
                  )
               )
            )
         ),
         conditionalPanel(
            condition = "input.donor_interest_checkbox == true",
            fileInput(
               "donor_of_interest_genotype",
               NULL,
               buttonLabel = "Donor HLA",
               multiple = FALSE
            )
         ),
         # Define the reset and submit buttons
         br(),
         fluidRow(
            column(6, actionButton(
               inputId = "submit",
               label = "Submit",
               class = "btn-success"
            )),
            column(6, actionButton(
               inputId = "resetButton",
               label = "Reset"
            ))
         ),
         br(),
         fluidRow(
            column(12, actionLink(
               inputId = "help",
               label = HTML('<span style="color: blue;">&quest;Help</span>')
            ))
         ),
         fluidRow(
            column(12, actionLink(
               inputId = "terms_of_use_link",
               label = HTML('<span style="color: blue;">&#128712;Terms of Use</span>')
            ))
         )
      ),
      ## Main panel ----
      mainPanel(
         fluid = TRUE,
         width = 10,
         position = c("left", "right"),
         padding = 0,
         tabsetPanel(
            id = "Dataset",
            # tabPanel("MFI", withSpinner(dataTableOutput("table_display"))),
            tabPanel("Genotype", withSpinner(dataTableOutput("geno_display"))),
            tabPanel("Cleaned MFI", withSpinner(dataTableOutput("cleaned_table"))),
            # tabPanel(
            #    "Organ Match",
            #    withSpinner(dataTableOutput("organ_match_tb")),
            #    conditionalPanel(
            #       condition = "output.organ_match_tb !== null && output.organ_match_tb !== undefined && $('#organ_match_tb').text().trim() !== ''",
            #       div(
            #          style = "text-align: right;",
            #          helpText("Click to download organ match"),
            #          downloadButton("downloadData0", "Download tsv")
            #       )
            #    )
            # ),
            # tabPanel(
            #    "CWD Match",
            #    withSpinner(dataTableOutput("cwd_table")),
            #    conditionalPanel(
            #       condition = "output.cwd_table !== null && output.cwd_table !== undefined && $('#cwd_table').text().trim() !== ''",
            #       div(
            #          style = "text-align: right;",
            #          helpText("Click to download CWD summary"),
            #          downloadButton("downloadData2", "Download tsv")
            #       )
            #    )
            # ),
            tabPanel(
               "Eplets (Donor & Bead)",
               withSpinner(dataTableOutput("eplets")),
               conditionalPanel(
                  condition = "output.eplets !== null && output.eplets !== undefined && $('#eplets').text().trim() !== ''",
                  div(
                     style = "text-align: right;",
                     helpText(
                        "Click to download shared donor and bead eplets minimal (no position)"
                     ),
                     downloadButton("downloadData4xls", "Download xlsx")
                  )
               )
            ),
            tabPanel(
               "Eplets (Recipient & Donor)",
               withSpinner(dataTableOutput("deduplicated")),
               conditionalPanel(
                  condition = "output.deduplicated !== null && output.deduplicated !== undefined && $('#deduplicated').text().trim() !== ''",
                  div(
                     style = "text-align: right;",
                     helpText(
                        "Click to download shared recipient and donor eplets minimal (no position)"
                     ),
                     downloadButton("downloadData7xls", "Download xlsx")
                  )
               )
            ),
            tabPanel(
               "Reactive Eplets",
               withSpinner(dataTableOutput("shortlisted_eplets")),
               conditionalPanel(
                  condition = "output.shortlisted_eplets !== null && output.shortlisted_eplets !== undefined && $('#shortlisted_eplets').text().trim() !== ''",
                  div(
                     style = "text-align: right;",
                     helpText(
                        "Click to download shortlisted reactive eplets (recipient and donor)"
                     ),
                     downloadButton("downloadData8xls", "Download xlsx")
                  )
               )
            ),
            tabPanel(
               "Donor of interest",
               withSpinner(dataTableOutput("donor_of_interest")),
               conditionalPanel(
                  condition = "output.donor_of_interest !== null && output.donor_of_interest !== undefined && $('#donor_of_interest').text().trim() !== 'Donor of interest data not provided'",
                  div(
                     style = "text-align: right;",
                     helpText(
                        "Download Donor of interest eplets"
                     ),
                     downloadButton("downloadData10xls", "Download xlsx")
                  )
               )
            ),
            tabPanel(
               "Downloads",
               div(
                  style = "text-align: right;",
                  helpText("Download CWD match"),
                  downloadButton("downloadData1", "Download tsv"),
                  br(),
                  br(),
                  helpText("Download Shared donor and bead eplets"),
                  downloadButton("downloadData3xls", "Download xlsx"),
                  br(),
                  br(),
                  helpText("Download Shared recipient and donor eplets"),
                  downloadButton("downloadData5xls", "Download xlsx"),
                  br(),
                  br(),
                  helpText(
                     "Download shared recipient and donor eplets unfiltered (no position)"
                  ),
                  downloadButton("downloadData6xls", "Download xlsx"),
                  br(),
                  br(),
                  helpText("Download Shortlisted reactive eplets (complete)"),
                  downloadButton("downloadData9xls", "Download xlsx"),
                  br(),
                  br()
               )
            )
         )
      )
   ),
   br(),

   ## App footer ----
   tags$div(
      class = "footer",
      HTML(
         '<div class="footer">
        Copyright &copy; Fredrick Mobegi | Department of Clinical Immunology; PathWest&#x00AE Laboratory Medicine;
        Government of Western Australia Department of Health. <br/>
        This application (HLA-EP-RESOLVER) is for research and reference purposes only and it may contain
        links to embargoed or legally privileged data.
        Except as permitted by the copyright law applicable to you,
        you may not reproduce or communicate any of the content produced on this page,
        including files downloadable from this page, without written permission
        of the copyright owner(s) or authorised PathWest personnel.<br/>
        The user acknowledges that they are using HLA-EP-RESOLVER at their own risk and they agree with the terms of use.
        <br/>
        This application is maintained by
        <a href="https://pathwest.health.wa.gov.au/Our-Services/Clinical-Services/Immunology"
        target="_blank">PathWest&#x00AE Immunology</a>.
        </div>'
      )
   )
)
