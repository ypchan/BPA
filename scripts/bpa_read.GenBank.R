library(shiny)
library(openxlsx)
library(ape)
library(DT)
library(shinyjs)
library(shinythemes)

# Define UI
ui <- fluidPage(
    theme = shinytheme("flatly"),  # Apply theme
    useShinyjs(),  # Use shinyjs
    titlePanel("Download barcode sequences"),
    sidebarLayout(
        sidebarPanel(
            fileInput("file", "Choose Excel File", accept = c(".xlsx")),
            helpText("Upload an Excel file with the taxa table."),
            uiOutput("barcode_ui"),
            helpText("Select the barcodes to process (multiple selection allowed)."),
            textInput("prefix", "Output Prefix", "Sulcatisporaceae"),
            helpText("Enter a prefix for the output file names."),
            actionButton("run", "Run"),
            actionButton("show_table", "Show Table"),
            actionButton("exit", "Exit"),
            verbatimTextOutput("log")
        ),
        mainPanel(
            DTOutput("table")
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    # Read and display Excel file columns
    observeEvent(input$file, {
        req(input$file)
        df <- read.xlsx(input$file$datapath, sheet = 1)
        output$barcode_ui <- renderUI({
            selectInput("barcode", "Select Barcodes", choices = colnames(df), selected = NULL, multiple = TRUE)
        })
        output$table <- renderDT({
            datatable(df, options = list(pageLength = 10, autoWidth = TRUE))
        }, server = FALSE)
    })
    
    # Processing logic
    observeEvent(input$run, {
        req(input$file, input$barcode, input$prefix)
        df <- read.xlsx(input$file$datapath, sheet = 1)
        
        logText <- "Downloaded FASTA flies:\n"
        
        for (barcode in input$barcode) {
            outfilename <- paste(input$prefix, barcode, sep = '_')
            outfilename <- paste(outfilename, 'fasta', sep = '.')
            
            tryCatch({
                df_barcode <- df[!is.na(df[[barcode]]), ]
                barcode_obj <- read.GenBank(unlist(df_barcode[[barcode]]), chunk.size = 50, species.names = FALSE)
                names(barcode_obj) <- df_barcode$longLabel
                write.FASTA(barcode_obj, outfilename)
                
                logText <- paste(logText, paste("    ", outfilename), sep = "\n")
            }, error = function(e) {
                logText <- paste(logText, paste("\nError processing barcode", barcode, ":", e$message), sep = "\n")
            })
        }
        
        output$log <- renderText(logText)
    })
    
    # Display table button
    observeEvent(input$show_table, {
        req(input$file)
        df <- read.xlsx(input$file$datapath, sheet = 1)
        output$table <- renderDT({
            datatable(df, options = list(pageLength = 10, autoWidth = TRUE))
        }, server = FALSE)
    })
    
    # Exit button
    observeEvent(input$exit, {
        stopApp()
    })
}

# Run the application
shinyApp(ui = ui, server = server)
