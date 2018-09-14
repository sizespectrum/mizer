library(shiny)
library(ggplot2)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="humboldt")
# library(mizer)
library(progress)

server <- function(input, output, session) {
    effort <- 1.4
    
    ## Load params object and store it as a reactive value ####
    params <- reactiveVal()
    params(readRDS("humboldt_params.rds"))
    
    ## Handle params object uploaded by user ####
    observeEvent(input$upload, {
        inFile <- input$upload
        tryCatch(
            pa <- readRDS(inFile$datapath),
            error = function(e) {stop(safeError(e))}
        )
        validate(
            need(is(pa, "MixerParams"),
                 "The file does not contain a mizer parameter object.")
        )
        params(pa)
    })
    
    output$plotSpectra <- renderPlot({
        plotSpectra(params(), total = TRUE)
    })
    
} #the server

#### User interface ####
ui <- fluidPage(
    
    titlePanel("Humboldt current ecosystem"),
    
    sidebarLayout(
        
        ## Sidebar ####
        sidebarPanel(
            fileInput("upload", "Upload new params object", accept = ".rds")
        ),  # endsidebarpanel
        
        ## Main panel ####
        mainPanel(
            plotOutput("plotSpectra")
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app


