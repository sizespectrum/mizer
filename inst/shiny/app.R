library(shiny)
library(devtools)
#install_github("gustavdelius/mizer")
#library(mizer)
library(progress)

server <- function(input, output, session) {

    sim <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        params <- set_trait_model(no_sp = input$no_sp,
                                  min_w_inf = 10, max_w_inf = 1e5,
                                  knife_edge_size = input$knife_edge_size)
        project(params, t_max = 200, effort = input$effort, 
                shiny_progress = progress)
        })

    output$plot <- renderPlot({plot(sim())})
    output$plotBiomass <- renderPlot({plotBiomass(sim())})
    output$plotYield <- renderPlot({plotYield(sim())})
    output$plotSpectra <- renderPlot({plotSpectra(sim())})
    output$plotM2 <- renderPlot({plotM2(sim())})

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Trait-based model simulation"),
    
    sidebarLayout(
        
        sidebarPanel(
            sliderInput("no_sp", "Number of species",
                        value=4, min=3, max=12, step=1, round=TRUE),
            sliderInput("effort", "Fishing effort",
                        value=0.4, min=0, max=2, step=0.1),
            sliderInput("knife_edge_size", "Minimum fishing size",
                        value=100, min=10, max=1000, ticks=FALSE)
        ), #endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("Biomass~Time", plotOutput("plotBiomass")),
                tabPanel("Yield~Time", plotOutput("plotYield")),
                tabPanel("Biomass~Size", plotOutput("plotSpectra")),
                tabPanel("Mortality~Size", plotOutput("plotM2")),
                tabPanel("Combo", plotOutput("plot"))
            )
        )#end mainpanel
    )# end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
