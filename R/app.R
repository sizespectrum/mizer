library(shiny)
library(mizer)
library(progress)

server <- function(input, output, session) {

    sim <- reactive({
        params <- set_trait_model(no_sp = input$no_sp,
                                  min_w_inf = 10, max_w_inf = 1e5,
                                  knife_edge_size = input$knife_edge_size)
        project(params, t_max = 20, effort = input$effort)
        })

    output$plot <- renderPlot({
        plot(sim())
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Trait-based model simulation"),
    
    sidebarLayout(
        
        sidebarPanel(
            sliderInput("no_sp", "Number of species",
                        value=4, min=1, max=20, step=1, round=TRUE),
            sliderInput("effort", "Fishing effort",
                        value=0.4, min=0, max=2, step=0.1),
            sliderInput("knife_edge_size", "Minimum fishing size",
                        value=100, min=10, max=1000, ticks=FALSE)
        ), #endsidebarpanel
        
        mainPanel(
            plotOutput("plot") 
        )#end mainpanel
    )# end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
