library(shiny)
library(devtools)
install_github("gustavdelius/mizer")
library(mizer)
library(progress)

server <- function(input, output, session) {

    sim <- reactive({
        q <- input$lambda - 2 + input$n
        params <- set_scaling_model(no_sp = input$no_sp,
                                    h = input$h, alpha = input$alpha,
                                    n = input$n, q = q
                                    )
        project(params, t_max = 5, t_save = 0.1, effort = 0)
        })

    output$plotBiomass <- renderPlot({plotBiomass(sim())})
    output$plotSpectra <- renderPlot({plotSpectra(sim(), total=TRUE, power=2)})
    output$plotGrowthCurves <- renderPlot({plotGrowthCurves(sim())})

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Trait-based model simulation"),
    
    sidebarLayout(
        
        sidebarPanel(
            sliderInput("no_sp", "Number of species",
                        value=11, min=8, max=20, step=1, round=TRUE),
            sliderInput("h", "h",
                        value=30, min=1, max=100, step=5, animate=TRUE),
            sliderInput("alpha", "conversion efficiency",
                        value=0.6, min=0.1, max=0.9, animate=TRUE),
            sliderInput("n", "n",
                        value=2/3, min=0.6, max=0.8, animate=TRUE),
            sliderInput("lambda", "Sheldon exponent",
                        value=2, min=1.8, max=2.2, step=0.02, animate=TRUE)
        ), #endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("Biomass~Time", plotOutput("plotBiomass")),
                tabPanel("Biomass~Size", plotOutput("plotSpectra")),
                tabPanel("Growth Curves", plotOutput("plotGrowthCurves"))
            )
        )#end mainpanel
    )# end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
