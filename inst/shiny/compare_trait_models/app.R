library(shiny)
library(devtools)
#install_github("gustavdelius/mizer")
#library(mizer)
library(progress)

server <- function(input, output, session) {

    sim_trait <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        params <- set_trait_model(min_w_inf = 10, max_w_inf = 1e5,
                                  knife_edge_size = input$knife_edge_size)
        sim_trait_initial <- project(params_trait, t_max = 50, effort = 0)
        project(params_trait, 
                initial_n = sim_trait_initial@n[dim(sim_trait_initial@n)[1],,],
                initial_n_pp = sim_trait_initial@n_pp[dim(sim_trait_initial@n_pp)[1],], 
                t_max = 15, effort = input$effort, 
                shiny_progress = progress)
        })
    
    sim <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        
        params <- set_scaling_model(
            min_w_inf = 10, max_w_inf = 1e5,
            knife_edge_size = input$knife_edge_size,
            min_egg = 1e-4, min_w_mat = 10^(0.4), rfac = Inf,
            kappa = 0.005)
        project(params, t_max = 15, effort = input$effort, 
                shiny_progress = progress)
    })

    output$plotBiomass <- renderPlot({
        b_trait <- getBiomassFrame(sim_trait())
        b <- getBiomassFrame(sim())
        display_frames(b_trait, b)})
    output$plotYield <- renderPlot({plotYield(sim_trait(), sim())})
    output$plotSpectra <- renderPlot({plotSpectra(sim())})

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Trait-based model simulation"),
    
    sidebarLayout(
        
        sidebarPanel(
            sliderInput("effort", "Fishing effort",
                        value=0.4, min=0, max=2, step=0.1),
            sliderInput("knife_edge_size", "Minimum fishing size",
                        value=100, min=10, max=1000, ticks=FALSE)
        ), #endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("Biomass", plotOutput("plotBiomass")),
                tabPanel("Yield", plotOutput("plotYield"))
            )
        )#end mainpanel
    )# end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
