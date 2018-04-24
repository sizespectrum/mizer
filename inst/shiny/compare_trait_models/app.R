library(shiny)
#library(devtools)
#install_github("gustavdelius/mizer", ref="scaling")
#library(mizer)
library(progress)

server <- function(input, output, session) {
    
    params_trait <- reactive({
        set_trait_model(no_sp = 11, min_w_inf = 10, max_w_inf = 1e4)
    })

    sim_trait_initial <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        project(params_trait(), t_max = 100, effort = 0)
    })
    
    sim_trait <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        sti <- sim_trait_initial()
        params <- params_trait()
        params@species_params$knife_edge_size = exp(input$log_knife_edge_size)
        project(params, 
                initial_n = sti@n[dim(sti@n)[1],,],
                initial_n_pp = sti@n_pp[dim(sti@n_pp)[1],], 
                t_max = 15, effort = input$effort, 
                shiny_progress = progress)
    })
    
    sim <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        params <- set_scaling_model(no_sp = 11, 
            min_w_inf = 10, max_w_inf = 1e4,
            knife_edge_size = exp(input$log_knife_edge_size),
            min_egg = 1e-4, min_w_mat = 10^(0.4), rfac = Inf,
            kappa = 0.01)
        
        project(params, t_max = 15, effort = input$effort, 
                shiny_progress = progress)
    })
    
    output$plotSSB <- renderPlot({
        b_trait <- getSSBFrame(sim_trait())
        b <- getSSBFrame(sim())
        display_frames(b_trait, b)
    })
    
    output$plotBiomass <- renderPlot({
        b_trait <- getBiomassFrame(sim_trait())
        b <- getBiomassFrame(sim())
        display_frames(b_trait, b)
    })
    
    output$plotYield <- renderPlot({
        if (input$effort > 0) {
            plotYield(sim_trait(), sim())
        } else {
            ggplot()
        }
    })
    
    output$plotSpectra1 <- renderPlot({
        plotSpectra(sim_trait(), time_range = input$year)
    })
    output$plotSpectra2 <- renderPlot({
        plotSpectra(sim(), time_range = input$year)
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Trait-based model simulation"),
    
    sidebarLayout(
        
        sidebarPanel(
            sliderInput("effort", "Fishing effort",
                        value=0, min=0, max=2, step=0.1),
            sliderInput("log_knife_edge_size", "Log10 of minimum fishing size",
                        value=2, min=1, max=3, ticks=FALSE)
        ), #endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("SSB", plotOutput("plotSSB")),
                tabPanel("Yield", plotOutput("plotYield")),
                tabPanel("Total Biomass", plotOutput("plotBiomass")),
                tabPanel("Spectra", 
                         wellPanel(
                             sliderInput("year", "Year",
                                         value = 0, min = 0, max = 15, step = 1, 
                                         animate = TRUE)
                         ),
                         fluidRow(column(6, plotOutput("plotSpectra1")),
                                  column(6, plotOutput("plotSpectra2"))
                         )
                )
            )
        )#end mainpanel
    )# end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
