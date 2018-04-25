library(shiny)
library(ggplot2)
## Uncomment the following 3 lines before publishing the app
#library(devtools)
#install_github("gustavdelius/mizer", ref="scaling")
#library(mizer)
library(progress)

server <- function(input, output, session) {
    
    p_bg  <- reactive({
        setBackground(set_scaling_model(no_sp = input$no_sp, 
                        min_w_inf = 10, max_w_inf = 1e5,
                        min_egg = 1e-4, min_w_mat = 10^(0.4),
                        knife_edge_size = 100, kappa = 500))
    })
    
    p <- reactive({
        
        ######### add mullet
        a_m <- 0.0085
        b_m <- 3.11
        L_inf_m <- 24.3
        # http://www.fishbase.org/summary/Mullus-barbatus+barbatus.html
        L_mat <- 11.1
        species_params <- data.frame(
            species = "Mullet",
            w_min = 0.001, # mizer's default egg weight, used in NS
            w_inf = a_m*L_inf_m^b_m, # from fishbase
            w_mat = a_m*L_mat^b_m, # from fishbase
            beta = 283, # = beta_gurnard from North sea. Silvia says gurnard is similar.
            sigma = 1.8, # = sigma_gurnard from North sea. Silvia says gurnard is similar.
            z0 = 0,
            alpha = 0.6, # unknown, mizer default=0.6
            erepro = 0.1, # unknown
            sel_func = "knife_edge", # not used but required
            knife_edge_size = 100, # we can choose
            gear = "knife_edge_gear",
            k = 0,
            r_max = 10^50,
            k_vb = 0.6,
            a = a_m,
            b = b_m,
            h = input$h_mullet
        )
        p <- addSpecies(p_bg(), species_params, SSB = input$SSB_mullet)
        
        ############# add hake 
        # Merluccius merluccius  (European hake)
        a <- 0.0046
        b <- 3.12
        L_inf <- 81.2
        L_mat <- 29.83
        
        species_params <- data.frame(
            species = "Hake",
            w_min = 0.001, # mizer default
            w_inf = a*L_inf^b, # from fishbase
            w_mat = a*L_mat^b, # from fishbase
            beta = exp(2.4), #RLD and Blanchard thesis p 88
            sigma = 1.1, #RLD and Blanchard thesis p 88
            z0 = 0,
            alpha = 0.6, # unknown, using mizer default=0.6
            erepro = 0.1, # unknown
            sel_func = "knife_edge", # not used but required
            knife_edge_size = 100, # can choose
            gear = "knife_edge_gear",
            k = 0,
            r_max = 10^50, #why do I need r_max after combining before
            k_vb = 0.1, # from FB website below
            a = a,
            b = b,
            h = input$h_hake
        )
        addSpecies(p, species_params, SSB = input$SSB_hake)
    })
    
    s <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(p(), t_max = 15, effort = 0, 
                shiny_progress = progress)
    })
    
    output$plotBiomass <- renderPlot({
        plotBiomass(s())
    })
    
    output$plotSpectra <- renderPlot({
        plotSpectra(p())
    })
    
    output$plotGrowthCurveMullet <- renderPlot({
        plotGrowthCurves(p(), species=c("Mullet"))
    })
    
    output$plotGrowthCurveHake <- renderPlot({
        plotGrowthCurves(p(), species=c("Hake"))
    })
    
    output$erepro <- renderPrint({
        p()@species_params$erepro
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Mullet and hake in background ecosystem"),
    
    sidebarLayout(
        
        sidebarPanel(
            tabsetPanel(type = "pills",
                tabPanel("Background",
                         sliderInput("no_sp", "Number of species",
                                     value=10, min=4, max=20, step=1,
                                     round = TRUE)
                        ),
                tabPanel("Mullet",
                         sliderInput("h_mullet", "max feeding rate",
                                     value=62, min=10, max=100, step=2),
                         sliderInput("SSB_mullet", "SSB",
                                     value=100, min=1, max=400)
                         ),
                tabPanel("Hake",
                         sliderInput("h_hake", "max feeding rate",
                                     value=30, min=10, max=100, step=2),
                         sliderInput("SSB_hake", "SSB",
                                     value=200, min=1, max=400)
                         )
            ),
            textOutput("erepro")
        ),  # endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("Spectra", plotOutput("plotSpectra")),
                tabPanel("Total Biomass", plotOutput("plotBiomass")),
                tabPanel("Growth", 
                    fluidRow(column(6, plotOutput("plotGrowthCurveMullet")),
                             column(6, plotOutput("plotGrowthCurveHake"))
                             )
                )
            )
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
