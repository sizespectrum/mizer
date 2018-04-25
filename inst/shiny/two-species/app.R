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
                        knife_edge_size = 100, kappa = 500,
                        lambda = input$lambda, gamma = input$gamma,
                        h = input$h))
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
            w_min = input$w_min_mullet,
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
            h = input$h_mullet,
            gamma = input$gamma_mullet
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
            w_min = input$w_min_hake,
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
            h = input$h_hake,
            gamma = input$gamma_hake
        )
        addSpecies(p, species_params, SSB = input$SSB_hake)
    })
    
    s <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(p(), t_max = 15, t_save = 0.2, effort = 0, 
                shiny_progress = progress)
    })
    
    output$plotBiomass <- renderPlot({
        plotBiomass(s())
    })
    
    output$plotSpectra <- renderPlot({
        plotSpectra(p(), total=TRUE)
    })
    
    output$plotGrowthCurveMullet <- renderPlot({
        plotGrowthCurves(p(), species=c("Mullet"))
    })
    
    output$plotGrowthCurveHake <- renderPlot({
        plotGrowthCurves(p(), species=c("Hake", max_age = 50))
    })
    
    output$plot_erepro <- renderPlot({
        ggplot(p()@species_params, aes(x = species, y = erepro)) + 
            geom_col() + geom_hline(yintercept = 1, color="red")
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Mullet and hake in background ecosystem"),
    
    sidebarLayout(
        
        sidebarPanel(
            tabsetPanel(type = "pills",
                tabPanel("Background",
                         sliderInput("lambda", "Sheldon exponent",
                                     value=2.14, min=1.9, max=2.2, step = 0.01),
                         sliderInput("gamma_mullet", "Predation rate coefficient",
                                     value=0.017, min=0.001, max=0.1),
                         sliderInput("h", "max feeding rate",
                                     value=30, min=10, max=100, step=2),
                         sliderInput("no_sp", "Number of species",
                                     value=10, min=4, max=20, step=1,
                                     round = TRUE)
                        ),
                tabPanel("Mullet",
                         sliderInput("SSB_mullet", "SSB",
                                     value=100, min=1, max=400),
                         sliderInput("w_min_mullet", "Egg weight",
                                     value=0.001, min=0.0001, max=0.01),
                         sliderInput("gamma_mullet", "Predation rate coefficient",
                                     value=0.017, min=0.001, max=0.1),
                         sliderInput("h_mullet", "max feeding rate",
                                     value=50, min=10, max=100, step=2)
                         ),
                tabPanel("Hake",
                         sliderInput("SSB_hake", "SSB",
                                     value=200, min=1, max=400),
                         sliderInput("w_min_hake", "Egg weight",
                                     value=0.001, min=0.0001, max=0.01),
                         sliderInput("gamma_hake", "Predation rate coefficient",
                                     value=0.025, min=0.001, max=0.1),
                         sliderInput("h_hake", "max feeding rate",
                                     value=30, min=10, max=100, step=2)
                         )
            ),
            plotOutput("plot_erepro")
        ),  # endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                tabPanel("Spectra", plotOutput("plotSpectra")),
                tabPanel("Growth Curves", 
                    fluidRow(column(6, plotOutput("plotGrowthCurveMullet")),
                             column(6, plotOutput("plotGrowthCurveHake"))
                             )
                ),
                tabPanel("Total Biomass", plotOutput("plotBiomass"))
            )
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
