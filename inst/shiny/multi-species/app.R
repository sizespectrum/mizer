library(shiny)
library(ggplot2)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="scaling")
# library(mizer)
library(progress)

server <- function(input, output, session) {
    
    params_data <- read.csv(system.file("extdata", "speciesNCME_edited2.csv", package = "mizer"))

    no_sp <- dim(params_data)[1]
    
    l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
    l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
    names(l25) <- as.character(params_data$species)
    names(l50) <- as.character(params_data$species)
    
    effort <- 1.4
    
    p_bg  <- reactive({
        setBackground(
            set_scaling_model(
                min_w_pp = 1e-12, no_sp = 10, no_w = 400, 
                min_w_inf = 2, max_w_inf = 6e5,
                min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
                lambda = 2.12, knife_edge_size = Inf
            )
        )
    })
    
    p <- reactive({
        p <- p_bg()
        for (i in (1:no_sp)) {
            a_m <- params_data$a2[i]
            b_m <- params_data$b2[i]
            L_inf_m <- params_data$Linf[i]
            L_mat <- params_data$Lmat[i]
            species_params <- data.frame(
                species = as.character(params_data$species[i]),
                w_min = params_data$Wegg[i],
                w_inf = params_data$w_inf[i],
                w_mat = params_data$w_mat[i],
                beta = params_data$beta[i],
                sigma = log(params_data$sigma[i]),
                z0 = 0,
                alpha = 0.6,
                erepro = 0.1, # unknown, determined later
                sel_func = "sigmoid_length",
                gear = "sigmoid_gear",
                l25 = l25[i],
                l50 = l50[i],
                k = 0,
                k_vb = params_data$k_vb[i],
                a = a_m,
                b = b_m
            )
            
            p <- addSpecies(p, species_params, effort = effort, rfac=Inf)
        }
        p
    })
    
    p_steady <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        steady(p(), effort = effort, t_max = 100, tol = 1e-2)
    })
    
    s <- reactive({
        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(p, t_max = 15, t_save = 0.1, effort = input$effort, 
                shiny_progress = progress)
    })
    
    output$params_sliders <- renderUI({
        sp <- p()@species_params[input$species_sel, ]
        list(sliderInput("gamma", "Predation rate coefficient gamma",
                    value = sp$gamma, 
                    min = sp$gamma/10, 
                    max = sp$gamma*2),
        sliderInput("h", "max feeding rate h",
                    value = sp$h, 
                    min = sp$h/10, 
                    max = sp$h*2))
    })
    output$species_sel <- renderUI({
        selectInput("species_sel", "Species:", as.character(params_data$species)) 
    })
    
    output$plotGrowthCurve <- renderPlot({
        plotGrowthCurves(p_steady(), species=input$species_sel)
    })
    
    output$plot_erepro <- renderPlot({
        ggplot(p_steady()@species_params, aes(x = species, y = erepro)) + 
            geom_col() + geom_hline(yintercept = 1, color="red")
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Humboldt current ecosystem"),
    
    sidebarLayout(
        
        sidebarPanel(
            tabsetPanel(
                tabPanel("General",
                    sliderInput("lambda", "Sheldon exponent",
                                value=2.08, min=1.9, max=2.2, step=0.005),
                    sliderInput("f0", "Feeding level",
                                value=0.6, min=0, max=1),
                    sliderInput("h", "max feeding rate",
                                value=34, min=10, max=100, step=2),
                    sliderInput("log_r_pp", "log10 Plankton replenishment rate",
                                value=-2, min=-4, max=0),
                    sliderInput("effort", "Fishing effort",
                                value=0.4, min=0, max=2, step=0.1),
                    sliderInput("no_sp", "Number of species",
                                value=10, min=4, max=20, step=1, round = TRUE)
                ),
                tabPanel("Species",
                    uiOutput("species_sel"),
                    uiOutput("params_sliders")
                    #plotOutput("plotGrowthCurve")
                ),
                tabPanel("Gear",
                    sliderInput("effort", "Effort",
                                value=1.4, min=0.3, max=0.5)
                )
            )
        ),  # endsidebarpanel
        
        mainPanel(
            plotOutput("plot_erepro", height = "150px"),
            tabsetPanel(type = "tabs",
                tabPanel("Growth", plotOutput("plotGrowthCurve")),
                tabPanel("Yield", plotOutput("plotYield")),
                tabPanel("SSB", plotOutput("plotSSB")),
                tabPanel("Total Biomass", plotOutput("plotBiomass"))
            )
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
