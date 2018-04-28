library(shiny)
library(ggplot2)
# Uncomment the following 3 lines before publishing the app
#library(devtools)
#install_github("gustavdelius/mizer", ref="scaling")
#library(mizer)
library(progress)

server <- function(input, output, session) {
    
    l50_mullet_old <- 16.6
    l25_mullet_old <- 10
    l50_hake_old <- 16.6
    l25_hake_old <- 10
    
    p_bg  <- reactive({
        setBackground(set_scaling_model(no_sp = input$no_sp, no_w = 400,
                        min_w_inf = 10, max_w_inf = 1e5,
                        min_egg = 1e-4, min_w_mat = 10^(0.4),
                        knife_edge_size = 10^5, kappa = 5000,
                        lambda = input$lambda, f0 = input$f0,
                        h = input$h, r_pp = 10^input$log_r_pp))
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
            sel_func = "sigmoid_length", # not used but required
            gear = "sigmoid_gear",
            l25 = l25_mullet_old,
            l50 = l50_mullet_old,
            k = 0,
            r_max = 10^50,
            k_vb = 0.6,
            a = a_m,
            b = b_m,
            h = input$h_mullet,
            gamma = input$gamma_mullet
        )
        p <- addSpecies(p_bg(), species_params, SSB = input$SSB_mullet, 
                        effort = input$effort, rfac = 1.01)
        
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
            sel_func = "sigmoid_length", # not used but required
            gear = "sigmoid_gear",
            l25 = l25_hake_old,
            l50 = l50_hake_old,
            k = 0,
            r_max = 10^50, #why do I need r_max after combining before
            k_vb = 0.1, # from FB website below
            a = a,
            b = b,
            h = input$h_hake,
            gamma = input$gamma_hake
        )
        p <- addSpecies(p, species_params, SSB = input$SSB_hake, 
                   effort = input$effort, rfac = 1.01)
        
        p
    })
    
    si <- reactive({
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(p(), t_max = 50, t_save = 5, effort = input$effort, 
                shiny_progress = progress)
    })
    
    s <- reactive({
        sim <- si()
        p <- sim@params
        no_sp <- length(p@species_params$species)
        no_t <- dim(sim@n)[1]
        p@initial_n <- sim@n[no_t, , ]
        p@initial_n_pp <- sim@n_pp[no_t, ]
        p@species_params$r_max <- Inf
        
        # Retune the values of erepro so that we get the correct level of
        # recruitment without stock-recruitment relationship
        mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = input$effort)
        gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
        rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
        for (i in (1:no_sp)) {
            gg0 <- gg[i, p@species_params$w_min_idx[i]]
            mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
            DW <- p@dw[p@species_params$w_min_idx[i]]
            p@species_params$erepro[i] <- p@species_params$erepro[i] *
                (p@initial_n[i, p@species_params$w_min_idx[i]] *
                     (gg0 + DW * mumu0)) / rdi[i]
        }
        
        # Set new gear for hake
        a <- p@species_params["Hake", "a"]
        b <- p@species_params["Hake", "b"]
        p@species_params["Hake", "l50"] <- input$l50_hake
        p@species_params["Hake", "l25"] <- input$l25_hake
        p@selectivity["sigmoid_gear", "Hake", ] <-
             sigmoid_length(p@w, input$l25_hake, input$l50_hake, a, b)
        # Set new gear for mullet
        a <- p@species_params["Mullet", "a"]
        b <- p@species_params["Mullet", "b"]
        p@species_params["Mullet", "l50"] <- input$l50_mullet
        p@species_params["Mullet", "l50"] <- input$l25_mullet
        p@selectivity["sigmoid_gear", "Mullet", ] <-
            sigmoid_length(p@w, input$l25_mullet, input$l50_mullet, a, b)
        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(p, t_max = 15, effort = input$effort_new, 
                shiny_progress = progress)
    })
    
    output$plotYield <- renderPlot({
        y <- getYield(s())
        y[1, ] <- getYield(si())[dim(si()@n)[1], ]
        ym <- reshape2::melt(y, varnames = c("Year", "Species"), 
                             value.name = "Yield")
        ym <- subset(ym, ym$Yield > 0)
        ggplot(ym) + 
            geom_line(aes(x=Year, y=Yield, colour=Species, linetype=Species)) +
            scale_y_continuous(name="Yield [g/year]", limits = c(0, NA))
    })
    
    output$plotSSB <- renderPlot({
        b <- getSSB(s())
        bm <- reshape2::melt(b, varnames = c("Year", "Species"), 
                             value.name = "SSB")
        bm <- subset(bm, bm$SSB > 0)
        ggplot(bm) + 
            geom_line(aes(x=Year, y=SSB, colour=Species, linetype=Species)) +
            scale_y_continuous(name="SSB [g]", limits = c(0, NA))
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
        plotGrowthCurves(p(), species=c("Hake"), max_age = 50)
    })
    
    output$plot_erepro <- renderPlot({
        p <- p()
        p@species_params$r_max <- Inf
        # Retune the values of erepro so that we get the correct level of
        # recruitment without stock-recruitment relationship
        mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = input$effort)
        gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
        rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
        # TODO: vectorise this
        no_sp <- length(p@species_params$species)
        erepro <- 1:no_sp  # set up vector of right dimension
        for (i in (1:no_sp)) {
            gg0 <- gg[i, p@species_params$w_min_idx[i]]
            mumu0 <- mumu[i, p@species_params$w_min_idx[i]]
            DW <- p@dw[p@species_params$w_min_idx[i]]
            erepro[i] <- p@species_params$erepro[i] *
                (p@initial_n[i, p@species_params$w_min_idx[i]] *
                     (gg0 + DW * mumu0)) / rdi[i]
        }
        p@species_params$erepro <- erepro
        ggplot(p@species_params, aes(x = species, y = erepro)) + 
            geom_col() + geom_hline(yintercept = 1, color="red")
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Mullet and hake in background ecosystem"),
    
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
                tabPanel("Mullet",
                    sliderInput("SSB_mullet", "SSB",
                                value=1400, min=10, max=2000),
                    sliderInput("w_min_mullet", "Egg weight",
                                value=0.001, min=0.0001, max=0.01),
                    sliderInput("gamma_mullet", "Predation rate coefficient",
                                value=0.0017, min=0.0001, max=0.01),
                    sliderInput("h_mullet", "max feeding rate",
                                value=50, min=10, max=100, step=2),
                    plotOutput("plotGrowthCurveMullet")
                ),
                tabPanel("Hake",
                    sliderInput("SSB_hake", "SSB",
                                value=600, min=10, max=2000),
                    sliderInput("w_min_hake", "Egg weight",
                                value=0.001, min=0.0001, max=0.01),
                    sliderInput("gamma_hake", "Predation rate coefficient",
                                value=0.003, min=0.0001, max=0.01),
                    sliderInput("h_hake", "max feeding rate",
                                value=20, min=10, max=100, step=2),
                         plotOutput("plotGrowthCurveHake")
                ),
                tabPanel("Gear",
                    sliderInput("effort_new", "Effort",
                                value=0.4, min=0.3, max=0.5),
                    h1("Hake selectivity"),
                    sliderInput("l50_hake", "L50",
                                value=16.6, min=10, max=20),
                    sliderInput("l25_hake", "L25",
                                value=10, min=5, max=15),
                    h1("Mullet selectivity"),
                    sliderInput("l50_mullet", "L50",
                                value=16.6, min=10, max=20),
                    sliderInput("l25_mullet", "L25",
                                value=10, min=5, max=15)
                )
            )
        ),  # endsidebarpanel
        
        mainPanel(
            plotOutput("plot_erepro", height = "150px"),
            tabsetPanel(type = "tabs",
                tabPanel("Spectra", plotOutput("plotSpectra")),
                tabPanel("Yield", plotOutput("plotYield")),
                tabPanel("SSB", plotOutput("plotSSB")),
                tabPanel("Total Biomass", plotOutput("plotBiomass"))
            )
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
