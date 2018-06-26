library(shiny)
library(ggplot2)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="scaling")
# library(mizer)
library(progress)

server <- function(input, output, session) {
    
    l50 <- c(12.5, 15.48, 15.48)
    names(l50) <- c("anchovy", "jmackerel", "jmackerelB")
    sd <- c(0.462, 2.10, 2.10)
    l25 = l50 - log(3) * sd
    
    
    # do we want to make ks variable too ?
    
    p_bg  <- reactive({
        setBackground(set_scaling_model(no_sp = input$no_sp, no_w = 400,
                        min_w_inf = 10, max_w_inf = 1e5,
                        min_egg = 1e-4, min_w_mat = 10^(0.4),
                        knife_edge_size = Inf, kappa = 10000,
                        lambda = input$lambda, f0 = input$f0,
                        h = input$h, r_pp = 10^input$log_r_pp))
    })
    
    p <- reactive({
        
        ######### add anchovy
        a_m <- 0.005
        b_m <- 3.17
        L_inf_m <- 20.25
        L_mat <- input$L_mat_anchovy
        species_params <- data.frame(
            species = "Anchovy",
            w_min = 10^-3,
            w_inf = a_m*L_inf_m^b_m,
            w_mat = a_m*L_mat^b_m,
            beta = exp(6.7),
            sigma = 2.32,
            z0 = 0,
            alpha = input$alpha_anchovy, # unknown, mizer default=0.6
            erepro = 0.1, # unknown, determined later
            sel_func = "sigmoid_length",
            gear = "sigmoid_gear",
            l25 = l25["anchovy"],
            l50 = l50["anchovy"],
            k = 0,
            k_vb = 0.88,
            a = a_m,
            b = b_m,
            gamma = input$gamma_anchovy,
            h = input$h_anchovy,
            linecolour = "red",
            linetype = "solid"
        )
        p <- addSpecies(p_bg(), species_params, SSB = input$SSB_anchovy, 
                        effort = input$effort, rfac = 1.01)
        
        ############# add Jack Mackerel 
        a <- 0.01
        b <- 3.058
        L_inf <- 70.80
        L_mat <- input$L_mat_jmackerel
        
        species_params <- data.frame(
            species = "Jack Mackerel",
            w_min = 10^-3,
            w_inf = a*L_inf^b, 
            w_mat = a*L_mat^b, 
            beta = 26.76,
            sigma = 2.83, 
            z0 = 0,
            alpha = input$alpha_jmackerel, # unknown, using mizer default=0.6
            erepro = 0.1, # unknown
            sel_func = "sigmoid_length", # not used but required
            gear = "sigmoid_gear",
            l25 = l25["jmackerel"],
            l50 = l50["jmackerel"],
            k = 0,
            k_vb = 0.094,
            a = a,
            b = b,
            gamma = input$gamma_jmackerel,
            h = input$h_jmackerel,
            linecolour = "blue",
            linetype = "solid"
        )
        p <- addSpecies(p, species_params, SSB = input$SSB_jmackerel, 
                   effort = input$effort, rfac = 1.01)

    ############# add Jack MackerelB 
    a <- 0.01
    b <- 3.058
    L_inf <- 70.80
    L_mat <- input$L_mat_jmackerelB
    
    species_params <- data.frame(
      species = "Jack MackerelB",
      w_min = 10^-3,
      w_inf = a*L_inf^b, 
      w_mat = a*L_mat^b, 
      beta = 26.76,
      sigma = 2.83, 
      z0 = 0,
      alpha = input$alpha_jmackerelB, # unknown, using mizer default=0.6
      erepro = 0.1, # unknown
      sel_func = "sigmoid_length", # not used but required
      gear = "sigmoid_gear",
      l25 = l25["jmackerelB"],
      l50 = l50["jmackerelB"],
      k = 0,
      k_vb = 0.094,
      a = a,
      b = b,
      gamma = input$gamma_jmackerelB,
      h = input$h_jmackerelB,
      linecolour = "blue",
      linetype = "solid"
    )
    p <- addSpecies(p, species_params, SSB = input$SSB_jmackerelB, 
                    effort = input$effort, rfac = 1.01)
    p
})


########################
    
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
        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(p, t_max = 15, t_save = 0.1, effort = input$effort_new, 
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
            scale_y_continuous(name="Yield [g/year]", limits = c(0, NA)) +
            scale_colour_manual(values = s()@params@linecolour) +
            scale_linetype_manual(values = s()@params@linetype)
    })
    
    output$plotSSB <- renderPlot({
        b <- getSSB(s())
        bm <- reshape2::melt(b, varnames = c("Year", "Species"), 
                             value.name = "SSB")
        bm <- subset(bm, bm$SSB > 0)
        ggplot(bm) + 
            geom_line(aes(x=Year, y=SSB, colour=Species, linetype=Species)) +
            scale_y_continuous(name="SSB [g]", limits = c(0, NA)) +
            scale_colour_manual(values = s()@params@linecolour) +
            scale_linetype_manual(values = s()@params@linetype)
    })
    
    output$plotBiomass <- renderPlot({
        plotBiomass(s())
    })
    
    output$plotSpectra <- renderPlot({
        plotSpectra(p(), total=TRUE)
    })
    
    output$plotGrowthCurveAnchovy <- renderPlot({
        plotGrowthCurves(p(), species=c("Anchovy"))
    })
    
    output$plotGrowthCurveJmackerel <- renderPlot({
        plotGrowthCurves(p(), species=c("Jack Mackerel"), max_age = 12)
    })
    
    output$plotGrowthCurveJmackerelB <- renderPlot({
      plotGrowthCurves(p(), species=c("Jack MackerelB"), max_age = 12)
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
                tabPanel("Anchovy",
                    sliderInput("SSB_anchovy", "SSB",
                                value=2800, min=100, max=10000),
                    sliderInput("gamma_anchovy", "Predation rate coefficient",
                                value=0.00085, min=0.00005, max=0.005),
                    sliderInput("h_anchovy", "max feeding rate",
                                value=50, min=10, max=100, step=2),
                    sliderInput("alpha_anchovy", "Conversion efficiency",
                                value=0.6, min=0.1, max=0.8, step=0.05),
                    sliderInput("L_mat_anchovy", "Maturity length",
                                value=12.15, min=8, max=18, step=0.05),
                    plotOutput("plotGrowthCurveAnchovy")
                ),
                tabPanel("Jack Mackerel",
                    sliderInput("SSB_jmackerel", "SSB",
                                value=1200, min=100, max=4000),
                    sliderInput("gamma_jmackerel", "Predation rate coefficient",
                                value=0.0015, min=0.00005, max=0.005),
                    sliderInput("h_jmackerel", "max feeding rate",
                                value=14, min=6, max=30, step=2),
                    sliderInput("alpha_jmackerel", "Conversion efficiency",
                                value=0.6, min=0.1, max=0.8, step=0.05),
                    sliderInput("L_mat_jmackerel", "Maturity length",
                                value=25.5, min=20, max=38, step=0.05),
                    plotOutput("plotGrowthCurveJmackerel")
                ),
                tabPanel("Jack MackerelB",
                         sliderInput("SSB_jmackerelB", "SSB",
                                     value=1200, min=100, max=4000),
                         sliderInput("gamma_jmackerelB", "Predation rate coefficient",
                                     value=0.0015, min=0.00005, max=0.005),
                         sliderInput("h_jmackerelB", "max feeding rate",
                                     value=14, min=6, max=30, step=2),
                         sliderInput("alpha_jmackerelB", "Conversion efficiency",
                                     value=0.6, min=0.1, max=0.8, step=0.05),
                         sliderInput("L_mat_jmackerelB", "Maturity length",
                                     value=25.5, min=20, max=38, step=0.05),
                         plotOutput("plotGrowthCurveJmackerelB")
                ),
                tabPanel("Gear",
                    sliderInput("effort_new", "Effort",
                                value=0.4, min=0.3, max=0.5)
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
