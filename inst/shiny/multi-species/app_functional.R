library(shiny)
library(ggplot2)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="scaling")
# library(mizer)
library(progress)

server <- function(input, output, session) {
    
    ## Load params object and store it as a reactive value
    data("humboldt_params")
    params <- reactiveVal()
    params(humboldt_params)
    
    ## Create dynamic ui for species parameters
    output$sp_sel <- renderUI({
        p <- isolate(params())
        species <- as.character(p@species_params$species[!is.na(p@A)])
        selectInput("sp_sel", "Species:", species) 
    })
    output$params_sliders <- renderUI({
        req(input$sp_sel)
        sp <- params()@species_params[input$sp_sel, ]
        list(numericInput("gamma", "Predation rate coefficient gamma",
                          value = sp$gamma,
                          min = sp$gamma/10,
                          max = sp$gamma*2,
                          step = 1),
             numericInput("h", "max feeding rate h",
                          value = sp$h,
                          min = sp$h/10,
                          max = sp$h*2),
             numericInput("alpha", "Assimilation efficiency alpha",
                          value = sp$alpha,
                          min = 0,
                          max = 1, step = 10^(-2)),
             numericInput("ks", "Coefficient of standard metabolism ks",
                          value = sp$ks,
                          min = sp$ks/10,
                          max = sp$ks*2, 
                          step = 10^(-2)),
             numericInput("beta", "Preferred predator-prey mass ratio beta",
                          value = sp$beta,
                          min = 1.01,
                          max = sp$beta*100, 
                          step = 10^(-2)),
             numericInput("sigma", "Width of size selection function sigma",
                          value = sp$sigma,
                          min = sp$sigma/10,
                          max = sp$sigma*100, 
                          step = 10^(-2)),
             numericInput("k_vb", "von Bertalanffy parameter k_vb",
                          value = sp$k_vb,
                          min = sp$k_vb/10,
                          max = sp$k_vb*100, 
                          step = 10^(-2)),
             numericInput("a", "Coefficient for length to weight conversion a",
                          value = sp$a,
                          min = sp$a/10,
                          max = sp$a*10, 
                          step = 10^(-4)),
             numericInput("b", "Exponent for length to weight conversion b",
                          value = sp$b,
                          min = sp$b/10,
                          max = sp$b*100, 
                          step = 10^(-2))
             )
    })
    ## When a species parameter input is changed, only change that
    # parameter value in the params object and run to steady state
    # starting at previous steady state
    observeEvent(input$sp_go, {
        req(input$gamma, input$h)

        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())

        p <- params()
        # Create updated species params data frame
        species_params <- p@species_params
        species_params[input$sp_sel, "gamma"] <- input$gamma
        species_params[input$sp_sel, "h"]     <- input$h
        species_params[input$sp_sel, "alpha"] <- input$alpha
        species_params[input$sp_sel, "ks"]    <- input$ks
        species_params[input$sp_sel, "beta"]  <- input$beta
        species_params[input$sp_sel, "sigma"] <- input$sigma
        species_params[input$sp_sel, "k_vb"]  <- input$k_vb
        species_params[input$sp_sel, "a"]     <- input$a
        species_params[input$sp_sel, "b"]     <- input$b
        # Create new params object identical to old one except for changed
        # species params 
        pc <- MizerParams(
            species_params,
            p = p@p,
            n = p@n,
            q = p@q,
            lambda = p@lambda,
            f0 = p@f0,
            kappa = p@kappa,
            min_w = min(p@w),
            max_w = max(p@w),
            no_w = length(p@w),
            min_w_pp = min(p@w_full),
            w_pp_cutoff = max(p@w_full),
            r_pp = (p@rr_pp / (p@w_full ^ (p@p - 1)))[1]
        )
        pc@linetype <- p@linetype
        pc@linecolour <- p@linecolour
        pc@A <- p@A
        pc@sc <- p@sc
        pc@cc_pp <- p@cc_pp
        pc@mu_b <- p@mu_b

        pc@initial_n <- p@initial_n
        pc@initial_n_pp <- p@initial_n_pp

        # Run to steady state
        #pc <- steady(pc, effort = input$effort, t_max = 100, tol = 1e-2,
        #             shiny_progress = progress)
        pc <- steady(pc, effort = c(knife_edge_gear = 0, sigmoid_gear = input$effort, sigmoid_gear_Anchovy = input$Anchovy_effort),
                     t_max = 100, tol = 1e-2,
                     shiny_progress = progress)
        # Update the reactive params object
        params(pc)
    })

    ## If a general parameter changes the entire params object is recalculated
    # from scratch
    observeEvent(input$bg_go, {

        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())

        p_old <- params()

        p <- setBackground(
            set_scaling_model(
                #min_w_pp = input$min_w_pp, no_sp = input$no_bg_sp, no_w = 400,
              min_w_pp = 1e-12, no_sp = input$no_bg_sp, no_w = 400,
                min_w_inf = 2, max_w_inf = 6e5,
                min_egg = 1e-4, min_w_mat = 2 / 10^0.6,
                lambda = input$lambda, knife_edge_size = Inf,
                f0 = input$f0, h = input$h_bkgd, r_pp = 10^input$log_r_pp
            )
        )
        # Loop over all foreground species and add them one-by-one to the new
        # background
        effort <- 0
        names(effort) <- "knife_edge_gear"
        all_efforts <- c(0, input$effort, input$Anchovy_effort)
        names(all_efforts) <- c("knife_edge_gear", "sigmoid_gear", "sigmoid_gear_Anchovy")
        no_sp <- length(p_old@A)
        for (sp in (1:no_sp)[!is.na(p_old@A)]) {
            if (!(p_old@species_params[sp, "gear"] %in% names(effort))) {
                effort <- c(effort, all_efforts[p_old@species_params[sp, "gear"]])
            }
            p <- addSpecies(p, p_old@species_params[sp, ],
                              effort = effort,
                          rfac=Inf)
        }

        # Run to steady state
        p <- steady(p, effort = c(knife_edge_gear = 0, sigmoid_gear = input$effort, sigmoid_gear_Anchovy = input$Anchovy_effort), 
                    t_max = 100, tol = 1e-2,
                    shiny_progress = progress)
        # Update the reactive params object
        params(p)
    })
    
    # Run a simulation
    sim <- reactive({
        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(params(), t_max = 15, t_save = 0.1, 
                effort = c(knife_edge_gear = 0, sigmoid_gear = input$effort, 
                           sigmoid_gear_Anchovy = input$Anchovy_effort), 
                shiny_progress = progress)
    })
    
    # Run a simulation with homogenous fishing efforts
    sim_new <- reactive({
      
      # Create a Progress object
      progress <- shiny::Progress$new(session)
      on.exit(progress$close())
      
      project(params(), t_max = 15, t_save = 0.1, 
              effort = c(knife_edge_gear = 0, sigmoid_gear = input$effort, 
                         sigmoid_gear_Anchovy = input$new_Anchovy_effort), 
              shiny_progress = progress)
    })
    
    output$plotGrowthCurve <- renderPlot({
        plotGrowthCurves(params(), species = input$sp_sel)
    })
    
    output$plotSpectra <- renderPlot({
        plotSpectra(params(), total=TRUE)
    })
    
    output$plotBiomass <- renderPlot({
        plotBiomass(sim())
    })
    
    output$plotSSB <- renderPlot({
      #b_new <- getSSB(sim_new())[, "Anchovy"]
      #b <- getSSB(sim())[, "Anchovy"]
      #plot((1:length(b_new))/10, b_new, sub="solid => Anchovy under old effort, dashed => Anchovy under new effort",xlab ="Time",
      #      ylab="SSB",type="l",lty=2)
      #lines((1:length(b))/10, b)
      
      
      b <- getSSB(sim())[, "Anchovy"]
      b_new <- getSSB(sim_new())[, "Anchovy"]
      b_df <- data.frame(
        "Year" = rep((1:length(b))/10,2),
        "Species" = rep(rep("Anchovy",length(b)),2),
        "SSB" = c(b,b_new),
        "Effort" = c(rep("Default",length(b)),rep("New",length(b)))
      )
      ggplot(b_df) + 
        geom_line(aes(x = Year, y = SSB, colour = Species, linetype = Effort)) +
            scale_y_continuous(name="SSB [tonnes]", limits = c(0, NA)) +
          #  scale_colour_manual(values = params()@linecolour) +
            scale_linetype_manual(values = c("New" = "dotted", "Default" = "solid")) +
            theme(text = element_text(size = 18))
    })
    
    output$plot_erepro <- renderPlot({
        ggplot(params()@species_params, aes(x = species, y = erepro)) + 
            geom_col() + geom_hline(yintercept = 1, color="red")
    })

} #the server

#### user interface
ui <- fluidPage(
    
    titlePanel("Humboldt current ecosystem"),
    
    sidebarLayout(
        
        sidebarPanel(
            tabsetPanel(
                tabPanel("Species",
                    uiOutput("sp_sel"),
                    actionButton("sp_go", "Go"),
                    uiOutput("params_sliders")
                ),
                tabPanel("General",
                    actionButton("bg_go", "Go"),
                    numericInput("lambda", "Sheldon exponent",
                                value=2.12, min=1.9, max=2.2, step=0.005),
                    sliderInput("f0", "Feeding level",
                                value=0.6, min=0, max=1),
                    sliderInput("h_bkgd", "max feeding rate",
                                value=30, min=10, max=100, step=2),
                    sliderInput("log_r_pp", "log10 Plankton replenishment rate",
                                value=-1, min=-4, max=0),
                    sliderInput("no_bg_sp", "Number of background species",
                                value=10, min=4, max=20, step=1, round = TRUE),
                    numericInput("min_w_pp", "Minimum plankton weight min_w_pp",
                                value=1e-12,  step=1e-13)
                ),
                tabPanel("Fishing @ steady state",
                    sliderInput("effort", "General Effort",
                                value=1.4, min=0.3, max=2),
                    sliderInput("Anchovy_effort", "Anchovy Effort",
                                value=1.1, min=0.3, max=2)
                ),
                tabPanel("Changed Fishing",
                         sliderInput("new_Anchovy_effort", "New Anchovy Effort",
                                     value=1.1, min=0.3, max=2)
                )
            )
        ),  # endsidebarpanel
        
        mainPanel(
            plotOutput("plot_erepro", height = "150px"),
            tabsetPanel(type = "tabs",
                tabPanel("Spectra", plotOutput("plotSpectra")),
                tabPanel("Growth", plotOutput("plotGrowthCurve")),
                tabPanel("Biomass", plotOutput("plotBiomass")),
                tabPanel("SSB Anchovy", plotOutput("plotSSB"))
            )
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
