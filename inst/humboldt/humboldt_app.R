library(shiny)
library(ggplot2)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="humboldt")
# library(mizer)
library(progress)

server <- function(input, output, session) {
  effort <- 1.4
  
  ## Load params object and store it as a reactive value ####
  params <- reactiveVal()
  params(readRDS("humboldt_params.rds"))
  
  ## Handle params object uploaded by user ####
  observeEvent(input$upload, {
    inFile <- input$upload
    tryCatch(
      p <- readRDS(inFile$datapath),
      error = function(e) {stop(safeError(e))}
    )
    validate(
      need(is(p, "MixerParams"),
           "The file does not contain a mizer parameter object.")
    )
    params(p)
  })
  
  ## Prepare for download of params object ####
  output$params <- downloadHandler(
    filename = "params.rds", 
    content = function(file) {
      saveRDS(params(), file = file)
    })
  
  ## Create dynamic ui for species parameters ####
  output$sp_sel <- renderUI({
    p <- isolate(params())
    species <- as.character(p@species_params$species[!is.na(p@A)])
    selectInput("sp_sel", "Species:", species) 
  })
  output$params_sliders <- renderUI({
    req(input$sp_sel)
    sp <- params()@species_params[input$sp_sel, ]
    list(
      # numericInput("rdd", "Reproductive rate",
      #              value = sp$SSB,
      #              min = 10^(-20),
      #              max = 10^4, 
      #              step = 10^(-5)),
      numericInput("gamma", "Predation rate coefficient gamma",
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
                   step = 10^(-2)),
      numericInput("l50", "L50",
                   value = sp$l50,
                   min = 0,
                   max = 100, 
                   step = 1),
      numericInput("l25", "L25",
                   value = sp$l25,
                   min = 0,
                   max = 100, 
                   step = 1)
    )
  })
  
  ## Handle species parameter change ####
  # observe({
  #   req(input$gamma, input$h)
  #   p <- isolate(params())
  #   sp <- isolate(input$sp_sel)
  # The version where a change in parameter automatically triggers this
  # observer does not work yet because it gets triggered also by the
  # rewriting of the input controls upon change of target species.
  # So for now require "Go" button.
  observeEvent(input$sp_go, {
    p <- params()
    sp <- input$sp_sel

    # Create updated species params data frame
    species_params <- p@species_params
    species_params[sp, "gamma"] <- input$gamma
    species_params[sp, "h"]     <- input$h
    species_params[sp, "alpha"] <- input$alpha
    species_params[sp, "ks"]    <- input$ks
    species_params[sp, "beta"]  <- input$beta
    species_params[sp, "sigma"] <- input$sigma
    species_params[sp, "k_vb"]  <- input$k_vb
    species_params[sp, "a"]     <- input$a
    species_params[sp, "b"]     <- input$b
    species_params[sp, "l50"]   <- input$l50
    species_params[sp, "l25"]   <- input$l25

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

    # The spectrum for the changed species is calculated with new
    # parameters but in the context of the original community
    # Compute death rate for changed species
    mumu <- getMort(pc, p@initial_n, p@initial_n_pp, effort = effort)[sp, ]
    # compute growth rate for changed species
    gg <- getEGrowth(pc, p@initial_n, p@initial_n_pp)[sp, ]
    # Compute solution for changed species
    w_inf_idx <- sum(pc@w < pc@species_params[sp, "w_inf"])
    idx <- p@w_min_idx[sp]:(w_inf_idx - 1)
    validate(
      need(!any(gg[idx] == 0),
           "Can not compute steady state due to zero growth rates")
    )
    n0 <- p@initial_n[sp, p@w_min_idx[sp]]
    pc@initial_n[sp, ] <- 0
    pc@initial_n[sp, pc@w_min_idx[sp]:w_inf_idx] <-
      c(1, cumprod(gg[idx] / ((gg + mumu * pc@dw)[idx + 1]))) *
      n0
    if (any(is.infinite(pc@initial_n))) {
      stop("Candidate steady state holds infinities")
    }
    if (any(is.na(pc@initial_n) || is.nan(pc@initial_n))) {
      stop("Candidate steady state holds none numeric values")
    }

    # Update the reactive params object
    params(pc)
  })

  
  ## Recompute all species ####
  # triggered by "Multi" button on "Species" tab
  observeEvent(input$sp_multi, {
    p <- params()
    
    # Recompute plankton
    plankton_mort <- getPlanktonMort(p, p@initial_n, p@initial_n_pp)
    p@initial_n_pp <- p@rr_pp * p@cc_pp / (p@rr_pp + plankton_mort)
    # Recompute all species
    mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = effort)
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    for (sp in 1:length(p@species_params$species)) {
      w_inf_idx <- sum(p@w < p@species_params[sp, "w_inf"])
      idx <- p@w_min_idx[sp]:(w_inf_idx - 1)
      validate(
        need(!any(gg[sp, idx] == 0),
             "Can not compute steady state due to zero growth rates")
      )
      n0 <- p@initial_n[sp, p@w_min_idx[sp]]
      p@initial_n[sp, ] <- 0
      p@initial_n[sp, p@w_min_idx[sp]:w_inf_idx] <- 
        c(1, cumprod(gg[sp, idx] / ((gg[sp, ] + mumu[sp, ] * p@dw)[idx + 1]))) *
        n0
    }    
    # Update the reactive params object
    params(p)
  })
    
    
  ## Find new steady state ####
  # triggered by "Steady" button on "species" tab
  observeEvent(input$sp_steady, {
    p <- params()
    
    # Create a Progress object
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    
    # Run to steady state
    p <- steady(p, effort = 1.4, t_max = 100, tol = 1e-2,
                 shiny_progress = progress)
    
    # Update the reactive params object
    params(p)
  })
  
  ## Reconstruct params object ####
  # This is triggered by the "Go" button on the "General" tab
  observeEvent(input$bg_go, {
    
    # Create a Progress object
    progress <- shiny::Progress$new(session)
    on.exit(progress$close())
    
    p_old <- params()
    
    p <- setBackground(
      set_scaling_model(
        #min_w_pp = input$min_w_pp, no_sp = input$no_bg_sp, no_w = 400,
        min_w_pp = 1e-12, no_sp = input$no_bg_sp, no_w = input$no_w,
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
      if (input$use_SSB) {
        p <- addSpecies(p, p_old@species_params[sp, ],
                        effort = effort,
                        rfac = Inf, SSB = input$SSB)
      } else {
        p <- addSpecies(p, p_old@species_params[sp, ],
                        effort = effort,
                        rfac = Inf)    
      }
    }
    
    # Run to steady state
    p <- steady(p, effort = c(knife_edge_gear = 0, sigmoid_gear = input$effort, 
                              sigmoid_gear_Anchovy = input$Anchovy_effort), 
                t_max = 100, tol = 1e-2,
                shiny_progress = progress)
    # Update the reactive params object
    params(p)
  })
  
  ## Create plots ####
  output$plotGrowthCurve <- renderPlot({
    plotGrowthCurves(params(), species = input$sp_sel)
  })
  
  output$plotSpectra <- renderPlot({
    plotSpectra(params())
  })
  
  output$plot_erepro <- renderPlot({
    p <- params()
    # Retune the values of erepro so that we get the correct level of
    # recruitment
    mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = effort)
    gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
    rdd <- getRDD(p, p@initial_n, p@initial_n_pp)
    # TODO: vectorise this
    for (i in (1:length(p@species_params$species))) {
      gg0 <- gg[i, p@w_min_idx[i]]
      mumu0 <- mumu[i, p@w_min_idx[i]]
      DW <- p@dw[p@w_min_idx[i]]
      p@species_params$erepro[i] <- p@species_params$erepro[i] *
        (p@initial_n[i, p@w_min_idx[i]] *
           (gg0 + DW * mumu0)) / rdd[i]
    }
    ggplot(p@species_params, aes(x = species, y = erepro)) + 
      geom_col() + geom_hline(yintercept = 1, color = "red")
  })
  
} #the server

#### User interface ####
ui <- fluidPage(
  
  titlePanel("Humboldt current ecosystem"),
  
  sidebarLayout(
    
    ## Sidebar ####
    sidebarPanel(
      tabsetPanel(
        tabPanel("Species",
                 uiOutput("sp_sel"),
                 actionButton("sp_go", "Go"),
                 actionButton("sp_multi", "Multi"),
                 actionButton("sp_steady", "Steady"),
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
                 sliderInput("no_w", "Number of weight brackets",
                             value=400, min=200, max=1200, step=50, round = TRUE),
                 numericInput("min_w_pp", "Minimum plankton weight min_w_pp",
                              value=1e-12,  step=1e-13),
                 checkboxInput("use_SSB", "Use target SSB", value = FALSE)
        ),
        tabPanel("File",
                 downloadButton("params", "Download current params object"),
                 fileInput("upload", "Upload new params object", accept = ".rds")
        )
      )
    ),  # endsidebarpanel
    
    ## Main panel ####
    mainPanel(
      plotOutput("plot_erepro", height = "150px"),
      tabsetPanel(type = "tabs",
                  tabPanel("Spectra", plotOutput("plotSpectra")),
                  tabPanel("Growth", plotOutput("plotGrowthCurve"))
      )
    )  # end mainpanel
  )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app


