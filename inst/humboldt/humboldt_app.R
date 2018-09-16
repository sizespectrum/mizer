library(shiny)
library(ggplot2)
library(plotly)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="humboldt")
# library(mizer)
library(progress)

server <- function(input, output, session) {
  effort <- 1.4
  
  ## Load params object and store it as a reactive value ####
  params <- reactiveVal()
  filename <- "humboldt_params.rds"
  params(readRDS(filename))
  output$filename <- renderText(paste0("Previously uploaded file: ", filename))
  
  ## Handle upload of params object ####
  observeEvent(input$upload, {
    inFile <- input$upload
    output$filename <- renderText(paste0("Previously uploaded file: ", 
                                         inFile$name))
    tryCatch(
      p <- readRDS(inFile$datapath),
      error = function(e) {stop(safeError(e))}
    )
    # Update the reactive params object
    params(p)
  })
  
  ## Prepare for download of params object ####
  output$params <- downloadHandler(
    filename = "params.rds", 
    content = function(file) {
      saveRDS(params(), file = file)
    })
  
  ## UI for species parameters ####
  output$sp_sel <- renderUI({
    p <- isolate(params())
    species <- as.character(p@species_params$species[!is.na(p@A)])
    selectInput("sp_sel", "Species:", species) 
  })
  output$params_sliders <- renderUI({
    # The parameter sliders get updated whenever the species selector changes
    req(input$sp_sel)
    # We do not want the updating of the slider triger any of the other
    # actions that trigger when a parameter value is changed, so we freeze
    # one of those inputs
    freezeReactiveValue(input, "h")
    
    p <- isolate(params())
    sp <- p@species_params[input$sp_sel, ]
    n0 <- p@initial_n[input$sp_sel, p@w_min_idx[input$sp_sel]]
    list(
      sliderInput("n0", "Egg density",
                  value = n0,
                  min = signif(n0/10, 3),
                  max = signif(n0*10, 3)),
      sliderInput("gamma", "Predation rate coefficient gamma",
                  value = sp$gamma,
                  min = signif(sp$gamma/2, 3),
                  max = signif(sp$gamma*2, 3)),
      sliderInput("h", "max feeding rate h",
                  value = sp$h,
                  min = signif(sp$h/2, 2),
                  max = signif(sp$h*2, 2)),
      sliderInput("alpha", "Assimilation efficiency alpha",
                  value = sp$alpha,
                  min = 0,
                  max = 1),
      sliderInput("ks", "Coefficient of standard metabolism ks",
                  value = sp$ks,
                  min = signif(sp$ks/2, 2),
                  max = signif(sp$ks*2, 2),
                  step = 0.05),
      sliderInput("beta", "Preferred predator-prey mass ratio beta",
                  value = round(sp$beta),
                  min = round(sp$beta/10),
                  max = round(sp$beta*10),
                  step = 1),
      sliderInput("sigma", "Width of size selection function sigma",
                  value = sp$sigma,
                  min = round(sp$sigma/2),
                  max = round(sp$sigma*2),
                  step = 0.05),
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
  
  ## UI for general parameters ####

  output$general_params <- renderUI({
    p <- isolate(params())
    i_bkgd <- which.max(is.na(p@A))
    bkgd_params <- p@species_params[i_bkgd, ]
    
    list(
      sliderInput("kappa", "kappa", value = p@kappa,
                  min = p@kappa / 2,
                  max = p@kappa * 1.5,
                  width = "80%"),
      numericInput("lambda", "Sheldon exponent",
                   value = p@lambda, min = 1.9, max = 2.2, step = 0.005),
      sliderInput("f0", "Feeding level",
                  value = p@f0, min = 0, max = 1),
      sliderInput("h_bkgd", "max feeding rate",
                  value = bkgd_params$h, min = 10, max = 100, step = 2),
      sliderInput("log_r_pp", "log10 Plankton replenishment rate",
                  value = -1, min = -4, max = 0),
      sliderInput("no_bg_sp", "Number of background species",
                  value = 10, min = 4, max = 20, step = 1, round = TRUE),
      sliderInput("no_w", "Number of weight brackets",
                  value = 400, min = 200, max = 1200, step = 50, round = TRUE),
      numericInput("min_w_pp", "Minimum plankton weight min_w_pp",
                   value = 1e-12,  step = 1e-13)
    )
  })
  
  ## Adjust kappa ####
  observe({
    req(input$kappa)
    p <- isolate(params())
    # We want a change in kappa to rescale all abundances by the same factor
    p@initial_n <- p@initial_n * input$kappa / p@kappa
    p@initial_n_pp <- p@initial_n_pp * input$kappa / p@kappa
    p@cc_pp <- p@cc_pp * input$kappa / p@kappa
    # To keep the same per-capity behaviour, we have to scale down the
    # search volume
    p@species_params$gamma <- p@species_params$gamma / (input$kappa / p@kappa)
    p@search_vol <- p@search_vol / (input$kappa / p@kappa)
    p@kappa <- input$kappa
    updateSliderInput(session, "kappa",
                      min = input$kappa / 2, 
                      max = input$kappa * 1.5)
    params(p)
  })
  
  ## Adjust k_vb ####
  observe({
    req(input$k_vb)
    p <- isolate(params())
    p@species_params[isolate(input$sp_sel), "k_vb"] <- input$k_vb
    params(p)
  })
  
  ## Handle species parameter change ####
  # This is triggered when any of the species inputs changes
  observe({
    # I do not want this to run at the start of the app, but don't know how
    # to avoid that. But at least I can make sure it does not run before
    # the last input value has been given its initial value.
    req(input$l25)
    
    p <- isolate(params())
    sp <- isolate(input$sp_sel)
    
    # wrap the code in trycatch so that when there is a problem we can
    # simply stay with the old parameters
    tryCatch({
      # rescale abundance to new egg density
      p@initial_n[sp, ] <- p@initial_n[sp, ] * input$n0 / 
        p@initial_n[sp, p@w_min_idx[sp]]
      
      # Create updated species params data frame
      species_params <- p@species_params
      species_params[sp, "gamma"] <- input$gamma
      species_params[sp, "h"]     <- input$h
      species_params[sp, "alpha"] <- input$alpha
      species_params[sp, "ks"]    <- input$ks
      species_params[sp, "beta"]  <- input$beta
      species_params[sp, "sigma"] <- input$sigma
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
      
      # Retune the value of erepro so that we get the correct level of
      # recruitment
      i <- which(pc@species_params$species == sp)
      rdd <- getRDD(pc, pc@initial_n, pc@initial_n_pp)[i]
      gg0 <- gg[pc@w_min_idx[i]]
      mumu0 <- mumu[pc@w_min_idx[i]]
      DW <- pc@dw[pc@w_min_idx[i]]
      pc@species_params$erepro[i] <- pc@species_params$erepro[i] *
        n0 * (gg0 + DW * mumu0) / rdd
      
      # Update slider min/max so that they are a fixed proportion of the 
      # parameter value
      updateSliderInput(session, "n0",
                        min = signif(input$n0 / 10, 3),
                        max = signif(input$n0 * 10, 3))
      updateSliderInput(session, "gamma",
                        min = signif(input$gamma/2, 3),
                        max = signif(input$gamma*2, 3))
      updateSliderInput(session, "h",
                        min = signif(input$h/2, 2),
                        max = signif(input$h*2, 2))
      updateSliderInput(session, "ks",
                        min = signif(input$ks/2, 2),
                        max = signif(input$ks*2, 2))
      updateSliderInput(session, "beta",
                        min = round(input$beta/10),
                        max = round(input$beta*10))
      updateSliderInput(session, "sigma",
                        min = round(input$sigma/2),
                        max = round(input$sigma*2))
      
      if (input$log_sp) {
        # Save new species params to disk
        time = format(Sys.time(), "_%Y_%m_%d_at_%H_%M_%S")
        file = paste0("species_params", time, ".rds")
        saveRDS(pc@species_params, file = file)
      }
      
      # Update the reactive params object
      params(pc)
    }, 
    error = function(e) {
      showModal(modalDialog(
        title = "Invalid parameters",
        HTML(paste0("These parameter values do not lead to an acceptable steady state. ",
                    "I will keep the previous values.<br>",
                    "The error message was:<br>", e)),
        easyClose = TRUE
      ))
      params(p)}
    )
  })

  
  ## Recompute all species ####
  # triggered by "Interact" button on "Species" tab
  observeEvent(input$sp_interact, {
    p <- params()
    
    tryCatch({
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
          p@initial_n[i, p@w_min_idx[i]] *
          (gg0 + DW * mumu0) / rdd[i]
      }
      
      # Update the reactive params object
      params(p)
    },
    error = function(e) {
      showModal(modalDialog(
        title = "Invalid parameters",
        HTML(paste0("These parameter do not lead to an acceptable steady state.",
                    "Please choose other values.<br>",
                    "The error message was:<br>", e)),
        easyClose = TRUE
      ))}
    )
  })
    
    
  ## Find new steady state ####
  # triggered by "Steady" button on "species" tab
  observeEvent(input$sp_steady, {
    p <- params()
    
    tryCatch({
      # Create a Progress object
      progress <- shiny::Progress$new(session)
      on.exit(progress$close())
      
      # Run to steady state
      p <- steady(p, effort = 1.4, t_max = 100, tol = 1e-2,
                  shiny_progress = progress)
      
      if (input$log_steady) {
        # Save new params object to disk
        time = format(Sys.time(), "_%Y_%m_%d_at_%H_%M_%S")
        file = paste0("mizer_params", time, ".rds")
        saveRDS(p, file = file)
      }
      
      # Update the reactive params object
      params(p)
    },
    error = function(e) {
      showModal(modalDialog(
        title = "Invalid parameters",
        HTML(paste0("These parameter do not lead to an acceptable steady state.",
                    "Please choose other values.<br>",
                    "The error message was:<br>", e)),
        easyClose = TRUE
      ))}
    )
  })
  
  ## Reconstruct params object ####
  # This is triggered by the "Go" button on the "General" tab
  observeEvent(input$bg_go, {
    
    tryCatch({
      # Create a Progress object
      progress <- shiny::Progress$new(session)
      on.exit(progress$close())
      
      p_old <- params()
      
      p <- setBackground(
        set_scaling_model(
          min_w_pp = input$min_w_pp, 
          no_sp = input$no_bg_sp, no_w = input$no_w,
          min_w_inf = 2, max_w_inf = 6e6,
          min_egg = 1e-4, min_w_mat = 2 / 10^0.8,
          lambda = input$lambda, knife_edge_size = Inf,
          beta = 500, sigma = 2,
          f0 = input$f0, h = input$h_bkgd, r_pp = 10^input$log_r_pp
        )
      )
      
      # Loop over all foreground species and add them one-by-one to the new
      # background
      gears <- "knife_edge"
      no_sp <- length(p_old@A)
      for (sp in (1:no_sp)[!is.na(p_old@A)]) {
        if (p_old@species_params$species[sp] == "Anchovy") {
          gear <- "sigmoid_gear_Anchovy"
        } else {
          gear <- "sigmoid_gear"
        }
        gears <- union(gears, gear)
        p <- addSpecies(p, p_old@species_params[sp, ],
                        effort = effort,
                        rfac = Inf)    
      }
      
      # Run to steady state
      p <- steady(p, effort = effort, 
                  t_max = 100, tol = 1e-2,
                  shiny_progress = progress)
      
      # Update the reactive params object
      params(p)
    },
    error = function(e) {
      showModal(modalDialog(
        title = "Invalid parameters",
        HTML(paste0("These parameter do not lead to an acceptable steady state.",
                    "Please choose other values.<br>",
                    "The error message was:<br>", e)),
        easyClose = TRUE
      ))}
    )
  })
  
  ## Growth curves ####
  output$k_vb_sel <- renderUI({
    req(input$sp_sel)
    k_vb <- params()@species_params[input$sp_sel, "k_vb"]
    numericInput("k_vb", "Von Bertalanffy k", value = k_vb)
  })
  output$plotGrowthCurve <- renderPlot({
    plotGrowthCurves(params(), species = input$sp_sel)
  })
  
  ## Spectra ####
  output$plotSpectra <- renderPlot({
    plotSpectra(params())
  })
  
  ## erepro plot ####
  output$plot_erepro <- renderPlot({
    p <- params()
    ggplot(p@species_params, aes(x = species, y = erepro)) + 
      geom_col() + geom_hline(yintercept = 1, color = "red") +
      scale_y_log10()
  })
  
  ## Biomass plot ####
  output$biomass_sel <- renderUI({
    sp <- input$sp_sel
    species_params <- params()@species_params[sp, ]
    list(
      div(style = "display:inline-block",
          numericInput("biomass_observed", 
                       paste0("Observed biomass for ", sp, 
                              "(megatonnes)"),
                       value = species_params$biomass_observed)),
      div(style = "display:inline-block",
          numericInput("biomass_cutoff", "Lower cutoff (grams)",
                       value = species_params$biomass_cutoff))
    )
  })
  output$plotBiomass <- renderPlot({
    req(input$sp_sel, input$biomass_cutoff, input$biomass_observed)
    sp <- input$sp_sel
    p <- params()
    biomass <- cumsum(p@initial_n[sp, ] * p@w * p@dw)
    
    cutoff_idx <- which.max(p@w >= input$biomass_cutoff)
    target <- input$biomass_observed + biomass[cutoff_idx]
    
    max_w <- p@species_params[sp, "w_inf"]
    min_w <- p@species_params[sp, "w_min"]
    sel <- p@w >= min_w & p@w <= max_w
    df <- data.frame(Size = p@w[sel], Biomass = biomass[sel])
    ggplot(df, aes(x = Size, y = Biomass)) + 
      geom_line(color = "blue") + scale_x_log10() +
      geom_hline(yintercept = biomass[cutoff_idx]) +
      geom_vline(xintercept = input$biomass_cutoff) +
      geom_hline(yintercept = target, color = "green")
  })
  output$plotObservedBiomass <- renderPlot({
    p <- params()
    no_sp <- length(p@species_params$species)
    cutoff <- p@species_params$biomass_cutoff
    observed <- p@species_params$biomass_observed
    # When no cutoff known, set it to maturity weight / 20
    cutoff[is.na(cutoff)] <- p@species_params$w_mat[is.na(cutoff)] / 20
    
    biomass_model <- 1:no_sp  # create vector of right length
    for (sp in 1:no_sp) {
      cum_biomass <- cumsum(p@initial_n[sp, ] * p@w * p@dw)
      cutoff_idx <- which.max(p@w >= cutoff[sp])
      biomass_model[sp] <- max(cum_biomass) - cum_biomass[cutoff_idx]
    }
    df <- rbind(
      data.frame(Species = p@species_params$species,
                 Type = "Observed",
                 Biomass = observed),
      data.frame(Species = p@species_params$species,
                 Type = "Model",
                 Biomass = biomass_model)
    )
    ggplot(df) +
      geom_col(aes(x = Species, y = Biomass, fill = Type),
               position = "dodge")
  })
  
  ## Plot catch ####
  
  # Catch by size for selected species
  output$plotCatch <- renderPlotly({
    req(input$sp_sel)
    p <- params()
    sp <- which.max(p@species_params$species == input$sp_sel)
    w_min_idx <- sum(p@w < (p@species_params$w_mat[sp] / 100))
    w_max_idx <- sum(p@w <= p@species_params$w_inf[sp])
    w_sel <- seq(w_min_idx, w_max_idx, by = 1)
    catch <- getFMort(p, effort = effort)[sp, w_sel] *
      p@w[w_sel] * p@dw[w_sel] * p@initial_n[sp, w_sel]
    df <- data.frame(Size = p@w[w_sel], Catch = catch)
    ggplot(df) +
      geom_line(aes(x = Size, y = Catch)) +
      scale_x_log10() +
      geom_vline(xintercept = p@species_params$w_mat[sp])
  })
  
  # Total catch by species
  output$plotObservedCatch <- renderPlot({
    p <- params()
    biomass <- sweep(p@initial_n, 2, p@w * p@dw, "*")
    catch <- rowSums(biomass * getFMort(p, effort = effort))
    df <- rbind(
      data.frame(Species = p@species_params$species,
                 Type = "Observed",
                 Catch = p@species_params$catch_observed),
      data.frame(Species = p@species_params$species,
                 Type = "Model",
                 Catch = catch)
    )
    ggplot(df) +
      geom_col(aes(x = Species, y = Catch, fill = Type),
               position = "dodge")
  })
  
  # Input field for observed catch
  output$catch_sel <- renderUI({
    sp <- input$sp_sel
    numericInput("catch_observed", 
                 paste0("Observed total catch for ", sp, "(megatonnes)"),
                 value = params()@species_params[sp, "catch_observed"])
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
                 tags$br(),
                 actionButton("sp_interact", "Interact"),
                 actionButton("sp_steady", "Steady"),
                 uiOutput("sp_sel"),
                 uiOutput("params_sliders")
        ),
        tabPanel("General",
                 tags$br(),
                 actionButton("bg_go", "Go"),
                 uiOutput("general_params")
        ),
        tabPanel("File",
                 tags$br(),
                 downloadButton("params", "Download current params object"),
                 checkboxInput("log_steady", "Log steady states",
                               value = FALSE),
                 checkboxInput("log_sp", "Log species parameters",
                               value = FALSE),
                 tags$hr(),
                 textOutput("filename"),
                 fileInput("upload", "Upload new params object", 
                           accept = ".rds")
        )
      ),
      width = 3
    ),  # endsidebarpanel
    
    ## Main panel ####
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel("Spectra", plotOutput("plotSpectra")),
        tabPanel("Biomass",
                 plotOutput("plotObservedBiomass"),
                 plotOutput("plotBiomass"),
                 uiOutput("biomass_sel")),
        tabPanel("Growth",
                 plotOutput("plotGrowthCurve"),
                 uiOutput("k_vb_sel")),
        tabPanel("Repro",
                 plotOutput("plot_erepro")),
        tabPanel("Catch",
                 plotOutput("plotObservedCatch"),
                 uiOutput("catch_sel"),
                 plotlyOutput("plotCatch"))
      )
    )  # end mainpanel
  )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app


