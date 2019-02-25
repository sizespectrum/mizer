library(shiny)
library(ggplot2)
library(plotly)
library(reshape2)
library(tidyr)
# # Uncomment the following 2 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="humboldt")
library(mizer)
library(progress)

server <- function(input, output, session) {
  
  ## Load params object and store it as a reactive value ####
  params <- reactiveVal()
  filename <- "humboldt_params.rds"
  params(readRDS(filename))
  output$filename <- renderText(paste0("Previously uploaded file: ", filename))
  
  # Define some globals to skip certain observers
  skip_update <- TRUE
  skip_update_n0 <- TRUE
  # Define a reactive value for triggering an update of species sliders
  trigger_update <- reactiveVal(0)
  
  # Load catch distribution
  catchdist <- readRDS("catchdistribution.rds")
  
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
    
    # Trigger an update of sliders
    trigger_update(runif(1))
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
    selectInput("sp", "Species:", species) 
  })
  output$sp_params <- renderUI({
    # The parameter sliders get updated whenever the species selector changes
    req(input$sp)
    # or when the trigger is set somewhere
    trigger_update()
    
    p <- isolate(params())
    sp <- p@species_params[input$sp, ]
    n0 <- p@initial_n[input$sp, p@w_min_idx[input$sp]]
    
    # We do not want the updating of the slider to triger an update of the
    # params object
    skip_update <<- TRUE
    skip_update_n0 <<- TRUE
    
    list(
      tags$h3("Predation"),
      sliderInput("gamma", "Predation rate coefficient gamma",
                  value = sp$gamma,
                  min = signif(sp$gamma / 2, 3),
                  max = signif(sp$gamma * 1.5, 3)),
      sliderInput("h", "max feeding rate h",
                  value = sp$h,
                  min = signif(sp$h / 2, 2),
                  max = signif(sp$h * 1.5, 2)),
      sliderInput("beta", "Preferred predator-prey mass ratio beta",
                  value = sp$beta,
                  min = signif(sp$beta / 2, 2),
                  max = signif(sp$beta * 1.5, 2)),
      sliderInput("sigma", "Width of size selection function sigma",
                  value = sp$sigma,
                  min = signif(sp$sigma / 2, 2),
                  max = signif(sp$sigma * 1.5, 2),
                  step = 0.05),
      tags$h3("Fishing"),
      sliderInput("catchability", "Catchability",
                   value = sp$catchability, min = 0, max = 1),
      sliderInput("l50", "L50",
                   value = sp$l50, 
                   min = 1, 
                   max = signif(sp$l50 * 2, 2),
                  step = 0.1),
      sliderInput("ldiff", "L50-L25",
                   value = sp$l50 - sp$l25, 
                   min = 0.1, 
                   max = signif(sp$l50 / 10, 2),
                  step = 0.1),
      numericInput("a", "Coefficient for length to weight conversion a",
                   value = sp$a),
      numericInput("b", "Exponent for length to weight conversion b",
                   value = sp$b),
      tags$h3("Maturity"),
      sliderInput("w_mat", "w_mat", value = sp$w_mat,
                  min = signif(sp$w_mat / 2, 2),
                  max = signif(sp$w_mat * 1.5, 2)),
      sliderInput("wfrac", "w25/w_mat", value = sp$w25/sp$w_mat,
                  min = 0.5,
                  max = 1,
                  step = 0.01),
      sliderInput("m", "m", value = sp$m,
                  min = 0,
                  max = 1),
      tags$h3("Others"),
      sliderInput("kappa", "kappa", value = p@kappa,
                  min = signif(p@kappa / 2, 2),
                  max = signif(p@kappa * 1.5, 2)),
      sliderInput("n0", "Egg density",
                  value = n0,
                  min = signif(n0 / 2, 3),
                  max = signif(n0 * 1.5, 3)),
      sliderInput("alpha", "Assimilation efficiency alpha",
                  value = sp$alpha,
                  min = 0,
                  max = 1),
      sliderInput("ks", "Coefficient of standard metabolism ks",
                  value = sp$ks,
                  min = signif(sp$ks / 2, 2),
                  max = signif(sp$ks * 1.5, 2),
                  step = 0.05),
      sliderInput("k", "Coefficient of activity k",
                  value = sp$k,
                  min = 0,
                  max = 1,
                  step = 0.01)
    )
  })
  
  ## UI for general parameters ####

  output$general_params <- renderUI({
    p <- isolate(params())
    
    list(
      numericInput("lambda", "Sheldon exponent",
                   value = p@lambda, min = 1.9, max = 2.2, step = 0.005),
      sliderInput("log_r_pp", "log10 Plankton replenishment rate",
                  value = -1, min = -4, max = 0),
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
                      min = signif(input$kappa / 2, 2),
                      max = signif(input$kappa * 1.5, 2))
    params(p)
    
    # Trigger an update of sliders
    trigger_update(runif(1))
  })
  
  ## Adjust k_vb ####
  observe({
    p <- isolate(params())
    p@species_params[isolate(input$sp), "k_vb"] <- req(input$k_vb)
    p@species_params[isolate(input$sp), "t0"] <- req(input$t0)
    params(p)
  })
  
  ## Adjust biomass observed ####
  observe({
    p <- isolate(params())
    p@species_params[isolate(input$sp), "biomass_observed"] <- 
      req(input$biomass_observed)
    p@species_params[isolate(input$sp), "biomass_cutoff"] <- 
      req(input$biomass_cutoff)
    params(p)
  })  
  
  ## Adjust catch observed ####
  observe({
    p <- isolate(params())
    p@species_params[isolate(input$sp), "catch_observed"] <- 
      req(input$catch_observed)
    params(p)
  })
  
  # Adjust egg density ####
  observe({
    n0 <- req(input$n0)
    p <- isolate(params())
    sp <- isolate(input$sp)
    
    if (skip_update_n0) {
      skip_update_n0 <<- FALSE
    } else {
      updateSliderInput(session, "n0",
                        min = signif(n0 / 2, 3),
                        max = signif(n0 * 1.5, 3))
      # rescale abundance to new egg density
      p@initial_n[sp, ] <- p@initial_n[sp, ] * n0 / 
        p@initial_n[sp, p@w_min_idx[sp]]
      
      # Update the reactive params object
      params(p)
    }
  })
  
  ## Adjust species parameters ####
  update_species <- function(sp, p, species_params) {
    
    # wrap the code in trycatch so that when there is a problem we can
    # simply stay with the old parameters
    tryCatch({
      
      # Create updated species params data frame
      
      
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
      mumu <- getMort(pc, p@initial_n, p@initial_n_pp, 
                      effort = 1)[sp, ]
      # compute growth rate for changed species
      gg <- getEGrowth(pc, p@initial_n, p@initial_n_pp)[sp, ]
      # Compute solution for changed species
      w_inf_idx <- sum(pc@w < pc@species_params[sp, "w_inf"])
      idx <- p@w_min_idx[sp]:(w_inf_idx - 1)
      if (any(gg[idx] == 0)) {
        weight <- p@w[which.max(gg[idx] == 0)]
        showModal(modalDialog(
          title = "Zero growth rate",
          paste0("With these parameter values the ", sp,
                 " does not have enough food to cover its metabolic cost"),
          easyClose = TRUE
        ))
      }
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
        HTML(paste0("These parameter values lead to an error.<br>",
                    "The error message was:<br>", e)),
        easyClose = TRUE
      ))
      params(p)}
    )
  }
  
  observe({
    req(input$sigma)
    p <- isolate(params())
    sp <- isolate(input$sp)
    
    species_params <- p@species_params
    species_params[sp, "gamma"] <- input$gamma
    species_params[sp, "h"]     <- input$h
    species_params[sp, "beta"]  <- input$beta
    species_params[sp, "sigma"] <- input$sigma
    species_params[sp, "catchability"]   <- input$catchability
    species_params[sp, "a"]     <- input$a
    species_params[sp, "b"]     <- input$b
    species_params[sp, "l50"]   <- input$l50
    species_params[sp, "l25"]   <- input$l50 - input$ldiff
    species_params[sp, "alpha"] <- input$alpha
    species_params[sp, "ks"]    <- input$ks
    species_params[sp, "k"]     <- input$k
    species_params[sp, "w25"]   <- input$w_mat * input$wfrac
    species_params[sp, "w_mat"]   <- input$w_mat
    species_params[sp, "m"]     <- input$m
    
    if (skip_update) {
      skip_update <<- FALSE
    } else {
      # Update slider min/max so that they are a fixed proportion of the 
      # parameter value
      updateSliderInput(session, "gamma",
                        min = signif(input$gamma / 2, 3),
                        max = signif(input$gamma * 1.5, 3))
      updateSliderInput(session, "h",
                        min = signif(input$h / 2, 2),
                        max = signif(input$h * 1.5, 2))
      updateSliderInput(session, "beta",
                        min = signif(input$beta / 2, 2),
                        max = signif(input$beta * 1.5, 2))
      updateSliderInput(session, "sigma",
                        min = signif(input$sigma / 2, 2),
                        max = signif(input$sigma * 1.5, 2))
      updateSliderInput(session, "l50",
                        max = signif(input$l50 * 2, 2))
      updateSliderInput(session, "ldiff",  
                        max = signif(input$l50 / 10, 2))
      updateSliderInput(session, "ks",
                        min = signif(input$ks / 2, 2),
                        max = signif(input$ks * 1.5, 2))
      updateSliderInput(session, "w_mat",
                        min = signif(input$w_mat / 2, 2),
                        max = signif(input$w_mat * 1.5, 2))
      
      update_species(sp, p, species_params)
    }
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
      mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = 1)
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
      mumu <- getMort(p, p@initial_n, p@initial_n_pp, effort = 1)
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
      p <- steady(p, effort = 1, t_max = 100, tol = 1e-2,
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
                        effort = 1, rfac = Inf)    
      }
      
      # Run to steady state
      p <- steady(p, effort = 1, 
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
    req(input$sp)
    k_vb <- params()@species_params[input$sp, "k_vb"]
    t0 <- params()@species_params[input$sp, "t0"]
    list(
      div(style = "display:inline-block",
          numericInput("k_vb", "Von Bertalanffy k", value = k_vb)),
      div(style = "display:inline-block",
          numericInput("t0", "t_0", value = t0))
    )
  })
  output$plotGrowthCurve <- renderPlot({
    plotGrowthCurves(params(), species = input$sp) +
      theme_grey(base_size = 18)
  })
  
  ## Spectra ####
  output$plotSpectra <- renderPlot({
    if (input$binning == "Logarithmic") {
      power <- 2
    } else {
      power <- 1
    }
    plotSpectra(params(), power = power) + theme_grey(base_size = 18)
  })
  
  ## erepro plot ####
  output$plot_erepro <- renderPlot({
    p <- params()
    ggplot(p@species_params, aes(x = species, y = erepro)) + 
      geom_col() + geom_hline(yintercept = 1, color = "red") +
      scale_y_log10() +
      theme_grey(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  })
  
  ## Plot feeding level ####
  output$plot_feeding_level <- renderPlotly({
    p <- params()
    fl <- getFeedingLevel(p, p@initial_n, p@initial_n_pp)
    df <- melt(fl)
    ggplot(df) +
      geom_line(aes(x = w, y = value, color = sp, linetype = sp)) +
      scale_x_log10()
  })
  
  ## Biomass plot ####
  output$biomass_sel <- renderUI({
    sp <- input$sp
    p <- isolate(params())
    species_params <- p@species_params[sp, ]
    if (is.na(species_params$biomass_observed)) {
      species_params$biomass_observed <- 0
    }
    if (is.na(species_params$biomass_cutoff)) {
      species_params$biomass_cutoff <- 0
    }
    list(
      div(style = "display:inline-block",
          numericInput("biomass_observed", 
                       paste0("Observed biomass for ", sp, " (megatonnes)"),
                       value = species_params$biomass_observed)),
      div(style = "display:inline-block",
          numericInput("biomass_cutoff", "Lower cutoff (grams)",
                       value = species_params$biomass_cutoff))
    )
  })
  output$plotBiomassDist <- renderPlot({
    req(input$sp, input$biomass_cutoff, input$biomass_observed)
    sp <- input$sp
    p <- params()
    biomass <- cumsum(p@initial_n[sp, ] * p@w * p@dw)
    
    max_w <- p@species_params[sp, "w_inf"]
    min_w <- p@species_params[sp, "w_min"]
    sel <- p@w >= min_w & p@w <= max_w
    df <- data.frame(Size = p@w[sel], Biomass = biomass[sel])
    pl <- ggplot(df, aes(x = Size, y = Biomass)) + 
      geom_line(color = "blue") + scale_x_log10() +
      geom_vline(xintercept = p@species_params[sp, "w_mat"], 
                 linetype = "dotted") +
      theme_grey(base_size = 18) +
      labs(x = "Size [g]", y = "Cummulative biomass [megatonnes]")  +
      geom_text(aes(x = p@species_params[sp, "w_mat"], 
                    y = max(Biomass * 0.2),
                    label = "\nMaturity"), 
                angle = 90)
    if (input$biomass_observed) {
      cutoff_idx <- which.max(p@w >= input$biomass_cutoff)
      target <- input$biomass_observed + biomass[cutoff_idx]
      pl <- pl +
        geom_hline(yintercept = biomass[cutoff_idx]) +
        geom_vline(xintercept = input$biomass_cutoff) +
        geom_hline(yintercept = target, color = "green")
    }
    pl
  })
  output$plotTotalBiomass <- renderPlot({
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
               position = "dodge") +
      theme_grey(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(x = "", y = "Biomass [megatonnes]")
  })
  
  ## Plot catch ####
  
  # Catch density for selected species
  output$plotCatchDist <- renderPlot({
    req(input$sp)
    p <- params()
    sp <- which.max(p@species_params$species == input$sp)
    a <- p@species_params$a[sp]
    b <- p@species_params$b[sp]
    
    # Check whether we have enough catch data for this species to plot it
    is_observed <- sum(catchdist$species == input$sp) > 3
    
    # To choose the range of sizes over which to plot we look at the range
    # of sizes for which a non-zero catch was observed. If no catch was
    # observed for the species, we use the range from w_mat/100 to w_inf.
    if (is_observed) {
      l_min = min(catchdist$length[catchdist$species == input$sp])
      w_min = a * l_min ^ b
      w_min_idx <- sum(p@w < w_min)
      l_max = max(catchdist$length[catchdist$species == input$sp])
      w_max = a * l_max ^ b
      w_max_idx <- sum(p@w <= w_max)
    } else {
      w_min_idx <- sum(p@w < (p@species_params$w_mat[sp] / 100))
      w_max_idx <- sum(p@w <= p@species_params$w_inf[sp])
    }
    w_sel <- seq(w_min_idx, w_max_idx, by = 1)
    w <- p@w[w_sel]
    l = (p@w[w_sel] / a) ^ (1 / b)
    
    catch_w <- getFMort(p, effort = 1)[sp, w_sel] * 
      p@initial_n[sp, w_sel]
    # We just want the distribution, so we rescale the density so its area is 1
    catch_w <- catch_w / sum(catch_w * p@dw[w_sel])
    # The catch density in l gets an extra factor of dw/dl
    catch_l <- catch_w * b * w / l
    df <- data.frame(w, l, catch_w, catch_l, type = "Model catch")
    
    # We also include the abundance density because that helps to understand
    # the catch density    
    catch_w <- p@initial_n[sp, w_sel]
    # We just want the distribution, so we rescale the density so its area is 1
    catch_w <- catch_w / sum(catch_w * p@dw[w_sel])
    # The catch density in l gets an extra factor of dw/dl
    catch_l <- catch_w * b * w / l
    df <- rbind(df, data.frame(w, l, catch_w, catch_l, type = "Abundance"))
    
    if (is_observed) {
      # The observed catch is binned in bins equally spaced in length.
      # We need that binsize to normalise the density
      l <- catchdist$length[catchdist$species == input$sp]
      binsize <- min(diff(l))
      catch_l <- catchdist$catch[catchdist$species == input$sp]
      catch_l <- catch_l / sum(catch_l * binsize)
      # To get the density in w we need to divide by dw/dl
      w <- a * l ^ b
      catch_w <- catch_l / b * l / w
      df <- rbind(df, data.frame(w, l, catch_w, catch_l, 
                                 type = "Observed catch"))
    }
    
    if (input$catch_x == "Weight") {
      mat  <- p@species_params$w_mat[sp]
      pl <- ggplot(df) +
        geom_line(aes(x = w, y = catch_w, color = type)) +
        geom_text(aes(x = mat, y = max(catch_w * 0.9), label = "\nMaturity"), 
                  angle = 90) +
        labs(x = "Size [g]", y = "Density")
    } else {
      mat <- (p@species_params$w_mat[sp] / a) ^ (1 / b)
      pl <- ggplot(df) +
        geom_line(aes(x = l, y = catch_l, color = type)) +
        geom_text(aes(x = mat, y = max(catch_l * 0.9), label = "\nMaturity"), 
                  angle = 90) +
        labs(x = "Size [cm]", y = "Density")
    }
    pl +
      scale_x_log10() +
      geom_vline(xintercept = mat, linetype = "dotted")  +
      theme_grey(base_size = 18)
  })
  
  # Total catch by species
  output$plotTotalCatch <- renderPlotly({
    p <- params()
    biomass <- sweep(p@initial_n, 2, p@w * p@dw, "*")
    catch <- rowSums(biomass * getFMort(p, effort = 1))
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
               position = "dodge") +
      theme_grey(base_size = 18) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(x = "", y = "Catch [megatonnes]")
  })
  
  # Input field for observed catch
  output$catch_sel <- renderUI({
    p <- isolate(params())
    sp <- input$sp
    numericInput("catch_observed", 
                 paste0("Observed total catch for ", sp, " (megatonnes)"),
                 value = p@species_params[sp, "catch_observed"])
  })
  
  ## Plot rates ####  
  output$plotGrowth <- renderPlot({
  req(input$sp)
  sp <- input$sp
  p <- params()
  
  max_w <- p@species_params[sp, "w_inf"]
  if (input$axis == "Logarithmic") {
    min_w <- p@species_params[sp, "w_min"]
  } else {
    min_w = p@species_params[sp, "w_mat"] / 10 # min(1, p@species_params[sp, "w_min"])
  }
  sel <- p@w >= min_w & p@w <= max_w
  len <- sum(sel)
  growth <- getEGrowth(p, p@initial_n, p@initial_n_pp)[sp,sel]
  growth_and_repro <- getEReproAndGrowth(p, p@initial_n, p@initial_n_pp)[sp,sel]
  metab <- p@metab[sp,sel]
  income <- growth_and_repro + metab
  repro <- growth_and_repro - growth
  df <- data.frame(
    w = rep(p@w[sel], 4),
    Type = c(rep("Growth", len),
             rep("Income", len),
             rep("Metabolic loss", len),
             rep("Reproduction", len)),
    value = c(growth, income, metab, repro)
  )
  pl <- ggplot(df, aes(x = w, y = value, color = Type)) + 
    geom_line() + 
    geom_vline(xintercept = p@species_params[sp, "w_mat"], 
               linetype = "dotted") +
    geom_vline(xintercept = p@species_params[sp, "w_inf"], 
               linetype = "dotted") +
    theme_grey(base_size = 18) +
    labs(x = "Size [g]", y = "Rate")  +
    geom_text(aes(x = p@species_params[sp, "w_mat"], 
                  y = max(value * 0.2),
                  label = "\nMaturity"), 
              angle = 90)  +
    geom_text(aes(x = p@species_params[sp, "w_inf"], 
                  y = max(value * 0.2),
                  label = "\nMaximum"), 
              angle = 90)
  if (input$axis == "Logarithmic") {
    pl <- pl + scale_x_log10()
  }
  pl
})
  output$plotDeath <- renderPlot({
    req(input$sp)
    sp <- input$sp
    p <- params()
    
    max_w <- p@species_params[sp, "w_inf"]
    if (input$axis == "Logarithmic") {
      min_w <- p@species_params[sp, "w_min"]
    } else {
      min_w = p@species_params[sp, "w_mat"] / 10# min(1, p@species_params[sp, "w_min"])
    }
    sel <- p@w >= min_w & p@w <= max_w
    len <- sum(sel)
    df <- data.frame(
      w = rep(p@w[sel], 4),
      Type = c(rep("Total", len),
               rep("Predation", len),
               rep("Fishing", len),
               rep("Background", len)),
      value = c(getMort(p, p@initial_n, p@initial_n_pp,
                        effort = 1)[sp,sel],
                getPredMort(p, p@initial_n, p@initial_n_pp)[sp, sel],
                getFMort(p, effort = 1)[sp, sel],
                p@mu_b[sp,sel])
    )
    pl <- ggplot(df, aes(x = w, y = value, color = Type)) + 
      geom_line() + 
      geom_vline(xintercept = p@species_params[sp, "w_mat"], 
                 linetype = "dotted") +
      geom_vline(xintercept = p@species_params[sp, "w_inf"], 
                 linetype = "dotted") +
      theme_grey(base_size = 18) +
      labs(x = "Size [g]", y = "Rate")  +
      geom_text(aes(x = p@species_params[sp, "w_mat"], 
                    y = max(value * 0.2),
                    label = "\nMaturity"), 
                angle = 90)  +
      geom_text(aes(x = p@species_params[sp, "w_inf"], 
                    y = max(value * 0.2),
                    label = "\nMaximum"), 
                angle = 90)
    if (input$axis == "Logarithmic") {
      pl <- pl + scale_x_log10()
    }
    pl
  })
  
  ## Plot prey ####
  output$pred_size_slider <- renderUI({
    p <- isolate(params())
    sp <- which.max(p@species_params$species == input$sp)
    sliderInput("pred_size", "log predator size",
                value = signif(log(p@species_params$w_mat[sp]), 2),
                min = signif(log(p@species_params$w_min[sp]), 2),
                max = signif(log(p@species_params$w_inf[sp]), 2),
                step = 0.2,
                width = "80%",
                animate = animationOptions(loop = TRUE))
  })
  output$plot_prey <- renderPlot({
    p <- params()
    sp <- which.max(p@species_params$species == input$sp)
    phi <- function(x, xp) {
      phi <- exp(-(x - xp + log(p@species_params$beta[sp])) ^ 2 /
                   (2 * p@species_params$sigma[sp] ^ 2))
      phi[x >= xp | x < (xp - log(p@species_params$beta[sp]) -
                           3 * p@species_params$sigma[sp])] <- 0
      return(phi)
    }
    x <- log(p@w_full)
    dx <- x[2] - x[1]
    xp <- req(input$pred_size)
    wp <- exp(xp)
    wp_idx <- sum(p@w <= wp)
    # Calculate total community abundance
    fish_idx <- (length(p@w_full) - length(p@w) + 1):length(p@w_full)
    total_n <- p@initial_n_pp
    total_n[fish_idx] <- total_n[fish_idx] + colSums(p@initial_n)
    totalx <- total_n * p@w_full
    totalx <- totalx / sum(totalx * dx)
    phix <- phi(x, xp)
    phix <- phix / sum(phix * dx)
    pr <- totalx * phix
    br <- pr * p@w_full
    pr <- pr / sum(pr * dx)
    br <- br / sum(br * dx)
    df <- tibble(x, Kernel = phix, Biomass = br, Abundance = pr) %>%
      gather(type, y, Kernel, Biomass, Abundance)
    ggplot(df) +
      geom_line(aes(x, y, color = type)) +
      labs(x = "log(w)", y = "Density") +
      geom_point(aes(x = xp, y = 0), size = 4, colour = "blue")
  })
  
  ## Plot predators ####
  output$plot_pred <- renderPlotly({
    req(input$sp)
    sp <- input$sp
    p <- params()
    fish_idx <- (length(p@w_full) - length(p@w) + 1):length(p@w_full)
    pred_rate <- p@interaction[, sp] * 
      getPredRate(p, p@initial_n, p@initial_n_pp)[, fish_idx]
    fishing <- getFMort(p, 1)[sp, ]
    total <- colSums(pred_rate) + p@mu_b[sp, ] + fishing
    pred_rate <- pred_rate / rep(total, each = dim(pred_rate)[[1]])
    background <- p@mu_b[sp, ] / total
    fishing <- fishing / total
    # Make data.frame for plot
    plot_dat <- 
      rbind(
        data.frame(value = c(pred_rate),
                           Cause = as.factor(dimnames(pred_rate)[[1]]),
                           w = rep(p@w, each = dim(pred_rate)[[1]])),
        data.frame(value = background,
                   Cause = "Background",
                   w = p@w),
        data.frame(value = fishing,
                   Cause = "Fishing",
                   w = p@w)
      )
    ggplot(plot_dat) +
      geom_line(aes(x = w, y = value, color = Cause)) +
      scale_x_log10() +
      labs(x = "Size [g]", y = "Proportion of all death")
  })
  
  ## Plot psi ####
  output$plot_psi <- renderPlot({
    p <- params()
    sp <- which.max(p@species_params$species == input$sp)
    w_min <- p@species_params$w_inf[sp] / 50
    sel <- p@w >= w_min & p@w <= p@species_params$w_inf[sp]
    df <- data.frame(Size = p@w[sel], psi = p@psi[sp, sel])
    ggplot(df, aes(x = Size, y = psi)) + 
      geom_line(color = "blue") +
      geom_vline(xintercept = p@species_params[sp, "w_mat"], 
                 linetype = "dotted") +
      theme_grey(base_size = 18) +
      labs(x = "Size [g]", y = "Proportion of energy for reproduction")  +
      geom_text(aes(x = p@species_params[sp, "w_mat"], 
                    y = max(psi * 0.8),
                    label = "\nMaturity"), 
                angle = 90)
  })
  
} #the server

#### User interface ####
ui <- fluidPage(
  
  # titlePanel("Humboldt current ecosystem"),
  
  sidebarLayout(
    
    ## Sidebar ####
    sidebarPanel(
      tabsetPanel(
        id = "sidebarTabs",
        tabPanel(
          "Species",
          tags$br(),
          actionButton("sp_interact", "Interact"),
          actionButton("sp_steady", "Steady"),
          tags$br(),
          uiOutput("sp_sel"),
          uiOutput("sp_params"),
          tags$head(tags$style(
            type = 'text/css',
            '#sp_params { max-height: 70vh; overflow-y: auto; }'
          ))
        ),
        # tabPanel(
        #   "General",
        #   tags$br(),
        #   actionButton("bg_go", "Go"),
        #   sliderInput("effort", "Effort",
        #               value = 1, min = 0, max = 2, step = 0.05),
        #   uiOutput("general_params")
        # ),
        tabPanel(
          "File",
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
      tabsetPanel(id = "mainTabs",
        type = "tabs",
        tabPanel("Spectra", plotOutput("plotSpectra"),
                 radioButtons("binning", "Binning:",
                              choices = c("Logarithmic", "Constant"), 
                              selected = "Logarithmic", inline = TRUE)
        ),
        tabPanel("Biomass",
                 plotOutput("plotTotalBiomass"),
                 uiOutput("biomass_sel"),
                 plotOutput("plotBiomassDist")),
        tabPanel("Growth",
                 plotOutput("plotGrowthCurve"),
                 uiOutput("k_vb_sel")),
        tabPanel("Repro",
                 plotOutput("plot_erepro")),
        tabPanel("Catch",
                 plotlyOutput("plotTotalCatch"),
                 uiOutput("catch_sel"),
                 plotOutput("plotCatchDist"),
                 radioButtons("catch_x", "Show size in:",
                              choices = c("Weight", "Length"), 
                              selected = "Length", inline = TRUE)),
        tabPanel("Rates",
                 radioButtons("axis", "x-axis scale:",
                              choices = c("Logarithmic", "Normal"), 
                              selected = "Logarithmic", inline = TRUE),
                 plotOutput("plotGrowth"),
                 plotOutput("plotDeath")),
        tabPanel("f",
                 plotlyOutput("plot_feeding_level")),
        tabPanel("Prey",
                 uiOutput("pred_size_slider"),
                 plotOutput("plot_prey")),
        tabPanel("Predators",
                 plotlyOutput("plot_pred")),
        tabPanel("psi",
                 plotOutput("plot_psi"))
      )
    )  # end mainpanel
  )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app


