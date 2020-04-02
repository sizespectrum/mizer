library(shiny)
library(shinyBS)
library(ggplot2)
library(mizer)
library(progress)

#### Server ####
server <- function(input, output, session) {
    
    # Show Introduction
    toggleModal(session, "Intro", toggle = "open")
    
    # Fishing parameters
    fixed_effort <- 0.4 
    # Gear parameters from Francesc Maynou
    l50 <- c(16.6, 15.48, 20.5, 15.85)
    names(l50) <- c("hake_old", "mullet_old", "hake_new", "mullet_new")
    sd <- c(0.462, 2.10, 0.331, 2.05)
    l25 = l50 - log(3) * sd
    
    # Load previously calculated simulation of hake_mullet system with old gear
    sim_old <- readRDS(file="hake_mullet.RDS")
    # Therefore we do not need the following:
    # p_bg  <- reactive({
    #     markBackground(set_scaling_model(no_sp = input$no_sp, no_w = 400,
    #                     min_w_inf = 10, max_w_inf = 1e5,
    #                     min_egg = 1e-4, min_w_mat = 10^(0.4),
    #                     knife_edge_size = Inf, kappa = 10000,
    #                     lambda = 2.08, f0 = 0.6, h = 34))
    # })
    # 
    # p_old <- reactive({
    #     rfac <- 1.01
    #     ######### add mullet
    #     a_m <- 0.0085
    #     b_m <- 3.11
    #     L_inf_m <- 24.3
    #     # http://www.fishbase.org/summary/Mullus-barbatus+barbatus.html
    #     L_mat <- 11.1
    #     species_params <- data.frame(
    #         species = "Mullet",
    #         w_min = 0.001,
    #         w_inf = a_m*L_inf_m^b_m, # from fishbase
    #         w_mat = a_m*L_mat^b_m, # from fishbase
    #         beta = 283, # = beta_gurnard from North sea. Silvia says gurnard is similar.
    #         sigma = 1.8, # = sigma_gurnard from North sea. Silvia says gurnard is similar.
    #         z0 = 0,
    #         alpha = 0.6, # unknown, mizer default=0.6
    #         erepro = 0.1, # unknown, determined later
    #         sel_func = "sigmoid_length",
    #         gear = "sigmoid_gear",
    #         l25 = l25["mullet_old"],
    #         l50 = l50["mullet_old"],
    #         k = 0,
    #         k_vb = 0.6,
    #         a = a_m,
    #         b = b_m,
    #         gamma = 0.0017,
    #         h = 50,
    #         linecolour = "red",
    #         linetype = "solid"
    #     )
    #     p <- addSpecies(p_bg(), species_params, SSB = 2800, 
    #                     effort = fixed_effort, rfac = rfac)
    #     
    #     ############# add hake 
    #     # Merluccius merluccius  (European hake)
    #     a <- 0.0046
    #     b <- 3.12
    #     L_inf <- 81.2
    #     L_mat <- 29.83
    #     
    #     species_params <- data.frame(
    #         species = "Hake",
    #         w_min = 0.001,
    #         w_inf = a*L_inf^b, # from fishbase
    #         w_mat = a*L_mat^b, # from fishbase
    #         beta = exp(2.4), #RLD and Blanchard thesis p 88
    #         sigma = 1.1, #RLD and Blanchard thesis p 88
    #         z0 = 0,
    #         alpha = 0.6, # unknown, using mizer default=0.6
    #         erepro = 0.1, # unknown
    #         sel_func = "sigmoid_length", # not used but required
    #         gear = "sigmoid_gear",
    #         l25 = l25["hake_old"],
    #         l50 = l50["hake_old"],
    #         k = 0,
    #         k_vb = 0.1, # from FB website below
    #         a = a,
    #         b = b,
    #         gamma = 0.003,
    #         h = 20,
    #         linecolour = "blue",
    #         linetype = "solid"
    #     )
    #     p <- addSpecies(p, species_params, SSB = input$SSB_hake, 
    #                effort = fixed_effort, rfac = 1.01)
    #     p
    # })
    # 
    # sim_old <- reactive({
    #     # Create a Progress object
    #     progress <- shiny::Progress$new(session)
    #     on.exit(progress$close())
    #     
    #     project(p(), t_max = 50, t_save = 5, effort = fixed_effort, 
    #             progress_bar = progress)
    # })
    
    # Data frame for yield plot
    no_sp <- length(sim_old@params@species_params$species)
    ym_old <- data.frame(
        "Year" = rep(c(2018, 2033), each = no_sp),
        "Species" = rep(sim_old@params@species_params$species, times = 2),
        "Yield" = rep(getYield(sim_old)[11, ], times = 2),
        "Gear" = "Current"
    )
    ym_old <- subset(ym_old, ym_old$Yield > 0)
    
    # Data frame for SSB plot
    bm_old <- data.frame(
        "Year" = rep(c(2018, 2033), each = 2),
        "Species" = rep(sim_old@params@species_params$species[11:12], times = 2),
        "SSB" = rep(getSSB(sim_old)[11, 11:12], times = 2),
        "Gear" = "Current"
    )
    
    # Set params ####
    params <- reactive({
        # sim <- si()
        p <- sim_old@params
        no_t <- dim(sim_old@n)[1]
        p@initial_n <- sim_old@n[no_t, , ]
        p@initial_n_pp <- sim_old@n_pp[no_t, ]
        p@species_params$R_max <- Inf
        
        # Retune the values of erepro so that we get the correct level of
        # reproduction without density dependence
        mumu <- getZ(p, p@initial_n, p@initial_n_pp, effort = fixed_effort)
        gg <- getEGrowth(p, p@initial_n, p@initial_n_pp)
        rdi <- getRDI(p, p@initial_n, p@initial_n_pp)
        for (i in (1:no_sp)) {
            gg0 <- gg[i, p@w_min_idx[i]]
            mumu0 <- mumu[i, p@w_min_idx[i]]
            DW <- p@dw[p@w_min_idx[i]]
            p@species_params$erepro[i] <- p@species_params$erepro[i] *
                (p@initial_n[i, p@w_min_idx[i]] *
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
        p
    })
    
    # Run simulation ####
    sim <- reactive({        
        # Create a Progress object
        progress <- shiny::Progress$new(session)
        on.exit(progress$close())
        
        project(params(), t_max = 15, t_save = 0.1, effort = fixed_effort, 
                progress_bar = progress)
    })
    
    # Plot yield ####
    output$plotYield <- renderPlot({
        y <- getYield(sim())
        ym <- reshape2::melt(y, varnames = c("Year", "Species"), 
                             value.name = "Yield")
        ym <- subset(ym, ym$Yield > 0)
        ym$Gear <- "Modified"
        ym$Year <- ym$Year + 2018
        ym <- rbind(ym_old, ym)
        ggplot(ym) + 
            geom_line(aes(x = Year, y = Yield, colour = Species, linetype = Gear)) +
            scale_y_continuous(name="Yield [tonnes/year]", limits = c(0, NA)) +
            scale_colour_manual(values = params()@linecolour) +
            scale_linetype_manual(values = c("Current" = "dotted", "Modified" = "solid")) +
            theme(text = element_text(size = 18))
    })
    
    # Plot SSB ####
    output$plotSSB <- renderPlot({
        b <- getSSB(sim())[, 11:12]
        bm <- reshape2::melt(b, varnames = c("Year", "Species"), 
                             value.name = "SSB")
        bm$Gear <- "Modified"
        bm$Year <- bm$Year + 2018
        bm <- rbind(bm_old, bm)
        ggplot(bm) + 
            geom_line(aes(x = Year, y = SSB, colour = Species, linetype = Gear)) +
            scale_y_continuous(name="SSB [tonnes]", limits = c(0, NA)) +
            scale_colour_manual(values = params()@linecolour) +
            scale_linetype_manual(values = c("Current" = "dotted", "Modified" = "solid")) +
            theme(text = element_text(size = 18))
    })
    
    # Plot percentage change ####
    output$plotChange <- renderPlot({
        # Plot changes in abundance
        s2 <- sim()
        p <- params()
        year <- input$change_year - 2018
        no_w <- length(p@w)
        w_sel <- seq.int(1, no_w, by = floor(no_w/50))
        w <- p@w[w_sel]
        change <- s2@n[10*year+1, ,w_sel]/s2@n[1, ,w_sel] - 1
        # change_total <- colSums(s2@n[10*year+1, ,w_sel], na.rm = TRUE) /
        #                      colSums(s2@n[1, ,w_sel], na.rm = TRUE) - 1
        # ch <- rbind(change, "Total" = change_total)
        # names(dimnames(ch)) <- names(dimnames(change))
        cf <- reshape2::melt(change)
        cf <- subset(cf, !is.nan(value))
        cf$Species <- as.character(cf$sp)
        cf$Species[cf$Species %in% 1:10] <- "Background"
        
        # data frame for special points
        w_mat <- p@species_params$w_mat[11:12]
        w50 <- p@species_params$a[11:12] * 
            (p@species_params$l50[11:12])^p@species_params$b[11:12]
        sp <- data.frame("w" = c(w_mat, w50),
                         "y" = c(change[11, which.min(w < w_mat[1])],
                                 change[12, which.min(w < w_mat[2])],
                                 change[11, which.min(w < w50[1])],
                                 change[12, which.min(w < w50[2])]),
                         "Points" = c("Maturity", "Maturity", "L50", "L50"),
                         "Species" = p@species_params$species[11:12])
        
        ggplot(cf, aes(x = w, y = value)) +
            geom_line(aes(colour = Species, linetype = Species, group = sp)) +
            geom_hline(yintercept = 0) +
            scale_x_log10(name = "Size [g]", labels = prettyNum,
                          breaks = 10^(-3:4)) +
            scale_y_continuous(name = "Percentage change", limits = c(-0.50, 0.60),
                               labels = scales::percent, breaks = (-7:9)/10) +
            scale_colour_manual(values = p@linecolour) +
            scale_linetype_manual(values = p@linetype) +
            theme(text = element_text(size = 14)) +
            geom_point(aes(x = w, y = y, colour = Species, shape = Points), 
                       data = sp, size=3)
    })
    
    # Plot selectivity curves ####
    output$plotSelectivity <- renderPlot({
        p <- params()
        w_min_idx <- sum(p@w < 0.5)
        w_max_idx <- which.min(p@w < 200)
        w_sel <- seq(w_min_idx, w_max_idx, by = floor((w_max_idx-w_min_idx)/50))
        w <- p@w[w_sel]
        selectivity <- p@selectivity[2, , w_sel]
        sf <- reshape2::melt(selectivity)
        sf$Gear <- "Modfied"
        selectivity_old <- sim_old@params@selectivity[2, , w_sel]
        sf_old <- reshape2::melt(selectivity_old)
        sf_old$Gear <- "Current"
        sf <- rbind(sf, sf_old)
        names(sf)[1] <- "Species"
        sf <- subset(sf, value > 0)
        if (input$selectivity_x == "Weight") {
            ggplot(sf, aes(x = w, y = value)) +
                geom_line(aes(colour = Species, linetype = Gear)) +
                scale_x_continuous(name = "Size [g]", labels = prettyNum) +
                scale_y_continuous(name = "Selectivity",
                                   labels = scales::percent) +
                scale_colour_manual(values = p@linecolour) +
                theme(text = element_text(size = 18))
        } else {
            a <- p@species_params$a
            names(a) <- p@species_params$species
            b <- p@species_params$b
            names(b) <- p@species_params$species
            ggplot(sf, aes(x = (w/a[Species])^(1/b[Species]), y = value)) +
                geom_line(aes(colour = Species, linetype = Gear)) +
                scale_x_continuous(name = "Length [cm]", labels = prettyNum) +
                scale_y_continuous(name = "Selectivity",
                                   labels = scales::percent) +
                scale_colour_manual(values = p@linecolour) +
                theme(text = element_text(size = 18))
        }
    })
    
    # Plot catch by size ####
    output$plotCatch <- renderPlot({
        year <- input$catch_year - 2018
        s2 <- sim()
        p <- params()
        w_min_idx <- sum(p@w < 4)
        w_max_idx <- which.min(p@w < 200)
        w_sel <- seq(w_min_idx, w_max_idx, by = 1)
        w <- p@w[w_sel]
        catch_old <- sim_old@params@selectivity[2, 11:12, w_sel] * 
            sim_old@n[11, 11:12,w_sel] * fixed_effort * rep(w, each = 2)
        catchf_old <- reshape2::melt(catch_old)
        catchf_old$Gear <- "Current"
        catch <- p@selectivity[2, 11:12, w_sel] * s2@n[10*year+1, 11:12,w_sel] * 
            fixed_effort * rep(w, each = 2)
        catchf <- reshape2::melt(catch)
        catchf$Gear <- "Modified"
        catchf <- rbind(catchf, catchf_old)
        names(catchf)[1] <- "Species"
        if (input$catch_x == "Weight") {
            ggplot(catchf, aes(x = w, y = value)) +
                geom_line(aes(colour = Species, linetype = Gear)) +
                scale_x_continuous(name = "Size [g]", labels = prettyNum) +
                scale_y_continuous(name = "Yield distribution") +
                scale_colour_manual(values = p@linecolour) +
                scale_linetype_manual(values = c("Current" = "dotted", 
                                                 "Modified" = "solid")) +
                theme(text = element_text(size = 18))
        } else {
            a <- p@species_params$a[11:12]
            names(a) <- p@species_params$species[11:12]
            b <- p@species_params$b[11:12]
            names(b) <- p@species_params$species[11:12]
            ggplot(catchf, aes(x = (w/a[Species])^(1/b[Species]), y = value)) +
                geom_line(aes(colour = Species, linetype = Gear)) +
                scale_x_continuous(name = "Length [cm]", labels = prettyNum) +
                scale_y_continuous(name = "Yield distribution") +
                scale_colour_manual(values = p@linecolour) +
                scale_linetype_manual(values = c("Current" = "dotted", 
                                                 "Modified" = "solid")) +
                theme(text = element_text(size = 18))
        }
    })
    
    # Plot size spectrum ####
    output$plotSpectrum <- renderPlot({
        plotSpectra(sim(), time_range = input$spectrum_year - 2018)
    })
    
    # Constraints on sliders ####
    # L25 must always be smaller than l50
    observe({
        if (input$l50_hake <= input$l25_hake) {
            updateSliderInput(session, "l25_hake", value = input$l50_hake-0.01)
            if (input$l50_hake <= 12) {
                updateSliderInput(session, "l50_hake", value = input$l50_hake+0.01)
            }
        }
    })
    observe({
        if (input$l50_mullet <= input$l25_mullet) {
            updateSliderInput(session, "l25_mullet", value = input$l50_mullet-0.01)
            if (input$l50_mullet <= 12) {
                updateSliderInput(session, "l50_mullet", value = input$l50_mullet+0.01)
            }
        }
    })
    
}

#### User interface ####
ui <- fluidPage(
    
    titlePanel("Gear modification: target species and the background ecosystem"),
    
    bsModal("Intro", "Welcome", "introBut", size = "large",
        p("When we change fishing on target species, we alter the ecosystems 
          in which the target species are embedded.  The ecosystem-based 
          approach to fisheries management recognises this, and tools are 
          needed to learn about these changes."),
        p("This web page is a tool for looking at these ecosystem effects.  
          It is an example based on Geographical SubArea GSA06 (Northern Spain), 
          an extension to a study by Maynou (2018). Specifically, we embed 
          chosen target species into a generic background ecosystem of multiple species, 
          allowing the target and background species to interact. The target 
          species of special importance in GSA06 are hake and red mullet.  
          A European discards ban on these two species is proposed, and we need 
          to understand the consequences of changing fishing gear in the 
          ecosystem context."),
        p("Calculations start in 2018 with the bottom-trawl fishing gear 
          currently in use, as given by Maynou (2018), and average abundances 
          given by the corresponding steady state.  This gear is called 
          'current' in the graphs.  At the start, a new fishing gear is chosen, 
          and an integration is run for a number of years.  As time goes on, 
          abundance and yield of the target species change: (a) directly through 
          the choice of fishing gear, and (b) indirectly through feedbacks with 
          other species."),
        p("Default settings of a new gear are those of a T90 extension net to 
          reduce the catches of undersize hake and red mullet (Maynou 2018). 
          Beyond this, there is flexibility for you to investigate the effect
          of gears with different selectivities."),
        p("The engine that drives the calculations is an updated version of 
          mizer, a software package for the dynamics of multispecies size 
          spectra.  This deals explicitly with the dynamics of marine species 
          that interact by feeding on each other.  We have modified mizer to 
          allow species with chosen life cycles to be inserted into a background 
          ecosystem comprising a an assemblage of multiple species."),
        p("This is a tool for you to explore effects of new kinds of fishing 
          gear.  The tool is in its early stages of development and is not, 
          at the moment, appropriate for making specific management decisions.  
          Like a weather forecast, prediction becomes more uncertain as you go 
          further into the future.  We welcome your feedback for improvements; 
          please contact: gustav.delius@york.ac.uk")
    ),
    
    sidebarLayout(
        
        sidebarPanel(
            h3("Modified fishing"),
            # sliderInput("effort", "Effort",
            #             value=0.4, min=0.3, max=0.5),
            # bsTooltip("effort", "Explanation of this slider", "right"),
            h4("Hake selectivity"),
            sliderInput("l50_hake", "L50", post = "cm",
                        value=20.5, min=12, max=22, step = 0.01),
            sliderInput("l25_hake", "L25", post = "cm",
                        value=20.1, min=12, max=22, step = 0.01),
            h4("Mullet selectivity"),
            sliderInput("l50_mullet", "L50", post = "cm",
                        value=15.8, min=12, max=22, step = 0.01),
            sliderInput("l25_mullet", "L25", post = "cm",
                        value=13.6, min=12, max=22, step = 0.01),
            img(src = "logo_minouw_blue.png", width = "200px"),
            actionButton("introBut", "View Introduction again")
        ),  # endsidebarpanel
        
        mainPanel(
            tabsetPanel(
                type = "tabs",
                tabPanel(
                    "Selectivity",
                    br(),
                    p("A comparison of the selectivity of the current and the 
                      modified gear on hake and mullet."),
                    p("You control the selectivity of the modified gear with 
                      the sliders on the left."),
                    plotOutput("plotSelectivity"),
                    radioButtons("selectivity_x", "Show size in:",
                                 choices = c("Weight", "Length"), 
                                 selected = "Length", inline = TRUE)
                ),
                tabPanel(
                    "Total Catch",
                    br(),
                    p("Yield from modified fishing compared with current fishing 
                      over time."),
                    plotOutput("plotYield"),
                    p("With the default settings, the yield of hake is initially 
                      below that of the current fishing.  As hake's size 
                      structure recovers, its yield increases.  Since the 
                      fishing mortality of red mullet remains close to that of 
                      current fishing, the decrease in red-mullet yield comes 
                      mostly from the change in abundance of hake and of 
                      background species.")
                ),
                tabPanel(
                    "Spectrum",
                    br(),
                    p("Abundance of species at a chosen time."),
                    wellPanel(
                        sliderInput("spectrum_year", "Year",
                                    value = 2023, min = 2018, 
                                    max = 2033, step = 1, 
                                    animate = TRUE, sep = "")
                    ),
                    plotOutput("plotSpectrum"),
                    p("This unpacks a bit of the underlying mizer model of 
                      multispecies dynamics.  The grey lines are background 
                      species, each with its own dynamics;  the green line is a 
                      resource resource spectrum without which the whole fish 
                      assemblage would go to extinction.  The target species 
                      both eat, and are eaten by, the background species."),
                    p("In 2018, the system is at the steady state under 
                      current fishing. At this time, your chosen fishing 
                      mortality is imposed on the target species 
                      (background species retain the old fishing).  As 
                      time goes on, the target species adjust to the new 
                      fishing. This adjustment leads to changes in the 
                      abundance of background species.  Since the species 
                      are coupled by predation, changes in the background 
                      species feed back to the target species.", 
                      "bottom")
                ),
                tabPanel(
                    "% Change",
                    br(),
                    p("Percentage change in abundance from 'current' to 
                      'modified' gear over time"),
                    wellPanel(
                        sliderInput("change_year", "Year",
                                    value = 2023, min = 2018, 
                                    max = 2033, step = 1, 
                                    animate = TRUE, sep = "")
                    ),
                    plotOutput("plotChange"),
                    p("With the default settings, the main effects on target 
                      species is to increase the abundance of hake, in the short 
                      term, and reduce the abundance of red mullet.  In 
                      addition, you can see how the change abundance of the 
                      target species, changes the background species, which in 
                      turn feeds back to the target species.")
                ),
                tabPanel(
                    "SSB",
                    br(),
                    p("Spawning stock biomass over time."),
                    plotOutput("plotSSB")
                ),
                tabPanel(
                    "Catch by Size",
                    br(),
                    p("Biomass catch as a function of size. This shows how the
                      improved selectivity leads to more of the caught biomass
                      to consist of larger fish."),
                    wellPanel(
                        sliderInput("catch_year", "Year",
                                    value = 2023, min = 2018, 
                                    max = 2033, step = 1, 
                                    animate = TRUE, sep = "")
                    ),
                    plotOutput("plotCatch"),
                    radioButtons("catch_x", "Show size in:",
                                 choices = c("Weight", "Length"), 
                                 selected = "Length") 
                )
            )
        
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app
