library(shiny)
library(shinyBS)
library(ggplot2)
# # Uncomment the following 3 lines before publishing the app
# library(devtools)
# install_github("gustavdelius/mizer", ref="scaling")
# library(mizer)
library(progress)

#### Server ####
server <- function(input, output, session) {
    
    # Show Introduction
    toggleModal(session, "Intro", toggle = "open")
    
    # Fishing parameters
    effort <- 0.4 
    # Gear parameters from Francesc Maynou
    l50 <- c(16.6, 15.48, 20.5, 15.85)
    names(l50) <- c("hake_old", "mullet_old", "hake_new", "mullet_new")
    sd <- c(0.462, 2.10, 0.331, 2.05)
    l25 = l50 - log(3) * sd
    
    # Load previously calculated simulation of hake_mullet system with old gear
    sim_old <- readRDS(file="hake_mullet.RData")
    # Therefore we do not need the following:
    # p_bg  <- reactive({
    #     setBackground(set_scaling_model(no_sp = input$no_sp, no_w = 400,
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
    #                     effort = effort, rfac = rfac)
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
    #                effort = input$effort, rfac = 1.01)
    #     p
    # })
    # 
    # sim_old <- reactive({
    #     # Create a Progress object
    #     progress <- shiny::Progress$new(session)
    #     on.exit(progress$close())
    #     
    #     project(p(), t_max = 50, t_save = 5, effort = effort, 
    #             shiny_progress = progress)
    # })
    
    sim <- reactive({
        # sim <- si()
        p <- sim_old@params
        no_sp <- length(p@species_params$species)
        no_t <- dim(sim_old@n)[1]
        p@initial_n <- sim_old@n[no_t, , ]
        p@initial_n_pp <- sim_old@n_pp[no_t, ]
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
        
        project(p, t_max = 15, t_save = 0.1, effort = input$effort, 
                shiny_progress = progress)
    })
    
    output$plotYield <- renderPlot({
        # sim_old <- si()
        y <- getYield(sim())
        y[1, ] <- getYield(sim_old)[dim(sim_old@n)[1], ]
        ym <- reshape2::melt(y, varnames = c("Year", "Species"), 
                             value.name = "Yield")
        ym <- subset(ym, ym$Yield > 0)
        ggplot(ym) + 
            geom_line(aes(x=Year, y=Yield, colour=Species, linetype=Species)) +
            scale_y_continuous(name="Yield [ton/year]", limits = c(0, NA)) +
            scale_colour_manual(values = sim()@params@linecolour) +
            scale_linetype_manual(values = sim()@params@linetype) +
            theme(text = element_text(size = 18))
    })
    
    output$plotChange <- renderPlot({
        # Plot changes in abundance
        s2 <- sim()
        p <- s2@params
        year <- input$change_year
        no_w <- length(p@w)
        w_sel <- seq.int(1, no_w, by = floor(no_w/50))
        w <- p@w[w_sel]
        change <- s2@n[10*year+1, ,w_sel]/s2@n[1, ,w_sel] - 1
        change_total <- colSums(s2@n[10*year+1, ,w_sel], na.rm = TRUE) /
            colSums(s2@n[1, ,w_sel], na.rm = TRUE) - 1
        ch <- rbind(change, "Total" = change_total)
        names(dimnames(ch)) <- names(dimnames(change))
        cf <- reshape2::melt(ch)
        cf <- subset(cf, !is.nan(value))
        names(cf)[1] <- "Species"
        ggplot(cf, aes(x = w, y = value)) +
            geom_line(aes(colour = Species, linetype = Species)) +
            geom_hline(yintercept = 0) +
            scale_x_log10(name = "Size [g]", labels = prettyNum,
                          breaks = 10^(-3:4)) +
            scale_y_continuous(name = "Percentage change", limits = c(-0.50, 0.60),
                               labels = scales::percent, breaks = (-7:9)/10) +
            scale_colour_manual(values = p@linecolour) +
            scale_linetype_manual(values = p@linetype) +
            # geom_vline(xintercept = l50) +
            # geom_vline(xintercept = l50_old) +
            theme(text = element_text(size = 18))
    })
    
    output$plotSelectivity <- renderPlot({
        p <- sim()@params
        w_min_idx <- sum(p@w < 0.5)
        w_max_idx <- which.min(p@w < 200)
        w_sel <- seq(w_min_idx, w_max_idx, by = floor((w_max_idx-w_min_idx)/50))
        w <- p@w[w_sel]
        selectivity <- p@selectivity[2, , w_sel]
        sf <- reshape2::melt(selectivity)
        sf$Gear <- "Modfied"
        selectivity_old <- sim_old@params@selectivity[2, , w_sel]
        sf_old <- reshape2::melt(selectivity_old)
        sf_old$Gear <- "Standard"
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
    
    output$plotCatch <- renderPlot({
        year <- input$catch_year
        s2 <- sim()
        p <- s2@params
        w_min_idx <- sum(p@w < 4)
        w_max_idx <- which.min(p@w < 200)
        w_sel <- seq(w_min_idx, w_max_idx, by = 1)
        w <- p@w[w_sel]
        catch <- p@selectivity[2, , w_sel] * s2@n[10*year+1, ,w_sel] * input$effort
        catchf <- reshape2::melt(catch)
        catchf$Gear <- "Modfied"
        catch_old <- sim_old@params@selectivity[2, , w_sel] * 
            sim_old@n[11, ,w_sel] * input$effort
        catchf_old <- reshape2::melt(catch_old)
        catchf_old$Gear <- "Standard"
        catchf <- rbind(catchf, catchf_old)
        names(catchf)[1] <- "Species"
        catchf <- subset(catchf, value > 0)
        if (input$catch_x == "Weight") {
            ggplot(catchf, aes(x = w, y = value)) +
                geom_line(aes(colour = Species, linetype = Gear)) +
                scale_x_continuous(name = "Size [g]", labels = prettyNum) +
                scale_y_continuous(name = "Yield distribution") +
                scale_colour_manual(values = s2@params@linecolour) +
                theme(text = element_text(size = 18))
        } else {
            a <- p@species_params$a
            names(a) <- p@species_params$species
            b <- p@species_params$b
            names(b) <- p@species_params$species
            ggplot(catchf, aes(x = (w/a[Species])^(1/b[Species]), y = value)) +
                geom_line(aes(colour = Species, linetype = Gear)) +
                scale_x_continuous(name = "Length [cm]", labels = prettyNum) +
                scale_y_continuous(name = "Yield distribution") +
                scale_colour_manual(values = p@linecolour) +
                theme(text = element_text(size = 18))
        }
    })
    
    output$plotSSB <- renderPlot({
        b <- getSSB(sim())[, 11:12]
        bm <- reshape2::melt(b, varnames = c("Year", "Species"), 
                             value.name = "SSB")
        bm <- subset(bm, bm$SSB > 0)
        ggplot(bm) + 
            geom_line(aes(x=Year, y=SSB, colour=Species, linetype=Species)) +
            scale_y_continuous(name="SSB [ton]", limits = c(0, NA)) +
            scale_colour_manual(values = sim()@params@linecolour) +
            scale_linetype_manual(values = sim()@params@linetype) +
            theme(text = element_text(size = 18))
    })
    
    output$plotSpectrum <- renderPlot({
        plotSpectra(sim(), time_range = input$spectrum_year)
    })
    
}

#### user interface ####
ui <- fluidPage(
    
    titlePanel("Effect of modfied gear"),
    
    bsModal("Intro", "Welcome", "introBut", size = "large",
            "Here comes some explanation"),
    
    sidebarLayout(
        
        sidebarPanel(
            h3("Modified fishing"),
            sliderInput("effort", "Effort",
                        value=0.4, min=0.3, max=0.5),
            bsTooltip("effort", "Explanation of this slider", "right"),
            h4("Hake selectivity"),
            sliderInput("l50_hake", "L50", post = "cm",
                        value=20.5, min=12, max=24, step = 0.01),
            sliderInput("l25_hake", "L25", post = "cm",
                        value=20.1, min=12, max=24, step = 0.01),
            h4("Mullet selectivity"),
            sliderInput("l50_mullet", "L50", post = "cm",
                        value=15.8, min=12, max=24, step = 0.01),
            sliderInput("l25_mullet", "L25", post = "cm",
                        value=13.6, min=12, max=24, step = 0.01),
            img(src = "logo_minouw_blue.png", width = "200px"),
            actionButton("introBut", "View Introduction again")
        ),  # endsidebarpanel
        
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Intro", "Introductory text"),
                        tabPanel("Yield", plotOutput("plotYield")),
                        tabPanel("Spectrum",
                                 wellPanel(
                                     sliderInput("spectrum_year", "Year",
                                                 value = 5, min = 0, max = 15, step = 1, 
                                                 animate = TRUE)
                                 ),
                                 plotOutput("plotSpectrum")
                        ),
                        tabPanel("Change",
                                 wellPanel(
                                     sliderInput("change_year", "Year",
                                                 value = 5, min = 0, max = 15, step = 1, 
                                                 animate = TRUE)
                                 ),
                                 plotOutput("plotChange")
                        ),
                        tabPanel("SSB", plotOutput("plotSSB")),
                        tabPanel("Selectivity",
                                 plotOutput("plotSelectivity"),
                                 radioButtons("selectivity_x", "Show size in:",
                                              choices = c("Weight", "Length"), 
                                              selected = "Length") 
                        ),
                        tabPanel("Catch",
                                 wellPanel(
                                     sliderInput("catch_year", "Year",
                                                 value = 5, min = 0, max = 15, step = 1, 
                                                 animate = TRUE)
                                 ),
                                 plotOutput("plotCatch"),
                                 radioButtons("catch_x", "Show size in:",
                                              choices = c("Weight", "Length"), 
                                              selected = "Length") 
                        )
                        #tabPanel("Spectra", plotOutput("plotSpectra"))
            )
        )  # end mainpanel
    )  # end sidebarlayout
)

shinyApp(ui = ui, server = server) # this launches your app