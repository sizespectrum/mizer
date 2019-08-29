
#' Launch shiny gadget for tuning parameters
#' 
#' The shiny gadget has sliders for the model parameters and tabs with a
#' variety of plots to visualise the steady-state
#' 
#' @param p MizerParams object to tune. If missing, the gadget tries to recover
#'   information from log files left over from aborted previous runs.
#' @param catch Data frame holding binned observed catch data. The data can
#'   be binned either into length bins or weight bins. In the former the data
#'   frame should have columns \code{length} and \code{dl} holding the start of
#'   the size bins in cm and the width of the size bins in cm respectively. In
#'   the latter case the data frame should have columns \code{weight} and
#'   \code{dw} holding the start of the size bins in frams and the width of the
#'   size bins in grams. The data frame also needs to have the columns
#'   \code{species} (the name of the species), \code{catch} (the number of
#'   individuals of a particular species caught in a size bin).
#' @param stomach Data frame holding observations of prey items in predator
#'   stomachs. The required columns are 
#'   \itemize{
#'   \item \code{species} holding the name of the predator species,
#'   \item \code{wpredator} with the weight in grams of the predator,
#'   \item \code{wprey} with the weight of the prey item.
#'   }
#'   In case prey items of the same weight have been aggregated in the data
#'   frame then there should be a column \code{Nprey} saying how many prey 
#'   items have been aggregated in each row.
#' 
#' @return The tuned MizerParams object
#' @export
tuneParams <- function(p, catch = NULL, stomach = NULL) {
    # Check arguments ----
    if (!is.null(catch)) {
        assert_that(
            is.data.frame(catch),
            "catch" %in% names(catch),
            "species" %in% names(catch),
            all(c("length", "dl") %in% names(catch)) |
                all(c("weight", "dw") %in% names(catch))
        )
    }
    if (!is.null(stomach)) {
        assert_that(
            is.data.frame(stomach),
            "wprey" %in% names(stomach),
            "wpredator" %in% names(stomach),
            "species" %in% names(stomach)
        )
        if (!("Nprey" %in% names(stomach))) stomach$Nprey <- 1
        stomach <- stomach %>% 
            mutate(logpredprey = log(wpredator / wprey),
                   weight = Nprey / sum(Nprey),
                   weight_biomass = Nprey * wprey / sum(Nprey * wprey),
                   weight_kernel = Nprey / wprey^(1 + alpha - lambda),
                   weight_kernel = weight_kernel / sum(weight_kernel))
    }
    
    # Define some globals to skip certain observers ----
    sp_old <- 1
    sp_old_kernel <- 1
    sp_old_predation <- 1
    sp_old_fishing <- 1
    sp_old_maturity <- 1
    sp_old_prey <- 1
    sp_old_n0 <- 1
    
    prepare_params <- function(p) {
        rownames(p@species_params) <- p@species_params$species
        p <- set_species_param_default(p, "a", 0.006)
        p <- set_species_param_default(p, "b", 3)
        p <- set_species_param_default(p, "t0", 0)
        p <- set_species_param_default(p, "sel_func", "sigmoid_length")
        rdi <- getRDI(p)
        rdd <- getRDD(p)
        p@species_params$erepro <- p@species_params$erepro * rdd / rdi
        p@species_params$r_max <- Inf
        return(p)
    }
    
    # Prepare logs for undo/redo functionality
    logs <- vector(mode = "character")
    log_idx <- 0
    add_to_logs <- function(p) {
        # Save params object to disk
        time = format(Sys.time(), "_%Y_%m_%d_at_%H_%M_%S")
        file = paste0(tempdir(), "/mizer_params", time, ".rds")
        saveRDS(p, file = file)
        # Update logs
        if (log_idx < length(logs)) {
            file.remove(logs[(log_idx + 1):length(logs)])
        }
        logs <<- append(logs[min(1, log_idx):log_idx], file)
        log_idx <<- log_idx + 1
        shinyjs::disable("redo")
        if (log_idx > 1) {
            shinyjs::enable("undo")
            shinyjs::enable("undo_all")
        }
    }
    
    if (missing(p)) {
        # Try to recover old log files ----
        logs <- sort(list.files(path = tempdir(), 
                                pattern = "mizer_params_...._.._.._at_.._.._..\\.rds",
                                full.names = TRUE))
        log_idx <- length(logs)
        if (log_idx == 0) {
            stop("You need to specify a MizerParams object. ",
                 "There are no temporary parameter files to recover.")
        }
        p <- readRDS(logs[log_idx])
    } else {
        validObject(p)
        p <- prepare_params(p)
    }
    
    # User interface ----
    ui <- fluidPage(
        shinyjs::useShinyjs(),
        # titlePanel("Humboldt current ecosystem"),
        
        sidebarLayout(
            
            ## Sidebar ####
            sidebarPanel(
                actionButton("sp_steady", "Steady"),
                actionButton("undo", "", icon = icon("undo")),
                actionButton("redo", "", icon = icon("redo")),
                actionButton("undo_all", "", icon = icon("fast-backward")),
                actionButton("done", "Done", icon = icon("check"),
                             onclick = "setTimeout(function(){window.close();},500);"),
                tags$br(),
                uiOutput("sp_sel"),
                "->",
                tags$a("Biomass", href = "#biomass"),
                "->",
                tags$a("Predation", href = "#predation"),
                "->",
                tags$a("Resources", href = "#resource"),
                "->",
                tags$a("Fishing", href = "#fishing"),
                "->",
                tags$a("Maturity", href = "#maturity"),
                "->",
                tags$a("Others", href = "#others"),
                "->",
                tags$a("Prey", href = "#interactions"),
                "->",
                tags$a("General", href = "#general"),
                "->",
                tags$a("Plankton", href = "#plankton"),
                "->",
                tags$a("File", href = "#file"),
                tags$br(),
                tags$div(id = "params",
                    uiOutput("sp_params"),
                    uiOutput("general_params")
                    ),
                tags$head(tags$style(
                    type = 'text/css',
                    '#params { max-height: 60vh; overflow-y: auto; }'
                )),
                width = 3
            ),  # endsidebarpanel
            
            ## Main panel ####
            mainPanel(
                tabsetPanel(id = "mainTabs",
                            type = "tabs",
                            tabPanel("Spectra", plotlyOutput("plotSpectra"),
                                     radioButtons("binning", "Binning:",
                                                  choices = c("Logarithmic", "Constant"), 
                                                  selected = "Logarithmic", inline = TRUE)
                            ),
                            tabPanel("Biomass",
                                     plotlyOutput("plotTotalBiomass"),
                                     plotlyOutput("plotTotalAbundance"),
                                     uiOutput("biomass_sel"),
                                     actionButton("tune_egg", "Tune egg density"),
                                     plotlyOutput("plotBiomassDist")),
                            tabPanel("Growth",
                                     radioButtons("all_growth", "Show:",
                                                  choices = c("All", "Selected species"), 
                                                  selected = "All", inline = TRUE),
                                     plotOutput("plotGrowthCurve",
                                                click = "growth_click"),
                                     textOutput("info"),
                                     uiOutput("k_vb_sel"),
                                     plotlyOutput("plot_feeding_level"),
                                     plotlyOutput("plot_psi")),
                            tabPanel("Repro",
                                     plotlyOutput("plot_erepro")),
                            tabPanel("Catch",
                                     actionButton("tune_catch", "Tune catchability"),
                                     uiOutput("catch_sel"),
                                     textOutput("catch_total"),
                                     plotlyOutput("plotCatchDist"),
                                     radioButtons("catch_x", "Show size in:",
                                                  choices = c("Weight", "Length"), 
                                                  selected = "Length", inline = TRUE),
                                     plotlyOutput("plotTotalCatch")
                            ),
                            tabPanel("Rates",
                                     radioButtons("axis", "x-axis scale:",
                                                  choices = c("Logarithmic", "Normal"), 
                                                  selected = "Logarithmic", inline = TRUE),
                                     plotlyOutput("plotGrowth"),
                                     plotlyOutput("plotDeath")),
                            tabPanel("Prey",
                                     uiOutput("pred_size_slider"),
                                     plotlyOutput("plot_prey")),
                            tabPanel("Diet",
                                     plotlyOutput("plot_diet")),
                            tabPanel("Death",
                                     radioButtons("death_prop", "Show",
                                                  choices = c("Proportion", "Rate"), 
                                                  selected = "Proportion", 
                                                  inline = TRUE),
                                     plotlyOutput("plot_pred")),
                            tabPanel("Plankton",
                                     plotOutput("plot_plankton", width = "84%"),
                                     plotlyOutput("plot_plankton_pred"),
                                     radioButtons("plankton_death_prop", "Show",
                                                  choices = c("Proportion", "Rate"), 
                                                  selected = "Proportion", 
                                                  inline = TRUE))
                            # tabPanel("Stomach",
                            #          plotOutput("plot_stomach"),
                            #          plotOutput("plot_kernel"))
                )
            )  # end mainpanel
        )  # end sidebarlayout
    )
    
    server <- function(input, output, session) {
        no_sp <- length(p@species_params$species)
        # Size of plot labels
        base_size <- 12
        # selector for foreground species
        foreground <- !is.na(p@A)
        foreground_indices <- (1:no_sp)[foreground]
        
        ## Store params object as a reactive value ####
        params <- reactiveVal(p)
        add_to_logs(p)
        if (log_idx == length(logs)) shinyjs::disable("redo")
        if (log_idx <= 1) {
            shinyjs::disable("undo")
            shinyjs::disable("undo_all")
        }
        output$filename <- renderText("")
        
        # Define a reactive value for triggering an update of species sliders
        trigger_update <- reactiveVal(0)
        # Fishing effort
        effort <- reactiveVal(1)
        
        ## Handle upload of params object ####
        observeEvent(input$upload, {
            inFile <- input$upload
            tryCatch({
                p <- readRDS(inFile$datapath)
                validObject(p)
                # Update species selector
                species <- as.character(p@species_params$species[foreground])
                updateSelectInput(session, "sp",
                                  choices = species,
                                  selected = species[1])
                
                # Update the reactive params object
                params(prepare_params(p))
                output$filename <- renderText(paste0("Previously uploaded file: ", 
                                                     inFile$name))
            },
            error = function(e) {
                showModal(modalDialog(
                    title = "Invalid parameter file",
                    HTML(paste0("Trying to load that file led to an error.<br>",
                                "The error message was:<br>", e)),
                    easyClose = TRUE
                ))
                p <- params()}
            )
        })
        
        ## Done ####
        observeEvent(input$done, {
            file.remove(logs)
            stopApp(params())
        })
        
        ## Prepare for download of params object ####
        output$params <- downloadHandler(
            filename = "params.rds", 
            content = function(file) {
                saveRDS(params(), file = file)
            })
        
        ## UI for side bar ####
        output$sp_sel <- renderUI({
            p <- isolate(params())
            species <- as.character(p@species_params$species[foreground])
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
            eff <- isolate(effort())
            
            l1 <- list(
                tags$h3(tags$a(id = "biomass"), "Biomass"),
                sliderInput("n0", "Egg density",
                            value = n0,
                            min = signif(n0 / 10, 3),
                            max = signif(n0 * 10, 3),
                            step = n0 / 50),
                # sliderInput("rescale", "Rescale all by", value = 1,
                #             min = 0.1,
                #             max = 5),
                tags$h3(tags$a(id = "predation"), "Predation"),
                sliderInput("gamma", "Predation rate coefficient gamma",
                            value = sp$gamma,
                            min = signif(sp$gamma / 2, 3),
                            max = signif(sp$gamma * 1.5, 3),
                            step = sp$gamma / 50, ticks = FALSE),
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
                            step = 0.05)
            )
            
            if (length(p@resource_dynamics) > 0) {
                l1 <- c(l1, list(tags$h3(tags$a(id = "resource"), "Resource encounter")))
                for (res in names(p@resource_dynamics)) {
                    res_var <- paste0("rho_", res)
                    l1 <- c(l1, list(
                        sliderInput(res_var, res,
                                    value = sp[[res_var]],
                                    min = signif(max(0, sp[[res_var]] - 0.1) / 2, 2),
                                    max = signif((sp[[res_var]] + 0.1) * 1.5, 2),
                                    step = 0.01)))
                }
            }
            
            l1 <- c(l1, list(tags$h3(tags$a(id = "fishing"), "Fishing"),
                             sliderInput("catchability", "Catchability",
                                         value = sp$catchability, 
                                         min = signif(max(0, sp$catchability / 2 - 1), 2), 
                                         max = signif(max(sp$catchability * 2, 2), 2), 
                                         step = 0.01)))
            
            if (sp$sel_func == "knife_edge") {
                l1 <- c(l1, list(
                             sliderInput("knife_edge_size", "knife_edge_size",
                                         value = sp$knife_edge_size, 
                                         min = 1, 
                                         max = signif(sp$knife_edge_size * 2, 2),
                                         step = 0.1)))
            } else if (sp$sel_func == "sigmoid_length") {
                l1 <- c(l1, list(
                             sliderInput("l50", "L50",
                                         value = sp$l50, 
                                         min = 1, 
                                         max = signif(sp$l50 * 2, 2),
                                         step = 0.1),
                             sliderInput("ldiff", "L50-L25",
                                         value = sp$l50 - sp$l25, 
                                         min = 0.1, 
                                         max = signif(sp$l50 / 10, 2),
                                         step = 0.1)))
            } else if (sp$sel_func == "double_sigmoid_length") {
                l1 <- c(l1, list(
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
                    sliderInput("l50_right", "L50 right",
                                value = sp$l50_right, 
                                min = 1, 
                                max = signif(sp$l50_right * 2, 2),
                                step = 0.1),
                    sliderInput("ldiff_right", "L50-L25 right",
                                value = sp$l25_right - sp$l50_right, 
                                min = 0.1, 
                                max = signif(sp$l50_right / 10, 2),
                                step = 0.1)
                ))
            }
            l1 <- c(l1, list(tags$h3(tags$a(id = "maturity"), "Maturity"),
                             sliderInput("w_mat", "w_mat", value = sp$w_mat,
                                         min = signif(sp$w_mat / 2, 2),
                                         max = signif(sp$w_mat * 1.5, 2)),
                             sliderInput("wfrac", "w_mat25/w_mat", value = sp$w_mat25/sp$w_mat,
                                         min = 0.5,
                                         max = 1,
                                         step = 0.01),
                             sliderInput("w_inf", "w_inf", value = sp$w_inf,
                                         min = signif(sp$w_inf / 2, 2),
                                         max = signif(sp$w_inf * 1.5, 2)),
                             sliderInput("m", "m", value = sp$m,
                                         min = 0,
                                         max = 2,
                                         step = 0.01),
                             tags$h3(tags$a(id = "others"), "Others"),
                             sliderInput("ks", "Coefficient of standard metabolism ks",
                                         value = sp$ks,
                                         min = signif(sp$ks / 2, 2),
                                         max = signif((sp$ks + 0.1) * 1.5, 2),
                                         step = 0.05),
                             sliderInput("k", "Coefficient of activity k",
                                         value = sp$k,
                                         min = signif(sp$k / 2, 2),
                                         max = signif((sp$k + 0.1) * 1.5, 2),
                                         step = 0.01),
                             # sliderInput("z0", "Mortality",
                             #             value = sp$z0,
                             #             min = signif(sp$z0 / 2, 2),
                             #             max = signif((sp$z0 + 0.1) * 1.5, 2),
                             #             step = 0.05),
                             sliderInput("alpha", "Assimilation efficiency alpha",
                                         value = sp$alpha,
                                         min = 0,
                                         max = 1),
                             tags$h3(tags$a(id = "interactions"), "Prey interactions"),
                             sliderInput("interaction_p", "Plankton",
                                         value = sp$interaction_p,
                                         min = 0,
                                         max = 1,
                                         step = 0.05)
            ))
            for (i in p@species_params$species) {
                inter_var <- paste0("inter_", i)
                l1 <- c(l1, list(
                    sliderInput(inter_var, i,
                                value = p@interaction[input$sp, i],
                                min = 0,
                                max = 1,
                                step = 0.05)
                ))
            }
            l1
        })
        
        output$general_params <- renderUI({
            
            p <- isolate(params())
            eff <- isolate(effort())
            
            l1 <- list(
                tags$h3(tags$a(id = "general"), "General"),
                sliderInput("effort", "Effort",
                            value = eff, 
                            min = 0,
                            max = signif(eff * 1.5, 2)),
                numericInput("p", "Exponent of metabolism",
                             value = p@p,
                             min = 0.6, max = 0.8, step = 0.005),
                numericInput("q", "Exponent of search volume",
                             value = p@q,
                             min = 0.6, max = 0.8, step = 0.005),
                numericInput("n", "Exponent of feeding rate",
                             value = p@n,
                             min = 0.6, max = 0.8, step = 0.005),
                tags$h3(tags$a(id = "plankton"), "Plankton"),
                numericInput("lambda", "Sheldon exponent lambda",
                             value = p@lambda, min = 1.9, max = 2.2, step = 0.005),
                numericInput("kappa", "Plankton coefficient kappa",
                             value = p@kappa),
                sliderInput("log_r_pp", "log10 Plankton replenishment rate",
                            value = 1, min = -1, max = 4, step = 0.05),
                numericInput("w_pp_cutoff", "Largest plankton",
                             value = p@w_full[which.min(p@cc_pp > 0)],
                             min = 1e-10,
                             max = 1e3),
                tags$h3(tags$a(id = "file"), "File management"),
                textOutput("filename"),
                fileInput("upload", "Upload new params", 
                          accept = ".rds"),
                downloadButton("params", "Download params")
            )
            l1
        })
        
        
        ## Adjust growth exponent ####
        observe({
            req(input$n)
            p <- isolate(params())
            sp <- isolate(input$sp)
            
            # change all h so that max intake rate at maturity stays the same
            p@species_params$h <- p@species_params$h * 
                p@species_params$w_mat^(p@n - input$n)
            h <- p@species_params[sp, "h"]
            updateSliderInput(session, "h",
                              value = h,
                              min = signif(h / 2, 2),
                              max = signif(h * 1.5, 2))
            
            # change all rho so that encounter rate at maturity stays the same
            for (res in names(p@resource_dynamics)) {
                res_var <- paste0("rho_", res)
                p@species_params[, res_var] <- p@species_params[, res_var] * 
                    p@species_params$w_mat^(p@n - input$n)
                new <- p@species_params[sp, res_var]
                updateSliderInput(session, res_var,
                                  value = new,
                                  min = signif(max(0, new - 0.1) / 2, 2),
                                  max = signif((new + 0.1) * 1.5, 2))
            }
            
            p <- setIntakeMax(p, n = input$n)
            p <- setResourceEncounter(p, n = input$n)
            params(p)
        })
        
        ## Adjust metabolism exponent ####
        observe({
            req(input$p)
            p <- isolate(params())
            sp <- isolate(input$sp)
            # change all ks so that metabolic rate at maturity stays the same
            p@species_params$ks <- p@species_params$ks * 
                p@species_params$w_mat^(p@p - input$p)
            p <- setMetab(p, p = input$p)
            params(p)
            ks <- p@species_params[sp, "ks"]
            updateSliderInput(session, "ks",
                              value = ks,
                              min = signif(ks / 2, 2),
                              max = signif((ks + 0.1) * 1.5, 2))
        })
        
        ## Rescale abundance ####
        # This is still buggy
        observe({
            req(input$rescale)
            p <- isolate(params())
            n0 <- isolate(input$n0)
            sp <- isolate(input$sp)
            # We want to rescale all abundances by the same factor
            p@initial_n <- p@initial_n * input$rescale
            p@initial_n_pp <- p@initial_n_pp * input$rescale
            p@initial_B <- p@initial_B * input$rescale
            p@cc_pp <- p@cc_pp * input$rescale
            p@kappa <- p@kappa * input$rescale
            n0 <- n0 * input$rescale
            # To keep the same per-capity behaviour, we have to scale down the
            # search volume
            p@species_params$gamma <- p@species_params$gamma / input$rescale
            p@search_vol <- p@search_vol / input$rescale
            params(p)
            
            gamma <- p@species_params[sp, "gamma"]
            updateSliderInput(session, "gamma",
                              value = gamma,
                              min = signif(gamma / 2, 3),
                              max = signif(gamma * 1.5, 3),
                              step = gamma / 50)
            updateSliderInput(session, "n0",
                              value = n0,
                              min = signif(n0 / 10, 3),
                              max = signif(n0 * 10, 3))
            updateSliderInput(session, "rescale", value = 1)
        })
        
        ## Adjust effort ####
        observe({
            req(input$effort)
            effort(input$effort)
            updateSliderInput(session, "effort",
                              max = signif(effort() * 1.5, 2))
        })
        
        ## Adjust von Bertalanffy ####
        observe({
            p <- isolate(params())
            p@species_params[isolate(input$sp), "k_vb"] <- req(input$k_vb)
            p@species_params[isolate(input$sp), "t0"] <- req(input$t0)
            p@species_params[isolate(input$sp), "a"] <- req(input$a)
            p@species_params[isolate(input$sp), "b"] <- req(input$b)
            params(p)
        })
        
        ## Adjust biomass observed ####
        observe({
            p <- isolate(params())
            p@species_params[isolate(input$sp), "biomass_observed"] <- 
                req(input$biomass_observed)
            p@species_params[isolate(input$sp), "abundance_observed"] <- 
                req(input$abundance_observed)
            p@species_params[isolate(input$sp), "cutoff_size"] <- 
                req(input$cutoff_size)
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
            
            if (sp != sp_old_n0) {
                # We came here after changing the species selector. This automatically
                # rewrote all sliders, so no need to update them. Also we do not want
                # update the params object.
                sp_old_n0 <<- sp
            } else {
                updateSliderInput(session, "n0",
                                  min = signif(n0 / 10, 3),
                                  max = signif(n0 * 10, 3))
                # rescale abundance to new egg density
                p@initial_n[sp, ] <- p@initial_n[sp, ] * n0 / 
                    p@initial_n[sp, p@w_min_idx[sp]]
                
                # Update the reactive params object
                params(p)
            }
        })
        
        
        # The "Tune egg density" button calculates the ratio of observed and
        # model biomass and then multiplies the egg density by that ratio. It
        # then runs the system to steady state.
        observeEvent(input$tune_egg, {
            p <- isolate(params())
            sp <- which.max(p@species_params$species == isolate(input$sp))
            if ("biomass_observed" %in% names(p@species_params) &&
                !is.na(p@species_params$biomass_observed[sp]) &&
                p@species_params$biomass_observed[sp] > 0) {
                total <- sum(p@initial_n[sp, ] * p@w * p@dw)
                n0 <- isolate(input$n0) * 
                    p@species_params$biomass_observed[sp] / total
                # rescale abundance to new egg density
                p@initial_n[sp, ] <- p@initial_n[sp, ] * n0 / 
                    p@initial_n[sp, p@w_min_idx[sp]]
                
                updateSliderInput(session, "n0",
                                  value = n0,
                                  min = signif(n0 / 10, 3),
                                  max = signif(n0 * 10, 3))
            }
        })
        
        ## Adjust plankton ####
        observe({
            req(input$kappa,
                input$lambda,
                input$log_r_pp,
                input$w_pp_cutoff)
            p <- isolate(params())
            p <- setPlankton(p, 
                             kappa = input$kappa, 
                             lambda = input$lambda,
                             r_pp = 10^input$log_r_pp,
                             w_pp_cutoff = input$w_pp_cutoff)
            params(p)
        })
        
        ## Adjust predation kernel ####
        observe({
            req(input$beta, input$sigma)
            p <- isolate(params())
            sp <- isolate(input$sp)
            
            if (sp != sp_old_kernel) {
                # We came here after changing the species selector. This automatically
                # rewrote all sliders, so no need to update them. Also we do not want
                # update the params object.
                sp_old_kernel <<- sp
            } else {
                # Update slider min/max so that they are a fixed proportion of the 
                # parameter value
                updateSliderInput(session, "beta",
                                  min = signif(input$beta / 2, 2),
                                  max = signif(input$beta * 1.5, 2))
                updateSliderInput(session, "sigma",
                                  min = signif(input$sigma / 2, 2),
                                  max = signif(input$sigma * 1.5, 2))
                p@species_params[sp, "beta"]  <- input$beta
                p@species_params[sp, "sigma"] <- input$sigma
                p <- setPredKernel(p)
                update_species(sp, p)
            }
        })
        
        ## Adjust predation ####
        observe({
            req(input$gamma, input$h, input$q)
            p <- isolate(params())
            sp <- isolate(input$sp)
            
            if (sp != sp_old_predation) {
                # We came here after changing the species selector. This automatically
                # rewrote all sliders, so no need to update them. Also we do not want
                # update the params object.
                sp_old_predation <<- sp
            } else {
                # Update slider min/max so that they are a fixed proportion of the 
                # parameter value
                updateSliderInput(session, "gamma",
                                  min = signif(input$gamma / 2, 3),
                                  max = signif(input$gamma * 1.5, 3))
                updateSliderInput(session, "h",
                                  min = signif(input$h / 2, 2),
                                  max = signif(input$h * 1.5, 2))
                p@species_params[sp, "gamma"] <- input$gamma
                p@species_params[sp, "h"]     <- input$h
                p@q <- input$q
                p <- setSearchVolume(p)
                p <- setIntakeMax(p)
                update_species(sp, p)
            }
        })
        
        ## Adjust fishing ####
        observe({
            req(input$catchability)
            p <- isolate(params())
            sp <- isolate(input$sp)
            
            if (sp != sp_old_fishing) {
                sp_old_fishing <<- sp
            } else {
                # Update slider min/max so that they are a fixed proportion of the 
                # parameter value
                updateSliderInput(session, "catchability",
                                  min = signif(max(input$catchability / 2 - 1, 0), 2),
                                  max = signif(max(input$catchability * 2, 2), 2))
                p@species_params[sp, "catchability"]  <- input$catchability
                
                if (p@species_params[sp, "sel_func"] == "knife_edge") {
                    updateSliderInput(session, "knife_edge_size",
                                      max = signif(input$knife_edge_size * 2, 2))
                    p@species_params[sp, "knife_edge_size"]   <- input$knife_edge_size
                }
                if (p@species_params[sp, "sel_func"] == "sigmoid_length" ||
                    p@species_params[sp, "sel_func"] == "double_sigmoid_length") {
                    updateSliderInput(session, "l50",
                                      max = signif(input$l50 * 2, 2))
                    updateSliderInput(session, "ldiff",  
                                      max = signif(input$l50 / 10, 2))
                    p@species_params[sp, "l50"]   <- input$l50
                    p@species_params[sp, "l25"]   <- input$l50 - input$ldiff
                }
                if (p@species_params[sp, "sel_func"] == "double_sigmoid_length") {
                    p@species_params[sp, "l50_right"]   <- input$l50_right
                    p@species_params[sp, "l25_right"]   <- input$l50_right + input$ldiff_right
                    updateSliderInput(session, "l50_right",
                                      max = signif(input$l50_right * 2, 2))
                    updateSliderInput(session, "ldiff_right",  
                                      max = signif(input$l50_right / 10, 2))
                }
                
                p <- setFishing(p)
                update_species(sp, p)
            }
        })
        
        # The Catch Tune button calculates the ratio of observed and
        # model catch and then multiplies the catchability by that ratio. It
        # then runs the system to steady state.
        observeEvent(input$tune_catch, {
            p <- isolate(params())
            sp <- which.max(p@species_params$species == isolate(input$sp))
            if ("catch_observed" %in% names(p@species_params) &&
                !is.na(p@species_params$catch_observed[sp]) &&
                p@species_params$catch_observed[sp] > 0) {
                total <- sum(p@initial_n[sp, ] * p@w * p@dw *
                                 getFMort(p, effort = effort())[sp, ])
                p@species_params$catchability[sp] <- 
                    p@species_params$catchability[sp] *
                    p@species_params$catch_observed[sp] / total
                updateSliderInput(session, "catchability",
                                  value = p@species_params$catchability[sp],
                                  min = signif(max(0, p@species_params$catchability[sp] / 2 - 1), 2),
                                  max = signif(max(p@species_params$catchability[sp] * 2, 2), 2))
                p <- setFishing(p)
                tryCatch({
                    # Create a Progress object
                    progress <- shiny::Progress$new(session)
                    on.exit(progress$close())
                    
                    # Run to steady state
                    p <- steady(p, effort = effort(), t_max = 100, tol = 1e-2,
                                progress_bar = progress)
                    
                    # Update the reactive params object
                    params(p)
                    add_to_logs(p)
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
                params(p)
            }
        })
        
        ## Adjust maturity ####
        observe({
            req(input$w_mat, input$wfrac, input$w_inf, input$m)
            p <- isolate(params())
            sp <- isolate(input$sp)
            
            if (sp != sp_old_maturity) {
                sp_old_maturity <<- sp
            } else {
                # Update slider min/max so that they are a fixed proportion of the 
                # parameter value
                updateSliderInput(session, "w_mat",
                                  min = signif(input$w_mat / 2, 2),
                                  max = signif(input$w_mat * 1.5, 2))
                updateSliderInput(session, "w_inf",
                                  min = signif(input$w_inf / 2, 2),
                                  max = signif(input$w_inf * 1.5, 2))
                
                p@species_params[sp, "w_mat25"]   <- input$w_mat * input$wfrac
                p@species_params[sp, "w_mat"]   <- input$w_mat
                p@species_params[sp, "w_inf"]   <- input$w_inf
                p@species_params[sp, "m"]     <- input$m
                
                p <- setReproduction(p)
                update_species(sp, p)
            }
        })
        
        update_species <- function(sp, p) {
            # wrap the code in trycatch so that when there is a problem we can
            # simply stay with the old parameters
            tryCatch({
                # The spectrum for the changed species is calculated with new
                # parameters but in the context of the original community
                # Compute death rate for changed species
                mumu <- getMort(p, effort = effort())[sp, ]
                # compute growth rate for changed species
                gg <- getEGrowth(p)[sp, ]
                # Compute solution for changed species
                w_inf_idx <- sum(p@w < p@species_params[sp, "w_inf"])
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
                p@initial_n[sp, ] <- 0
                p@initial_n[sp, p@w_min_idx[sp]:w_inf_idx] <-
                    c(1, cumprod(gg[idx] / ((gg + mumu * p@dw)[idx + 1]))) *
                    n0
                if (any(is.infinite(p@initial_n))) {
                    stop("Candidate steady state holds infinities")
                }
                if (any(is.na(p@initial_n) || is.nan(p@initial_n))) {
                    stop("Candidate steady state holds non-numeric values")
                }
                
                # Retune the value of erepro so that we get the correct level of
                # recruitment
                i <- which(p@species_params$species == sp)
                rdd <- getRDD(p, p@initial_n, p@initial_n_pp)[i]
                gg0 <- gg[p@w_min_idx[i]]
                mumu0 <- mumu[p@w_min_idx[i]]
                DW <- p@dw[p@w_min_idx[i]]
                p@species_params$erepro[i] <- p@species_params$erepro[i] *
                    n0 * (gg0 + DW * mumu0) / rdd
                
                # Update the reactive params object
                params(p)
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
        
        ## Adjust species parameters ####
        observe({
            req(input$interaction_p,
                input$alpha, input$ks, input$k)
            p <- isolate(params())
            sp <- isolate(input$sp)
            
            # Create updated species params data frame
            species_params <- p@species_params
            species_params[sp, "interaction_p"] <- input$interaction_p
            species_params[sp, "alpha"] <- input$alpha
            species_params[sp, "ks"]    <- input$ks
            species_params[sp, "k"]     <- input$k
            #species_params[sp, "z0"]     <- input$z0
            if (length(p@resource_dynamics) > 0) {
                for (res in names(p@resource_dynamics)) {
                    res_var <- paste0("rho_", res)
                    species_params[sp, res_var] <- input[[res_var]]
                }
            }
            for (i in p@species_params$species) {
                inter_var <- paste0("inter_", i)
                p@interaction[sp, i] <- input[[inter_var]]
            }
            
            if (sp != sp_old) {
                # We came here after changing the species selector. This automatically
                # rewrote all sliders, so no need to update them. Also we do not want
                # update the params object.
                sp_old <<- sp
            } else {
                # Update slider min/max so that they are a fixed proportion of the 
                # parameter value
                updateSliderInput(session, "ks",
                                  min = signif(input$ks / 2, 2),
                                  max = signif((input$ks + 0.1) * 1.5, 2))
                updateSliderInput(session, "k",
                                  min = signif(input$k / 2, 2),
                                  max = signif((input$k + 0.1) * 1.5, 2))
                # updateSliderInput(session, "z0",
                #                   min = signif(input$z0 / 2, 2),
                #                   max = signif((input$z0 + 0.1) * 1.5, 2))
                
                if (length(p@resource_dynamics) > 0) {
                    for (res in names(p@resource_dynamics)) {
                        res_var <- paste0("rho_", res)
                        updateSliderInput(session, res_var,
                                          min = signif(max(0, input[[res_var]] - 0.1) / 2, 2),
                                          max = signif((input[[res_var]] + 0.1) * 1.5, 2))
                    }
                }
                # Create new params object identical to old one except for changed
                # species params
                p@species_params <- species_params
                p <- setInteraction(p)
                p <- setMetab(p)
                #p <- setBMort(p)
                p <- setReproduction(p)
                p <- setFishing(p)
                p <- setResourceEncounter(p)
                update_species(sp, p)
            }
        })
        
        ## Undo ####
        observeEvent(input$undo, {
            if (log_idx <= 1) return()
            p_new <- readRDS(logs[log_idx])
            p_old <- params()
            # if the params have not changed, go to the previous one
            if (all(p_old@species_params == p_new@species_params, na.rm = TRUE)) {
                log_idx <<- log_idx - 1
                shinyjs::enable("redo")
                p_new <- readRDS(logs[log_idx])
                if (log_idx == 1) {
                    shinyjs::disable("undo")
                    shinyjs::disable("undo_all")
                }
            }
            params(p_new)
            # Trigger an update of sliders
            trigger_update(runif(1))
        })
        ## Redo ####
        observeEvent(input$redo, {
            if (log_idx >= length(logs)) return()
            log_idx <<- log_idx + 1
            params(readRDS(logs[log_idx]))
            # Trigger an update of sliders
            trigger_update(runif(1))
            shinyjs::enable("undo")
            shinyjs::enable("undo_all")
            if (log_idx == length(logs)) shinyjs::disable("redo")
        })
        ## Cancel ####
        observeEvent(input$undo_all, {
            if (log_idx > 1) shinyjs::enable("redo")
            shinyjs::disable("undo")
            shinyjs::disable("undo_all")
            log_idx <- 1
            params(readRDS(logs[log_idx]))
            # Trigger an update of sliders
            trigger_update(runif(1))
        })
        
        observeEvent(input$growth_click, {
            if (input$growth_click$panelvar1 != input$sp) {
                updateSelectInput(session, "sp",
                                  selected = input$growth_click$panelvar1)
                updateRadioButtons(session, "all_growth",
                                   selected = "Selected species")
            }
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
                p <- steady(p, effort = effort(), t_max = 100, tol = 1e-2,
                            progress_bar = progress)
                
                # Update the reactive params object
                params(p)
                add_to_logs(p)
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
            if (!input$all_growth == "All") {
                k_vb <- params()@species_params[input$sp, "k_vb"]
                t0 <- params()@species_params[input$sp, "t0"]
                a <- params()@species_params[input$sp, "a"]
                b <- params()@species_params[input$sp, "b"]
                list(
                    div(style = "display:inline-block",
                        numericInput("k_vb", "Von Bertalanffy k", value = k_vb, width = "9em")),
                    div(style = "display:inline-block",
                        numericInput("t0", "t_0", value = t0, width = "6em")),
                    p("Parameters for length-weight relationship l = a w^b"),
                    div(style = "display:inline-block",
                        numericInput("a", "a", value = a, width = "8em")),
                    div(style = "display:inline-block",
                        numericInput("b", "b", value = b, width = "8em"))
                )
            }
        })
        output$plotGrowthCurve <- renderPlot({
            p <- params()
            if (input$all_growth == "All") {
                gc <- getGrowthCurves(p)[foreground, ] %>% 
                    as.tbl_cube(met_name = "Size") %>% 
                    as_tibble() %>%
                    mutate(Legend = "Model")
                
                vb <- gc %>% 
                    mutate(Legend = "von Bertalanffy") %>% 
                    mutate(a = p@species_params[Species, "a"],
                           b = p@species_params[Species, "b"],
                           k_vb = p@species_params[Species, "k_vb"],
                           t0 = p@species_params[Species, "t0"],
                           L_inf = (p@species_params[Species, "w_inf"] / a)^(1 / b),
                           Size = a * (L_inf * (1 - exp(-k_vb * (Age - t0))))^b) %>% 
                    select(names(gc))
                
                ggplot(bind_rows(gc, vb)) +
                    geom_line(aes(x = Age, y = Size, colour = Legend)) +
                    scale_x_continuous(name = "Age [years]") +
                    scale_y_continuous(name = "Size [g]") +
                    geom_hline(aes(yintercept = w_mat), 
                               data = tibble(Species = p@species_params$species[foreground],
                                             w_mat = p@species_params$w_mat[foreground]), 
                               linetype = "dashed",
                               colour = "grey") +
                    facet_wrap(~Species, scales = "free_y")
            } else {
                plotGrowthCurves(p, species = input$sp) +
                    theme_grey(base_size = base_size)
            }
        })
        
        ## Spectra ####
        output$plotSpectra <- renderPlotly({
            if (input$binning == "Logarithmic") {
                power <- 2
            } else {
                power <- 1
            }
            plotSpectra(params(), power = power, highlight = input$sp, total = TRUE) + 
                theme_grey(base_size = base_size)
        })
        
        ## erepro plot ####
        output$plot_erepro <- renderPlotly({
            p <- params()
            df <- data.frame(Species = factor(p@species_params$species[foreground],
                                              levels = p@species_params$species),
                             erepro = p@species_params$erepro[foreground])
            ggplot(df, aes(x = Species, y = erepro)) + 
                geom_col() + geom_hline(yintercept = 1, color = "red") +
                scale_y_log10() +
                theme_grey(base_size = base_size) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        })
        
        ## Plot feeding level ####
        output$plot_feeding_level <- renderPlotly({
            if (input$all_growth == "All") {
                plotFeedingLevel(params(), highlight = input$sp) + 
                    theme_grey(base_size = base_size)
            } else {
                plotFeedingLevel(params(), species = input$sp) + 
                    theme_grey(base_size = base_size)
            }
        })
        
        ## Biomass plot ####
        output$biomass_sel <- renderUI({
            sp <- input$sp
            p <- isolate(params())
            species_params <- p@species_params[sp, ]
            if (is.null(species_params$biomass_observed) ||
                is.na(species_params$biomass_observed)) {
                species_params$biomass_observed <- 0
            }
            if (is.null(species_params$abundance_observed) ||
                is.na(species_params$abundance_observed)) {
                species_params$abundance_observed <- 0
            }
            if (is.null(species_params$cutoff_size) ||
                is.na(species_params$cutoff_size)) {
                species_params$cutoff_size <- 0
            }
            list(
                div(style = "display:inline-block",
                    numericInput("biomass_observed", 
                                 paste0("Observed biomass for ", sp),
                                 value = species_params$biomass_observed)),
                div(style = "display:inline-block",
                    numericInput("abundance_observed", 
                                 paste0("Observed abundance for ", sp),
                                 value = species_params$abundance_observed)),
                div(style = "display:inline-block",
                    numericInput("cutoff_size", "Lower cutoff",
                                 value = species_params$cutoff_size))
            )
        })
        output$plotBiomassDist <- renderPlotly({
            req(input$sp, input$cutoff_size, input$biomass_observed)
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
                theme_grey(base_size = base_size) +
                labs(x = "Size [g]", y = "Cummulative biomass")  +
                geom_text(aes(x = p@species_params[sp, "w_mat"], 
                              y = max(Biomass * 0.2),
                              label = "\nMaturity"), 
                          angle = 90)
            if (input$biomass_observed) {
                cutoff_idx <- which.max(p@w >= input$cutoff_size)
                target <- input$biomass_observed + biomass[cutoff_idx]
                pl <- pl +
                    geom_hline(yintercept = biomass[cutoff_idx]) +
                    geom_vline(xintercept = input$cutoff_size) +
                    geom_hline(yintercept = target, color = "green")
            }
            pl
        })
        output$plotTotalBiomass <- renderPlotly({
            p <- params()
            cutoff <- p@species_params$cutoff_size
            # When no cutoff known, set it to maturity weight / 20
            if (is.null(cutoff)) cutoff <- p@species_params$w_mat / 20
            cutoff[is.na(cutoff)] <- p@species_params$w_mat[is.na(cutoff)] / 20
            observed <- p@species_params$biomass_observed
            if (is.null(observed)) observed <- 0
            
            biomass_model <- foreground_indices  # create vector of right length
            for (i in seq_along(foreground_indices)) {
                sp <- foreground_indices[i]
                cum_biomass <- cumsum(p@initial_n[sp, ] * p@w * p@dw)
                cutoff_idx <- which.max(p@w >= cutoff[sp])
                biomass_model[i] <- max(cum_biomass) - cum_biomass[cutoff_idx]
            }
            species <- factor(p@species_params$species[foreground],
                              levels = p@species_params$species[foreground])
            df <- rbind(
                data.frame(Species = species,
                           Type = "Observed",
                           Biomass = observed[foreground]),
                data.frame(Species = species,
                           Type = "Model",
                           Biomass = biomass_model)
            )
            ggplot(df) +
                geom_col(aes(x = Species, y = Biomass, fill = Type),
                         position = "dodge") +
                scale_y_continuous(name = "Biomass [g]", trans = "log10",
                                   breaks = log_breaks()) +
                theme_grey(base_size = base_size) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        })
        output$plotTotalAbundance <- renderPlotly({
            p <- params()
            cutoff <- p@species_params$cutoff_size
            # When no cutoff known, set it to maturity weight / 20
            if (is.null(cutoff)) cutoff <- p@species_params$w_mat / 20
            cutoff[is.na(cutoff)] <- p@species_params$w_mat[is.na(cutoff)] / 20
            observed <- p@species_params$abundance_observed
            if (is.null(observed)) observed <- 0
            
            abundance_model <- foreground_indices  # create vector of right length
            for (i in seq_along(foreground_indices)) {
                sp <- foreground_indices[i]
                cum_abundance <- cumsum(p@initial_n[sp, ] * p@dw)
                cutoff_idx <- which.max(p@w >= cutoff[sp])
                abundance_model[i] <- max(cum_abundance) - cum_abundance[cutoff_idx]
            }
            species <- factor(p@species_params$species[foreground],
                              levels = p@species_params$species[foreground])
            df <- rbind(
                data.frame(Species = species,
                           Type = "Observed",
                           Abundance = observed[foreground]),
                data.frame(Species = species,
                           Type = "Model",
                           Abundance = abundance_model)
            )
            ggplot(df) +
                geom_col(aes(x = Species, y = Abundance, fill = Type),
                         position = "dodge") +
                scale_y_continuous(name = "Abundance", trans = "log10",
                                   breaks = log_breaks()) +
                theme_grey(base_size = base_size) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        })
        
        ## Plot catch ####
        
        # Catch density for selected species
        output$plotCatchDist <- renderPlotly({
            req(input$sp)
            p <- params()
            sp <- which.max(p@species_params$species == input$sp)
            p <- set_species_param_default(p, "a", 0.006)
            p <- set_species_param_default(p, "b", 3)
            a <- p@species_params[sp, "a"]
            b <- p@species_params[sp, "b"]
            
            # Check whether we have enough catch data for this species to plot it
            is_observed <- sum(catch$species == input$sp) > 3
            
            # To choose the range of sizes over which to plot we look at the range
            # of sizes for which a non-zero catch was observed. If no catch was
            # observed for the species, we use the range from w_mat/100 to w_inf.
            if (is_observed) {
                if ("length" %in% names(catch)) {
                    l_min = min(catch$length[catch$species == input$sp])
                    w_min = a * l_min ^ b
                    l_max = max(catch$length[catch$species == input$sp])
                    w_max = a * l_max ^ b
                } else {
                    w_min = min(catch$weight[catch$species == input$sp])
                    w_max = max(catch$weight[catch$species == input$sp])
                }
                w_min_idx <- sum(p@w < w_min)
                w_max_idx <- sum(p@w <= w_max)
            } else {
                w_min_idx <- sum(p@w < (p@species_params$w_mat[sp] / 100))
                w_max_idx <- sum(p@w <= p@species_params$w_inf[sp])
            }
            w_sel <- seq(w_min_idx, w_max_idx, by = 1)
            w <- p@w[w_sel]
            l = (p@w[w_sel] / a) ^ (1 / b)
            
            catch_w <- getFMort(p, effort = effort())[sp, w_sel] * 
                p@initial_n[sp, w_sel]
            # We just want the distribution, so we rescale the density so its area is 1
            if (sum(catch_w) > 0) catch_w <- catch_w / sum(catch_w * p@dw[w_sel])
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
            abundance <- data.frame(w, l, catch_w, catch_l, type = "Abundance")
            
            if (is_observed) {
                sel <- (catch$species == input$sp)
                if ("length" %in% names(catch)) {
                    l <- catch$length[sel]
                    dl <- catch$dl[sel]
                    catch_l <- catch$catch[sel]
                    # normalise to a density in l
                    catch_l <- catch_l / sum(catch_l * dl)
                    # To get the density in w we need to divide by dw/dl
                    w <- a * l ^ b
                    catch_w <- catch_l / b * l / w
                } else {
                    w <- catch$weight[sel]
                    dw <- catch$dw[sel]
                    catch_w <- catch$catch[sel]
                    # normalise to a density in w
                    catch_w <- catch_w / sum(catch_w * dw)
                    # To get the density in l we need to divide by dl/dw
                    l <- (w / a)^(1/b)
                    catch_l <- catch_w * b / l * w
                }
                df <- rbind(df, data.frame(w, l, catch_w, catch_l, 
                                           type = "Observed catch"))
            }
            # From the abundance only keep values that are no larger than
            # the maximum of the other shown densities.
            if (input$catch_x == "Weight") {
                abundance <- subset(abundance, catch_w < max(df$catch_w))
            } else {
                abundance <- subset(abundance, catch_l < max(df$catch_l))
            }
            # Add the abundance to the data frame last so that it shows up
            # last also in legend
            df <- rbind(df, abundance)
            
            if (input$catch_x == "Weight") {
                mat  <- p@species_params$w_mat[sp]
                pl <- ggplot(df) +
                    geom_line(aes(x = w, y = catch_w, color = type)) +
                    geom_text(aes(x = mat, y = max(catch_w * 0.9), label = "\nMaturity"), 
                              angle = 90) +
                    labs(x = "Size [g]", y = "Normalised number density [1/g]")
            } else {
                mat <- (p@species_params$w_mat[sp] / a) ^ (1 / b)
                pl <- ggplot(df) +
                    geom_line(aes(x = l, y = catch_l, color = type)) +
                    geom_text(aes(x = mat, y = max(catch_l * 0.9), label = "\nMaturity"), 
                              angle = 90) +
                    labs(x = "Size [cm]", y = "Normalised number density [1/cm]")
            }
            pl +
                geom_vline(xintercept = mat, linetype = "dotted")  +
                theme_grey(base_size = base_size) +
                scale_colour_manual(values = c("Model catch" = "blue",
                                               "Observed catch" = "red",
                                               "Abundance" = "grey"))
        })
        
        # Total catch by species
        output$plotTotalCatch <- renderPlotly({
            p <- params()
            if (is.null(p@species_params$catch_observed)) {
                p@species_params$catch_observed <- 0
            }
            biomass <- sweep(p@initial_n, 2, p@w * p@dw, "*")
            total <- rowSums(biomass * getFMort(p, effort = effort()))
            df <- rbind(
                data.frame(Species = p@species_params$species[foreground],
                           Type = "Model",
                           Catch = total[foreground]),
                data.frame(Species = p@species_params$species[foreground],
                           Type = "Observed",
                           Catch = p@species_params$catch_observed[foreground])
            )
            ggplot(df) +
                geom_col(aes(x = Species, y = Catch, fill = Type),
                         position = "dodge") +
                theme_grey(base_size = base_size) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
                labs(x = "", y = "Catch [megatonnes]") +
                scale_fill_manual(values = c("Model" = "blue",
                                             "Observed" = "red"))
        })
        
        # Input field for observed catch
        output$catch_sel <- renderUI({
            p <- isolate(params())
            sp <- input$sp
            numericInput("catch_observed", 
                         paste0("Observed total catch for ", sp, " (megatonnes)"),
                         value = p@species_params[sp, "catch_observed"])
        })
        
        # Output of model catch
        output$catch_total <- renderText({
            p <- params()
            sp <- which.max(p@species_params$species == input$sp)
            total <- sum(p@initial_n[sp, ] * p@w * p@dw *
                             getFMort(p, effort = effort())[sp, ])
            paste("Model catch:", total)
        })
        
        ## Plot rates ####  
        output$plotGrowth <- renderPlotly({
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
                theme_grey(base_size = base_size) +
                labs(x = "Size [g]", y = "Rate [g/year]")  +
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
        output$plotDeath <- renderPlotly({
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
                                  effort = effort())[sp,sel],
                          getPredMort(p, p@initial_n, p@initial_n_pp)[sp, sel],
                          getFMort(p, effort = effort())[sp, sel],
                          p@mu_b[sp,sel])
            )
            pl <- ggplot(df, aes(x = w, y = value, color = Type)) + 
                geom_line() + 
                geom_vline(xintercept = p@species_params[sp, "w_mat"], 
                           linetype = "dotted") +
                geom_vline(xintercept = p@species_params[sp, "w_inf"], 
                           linetype = "dotted") +
                theme_grey(base_size = base_size) +
                labs(x = "Size [g]", y = "Rate [1/year]")  +
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
        output$plot_prey <- renderPlotly({
            p <- params()
            sp <- which.max(p@species_params$species == input$sp)
            x <- log(p@w_full)
            dx <- x[2] - x[1]
            xp <- req(input$pred_size)
            wp <- exp(xp)
            wp_idx <- sum(p@w <= wp)
            # Calculate total community abundance
            # Todo: take interaction matrix into account
            fish_idx <- (length(p@w_full) - length(p@w) + 1):length(p@w_full)
            total_n <- p@initial_n_pp
            total_n[fish_idx] <- total_n[fish_idx] + 
                p@interaction[sp, ] %*% p@initial_n
            totalx <- total_n * p@w_full
            #totalx <- totalx / sum(totalx * dx)
            phix <- getPredKernel(p)[sp, wp_idx, ]
            pr <- totalx * phix
            br <- pr * p@w_full
            # convert to proportions
            phix <- phix / sum(phix * dx)
            pr <- pr / sum(pr * dx)
            br <- br / sum(br * dx)
            df <- tibble::tibble(w = p@w_full, 
                                 Kernel = phix, 
                                 Biomass = br, 
                                 Numbers = pr) %>%
                tidyr::gather(type, y, Kernel, Biomass, Numbers)
            ggplot(df) +
                geom_line(aes(w, y, color = type)) +
                labs(x = "Weight [g]", y = "Density") +
                geom_point(aes(x = wp, y = 0), size = 4, colour = "blue") +
                scale_x_log10()
        })
        
        ## Plot diet ####
        output$plot_diet <- renderPlotly({
            req(input$sp)
            plotDiet(params(), input$sp)
        })
        
        ## Plot predators ####
        output$plot_pred <- renderPlotly({
            req(input$sp)
            sp <- input$sp
            p <- params()
            species <- factor(p@species_params$species,
                              levels = p@species_params$species)
            fish_idx_full <- (p@w_full >= p@species_params[sp, "w_min"]) &
                (p@w_full <= p@species_params[sp, "w_inf"])
            fish_idx <- (p@w >= p@species_params[sp, "w_min"]) &
                (p@w <= p@species_params[sp, "w_inf"])
            pred_rate <- p@interaction[, sp] * 
                getPredRate(p, p@initial_n, p@initial_n_pp, p@initial_B)[, fish_idx_full]
            fishing <- getFMort(p, effort = effort())[sp, fish_idx]
            total <- colSums(pred_rate) + p@mu_b[sp, fish_idx] + fishing
            ylab <- "Death rate [1/year]"
            background <- p@mu_b[sp, fish_idx]
            if (input$death_prop == "Proportion") {
                pred_rate <- pred_rate / rep(total, each = dim(pred_rate)[[1]])
                background <- background / total
                fishing <- fishing / total
                ylab <- "Proportion of all death"
            }
            # Make data.frame for plot
            plot_dat <- 
                rbind(
                    data.frame(value = background,
                               Cause = "Background",
                               w = p@w[fish_idx]),
                    data.frame(value = fishing,
                               Cause = "Fishing",
                               w = p@w[fish_idx]),
                    data.frame(value = c(pred_rate),
                               Cause = species,
                               w = rep(p@w[fish_idx], each = dim(pred_rate)[[1]]))
                )
            ggplot(plot_dat) +
                geom_area(aes(x = w, y = value, fill = Cause)) +
                scale_x_log10() +
                labs(x = "Size [g]", y = ylab) +
                scale_fill_manual(values = p@linecolour)
        })
        
        ## Plot psi ####
        output$plot_psi <- renderPlotly({
            if (!input$all_growth == "All") {
                p <- params()
                sp <- which.max(p@species_params$species == input$sp)
                w_min <- p@species_params$w_inf[sp] / 50
                sel <- p@w >= w_min & p@w <= p@species_params$w_inf[sp]
                df <- data.frame(Size = p@w[sel], psi = p@psi[sp, sel])
                ggplot(df, aes(x = Size, y = psi)) + 
                    geom_line(color = "blue") +
                    geom_vline(xintercept = p@species_params[sp, "w_mat"], 
                               linetype = "dotted") +
                    theme_grey(base_size = base_size) +
                    labs(x = "Size [g]", y = "Proportion of energy for reproduction")  +
                    geom_text(aes(x = p@species_params[sp, "w_mat"], 
                                  y = max(psi * 0.8),
                                  label = "\nMaturity"), 
                              angle = 90)
            }
        })
        
        ## Plot plankton ####
        output$plot_plankton <- renderPlot({
            p <- params()
            select <- (p@cc_pp > 0)
            plot_dat <- data.frame(
                x = p@w_full[select],
                y = p@initial_n_pp[select] / p@cc_pp[select]
            )
            ggplot(plot_dat) +
                geom_line(aes(x, y)) +
                scale_x_log10("Plankton size [g]") +
                ylab("Proportion of carrying capacity") +
                theme_grey(base_size = 16)
        })
        
        ## Plot plankton predators ####
        output$plot_plankton_pred <- renderPlotly({
            p <- params()
            species <- factor(p@species_params$species,
                              levels = p@species_params$species)
            select <- (p@cc_pp > 0)
            pred_rate <- p@species_params$interaction_p * 
                getPredRate(p)[, select]
            total <- colSums(pred_rate)
            ylab <- "Death rate [1/year]"
            if (input$plankton_death_prop == "Proportion") {
                pred_rate <- pred_rate / rep(total, each = dim(pred_rate)[[1]])
                ylab = "Proportion of predation"
            }
            # Make data.frame for plot
            plot_dat <- data.frame(
                value = c(pred_rate),
                Predator = species,
                w = rep(p@w_full[select], each = dim(pred_rate)[[1]]))
            ggplot(plot_dat) +
                geom_area(aes(x = w, y = value, fill = Predator)) +
                scale_x_log10("Plankton size [g]") +
                ylab(ylab) +
                scale_fill_manual(values = p@linecolour) +
                theme_grey(base_size = base_size)
        })
        
        ## Plot stomach content ####
        output$plot_stomach <- renderPlot({
            req(input$sp)
            p <- params()
            sp <- which.max(p@species_params$species == input$sp)
            
            df <- tibble(
                x = x,
                Kernel = double_sigmoid(
                    x, 
                    p_l = input$p_l, 
                    s_l = input$s_l, 
                    p_r = input$p_r, 
                    s_r = input$s_r, 
                    ex = input$ex)) %>% 
                mutate(Numbers = Kernel / exp((1 + alpha - lambda) * x),
                       Biomass = Numbers / exp(x),
                       Kernel = Kernel / sum(Kernel) / dx,
                       Numbers = Numbers / sum(Numbers) / dx,
                       Biomass = Biomass / sum(Biomass) / dx) %>% 
                gather(key = Type, value = "Density",
                       Numbers, Biomass)
            
            pl + geom_line(data = df,
                           aes(x, Density, colour = Type),
                           size = 3) 
            st <- stomach %>% 
                filter(species == input$sp)
        })
        
    } #the server
    
    runGadget(ui, server, viewer = browserViewer())
}