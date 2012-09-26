# Class outline for sbm base parameters class
# Class has members to store parameters of size based model

# Copyright 2012 Finlay Scott and Julia Blanchard. 
# Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, CEFAS. finlay.scott@cefas.co.uk

#Naming conventions:
#
#S4 classes and constructors: AClass
#S4 methods: aMethod
#functions a_function


# Validity function - pretty long...
# Not documented as removed later on
valid_MizerParams <- function(object) {

    errors <- character()
    # grab some dims
    length_w <- length(object@w)
    length_w_full <- length(object@w_full)
    #no_species <- 

    # Check dw and dw_full are correct length
    if(length(object@dw) != length_w){
	msg <- paste("dw is length ", length(object@dw), " and w is length ", length_w, ". These should be the same length", sep="")
	errors <- c(errors, msg)
    }

    if(length(object@dw_full) != length_w_full){
	msg <- paste("dw_full is length ", length(object@dw_full), " and w_full is length ", length_w_full, ". These should be the same length", sep="")
	errors <- c(errors, msg)
    }
    # Check the array dimensions are good
    # 2D arrays
    if(!all(c(
	length(dim(object@psi)),
	length(dim(object@intake_max)),
	length(dim(object@search_vol)),
	length(dim(object@activity)),
	length(dim(object@std_metab)),
	length(dim(object@interaction)),
	length(dim(object@catchability))) == 2)){
	    msg <- "psi, intake_max, search_vol, activity, std_metab, interaction and catchability must all be two dimensions"
	    errors <- c(errors, msg)
    }
    # 3D arrays
    if(!all(c(
	length(dim(object@pred_kernel)),
	length(dim(object@selectivity))) == 3)){
	    msg <- "pred_kernel and selectivity must be three dimensions"
	    errors <- c(errors, msg)
    }
    # Check number of species is equal across relevant slots
    if(!all(c(
	dim(object@psi)[1],
	dim(object@intake_max)[1],
	dim(object@search_vol)[1],
	dim(object@std_metab)[1],
	dim(object@activity)[1],
	dim(object@pred_kernel)[1],
	dim(object@selectivity)[2],
	dim(object@catchability)[2],
	dim(object@interaction)[1],
	dim(object@interaction)[2]) == 
	    dim(object@species_params)[1])){
	    msg <- "The number of species in the model must be consistent across the species_params, psi, intake_max, search_vol, activity, pred_kernel, interaction (dim 1), selectivity, catchability and interaction (dim 2) slots"
	    errors <- c(errors, msg)
    }
    # Check number of size groups
    if(!all(c(
	dim(object@psi)[2],
	dim(object@intake_max)[2],
	dim(object@search_vol)[2],
	dim(object@activity)[2],
	dim(object@std_metab)[2],
	dim(object@pred_kernel)[2],
	dim(object@selectivity)[3]) ==
	    length_w)){
	    msg <- "The number of size bins in the model must be consistent across the w, psi, intake_max, search_vol, activity, pred_kernel (dim 2) and selectivity (dim 3) slots"
	    errors <- c(errors, msg)
    }
    # Check number of full spectrum size groups
    if(!isTRUE(all.equal(dim(object@pred_kernel)[3],length(object@w_full)))){
	msg <- "The length of the full size spectrum in the third dimension of the pred_kernel slot must be the same length as the w_full slot"
	errors <- c(errors, msg)
    }
    # Check numbe of gears
    if(!isTRUE(all.equal(dim(object@selectivity)[1], dim(object@catchability)[1]))){
	msg <- "The number of fishing gears must be consistent across the catchability and selectivity (dim 1) slots"
	errors <- c(errors, msg)
    }
    # Check names of dimnames of arrays
    # sp dimension
    if(!all(c(
	names(dimnames(object@psi))[1],
	names(dimnames(object@intake_max))[1],
	names(dimnames(object@search_vol))[1],
	names(dimnames(object@activity))[1],
	names(dimnames(object@std_metab))[1],
	names(dimnames(object@pred_kernel))[1],
	names(dimnames(object@selectivity))[2],
	names(dimnames(object@catchability))[2]) == "sp")){
	    msg <- "Name of first dimension of psi, intake_max, search_vol, std_metab, activity and pred_kernel and the second dimension of selectivity and catchability must be 'sp'"
	    errors <- c(errors, msg)
	}
    #interaction dimension names
    if(names(dimnames(object@interaction))[1] != "predator"){
	msg <- "The first dimension of interaction must be called 'predator'"
	errors <- c(errors, msg)
    }
    if(names(dimnames(object@interaction))[2] != "prey"){
	msg <- "The first dimension of interaction must be called 'prey'"
	errors <- c(errors, msg)
    }
    # w dimension
    if(!all(c(
	names(dimnames(object@psi))[2],
	names(dimnames(object@intake_max))[2],
	names(dimnames(object@search_vol))[2],
	names(dimnames(object@std_metab))[2],
	names(dimnames(object@activity))[2],
	names(dimnames(object@selectivity))[3]) == "w")){
	    msg <- "Name of second dimension of psi, intake_max, search_vol, std_metab, activity and third dimension of selectivity must be 'w'"
	    errors <- c(errors, msg)
	}
    if(names(dimnames(object@pred_kernel))[2] != "w_pred"){
	msg <- "Name of second dimension of pred_kernel must be 'w_pred'"
	errors <- c(errors, msg)
    }
    if(names(dimnames(object@pred_kernel))[3] != "w_prey"){
	msg <- "Name of third dimension of pred_kernel must be 'w_prey'"
	errors <- c(errors, msg)
    }
    if(!all(c(
	  names(dimnames(object@selectivity))[1],
	  names(dimnames(object@catchability))[1]) == "gear")){
	msg <- "Name of first dimension of selectivity and catchability must be 'gear'"
	errors <- c(errors, msg)
    }

    # Check dimnames of species are identical
    # Bit tricky this one as I don't know of a way to compare lots of vectors at the same time. Just use == and the recycling rule
    if(!all(c(
	dimnames(object@psi)[[1]],
	dimnames(object@intake_max)[[1]],
	dimnames(object@search_vol)[[1]],
	dimnames(object@std_metab)[[1]],
	dimnames(object@activity)[[1]],
	dimnames(object@pred_kernel)[[1]],
	dimnames(object@selectivity)[[2]],
	dimnames(object@catchability)[[2]],
	dimnames(object@interaction)[[1]],
	dimnames(object@interaction)[[2]]) ==
	    object@species_params$species)){
	    msg <- "The species names of species_params, psi, intake_max, search_vol, std_metab, activity, pred_kernel, selectivity, catchability and interaction must all be the same"
	    errors <- c(errors, msg)
    }
    # Check dimnames of w
    if(!all(c(
	dimnames(object@psi)[[2]],
	dimnames(object@intake_max)[[2]],
	dimnames(object@search_vol)[[2]],
	dimnames(object@std_metab)[[2]],
	dimnames(object@activity)[[2]],
	dimnames(object@pred_kernel)[[2]]) == 
	    dimnames(object@selectivity)[[3]])){
	    msg <- "The size names of psi, intake_max, search_vol, std_metab, activity, pred_kernel and selectivity must all be the same"
	    errors <- c(errors, msg)
    }
    # Check dimnames of gear
    if(!isTRUE(all.equal(
	dimnames(object@catchability)[[1]],
	dimnames(object@selectivity)[[1]]))){
	    msg <- "The gear names of selectivity and catchability must all be the same"
	    errors <- c(errors, msg)
    }
    # Check the vector slots
    if(length(object@rr_pp) != length(object@w_full)){
	msg <- "rr_pp must be the same length as w_full"
	errors <- c(errors, msg)
    }
    if(!isTRUE(all.equal(names(object@rr_pp),dimnames(object@pred_kernel)[[3]]))){
	msg <- "Names of rr_pp and third dimension of pred_kernel must be consistent"
	errors <- c(errors, msg)
    }
    if(length(object@cc_pp) != length(object@w_full)){
	msg <- "cc_pp must be the same length as w_full"
	errors <- c(errors, msg)
    }
    if(!isTRUE(all.equal(names(object@cc_pp),dimnames(object@pred_kernel)[[3]]))){
	msg <- "Names of cc_pp and third dimension of pred_kernel must be consistent"
	errors <- c(errors, msg)
    }

    # SRR
    # Must have two arguments: rdi amd species_params
    if(!isTRUE(all.equal(names(formals(object@srr)), c("rdi", "species_params")))){
	msg <- "Arguments of srr function must be 'rdi' and 'species_params'"
	errors <- c(errors, msg)
    }

    # species_params
    # Column check done in constructor
    # If everything is OK
    if (length(errors) == 0) TRUE else errors
}

# Try some documentation with roxygenize
# All lines start with single quote
# First line is title of the class
# Second line is details about the class (description)
# @name (override default topic name), @rdname (output filename), @export (for namespace)
# @slot doesn't work so need ugliness below
# @param w a numeric vector of size bins - can use this style but not really right
# Add Methods section?
# And See Also when there is something also to see

#' MizerParams
#'
#' A class to hold the parameters for the size based model. Dynamic simulations are performed using the project() method on this class. Objects of this class are created using the MizerParams() methods
#'
#' \section{Slots}{
#'     \describe{
#'         \item{w}{A numeric vector of size bins used for the community (i.e. fish) part of the model. These are usually spaced on a log10 scale}
#'         \item{dw}{The absolute difference between the size bins specified in the w slot. A vector the same length as the w slot. The final value is the same as the second to last value}
#'         \item{w_full}{A numeric vector of size bins used for the whole (i.e. community and background) model. These are usually spaced on a log10 scale}
#'         \item{dw_full}{The absolute difference between the size bins specified in the w_full slot. A vector the same length as the w_full slot. The final value is the same as the second to last value}
#'         \item{psi}{An array (species x size) that holds the allocation to reproduction for each species at size}
#'         \item{intake_max}{An array (species x size) that holds the maximum intake for each species at size}
#'         \item{search_vol}{An array (species x size) that holds the search volume for each species at size}
#'         \item{activity}{An array (species x size) that holds the activity for each species at size}
#'         \item{std_metab}{An array (species x size) that holds the standard metabolism for each species at size}
#'         \item{pred_kernel}{An array (species x predator size x prey size) that holds the predation coefficient of each predator at size on each prey size}
#'         \item{rr_pp}{A vector the same length as the w_full slot. The size specific growth rate of the background spectrum}
#'         \item{cc_pp}{A vector the same length as the w_full slot. The size specific carrying capacity of the background spectrum}
#'         \item{species_params}{A data.frame to hold the species specific parameters (see Details)}
#'         \item{interaction}{The species specific interation matrix PREDATOR OR PREY ON WHICH AXIS}
#'         \item{srr}{Function to calculate the realised (density dependent) recruitment. Has two arguments which are rdi and species_params}
#'         \item{selectivity}{An array (gear x species x w) that holds the selectivity of each species by gear and species size}
#'         \item{catchability}{An array (gear x species) that holds the catchability of each species by each gear}
#'     }
#' }
#' @name MizerParams-class
#' @rdname MizerParams-class
#' @docType class
#' @export
setClass("MizerParams",
    representation(
	w = "numeric",
	dw = "numeric",
	w_full = "numeric",
	dw_full = "numeric",
	psi = "array", 
	intake_max = "array",
	search_vol = "array",
	activity = "array",
	std_metab = "array",
	pred_kernel = "array",
	#z0 = "numeric",
	rr_pp = "numeric",
	cc_pp = "numeric", # was NinPP, carrying capacity of background
	species_params = "data.frame",
	interaction = "array",
	srr  = "function",
	selectivity = "array",
	catchability = "array"
    ),
    prototype = prototype(
	w = NA_real_,
	dw = NA_real_,
	w_full = NA_real_,
	dw_full = NA_real_,
	psi = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
	intake_max = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
	search_vol = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
	activity = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
	std_metab = array(NA,dim=c(1,1), dimnames = list(sp=NULL,w=NULL)),
	pred_kernel = array(NA,dim=c(1,1,1), dimnames = list(sp=NULL,w_pred=NULL,w_prey=NULL)),
	#z0 = NA_real_,
	rr_pp = NA_real_,
	cc_pp = NA_real_,
	#speciesParams = data.frame(),
	interaction = array(NA,dim=c(1,1), dimnames=list(predator=NULL, prey=NULL)), # which dimension is prey and which is prey?
	selectivity = array(NA, dim=c(1,1,1), dimnames = list(gear=NULL, sp=NULL, w=NULL)),
	catchability = array(NA, dim=c(1,1), dimnames = list(gear=NULL, sp=NULL))
    ),
    validity = valid_MizerParams
)


# Generic constructor

#' Constructors for objects of \code{MizerParams} class
#'
#' A range of constructors for the \code{MizerParams} class. Some are more user friendly than others and some are really only used 
#' internally to help construct \code{sbm_param} objects of the right dimension
#'
#' @param object Can either be an integer whose value determines the number of species in the community or a data.frame of species specific parameter values
#' @param interaction Optional argument to specify the interaction matrix of the species (predator by prey). If missing the a default interaction is used where all interactions are set to 1. 
#' @param ... Additional arguments used to specify the dimensions and sizes of the model. These include:
#'
#' \itemize{
#'     \item{\code{min_w} The smallest size of the community spectrum}
#'     \item{\code{max_w} The largest size of the community spectrum. Default value is the largest w_inf in the community x 1.1}
#'     \item{\code{no_w} The number of size bins in the community spectrum}
#'     \item{\code{min_w_pp} The smallest size of the background spectrum}
#'     \item{\code{no_w_pp} The number of the extra size bins in the background spectrum (i.e. the difference between the number of sizes bins in the community spectrum and the full spectrum)}
#'     \item{\code{species_names} A vector of the species names in the community}
#'     \item{\code{gear_names} A vector of the names of the fishing gears} 
#'     \item{\code{n} Scaling of the intake. Default value is 2/3} 
#'     \item{\code{p} Scaling of the standard metabolism. Default value is 0.7} 
#'     \item{\code{q} Exponent of the search volume. Default value is 0.8} 
#'     \item{\code{r_pp} Growth rate of the primary productivity. Default value is 10} 
#'     \item{\code{kappa} Carrying capacity of the resource spectrum. Default value is 1e11} 
#'     \item{\code{lambda} Exponent of the resource spectrum. Default value is (2+q-n)} 
#'     \item{\code{w_pp_cutoff} The cut off size of the background spectrum. Default value is 10} 
#' }
#'
#' @return An object of type \code{sbm_param}
#' @note or details Something on the default formulations of setting psi, intakeMax etc
#' @export
#' @docType methods
#' @rdname MizerParams-methods
#'
#' @examples
#' params <- MizerParams(object=3, species_names = c("cod", "haddock", "whiting"))
setGeneric('MizerParams', function(object, interaction, ...)
    standardGeneric('MizerParams'))

# Basic constructor with only the number of species as dispatching argument
# Only really used to make MizerParams of the right size

#' @rdname MizerParams-methods
#' @aliases MizerParams,numeric,missing-method
setMethod('MizerParams', signature(object='numeric', interaction='missing'),
    function(object, min_w = 0.001, max_w = 1000, no_w = 100,  min_w_pp = 1e-10, no_w_pp = round(no_w)*0.3, species_names=1:object, gear_names=species_names){
	#args <- list(...)

	# Some checks
	if (length(species_names) != object)
	    stop("species_names must be the same length as the value of object argument")

	# Set up grids:
	# Community grid
	w <- 10^(seq(from=log10(min_w),to=log10(max_w),length.out=no_w))
	dw <- diff(w)
	dw[no_w] <- dw[no_w-1] # Set final dw as same as one before

	# Set up full grid - background + community
	# ERROR if dw > w, nw must be at least... depends on minw, maxw and nw
	if(w[1] <= dw[1])
	    stop("Your size bins are too close together. You should consider increasing the number of bins, or changing the size range")
	w_full <- c(10^seq(from=log10(min_w_pp), to =  log10(w[1]-dw[1]),length.out=no_w_pp),w)
	no_w_full <- length(w_full)
	dw_full <- diff(w_full)
	dw_full[no_w_full] <- dw_full[no_w_full-1]

	# Basic arrays for templates
	mat1 <- array(NA, dim=c(object,no_w), dimnames = list(sp=species_names,w=signif(w,3)))
	mat2 <- array(NA, dim=c(object,no_w,no_w_full), dimnames = list(sp=species_names,w_pred=signif(w,3), w_prey=signif(w_full,3)))
	selectivity <- array(0, dim=c(length(gear_names), object, no_w), dimnames=list(gear=gear_names, sp=species_names, w=signif(w,3)))
	catchability <- array(0, dim=c(length(gear_names), object), dimnames = list(gear=gear_names, sp=species_names))
	interaction <- array(1, dim=c(object,object), dimnames = list(predator = species_names, prey = species_names))
	vec1 <- as.numeric(rep(NA, no_w_full))
	names(vec1) <- signif(w_full,3)
	
	# Make an empty data.frame for species_params
	# This is just to pass validity check. There should be a seperate function to check if the species_params data.frame has all the right columns and dims
	species_params <- data.frame(species = species_names)

	# Make an empty srr function, just to pass validity check
	srr <- function(rdi, species_params) return(0)

	# Make the new object
	# Should Z0, rrPP and ccPP have names (species names etc)?
	res <- new("MizerParams",
	    w = w, dw = dw, w_full = w_full, dw_full = dw_full,
	    psi = mat1, intake_max = mat1, search_vol = mat1, activity = mat1, std_metab = mat1, pred_kernel = mat2,
	    selectivity=selectivity, catchability=catchability,
	    rr_pp = vec1, cc_pp = vec1, species_params = species_params,
	    interaction = interaction, srr = srr) 
	return(res)
    }
)

# Constructor that takes the species_params data.frame and the interaction matrix

#' @rdname MizerParams-methods
#' @aliases MizerParams,data.frame,matrix-method
setMethod('MizerParams', signature(object='data.frame', interaction='matrix'),
    function(object, interaction,  n = 2/3, p = 0.7, q = 0.8, r_pp = 10, kappa = 1e11, lambda = (2+q-n), w_pp_cutoff = 10, ...){
	args <- list(...)

	# Check species_params dataframe (with a function) for essential cols
	check_species_params_dataframe(object)
	# Essential columns: species (name) # wInf # wMat # h # gamma - search Volume # k # ks # beta # sigma # alpha
	# If no gear_name column in object, then named after species
	if(!("gear" %in% colnames(object)))
	    object$gear <- object$species
	no_gear <- length(unique(object$gear))

	no_sp <- nrow(object)

	# if max_w not passed in by user , set from w_inf in data_fram
	if ("max_w" %in% names(args))
	    max_w <- args[["max_w"]]
	else max_w <- max(object$w_inf)*1.1

	# Make an empty object of the right dimensions
	res <- MizerParams(no_sp, species_names=object$species, gear_names=unique(object$gear), max_w=max_w,...)

	# If not w_min column in species_params, set to w_min of community
	# Check min_w argument is not > w_min in species_params
	if (!("w_min" %in% colnames(object)))
	    object$w_min <- min(res@w)
	if(any(object$w_min < min(res@w)))
	    stop("One or more of your w_min values is less than the smallest size of the community spectrum")

	# Add w_min_idx column which has the reference index of the size class closest to w_min - this is a short cut for later on and prevents repetition
	object$w_min_idx <- as.vector(tapply(object$w_min,1:length(object$w_min),function(w_min,wx) max(which(wx<=w_min)),wx=res@w))

	# If no sel_func column in species_params, set to 'sigmoid_length'
	if(!("sel_func" %in% colnames(object)))
	    object$sel_func <- 'sigmoid_length'
	# Start filling the slots
	res@species_params <- object
	# Check dims of interaction argument - make sure it's right
	if (!isTRUE(all.equal(dim(res@interaction), dim(interaction))))
	    stop("interaction matrix is not of the right dimensions. Must be number of species x number of species")
	res@interaction[] <- interaction

	# Now fill up the slots using default formulations:
	# psi - allocation to reproduction - from original Setup() function
	res@psi[] <- unlist(tapply(res@w,1:length(res@w),function(wx,w_inf,w_mat,n){
	    ((1 + (wx/(w_mat))^-10)^-1) * (wx/w_inf)^(1-n)},w_inf=object$w_inf,w_mat=object$w_mat,n=n))
	# Set w < 10% of w_mat to 0
	res@psi[unlist(tapply(res@w,1:length(res@w),function(wx,w_mat)wx<(w_mat*0.1)  ,w_mat=object$w_mat))] <- 0
	# Set all w > w_inf to 1 # Check this is right...
	res@psi[unlist(tapply(res@w,1:length(res@w),function(wx,w_inf)(wx/w_inf)>1,w_inf=object$w_inf))] <- 1

	res@intake_max[] <- unlist(tapply(res@w,1:length(res@w),function(wx,h,n)h * wx^n,h=object$h,n=n))
	res@search_vol[] <- unlist(tapply(res@w,1:length(res@w),function(wx,gamma,q)gamma * wx^q, gamma=object$gamma, q=q))
	res@activity[] <-  unlist(tapply(res@w,1:length(res@w),function(wx,k)k * wx,k=object$k))
	res@std_metab[] <-  unlist(tapply(res@w,1:length(res@w),function(wx,ks,p)ks * wx^p, ks=object$ks,p=p))
	# Could maybe improve this. Pretty ugly at the moment
	res@pred_kernel[] <- object$beta
	res@pred_kernel <- exp(-0.5*sweep(log(sweep(sweep(res@pred_kernel,3,res@w_full,"*")^-1,2,res@w,"*")),1,object$sigma,"/")^2)
	res@pred_kernel <- sweep(res@pred_kernel,c(2,3),combn(res@w_full,1,function(x,w)x<w,w=res@w),"*") # find out the untrues and then multiply

	# How to choose which z0 function? z0 set before
	#res@z0 = z0pre*object$wInf^z0exp    # background natural mortality
	#Z0= (param$Z0pre*sp$k_vb)   # background mortality of Pope et al. 2006

	# Background spectrum
	res@rr_pp[] <- r_pp * res@w_full^(n-1) #weight specific plankton growth rate ##
	res@cc_pp[] <- kappa*res@w_full^(-lambda) # the resource carrying capacity - one for each mp and m (130 of them)
	res@cc_pp[res@w_full>w_pp_cutoff] <- 0      #set density of sizes < plankton cutoff size
	# Set the SRR to be a Beverton Holt esque relationship
	# Can add more functional forms or user specifies own
	res@srr <- function(rdi, species_params){
	    return(species_params$r_max * rdi / (species_params$r_max+rdi))
	}

	# Set fishing parameters: selectivity and catchability
	# At the moment, each species is only caught by 1 gear so in species_params there are the columns: gear_name and sel_func
	# BEWARE! This routine assumes that each species has only one gear operating on it
	# So we can just go row by row through the species parameters
	# However, I really hope we can do something better soon
	for (g in 1:nrow(object)){
	    # Do selectivity first
	    # get args
	    # These as.characters are annoying - but factors everywhere
	    arg <- names(formals(as.character(object[g,'sel_func'])))
	    # lop off w as that is always the first argument of the selectivity functions
	    arg <- arg[!(arg %in% "w")]
	    if(!all(arg %in% colnames(object)))
		stop("The arguments needed for the selectivity function are notin the object argement")
	    # Check that there is only one column in object with the same name
	    # Check that column of arguments exists
	    par <- c(w=list(res@w),as.list(object[g,arg]))
	    sel <- do.call(as.character(object[g,'sel_func']), args=par)
	    # Dump Sel in the right place
	    res@selectivity[as.character(object[g,'gear']), g, ] <- sel
	    # Now do catchability
	    res@catchability[as.character(object[g,'gear']), g] <- object[g,"catchability"]
	}

	return(res)
    }
)

# If interaction is missing, make one of the right size and fill with 1s
#' @rdname MizerParams-methods
#' @aliases MizerParams,data.frame,missing-method
setMethod('MizerParams', signature(object='data.frame', interaction='missing'),
    function(object, ...){
	interaction <- matrix(1,nrow=nrow(object), ncol=nrow(object))
	res <- MizerParams(object,interaction, ...)
	return(res)
})

# Check that the species_params dataset is OK
# internal only
check_species_params_dataframe <- function(species_params){
    # Check species_params dataframe (with a function) for essential cols
    # Essential columns: species (name) # wInf # wMat # h # gamma - search Volume # k # ks # beta # sigma # alpha
    essential_cols <- c("species","w_inf","w_mat","h","gamma","k","ks","beta","sigma","alpha", "erepro")
    missing_cols <- !(essential_cols %in% colnames(species_params))
    if(any(missing_cols))
    {
	errors <- character()
	for (i in essential_cols[missing_cols])
	    errors <- paste(errors, i, sep=" ")
	stop("You are missing these columns from the input dataframe: ", errors)
    }
    return(TRUE)
}
