\docType{methods}
\name{MizerParams}
\alias{MizerParams}
\alias{MizerParams,data.frame,matrix-method}
\alias{MizerParams,data.frame,missing-method}
\alias{MizerParams,numeric,missing-method}
\title{Constructors for objects of \code{MizerParams} class}
\arguments{
  \item{object}{Can either be an integer whose value
  determines the number of species in the community or a
  data.frame of species specific parameter values}

  \item{interaction}{Optional argument to specify the
  interaction matrix of the species (predator by prey). If
  missing the a default interaction is used where all
  interactions are set to 1.}

  \item{...}{Additional arguments used to specify the
  dimensions and sizes of the model. These include:

  \itemize{ \item{\code{min_w} The smallest size of the
  community spectrum} \item{\code{max_w} The largest size
  of the community spectrum. Default value is the largest
  w_inf in the community x 1.1} \item{\code{no_w} The
  number of size bins in the community spectrum}
  \item{\code{min_w_pp} The smallest size of the background
  spectrum} \item{\code{no_w_pp} The number of the extra
  size bins in the background spectrum (i.e. the difference
  between the number of sizes bins in the community
  spectrum and the full spectrum)}
  \item{\code{species_names} A vector of the species names
  in the community} \item{\code{gear_names} A vector of the
  names of the fishing gears} \item{\code{n} Scaling of the
  intake. Default value is 2/3} \item{\code{p} Scaling of
  the standard metabolism. Default value is 0.7}
  \item{\code{q} Exponent of the search volume. Default
  value is 0.8} \item{\code{r_pp} Growth rate of the
  primary productivity. Default value is 10}
  \item{\code{kappa} Carrying capacity of the resource
  spectrum. Default value is 1e11} \item{\code{lambda}
  Exponent of the resource spectrum. Default value is
  (2+q-n)} \item{\code{w_pp_cutoff} The cut off size of the
  background spectrum. Default value is 10} }}
}
\value{
  An object of type \code{MizerParams}
}
\description{
  A range of constructors for the \code{MizerParams} class.
  Some are more user friendly than others and some are
  really only used internally to help construct
  \code{sbm_param} objects of the right dimension
}
\note{
  Something on the default formulations of setting psi,
  intakeMax etc
}

