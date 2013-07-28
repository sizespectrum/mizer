\name{set_trait_model}
\alias{set_trait_model}
\title{Sets up parameters for a trait-based model}
\usage{
  set_trait_model(no_sp = 10, min_w_inf = 10,
    max_w_inf = 1e+06, no_w = 100, min_w = 0.001,
    max_w = max_w_inf * 1.1, min_w_pp = 1e-10,
    no_w_pp = round(no_w) * 0.3, w_pp_cutoff = 1, k0 = 100,
    n = 2/3, p = 0.75, q = 0.9, eta = 0.25, r_pp = 4,
    kappa = 0.005, lambda = 2 + q - n, alpha = 0.6, ks = 4,
    z0 = 0.6, h = 30, beta = 100, sigma = 1.3, f0 = 0.5,
    knife_edge_size = 1000, knife_is_min = TRUE, ...)
}
\arguments{
  \item{no_sp}{The number of species in the model. The
  default value is 10. The more species, the longer takes
  to run.}

  \item{min_w_inf}{The asymptotic size of the smallest
  species in the community.}

  \item{max_w_inf}{The asymptotic size of the largest
  species in the community.}

  \item{no_w}{The number of size bins in the community
  spectrum.}

  \item{min_w}{The smallest size of the community
  spectrum.}

  \item{max_w}{The largest size of the community spectrum.
  Default value is the largest w_inf in the community x
  1.1.}

  \item{min_w_pp}{The smallest size of the background
  spectrum.}

  \item{no_w_pp}{The number of the extra size bins in the
  background spectrum (i.e. the difference between the
  number of sizes bins in the community spectrum and the
  full spectrum).}

  \item{w_pp_cutoff}{The cut off size of the background
  spectrum. Default value is 1.}

  \item{k0}{Scaling parameter for the maximum recruitment.}

  \item{n}{Scaling of the intake. Default value is 2/3.}

  \item{p}{Scaling of the standard metabolism. Default
  value is 0.7.}

  \item{q}{Exponent of the search volume. Default value is
  0.8.}

  \item{eta}{Factor to calculate \code{w_mat} from
  asymptotic size.}

  \item{r_pp}{Growth rate of the primary productivity.
  Default value is 10.}

  \item{kappa}{Carrying capacity of the resource spectrum.
  Default value is 1e11.}

  \item{lambda}{Exponent of the resource spectrum. Default
  value is (2+q-n).}

  \item{alpha}{The assimilation efficiency of the
  community. The default value is 0.2 (from Andersen et.
  al., 2010).}

  \item{ks}{Standard metabolism coefficient. Default value
  is 4.}

  \item{z0}{The background mortality of the community. The
  default value is 0.6.}

  \item{h}{Maximum food intake rate. Default value is 30.}

  \item{beta}{Preferred predator prey mass ratio. Default
  value is 100.}

  \item{sigma}{Width of prey size preference. Default value
  is 1.3.}

  \item{f0}{Expected average feeding level. Used to set
  \code{gamma}, the factor for the search volume.}

  \item{knife_edge_size}{The size at the edge of the
  knife-selectivity function.}

  \item{knife_is_min}{Is the knife-edge selectivity
  function selecting above (TRUE) or below (FALSE) the
  edge.}
}
\value{
  An object of type \code{MizerParams}
}
\description{
  This functions creates a \code{MizerParams} object so
  that trait-based-type models can be easily set up and
  run. The trait-based size spectrum model can be derived
  as a simplification of the general size-based model used
  in \code{mizer}. All the species-specific parameters are
  the same, except for the asymptotic size, which is
  considered the most important trait characterizing a
  species. Other parameters are related to the asymptotic
  size. For example, the size at maturity is given by w_inf
  * eta, where eta is the same for all species. For the
  trait-based model the number of species is not important.
  For applications of the trait-based model see Andersen &
  Pedersen (2010). See the \code{mizer} vignette for more
  details and examples of the trait-based model.
}
\details{
  The function has many arguments, all of which have
  default values. Of particular interest to the user are
  the number of species in the model and the minimum and
  maximum asymptotic sizes. The asymptotic sizes of the
  species are spread evenly on a logarithmic scale within
  this range.

  The stock recruitment relationship in the trait-based
  model is a 'hockey-stick' (following Andersen and
  Pedersen, 2009). The maximum value of the recruitment can
  be controlled using the \code{k0} parameter. Users should
  adjust this value to get the spectra they want.

  The factor for the search volume, \code{gamma}, is
  calculated using the expected feeding level, \code{f0}.

  Fishing selectivity is modelled as a sigmoid function
  with two parameters: \code{l25} and \code{l50} which
  determine the lengths at which selectivity is 0.25 and
  0.5 respectively. Lengths are converted to weights using
  the default parameters a = 0.001 and b = 3.0.

  The resulting \code{MizerParams} object can be projected
  forward using \code{project()} like any other
  \code{MizerParams} object. When projecting the community
  model it may be necessary to reduce \code{dt} to 0.1 to
  avoid any instabilities with the solver. You can check
  this by plotting the biomass or abundance through time
  after the projection.
}
\examples{
trait_params <- set_trait_model(no_sp = 15)
init_pop <- get_initial_n(trait_params, n0_mult = 0.001)
sim <- project(trait_params, effort = 0, t_max = 50, dt=0.2, initial_n = init_pop, t_save = 1)
plot(sim)
}
\references{
  K. H. Andersen and M. Pedersen, 2010, Damped trophic
  cascades driven by fishing in model marine ecosystems.
  Proceedings of the Royal Society V, Biological Sciences,
  1682, 795-802.
}
\seealso{
  \code{\link{MizerParams}}
}

