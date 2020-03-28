#' mizer: Multi-species size-based modelling in R
#'
#' The mizer package implements multi-species size-based modelling in R. It has 
#' been designed for modelling marine ecosystems.
#'
#' Using \pkg{mizer} is relatively simple.  There are three main stages: 
#' \enumerate{
#'
#' \item Setting the model parameters. This is done by creating an object of
#' class \linkS4class{MizerParams}. This includes model parameters such as the
#' life history parameters of each species, and the range of the size spectrum.
#' There are several setup functions that help to create a MizerParams objects
#' for particular types of models:
#' \itemize{
#'   \item [newCommunityParams()]
#'   \item [newTraitParams()]
#'   \item [newMultispeciesParams()]
#' }
#' \item Running a simulation. This is done by calling the
#' [project()] function with the model parameters. This produces an
#' object of \linkS4class{MizerSim} that contains the results of the simulation.
#'
#' \item Exploring results. After a simulation has been run, the results can be
#' explored using a range of [plotting_functions()] and
#' [summary_functions()].
#' }
#'
#' See the mizer website and vignettes for full details of the principles behind
#' mizer and how the package can be used to perform size-based modelling.
#'
#' @import ggplot2 methods assertthat shiny dplyr
#' @importFrom plotly ggplotly plotlyOutput renderPlotly
#' @importFrom stats fft mvfft lm pnorm runif complete.cases
#' @importFrom grDevices col2rgb
#' @importFrom utils modifyList
#' @docType package
#' @name mizer
#' @aliases mizer-package
NULL

#' @importFrom reshape2 melt
#' @export
reshape2::melt

#' Description of the general mizer size-spectrum model
#' 
#' On this page we will present the general mizer model
#' We will provide links to the functions that can be used to set
#' or change the various model parameters (all collected together in `setParams()`)
#' as well as to the functions that calculate the various ecological rates in the
#' model (all collected together in `getRates()`) because the help pages of these
#' functions will provide useful additional details.
#' 
#' @section Consumer size spectra:
#' The model assumes that to a first approximation an individual can be characterized by its 
#' weight \eqn{w} and its species \eqn{i} only. The aim of the model is to calculate 
#' the size spectrum \eqn{N_i(w)}, which is the *density* of individuals of
#' species \eqn{i} such that \eqn{\int_w^{w+dw}N_i(w)dw} is the *number* of individuals of
#' species \eqn{i} in the size interval \eqn{[w,w+dw]}. In other word: the number of 
#' individuals in a size range is the area under the number density \eqn{N_i(w)}.
# 
# Here is a plot of an example size spectrum for two species with
# \eqn{N_i(w)} on the vertical axis for \eqn{i=1,2} and \eqn{w} on the horizontal axis.
# ```{r, message=FALSE}
# library(mizer)
# params <- newTraitParams(no_sp = 2, min_w = 1e-3)
# plotSpectra(params, plankton = FALSE, power = 0)
# ```
#' 
#' To represent this continuous size spectrum in the computer, the size
#' variable \eqn{w} is discretized into a vector `w` of discrete weights,
#' providing a grid of sizes spanning the range from the smallest egg size
#' to the largest asymptotic size. These grid values divide the full size
#' range into a finite number of size bins. The size bins should be chosen
#' small enough to avoid the discretisation errors from becoming too big.
#' 
#' The weight grid is set up to be logarithmically spaced, so that
#' `w[j]=w[1]*exp(j*dx)` for some fixed `dx`. 
#' This grid is set up automatically when creating a MizerParams object.
#' 
#' In the code the size spectrum is stored as an array such that `n[i, a]`
#' holds the density \eqn{N_i(w_a)} at weights \eqn{w_a=}`w[a]`, or, if time
#' dependence is included, an array such that `n[i, a, u]`
#' holds \eqn{N_i(w_a,t_u)}.
#' 
#' Note that, contrary to what one might have expected,
#' `n[i, a]` is not the *number* of individuals in a size bin
#' but the *density* at a grid point.
#' The number of individuals in the size bin between `w[a]` and
#' `w[a+1]=w[a]+dw[a]` is only approximately given as `n[i, a]*dw[a]`,
#' where `dw[a]= w[a+1]-w[a]`.
#'
#'     
#' The time evolution of the size spectrum is described by the 
#' McKendrik-von Foerster equation, which is a continuity equation with loss:
#' \deqn{
#' \frac{\partial N_i(w)}{\partial t} + \frac{\partial g_i(w) N_i(w)}{\partial w} 
#' = -\mu_i(w) N_i(w),
#' }
#' where individual growth \eqn{g_i(w)} and mortality \eqn{\mu_i(w)} will be described
# below.
# 
# <details>
#     This McKendrik-von Foerster equation is approximated in mizer by a
# finite-difference method (to be described in section ...). This allows the
# `project()` function in mizer to project the size
# spectrum forwards in time: Given the spectrum at one time the `project()`
# function calculates it at a set of later times.
# 
# </details>
#     
#' The McKendrik-von Foerster equation is supplemented by a boundary condition at
#' the egg weight \eqn{w_0} where the flux of individuals (numbers per time) 
#' \eqn{g_i(w_0)N_i(w_0)} is determined by the rate \eqn{R_i} of production of offspring by
#' mature individuals in the population:
#' \deqn{
#' g_i(w_0)N_i(w_0) = R_i.
#' }
# 
# In the code this boundary condition is actually implemented as an equation for
# the rate of change of the number of individuals in the smallest size bin:
#     \begin{equation}
# \frac{d N_i(w_0)}{dt}=...
# \end{equation}
#     
#' @section Plankton size spectrum:
#' Besides the fish spectrum there is also a plankton spectrum \eqn{N_R(w)},
#' representing for example the phytoplankton. This spectrum starts at a smaller
#' size than the fish spectrum, in order to provide food also for the smallest
#' individuals (larvae) of the fish spectrum.
#' 
#' By default the time evolution of the plankton spectrum is described by a 
#' semi-chemostat equation:
#' \deqn{
#' \frac{\partial N_R(w,t)}{\partial t} 
#' = r_p(w) \Big[ c_p (w) - N_R(w,t) \Big] - \mu_p(w) N_R(w,t).
#' }
#' Here \eqn{r_p(w)} is the plankton regeneration rate and \eqn{c_p(w)} is the carrying
#' capacity in the absence of predation. These parameters are changed with
#' `setPlankton()`. By default mizer assumes allometric forms
#' \deqn{r_p(w)= r_p\, w^{n-1}.}
#' \deqn{c_p(w)=\kappa\, w^{-\lambda}.}
#' It is also possible to implement other plankton dynamics, as
#' described in the help page for `setPlankton()`. The death \eqn{\mu_p(w)} is
#' described in the subsection [Plankton mortality].
#' 
#' Because the plankton spectrum spans a different range of sizes these sizes
#' are discretized into a different vector of weights
#' `w_full`. The last entries of `w_full` have to coincide with the entries
#' of `w`. The plankton spectrum is then stored in a vector `n_pp` such that
#' `n_pp[c]` =\eqn{N_R(}`w_full[c]`\eqn{)}.
# 
# The plankton regeneration rate is stored as a vector
# `params@rr_pp[c]`\eqn{=r_p(}`w_full[c]`\eqn{)}, and the carrying capacity
# is stored as a vector
# `params@cc_pp[c]`\eqn{=c_p(}`w_full[c]`\eqn{)}.
#     
#     
#     
#' @section Predator-prey encounter rate:
#' The rate at which a predator of species \eqn{i} and weight \eqn{w} encounters food 
#' (mass per time) is
#' determined by summing over all prey species and the plankton spectrum and
#' integrating over all prey sizes \eqn{w_p}, weighted by the selectivity factors:
#     \begin{equation}
# \label{eq:1}
# E_{i}(w) = \gamma_i(w) \int \left(\sum_{j} \theta_{ij} N_j(w_p) +
#                                        \theta_{ip} N_R(w_p) + \right) 
# \phi_i(w,w_p) w_p \, dw_p.
# \end{equation}
# This is calculated by `getEncounter()`.
# The overall prefactor \eqn{\gamma_i(w)} sets the predation power of the predator. It
# could be interpreted as a search volume. It is set by `setSearchVolume()`. By 
# default it is assumed to scale allometrically as
# \eqn{\gamma_i(w) = \gamma_i\, w^q.}
#     
#     The \eqn{\theta} matrix sets the interaction strength between predators and the
# various prey species and plankton. It is changed with `setInteraction()`.
# 
# The size selectivity is encoded in the predation kernel \eqn{\phi_i(w,w_p)}. This is
# changed with `setPredKernel()`.
# 
# <details>
#     An important simplification occurs when the predation kernel \eqn{\phi_i(w,w_p)}
#     depends on the size of the prey **only** through the predator/prey size ratio
# \eqn{w_p/w},
# \[\phi_i(w, w_p)=\tilde{\phi}_i(w/w_p).\]
# This is assumed by default but can be overruled.
# The default for the predation kernel is the truncated log-normal function
# \[
#     \label{eq:4}
#     \tilde{\phi}_i(x) = \begin{cases}
#     \exp \left[ \dfrac{-(\ln(x / \beta_i))^2}{2\sigma_i^2} \right]
#     &\text{ if }x\in\left[0,\beta_i\exp(3\sigma_i)\right]\\
#     0&\text{ otherwise,}
#     \end{cases}
#     \]
# where \eqn{\beta_i} is the preferred predator-prey mass ratio and \eqn{\sigma_i} sets
# the width of the predation kernel.
# 
# The integral in the expression for the encounter rate is approximated by a 
# Riemann sum over all weight brackets:
#     \[
#         {\tt encounter}[i,a] = {\tt search\_vol}[i,a]\sum_{k}
#         \left( n_{pp}[k] + \sum_{j} \theta[i,j] n[j,k] \right) 
#         \phi_i\left(w[a],w[k]\right) w[k]\, dw[k].
#         \]
# In the case of a predation kernel that depends on \eqn{w/w_p} only, this becomes
# a convolution sum and can be evaluated efficiently via fast Fourier transform.
# 
# </details>
#     
#     
#     ## Consumption
#     
#     The encountered food is consumed subjected to a standard Holling functional 
# response type II to represent satiation. This determines the 
# *feeding level* \eqn{f_i(w)}, which is a dimensionless number between 0 
# (no food) and 1 (fully satiated) so that \eqn{1-f_i(w)} is the proportion of the
# encountered food that is consumed. The feeding level is given by
# 
# \begin{equation}
# \label{eq:f}
# f_i(w) = \frac{E_{i}(w)}{E_{e.i}(w) + h_i(w)},
# \end{equation}
# 
# where \eqn{h_i(w)} is the maximum consumption rate. This is changed with
# `setMaxIntakeRate()`. By default mizer assumes an allometric form
# \eqn{h_i(w) = h_i\, w^n.}
#     The feeding level is calculated with the function `getFeedingLevel()`.
# 
# The rate at which food is consumed is then 
# \begin{equation}
# (1-f_i(w))E_{i}(w)=f_i(w)\, h_i(w).
# \end{equation}
# 
# ## Metabolic losses
# 
# Some of the consumed food is used to fuel the needs for metabolism and 
# activity and movement, at a rate \eqn{{\tt metab}_i(w)}. By default 
# this is made up out of standard metabolism, scaling with exponent \eqn{p}, and
# loss due to activity and movement, scaling with exponent \eqn{1}:
#     \[{\tt metab}_i(w) = k_{s.i}\,w^p + k_i\,w.\]
# See the help page for `setMetabolicRate()`. 
# 
# The remaining rate, if any, is assimilated with an efficiency \eqn{\alpha_i} and 
# is then available for growth and reproduction. So the rate at which energy
# becomes available for growth and reproduction is
# \begin{equation}
# \label{eq:Er}
# E_{r.i}(w) = \max(0, \alpha_i f_i(w)\, h_i(w) - {\tt metab}_i(w))
# \end{equation}
# This is calculated with the `getEReproAndGrowth()` function.
# 
# 
# ## Investment into reproduction  {#sec:repro}
# A proportion \eqn{\psi_i(w)} of the energy available for growth and reproduction is
# used for reproduction. This proportion should change from zero below the weight
# \eqn{w_{m.i}} of maturation to one at the asymptotic weight \eqn{w_{\infty.i}}, where
# all available energy is used for reproduction. This function is changed with
# `setReproduction()`. Mizer provides a default form for the function which you
# can however overrule.
# 
# The total production rate of egg production \eqn{R_{p.i}} (numbers per year) is
# found by integrating the contribution from all individuals of species \eqn{i}:
#     \begin{equation}
# \label{eq:Rp}
# R_{p.i} = \frac{\epsilon}{2 w_0} \int N_i(w)  E_{r.i}(w) \psi_i(w) \, dw,
# \end{equation}
# where the individual contribution is obtained by multiplying the rate at which the
# individual allocates energy to reproduction by an efficiency factor \eqn{\epsilon} 
#     and then dividing by the egg weight \eqn{w_0} to convert the energy into number of eggs.
# The result is multiplied by a factor \eqn{1/2} to take into account that only 
# females reproduce.
# 
# 
# ## Density-dependence in reproduction
# In mizer, density dependence is modelled as a compensation on the egg 
# production. This gives an additional contribution to the stock-recruitment 
# relationship. The default functional form of this density dependence is such 
# that the reproduction rate \eqn{R_i} (numbers per time) approaches a maximum as the 
# energy invested in reproduction increases, modelled mathematically analogous to 
# the Holling type II function response as a *Beverton-Holt* type function:
#     
#     \begin{equation}
# \label{eq:R}
# R_i = R_{\max.i} \frac{R_{p.i}}{R_{p.i} + R_{\max.i}},
# \end{equation}
# where \eqn{R_{\max.i}} is the maximum reproduction rate of each trait class (Figure~\ref{fig:recruitment}). 
# 
# The *Beverton-Holt* type is not the only density dependence model that mizer can
# use. Users are able to write their own functions, e.g. fixed reproduction (as
#                                                                            used in the community-type model) or *hockey-stick*.
# 
# 
# ## Growth
# What is left over after metabolism and reproduction is taken into account
# is invested in somatic growth. Thus the growth rate is
# \begin{equation}
# \label{eq:growth}
# g_i(w) = E_{r.i}(w)\left(1-\psi_i(w)\right).
# \end{equation}
# It is calculated by the `getEGrowth()` function.
# 
# When food supply does not cover the requirements of metabolism and activity, 
# growth and reproduction stops, i.e. there is no negative growth.
# The individual should then be subjected to a starvation mortality, but starvation
# mortality is not implemented in mizer at the moment.
# 
# 
# ## Mortality
# The mortality rate of an individual \eqn{\mu_i(w)} has three sources: 
#     predation mortality \eqn{\mu_{p.i}(w)}, background mortality \eqn{\mu_{b.i}(w)} and 
# fishing mortality \eqn{\mu_{f.i}(w)}. 
# 
# Predation mortality is calculated such that all that is eaten translates into 
# corresponding predation mortalities on the ingested prey individuals. 
# Recalling that \eqn{1-f_j(w)} is the proportion of the food encountered by a 
# predator of species \eqn{j} and weight \eqn{w} that is actually consumed, the
# rate at which all predators of species \eqn{j} consume prey of size \eqn{w_p} is
# \begin{equation}
# \label{eq:pred_rated}
# {\tt pred\_rate}_j(w_p) = \int \phi_j(w,w_p) (1-f_j(w))
# \gamma_j(w) N_j(w) \, dw.
# \end{equation}
# This predation rate is calculated by the function `getPredRate()`.
# 
# <details>
#     The integral is approximated by a Riemann sum over all fish weight brackets.
# \[
#     {\tt pred\_rate}[j,c] = \sum_{a}
#     {\tt pred_kernel}[j,a,c]\,(1-{\tt feeding_level}[j,a])\,
#     \gamma[j,a]\,n[j,a]\,dw[a].
#     \]
# 
# </details>
#     
#     The mortality rate due to predation is then obtained as
# \begin{equation}
# \label{eq:mup}
# \mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, \theta_{ji}.
# \end{equation}
# This predation mortality rate is calculated by the function `getPredMort()`.
# 
# External mortality \eqn{z0_i(w)} is independent of the abundances and is
# changed with `setExtMort()`. By default mizer assumes an allometric form
# \[z0_i(w) = z0_{pre} w_{\infty.i}^{1-n},\]
# where \eqn{w_{\infty.i}} is the asymptotic size of species \eqn{i}.
# 
# Fishing mortality \eqn{\mu_{f.i}(w)} will be discussed in section [Fishing].
# 
# The total mortality rate
# \[\mu_i(w)=\mu_{p.i}(w)+z0_i(w)+\mu_{f.i}(w)\]
# is calculated with the function `getMort()`.
# 
# 
# ## Plankton Mortality
# 
# The predation mortality rate on plankton is given by a similar expression
# as the predation mortality on fish:
#     \begin{equation}
# \label{eq:mupp}
# \mu_{p}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, \theta_{jp}.
# \end{equation}
# This is the only mortality on plankton currently implemented in mizer.
# It is calculated with the function `getPlanktonMort()`.
# 
# 
# ## Fishing
# This section still to be written.
# 

#' @name mizer_model
NULL