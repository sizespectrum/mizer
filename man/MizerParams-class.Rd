\docType{class}
\name{MizerParams-class}
\alias{MizerParams-class}
\title{MizerParams}
\description{
  A class to hold the parameters for the size based model.
  Dynamic simulations are performed using the project()
  method on this class. Objects of this class are created
  using the MizerParams() methods
}
\details{
  \section{Slots}{ \describe{ \item{w}{A numeric vector of
  size bins used for the community (i.e. fish) part of the
  model. These are usually spaced on a log10 scale}
  \item{dw}{The absolute difference between the size bins
  specified in the w slot. A vector the same length as the
  w slot. The final value is the same as the second to last
  value} \item{w_full}{A numeric vector of size bins used
  for the whole (i.e. community and background) model.
  These are usually spaced on a log10 scale}
  \item{dw_full}{The absolute difference between the size
  bins specified in the w_full slot. A vector the same
  length as the w_full slot. The final value is the same as
  the second to last value} \item{psi}{An array (species x
  size) that holds the allocation to reproduction for each
  species at size} \item{intake_max}{An array (species x
  size) that holds the maximum intake for each species at
  size} \item{search_vol}{An array (species x size) that
  holds the search volume for each species at size}
  \item{activity}{An array (species x size) that holds the
  activity for each species at size} \item{std_metab}{An
  array (species x size) that holds the standard metabolism
  for each species at size} \item{pred_kernel}{An array
  (species x predator size x prey size) that holds the
  predation coefficient of each predator at size on each
  prey size} \item{rr_pp}{A vector the same length as the
  w_full slot. The size specific growth rate of the
  background spectrum} \item{cc_pp}{A vector the same
  length as the w_full slot. The size specific carrying
  capacity of the background spectrum}
  \item{species_params}{A data.frame to hold the species
  specific parameters (see Details)} \item{interaction}{The
  species specific interation matrix PREDATOR OR PREY ON
  WHICH AXIS} \item{srr}{Function to calculate the realised
  (density dependent) recruitment. Has two arguments which
  are rdi and species_params} \item{selectivity}{An array
  (gear x species x w) that holds the selectivity of each
  species by gear and species size} \item{catchability}{An
  array (gear x species) that holds the catchability of
  each species by each gear} } }
}

