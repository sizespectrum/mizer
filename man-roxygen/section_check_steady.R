#' @section What you get back may not be a steady state:
#'
#'   The stopping criterion is a proxy. It says that two states `t_per` years
#'   apart differ by less than `tol` on whatever scale the criterion is measured
#'   on; it does not say that the state reached is a fixed point. There are four
#'   ways the returned object can fail to be one:
#'
#'   - the run converged on its own scale while the biomasses are still
#'     visibly drifting, because `tol` was too loose (`type =
#'     "below_tolerance"`, which is why that type is not called `"steady"`);
#'   - the run reached `t_max` without converging (`type = "not_converged"`);
#'   - the run settled on a limit cycle (`type = "cycle"`), in which case the
#'     state stored is one point on that cycle;
#'   - the run stopped because a species was going extinct
#'     (`type = "extinction"`).
#'
#'   So treat the result as a claim to be checked rather than as a guarantee:
#'
#'   ```r
#'   attr(params, "convergence")$type      # "below_tolerance", "cycle", "not_converged", "extinction"
#'   attr(params, "convergence")$residual  # largest biomass drift, in 1/year
#'   isSteady(params)                      # TRUE if within tolerance
#'   summary(params)                       # includes the biomass-drift verdict
#'   plot(getSteadyResidual(params))       # which species, and at which sizes
#'   ```
#'
#'   The first four say *whether* the model is steady, the last says *where* it
#'   is not, which is the one to reach for when it is not: a model that is off
#'   steady state is usually off in one species or one part of the size range,
#'   and the plot names it. See [getSteadyResidual()] for why the verdict is
#'   phrased in terms of biomass drift rather than the largest per-capita rate.
#'
#'   The messages this function prints say the same thing — a converged run
#'   whose biomasses are still moving reports the drift and adds "Reduce `tol`
#'   to converge further." — but they are suppressed by `info_level = 0`, so in
#'   a script the `"convergence"` attribute is the reliable check.
#'
#'   Finally, a genuine fixed point need not be a *stable* one. Use
#'   [getStability()] to find out, and `solver = "newton"` to converge onto a
#'   fixed point that the dynamics themselves would run away from.
