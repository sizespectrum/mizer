context("Wrapper functions for trait and community models")

test_that("scaling model is set up correctly", {
  p <- set_scaling_model()
  sim <- project(p, t_max=5, effort = 0)
  # Check that total biomass changes little (relatively)
  bm <- getBiomass(sim)
  expect_lt(abs(max(bm[1, ]-bm[6, ])), 2.6*10^(-6))
})

test_that("trait-based model multiple gears",{
    # Check multiple gears are working properly
    min_w_inf <- 10
    max_w_inf <- 1e5
    no_sp <- 10
    w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
    knife_edges <- w_inf * 0.05
    params <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf, max_w_inf = max_w_inf, knife_edge_size = knife_edges)
    expect_that(params@species_params$knife_edge_size, is_identical_to(knife_edges))
    # All gears fire
    sim1 <- project(params, t_max = 10, effort = 1)
    fmg <- getFMortGear(sim1)
    for (i in 1:no_sp){
        expect_that(all(fmg[10,1,i,params@w < knife_edges[i]] == 0), is_true())
        expect_that(all(fmg[10,1,i,params@w >= knife_edges[i]] == 1), is_true())
    }
    # Only the 4th gear fires
    params <- set_trait_model(no_sp = no_sp, min_w_inf = min_w_inf, max_w_inf = max_w_inf, knife_edge_size = knife_edges, gear_names = 1:no_sp)
    effort <- c(0,0,0,1,0,0,0,0,0,0)
    names(effort) = 1:no_sp
    sim2 <- project(params, t_max = 10, effort = effort)
    fmg <- getFMortGear(sim2)
    expect_that(all(fmg[10,c(1:3,5:10),c(1:3,5:10),] == 0), is_true())
        expect_that(all(fmg[10,4,4,params@w < knife_edges[4]] == 0), is_true())
        expect_that(all(fmg[10,4,4,params@w >= knife_edges[4]] == 1), is_true())

})
