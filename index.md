[![CRAN Status](https://www.r-pkg.org/badges/version/mizer)](https://cran.r-project.org/package=mizer)
[![Travis-CI Build Status](https://travis-ci.org/sizespectrum/mizer.svg?branch=master)](https://travis-ci.org/sizespectrum/mizer)
[![Coverage status](https://codecov.io/gh/gustavdelius/mizer/branch/master/graph/badge.svg)](https://codecov.io/github/gustavdelius/mizer?branch=master)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mizer)](https://cran.r-project.org/package=mizer)
[![Rdoc](http://www.rdocumentation.org/badges/version/mizer)](http://www.rdocumentation.org/packages/mizer)

# mizer
Mizer is an R package to run multi-species size-spectrum models of fish
communities. The package has been developed to model marine ecosystems that are
subject to fishing. However, it may also be appropriate for other ecosystems.

The package contains routines and methods to allow users to set up an ecosystem 
mode, and then project it through time under different fishing strategies.
Methods are included to explore the results, including plots and calculation of
community indicators such as the slope of the size spectrum. Size-based models
can be complicated so mizer contains many default options that can be easily
changed by the user.

<iframe width="560" height="315" src="https://www.youtube.com/embed/0RlXqLbFbWc"
frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; 
picture-in-picture" allowfullscreen></iframe>

# Learning resources
No matter whether you are a new to ecological modelling, a fisheries scientist, a policy maker, a mathematician, or someone interested in further developing mizer, you will find useful resources on this website. The different learning resources are further described in the [Getting started page](articles/mizer.html), although we also also discuss them below.

For those new to mizer, perhaps the easiest way to start learning about mizer is to watch the [online video tutorials,](https://www.youtube.com/watch?v=zh0PDyTUssw&list=PLCTMeyjMKRkqR7uohI3p-61P7ZJj8sd5B) which explain how to make use of the different features that the mizer software provides. A gentle way to get into actually using mizer is to try out the [online mizer simulation app.](https://mizer.shinyapps.io/selectivity/) This app allows you to experiment with running mizer simulations under different conditions, and investigate the outcome of changing different fishing policies on the system, without having to learn how to use R. 

 This website holds complete documentation about the different aspects of mizer. The articles on [the general model description](articles/model_description.html), [the community model](articles/community_model.html), [the trait-based model](articles/trait_model.html), [the multispecies model](articles/multispecies_model.html), on [running](articles/running_a_simulation.html) and [exploring](articles/exploring_the_simulation_results.html) simulations, on [the North sea model](articles/a_multispecies_model_of_the_north_sea.html) and on [the scale-invariant model.](articles/scale_invariant_trait_based_model.html) These articles provide a firm background in using mizer within R. We recommend reading these articles to users who wish to enjoy the full functionality of mizer by using it in R. Although you should feel free to skip around to whichever articles suit you best. For example, if you are not mathematically inclined you may wish to skip the article on the model description, and if you want to get straight into setting up your own ecosystem model, you may want to go straight to the article on [the multispecies model.](articles/multispecies_model.html)

If you are mathematically inclined, you may wish to read the article on [the model description](articles/model_description.html), and the article on [Mathematical Details behind Mizer.](articles/mathematical_details.html) People wishing to develop and expand the mizer software further will find [the Developer Guide](articles/developer_vignette.html) helpful.

In addition to this documentation there are [help pages](reference/index.html) explaining the details of implementing all of the functions and methods of mizer.

Mizer can also be used to create web apps that allow users to explore models
without the need to install R. An [example of such an
app](https://mizer.shinyapps.io/selectivity/) investigates the effect of
switching to a gear with a T90 extension net to reduce the catches of undersized
hake and red mullet

Mizer is still under active development, currently funded by the European
Commission Horizon 2020 Research and Innovation Programme under Grant Agreement
No. 634495 for the project MINOUW (http://minouw-project.eu/) and the Australian
Research Council Discovery Project ("Rewiring Marine Food Webs").

Does your project or publication use mizer? If so, we'd love to know. You can
also join our Google Discussion group here:
https://groups.google.com/forum/#!forum/size-spectrum-models

The package is on [CRAN](https://cran.r-project.org/package=mizer) and therefore
available from R's build-it package manager.
See the accompanying
[vignette](https://cran.r-project.org/web/packages/mizer/vignettes/mizer_vignette.pdf)
for more details on how the package works, including detailed examples.


