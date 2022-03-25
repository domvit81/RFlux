# RFlux

An R graphical user interface for processing eddy covariance raw data and release high quality fluxes of the main GHGs exchanged by ecosystems and agricultural fields. Fluxes are estimated through a call to the open source EddyPro software (registered trademark, LI-COR, Biosciences, 2021). 'RFlux' provides tools for the metadata management as well as for the implementation of the robust data cleaning procedure described by Vitale et al (2020, https://doi.org/10.5194/bg-17-1367-2020).

'RFlux' has been developed in the context of the ICOS Ecosystem Thematic Centre. The authors thanks the ENVRIPLUS H2020 European project (Grant Agreement 654182) and the ENVRIFAIR H2020 European project (Grant Agreement 824068) for the support.

To install the RFlux package
```{r, eval = F}
install.packages("devtools")
devtools::install_github("icos-etc/RFlux")
help(package=RFlux)
citation("RFlux")

Domenico Vitale, Dario Papale and ICOS-ETC Team (2021)
RFlux: An R Package for Processing and Cleaning Eddy Covariance Flux Measurements
R package version 3.2.0

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {{RFlux}: An R Package for Processing and Cleaning Eddy Covariance Flux Measurements},
    author = {Domenico Vitale and Dario Papale and {ICOS-ETC Team}},
    organization = {ICOS-ETC},
    address = {Viterbo, Italy},
    year = {2022},
    note = {R package version 3.2.0},
    url = {https://github.com/icos-etc/RFlux},
  }

```


