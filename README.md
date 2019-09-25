# RFlux

An R graphical user interface for processing eddy covariance raw data and release high quality fluxes of the main GHGs exchanged by ecosystems and agricultural fields. Fluxes are estimated through a call to the open source EddyPro\textregistered software (LI-COR, Biosciences, 2019). 'RFlux' provides tools for the metadata management as well as for the implementation of the robust data cleaning procedure described by Vitale et al (2019, https://doi.org/10.5194/bg-2019-270).

'RFlux' has been developed in the context of the ICOS Ecosystem Thematic Centre. The authors thanks the ENVRIPLUS H2020 European project (Grant Agreement 654182) and the ENVRIFAIR H2020 European project (Grant Agreement 824068) for the support.

To install the RFlux package
```{r, eval = F}
# install.packages("devtools")
devtools::install_github("domvit81/RFlux")
```


