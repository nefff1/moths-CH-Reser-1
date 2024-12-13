# Five decades of trapping reveal effects of season, weather and underlying environmental variables on moth abundance, richness and biomass

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14439275.svg)](https://doi.org/10.5281/zenodo.14439275)

This repository contains codes and some data that were used in the analyses for the manuscript:

Neff F, Chittaro Y, Korner-Nievergelt F, Litsios G, Rey E, Knop E. **Five decades of trapping reveal effects of season, weather and underlying environmental variables on moth abundance, richness and biomass.**

The study is based on a vast moth community dataset from Switzerland collected by Dr. Ladislaus Rezbanyai-Reser.

The following R files are included in the folder *R_Code*:

-   **R_prepare_data.R**: File used to prepare all data frames used for the analyses.

-   **R_analyses.R**: Code used for the analyses and creation of outputs.

The following Stan code files are included in the folder *Stan_Code*:

-   **Stan_hg_spline_s2p1_r4.stan**: Stan model code of a regression model with 2 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a hurdle gamma distribution.

-   **Stan_nb_spline_s2p1_r4.stan**: Stan model code of a regression model with 2 smoothing term, an additional smoothing term with restricted frame (hours active), 4 random terms and a zero-inflated negative binomial distribution.

Besides data available from other sources indicated in the code (e.g. GBIF), the folder *Data* contains:

-   **d_mass.txt**: Estimated species-level biomass used to estimate community-wide biomass. Based on a set of allometric relationships.
-   **d_nullnights.txt**: List of nights in which nothing was caught (missing from GBIF dataset).
-   **d_overwintering_stage.txt**: Overwintering stages for the species contained in the studied dataset.
-   **d_samplings.txt**: Details on sampling site-year combinations used in the analyses: Spatio-temporal clusters, sampling pairs, different land-use proportions.
-   **d_taxonomy.txt**: Moth species names according to the taxonomy used in the current analyses. Can be joint to GBIF data through the *taxonID* variable.
