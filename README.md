# qsp-model-steroid-afp
Quantitative systems pharmacology model of a glucocorticoid-binding antibody fusion protein

This is a quantitative systems pharmacology model simulating a glucocorticoid-binding antibody fusion protein in the human body. The drug is intended to create immunosuppressive effects in selected cell types or tissues, preventing damaging steroid side effects. The protein drug consists of a leukocyte-binding minibody and a domain that non-covalently binds a glucocorticoid. This model predicts the occupancy of the glucocorticoid receptor over time in theoretical on- and off-target tissues under various conditions. The R script uses dMod to organize and solve the system of mass-action-based differential equations. 

This model was used to determine several critical features to design into the protein drug. First, I found that the antibody’s target antigen must endocytose quickly upon protein binding. However, endosomal release of steroid or antibody using histidine switching is unlikely to be helpful. My model also showed that the fusion protein could direct endogenous cortisol to leukocytes, achieving immunosuppression without any synthetic glucocorticoid.

For a detailed description of the model, please refer to the paper [Pharmacokinetic Model of a Glucocorticoid-Binding Antibody Fusion Protein](thesis-chapter-pk-model-steroid-afp.pdf). This document covers the model setup, shows the predicted glucocorticoid receptor occupancy over time, and elaborates on the results. It also contains citations and derivations for compartment volumes and rate constants. This version of the model and accompanying paper was written in March 2019. 

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
