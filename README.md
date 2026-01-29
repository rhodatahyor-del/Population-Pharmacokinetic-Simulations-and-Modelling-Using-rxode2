# Population-Pharmacokinetic-Simulations-and-Modelling-Using-rxode2

Population pharmacokinetic simulations in R using rxode2, illustrating one- and two-compartment oral and IV models with inter-individual variability

**Overview**

This repository contains a set of population pharmacokinetic (PopPK) simulations implemented using ODE-based models in rxode2. The simulations are designed to demonstrate how model structure, route of administration, and inter-individual variability (IIV) influence concentration–time profiles.

Rather than focusing on parameter estimation, this work emphasizes mechanistic interpretation of drug exposure and variability arising from pharmacokinetic processes.

**Background**

Many apparent failures in pharmacotherapy can be traced to limitations in drug exposure at the site of action, rather than insufficient intrinsic drug potency. Pharmacokinetic modeling provides a quantitative framework to examine how absorption, distribution, and elimination govern systemic drug concentrations.
Population PK models extend this framework by incorporating inter-individual variability, enabling exploration of heterogeneity in drug response under otherwise identical dosing conditions.

**Modeling Framework**

All models are implemented in rxode2, which numerically solves systems of ordinary differential equations (ODEs). Log-normal inter-individual variability is introduced using lotri variance–covariance structures.
Simulations are performed deterministically (no residual error) to isolate the effects of structural and population variability.

**Implemented Models**

****Model 1**: One-Compartment Oral (No IIV)**

•	First-order absorption and elimination

•	Identical PK parameters across individuals

•	Serves as a deterministic baseline

Purpose:

To establish a reference concentration–time profile without population variability.

**Model 2: One-Compartment Oral (With IIV)**

•	Log-normal variability on bioavailability, absorption rate, clearance, and volume

•	Multiple individuals simulated under identical dosing

Purpose:

To illustrate how inter-individual variability alone can generate divergent exposure profiles.

**Model 3: One-Compartment IV Bolus (With IIV)**

•	No absorption phase

•	Variability on clearance and volume of distribution only

Purpose:
To isolate the effect of systemic disposition variability in the absence of absorption-related processes.

**Model 4: Two-Compartment Oral (With IIV)**

•	First-order absorption

•	Central and peripheral compartments with inter-compartmental clearance

•	Variability on absorption, distribution, and elimination parameters

Purpose:
To capture the combined influence of absorption kinetics and tissue distribution on population exposure.

**Simulation Design**
•	Dose: Single 100 mg administration
•	Routes: Oral and intravenous
•	Simulation horizon: 48 hours
•	Output resolution: 0.1 h
•	Solver: Explicit Runge–Kutta (dop853)
•	Variability model: Log-normal IIV

**Outputs**
Each simulation produces:

•	Individual concentration–time profiles

•	Visual comparisons across model structures and routes of administration
Plots are generated using ggplot2 with consistent scaling to facilitate comparison.

#Model 1

![Concentration–Time Profile](outputs/model1_plot.png)

## Model 2

![Concentration–Time Profile](outputs/model2_plot.png)

#Model 3

![Concentration–Time Profile](outputs/model3_plot.png)

#Model 4

![Concentration–Time Profile](outputs/model4_plot.png)

**Intended Use**

This repository is intended for:

•	Educational demonstrations of PopPK concepts

•	Methodological illustration of variability effects

•	Conceptual support for drug delivery and exposure-driven discussions

It is not intended for clinical decision-making or parameter inference.

**Dependencies**

•	rxode2

•	lotri

•	ggplot2, dplyr, tidyr, purrr, readr

