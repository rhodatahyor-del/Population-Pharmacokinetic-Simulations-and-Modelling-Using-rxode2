
############################################################
# Population Pharmacokinetic Simulations using rxode2
#
# This script implements a series of population
# pharmacokinetic (PopPK) simulations using ODE-based
# models defined in rxode2.
#
# The simulations illustrate the effects of:
# - route of administration (oral vs intravenous),
# - model structure (one- vs two-compartment),
# - inter-individual variability (IIV)
# on concentration–time profiles.

rm(list = ls(all.names = TRUE))  # Remove all objects, including hidden ones
gc()    

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(purrr)
library(rxode2)
library(lotri)

options(scipen = 999)

global_font      <- "Segoe UI"
themeing_pattern <- theme(legend.position       = c(.65, .08),
                          legend.title          = element_blank(),
                          legend.text           = element_text(size = 16),
                          legend.key            = element_rect(fill = "transparent", colour = NA_character_),
                          text                  = element_text(size = 14, family = global_font, colour = "black"),
                          plot.title            = element_text(hjust = 0.5, size = 20, face = "bold"),
                          axis.text             = element_text(size = 16, colour = "black"),
                          axis.title            = element_text(size = 20),
                          axis.ticks            = element_line(colour = "black"),
                          panel.background      = element_rect(fill = "transparent", colour = "black"),
                          panel.border          = element_rect(colour = "black"),
                          plot.background       = element_rect(fill = "transparent", colour = NA_character_),
                          legend.background     = element_rect(fill = "transparent"),
                          legend.box.background = element_rect(fill = "transparent", colour = NA_character_),
                          panel.grid.major      = element_blank(),
                          panel.grid.minor      = element_blank())

# Simulations with rxODE

# Model 1 (one compartment)
############################################################
# One-compartment oral PK model without inter-individual variability.
# Purpose:
# To establish a baseline deterministic concentration–time
# profile assuming identical pharmacokinetic parameters across individuals.

my_model_1_cmpt_no_iiv <- rxode2({

  Bioav      = F_pop;
  ka         = ka_pop;
  Cl         = Cl_pop;
  V          = V_pop;

  f(Ad)      =   Bioav;
  d/dt(Ad)   = - ka*Ad;
  d/dt(Ac)   =   ka*Ad - (Cl/V)*Ac; 
  
  Cc         =   Ac/V;
  d/dt(AUC)  =   Cc;

})

# closed form (?; untested)
# my_model_1_cmpt_no_iiv_closed_form <- rxode2({

#   # Parameters (no IIV)
#   Bioav = F_pop;
#   ka    = ka_pop;
#   Cl    = Cl_pop;
#   V     = V_pop;

#   # Convenience
#   k = Cl / V;              # elimination rate

#   # Bioavailability applied to depot dosing
#   f(Ad) = Bioav;

#   # ---- Closed-form expressions (single oral dose at t=0) ----
#   # Notes:
#   # 1) These expressions assume a single dose placed into Ad at t = 0 (handled by your event table),
#   #    and bioavailability via f(Ad) = Bioav so that the effective depot "initial amount" is Bioav * Dose.
#   # 2) If you have multiple doses, use the same formulas per-dose and superpose (sum over doses),
#   #    with t replaced by (t - t_dose) for each dose time.

#   # Effective absorbed dose for a single dose at t=0
#   FD = Bioav * Dose;

#   # Depot amount
#   Ad = FD * exp(-ka * t);

#   # Central amount, concentration, and AUC
#   if (abs(ka - k) > 1e-12) {
#     Ac  = (ka * FD / (ka - k)) * (exp(-k * t) - exp(-ka * t));
#     Cc  = Ac / V;

#     AUC = (ka * FD / (V * (ka - k))) *
#       ( (1 - exp(-k * t)) / k - (1 - exp(-ka * t)) / ka );
#   } else {
#     # Degenerate limit ka == k
#     Ac  = FD * k * t * exp(-k * t);
#     Cc  = Ac / V;

#     AUC = (FD / (V * k)) * (1 - exp(-k * t) * (1 + k * t));
#   }

# })

# thetas
theta_no_iiv             <- c(F_pop  = 1, 
                              ka_pop = 0.45, 
                              Cl_pop = 7.5, 
                              V_pop  = 25)

# event function
event_function           <- function(dose, 
                                     sim_len) {dose_times <- 0
                                               event      <- et(amountUnits = "mg", timeUnits = "hours") %>% et(time = dose_times,
                                                                                                                amt  = as.numeric(dose),
                                                                                                                cmt  = "Ad") %>% et(time = seq(0, sim_len, by = 0.1))
                                              } 

event_no_iiv              <- event_function(100, 48) # 100 mg dose
                                                                            
simulation_rxode_no_iiv   <- rxSolve(my_model_1_cmpt_no_iiv,
                                     params   = theta_no_iiv,
                                     events   = event_no_iiv,
                                     method   = "dop853",   # <- explicit Runge–Kutta solver
                                     maxsteps = 200000,
                                     rtol     = 1e-6,
                                     atol     = 1e-8)

sim_rxode_dataframe_no_iiv <- as.data.frame(simulation_rxode_no_iiv)

# Plot 1
themeing_pattern <- theme(panel.grid.major = element_line(linetype = "dashed"), panel.grid.minor = element_blank(), legend.position = "none") #chatgptgenerated
data_plot_1                <- ggplot() +
                              geom_line(data = sim_rxode_dataframe_no_iiv, aes(x = as.numeric(time), y = as.numeric(Cc))) +
                              geom_point(data = sim_rxode_dataframe_no_iiv, aes(x = as.numeric(time), y = as.numeric(Cc)), size = 1.5) +
                              scale_x_continuous(breaks = seq(0, 48, by = 2), labels = seq(0, 48, by = 2)) +
                              scale_y_continuous(breaks = seq(0, 4, by = 1), labels = scales::number_format(accuracy = 1)) +
                              ggtitle("Concentration vs Time plot (1 CPT)") +
                              labs(x = "Time",
                                   y = "Concentration, Drug (Oral dose, No Variability)") +
                              coord_cartesian(xlim = c(0, 12), ylim = c(0, 4)) +
                              theme_light() + themeing_pattern +
                              theme(panel.grid.major = element_line(linetype = "dashed"), legend.position = "none")
data_plot_1



# Model 2
############################################################
# One-compartment oral PK model with inter-individual
# variability (IIV) on F, ka, CL, and V.
#
# Purpose:
# To illustrate the impact of log-normal parameter
# variability on population concentration–time profiles.
############################################################

my_model_1_cmpt_with_iiv <- rxode2({

  # IIV
  Bioav    = F_pop  * exp(omega_F);
  ka       = ka_pop * exp(omega_ka);
  Cl       = Cl_pop * exp(omega_Cl);
  V        = V_pop  * exp(omega_V);

  # Bioavailability
  f(Ad)      =   Bioav;
  d/dt(Ad)   = - ka*Ad;
  d/dt(Ac)   =   ka*Ad - (Cl/V)*Ac; 
  
  Cc         =   Ac/V;
  d/dt(AUC)  =   Ac/V;

})

# thetas
theta_with_iiv             <- c(F_pop  = 1, 
                                ka_pop = 0.45, 
                                Cl_pop = 7.5, 
                                V_pop  = 25)

# IIV matrix
omega_with_iiv             <- lotri(omega_F   ~  0.1,
                                    omega_ka  ~  0.1,  
                                    omega_Cl  ~  0.1,
                                    omega_V   ~  0.1)

# evid QD                        
event_function_with_iiv    <- function(dose_mg, 
                                       sim_len) {dose_times <- 0
                                                 pre_event  <- et(amountUnits = "mg", timeUnits = "hours") %>% et(time = dose_times,
                                                                                                                  amt  = as.numeric(dose_mg),
                                                                                                                  cmt  = "Ad") %>% et(time = seq(0, sim_len, by = 0.1))

                                                 event <- pre_event %>%
                                                   et(id = seq(1, 10, by = 1)) %>%
                                                   as_tibble() %>%
                                                   mutate(
                                                     time = as.numeric(time),
                                                     amt  = as.numeric(amt)
                                                   ) 
                                                 event
                                                            } 

event_with_iiv             <- event_function_with_iiv(100, 48) # 100 mg dose

simulation_rxode_with_iiv  <- rxSolve(my_model_1_cmpt_with_iiv,
                                      params   = theta_with_iiv,
                                      events   = event_with_iiv,
                                      omega    = omega_with_iiv,
                                      method   = "dop853",   # <- explicit Runge–Kutta solver
                                      maxsteps = 200000,
                                      rtol     = 1e-6,
                                      atol     = 1e-8)

sim_rxode_dataframe_with_iiv <- as.data.frame(simulation_rxode_with_iiv)

# Plot 2
data_plot_2                <- ggplot() +
                              geom_line(data = sim_rxode_dataframe_with_iiv, aes(x = as.numeric(time), y = as.numeric(Cc), colour = as.factor(id))) +
                              geom_point(data = sim_rxode_dataframe_with_iiv, aes(x = as.numeric(time), y = as.numeric(Cc), colour = as.factor(id)), size = 1.5) +
                              scale_x_continuous(breaks = seq(0, 48, by = 2), labels = seq(0, 48, by = 2)) +
                              scale_y_continuous(breaks = seq(0, 4, by = 1), labels = scales::number_format(accuracy = 1)) +
                              ggtitle("Concentration vs Time plot (1 CPT)") +
                              labs(x = "Time",
                                   y = "Concentration, Drug (Oral dose, Variability)") +
                              coord_cartesian(xlim = c(0, 12), ylim = c(0, 4)) +
                              theme_light() + themeing_pattern +
                              theme(panel.grid.major = element_line(linetype = "dashed"), legend.position = "none")
data_plot_2

# Model 3 (IV bolus)

############################################################
# One-compartment IV bolus PK model with inter-individual
# variability (IIV) on clearance and volume of distribution.
#
# Purpose:
# To isolate the impact of systemic disposition variability
# on concentration–time profiles in the absence of
# absorption phase
############################################################


my_model_1_cmpt_with_iiv_iv <- rxode2({

  # IIV
  Cl         = Cl_pop * exp(omega_Cl);
  V          = V_pop  * exp(omega_V);

  d/dt(Ac)   = - (Cl/V)*Ac; 
  
  Cc         =   Ac/V;
  d/dt(AUC)  =   Ac/V;

})

# thetas
theta_with_iiv_iv          <- c(Cl_pop = 7.5, 
                                V_pop = 25)

# IIV matrix
omega_with_iiv_iv          <- lotri(omega_Cl  ~  0.1,
                                    omega_V   ~  0.1)

# evid QD                        
event_function_with_iiv_iv <- function(dose_mg, 
                                       sim_len) {dose_times <- 0
                                                 pre_event  <- et(amountUnits = "mg", timeUnits = "hours") %>% et(time = dose_times,
                                                                                                                  amt  = as.numeric(dose_mg),
                                                                                                                  cmt  = "Ac") %>% et(time = seq(0, sim_len, by = 0.1))

                                                 event <- pre_event %>%
                                                   et(id = seq(1, 10, by = 1)) %>%
                                                   as_tibble() %>%
                                                   mutate(
                                                     time = as.numeric(time),
                                                     amt  = as.numeric(amt)
                                                   ) 
                                                 event
                                                            } 

event_with_iiv_iv            <- event_function_with_iiv_iv(100, 48) # 100 mg dose

simulation_rxode_with_iiv_iv <- rxSolve(my_model_1_cmpt_with_iiv_iv,
                                        params   = theta_with_iiv_iv,
                                        events   = event_with_iiv_iv,
                                        omega    = omega_with_iiv_iv,
                                        method   = "dop853",   # <- explicit Runge–Kutta solver
                                        maxsteps = 200000,
                                        rtol     = 1e-6,
                                        atol     = 1e-8)

sim_rxode_dataframe_with_iiv_iv <- as.data.frame(simulation_rxode_with_iiv_iv)

# Plot 3
data_plot_3                <- ggplot() +
                              geom_line(data = sim_rxode_dataframe_with_iiv_iv, aes(x = as.numeric(time), y = as.numeric(Cc), colour = as.factor(id))) +
                              geom_point(data = sim_rxode_dataframe_with_iiv_iv, aes(x = as.numeric(time), y = as.numeric(Cc), colour = as.factor(id)), size = 1.5) +
                              scale_x_continuous(breaks = seq(0, 48, by = 2), labels = seq(0, 48, by = 2)) +
                              scale_y_continuous(breaks = seq(0, 10, by = 2), labels = scales::number_format(accuracy = 1)) +
                              ggtitle("Concentration vs Time plot (1 CPT)") +
                              labs(x = "Time",
                                   y = "Concentration, Drug (IV dose, Variability)") +
                              coord_cartesian(xlim = c(0, 12), ylim = c(0, 10)) +
                              theme_light() + themeing_pattern +
                              theme(panel.grid.major = element_line(linetype = "dashed"), legend.position = "none")
data_plot_3

# Model 4, Oral dose, 2 CPT
############################################################
# Two-compartment oral PK model with inter-individual
# variability (IIV) on absorption, distribution, and
# elimination parameters.
#
# Purpose:
# To capture the combined effects of absorption kinetics
# and peripheral tissue distribution on population
# concentration–time profiles.
############################################################

my_model_2_cmpt_with_iiv <- rxode2({

  # IIV
  Bioav    = F_pop  * exp(omega_F);
  ka       = ka_pop * exp(omega_ka);
  Cl       = Cl_pop * exp(omega_Cl);
  V        = V_pop  * exp(omega_V);
  Q        = Q_pop  * exp(omega_Q);
  Vp       = Vp_pop * exp(omega_Vp);

  # Bioavailability
  f(Ad)      =   Bioav;
  d/dt(Ad)   = - ka*Ad;
  d/dt(Ac)   =   ka*Ad - (Cl/V)*Ac - (Q/V)*Ac + (Q/Vp)*Ap;
  d/dt(Ap)   =                       (Q/V)*Ac - (Q/Vp)*Ap;
  
  Cc         =   Ac/V;
  d/dt(AUC)  =   Ac/V;

})

# thetas
theta_with_iiv_2cpt        <- c(F_pop  = 1, 
                                ka_pop = 0.75, 
                                Cl_pop = 14.5, 
                                V_pop  = 15,
                                Q_pop  = 10,
                                Vp_pop = 115)

# IIV matrix
omega_with_iiv_2cpt        <- lotri(omega_F   ~  0.1,
                                    omega_ka  ~  0.1,  
                                    omega_Cl  ~  0.1,
                                    omega_V   ~  0.1,
                                    omega_Q   ~  0.1,
                                    omega_Vp  ~  0.1)

# evid QD                        
event_function_with_iiv_2cpt <- function(dose_mg, 
                                         sim_len) {dose_times <- 0
                                                   pre_event  <- et(amountUnits = "mg", timeUnits = "hours") %>% et(time = dose_times,
                                                                                                                    amt  = as.numeric(dose_mg),
                                                                                                                    cmt  = "Ad") %>% et(time = seq(0, sim_len, by = 0.1))

                                                  event <- pre_event %>%
                                                    et(id = seq(1, 10, by = 1)) %>%
                                                    as_tibble() %>%
                                                    mutate(
                                                      time = as.numeric(time),
                                                      amt  = as.numeric(amt)
                                                    ) 
                                                  event
                                                              } 

event_with_iiv_2cpt            <- event_function_with_iiv_2cpt(100, 48) # 100 mg dose

simulation_rxode_with_iiv_2cpt <- rxSolve(my_model_2_cmpt_with_iiv,
                                          params   = theta_with_iiv_2cpt,
                                          events   = event_with_iiv_2cpt,
                                          omega    = omega_with_iiv_2cpt,
                                          method   = "dop853",   # <- explicit Runge–Kutta solver
                                          maxsteps = 200000,
                                          rtol     = 1e-6,
                                          atol     = 1e-8)

sim_rxode_dataframe_with_iiv_2cpt <- as.data.frame(simulation_rxode_with_iiv_2cpt)

# Plot 4
data_plot_4                <- ggplot() +
                              geom_line(data = sim_rxode_dataframe_with_iiv_2cpt, aes(x = as.numeric(time), y = as.numeric(Cc), colour = as.factor(id))) +
                              geom_point(data = sim_rxode_dataframe_with_iiv_2cpt, aes(x = as.numeric(time), y = as.numeric(Cc), colour = as.factor(id)), size = 1.5) +
                              scale_x_continuous(breaks = seq(0, 48, by = 2), labels = seq(0, 48, by = 2)) +
                              scale_y_continuous(breaks = seq(0, 4, by = 1), labels = scales::number_format(accuracy = 1)) +
                              ggtitle("Concentration vs Time plot (2 CPT)") +
                              labs(x = "Time",
                                   y = "Concentration, Drug (IV dose, Variability)") +
                              coord_cartesian(xlim = c(0, 12), ylim = c(0, 4)) +
                              theme_light() + themeing_pattern +
                              theme(panel.grid.major = element_line(linetype = "dashed"), legend.position = "none")
data_plot_4

