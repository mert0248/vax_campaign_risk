
# --------------------------------------------------------------------------------------------------
# estimate excess infection risk due to vaccination campaigns during covid-19 pandemic
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# load libraries
# --------------------------------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(gridExtra)

options(dplyr.summarise.inform = FALSE)

#--------------------------------------------------------------------------------------------------
# Simulate epidemics using CovidM (must have CovidM installed)
#--------------------------------------------------------------------------------------------------

print("Running CovidM scenarios")

# CovidM options
cm_path = "~/Documents/GitHub/covidm/"; ### TO UPDATE: path to CovidM
cm_force_rebuild = F;
cm_build_verbose = T;
cm_version = 1;
source(paste0(cm_path, "/R/covidm.R"))

# Load data on age-varying symptomatic rate
covid_scenario = qread(paste0( "./data/2-linelist_symp_fit_fIa0.5.qs"));
covy = unname(unlist(covid_scenario[sample.int(nrow(covid_scenario), 1), f_00:f_70]));
covy = rep(covy, each = 2);

# start and end dates for simulations and for intervention
model_start = "2020-01-01"
model_end = "2022-01-01"
iv_start = "2020-02-01"
iv_end = "2022-01-01"

# countries to model
countries <- c("Burkina Faso","Ethiopia","Brazil")

# interventions - reduces home, work, school and other contacts
interventions <- list(soc_dist_0.4=c(1,0.4,0.4,0.4),
                      soc_dist_0.6=c(1,0.6,0.6,0.6),
                      soc_dist_0.7=c(1,0.7,0.7,0.7),
                      soc_dist_0.8=c(1,0.8,0.8,0.8))
# R0 assumptions
r0s <- c(1.8,2.0,2.2)

# somewhere to keep model outputs
m <- list() 

# run simulations
for (country in countries) {
  for (R0 in r0s){
    for (i in 1:length(interventions)){
      
      # build parameters
      params = cm_parameters_SEI3R(country, deterministic = T, date_start = model_start, date_end = model_end)
      params$pop[[1]]$y = covy
      params$pop[[1]]$fIa = rep(1,16)
      
      # calculate u to give target R0
      current_R0 = cm_calc_R0(params, 1) 
      target_R0 = R0
      params$pop[[1]]$u = params$pop[[1]]$u * target_R0 / current_R0
      
      # add impact of interventions
      iv_params <- interventions[[i]]
      iv = cm_iv_build(params) # this sets up a data structure for doing interventions
      cm_iv_contact(iv, iv_start, iv_end, iv_params) # changes contact matrix components: home, work, school and other
      params = cm_iv_apply(params, iv) # sets the "schedule" parameter to follow interventions in iv.
      
      # Seed outbreak with infections equal to 1 in 10^5 of total population seeded
      # over 7 days across age groups 20 to 50 
      n_to_seed = floor(sum(params$pop[[1]]$size)/(6*7*10^5)) 
      params$pop[[1]]$seed_times = rep(0:6, each = n_to_seed) # seed over 7 days
      params$pop[[1]]$dist_seed_ages = cm_age_coefficients(20, 50, 5 * (0:length(params$pop[[1]]$size))) # infections in ages 20-50
      
      # run the model
      run = cm_simulate(params, 1)
      
      # Collapse age-structure and reformat data
      seir_tot <- run$dynamics %>%
        filter(run == 1) %>%
        filter(compartment %in% c("S","E","Ia","Is","Ip","R")) %>%
        group_by(t, compartment) %>%
        summarise(n=sum(value)) %>%
        spread(compartment,n) %>%
        mutate(I=Ia+Is+Ip) %>%
        mutate(N=S+E+I+R)
      
      # Normalised output
      seir_norm <- seir_tot %>%
        mutate_at(vars(-t),list(~./N))
      
      # Add incidence based on change in S
      seir_norm$inc = 0
      for (k in 2:length(seir_norm$S)){
        seir_norm$inc[k] = seir_norm$S[k-1]-seir_norm$S[k]
      }
      
      # Store results
      n <- paste0(country,"-R0_",R0,"-",names(interventions[i]))
      
      m[[n]]$params       <- params   
      m[[n]]$seir_norm    <- seir_norm
      m[[n]]$seir_name    <- n
      m[[n]]$country      <- country
      m[[n]]$r0           <- R0
      m[[n]]$intervention <- names(interventions[i])
    }
  }
}

# --------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# covid sia risk model 3 - linking increased infection to vaccinators and community members
# --------------------------------------------------------------------------------------------------

sia_risk_3 <- function(t_c,           # time index of first day of campaign
                       S,             # vector of susceptibles in community
                       I,             # vector of (all) infecteds in community
                       Is,            # vector of symptomatic infecteds in community
                       Ia,            # vector of asymptomatic infecteds in community
                       Ip,            # vector of presymptomatic infecteds in community
                       S_v,           # vector of susceptibles among vaccinators
                       I_v,           # vector of (all) infecteds among vaccinators
                       Is_v,          # vector of symptomatic infecteds among vaccinators
                       Ia_v,          # vector of asymptomatic infecteds among vaccinators
                       Ip_v,          # vector of pre-symptomatic infecteds among vaccinators
                       n_cj,          # no. community contacts by child/caregiver during vaccine visit
                       n_cv,          # no. contacts with vaccinators by child during vaccine visit
                       n_vj,          # no. additional community contacts by vaccinator per hh visited
                       u,             # transmissibility per contact
                       d_c,           # duration of campaign (days)
                       v_per_d,       # no. of children vaccinated per day by vaccinator
                       ppe_eff_v,     # effect of ppe on transmission to vaccinator
                       ppe_eff_c,     # effect of ppe on transmission to vaccinee
                       d_e,           # latent period before infectious
                       d_r,           # duration from infectious until recovered/removed
                       scrn_sympt = F # whether screening of symptomatics is employed
){
  
  # somewhere to store results for each campaign start date
  
  ind_risk_t <- c()       # risk to susceptible individuals if campaign starts at time t
  ind_ex_risk_t <- c()    # excess risk to individuals accounting for proportion susceptible
  vac_risk_t <- c()       # risk to susceptible vaccinators if campaign starts at time t
  vac_ex_risk_t <- c()    # excess risk to vaccinators accounting for proportion susceptible
  
  for (t in t_c[1:(length(t_c)-d_c)]){  # loop of days of the epidemic on which campaign could start
    
    # somewhere to store results by day of campaign
    ind_risk_d <- c()
    ind_ex_risk_d <- c()
    vac_risk_d <- c()
    vac_ex_risk_d <- c()
    
    for (d in 1:d_c){ # for each start day loop over days of the campaign
      
      # ***risk that a susceptible vaccinator is infected ON day d***
    
      if (scrn_sympt) { # assume child/caregiver symptomatics are removed so prevalence is Ia + Ip
        vac_risk_d[d] <- (
          1 - ( 
            (1 - (u * (1 - ppe_eff_v) * (Ip[[t+d]] + Ia[[t+d]])) )^(2*v_per_d) * # contacts with child/caregiver
            (1 - (u * (1 - ppe_eff_v) * I[[t+d]]))^(n_vj*v_per_d) # community contacts
          )
        )
      } else { # assume community prevalence is I
        vac_risk_d[d] <- (
          1 - (
            (1 - (u * (1 - ppe_eff_v) * I[[t+d]]) )^(2*v_per_d) * # contacts with child/caregiver
            (1 - (u * (1 - ppe_eff_v) * I[[t+d]]) )^(n_vj*v_per_d) # community contacts
          )
        )
      }
      
      # cumulative excess risk on day d accounting for proportion susceptible
      vac_ex_risk_d[d] <- (1 - prod(1 - vac_risk_d[1:d])) * S_v[[t+d]] 
      
      # ***risk that child/caregiver being vaccinated is infected on day d***
      
      # calculate prevalence prev_d among vaccinators on day d of campaign
      if (scrn_sympt) { # symptomatic screening
        
        # among vaccinators who become infected during vaccine delivery assume 
        # that proportion of infectious who are asymptomatic is the same as 
        # Ia[t]/I[t] and proportion pre-symptomatic is Ip[t]/I[t]
        
        if (d <= d_e) { # no increase in infectious until after latent period d_e
          prev_d <- Ia_v[[t+d]] + Ip_v[[t+d]] 
        } else if (d <= (d_e + d_r)) { # increased infectious but none yet recovered
          prev_d <- Ia_v[[t+d]] + Ip_v[[t+d]] + (vac_ex_risk_d[[d-d_e]] * (Ia[[t+d]]+Ip[[t+d]])/I[[t+d]])
        } else { # increased infectious, some recovered
          prev_d <- Ia_v[[t+d]] + Ip_v[[t+d]] + ((vac_ex_risk_d[[d-d_e]] - vac_ex_risk_d[[d-(d_e+d_r)]]) * (Ia[[t+d]]+Ip[[t+d]])/I[[t+d]])
        }
        ind_risk_d[d] <- (
          1 - (
            # risk from contact with vaccinators
            (1 - (u * (1 - ppe_eff_c) * prev_d))^(2*n_cv) *
            # risk from contact with community members
            (1 - (u * I[[t+d]]))^(2*n_cj)
          )
        )
      } else { # no symptomatic screening
        if (d <= d_e) { # no increase in infectious until after latent period d_e
          prev_d <- I_v[[t+d]]
        } else if (d <= (d_e + d_r)) { # increased infectious but none yet recovered
          prev_d <- I_v[[t+d]] + vac_ex_risk_d[[d-d_e]] 
        } else { # increased infectious, some recovered
          prev_d <- I_v[[t+d]] + vac_ex_risk_d[[d-d_e]] - vac_ex_risk_d[[d-(d_e+d_r)]]
        }
        ind_risk_d[d] <- (
          1 - (
            # risk from contact with vaccinators
            (1 - (u * (1 - ppe_eff_c) * prev_d)) ^ (2*n_cv) *
            # risk from contact with community members
            (1 - (u * I[[t+d]]) ) ^ (2*n_cj)
          )
        )
      }
   
      # excess risk accounting for proportion susceptible
      ind_ex_risk_d[d] <- ind_risk_d[[d]] * S[[t+d]] 
    } 
    
    # vaccinators - store cumulative risk at end of campaign
    vac_risk_t[t+1]     <- vac_risk_d[[d_c]]
    vac_ex_risk_t[t+1]  <- vac_ex_risk_d[[d_c]]
    
    # individuals - store average risk across days of the campaign
    ind_risk_t[t+1]    <- mean(ind_risk_d)
    ind_ex_risk_t[t+1]  <- mean(ind_ex_risk_d)
    
  }
  
  # pad results for last d_c time-points for which results are undefined
  for (t in t_c[(length(t_c)-d_c+1):length(t_c)]){
    vac_risk_t[t+1]     <- NA
    vac_ex_risk_t[t+1]  <- NA
    ind_risk_t[t+1]     <- NA
    ind_ex_risk_t[t+1]  <- NA
  }
  
  # return results
  result <- data.table(
    t = t_c,
    ind_ex_risk = ind_ex_risk_t,
    vac_ex_risk = vac_ex_risk_t
    )
  
  return(result)
}
# --------------------------------------------------------------------------------------------------


#

# --------------------------------------------------------------------------------------------------
# Summarise some SEIR output
# --------------------------------------------------------------------------------------------------

# extract seir data for plotting each model run

seir_dat <- list()
for (mod in m){
  seir_dat <- rbind(seir_dat,
                    cbind(seir_name = mod$seir_name,
                          country = mod$country,
                          r0 = mod$r0,
                          intervention = mod$intervention,
                          as.data.frame(mod$seir_norm)
                          )
                    )
}
seir_dat <- as.data.table(seir_dat)

# --------------------------------------------------------------------------------------------------
# Run covid vax risk model scenarios for each modelled seir epidemic curve
# --------------------------------------------------------------------------------------------------

scenarios <- as.data.table(fread("./data/vax_scenarios.csv"))

sia_model_run = 0
total_runs <- length(countries) * max(scenarios[,scen_id]) * 5 
sia_dat <- list()

for (i in 1:length(countries)){
  for (s in 1:max(scenarios[,scen_id])){
    for (ppe_effect in c(0,0.5,0.75,0.9,1)){

      sia_model_run = sia_model_run + 1
      print(paste("Vaccine campaign model run:",sia_model_run, " of ", total_runs))

      name = paste0(countries[i], 
                       "-R0_", scenarios[s, r0],
                       "-", scenarios[s,soc_dist]
                       )
      v_name = paste0(countries[i], 
                       "-R0_", scenarios[s, r0],
                       "-", scenarios[s,v_soc_dist]
      )
      
      result <- sia_risk_3(
        t_c =  m[[name]]$seir_norm$t,
        S =    m[[name]]$seir_norm$S,
        I =    m[[name]]$seir_norm$I,
        Is =   m[[name]]$seir_norm$Is,
        Ia =   m[[name]]$seir_norm$Ia,
        Ip =   m[[name]]$seir_norm$Ip,
        S_v =  m[[v_name]]$seir_norm$S,
        I_v =  m[[v_name]]$seir_norm$I,
        Is_v = m[[v_name]]$seir_norm$Is,
        Ia_v = m[[v_name]]$seir_norm$Ia,
        Ip_v = m[[v_name]]$seir_norm$Ip,
        u = m[[name]]$params$pop[[1]]$u[[1]],
        d_c = scenarios[s, d_c],
        d_e = scenarios[s, d_e],
        d_r = scenarios[s, d_r],
        n_cv = scenarios[s, n_cv],
        n_cj = scenarios[s, n_cj],
        n_vj = scenarios[s, n_vj],
        v_per_d = scenarios[s, v_per_d],
        scrn_sympt = scenarios[s, scrn_sympt],
        ppe_eff_v = ppe_effect * scenarios[s, diff_ppe_v],
        ppe_eff_c = ppe_effect * scenarios[s, diff_ppe_c]
      )
      
      result[,scen_id:=s]
      result[,country:=eval(countries[i])]
      result[,ppe_effect:=ppe_effect]
      sia_dat[[sia_model_run]] <- merge(result,scenarios[s,],by="scen_id",all=T)
      
    }
  }
}

sia_dat <- rbindlist(sia_dat, fill=T)

sia_dat[, country := factor(country,levels=c("Burkina Faso","Ethiopia","Brazil"))]

# --------------------------------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------------------------------
max_sia_dat <- sia_dat[
  !is.na(vac_ex_risk) & !is.na(ind_ex_risk),
  .(vac_ex_risk=max(vac_ex_risk),ind_ex_risk=max(ind_ex_risk)),
  by=setdiff(names(sia_dat),c("t","vac_ex_risk","ind_ex_risk"))
]

base_dat <- max_sia_dat[
  scen_id==1, 
  .(country=country, 
    ppe_effect=ppe_effect, 
    base_ind_ex_risk=ind_ex_risk, 
    base_vac_ex_risk=vac_ex_risk)
]
max_sia_dat <- max_sia_dat[scen_id!=1,]
max_sia_dat <- merge(max_sia_dat,base_dat,by=c("country","ppe_effect"))

max_sia_dat <- dcast(
  max_sia_dat[sens_type %in% c("max","min")],
  name + country + ppe_effect + base_ind_ex_risk + base_vac_ex_risk ~ sens_type,
  value.var = c("ind_ex_risk","vac_ex_risk")
)
# replace NAs with baseline value if a scenario has only a min or max, but not both
max_sia_dat <- max_sia_dat[is.na(vac_ex_risk_min),vac_ex_risk_min:=base_vac_ex_risk]
max_sia_dat <- max_sia_dat[is.na(vac_ex_risk_max),vac_ex_risk_max:=base_vac_ex_risk]
max_sia_dat <- max_sia_dat[is.na(ind_ex_risk_min),ind_ex_risk_min:=base_ind_ex_risk]
max_sia_dat <- max_sia_dat[is.na(ind_ex_risk_max),ind_ex_risk_max:=base_ind_ex_risk]

max_sia_dat <- max_sia_dat[,name := factor(name, levels = rev(c(
  "r0",
  # 'soc_dist',
  "v_per_d",
  "d_c_fixed_vpd",
  "d_c_fixed_tot_v",
  "house-to-house",
  "v_hi_contact",
  "screening",
  "dif_ppe",
  "n_cj",
  "n_cv",
  "d_r",
  "d_e"
)))]


# --------------------------------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# fig 1 - epidemic curves by R0 and country
# --------------------------------------------------------------------------------------------------

country.lbl <- c('Brazil'='Brazil',
                'Ethiopia' = 'Ethiopia', 'Burkina Faso'='Demography: Burkina Faso')

seir_dat[, country := factor(country,levels=c("Burkina Faso","Ethiopia","Brazil"))]

# fig 1a
fig1a <- ggplot(data = seir_dat[intervention=="soc_dist_0.6"],) +
  geom_line(aes(x=t, y=inc, 
                color=as.factor(r0)
  )) +
  facet_grid(
    country~.,
    labeller = labeller(country = country.lbl)
  ) +
  ggtitle("A") +
  scale_color_viridis_d(option="plasma") +
  ylab("Incidence of infection") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  labs(color=expression(R[0])) +
  scale_x_continuous(
    expand = c(0,0),
    name = " ",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="top") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# fig 1b
fig1b <- ggplot(data = seir_dat[intervention=="soc_dist_0.6"],) +
  geom_line(aes(x=t, y=I, 
                color=as.factor(r0)
  )) +
  facet_grid(
    country~.,
    labeller = labeller(country = country.lbl)
  ) +
  ggtitle("B") +
  scale_color_viridis_d(option="plasma") +
  ylab("Prevalence of infection") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  labs(color=expression(R[0])) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="top") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# fig 1c
fig1c <- ggplot(data = seir_dat[intervention=="soc_dist_0.6"],) +
  geom_line(aes(x=t, y=R, 
                color=as.factor(r0)
  )) +
  facet_grid(
    country~.,
    labeller = labeller(country = country.lbl)
  ) +
  scale_color_viridis_d(option="plasma") +
  ggtitle("C") +
  ylab("Cumulative proportion infected") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  labs(color=expression(R[0])) +
  scale_x_continuous(
    expand = c(0,0),
    name = " ",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="top") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

leg <- get_legend(fig1a)
fig1a <- fig1a + theme(legend.position="none")
fig1b <- fig1b + theme(legend.position="none")
fig1c <- fig1c + theme(legend.position="none")

fig1 <- grid.arrange(fig1a, fig1b, fig1c, leg, ncol=3, nrow = 2,
             layout_matrix = rbind(c(1,2,3), c(4,4,4)),
             widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))

ggsave("./figures/fig1 - epi curves.png",fig1,width=7.5, height=6, units="in")

# --------------------------------------------------------------------------------------------------
# fig S3 - epidemic curves by % contact reduction (used for sensitivity analyses)
# --------------------------------------------------------------------------------------------------

# figS3 a
figS3a <- ggplot(data = seir_dat[r0==2],) +
  geom_line(aes(x=t, y=inc, 
                color=as.factor(intervention)
  )) +
  facet_grid(
    country~.,
    labeller = labeller(country = country.lbl)
  ) +
  ggtitle("A") +
  ylab("Incidence") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_color_viridis_d(
    option="plasma",
    # name="Vaccinator baseline \ncontact reduction",
    # labels = c("40%","30%","20%")
    name="Baseline contact reduction",
    labels = c("60%","40%","30%","20%")
  ) +
  scale_x_continuous(
    expand = c(0,0),
    name = " ",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="top") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# figS3 b
figS3b <- ggplot(data = seir_dat[r0==2],) +
  geom_line(aes(x=t, y=I, 
                color=as.factor(intervention)
  )) +
  facet_grid(
    country~.,
    labeller = labeller(country = country.lbl)
  ) +
  ggtitle("B") +
  ylab("Prevalence") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_color_viridis_d(
    option="plasma",
    # name="Vaccinator baseline \ncontact reduction",
    # labels = c("40%","30%","20%")
    name="Baseline contact reduction",
    labels = c("60%","40%","30%","20%")
  ) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="top") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# figS3 c
figS3c <- ggplot(data = seir_dat[r0==2],) +
  geom_line(aes(x=t, y=R, 
                color=as.factor(intervention)
  )) +
  facet_grid(
    country~.,
    labeller = labeller(country = country.lbl)
  ) +
  ggtitle("C") +
  ylab("Cumulative proportion infected") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_color_viridis_d(
    option="plasma",
    # name="Vaccinator baseline \ncontact reduction",
    # labels = c("40%","30%","20%")
    name="Baseline contact reduction",
    labels = c("60%","40%","30%","20%")
  ) +
  scale_x_continuous(
    expand = c(0,0),
    name = " ",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="top") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

leg <- get_legend(figS3a)
figS3a <- figS3a + theme(legend.position="none")
figS3b <- figS3b + theme(legend.position="none")
figS3c <- figS3c + theme(legend.position="none")

figS3 <- grid.arrange(figS3a, figS3b, figS3c, leg, ncol=3, nrow = 2,
                     layout_matrix = rbind(c(1,2,3), c(4,4,4)),
                     widths = c(2.7, 2.7, 2.7), heights = c(2.5, 0.2))

ggsave("./figures/figS3 - epi curves.png",figS3,width=7.5, height=6, units="in")

# --------------------------------------------------------------------------------------------------
# fig 2 - Excess risk by R0 and country
# --------------------------------------------------------------------------------------------------

country.lbl <- c('Brazil'='Brazil',
                 'Ethiopia' = 'Ethiopia', 'Burkina Faso'='Demography:\nBurkina Faso')

# fig 2a - HCW risk by R0 and country
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs <- c('1.8'="R0 = 1.8",
             '2'="R0 = 2.0",
             '2.2'="R0 = 2.2")

fig2a <- ggplot(data = sia_dat[
  name=="base" | name=="screening" | name=="r0" | name=="r0_screening",
]) +
  geom_line(aes(
    x=t,
    y=vac_ex_risk, 
    color=factor(ppe_effect),
    linetype=scrn_sympt
  )) +
  facet_grid(
    country ~ r0,
    labeller = labeller(
      r0 = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option="viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  scale_linetype(
    name = "Symptomatic \nscreening",
    labels = line.labs
  ) +
  ggtitle("A - Excess risk to vaccinators") +
  ylab("Excess infection risk to vaccinators") +
  scale_y_continuous(
    labels = function(x) sprintf("%.2g%%",x*100)
    ) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# fig 2b - individual excess risk by R0 and country
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs <- c('1.8'="R0 = 1.8",
              '2'="R0 = 2.0",
              '2.2'="R0 = 2.2")

fig2b <- ggplot(data = sia_dat[
  name=="base" | name=="screening" | name=="r0" | name=="r0_screening",
]) +
  geom_line(aes(
    x=t,
    y=ind_ex_risk, 
    color=factor(ppe_effect),
    linetype=scrn_sympt
  )) +
  facet_grid(
    country ~ r0,
    labeller = labeller(
      r0 = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option="viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  scale_linetype(
    name = "Symptomatic \nscreening",
    labels = line.labs
  ) +
  ggtitle("B - Excess risk to vaccinees/caregivers") +
  ylab("Excess infection risk to vaccinees/caregivers") + 
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

fig2 <- grid.arrange(fig2a, fig2b, ncol=1, nrow = 2)

ggsave("./figures/fig2 - hcw & ind ex risk.png",fig2,width=7.5, height=10, units="in")

# --------------------------------------------------------------------------------------------------
# fig S4 - excess risk assuming faster epidemic in hcws
# --------------------------------------------------------------------------------------------------

# figS4 a - hcw excess risk assuming faster epidemic in hcws
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs  <- c('soc_dist_0.6'="Reduction in vaccinator \nbaseline contacts: 40%",
               'soc_dist_0.7'="30%",
               'soc_dist_0.8'="20%"
)

figS4a <- ggplot(data = sia_dat[
  name=="base" |name=="v_hi_contact",
]) +
  geom_line(aes(
    x=t,
    y=vac_ex_risk, 
    color=factor(ppe_effect)
  )) +
  facet_grid(
    country ~ v_soc_dist,
    labeller = labeller(
      v_soc_dist = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option = "viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  ggtitle("A - Excess risk to vaccinators") +
  ylab("Excess infection risk to vaccinators") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# figS4 b - individual excess risk assuming faster epidemic in hcws
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs  <- c('soc_dist_0.6'="Reduction in vaccinator \nbaseline contacts: 40%",
               'soc_dist_0.7'="30%",
               'soc_dist_0.8'="20%"
)

figS4b <- ggplot(data = sia_dat[
  name=="base" |name=="v_hi_contact",
]) +
  geom_line(aes(
    x=t,
    y=ind_ex_risk, 
    color=factor(ppe_effect)
  )) +
  facet_grid(
    country ~ v_soc_dist,
    labeller = labeller(
      v_soc_dist = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option = "viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  ggtitle("B - Excess risk to vaccinees/caregivers") +
  ylab("Excess infection risk to vaccinees/caregivers") + 
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

figS4 <- grid.arrange(figS4a, figS4b, ncol=1, nrow = 2)
ggsave("./figures/figS4 - hcw & ind ex risk with faster hcw epidemic.png",figS4,width=7.5, height=10, units="in")

# --------------------------------------------------------------------------------------------------
# fig S5 - risk - house to house delivery
# --------------------------------------------------------------------------------------------------

# figS5 a - hcw risk - house to house delivery
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs <- c(
  '0'="Extra vaccinator community \ncontacts: None",
  '1'="One per household",
  '2'="Two per household"
)

figS5a <- ggplot(data = sia_dat[
  name %in% c("house-to-house","n_vj") & r0==2 & scrn_sympt==FALSE
]) +
  geom_line(aes(
    x=t,
    y=vac_ex_risk, 
    color=factor(ppe_effect)
  )) +
  facet_grid(
    country ~ n_vj,
    labeller = labeller(
      n_vj = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option="viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  ggtitle("A - Excess risk to vaccinators") +
  ylab("Excess infection risk to vaccinators") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# figS5 b - individual risk - house to house delivery
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs <- c(
  '0'="Extra vaccinator community \ncontacts: None",
  '1'="One per household",
  '2'="Two per household"
)


figS5b <- ggplot(data = sia_dat[
  name %in% c("house-to-house","n_vj") & r0==2 & scrn_sympt==FALSE
]) +
  geom_line(aes(
    x=t,
    y=ind_ex_risk, 
    color=factor(ppe_effect)
  )) +
  facet_grid(
    country ~ n_vj,
    labeller = labeller(
      n_vj = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option="viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  ggtitle("B - Excess risk to vaccinees/caregivers") +
  ylab("Excess infection risk to vaccinees/caregivers") + 
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL


figS5 <- grid.arrange(figS5a, figS5b, ncol=1, nrow = 2)
ggsave("./figures/figS5 - hcw & ind ex risk h-2-h delivery.png",figS5,width=7.5, height=10, units="in")

# --------------------------------------------------------------------------------------------------
# fig S6 - risk - contact reduction
# --------------------------------------------------------------------------------------------------

# figS6 a - hcw risk - under different levels of contact reduction
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs <- c(
  'soc_dist_0.8'="Reduction in baseline\ncontacts: 20%",
  'soc_dist_0.6'="40%",
  'soc_dist_0.4'="60%"
)

figS6a <- ggplot(data = sia_dat[
  name %in% c("base","soc_dist") & r0==2 & scrn_sympt==FALSE
]) +
  geom_line(aes(
    x=t,
    y=vac_ex_risk, 
    color=factor(ppe_effect)
  )) +
  facet_grid(
    country ~ soc_dist,
    labeller = labeller(
      soc_dist = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option="viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  ggtitle("A - Excess risk to vaccinators") +
  ylab("Excess infection risk to vaccinators") +
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

# figS6 b - individual risk - house to house delivery
col.labs  <- c("0%","50%","75%","90%","100%")
line.labs <- c("No","Yes")
fac.labs <- c(
  'soc_dist_0.8'="Reduction in baseline\ncontacts: 20%",
  'soc_dist_0.6'="40%",
  'soc_dist_0.4'="60%"
)


figS6b <- ggplot(data = sia_dat[
  name %in% c("base","soc_dist") & r0==2 & scrn_sympt==FALSE
]) +
  geom_line(aes(
    x=t,
    y=ind_ex_risk, 
    color=factor(ppe_effect)
  )) +
  facet_grid(
    country ~ soc_dist,
    labeller = labeller(
      soc_dist = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_color_viridis_d(
    option="viridis",
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  ggtitle("B - Excess risk to vaccinees/caregivers") +
  ylab("Excess infection risk to vaccinees/caregivers") + 
  scale_y_continuous(labels = function(x) sprintf("%.2g%%",x*100)) +
  scale_x_continuous(
    expand = c(0,0),
    name = "Days since introduction",
    breaks = seq(0,700,100)
  ) +
  theme_minimal() + 
  theme(legend.position="right") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL


figS6 <- grid.arrange(figS6a, figS6b, ncol=1, nrow = 2)
ggsave("./figures/figS6 - hcw & ind ex risk by contact reduction.png",figS6,width=7.5, height=10, units="in")

# --------------------------------------------------------------------------------------------------
# fig 3 - vaccinator tornado plot
# --------------------------------------------------------------------------------------------------

scen.labs = rev(c(
  expression(R[0]),
  "Reduction in community contacts due to NPIs",
  # "Number of vaccinations by vaccinator per day",
  "Campaign duration (fixed daily rate per vaccinator)",
  "Campaign duration (fixed total vaccinations per vaccinator)",
  "House-to-house campaign",
  "Higher baseline transmission among vaccinators",
  "Including symptomatic screening",
  "Lower PPE protection for vaccinator or vaccinee",
  "No. of contacts during journey to vaccine clinic",
  "Number of contacts with healthcare workers during visit",
  "Infectious duration",
  "Latent duration")
)

pal = c("#440154FF","#21908CFF") # to match 5 level discrete viridis

fac.labs = c('0'="",'0.75'="")
col.labs = c("0%","75%")
width = 0.8 # width

fig3 <- ggplot(data = max_sia_dat[ppe_effect %in% c(0,0.75)]) +

  geom_rect(
    aes(xmin = as.numeric(name) - width/2,
        xmax = as.numeric(name) + width/2,
        ymin = pmin(base_vac_ex_risk,vac_ex_risk_min),
        ymax = pmax(base_vac_ex_risk,vac_ex_risk_min),
        fill = factor(ppe_effect)
    ),
    alpha=0.5
  ) +
  geom_rect(
    aes(xmin = as.numeric(name) - width/2,
        xmax = as.numeric(name) + width/2,
        ymin = pmin(base_vac_ex_risk,vac_ex_risk_max),
        ymax = pmax(base_vac_ex_risk,vac_ex_risk_max),
        fill = factor(ppe_effect)
    ),
    alpha=0.8
  ) +
  geom_hline(aes(yintercept = base_vac_ex_risk)) +
  expand_limits(y=0) +
  coord_flip() +
  facet_grid(
    ppe_effect ~ country,
    labeller = labeller(
      ppe_effect = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_fill_discrete(
    type = pal,
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  scale_x_continuous(
    name="Parameter / Scenario",
    breaks = c(1:length(scen.labs)),
    labels = scen.labs
  ) +
  scale_y_continuous(
    name="Excess infection risk to vaccinators",
    labels = function(x) sprintf("%.2g%%",x*100)
  ) +
  ggtitle("Excess risk to vaccinators") +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  NULL

ggsave("./figures/fig3 - vaccinators sensitivity.png",fig3, width=10, height=6, units="in")

# --------------------------------------------------------------------------------------------------
# fig 4 - vaccinee / caregiver tonrnado plot
# --------------------------------------------------------------------------------------------------

pal = c("#440154FF","#21908CFF")

fac.labs = c('0'="",'0.75'="")
col.labs = c("0%","75%")
width = 0.8 # width

fig4 <- ggplot(data = max_sia_dat[ppe_effect %in% c(0,0.75)]) +
  geom_rect(
    aes(xmin = as.numeric(name) - width/2,
        xmax = as.numeric(name) + width/2,
        ymin = pmin(base_ind_ex_risk,ind_ex_risk_min),
        ymax = pmax(base_ind_ex_risk,ind_ex_risk_min),
        fill = factor(ppe_effect)
    ),
    alpha=0.5
  ) +
  geom_rect(
    aes(xmin = as.numeric(name) - width/2,
        xmax = as.numeric(name) + width/2,
        ymin = pmin(base_ind_ex_risk,ind_ex_risk_max),
        ymax = pmax(base_ind_ex_risk,ind_ex_risk_max),
        fill = factor(ppe_effect)
    ),
    alpha=0.8
  ) +
  geom_hline(aes(yintercept = base_ind_ex_risk)) +
  expand_limits(y=0) +
  coord_flip() +
  facet_grid(
    ppe_effect ~ country,
    labeller = labeller(
      ppe_effect = fac.labs,
      country = country.lbl
    )
  ) + 
  scale_fill_discrete(
    type = pal,
    name="PPE \neffectiveness",
    labels = col.labs
  ) +
  scale_x_continuous(
    name="Parameter / Scenario",
    breaks = c(1:length(scen.labs)),
    labels = scen.labs
  ) +
  scale_y_continuous(
    name="Excess infection risk to vaccinees/caregivers",
    labels = function(x) sprintf("%.2g%%",x*100)
  ) +
  ggtitle("Excess risk to vaccinees/caregivers") +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  NULL

ggsave("./figures/fig4 - ind sensitivity.png",fig4,width=10, height=6, units="in")


# --------------------------------------------------------------------------------------------------
# Plot population demographics and contact matrices
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# FigS1 population demographics
# --------------------------------------------------------------------------------------------------
pop_bfa <-  data.table(
      age =  factor(m$`Burkina Faso-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names,m$`Burkina Faso-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names),
      pop = m$`Burkina Faso-R0_1.8-soc_dist_0.6`$params$pop[[1]]$size,
      country = rep("Burkina Faso",length(m$`Burkina Faso-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names))
    )
pop_eth <- data.table(
      age =  factor(m$`Ethiopia-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names,m$`Ethiopia-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names),
      pop = m$`Ethiopia-R0_1.8-soc_dist_0.6`$params$pop[[1]]$size,
      country = rep("Ethiopia",length(m$`Ethiopia-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names))
    )
pop_chl <-  data.table(
      age =  factor(m$`Brazil-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names,m$`Brazil-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names),
      pop = m$`Brazil-R0_1.8-soc_dist_0.6`$params$pop[[1]]$size,
      country = rep("Brazil",length(m$`Brazil-R0_1.8-soc_dist_0.6`$params$pop[[1]]$group_names))
    )

figS1a <- ggplot(data=pop_bfa) +
  geom_bar(aes(x=age,y=100*pop/sum(pop)),stat="identity") +
  ggtitle("Population distribution by age: Burkina Faso") +
  xlab("Age band") + ylab("% of population") +
  ylim(0,20) +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

figS1b <- ggplot(data=pop_eth) +
  geom_bar(aes(x=age,y=100*pop/sum(pop)),stat="identity") +
  ggtitle("Population distribution by age: Ethiopia") +
  xlab("Age band") + ylab("% of population") +
  ylim(0,20) +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL

figS1c <- ggplot(data=pop_chl) +
  geom_bar(aes(x=age,y=100*pop/sum(pop)),stat="identity") +
  ggtitle("Population distribution by age: Brazil") +
  xlab("Age band") + ylab("% of population") +
  ylim(0,20) +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(panel.border=element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(angle = 45)) +
  NULL


figS1 <- grid.arrange(figS1a, figS1b, figS1c, ncol=3, nrow = 1)
ggsave("./figures/figS1 - population by age.png",figS1, width=14, height=5, units="in")

# --------------------------------------------------------------------------------------------------
# FigS2 contact matrices
# --------------------------------------------------------------------------------------------------
ml <- list(
  'Burkina Faso' =  m$`Burkina Faso-R0_1.8-soc_dist_0.6`$params$pop[[1]]$matrices,
  'Ethiopia' =  m$`Ethiopia-R0_1.8-soc_dist_0.6`$params$pop[[1]]$matrices,
  'Brazil' =  m$`Brazil-R0_1.8-soc_dist_0.6`$params$pop[[1]]$matrices
)

mat_dat <- NULL

for (c in 1:length(ml)) {
  for (s in 1:length(ml[[c]])) {
  country = names(ml[c])
  mat_dat <- rbind(
    mat_dat,
    cbind(
      country,
      reshape2::melt(ml[c][[1]][s])
    )
  )
  }
}

mat_dat <- as.data.table(mat_dat)
mat_dat[value<0,value:=0]
mat_dat[,country:=factor(country,levels=c("Burkina Faso","Ethiopia","Brazil"))]

figS2 <- ggplot(mat_dat, aes(x = Var1, y = Var2)) + 
  geom_raster(aes(fill=cut(value,c(-1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,3.0,6.0,max(value))))) + 
  facet_grid(L1~country) +
  scale_fill_viridis_d(
  ) +
  labs(
    title = "Country-specific contact matrices",
    x = "Age of individual",
    y = "Age of contact",
    fill = "Average daily number \nof contacts"
  )+
  theme_minimal() + 
  theme(
    legend.position="bottom", 
    ) +
  theme(axis.text.x = element_text(angle = 90)) +
  NULL



ggsave("./figures/figS2 - contact matrices.png",figS2, width=8, height=11, units="in")

# End
# --------------------------------------------------------------------------------------------------



