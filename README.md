# winter_covid_infection_study-parameters

This repository contains code associated with the paper "Epidemiological Parameters of SARS-CoV-2 in the UK during the 2023/2024 Winter: A Cohort Study" by Overton C. E. et al. (2024). This paper concerns estimating paramters from the Winter Coronavirus Infection Study (https://www.gov.uk/government/statistics/winter-coronavirus-covid-19-infection-study-estimates-of-epidemiological-characteristics-england-and-scotland-2023-to-2024)

Stan files are provided for estimating the sensitivity and positivity and the incubation period. The sensitivity and positivity models are combined within a single stan file, to allow integration with Bayesian models for estimating incidence and prevalence as used in the wider study. Three incubation period models for different probability distrubutions, namely gamma, lognormal, and Weibull. 

Full posterior distributions for the fitted models will be uploaded during July 2024. Wrapper functions for calling the stan files and example input data formats will be added in late July 2024.
