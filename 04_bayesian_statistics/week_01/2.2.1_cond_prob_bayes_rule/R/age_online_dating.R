d <- read.csv("June 10-July 12, 2015 - Gaming, Jobs and Broadband/June 10-July 12, 2015 - Gaming, Jobs and Broadband - CSV.csv")

library(dplyr)

d <- d %>%
  mutate(age_cat = ifelse(age <= 29, "18-29",
                          ifelse(age >= 30 & age <= 49, "30-49",
                                 ifelse(age >= 50 & age <= 64, "50-64",
                                        "65+")))) %>%
  filter(date1a %in% c(1, 2))

addmargins(table(d$date1a, d$age_cat))

table(d$date1a)/sum(table(d$date1a))
