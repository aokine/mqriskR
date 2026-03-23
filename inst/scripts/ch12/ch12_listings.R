## Listing 12.1: Joint-life survival probability
library(mqriskR)

x <- 50; y <- 60; t <- 10

px <- tpx(x = x,t = t,model = "uniform",omega = 100)

py <- tpx(x = y,t = t,model = "uniform",omega = 100)

joint_from_product <- px * py

joint_from_package <- tpxy(x = x,y = y,t = t,model = "uniform",omega = 100)

c(px = px,py = py,
  joint_from_product = joint_from_product,
  joint_from_package = joint_from_package)




## Listing 12.2: Joint-life vs last-survivor survival
library(mqriskR)

x <- 50; y <- 60; t <- 10

joint_survival <- tpxy(x, y, t, model = "uniform", omega = 100)
last_survival  <- tpxybar(x, y, t, model = "uniform", omega = 100)

c(joint = joint_survival,
  last_survivor = last_survival)



##Listing 12.3: Joint-life and last-survivor annuity values
library(mqriskR)

x <- 50; y <- 60; i <- 0.05

joint_annuity <- adotxy(x = x,y = y,i = i,model = "uniform",omega = 100)

last_annuity <- adotxybar(x = x,y = y,i = i,model = "uniform",omega = 100)

c(joint_annuity = joint_annuity,
  last_survivor_annuity = last_annuity)




##Listing 12.4: Reversionary annuity to $(y)$ after death of $(x)$
library(mqriskR)

x <- 50; y <- 60; i <- 0.05

reversionary_annuity <- ax_y(x = x,y = y,i = i,model = "uniform", omega = 100)

check_value <- ax(x = y,i = i,model = "uniform",omega = 100) - axy(x = x, y = y, i = i, model = "uniform", omega = 100)

c(reversionary_annuity = reversionary_annuity,
  check_value = check_value)


##Listing 12.5: Continuous contingent insurance values}, label={lst:ch12_contingent_insurance}]
library(mqriskR)

x <- 50; y <- 60; i <- 0.05

benefit_if_x_before_y <- Abarxy1(x = x, y = y, i = i, model = "uniform", omega = 100)

benefit_if_x_after_y <- Abarxy2(x = x, y = y, i = i, model = "uniform", omega = 100)

c(x_before_y = benefit_if_x_before_y,
  x_after_y = benefit_if_x_after_y)

