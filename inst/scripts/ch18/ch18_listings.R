
##Listing 18.1:Target contribution rate in a defined contribution plan
library(mqriskR)

c_target <- contribution_rate_target(
  x = 30, z = 65, Sx = 60000,
  RR_target = 0.60,
  i = 0.06, adue_z = 11,
  g = 0.04)

RR_check <- replacement_ratio_dc(
  x = 30, z = 65, Sx = 60000,
  c = c_target, i = 0.06,
  adue_z = 11, g = 0.04)

c(c_target = c_target,
  RR_check = RR_check)


##Listing 18.2: Defined contribution accumulation and replacement ratio
library(mqriskR)

AV_65 <- AVz_dc(
  x = 30, z = 65, Sx = 50000,
  c = 0.10, i = 0.05, g = 0.04
)

Income_65 <- Income_dc(
  AVz = AV_65,
  adue_z = 12
)

RR_65 <- replacement_ratio_dc(
  x = 30, z = 65, Sx = 50000,
  c = 0.10, i = 0.05,
  adue_z = 12, g = 0.04
)

c(
  AV_65 = AV_65,
  Income_65 = Income_65,
  RR_65 = RR_65
)






##Listing 18.3: Projected annual benefit under a final average salary plan
library(mqriskR)

PAB_65 <- PAB_fas(
  x = 35, z = 65, CASx = 60000,
  p = 2, fas_years = 3, g = 0.04
)

final_salary <- 60000 * (1.04)^29

RR_final_salary <- replacement_ratio_db(
  benefit = PAB_65,
  salary = final_salary
)

c(PAB_65 = PAB_65,
  final_salary = final_salary,
  RR_final_salary = RR_final_salary)



