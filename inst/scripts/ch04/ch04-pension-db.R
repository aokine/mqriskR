# Chapter 4 — Defined benefit formula illustration (Section 4.3.1)

# NYC Teachers example: 1 2/3% * FAS * years
db_pension <- function(final_avg_salary, years_service, accrual = 1 + 2/3) {
  rate <- accrual / 100  # percent to decimal
  rate * final_avg_salary * years_service
}

db_pension(final_avg_salary = 90000, years_service = 25)
