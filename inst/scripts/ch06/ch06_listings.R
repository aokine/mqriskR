## Listing 6.1: Converting between $S_0(x)$ and $l_x$
library(mqriskR)

S0_vals <- c(1.00000, 0.97408, 0.97259, 0.97160, 0.97082)

lx_vals <- S0_to_lx(S0_vals, radix = 100000)
lx_vals

lx_to_S0(lx_vals)


## Listing 6.2: Constructing a life table from one-year survival probabilities
library(mqriskR)

tbl <- life_table(
  x = 0:4,
  px = c(0.90, 0.80, 0.60, 0.30, 0.00),
  radix = 10000
)


tbl

lx(tbl, 2)        # l_2
dx(tbl, 2)        # d_2
ndx(tbl, 1, 3)    # {}_3 d_1
npx(tbl, 1, 3)    # {}_3 p_1
nqx(tbl, 2, 3)    # {}_3 q_2


## Listing 6.3: Curtate and complete expectations from a life table
library(mqriskR)

tbl <- life_table(
  x = 0:4,
  px = c(0.90, 0.80, 0.60, 0.30, 0.00),
  radix = 10000
)

ex_curtate_tab(tbl, 2)                      # e_2
ex_complete_tab(tbl, 2, assumption = "udd") # e^o_2 under UDD

ex_temp_curtate_tab(tbl, x = 1, n = 3)
ex_temp_complete_tab(tbl, x = 2, n = 1.5, assumption = "cf")




## Listing 6.4: Fractional-age calculations under UDD, constant force, and Balducci
library(mqriskR)

tbl <- life_table(
  x = 0:4,
  px = c(0.90, 0.80, 0.60, 0.30, 0.00),
  radix = 10000
)

tpx_tab(tbl, x = 2, t = 0.25, assumption = "udd")
tpx_tab(tbl, x = 2, t = 0.25, assumption = "cf")
tpx_tab(tbl, x = 2, t = 0.25, assumption = "balducci")

mux_tab(tbl, x = 2, t = 0.25, assumption = "udd")
mux_tab(tbl, x = 2, t = 0.25, assumption = "cf")
mux_tab(tbl, x = 2, t = 0.25, assumption = "balducci")




## Listing 6.5: Working with a select life table
library(mqriskR)

sel_tbl <- select_life_table(
  x_sel = c(20,20,20,
            21,21,21,
            22,22,22),
  duration = c(0,1,2,
               0,1,2,
               0,1,2),
  attained_age = c(20,21,22,
                   21,22,23,
                   22,23,24),
  lx = c(100000, 99000, 97000,
         99500, 97500, 95000,
         98000, 95500, 92000)
)

lx_select(sel_tbl, 20, 1)                  # l_[20]+1
npx_select(sel_tbl, x_sel = 20, t = 0, n = 2)
nqx_select(sel_tbl, x_sel = 20, t = 0, n = 2)
nmxq_select(sel_tbl, x_sel = 20, t = 0, n = 1, m = 1)

