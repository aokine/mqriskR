# Package index

## Interest Theory

- [`interest_convert()`](https://aokine.github.io/mqriskR/reference/interest_convert.md)
  : Convert between compound-interest quantities
- [`discount()`](https://aokine.github.io/mqriskR/reference/discount.md)
  : Discount factor for compound interest
- [`pv_cashflows()`](https://aokine.github.io/mqriskR/reference/pv_cashflows.md)
  : Present value of cash flows at time 0
- [`solve_yield()`](https://aokine.github.io/mqriskR/reference/solve_yield.md)
  : Solve the yield rate by the equation of value
- [`annuity_certain()`](https://aokine.github.io/mqriskR/reference/annuity_certain.md)
  : Present value of a level annuity-certain

## Survival Models

- [`S0()`](https://aokine.github.io/mqriskR/reference/S0.md) : Survival
  function for age-at-failure T0
- [`F0()`](https://aokine.github.io/mqriskR/reference/dist0.md)
  [`f0()`](https://aokine.github.io/mqriskR/reference/dist0.md) :
  Distribution functions for age-at-failure T0
- [`hazard0()`](https://aokine.github.io/mqriskR/reference/hazard0.md) :
  Hazard / force for age-at-failure T0
- [`cumhaz0()`](https://aokine.github.io/mqriskR/reference/cumhaz0.md) :
  Cumulative hazard for age-at-failure T0
- [`tpx()`](https://aokine.github.io/mqriskR/reference/tpx.md) :
  Conditional survival probability for Tx
- [`tqx()`](https://aokine.github.io/mqriskR/reference/tqx.md) :
  Conditional failure probability for Tx
- [`fx()`](https://aokine.github.io/mqriskR/reference/fx.md) :
  Conditional density for Tx
- [`ex_complete()`](https://aokine.github.io/mqriskR/reference/ex_complete.md)
  : Complete expectation of life at age x
- [`ex_curtate()`](https://aokine.github.io/mqriskR/reference/ex_curtate.md)
  : Curtate expectation of life at age x

## Life Tables

- [`S0_to_lx()`](https://aokine.github.io/mqriskR/reference/S0_to_lx.md)
  : Convert survival probabilities to life-table values
- [`lx_to_S0()`](https://aokine.github.io/mqriskR/reference/lx_to_S0.md)
  : Convert life-table values to survival probabilities
- [`px_to_lx()`](https://aokine.github.io/mqriskR/reference/px_to_lx.md)
  : Construct life-table values from p_x values
- [`qx_to_lx()`](https://aokine.github.io/mqriskR/reference/qx_to_lx.md)
  : Construct life-table values from q_x values
- [`life_table()`](https://aokine.github.io/mqriskR/reference/life_table.md)
  : Construct a life table
- [`lx()`](https://aokine.github.io/mqriskR/reference/lx.md) : Extract
  life-table survivor values
- [`dx()`](https://aokine.github.io/mqriskR/reference/dx.md) : Compute
  deaths between ages x and x+1
- [`ndx()`](https://aokine.github.io/mqriskR/reference/ndx.md) : Compute
  deaths over an n-year interval from a life table
- [`qx_tab()`](https://aokine.github.io/mqriskR/reference/qx_tab.md) :
  Compute one-year death probability from a life table
- [`npx()`](https://aokine.github.io/mqriskR/reference/npx.md) : Compute
  n-year survival probability from a life table
- [`nqx()`](https://aokine.github.io/mqriskR/reference/nqx.md) : Compute
  n-year death probability from a life table
- [`select_life_table()`](https://aokine.github.io/mqriskR/reference/select_life_table.md)
  : Construct a select life table
- [`lx_select()`](https://aokine.github.io/mqriskR/reference/lx_select.md)
  : Extract select-table survivor value
- [`npx_select()`](https://aokine.github.io/mqriskR/reference/npx_select.md)
  : Select-life survival probability
- [`nqx_select()`](https://aokine.github.io/mqriskR/reference/nqx_select.md)
  : Select-life death probability
- [`nmxq_select()`](https://aokine.github.io/mqriskR/reference/nmxq_select.md)
  : Deferred select-life death probability
- [`tpx_tab()`](https://aokine.github.io/mqriskR/reference/tpx_tab.md) :
  Fractional survival probability from a life table
- [`tqx_tab()`](https://aokine.github.io/mqriskR/reference/tqx_tab.md) :
  Fractional failure probability from a life table
- [`mux_tab()`](https://aokine.github.io/mqriskR/reference/mux_tab.md) :
  Fractional force of mortality from a life table
- [`fx_tab()`](https://aokine.github.io/mqriskR/reference/fx_tab.md) :
  Fractional conditional density from a life table
- [`nmxq()`](https://aokine.github.io/mqriskR/reference/nmxq.md) :
  Deferred death probability from a life table
- [`nkqx()`](https://aokine.github.io/mqriskR/reference/nkqx.md) :
  Curtate death probability from a life table
- [`ex_curtate_tab()`](https://aokine.github.io/mqriskR/reference/ex_curtate_tab.md)
  : Curtate expectation of life from a life table
- [`ex_temp_curtate_tab()`](https://aokine.github.io/mqriskR/reference/ex_temp_curtate_tab.md)
  : Temporary curtate expectation of life from a life table
- [`ex_temp_complete_tab()`](https://aokine.github.io/mqriskR/reference/ex_temp_complete_tab.md)
  : Temporary complete expectation of life from a life table
- [`ex_complete_tab()`](https://aokine.github.io/mqriskR/reference/ex_complete_tab.md)
  : Complete expectation of life from a life table

## Contingent Payment Models

- [`double_force_i()`](https://aokine.github.io/mqriskR/reference/double_force_i.md)
  : Effective annual interest at doubled force
- [`double_force_delta()`](https://aokine.github.io/mqriskR/reference/double_force_delta.md)
  : Doubled force of interest
- [`udd_continuous_multiplier()`](https://aokine.github.io/mqriskR/reference/udd_continuous_multiplier.md)
  : UDD multiplier for continuous insurance approximations
- [`udd_mthly_multiplier()`](https://aokine.github.io/mqriskR/reference/udd_mthly_multiplier.md)
  : UDD multiplier for m-thly insurance approximations
- [`Abarx_udd()`](https://aokine.github.io/mqriskR/reference/Abarx_udd.md)
  : UDD approximation of continuous whole life insurance
- [`Abarxn1_udd()`](https://aokine.github.io/mqriskR/reference/Abarxn1_udd.md)
  : UDD approximation of continuous term insurance
- [`nAbarx_udd()`](https://aokine.github.io/mqriskR/reference/nAbarx_udd.md)
  : UDD approximation of continuous deferred insurance
- [`Abarxn_udd()`](https://aokine.github.io/mqriskR/reference/Abarxn_udd.md)
  : UDD approximation of continuous endowment insurance
- [`Ax_m_udd()`](https://aokine.github.io/mqriskR/reference/Ax_m_udd.md)
  : UDD approximation of m-thly whole life insurance
- [`Axn1_m_udd()`](https://aokine.github.io/mqriskR/reference/Axn1_m_udd.md)
  : UDD approximation of m-thly term insurance
- [`nAx_m_udd()`](https://aokine.github.io/mqriskR/reference/nAx_m_udd.md)
  : UDD approximation of m-thly deferred insurance
- [`Axn_m_udd()`](https://aokine.github.io/mqriskR/reference/Axn_m_udd.md)
  : UDD approximation of m-thly endowment insurance
- [`Ax()`](https://aokine.github.io/mqriskR/reference/Ax.md) : Whole
  life insurance APV
- [`Axn1()`](https://aokine.github.io/mqriskR/reference/Axn1.md) : Term
  insurance APV
- [`nEx()`](https://aokine.github.io/mqriskR/reference/nEx.md) : Pure
  endowment APV
- [`nAx()`](https://aokine.github.io/mqriskR/reference/nAx.md) :
  Deferred insurance APV
- [`Axn()`](https://aokine.github.io/mqriskR/reference/Axn.md) :
  Endowment insurance APV
- [`A2x()`](https://aokine.github.io/mqriskR/reference/A2x.md) : Second
  moment of whole life insurance PV
- [`A2xn1()`](https://aokine.github.io/mqriskR/reference/A2xn1.md) :
  Second moment of term insurance PV
- [`A2nEx()`](https://aokine.github.io/mqriskR/reference/A2nEx.md) :
  Second moment of pure endowment PV
- [`A2nAx()`](https://aokine.github.io/mqriskR/reference/A2nAx.md) :
  Second moment of deferred insurance PV
- [`A2xn()`](https://aokine.github.io/mqriskR/reference/A2xn.md) :
  Second moment of endowment insurance PV
- [`var_Ax()`](https://aokine.github.io/mqriskR/reference/var_Ax.md) :
  Variance of whole life insurance PV
- [`var_Axn1()`](https://aokine.github.io/mqriskR/reference/var_Axn1.md)
  : Variance of term insurance PV
- [`var_nEx()`](https://aokine.github.io/mqriskR/reference/var_nEx.md) :
  Variance of pure endowment PV
- [`var_nAx()`](https://aokine.github.io/mqriskR/reference/var_nAx.md) :
  Variance of deferred insurance PV
- [`var_Axn()`](https://aokine.github.io/mqriskR/reference/var_Axn.md) :
  Variance of endowment insurance PV
- [`cov_term_deferred()`](https://aokine.github.io/mqriskR/reference/cov_term_deferred.md)
  : Covariance of term and deferred insurance PVs
- [`cov_term_endow()`](https://aokine.github.io/mqriskR/reference/cov_term_endow.md)
  : Covariance of term insurance and pure endowment PVs
- [`Abarx()`](https://aokine.github.io/mqriskR/reference/Abarx.md) :
  Continuous whole life insurance APV
- [`Abarxn1()`](https://aokine.github.io/mqriskR/reference/Abarxn1.md) :
  Continuous term insurance APV
- [`nAbarx()`](https://aokine.github.io/mqriskR/reference/nAbarx.md) :
  Continuous deferred insurance APV
- [`Abarxn()`](https://aokine.github.io/mqriskR/reference/Abarxn.md) :
  Continuous endowment insurance APV
- [`A2barx()`](https://aokine.github.io/mqriskR/reference/A2barx.md) :
  Second moment of continuous whole life insurance PV
- [`A2barxn1()`](https://aokine.github.io/mqriskR/reference/A2barxn1.md)
  : Second moment of continuous term insurance PV
- [`A2nAbarx()`](https://aokine.github.io/mqriskR/reference/A2nAbarx.md)
  : Second moment of continuous deferred insurance PV
- [`A2barxn()`](https://aokine.github.io/mqriskR/reference/A2barxn.md) :
  Second moment of continuous endowment insurance PV
- [`var_Abarx()`](https://aokine.github.io/mqriskR/reference/var_Abarx.md)
  : Variance of continuous whole life insurance PV
- [`var_Abarxn1()`](https://aokine.github.io/mqriskR/reference/var_Abarxn1.md)
  : Variance of continuous term insurance PV
- [`var_nAbarx()`](https://aokine.github.io/mqriskR/reference/var_nAbarx.md)
  : Variance of continuous deferred insurance PV
- [`var_Abarxn()`](https://aokine.github.io/mqriskR/reference/var_Abarxn.md)
  : Variance of continuous endowment insurance PV
- [`Ax_m()`](https://aokine.github.io/mqriskR/reference/Ax_m.md) :
  m-thly whole life insurance APV
- [`Axn1_m()`](https://aokine.github.io/mqriskR/reference/Axn1_m.md) :
  m-thly term insurance APV
- [`nAx_m()`](https://aokine.github.io/mqriskR/reference/nAx_m.md) :
  m-thly deferred insurance APV
- [`Axn_m()`](https://aokine.github.io/mqriskR/reference/Axn_m.md) :
  m-thly endowment insurance APV
- [`A2x_m()`](https://aokine.github.io/mqriskR/reference/A2x_m.md) :
  Second moment of m-thly whole life insurance PV
- [`A2xn1_m()`](https://aokine.github.io/mqriskR/reference/A2xn1_m.md) :
  Second moment of m-thly term insurance PV
- [`A2nAx_m()`](https://aokine.github.io/mqriskR/reference/A2nAx_m.md) :
  Second moment of m-thly deferred insurance PV
- [`A2xn_m()`](https://aokine.github.io/mqriskR/reference/A2xn_m.md) :
  Second moment of m-thly endowment insurance PV
- [`var_Ax_m()`](https://aokine.github.io/mqriskR/reference/var_Ax_m.md)
  : Variance of m-thly whole life insurance PV
- [`var_Axn1_m()`](https://aokine.github.io/mqriskR/reference/var_Axn1_m.md)
  : Variance of m-thly term insurance PV
- [`var_nAx_m()`](https://aokine.github.io/mqriskR/reference/var_nAx_m.md)
  : Variance of m-thly deferred insurance PV
- [`var_Axn_m()`](https://aokine.github.io/mqriskR/reference/var_Axn_m.md)
  : Variance of m-thly endowment insurance PV
- [`IAx()`](https://aokine.github.io/mqriskR/reference/IAx.md) :
  Increasing whole life insurance
- [`IAxn1()`](https://aokine.github.io/mqriskR/reference/IAxn1.md) :
  Increasing n-year term insurance
- [`DAxn1()`](https://aokine.github.io/mqriskR/reference/DAxn1.md) :
  Decreasing n-year term insurance
- [`IbarAbarx()`](https://aokine.github.io/mqriskR/reference/IbarAbarx.md)
  : Fully continuous increasing whole life insurance
- [`IAbarx()`](https://aokine.github.io/mqriskR/reference/IAbarx.md) :
  Piecewise-continuous increasing whole life insurance
- [`IbarAbarxn1()`](https://aokine.github.io/mqriskR/reference/IbarAbarxn1.md)
  : Fully continuous increasing n-year term insurance
- [`DbarAbarxn1()`](https://aokine.github.io/mqriskR/reference/DbarAbarxn1.md)
  : Fully continuous decreasing n-year term insurance
- [`DAbarxn1()`](https://aokine.github.io/mqriskR/reference/DAbarxn1.md)
  : Piecewise-continuous decreasing n-year term insurance

## Contingent Annuity Models

- [`ax()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`adotx()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`abarx()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`axn()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`adotxn()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`abarxn()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`nax()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`nadotx()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`nabarx()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`sxn()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`sdotxn()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  [`sbarxn()`](https://aokine.github.io/mqriskR/reference/annuity_annual.md)
  : Annual annuity functions (Chapter 8)
- [`ax_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_whole_immediate.md)
  : Whole life m-thly annuity-immediate
- [`adotx_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_whole_due.md)
  : Whole life m-thly annuity-due
- [`axn_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_temp_immediate.md)
  : Temporary m-thly annuity-immediate
- [`adotxn_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_temp_due.md)
  : Temporary m-thly annuity-due
- [`nax_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_deferred_immediate.md)
  : Deferred whole life m-thly annuity-immediate
- [`nadotx_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_deferred_due.md)
  : Deferred whole life m-thly annuity-due
- [`sxn_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_accum_immediate.md)
  : Temporary m-thly annuity-immediate actuarial accumulated value
- [`sdotxn_m()`](https://aokine.github.io/mqriskR/reference/annuity_mthly_accum_due.md)
  : Temporary m-thly annuity-due actuarial accumulated value
- [`annuity_identity_ax()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_adotx()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_abarx()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_axn()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_adotxn()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_abarxn()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_nax()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_nadotx()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  [`annuity_identity_nabarx()`](https://aokine.github.io/mqriskR/reference/annuity_relationships.md)
  : Annuity-insurance relationships (Chapter 8)
- [`Iax()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Iaxn()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Daxn()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Iadotx()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Iadotxn()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Dadotxn()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Iabarx()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Iabarxn()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  [`Dabarxn()`](https://aokine.github.io/mqriskR/reference/annuity_varying_payments.md)
  : Varying-payment annuity functions (Chapter 8)
- [`qx_proj()`](https://aokine.github.io/mqriskR/reference/qx_proj.md) :
  Project one-year death probability under mortality improvement
- [`px_proj()`](https://aokine.github.io/mqriskR/reference/px_proj.md) :
  Project one-year survival probability under mortality improvement
- [`tpx_improved()`](https://aokine.github.io/mqriskR/reference/tpx_improved.md)
  : Multi-year survival probability under mortality improvement
- [`axn_improved()`](https://aokine.github.io/mqriskR/reference/axn_improved.md)
  : Temporary annuity-immediate under mortality improvement
- [`naxn_improved()`](https://aokine.github.io/mqriskR/reference/naxn_improved.md)
  : Deferred temporary annuity-immediate under mortality improvement
- [`ax_improved()`](https://aokine.github.io/mqriskR/reference/ax_improved.md)
  : Whole life annuity-immediate under mortality improvement
- [`adotx_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`adotxn_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`nadotx_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`ax_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`axn_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`nax_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`sdotxn_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`sxn_m_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`abarx_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`abarxn_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  [`nabarx_udd()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_udd.md)
  : UDD annuity approximations
- [`ax_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`adotx_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`nax_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`nadotx_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`axn_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`adotxn_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`sxn_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`sdotxn_m_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  [`abarx_woolhouse2()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse2.md)
  : Woolhouse 2-term annuity approximations
- [`ax_m_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  [`adotx_m_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  [`nax_m_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  [`nadotx_m_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  [`axn_m_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  [`adotxn_m_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  [`abarx_woolhouse3()`](https://aokine.github.io/mqriskR/reference/annuity_approximations_woolhouse3.md)
  : Woolhouse 3-term annuity approximations

## Funding Plans for Contingent Contracts

- [`Px()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pxn1()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`PnEx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pxn()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`tPx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`tPxn1()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`tPnEx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`tPxn()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`PnAx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`tPnAx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pbarx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pbarxn1()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pbarxn()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`PbarAbarx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`PbarAbarxn1()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`PbarAbarxn()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Px_m()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pxn1_m()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Pxn_m()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`PnAx_m()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`EL0x()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`varL0x()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`EL0xn1()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`varL0xn1()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`EL0xn()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`varL0xn()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`EL0barAbarx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`varL0barAbarx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md)
  [`Gx()`](https://aokine.github.io/mqriskR/reference/premium_ch9.md) :
  Chapter 9 premium, loss, and expense functions

## Contingent Contract Reserves

- [`tVx()`](https://aokine.github.io/mqriskR/reference/tVx.md) : Whole
  life net level premium reserve
- [`tVxn1()`](https://aokine.github.io/mqriskR/reference/tVxn1.md) :
  Term insurance net level premium reserve
- [`tVnEx()`](https://aokine.github.io/mqriskR/reference/tVnEx.md) :
  Pure endowment net level premium reserve
- [`tVxn()`](https://aokine.github.io/mqriskR/reference/tVxn.md) :
  Endowment insurance net level premium reserve
- [`htVx()`](https://aokine.github.io/mqriskR/reference/htVx.md) : h-pay
  whole life net level premium reserve
- [`ELtx()`](https://aokine.github.io/mqriskR/reference/ELtx.md) : Mean
  present value of loss at duration t for whole life insurance
- [`varLtx()`](https://aokine.github.io/mqriskR/reference/varLtx.md) :
  Variance of present value of loss at duration t for whole life
  insurance
- [`tVbarx()`](https://aokine.github.io/mqriskR/reference/tVbarx.md) :
  Whole life reserve with continuous premiums
- [`tVbarAbarx()`](https://aokine.github.io/mqriskR/reference/tVbarAbarx.md)
  : Fully continuous whole life reserve
- [`tVx_m()`](https://aokine.github.io/mqriskR/reference/tVx_m.md) :
  Whole life reserve with m-thly premiums
- [`GT_disc()`](https://aokine.github.io/mqriskR/reference/GT_disc.md) :
  Total gain for a discrete insurance contract
- [`GM_disc()`](https://aokine.github.io/mqriskR/reference/GM_disc.md) :
  Mortality gain for a discrete insurance contract
- [`GI_disc()`](https://aokine.github.io/mqriskR/reference/GI_disc.md) :
  Interest gain for a discrete insurance contract
- [`GT_cont()`](https://aokine.github.io/mqriskR/reference/GT_cont.md) :
  Total gain for a continuous-style one-step recursion
- [`GM_cont()`](https://aokine.github.io/mqriskR/reference/GM_cont.md) :
  Mortality gain helper for continuous-style recursion
- [`GI_cont()`](https://aokine.github.io/mqriskR/reference/GI_cont.md) :
  Interest gain helper for continuous-style recursion
- [`tVx_ret()`](https://aokine.github.io/mqriskR/reference/tVx_ret.md) :
  Whole life net level premium reserve by retrospective method
- [`tVxn_ret()`](https://aokine.github.io/mqriskR/reference/tVxn_ret.md)
  : Endowment insurance reserve by retrospective method
- [`tVxn1_ret()`](https://aokine.github.io/mqriskR/reference/tVxn1_ret.md)
  : Term insurance reserve by retrospective method
- [`tVnAx()`](https://aokine.github.io/mqriskR/reference/reserve_deferred_insurance_ch10.md)
  [`htVnAx()`](https://aokine.github.io/mqriskR/reference/reserve_deferred_insurance_ch10.md)
  : Deferred insurance reserve functions
- [`PnAdotx()`](https://aokine.github.io/mqriskR/reference/PnAdotx.md) :
  Deferred annuity-due premium
- [`tVnAdotx()`](https://aokine.github.io/mqriskR/reference/tVnAdotx.md)
  : Deferred annuity-due reserve
- [`Pnax()`](https://aokine.github.io/mqriskR/reference/Pnax.md) :
  Deferred annuity-immediate premium
- [`tVnax()`](https://aokine.github.io/mqriskR/reference/tVnax.md) :
  Deferred annuity-immediate reserve
- [`thiele_backward_step()`](https://aokine.github.io/mqriskR/reference/thiele_backward_step.md)
  : One backward Euler-style Thiele step
- [`thiele_dVdt()`](https://aokine.github.io/mqriskR/reference/thiele_dVdt.md)
  : Reserve derivative from Thiele's equation
- [`thiele_backward_path()`](https://aokine.github.io/mqriskR/reference/thiele_backward_path.md)
  : Backward Euler reserve path from maturity

## Reserves as Financial Liabilities

- [`alphaF()`](https://aokine.github.io/mqriskR/reference/alphaF.md) :
  Full preliminary term first-year modified premium
- [`betaF()`](https://aokine.github.io/mqriskR/reference/betaF.md) :
  Full preliminary term renewal modified premium
- [`tVFx()`](https://aokine.github.io/mqriskR/reference/tVFx.md) : Full
  preliminary term reserve for whole life insurance
- [`tsVx()`](https://aokine.github.io/mqriskR/reference/tsVx.md) :
  Fractional-duration whole life reserve
- [`meanVx()`](https://aokine.github.io/mqriskR/reference/meanVx.md) :
  Mean reserve for whole life insurance
- [`tsVxn()`](https://aokine.github.io/mqriskR/reference/tsVxn.md) :
  Fractional-duration endowment reserve
- [`tsVxn1()`](https://aokine.github.io/mqriskR/reference/tsVxn1.md) :
  Fractional-duration term reserve
- [`tVGx()`](https://aokine.github.io/mqriskR/reference/tVGx.md) : Whole
  life gross premium reserve
- [`tVEx()`](https://aokine.github.io/mqriskR/reference/tVEx.md) : Whole
  life expense reserve
- [`GTg_disc()`](https://aokine.github.io/mqriskR/reference/GTg_disc.md)
  : Total gross gain for a discrete insurance contract
- [`decompGg_disc()`](https://aokine.github.io/mqriskR/reference/decompGg_disc.md)
  : Ordered gross gain decomposition

## Multi-Life Models

- [`tpxy()`](https://aokine.github.io/mqriskR/reference/tpxy.md) :
  Joint-life survival probability
- [`tqxy()`](https://aokine.github.io/mqriskR/reference/tqxy.md) :
  Joint-life failure probability
- [`tpxybar()`](https://aokine.github.io/mqriskR/reference/tpxybar.md) :
  Last-survivor survival probability
- [`tqxybar()`](https://aokine.github.io/mqriskR/reference/tqxybar.md) :
  Last-survivor failure probability
- [`tqxy1()`](https://aokine.github.io/mqriskR/reference/tqxy1.md) :
  Probability that (x) fails before (y) within n years
- [`tqyx1()`](https://aokine.github.io/mqriskR/reference/tqyx1.md) :
  Probability that (y) fails before (x) within n years
- [`tqxy2()`](https://aokine.github.io/mqriskR/reference/tqxy2.md) :
  Probability that (x) fails after (y) within n years
- [`tqyx2()`](https://aokine.github.io/mqriskR/reference/tqyx2.md) :
  Probability that (y) fails after (x) within n years
- [`nExy()`](https://aokine.github.io/mqriskR/reference/nExy.md) :
  Joint-life pure endowment
- [`nExybar()`](https://aokine.github.io/mqriskR/reference/nExybar.md) :
  Last-survivor pure endowment
- [`adotxyn()`](https://aokine.github.io/mqriskR/reference/adotxyn.md) :
  Joint-life temporary annuity-due
- [`axyn()`](https://aokine.github.io/mqriskR/reference/axyn_ch12.md) :
  Joint-life temporary annuity-immediate
- [`adotxy()`](https://aokine.github.io/mqriskR/reference/adotxy.md) :
  Joint-life whole life annuity-due
- [`axy()`](https://aokine.github.io/mqriskR/reference/axy_ch12.md) :
  Joint-life whole life annuity-immediate
- [`Axyn1()`](https://aokine.github.io/mqriskR/reference/Axyn1.md) :
  Joint-life term insurance
- [`Axyn()`](https://aokine.github.io/mqriskR/reference/Axyn.md) :
  Joint-life endowment insurance
- [`Axy()`](https://aokine.github.io/mqriskR/reference/Axy.md) :
  Joint-life whole life insurance
- [`adotxybarn()`](https://aokine.github.io/mqriskR/reference/adotxybarn.md)
  : Last-survivor temporary annuity-due
- [`axybarn()`](https://aokine.github.io/mqriskR/reference/axybarn_ch12.md)
  : Last-survivor temporary annuity-immediate
- [`adotxybar()`](https://aokine.github.io/mqriskR/reference/adotxybar.md)
  : Last-survivor whole life annuity-due
- [`axybar()`](https://aokine.github.io/mqriskR/reference/axybar_ch12.md)
  : Last-survivor whole life annuity-immediate
- [`Axybarn1()`](https://aokine.github.io/mqriskR/reference/Axybarn1.md)
  : Last-survivor term insurance
- [`Axybarn()`](https://aokine.github.io/mqriskR/reference/Axybarn.md) :
  Last-survivor endowment insurance
- [`Axybar()`](https://aokine.github.io/mqriskR/reference/Axybar.md) :
  Last-survivor whole life insurance
- [`ax_y()`](https://aokine.github.io/mqriskR/reference/ax_y.md) :
  Reversionary annuity to (y) after death of (x)
- [`ay_x()`](https://aokine.github.io/mqriskR/reference/ay_x.md) :
  Reversionary annuity to (x) after death of (y)
- [`abarx_y()`](https://aokine.github.io/mqriskR/reference/abarx_y.md) :
  Continuous reversionary annuity to (y) after death of (x)
- [`abary_x()`](https://aokine.github.io/mqriskR/reference/abary_x.md) :
  Continuous reversionary annuity to (x) after death of (y)
- [`abarxy()`](https://aokine.github.io/mqriskR/reference/abarxy_ch12.md)
  : Continuous joint-life whole life annuity
- [`abarxybar()`](https://aokine.github.io/mqriskR/reference/abarxybar_ch12.md)
  : Continuous last-survivor whole life annuity
- [`Abarxy()`](https://aokine.github.io/mqriskR/reference/Abarxy.md) :
  Continuous joint-life whole life insurance
- [`Abarxybar()`](https://aokine.github.io/mqriskR/reference/Abarxybar.md)
  : Continuous last-survivor whole life insurance
- [`Abarxy1()`](https://aokine.github.io/mqriskR/reference/Abarxy1.md) :
  Continuous contingent insurance: benefit on death of (x) if before (y)
- [`Abaryx1()`](https://aokine.github.io/mqriskR/reference/Abaryx1.md) :
  Continuous contingent insurance: benefit on death of (y) if before (x)
- [`Abarxy2()`](https://aokine.github.io/mqriskR/reference/Abarxy2.md) :
  Continuous contingent insurance: benefit on death of (x) if after (y)
- [`Abaryx2()`](https://aokine.github.io/mqriskR/reference/Abaryx2.md) :
  Continuous contingent insurance: benefit on death of (y) if after (x)

## Multiple-Decrement Models

- [`qxtau()`](https://aokine.github.io/mqriskR/reference/qxtau.md) :
  Total probability of decrement \\q_x^{(\tau)}\\
- [`pxtau()`](https://aokine.github.io/mqriskR/reference/pxtau.md) :
  Survival probability \\p_x^{(\tau)}\\
- [`dxj()`](https://aokine.github.io/mqriskR/reference/dxj.md) :
  Cause-specific decrements \\d_x^{(j)}\\
- [`dxtau()`](https://aokine.github.io/mqriskR/reference/dxtau.md) :
  Total decrements \\d_x^{(\tau)}\\
- [`md_table()`](https://aokine.github.io/mqriskR/reference/md_table.md)
  : Build a multiple-decrement table
- [`npxtau_md()`](https://aokine.github.io/mqriskR/reference/npxtau_md.md)
  : \\{}\_n p_x^{(\tau)}\\ from a multiple-decrement table
- [`nqxj_md()`](https://aokine.github.io/mqriskR/reference/nqxj_md.md) :
  \\{}\_n q_x^{(j)}\\ from a multiple-decrement table
- [`nqxtau_md()`](https://aokine.github.io/mqriskR/reference/nqxtau_md.md)
  : \\{}\_n q_x^{(\tau)}\\ from a multiple-decrement table
- [`tpxprimej_cf()`](https://aokine.github.io/mqriskR/reference/tpxprimej_cf.md)
  : Single-decrement survival probability \\{}\_t p_x^{\prime(j)}\\
  under constant force
- [`tqxprimej_cf()`](https://aokine.github.io/mqriskR/reference/tqxprimej_cf.md)
  : Single-decrement failure probability \\{}\_t q_x^{\prime(j)}\\ under
  constant force
- [`tpx_tau_cf()`](https://aokine.github.io/mqriskR/reference/tpx_tau_cf.md)
  : Total survival probability \\{}\_t p_x^{(\tau)}\\ under constant
  forces
- [`tqxj_cf()`](https://aokine.github.io/mqriskR/reference/tqxj_cf.md) :
  Cause-specific probability \\{}\_t q_x^{(j)}\\ under constant forces
- [`qx_dep_cf()`](https://aokine.github.io/mqriskR/reference/qx_dep_cf.md)
  : Dependent probabilities \\q_x^{(j)}\\ from independent probabilities
  \\q_x^{\prime(j)}\\ under constant force
- [`qxprime_mudd()`](https://aokine.github.io/mqriskR/reference/qxprime_mudd.md)
  : Independent probabilities \\q_x^{\prime(j)}\\ from dependent
  probabilities \\q_x^{(j)}\\ under MUDD
- [`tqxprime_mudd()`](https://aokine.github.io/mqriskR/reference/tqxprime_mudd.md)
  : Independent probabilities \\{}\_t q_x^{\prime(j)}\\ under MUDD
- [`qx_dep_sudd()`](https://aokine.github.io/mqriskR/reference/qx_dep_sudd.md)
  : Dependent probabilities \\q_x^{(j)}\\ from independent probabilities
  \\q_x^{\prime(j)}\\ under SUDD
- [`qxprime_sudd()`](https://aokine.github.io/mqriskR/reference/qxprime_sudd.md)
  : Independent probabilities \\q_x^{\prime(j)}\\ from dependent
  probabilities \\q_x^{(j)}\\ under SUDD

## Multiple-Decrement Applications

- [`Axj_md()`](https://aokine.github.io/mqriskR/reference/Axj_md.md) :
  Discrete multiple-decrement insurance APV \\A\_{x}^{(j)}\\
- [`Abarxj_md()`](https://aokine.github.io/mqriskR/reference/Abarxj_md.md)
  : Continuous multiple-decrement insurance APV
  \\\overline{A}\_{x}^{(j)}\\
- [`AS_path()`](https://aokine.github.io/mqriskR/reference/AS_path.md) :
  Projected asset share path \\{}\_{k}AS\\
- [`AS_path_md()`](https://aokine.github.io/mqriskR/reference/AS_path_md.md)
  : General projected asset share path (multiple decrements)
- [`tp00_tp01_euler()`](https://aokine.github.io/mqriskR/reference/tp00_tp01_euler.md)
  : Euler approximation for \\{}\_{t}p\_{x}^{00}\\ and
  \\{}\_{t}p\_{x}^{01}\\
- [`Pbar_trapz_ms()`](https://aokine.github.io/mqriskR/reference/Pbar_trapz_ms.md)
  : Continuous premium approximation \\\overline{P}\\ by trapezoidal
  rule
- [`thiele_dVdt_01()`](https://aokine.github.io/mqriskR/reference/thiele_dVdt_01.md)
  : Reserve derivatives for the disability model with recovery
- [`thiele_path_01()`](https://aokine.github.io/mqriskR/reference/thiele_path_01.md)
  : Backward reserve path for the disability model with recovery
- [`markov_nstep_prob()`](https://aokine.github.io/mqriskR/reference/markov_nstep_prob.md)
  : n-step transition probability for a discrete-time Markov chain
- [`gain_loss_md()`](https://aokine.github.io/mqriskR/reference/gain_loss_md.md)
  : Gain or loss in a multiple-decrement model

## Variable Interest Models

- [`nEx_var()`](https://aokine.github.io/mqriskR/reference/chapter15_variable_interest_apv.md)
  [`Axn1_var()`](https://aokine.github.io/mqriskR/reference/chapter15_variable_interest_apv.md)
  [`Axn_var()`](https://aokine.github.io/mqriskR/reference/chapter15_variable_interest_apv.md)
  [`axn_var()`](https://aokine.github.io/mqriskR/reference/chapter15_variable_interest_apv.md)
  : Variable-interest actuarial present value functions
- [`nEx_spot()`](https://aokine.github.io/mqriskR/reference/chapter15_spot_interest_apv.md)
  [`Axn1_spot()`](https://aokine.github.io/mqriskR/reference/chapter15_spot_interest_apv.md)
  [`Axn_spot()`](https://aokine.github.io/mqriskR/reference/chapter15_spot_interest_apv.md)
  [`axn_spot()`](https://aokine.github.io/mqriskR/reference/chapter15_spot_interest_apv.md)
  : Spot-rate actuarial present value functions
- [`vt_var()`](https://aokine.github.io/mqriskR/reference/vt_var.md) :
  Discount factors under a variable annual interest scenario
- [`pv_spot_cashflows()`](https://aokine.github.io/mqriskR/reference/pv_spot_cashflows.md)
  : Present value of cash flows using spot rates
- [`z_from_coupon_semi()`](https://aokine.github.io/mqriskR/reference/z_from_coupon_semi.md)
  : Bootstrap semiannual nominal spot rates from coupon-bond yields
- [`z_from_coupon_annual()`](https://aokine.github.io/mqriskR/reference/z_from_coupon_annual.md)
  : Bootstrap annual spot rates from annual coupon-bond yields
- [`fnk_from_z()`](https://aokine.github.io/mqriskR/reference/fnk_from_z.md)
  : Forward rate \\f\_{n,k}\\ from spot rates
- [`forward_matrix_from_z()`](https://aokine.github.io/mqriskR/reference/forward_matrix_from_z.md)
  : Matrix of all determinable forward rates from spot rates
- [`z_from_fn1()`](https://aokine.github.io/mqriskR/reference/z_from_fn1.md)
  : Spot rates from forward one-year rates

## Universal Life Insurance

- [`coi_ul_typeB()`](https://aokine.github.io/mqriskR/reference/coi_ul_typeB.md)
  : Cost of insurance for Type B universal life
- [`AV_path_ul_typeB()`](https://aokine.github.io/mqriskR/reference/AV_path_ul_typeB.md)
  : Account-value path for Type B universal life
- [`AV_path_ul_typeA()`](https://aokine.github.io/mqriskR/reference/AV_path_ul_typeA.md)
  : Account-value path for Type A universal life
- [`iP_eiul()`](https://aokine.github.io/mqriskR/reference/iP_eiul.md) :
  Point-to-point index growth rates
- [`iMA_eiul()`](https://aokine.github.io/mqriskR/reference/iMA_eiul.md)
  : Monthly-average index growth rate
- [`i_credit_eiul()`](https://aokine.github.io/mqriskR/reference/i_credit_eiul.md)
  : Credited rates from raw index growth rates
- [`pxtau_ul()`](https://aokine.github.io/mqriskR/reference/pxtau_ul.md)
  : One-year persistency rates for universal life
- [`tpxtau_ul()`](https://aokine.github.io/mqriskR/reference/tpxtau_ul.md)
  : Cumulative persistency to the end of each policy year
- [`GMF_rollforward_ul()`](https://aokine.github.io/mqriskR/reference/GMF_rollforward_ul.md)
  : Guaranteed maturity fund roll-forward
- [`rt_ul()`](https://aokine.github.io/mqriskR/reference/rt_ul.md) :
  Ratio \\r_t = AV_t / GMF_t\\ capped at 1
- [`Vprefloor_crvm_ul()`](https://aokine.github.io/mqriskR/reference/Vprefloor_crvm_ul.md)
  : Pre-floor CRVM reserve for universal life
- [`ag38_prefunding_ratio()`](https://aokine.github.io/mqriskR/reference/ag38_prefunding_ratio.md)
  : AG 38 prefunding ratio
- [`ag38_reserve_ul()`](https://aokine.github.io/mqriskR/reference/ag38_reserve_ul.md)
  : AG 38 reserve calculation

## Profit Analysis

- [`Pr_vector_disc()`](https://aokine.github.io/mqriskR/reference/Pr_vector_disc.md)
  : Profit vector for a discrete profit-analysis model
- [`Pi_signature()`](https://aokine.github.io/mqriskR/reference/Pi_signature.md)
  : Profit signature from a profit vector
- [`NPV_profit()`](https://aokine.github.io/mqriskR/reference/NPV_profit.md)
  : Net present value of a profit signature
- [`NPV_partial()`](https://aokine.github.io/mqriskR/reference/NPV_partial.md)
  : Partial net present values
- [`discounted_payback_period()`](https://aokine.github.io/mqriskR/reference/discounted_payback_period.md)
  : Discounted payback period
- [`IRR_profit()`](https://aokine.github.io/mqriskR/reference/IRR_profit.md)
  : Internal rate of return of a profit signature
- [`APV_gross_premiums()`](https://aokine.github.io/mqriskR/reference/APV_gross_premiums.md)
  : APV of gross premiums under a risk discount rate
- [`profit_margin()`](https://aokine.github.io/mqriskR/reference/profit_margin.md)
  : Profit margin
- [`V_zeroized()`](https://aokine.github.io/mqriskR/reference/V_zeroized.md)
  : Zeroized reserves for a discrete death-only contract

## Pension Benefits

- [`salary_scale()`](https://aokine.github.io/mqriskR/reference/salary_scale.md)
  : Salary scale under constant annual growth
- [`AVz_dc()`](https://aokine.github.io/mqriskR/reference/AVz_dc.md) :
  Accumulated value of defined contribution plan contributions
- [`Income_dc()`](https://aokine.github.io/mqriskR/reference/Income_dc.md)
  : Retirement income from a defined contribution accumulation
- [`replacement_ratio_dc()`](https://aokine.github.io/mqriskR/reference/replacement_ratio_dc.md)
  : Replacement ratio for a defined contribution plan
- [`contribution_rate_target()`](https://aokine.github.io/mqriskR/reference/contribution_rate_target.md)
  : Target contribution rate for a defined contribution plan
- [`PAB_fas()`](https://aokine.github.io/mqriskR/reference/PAB_fas.md) :
  Projected annual benefit under a final average salary DB plan
- [`PAB_cae()`](https://aokine.github.io/mqriskR/reference/PAB_cae.md) :
  Projected annual benefit under a career average earnings DB plan
- [`replacement_ratio_db()`](https://aokine.github.io/mqriskR/reference/replacement_ratio_db.md)
  : Replacement ratio for a defined benefit plan
- [`AB_fas()`](https://aokine.github.io/mqriskR/reference/AB_fas.md) :
  Accrued benefit for a final average salary plan
- [`AB_cae()`](https://aokine.github.io/mqriskR/reference/AB_cae.md) :
  Accrued benefit for a career average earnings plan
- [`APV_NR_db()`](https://aokine.github.io/mqriskR/reference/APV_NR_db.md)
  : APV of normal retirement benefit for a DB plan
- [`NC_EAN_db()`](https://aokine.github.io/mqriskR/reference/NC_EAN_db.md)
  : Entry Age Normal normal cost for a DB plan
- [`NC_TUC_db()`](https://aokine.github.io/mqriskR/reference/NC_TUC_db.md)
  : Traditional Unit Credit normal cost for a DB plan
- [`AAL_TUC_db()`](https://aokine.github.io/mqriskR/reference/AAL_TUC_db.md)
  : Traditional Unit Credit accrued liability for a DB plan
- [`NC_PUC_db()`](https://aokine.github.io/mqriskR/reference/NC_PUC_db.md)
  : Projected Unit Credit normal cost for a DB plan
- [`AAL_PUC_db()`](https://aokine.github.io/mqriskR/reference/AAL_PUC_db.md)
  : Projected Unit Credit accrued liability for a DB plan
