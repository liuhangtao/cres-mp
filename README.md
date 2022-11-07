# cres-mp
A MATLAB implementation of the crossover residual entropy scaling model incorporating crossover multiparameter equation of state for the viscosity and thermal conductivity of carbon dioxide.

Example:
>> [vis, tcx] = co2prop('vl','t',300,'p',20e6)
vis =
   9.1026e-05
tcx =
   0.1001

Refrence:
[1] H. Liu, F. Yang, Z. Yang, Y. Duan, Crossover residual entropy scaling of the viscosity and thermal conductivity of 
carbon dioxide, J Mol Liq (submitted).
[2] F. Yang, Q. Liu, Y. Duan, Z. Yang, Crossover multiparameter equation of state: General procedure and demonstration with carbon dioxide, Fluid Phase Equilibr 494 (2019) 161171. URL: http://www.sciencedirect.com/science/article/pii/S0378381219301992. doi:10.1016/j.fluid.2019.04.035.
