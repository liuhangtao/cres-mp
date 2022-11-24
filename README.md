# cres-mp
A MATLAB implementation of the crossover residual entropy scaling model incorporating crossover multiparameter equation of state for the viscosity and thermal conductivity of carbon dioxide.

Example:
> [v, l] = co2prop('VL', 'T', 300, 'P', 20e6)
>
>> v = 9.1026e-05
>>   
>> l = 0.1001

Refrence:

[1] H. Liu, F. Yang, Z. Yang, Y. Duan, Crossover residual entropy scaling of the viscosity and thermal conductivity of carbon dioxide, Journal of Molecular Liquids 368 (2022) 120799. doi:10.1016/j.molliq.2022.120799.

[2] F. Yang, Q. Liu, Y. Duan, Z. Yang, Crossover multiparameter equation of state: General procedure and demonstration with carbon dioxide, Fluid Phase Equilibria 494 (2019) 161171. doi:10.1016/j.fluid.2019.04.035.
