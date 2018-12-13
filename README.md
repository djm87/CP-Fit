# CP-Fit
Automated fitting of stress-strain data and texture for crystal plasticity codes such as VPSC, EPSC, CPFE, and VPFFT

## Overview
Crystal plasticity and more specifically the implementations that account for twinning and dislocation evolution have several dozen parameters. For a skilled and knowledgable person, a robust dataset that constrains these parameters requires weeks to months to configure. 

This project couples Matlab's parallel environement, global optimization tool kit, and MTEX compatibility to quickly fit large datasets. The cost function currently handles fitting of the stress-strain response and will be extended in the future to include volume fraction of main texture components, general texture difference, and quantities such as twin fractions when determining the global optimum parameters.
