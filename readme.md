# STEADY STATE TURBULENT SOLVER (WIP)

This code models channel flow with present turbulence phenomena. The turbulence model used is the k-epsilon one, and first order upwinding is used for all transport phenomena. Moreover, wall functions are inluded and the solver is based on a collocatedgrid

- Written in Fortran 90
- Valid for uniform, orthonomal grid.
- Use of k-epsilon model to simulate mean flow characteristics.
- Postprocessing in Paraview

The example code has been setup, tested, and validated for a simple channel flow. However, it is (propably!) able to model different problems such as lid-driven cavity.

