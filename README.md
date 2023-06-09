# MultiState@GWEV
The MultiState@GWEV (or MS@GWEV) program is an interface to Gaussian 16 (through the keyword `External`) to calculate spin-forbidden reactions involving multiple spin states.

<img src="https://raw.githubusercontent.com/zorkzou/MultiState/master/mssm-logo.png" />

## Recent Changes

04/20/2023

1. Bug fix for the default `chi` value.

2. Bug fix for the `ONIOM` calculation by Gaussian 16.c.

3. Single spin state calculation may be performed through `-nst 1`.

## Citation

Please cite the following paper if you use this program.

* L. Zhao and W. Zou, A general method for locating stationary points on the mixed-spin surface of spin-forbidden reaction with multiple spin states, J. Chem. Phys. 158, 224110 (2023). [doi](https://doi.org/10.1063/5.0151630)

The closely related two-state spin-mixing model may be found in the following papers.

* B. Yang, L. Gagliardi, and D. G. Truhlar, Transition states of spin-forbidden reactions, Phys. Chem. Chem. Phys. 20, 4129 (2018). [doi](https://doi.org/10.1063/5.0151630)
* T. Takayanagi and T. Nakatomi, Automated reaction path searches for spin-forbidden reactions, J. Comput. Chem. 39, 1319 (2018). [doi](https://doi.org/10.1002/jcc.25202)

