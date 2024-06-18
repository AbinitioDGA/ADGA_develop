# AbinitioDGA development

This repository includes the addition of the calculation of the current-current susceptibility and the calculation of the resolved self-energy
to perform fluctuation diagonstics.

The development of these branches were never inteneded for public use and are hence slightly cumbersome to use.


## Current-current susceptibility calculation

In order to calculate the optical conductivity, first a full ADGA run has to be performed (either oneshot or via the self-consistent extension -- which is not part of this package).

From these new k-resolved self-energies we compute the DGA Green's function with the auxiliary script in scripts/compute\_GDGA.py.

Further, the derivative of the Hamiltonian has be calculated, whose format is defined in scripts/derive.py

For the optical run, the config file needs to be extended to include the derivative of the Hamiltonian, the newly created Green's function file and additional frequency parameters that determine the internal frequency sums.
(Note that some of these flags in the updated configspec file are outdated)

The optical run is configured by the following settings:

`calc-cond = T`

`N1iwbc = X` (where X is the number of positive bosonic frequencies for chi\_jj)

`HkdkFile = X` (derivative of the Hamiltonian -- created by e.g. the derive.py file)

`cond-phhor = T` (if the ph-contribution of vertex corrections should be calculated)

`cond-phver = T` (if the transverse ph-contribution of vertex corrections should be calculated)

`cond-legs = X` (where X is the output file from the compute\_GDGA.py script)

`extend-cond-bubble = T` (if the bubble should be calculated on the largest possible frequency grid)

`QDataPHBAR = qdata` (a text file with '0 0 0' in it to signal that we want the zero momentum transfer response)

## Nota bene
The computation for q=0 has been validated with an independent code package developed in the group.
I, however, cannot guarantee the validity of any calculations of chi\_jj at finite q.

The optical conducitivity calculation requires (contrary to susceptibility calculations) a 'full run', i.e. the full q-grid and the full w-box should be used.

The output of the current-current susceptibility will be in the HDF5 under the path /conductivity
and will only include the positive frequencies.

For more details please refer to the Dissertation "Electronic correlation and transport phenomena in Mott and Kondo materials" by Matthias Pickem
