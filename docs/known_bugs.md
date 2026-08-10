# Known Bugs

This is a list of known issues that are still in the Nebula or cstool code.

## Nebula:

- Energy deposition mode double counts SE energies, meaning that if you sum up the energy of all energy deposition events dE you get a number that is
higher than the number of primaries $\times$ the primary landing energy

- Incorrectly sampled q values for inelastic scattering cause a very small number of scattering angle costheta values to be larger than 1. This is
likely due to an interpolation error near 0 of the icdf table.

- The trajectories function does not output all SE generation events, for some reason.

- Setting the minimum simulation energy to a negative number in order to lower the simulation cutoff threshold energy results in an unexpectedly large
number of 0 energy loss edep events

## cstool: