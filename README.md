# NebulaLE (Nebula Low Energy)

This is a follow up project that continues development of Nebula, Monte Carlo simulator for electron-matter interaction. 
The main simulator repository can be found [here](https://github.com/Nebula-simulator/Nebula).
The documentation for the simulator can be found [here](https://nebula-simulator.github.io).

The focus of this version is to improve the physics models of the low energy range < 500 eV, with special focus between 0 and 100 eV.


## Some Documentation Points:

Energy deposition function has dtype:

energy_dtype = np.dtype([
    ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
    ('E',  '=f'), # Energy of electron before collision at position x, y, z (in eV)
    ('dE',  '=f'), # Energy loss (deV)
    ('px', '=i'), ('py', '=i')]) # tags of primary electron
	
Trajectory function has dtype:

trajectory_dtype = np.dtype([
    ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
    ('E',  '=f'), # Energy of electron before collision at position x, y, z (in eV)
    ('dE',  '=f'), # Energy loss (deV)
    ('px', '=i'), ('py', '=i'), # Pixel tag of primary
    ('tag', '=i'), # Tag of primary
    ('ns', '=i')]) # Consecutive scattering event number
