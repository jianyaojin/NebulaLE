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

# Update Log:
### Update 14-05-2025
Added Nebula trajectories function. This function views the electron scattering path as a tree, where each energy deposition event is a node and the path traversed by the electron between nodes are edges.

This function can be operated via "nebula_cpu_traj", which tracks the electron trajectories by keeping track for each node:
- The parent edge ID
- Child edge 1's ID
- Child edge 2's ID

The edge ID tracking piggybacks off of the electron particle class/type. This means that while each edge in the tree has a unique ID, the same electron may (and will) carry multiple different edge IDs across its lifetime.
The data type outputted by this function is:
trajectory_dtype = np.dtype([
        ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
        ('E',  '=f'), # Energy of electron before collision at position x, y, z (in eV)
        ('dE',  '=f'), # Energy loss (deV)
        ('tag', '=i'), # Tag of primary
        ('pe', '=i'),  # Parent edge id
        ('ce1', '=i'), # Child edge one
        ('ce2', '=i')])# Child edge two 

Currently the function does not track the following events:
- Boundary crossing events
- Electron termination events

Also included is a complementary script which processes the trajectories output data and plots them using matplotlib. The matplotlib plot is very slow since there are a lot of edges.