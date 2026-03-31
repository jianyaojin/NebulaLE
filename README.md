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
`
trajectory_dtype = np.dtype([
        ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
        ('E',  '=f'), # Energy of electron before collision at position x, y, z (in eV)
        ('dE',  '=f'), # Energy loss (deV)
        ('tag', '=i'), # Tag of primary
        ('pe', '=i'),  # Parent edge id
        ('ce1', '=i'), # Child edge one
        ('ce2', '=i')])# Child edge two 
`

Also included is a complementary script which processes the trajectories output data and plots them using matplotlib. The matplotlib plot is very slow since there are a lot of edges.

### Update 20-05-2025

Added functionality to track termination events (output from Nebula + analyze_trajectories.py). Termination events have no child edges and so have both ce1 and ce2 equal to -999 (placeholder value, may change).

### Update 27-05-2025

Added functionality to track both boundary crossing and detection events. 

- Boundary crossing events are similar to elastic scattring events. They can be identified by child tags with opposite signs, so for instance: ce1 = a, ce2 = -a. The current implementation assumes that it is impossible for electrons to be absorbed at a material interface. In other words, it assumes there is always a child edge and thus at least one more logged event after it.
- Detection events are logged with the child edge ID -998 for both ce1 and ce2. Detection events here differ from those logged by the detector. Here we specifically log the boundary intersection event with a detector, after which the particle would normally be detected and saved.

The python trajectory analysis script has been updated accordingly as well.


### Update 23-09-2025

Implemented the `nebula_cpu_SE_mode` functionality. The datatype here is given by:
SE_mode_dtype = np.dtype([
            ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
            ('dx', '=f'), ('dy', '=f'), ('dz', '=f'), # Direction
            ('E',  '=f'),                             # Energy
            ('sc', '=i'),                             # Scatter counter
            ('px', '=i'), ('py', '=i')])              # Pixel index

The only difference between SE mode and `nebula_cpu_mt` is that SE mode also outputs the "scatter counter", which counts the number of scattering events undergone by the electron before being detected (elastic and inelastic both). When a secondary is generated, the secondary inherits the scatter number of its parent electron (This is after it has been updated due to the scattering event)

The use of this functionality is to be able to separate between SE1 electrons and SE2 electrons. The exact difference between the two is not clearly defined so that is left to the user to decide. 