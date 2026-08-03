# Usage

There are currently 5 commands one can run in order to use NebulaLE:

- `nebula_gpu`

- `nebula_cpu_mt`

- `nebula_cpu_edep`

- `nebula_cpu_traj`

- `nebula_cpu_SE_mode`

## Commands & Features

### GPU version: `nebula_gpu`
Run Nebula on the GPU. Not all changes in NebulaLE are guaranteed to work for GPU, so use of this version is not recommended.

### CPU main version: `nebula_cpu_mt`
Run Nebula on the CPU. This is the standard version of Nebula.

### CPU edep version: `nebula_cpu_edep`
Run Nebula on the CPU, with a custom output datatype that logs energy deposition events. Detection events (if detector is present)
are stored in a .bin file instead of a .det file. 

The output datatype is given by:

```python
import numpy as np

energy_dtype = np.dtype([ 
('x', '=f'), 
('y', '=f'), 
('z', '=f'), # Position 
('E', '=f'), # Energy of electron before collision at position x, y, z (in eV) 
('dE', '=f'), # Energy loss (deV) 
('px', '=i'), # Pixel x of initial primary beam electron 
('py', '=i')]) # Pixel y of initial primary beam electron 
```

(Current bug: total energy outputted is too high: Secondary electron energy is double counted.)

### CPU trajectories version: `nebula_cpu_traj`
Run Nebula on the CPU, with a custom output datatype that logs every single event that the electron encounters (inelastic,
elastic, boundary crossing (transmission and reflection), mirroring (against a mirror plane), detection, termination. IDs of each event type are
custom, and require a special script to interpret.

The output datatype is given by:

```python
import numpy as np

trajectory_dtype = np.dtype([ 
('x', '=f'), 
('y', '=f'), 
('z', '=f'), # Position 
('E', '=f'), # Energy of electron before collision at position x, y, z (in eV) 
('dE', '=f'), # Energy loss (deV) 
('px', '=i'), # Pixel x of initial primary beam electron 
('py', '=i'), # Pixel y of initial primary beam electron 
('tag', '=i'), # Unique ID tag of initial primary beam electron  
('ns', '=i')]) # Consecutive scattering event number
```

(Current bug: blocks GPU compilation, there seem to be some inelastic events missing in the list of output events.)


### CPU SE version: `nebula_cpu_SE_mode`
Run Nebula on the CPU, with an additional output of the scatter counter. This is distinct from the trajectories version because it focuses on 
detected electrons and not "events" in general.

The output datatype is given by:

```python
import numpy as np

SE_mode_dtype = np.dtype([ 
('x', '=f'), 
('y', '=f'), 
('z', '=f'), # Position 
('dx', '=f'), 
('dy', '=f'), 
('dz', '=f'), # Direction 
('E', '=f'), # Energy 
('sc', '=i'), # Scatter counter 
('px', '=i'), 
('py', '=i')]) # Pixel index
```

(Current bug: blocks GPU compilation)