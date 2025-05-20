import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D 

"""
This script is a testing bench for the trajectories function. The goal is to
plot the trajectories that are outputted from Nebula.
"""

def construct_tree(dat):
    """
    Function takes data from Nebula output, plots trajectories by going down the list of
    energy deposition nodes. It does this by:
    - Plotting each node first
    - Looking at the node parent ID, finding the corresponding node with that parent ID as
    a child ID
    - Plotting the edge between parent and child.

    Outputs a plot of the electron trajectories. This functionality does not visualize:
    - boundary crossing events
    - Electron termination events
    """
    print(dat)
    
    # Set up plot environment, emphasize root node
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dat[0]['x'],dat[0]['y'],dat[0]['z'], marker='o', color='lime')

    # Plot the rest
    for i in range(len(dat)):
        
        # Distinguish between inelastic, elastic and termination events
        if dat[i]['ce1'] == dat[i]['ce2'] and not dat[i]['ce1'] == -999:
            # Elastic is blue-ish
            ax.scatter(dat[i]['x'],dat[i]['y'],dat[i]['z'], marker='.', color='mediumturquoise')
        elif dat[i]['ce1'] == -999:
            # Termination is red
            ax.scatter(dat[i]['x'],dat[i]['y'],dat[i]['z'], marker='.', color='red')
        else:
            # Inelastic is orange
            ax.scatter(dat[i]['x'],dat[i]['y'],dat[i]['z'], marker='.', color='orange')

        # Plot the parent edge
        p_edge = dat[i]['pe']
        
        if p_edge != 0:
            # Find the node with p_edge as the child edge:
            mask = (dat['ce1'] == p_edge) | (dat['ce2'] == p_edge)
            matches = dat[mask]
            index = np.nonzero(mask)[0]
            #print(f'Found that node {i} is the child of node {index[0]}')
            
            x_parent, y_parent, z_parent = dat[index[0]][['x','y','z']]
            x_child, y_child, z_child = dat[i][['x','y','z']]

            #print(f'Coords check: {x_parent, y_parent, z_parent} and child {x_child, y_child, z_child}')
            #raise NotImplementedError
            ax.plot([x_parent, x_child], [y_parent, y_child], [z_parent, z_child], color='teal', linewidth=1)

    # Optional visualization choices:
    # Draw material surface at z = 0
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Build a grid in the x–y plane
    xx = np.linspace(xlim[0], xlim[1], 10)
    yy = np.linspace(ylim[0], ylim[1], 10)
    XX, YY = np.meshgrid(xx, yy)
    ZZ = np.zeros_like(XX)    # z = 0 plane

    # Plot a semitransparent brown surface
    ax.plot_surface(XX, YY, ZZ,
                    color='saddlebrown',   # or simply 'brown'
                    alpha=0.4,             # adjust transparency (0=fully transparent, 1=opaque)
                    linewidth=0,           # no mesh lines
                    antialiased=True)
    
    # Clean up and show plot    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')    
    plt.show()


def filter_one(dat, prim=0):
    """
    Function filters the input data and returns the dataset of a single primary
    electron's cascade. Can specify which electron you specifically want.
    """
    return dat[dat['tag']==prim]

def filter_primary(dat):
    """
    Filters the data so that only the primary trajectory is shown. Since Nebula
    fully traverses the primary path before moving onto secondary paths, we just
    need to cut off the data at the point where we see a parent edge ID of 2
    (assuming that an inelastic scattering event occurs first)
    """
    for v in range(len(dat)):
        if dat[v]['pe'] == 2:
            break
    return dat[:v+1]  

    
if __name__ == "__main__":

    ### Define the two datatypes: Trajectory datatype and edep datatype.
    trajectory_dtype = np.dtype([
        ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
        ('E',  '=f'), # Energy of electron before collision at position x, y, z (in eV)
        ('dE',  '=f'), # Energy loss (deV)
        ('tag', '=i'), # Tag of primary
        ('pe', '=i'),  # Parent edge id
        ('ce1', '=i'), # Child edge one
        ('ce2', '=i')])# Child edge two 

    energy_dtype = np.dtype([
        ('x',  '=f'), ('y',  '=f'), ('z',  '=f'), # Position
        ('E',  '=f'), # Energy of electron before collision at position x, y, z (in eV)
        ('dE',  '=f'), # Energy loss (deV)
        ('px', '=i'), ('py', '=i')]) # tags of primary electron

    ### Read output file
    data = np.fromfile("output_with_terminated.det", dtype=trajectory_dtype)

    ### Process and prepare the data
    one_traj = filter_one(data)
    print(one_traj)
    #one_traj = filter_primary(one_traj)
    #coordinates = one_traj[['x', 'y', 'z']]
    #print(f'How many scattering events do we have? {len(one_traj)}')
    #print(f'One trajectory data {one_traj}')
    #print(f'What are the coordinates? {coordinates}')

    construct_tree(one_traj)
    """    
    edges = prepare_edges(one_traj[['pe','ce1','ce2']])
    print(f'Edges output from the prepare_edges function {edges}')

    edges = edges[:4]
    print(f'Raw data {data[:100]}')
    print(f'Coordinates at positions 0, 71 in data {coordinates[0], coordinates[71]}')
    print(f'Coordinates are {coordinates[:4]}')
    print(f'Edges are {edges}')

        
    plot_trajectories(edges, coordinates)
    """
    raise NotImplementedError



"""
Legacy Code:

def prepare_edges(dat):
    """"""
    Convert raw node–edge data into a list of (parent, child) index pairs.

    Parameters
    ----------
    data : numpy.ndarray
        Structured array of length N, with integer fields:
          - 'pe'  : parent‐edge ID for node i
          - 'ce1' : first child‐edge ID for node i
          - 'ce2' : second child‐edge ID for node i
        Each entry’s index (0..N-1) is the node’s numeric ID.

    Returns
    -------
    edges : list of tuple(int, int)
        All (parent_node, child_node) pairs, deduplicated.
    """"""
    # Inspect data
    print(dat)

    # Edges array is an array of lists of format [A,B], which tells us that
    # node A connects to node B with an edge (A and B are some integer that
    # represents a node ID). We generate this array by checking each entry
    # in the input data and leaving the node pair (e.g. B in the above
    # example) empty, to be filled later.
    edges = []

    # Matchmaker list is an ordered list of edge IDs. The index of a value
    # corresponds to the parent
    # Note to self: Matchmaker could maybe be replaced by np.arange!
    matchmaker = []

    # Traverse the data, fill in the two lists (edges only half filled)
    for i in range(len(dat)):
        if dat[i]['ce1'] != dat[i]['ce2']:
            edges.append([i,None])
            matchmaker.append(dat[i]['ce1'])
            edges.append([i,None])
            matchmaker.append(dat[i]['ce2'])
        else:
            edges.append([i,None])
            matchmaker.append(dat[i]['ce1'])
    print(f'-------Edges looks like (Node1, Node2) {edges}')
    print(f'-------The Matchmaker looks like {matchmaker}')

    # Find the matching node for each edge in the edges list by looking up
    # the matchmaker edge ID

    for i in range(len(matchmaker)):
        # Every node only has one parent, we just need to find the correct one
        for j in range(len(dat)):
            if dat[j]['pe'] == matchmaker[i]:
                print(f'Comparing {dat[j]} and {matchmaker[i]}')
                #print(f'========== We have found an edge between node {matchmaker[i]} and {j}')
                break
        edges[i][1] = j

    #print(f'After filling, the edges are {edges}')
    #print(f'Entry 911 is {dat[911]}')
    
    return edges
"""

