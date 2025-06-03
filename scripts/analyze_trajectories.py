import numpy as np
import mpl_toolkits
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

"""
This script is a testing bench for the trajectories function. The goal is to
plot the trajectories that are outputted from Nebula.
"""

def read_tri_file(filename):
    """
    Read a .tri file and return triangles whose first two columns are >= 0.
    Each line: m1 m2 x1 y1 z1 x2 y2 z2 x3 y3 z3
    """
    triangles_geom = []
    triangles_detplane = []
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 11:
                continue
            m1, m2 = int(parts[0]), int(parts[1])
            if m1 >= 0 or m2 >= 0:
                coords = list(map(float, parts[2:]))
                tri = np.array(coords).reshape(3, 3)
                triangles_geom.append(tri)
            elif m1 == -126 and m2 == -126:
                coords = list(map(float, parts[2:]))
                # reshape into 3 points of (x, y, z)
                tri = np.array(coords).reshape(3, 3)
                triangles_detplane.append(tri)
                
    return triangles_geom, triangles_detplane

def plot_triangles(triangles_geom, triangles_detplane, x_range=None, y_range=None, z_range=None):
    """
    Plot a list of triangles in 3D, optionally restricting axis ranges.
    x_range, y_range, z_range: tuple (min, max) or None to auto-scale.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    poly = Poly3DCollection(triangles_geom, alpha=0.7)
    poly.set_edgecolor('k')
    ax.add_collection3d(poly)

    poly = Poly3DCollection(triangles_detplane, alpha=0.5)
    poly.set_edgecolor('r')
    ax.add_collection3d(poly)

    # Compute data bounds
    all_pts = np.vstack(triangles_geom + triangles_detplane)
    xmin, ymin, zmin = all_pts.min(axis=0)
    xmax, ymax, zmax = all_pts.max(axis=0)

    # Apply user-specified ranges or auto-scale
    if x_range:
        ax.set_xlim(x_range)
    else:
        ax.set_xlim(xmin, xmax)
    if y_range:
        ax.set_ylim(y_range)
    else:
        ax.set_ylim(ymin, ymax)
    if z_range:
        ax.set_zlim(z_range)
    else:
        ax.set_zlim(zmin, zmax)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show()

def construct_tree(dat, triangles_geom, triangles_detplane, x_range=None, y_range=None, z_range=None):
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
    ax.scatter(dat[0]['x'],dat[0]['y'],dat[0]['z'], marker='o', color='purple')

    inel_counter = 0
    el_counter = 0
    bc_counter = 0
    term_counter = 0
    
    # Arrays to keep track of edge segments and colors
    segs = []
    colorsVertex = []
    colorsEdge = []

    # Plot the rest
    for i in range(len(dat)):

        # Distinguish between inelastic, elastic, boundary crossing and termination events
        if dat[i]['ce1'] == dat[i]['ce2'] and not dat[i]['ce1'] == -999 and not dat[i]['ce1'] == -998:
            # Elastic is blue-ish
            el_counter += 1
            colorsVertex.append('mediumturquoise')
        elif dat[i]['ce1'] == -dat[i]['ce2'] and not dat[i]['ce1'] == -999 and not dat[i]['ce1'] == -998:
            # Boundary crossing is green
            bc_counter += 1
            colorsVertex.append('lime')
        elif dat[i]['ce1'] == -999:
            # Termination is red
            term_counter += 1
            colorsVertex.append('red')
        elif dat[i]['ce1'] == -998:
            # Detection is yellow
            bc_counter += 1
            colorsVertex.append('forestgreen')
        else:
            # Inelastic is orange
            inel_counter += 1
            colorsVertex.append('orange')

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
            # Plot the detection events with a different colour line
            cl = None
            if dat[i]['ce1'] == -998:
                cl = 'darkgreen'
            else:
                cl = 'teal'

            segs.append(((x_parent, y_parent, z_parent), (x_child, y_child, z_child)))
            colorsEdge.append(cl)

    ln_coll = mpl_toolkits.mplot3d.art3d.Line3DCollection(segs,colors=colorsEdge,linewidth=1)
    ax.add_collection(ln_coll)
    ax.scatter(dat[:]['x'],dat[:]['y'],dat[:]['z'], marker='.', color=colorsVertex)

    # Optional visualization choices:
    # Draw material surface at z = 0
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Build a grid in the x–y plane
    xx = np.linspace(xlim[0], xlim[1], 10)
    yy = np.linspace(ylim[0], ylim[1], 10)
    XX, YY = np.meshgrid(xx, yy)
    ZZ = np.zeros_like(XX)    # z = 0 plane

    ###### Plot the sample geometry from the .tri file
    poly = Poly3DCollection(triangles_geom, alpha=0.2)
    poly.set_edgecolor('k')
    ax.add_collection3d(poly)

    poly = Poly3DCollection(triangles_detplane, alpha=0.1)
    poly.set_edgecolor('r')
    ax.add_collection3d(poly)

    # Compute data bounds
    all_pts = np.vstack(triangles_geom + triangles_detplane)
    xmin, ymin, zmin = all_pts.min(axis=0)
    xmax, ymax, zmax = all_pts.max(axis=0)

    # Apply user-specified ranges or auto-scale
    if x_range:
        ax.set_xlim(x_range)
    else:
        ax.set_xlim(xmin, xmax)
    if y_range:
        ax.set_ylim(y_range)
    else:
        ax.set_ylim(ymin, ymax)
    if z_range:
        ax.set_zlim(z_range)
    else:
        ax.set_zlim(zmin, zmax)

    print(f'Number of inelastic scattering events: {inel_counter}')
    print(f'Number of elastic scattering events: {el_counter}')
    print(f'Number of termination scattering events: {term_counter}')
    
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
    data = np.fromfile("singleline_primaries_at_edge_19nm.det", dtype=trajectory_dtype)

    ### Process and prepare the data
    one_traj = filter_one(data, 1)
    print(one_traj)
    #print(f'How many Lower than 0 nm? {one_traj[one_traj["z"] < 0]}')
    #print(f'Any detected electrons: {one_traj[one_traj["ce1"] == -998]}')
    #print(f'BC events: {one_traj[one_traj["ce1"] == -one_traj["ce2"]]}')
    #print(f'Detection point {one_traj[one_traj["ce1"] == 455]}')

    t_geom, t_detplane = read_tri_file("single_40x20_line.tri")
    construct_tree(one_traj, t_geom, t_detplane, y_range=[-40,40])


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
