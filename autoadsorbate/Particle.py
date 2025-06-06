import numpy as np
from ase import Atoms

def grid_sphere(center, radius, num_theta=36, num_phi=18):
    """
    Generate a grid of points evenly spaced along a sphere surface.

    Parameters:
    - center: tuple or array of shape (3,), the (x, y, z) center of the sphere
    - radius: float, the radius of the sphere
    - num_theta: int, number of divisions along azimuthal angle (longitude)
    - num_phi: int, number of divisions along polar angle (latitude)

    Returns:
    - numpy array of shape (num_phi * num_theta, 3) with Cartesian coordinates
    """
    theta = np.linspace(0, 2 * np.pi, num_theta, endpoint=False)
    phi = np.linspace(0, np.pi, num_phi, endpoint=True)

    theta, phi = np.meshgrid(theta, phi)
    theta = theta.flatten()
    phi = phi.flatten()

    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]

    points = np.stack([x, y, z], axis=-1)
    return points

def fibonacci_sphere(center, radius, point_distance):
    """
    Generate spiral grid of points over a sphere using Fibonacci lattice with approximately even spacing.

    Parameters:
    - center: tuple of (x, y, z), center of the sphere
    - radius: float, radius of the sphere
    - point_distance: float, approximate desired spacing between points on the surface

    Returns:
    - numpy array of shape (n_points, 3) with Cartesian coordinates
    """
    surface_area = 4 * np.pi * radius**2
    approx_area_per_point = point_distance**2
    n_points = int(surface_area / approx_area_per_point)

    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5**0.5) * indices

    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]

    return np.stack((x, y, z), axis=-1)

def move_sphere_points_toward_center(sphere_points, center, point_cloud, touch_criteria, step_size=0.01, max_steps=1000):
    """
    Move each point on the sphere toward the center until it's within `touch_criteria`
    of any point in the point_cloud.

    Parameters:
    - sphere_points: (N, 3) array of points on the sphere
    - center: (3,) tuple or array, the center of the sphere
    - point_cloud: (M, 3) array of reference points
    - touch_criteria: float, distance threshold to stop movement
    - step_size: float, how far to move each step
    - max_steps: int, maximum number of steps to avoid infinite loops

    Returns:
    - (N, 3) numpy array of adjusted sphere points
    """
    from scipy.spatial import KDTree
    sphere_points = np.array(sphere_points, dtype=float)
    center = np.array(center, dtype=float)
    directions = center - sphere_points
    directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]

    tree = KDTree(point_cloud)
    updated_points = sphere_points.copy()

    for i in range(len(updated_points)):
        for _ in range(max_steps):
            dist, _ = tree.query(updated_points[i], k=1)
            if dist < touch_criteria:
                break
            updated_points[i] += directions[i] * step_size

    return updated_points

def random_point_in_sphere(radius, center=(0, 0, 0)):
    """
    Generate a single random 3D point uniformly inside a sphere.

    Parameters:
    - radius: float, radius of the sphere
    - center: tuple of 3 floats, center of the sphere

    Returns:
    - numpy array of shape (3,)
    """
    # Random direction
    vec = np.random.normal(0, 1, 3)
    vec /= np.linalg.norm(vec)

    # Random radius with cube root to ensure uniform volume distribution
    r = radius * np.random.uniform(0, 1) ** (1/3)

    return np.array(center) + vec * r

def random_points_in_sphere(radius, n_points, center=(0, 0, 0)):
    vecs = np.random.normal(0, 1, (n_points, 3))
    vecs /= np.linalg.norm(vecs, axis=1)[:, None]
    rs = radius * np.random.uniform(0, 1, n_points) ** (1/3)
    return np.array(center) + vecs * rs[:, None]


def get_nearby_point_indices(sphere_points, point_cloud, threshold):
    """
    For each sphere point, return the unique indices of point_cloud points within the given threshold.

    Parameters:
    - sphere_points: (N, 3) numpy array of query points (e.g., sphere surface)
    - point_cloud: (M, 3) numpy array of atoms or points in space
    - threshold: float, distance threshold

    Returns:
    - A set of unique indices from point_cloud within threshold of any sphere point
    """
    from scipy.spatial import KDTree
    
    tree = KDTree(point_cloud)
    unique_indices = set()

    for point in sphere_points:
        indices = tree.query_ball_point(point, r=threshold)
        unique_indices.update(indices)

    return np.array(sorted(unique_indices))

def get_all_nearby_indices_per_point(sphere_points, point_cloud, threshold):
    from scipy.spatial import KDTree
    tree = KDTree(point_cloud)
    return [tree.query_ball_point(p, r=threshold) for p in sphere_points]

def keep_unique_inds(inds):
    unique=[]
    for lst in inds:
        
        lst.sort()
        string = ''
        for i in lst:
            string+=f'{i}#'
        string = string[:-1]
        unique.append(string)
    unique = list(set(unique))
    for i, u in enumerate(unique):
        unique[i] = [int(x) for x in u.split('#')] 
    return unique


def mean_of_close_sphere_points(query_point, sphere_points, threshold):
    """
    Find all points on the sphere within `threshold` distance to `query_point`,
    and return their mean position.

    Parameters:
    - query_point: array-like of shape (3,)
    - sphere_points: numpy array of shape (N, 3)
    - threshold: float, distance threshold

    Returns:
    - mean_point: numpy array of shape (3,), or None if no points are within threshold
    """
    from scipy.spatial import KDTree
    tree = KDTree(sphere_points)
    indices = tree.query_ball_point(query_point, r=threshold)
    
    if not indices:
        return None  # No points within threshold

    nearby_points = sphere_points[indices]
    return np.mean(nearby_points, axis=0)


def calculate_sites(inds, particle_atoms, shrinkwrap, threshold=2.7):
    particle = particle_atoms.positions
    
    site_dict = {}
    for k in ['coordinates', 'connectivity', 'topology', 'n_vector', 'h_vector', 'site_formula']:
        site_dict[k]=[]
        
    for i in inds:
        site_dict['topology'].append(i)
        site_dict['connectivity'].append(len(i))
        p1 = np.mean(particle[i], axis=0)
        site_dict['coordinates'].append(p1)
        p2 = mean_of_close_sphere_points(query_point=p1, sphere_points=shrinkwrap, threshold=threshold)
        n_vector = p2-p1
        n_vector /= np.linalg.norm(n_vector)
        site_dict['n_vector'].append(n_vector)
        
        if len(i)==1:
            h_vector = [1.,0,0]
        else:
            h_vector = particle[i[0]] - particle[i[1]]
            h_vector /= np.linalg.norm(h_vector)
        site_dict['h_vector'].append(h_vector)
        
        site_dict['site_formula'].append(particle_atoms[i].symbols.formula.count())
    
    return site_dict

def get_shrinkwrap_particle_ads_sites(
    particle_atoms: Atoms,
    grid_mode = 'fibonacci',
    precision: float = 1.,
    touch_sphere_size: float = 3.
):
    """Identifies adsorption sites on a surface using a shrinkwrap grid.

    Args:
        particle_atoms (Atoms): Atoms slab.
        grid_mode (str): 'fibonacci' or 'grid'
        precision (float): Precision for the shrinkwrap grid.
        touch_sphere_size (float): Radius to consider for grid points.
        return_trj (bool): Whether to return the trajectory for demo mode.

    Returns:
        dict: Dictionary containing site information.
    """
    touch_buffer = 0.2
     
    center = np.mean(particle_atoms.positions, axis=0)
    
    diffs = particle_atoms.positions - center
    dists = np.linalg.norm(diffs, axis=1)
    index = np.argmax(dists)
    particle_radius = dists[index]
    
    grid_radius = particle_radius + touch_sphere_size + 0.5 # 0.5 is safety buffer
    
    if grid_mode == 'fibonacci':
        grid = fibonacci_sphere(center=center, radius=grid_radius, point_distance=precision)
    elif grid_mode == 'grid':
        grid = grid_sphere(center=center, radius=grid_radius, num_theta=36, num_phi=18) # never useful?
    else:
        raise ValueError('grid_mode supported: fibonacci, grid')

    shrinkwrap = move_sphere_points_toward_center(
        sphere_points = grid,
        center = center,
        point_cloud = particle_atoms.positions,
        touch_criteria=touch_sphere_size,
        step_size=0.05,
        max_steps=1000
        )
    
    
    inds = get_all_nearby_indices_per_point(shrinkwrap, particle_atoms.positions, touch_sphere_size+touch_buffer)
    inds = keep_unique_inds(inds)

    sites_dict = calculate_sites(inds,
                                 particle_atoms=particle_atoms,
                                 shrinkwrap=shrinkwrap,
                                 threshold=touch_sphere_size+touch_buffer)
    
    return sites_dict