# helper_functions.py
import numpy as np
import io
#from scipy.spatial.distance import pdist, squareform
import numpy as np

def pairwise_distances(coords):
   
    n = len(coords)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            distances[i, j] = np.linalg.norm(coords[i] - coords[j])
    return distances
# Replace `pdist` and `squareform` in `infer_bonds` with:
distances = pairwise_distances(atomic_coordinates)


def read_xyz(file):
    content = file.getvalue().decode("utf-8").splitlines()
    lines = content[2:]  # Skipping the first two lines
    atomic_symbols = [line.split()[0] for line in lines]
    atomic_coordinates = np.array(
        [line.split()[1:4] for line in lines], dtype=np.float64
    )
    return atomic_symbols, atomic_coordinates
# Replace `pdist` and `squareform` in `infer_bonds` with:
distances = pairwise_distances(atomic_coordinates)



def write_xyz(atomic_symbols, atomic_coordinates):
    output = io.StringIO()
    output.write(f"{len(atomic_symbols)}\n")
    output.write("Modified molecule\n")
    for symbol, (x, y, z) in zip(atomic_symbols, atomic_coordinates):
        output.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")
    return output.getvalue()


def replace_atom_with_group(
    atomic_symbols, atomic_coordinates, atom_index, new_group, bond_length=1.0
):
    # Remove the selected atom
    del atomic_symbols[atom_index]
    removed_atom_coords = atomic_coordinates[atom_index]
    atomic_coordinates = np.delete(atomic_coordinates, atom_index, axis=0)

    # Find the nearest neighbor (assuming it's bonded)
    distances = np.linalg.norm(atomic_coordinates - removed_atom_coords, axis=1)
    nearest_neighbor_index = np.argmin(distances)
    nearest_neighbor_coords = atomic_coordinates[nearest_neighbor_index]

    # Calculate the bond vector
    bond_vector = removed_atom_coords - nearest_neighbor_coords
    bond_vector /= np.linalg.norm(bond_vector)


    # Add new group atoms
    for i, new_atom in enumerate(new_group):
        new_coords = nearest_neighbor_coords + (i + 1) * bond_length * bond_vector
        atomic_symbols.append(new_atom)
        atomic_coordinates = np.vstack([atomic_coordinates, new_coords])

    return atomic_symbols, atomic_coordinates


def create_xyz_string(atomic_symbols, atomic_coordinates):
    xyz_string = f"{len(atomic_symbols)}\nModified molecule\n"
    for symbol, coords in zip(atomic_symbols, atomic_coordinates):
        xyz_string += f"{symbol} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}\n"
    return xyz_string


def add_group_to_atom(
    atomic_symbols, atomic_coordinates, atom_index, new_group, bond_length=1.0
):
    base_atom_coords = atomic_coordinates[atom_index]

    # Find the nearest neighbor (assuming it's bonded)
    distances = np.linalg.norm(atomic_coordinates - base_atom_coords, axis=1)
    distances[atom_index] = np.inf  # Exclude the atom itself
    nearest_neighbor_index = np.argmin(distances)
    nearest_neighbor_coords = atomic_coordinates[nearest_neighbor_index]

    # Calculate the bond vector (pointing away from the nearest neighbor)
    bond_vector = base_atom_coords - nearest_neighbor_coords
    bond_vector /= np.linalg.norm(bond_vector)

    # Add new group atoms
    for i, new_atom in enumerate(new_group):
        new_coords = base_atom_coords + (i + 1) * bond_length * bond_vector
        atomic_symbols.append(new_atom)
        atomic_coordinates = np.vstack([atomic_coordinates, new_coords])

    return atomic_symbols, atomic_coordinates


def delete_atoms(atomic_symbols, atomic_coordinates, atom_indices):
    # Sort indices in reverse order to avoid shifting problems
    for index in sorted(atom_indices, reverse=True):
        del atomic_symbols[index]
        atomic_coordinates = np.delete(atomic_coordinates, index, axis=0)
    return atomic_symbols, atomic_coordinates

def infer_bonds(atomic_symbols, atomic_coordinates, bond_threshold=1.6):
    """
    Infers bonds between atoms based on distances.
    Returns a list of tuples representing bonds and their inferred order.
    """
    bond_list = []
    bond_orders = []
    
    # Distance matrix for atomic coordinates
    distances = squareform(pdist(atomic_coordinates))
    
    # Covalent radii approximation (Angstroms, sourced from standard values)
    covalent_radii = {
        "H": 0.31, "C": 0.76, "O": 0.66, "N": 0.71, "F": 0.57,
        "Cl": 0.99, "Br": 1.14, "I": 1.33, "S": 1.02, "P": 1.07,
    }
    
    # Iterate over atom pairs to identify bonds
    for i in range(len(atomic_coordinates)):
        for j in range(i + 1, len(atomic_coordinates)):
            expected_bond_length = covalent_radii.get(atomic_symbols[i], 0.5) + covalent_radii.get(
                atomic_symbols[j], 0.5)
            margin = 0.2  # Adjust margin depending on desired accuracy
            if distances[i, j] <= expected_bond_length + margin:
                bond_list.append((i, j))
                # Infer bond order based on bond length (simplified logic)
                if distances[i, j] < expected_bond_length * 0.9:
                    bond_orders.append(2)  # Double bond
                else:
                    bond_orders.append(1)  # Single bond
                
    return bond_list, bond_orders
