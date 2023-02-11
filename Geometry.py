import numpy as np

vanderwaals_radii = {  # Angstrom
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'P': 1.8,
    'S': 1.8,
    'Cl': 1.75,
    'Br': 1.85,
    'I': 1.98,
    'B': 1.85,
    'Si': 2.1,
    'Se': 1.9,
    'Te': 2.06,
    'Zn': 1.39,
    'Cu': 1.4,
}


def distance_matrix(atoms: np.ndarray) -> np.ndarray:
    """Calculate the distance matrix between atoms."""
    dm = np.zeros((len(atoms), len(atoms)))
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            dm[i, j] = np.linalg.norm(atoms[i] - atoms[j])
    return dm


def distance_cutoff(atoms: np.ndarray, cutoff: float) -> np.ndarray:
    """Calculate the distance matrix between atoms and apply a cutoff."""
    dm = distance_matrix(atoms)
    dm[dm > cutoff] = 0
    return dm


def get_vanderwaals_radii(atoms: list) -> np.ndarray:
    """Get the van der Waals radii for the atoms."""
    radii = np.zeros(len(atoms))
    for i, atom in enumerate(atoms):
        radii[i] = vanderwaals_radii[atom]
    return radii


def get_edges_from_distance_matrix(
        atoms: np.ndarray,
        atom_names: list) \
        -> list:
    """Get the edges from the distance matrix."""
    dm = distance_matrix(atoms)
    radii = get_vanderwaals_radii(atom_names)
    edges = []
    for i in range(len(dm)):
        for j in range(len(dm)):
            if i < j and dm[i, j] < (radii[i] + radii[j])*0.6:
                edges.append((i, j))
    return edges