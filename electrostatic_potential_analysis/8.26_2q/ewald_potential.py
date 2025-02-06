import numpy as np
import pandas as pd
from math import sqrt, pi, exp, log
from tqdm import tqdm
from scipy.special import erfc

# Constants
DELTA = 0.125  # Grid spacing in Å
AXIAL_LENGTH = 90.0  # May not be used with Ewald, but kept for consistency
EPSILON = 8.8541878176e-12
ANG_TO_METER_CONV_FACTOR = 1e-10
COULOMB_PER_UNIT_CHARGE = 1.60217733e-19
K = (1/(4 * pi * EPSILON * ANG_TO_METER_CONV_FACTOR)) * COULOMB_PER_UNIT_CHARGE

def write_results_to_csv(potentials, filename):
    """Write grid potentials to a CSV file."""
    df = pd.DataFrame(potentials, columns=["X", "Y", "Z", "Potential","Ex","Ey","Ez"])
    df.to_csv(filename, index=False)

def parse_cif(cif_file):
    """Parse atomic data from a CIF file."""
    with open(cif_file, 'r') as f:
        lines = f.readlines()

    a, b, c = 0, 0, 0
    atomic_data = []

    for line in lines:
        line = line.strip()
        if line.startswith("_cell_length_a"):
            a = float(line.split()[-1])
        elif line.startswith("_cell_length_b"):
            b = float(line.split()[-1])
        elif line.startswith("_cell_length_c"):
            c = float(line.split()[-1])

    start_idx = None
    for i, line in enumerate(lines):
        if "_atom_site_label" in line:
            start_idx = i + 1
            break

    if start_idx is None:
        raise ValueError("_atom_site_label not found in the CIF file.")

    for line in lines[start_idx:]:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("_"):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue

        label = parts[0]
        atom_type = parts[1]
        fract_x = float(parts[2])
        fract_y = float(parts[3])
        fract_z = float(parts[4])
        charge = float(parts[5])

        atomic_data.append({
            "label": label,
            "atom_type": atom_type,
            "fract_x": fract_x,
            "fract_y": fract_y,
            "fract_z": fract_z,
            "charge": charge
        })

    return atomic_data, a, b, c

def fractional_to_cartesian_with_images(atomic_data, a, b, c):
    """Convert fractional coordinates to Cartesian coordinates and include ±1 periodic images."""
    
    # Convert fractional to cartesian for the original cell
    cartesian_data = []
    for atom in atomic_data:
        x = atom["fract_x"] * a
        y = atom["fract_y"] * b
        z = atom["fract_z"] * c
        cartesian_data.append({
            "label": atom["label"],
            "atom_type": atom["atom_type"],
            "charge": atom["charge"],
            "x": x,
            "y": y,
            "z": z
        })
    
    # Create periodic images in ±1 directions
    # For each atom, we will generate images by adding or subtracting a, b, c for each coordinate
    expanded_data = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                for atom in cartesian_data:
                    # Translate atom coordinates by dx*a, dy*b, dz*c
                    x_image = atom["x"] + dx * a
                    y_image = atom["y"] + dy * b
                    z_image = atom["z"] + dz * c
                    expanded_data.append({
                        "label": atom["label"],
                        "atom_type": atom["atom_type"],
                        "charge": atom["charge"],
                        "x": x_image,
                        "y": y_image,
                        "z": z_image
                    })
    
    return expanded_data

def fractional_to_cartesian(atomic_data, a, b, c):
    """Convert fractional coordinates to Cartesian coordinates."""
    for atom in atomic_data:
        atom["x"] = atom["fract_x"] * a
        atom["y"] = atom["fract_y"] * b
        atom["z"] = atom["fract_z"] * c
    return atomic_data

def generate_grid(a, b, diameter, starting_z, z_delta, z_layers):
    """Generate a 3D grid."""
    center_x = 0.5 * a
    center_y = 0.5 * b
    radius = diameter / 2

    x_range = np.arange(center_x - radius, center_x + radius, DELTA)
    y_range = np.arange(center_y - radius, center_y + radius, DELTA)
    z_range = [starting_z + i * z_delta for i in range(z_layers)]

    grid = []
    for x in x_range:
        for y in y_range:
            for z in z_range:
                grid.append((x, y, z))
    return np.array(grid)

def compute_alpha(r_cut, epsilon=1e-6):
    """Compute the Ewald parameter alpha from cutoff and precision."""
    return np.sqrt(-np.log(epsilon)) / r_cut

def compute_structure_factors(atomic_data, k_vectors):
    """Compute S(k) = sum_j q_j e^{-i k·r_j}."""
    charges = np.array([atom["charge"] for atom in atomic_data])
    positions = np.array([[atom["x"], atom["y"], atom["z"]] for atom in atomic_data])

    Sk = np.zeros(len(k_vectors), dtype=complex)
    for idx, k_vec in enumerate(k_vectors):
        phase = np.exp(-1j * np.dot(positions, k_vec))
        Sk[idx] = np.sum(charges * phase)
    return Sk

def calculate_potential_ewald(grid, atomic_data, a, b, c, k_vectors, r_cut, epsilon=1e-6):
    """
    Calculate electric potential at each grid point using Ewald summation.
    Includes K in both real-space and reciprocal-space terms.
    """
    # Compute alpha
    alpha = compute_alpha(r_cut, epsilon)
    print(alpha)
    V = a * b * c  # Volume of the cell

    # Precompute structure factors for reciprocal space
    Sk = compute_structure_factors(atomic_data, k_vectors)

    positions = np.array([[atom["x"], atom["y"], atom["z"]] for atom in atomic_data])
    charges = np.array([atom["charge"] for atom in atomic_data])

    # Real-space: consider only atoms within r_cut of each grid point
    # This requires no explicit periodic expansions if r_cut < min(a,b,c),
    # or if you include nearest neighbor images if needed.

    potentials = []
    total_points = len(grid)
    with tqdm(total=total_points, desc="Calculating potentials (Ewald)") as pbar:
        for gx, gy, gz in grid:
            r_grid = np.array([gx, gy, gz])

            # --- Real-space contribution ---
            diff = positions - r_grid
            dist = np.sqrt(np.sum(diff**2, axis=1))
            mask = (dist > 1e-12) & (dist <= r_cut)
            dist_masked = dist[mask]
            q_masked = charges[mask]
            diff_masked = diff[mask]

            # Potential contribution
            # phi_real = sum_j K * q_j * erfc(alpha*r)/r
            real_contrib = np.sum(q_masked * erfc(alpha * dist_masked) / dist_masked) * K

            # Field contributions
            # For simplicity, we use the approximate factor from the original code:
            # factor = K * q_j * erfc(alpha*r)/r^3
            factor = (K * q_masked * erfc(alpha * dist_masked)) / (dist_masked**3)
            ex = np.sum(factor * diff_masked[:,0])
            ey = np.sum(factor * diff_masked[:,1])
            ez = np.sum(factor * diff_masked[:,2])

            # --- Reciprocal-space contribution ---
            # phi_rec(r) = K * (4*pi/V)* sum_{k!=0} [exp(-k^2/(4 alpha^2))/k^2 * S(k)*exp(i k·r)]
            phi_rec = 0.0
            for idx, k_vec in enumerate(k_vectors):
                k_sq = np.dot(k_vec, k_vec)
                if k_sq < 1e-14:
                    continue
                phase_grid = np.exp(1j * np.dot(k_vec, r_grid))
                factor_rec = (4.0 * pi / V) * exp(-k_sq/(4.0 * alpha**2)) / k_sq
                phi_rec += factor_rec * Sk[idx] * phase_grid

            phi_rec = phi_rec.real * K

            potential = real_contrib + phi_rec
            potentials.append((gx, gy, gz, potential, ex, ey, ez))
            pbar.update(1)

    return potentials

def main(cif_file, output_file, cylindrical_diameter, starting_z, z_delta, z_layers, k_vectors, r_cut):
    # Parse CIF file
    atomic_data, a, b, c = parse_cif(cif_file)
    atomic_data = fractional_to_cartesian_with_images(atomic_data, a, b, c)

    print(f"Cylindrical diameter: {cylindrical_diameter} Å")
    print(f"z_delta: {z_delta}, z_layers: {z_layers}")

    # Generate grid
    grid = generate_grid(a, b, cylindrical_diameter, starting_z, z_delta, z_layers)
    print(f"Grid points: {len(grid)}")

    # Calculate potentials using Ewald
    potentials = calculate_potential_ewald(grid, atomic_data, a, b, c, k_vectors, r_cut)

    # Write results to file
    write_results_to_csv(potentials, output_file)
    print(f"Potential data written to {output_file}")

if __name__ == "__main__":
    CIF_FILE = "AR/IC_24_AA.cif"
    OUTPUT_FILE = "potential_data.csv"
    PROBE_LENGTH = 16
    STARTING_Z = 0.55
    Z_DELTA = 0.11
    Z_LAYERS = 20
    # Example k_vectors and cutoff:
    a, b, c = 30.8, 30.8, 30.8  
    KX, KY, KZ = 0, 0, 9
    k_vectors = []
    for i in range(-KX, KX+1):
        for j in range(-KY, KY+1):
            for k in range(-KZ, KZ+1):
                if i == 0 and j == 0 and k == 0:
                    continue
                kx = (2.0 * pi * i) / a
                ky = (2.0 * pi * j) / b
                kz = (2.0 * pi * k) / c
                k_vectors.append([kx, ky, kz])
    k_vectors = np.array(k_vectors)

    r_cut = 14.0  # Example real-space cutoff in Å

    main(CIF_FILE, OUTPUT_FILE, PROBE_LENGTH, STARTING_Z, Z_DELTA, Z_LAYERS, k_vectors, r_cut)

