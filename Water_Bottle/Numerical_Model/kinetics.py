import numpy as np
from constants import *

# Rotational inertia of water 
def rotational_water(positions, center_of_mass, water_mass, n_slices):
    positions = np.asarray(positions, dtype=float)
    rel = positions - center_of_mass
    slice_mass = water_mass / n_slices
    return slice_mass * np.sum(rel ** 2)

# Rotational inertia of bottle 
def rotational_bottle(center_of_mass,
                              mass_bottle=BOTTLE_MASS,
                              radius=BOTTLE_RADIUS,
                              length=BOTTLE_HEIGHT):
    j_cm = mass_bottle * (radius ** 2 / 2.0 + length ** 2 / 12.0)
    d = length / 2.0 - center_of_mass
    return j_cm + mass_bottle * d ** 2

# Center of mass of bottle + water (eps = m_bottle / (m_bottle + m_water))
def find_center_of_mass(positions, eps, length_bottle=BOTTLE_HEIGHT):
    positions = np.asarray(positions, dtype=float)
    z_bottle = length_bottle / 2.0
    z_water = np.mean(positions)
    return eps * z_bottle + (1.0 - eps) * z_water

# Time integration, Velocity–Verlet
def update_slice_positions(positions, center_of_mass, angular_velocity,
                           angle, velocities, alpha, mass, dt, gravity):
    positions = np.asarray(positions, dtype=float)
    velocities = np.asarray(velocities, dtype=float)

    drag = alpha / mass
    cos_theta = np.cos(angle)
    omega2 = angular_velocity ** 2
                             
    acc = -drag * velocities - omega2 * (positions - center_of_mass) - gravity * cos_theta
    v_half = velocities + 0.5 * acc * dt
    new_positions = positions + v_half * dt
    acc_new = -drag * v_half - omega2 * (new_positions - center_of_mass) - gravity * cos_theta
    new_velocities = v_half + 0.5 * acc_new * dt

    return new_positions, new_velocities

# Bounce on boundaries with restitution
def check_boundary_conditions(positions, previous_positions, velocities,
                              restitution, l_min, l_max):
    positions = np.asarray(positions, dtype=float)
    previous_positions = np.asarray(previous_positions, dtype=float)
    velocities = np.asarray(velocities, dtype=float)

    # Lower boundary
    below = positions < l_min
    if np.any(below):
        # Reflect inside and reverse direction
        positions[below] = l_min + (l_min - positions[below])
        velocities[below] = -restitution * velocities[below]

    # Upper boundary
    above = positions > l_max
    if np.any(above):
        positions[above] = l_max - (positions[above] - l_max)
        velocities[above] = -restitution * velocities[above]

    return positions, velocities
