import numpy as np
import matplotlib.pyplot as plt

from constants import *
from kinetics import (
    rotational_water,
    rotational_bottle,
    find_center_of_mass,
    update_slice_positions,
    check_boundary_conditions,
)


def simulate_flip(
    filling_fraction=0.4,
    n_slices=50,
    t_max=1.2,
    dt=1e-3,
    omega0=30.0,
    theta0=np.deg2rad(75.0),
    restitution=0.2,
    alpha=0.02,
):
    water_mass = filling_fraction * WATER_MASS_MAX
    bottle_mass = BOTTLE_MASS
    total_mass = bottle_mass + water_mass
    eps = bottle_mass / total_mass
    slice_mass = water_mass / n_slices

    z0 = np.linspace(0.0, filling_fraction * BOTTLE_HEIGHT, n_slices)
    positions = z0.copy()
    velocities = np.zeros_like(positions)

    angle = theta0
    omega = omega0

    n_steps = int(t_max / dt) + 1
    times = np.linspace(0.0, t_max, n_steps)
    angles = np.zeros(n_steps)
    omegas = np.zeros(n_steps)
    coms = np.zeros(n_steps)

    center_of_mass = find_com(positions, eps)
    J_bottle = rotational_bottle(center_of_mass)
    J_water = rotational_water(positions, center_of_mass, water_mass, n_slices)
    J_tot = J_bottle + J_water
    L = J_tot * omega  # conserved angular momentum

    for i, _ in enumerate(times):
        center_of_mass = find_com(positions, eps)
        J_bottle = rotational_bottle(center_of_mass)
        J_water = rotational_water(positions, center_of_mass, water_mass, n_slices)
        J_tot = J_bottle + J_water

        omega = L / J_tot

        prev_positions = positions.copy()
        positions, velocities = update_slice_positions(
            positions,
            center_of_mass,
            omega,
            angle,
            velocities,
            alpha,
            slice_mass,
            dt,
            G,
        )

        positions, velocities = check_boundary_conditions(
            positions,
            prev_positions,
            velocities,
            restitution,
            l_min=0.0,
            l_max=BOTTLE_HEIGHT,
        )

        angle += omega * dt

        angles[i] = angle
        omegas[i] = omega
        coms[i] = center_of_mass

    return times, angles, omegas, coms


if __name__ == "__main__":
    t, theta, omega, com = simulate_flip()

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(t, omega)
    plt.ylabel(r'$\omega$ (rad/s)')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(t, com)
    plt.xlabel('t (s)')
    plt.ylabel('COM (m)')
    plt.grid(True)

    plt.tight_layout()
    plt.show()
