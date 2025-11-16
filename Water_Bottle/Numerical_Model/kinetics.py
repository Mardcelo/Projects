import numpy as np
from constants import *


def rotational_water(positions, center_of_mass, water_mass, n_slices):
    positions = np.asarray(positions, float)
    return water_mass / n_slices * np.sum((positions - center_of_mass) ** 2)


def rotational_bottle(center_of_mass,
                              mass_bottle=BOTTLE_MASS,
                              radius=BOTTLE_RADIUS,
                              length=BOTTLE_HEIGHT):
    j_bottle_axis = mass_bottle * (radius ** 2 / 2.0 + length ** 2 / 12.0)
    return j_bottle_axis + mass_bottle * (length / 2.0 - center_of_mass) ** 2

# Find Center Of Mass (COM) - I made it lol 
def find_com(positions, eps, length_bottle=BOTTLE_HEIGHT):
    positions = np.asarray(positions, float)
    return eps * length_bottle / 2.0 + (1.0 - eps) * np.mean(positions)


def update_slice_positions(positions, center_of_mass, angular_velocity,
                           angle, velocities, alpha, mass, dt, gravity):
    positions = np.asarray(positions, float)
    velocities = np.asarray(velocities, float)

    if alpha == 0.0:
        # No friction
        acc = -angular_velocity**2 * (positions - center_of_mass) - gravity * np.cos(angle)
        v_half = velocities + 0.5 * acc * dt
        new_positions = positions + v_half * dt
        acc_new = -angular_velocity**2 * (new_positions - center_of_mass) - gravity * np.cos(angle)
        new_velocities = v_half + 0.5 * acc_new * dt
        return new_positions, new_velocities

    tau = (2.0 * mass) / alpha
    omega = float(angular_velocity)
    gcos = gravity * np.cos(angle)

    # \omega -> 0 limit of the analytic solution
    if abs(omega) < 1e-6:
        exp2 = np.exp(-2.0 * dt / tau)
        v0 = velocities
        z0 = positions
        v_new = -gcos * tau / 2.0 + (gcos * tau / 2.0 + v0) * exp2
        z_new = (
            z0
            + tau * v0 / 2.0 * (1.0 - exp2)
            - gcos * (0.5 * dt * tau + 0.25 * tau**2 * (1.0 - exp2))
        )
        return z_new, v_new

    factor_1 = positions - center_of_mass + (gcos / (omega ** 2))
    factor_2 = velocities / omega + (1.0 / (tau * omega)) * factor_1

    exp = np.exp(-dt / tau)
    cosh = np.cosh(omega * dt)
    sinh = np.sinh(omega * dt)

    new_positions = (
        center_of_mass
        - (gcos / (omega ** 2))
        + factor_1 * exp * cosh
        + factor_2 * exp * sinh
    )

    factor_1_v = ((omega ** 2 - 1.0 / (tau ** 2)) / omega) * factor_1 + (1.0 / (tau * omega)) * velocities
    factor_2_v = velocities
    new_velocities = factor_1_v * exp * sinh + factor_2_v * exp * cosh
    return new_positions, new_velocities


def check_boundary(positions, previous_positions, velocities,
                              restitution, l_min, l_max):
    positions = np.asarray(positions, float)
    previous_positions = np.asarray(previous_positions, float)
    velocities = np.asarray(velocities, float)

    below = positions < l_min
    if np.any(below):
        positions[below] = l_min + (l_min - positions[below])
        velocities[below] = -restitution * velocities[below]

    above = positions > l_max
    if np.any(above):
        positions[above] = l_max - (positions[above] - l_max)
        velocities[above] = -restitution * velocities[above]

    return positions, velocities
