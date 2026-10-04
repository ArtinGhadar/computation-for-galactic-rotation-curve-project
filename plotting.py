import matplotlib.pyplot as plt
import numpy as np

import constants
import data_format


def plot_curve_observed(galaxy):
    radius = data_format.data_format(galaxy, "Radius")
    observed_velocity = data_format.data_format(galaxy, "Vobs")
    velocity_error = data_format.data_format(galaxy, "errV")

    plt.errorbar(
        radius,
        observed_velocity,
        yerr=velocity_error,
        fmt="o",
        capsize=3,
        label=galaxy,
    )


def plot_curve_constant_density(galaxy):
    radius_kpc = np.array(data_format.data_format(galaxy, "Radius"))
    radius_m = radius_kpc * constants.kpc
    velocity_m_s = radius_m * np.sqrt(
        4 * np.pi * constants.G * constants.rho_0 / 3
    )
    velocity_km_s = velocity_m_s / 1000
    plt.plot(
        radius_kpc,
        velocity_km_s,
        linestyle=":",
        label=f"{galaxy} (constant density)",
    )


def plot_curve_nfw(galaxy):
    radius_kpc = np.array(data_format.data_format(galaxy, "Radius"))
    radius_m = radius_kpc * constants.kpc
    radius_ratio = radius_m / constants.halo_scale_radius
    enclosed_mass = 4 * np.pi * constants.halo_scale_density * (
        constants.halo_scale_radius**3
    ) * (
        np.log1p(radius_ratio) - radius_ratio / (1 + radius_ratio)
    )
    velocity_m_s = np.sqrt(constants.G * enclosed_mass / radius_m)
    velocity_km_s = velocity_m_s / 1000
    plt.plot(
        radius_kpc,
        velocity_km_s,
        linestyle="--",
        label=f"{galaxy} (NFW halo)",
    )


def plot_curve_pseudo_isothermal(galaxy):
    radius_kpc = np.array(data_format.data_format(galaxy, "Radius"))
    radius_m = radius_kpc * constants.kpc
    core_radius = constants.pseudo_iso_core_radius
    radius_ratio = radius_m / core_radius

    velocity_squared = (
        4
        * np.pi
        * constants.G
        * constants.pseudo_iso_central_density
        * core_radius**2
        * (1 - np.arctan(radius_ratio) / radius_ratio)
    )
    velocity_km_s = np.sqrt(velocity_squared) / 1000

    plt.plot(
        radius_kpc,
        velocity_km_s,
        linestyle="-.",
        label=(
            f"{galaxy} (pseudo-isothermal halo; "
            f"rho_c=0.1 Msun/pc^3, r_c=5 kpc)"
        ),
    )


def plot_curve_mass_distribution(galaxy):
    radius_kpc = np.array(data_format.data_format(galaxy, "Radius"))
    gas_velocity = np.array(data_format.data_format(galaxy, "Vgas"))
    disk_velocity = np.array(data_format.data_format(galaxy, "Vdisk"))
    bulge_velocity = np.array(data_format.data_format(galaxy, "Vbul"))

    baryonic_velocity_squared = (
        gas_velocity * np.abs(gas_velocity)
        + constants.stellar_mass_to_light_disk * disk_velocity**2
        + constants.stellar_mass_to_light_bulge * bulge_velocity**2
    )
    baryonic_velocity = np.sqrt(np.maximum(baryonic_velocity_squared, 0))

    radius_m = radius_kpc * constants.kpc
    radius_ratio = radius_m / constants.halo_scale_radius
    enclosed_halo_mass = 4 * np.pi * constants.halo_scale_density * (
        constants.halo_scale_radius**3
    ) * (np.log1p(radius_ratio) - radius_ratio / (1 + radius_ratio))
    halo_velocity_squared = constants.G * enclosed_halo_mass / radius_m / 1e6
    total_velocity = np.sqrt(
        np.maximum(baryonic_velocity_squared + halo_velocity_squared, 0)
    )

    plt.plot(
        radius_kpc,
        baryonic_velocity,
        linestyle=":",
        label=f"{galaxy} (baryons; disk M/L=0.5, bulge M/L=0.7)",
    )
    plt.plot(
        radius_kpc,
        total_velocity,
        linestyle="--",
        label=f"{galaxy} (baryons + fixed NFW halo)",
    )


def plot_curve_mond(galaxy):
    radius_kpc = np.array(data_format.data_format(galaxy, "Radius"))
    radius_m = radius_kpc * constants.kpc
    gas_velocity = np.array(data_format.data_format(galaxy, "Vgas"))
    disk_velocity = np.array(data_format.data_format(galaxy, "Vdisk"))
    bulge_velocity = np.array(data_format.data_format(galaxy, "Vbul"))

    baryonic_velocity_squared = (
        gas_velocity * np.abs(gas_velocity)
        + constants.stellar_mass_to_light_disk * disk_velocity**2
        + constants.stellar_mass_to_light_bulge * bulge_velocity**2
    )
    newtonian_acceleration = (
        np.maximum(baryonic_velocity_squared, 0) * 1e6 / radius_m
    )
    mond_acceleration = 0.5 * (
        newtonian_acceleration
        + np.sqrt(
            newtonian_acceleration**2
            + 4 * newtonian_acceleration * constants.mond_acceleration_scale
        )
    )
    mond_velocity = np.sqrt(mond_acceleration * radius_m) / 1000

    plt.plot(
        radius_kpc,
        mond_velocity,
        linestyle="--",
        label=f"{galaxy} (MOND; a0={constants.mond_acceleration_scale:.1e} m/s2)",
    )


galaxy = "NGC3198"
plot_curve_observed(galaxy)
plot_curve_mass_distribution(galaxy)
# plot_curve_constant_density(galaxy)
plot_curve_mond(galaxy)
plot_curve_nfw(galaxy)
# plot_curve_pseudo_isothermal(galaxy)



plt.title("Galaxy rotation curves")
plt.xlabel("Radius (kpc)")
plt.ylabel("Velocity (km/s)")
plt.grid(True, alpha=0.3)
plt.legend(title="Galaxy")
plt.tight_layout()
plt.show()


