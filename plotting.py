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
        label=f"{galaxy} (observed)",
    )


def plot_curve_nfw(galaxy):
    radius_kpc = np.array(data_format.data_format(galaxy, "Radius"))
    radius_m = radius_kpc * constants.kpc

    radius_ratio = radius_m / constants.halo_scale_radius

    enclosed_mass = (
        4
        * np.pi
        * constants.halo_scale_density
        * constants.halo_scale_radius**3
        * (
            np.log1p(radius_ratio)
            - radius_ratio / (1 + radius_ratio)
        )
    )

    velocity_m_s = np.sqrt(
        constants.G * enclosed_mass / radius_m
    )

    velocity_km_s = velocity_m_s / 1000

    plt.plot(
        radius_kpc,
        velocity_km_s,
        linestyle="--",
        label=f"{galaxy} (NFW halo contribution)",
    )



def plot_curve_mass_distribution(galaxy):
    radius_kpc = np.array(
        data_format.data_format(galaxy, "Radius")
    )

    gas_velocity = np.array(
        data_format.data_format(galaxy, "Vgas")
    )

    disk_velocity = np.array(
        data_format.data_format(galaxy, "Vdisk")
    )

    bulge_velocity = np.array(
        data_format.data_format(galaxy, "Vbul")
    )

    # ---------------------------------------------------------
    # BARYONIC CONTRIBUTION
    # ---------------------------------------------------------
    #
    # SPARC's contributions combine in quadrature:
    #
    # V_baryon^2 =
    # V_gas^2
    # + Upsilon_disk * V_disk^2
    # + Upsilon_bulge * V_bulge^2
    #
    baryonic_velocity_squared = (
        gas_velocity**2
        + constants.stellar_mass_to_light_disk
        * disk_velocity**2
        + constants.stellar_mass_to_light_bulge
        * bulge_velocity**2
    )

    baryonic_velocity = np.sqrt(
        np.maximum(baryonic_velocity_squared, 0)
    )

    # ---------------------------------------------------------
    # NFW HALO CONTRIBUTION
    # ---------------------------------------------------------

    radius_m = radius_kpc * constants.kpc
    radius_ratio = radius_m / constants.halo_scale_radius

    enclosed_halo_mass = (
        4
        * np.pi
        * constants.halo_scale_density
        * constants.halo_scale_radius**3
        * (
            np.log1p(radius_ratio)
            - radius_ratio / (1 + radius_ratio)
        )
    )

    # G*M/r gives m^2/s^2.
    # Divide by 1e6 to convert to (km/s)^2.
    halo_velocity_squared = (
        constants.G
        * enclosed_halo_mass
        / radius_m
        / 1e6
    )

    # ---------------------------------------------------------
    # TOTAL BARYONS + NFW
    # ---------------------------------------------------------

    total_velocity = np.sqrt(
        np.maximum(
            baryonic_velocity_squared
            + halo_velocity_squared,
            0
        )
    )

    # Plot baryons
    plt.plot(
        radius_kpc,
        baryonic_velocity,
        linestyle=":",
        label=(
            f"{galaxy} (baryons; "
            f"disk M/L={constants.stellar_mass_to_light_disk}, "
            f"bulge M/L={constants.stellar_mass_to_light_bulge})"
        ),
    )

    # Plot baryons + NFW
    plt.plot(
        radius_kpc,
        total_velocity,
        linestyle="--",
        label=f"{galaxy} (baryons + fixed NFW halo)",
    )


def plot_curve_mond(galaxy):
    radius_kpc = np.array(
        data_format.data_format(galaxy, "Radius")
    )

    radius_m = radius_kpc * constants.kpc

    gas_velocity = np.array(
        data_format.data_format(galaxy, "Vgas")
    )

    disk_velocity = np.array(
        data_format.data_format(galaxy, "Vdisk")
    )

    bulge_velocity = np.array(
        data_format.data_format(galaxy, "Vbul")
    )

    # ---------------------------------------------------------
    # NEWTONIAN BARYONIC ROTATION CURVE
    # ---------------------------------------------------------

    baryonic_velocity_squared = (
        gas_velocity**2
        + constants.stellar_mass_to_light_disk
        * disk_velocity**2
        + constants.stellar_mass_to_light_bulge
        * bulge_velocity**2
    )

    # Convert (km/s)^2 to (m/s)^2 and divide by radius
    # to obtain Newtonian acceleration.
    newtonian_acceleration = (
        np.maximum(baryonic_velocity_squared, 0)
        * 1e6
        / radius_m
    )

    # ---------------------------------------------------------
    # MOND
    # ---------------------------------------------------------

    # Simple MOND interpolation function:
    #
    # g = 1/2 * [g_N + sqrt(g_N^2 + 4*g_N*a_0)]
    #
    mond_acceleration = 0.5 * (
        newtonian_acceleration
        + np.sqrt(
            newtonian_acceleration**2
            + 4
            * newtonian_acceleration
            * constants.mond_acceleration_scale
        )
    )

    # V = sqrt(gR)
    mond_velocity = (
        np.sqrt(mond_acceleration * radius_m)
        / 1000
    )

    plt.plot(
        radius_kpc,
        mond_velocity,
        linestyle="--",
        label=(
            f"{galaxy} (MOND; "
            f"a0={constants.mond_acceleration_scale:.1e} m/s²)"
        ),
    )


# =============================================================
# PLOT
# =============================================================

galaxy = "NGC3198"

plot_curve_observed(galaxy)
plot_curve_mass_distribution(galaxy)
plot_curve_mond(galaxy) 
plot_curve_nfw(galaxy)


plt.title("Galaxy rotation curves")
plt.xlabel("Radius (kpc)")
plt.ylabel("Velocity (km/s)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()



