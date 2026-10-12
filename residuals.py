import numpy as np
import matplotlib.pyplot as plt

import constants
import data_format


def baryonic_velocity(galaxy):
    gas_velocity = np.asarray(data_format.data_format(galaxy, "Vgas"), dtype=float)
    disk_velocity = np.asarray(data_format.data_format(galaxy, "Vdisk"), dtype=float)
    bulge_velocity = np.asarray(data_format.data_format(galaxy, "Vbul"), dtype=float)

    baryonic_velocity_squared = (
        gas_velocity**2
        + constants.stellar_mass_to_light_disk * disk_velocity**2
        + constants.stellar_mass_to_light_bulge * bulge_velocity**2
    )

    return np.sqrt(np.maximum(baryonic_velocity_squared, 0.0))


def nfw_total_velocity(galaxy):
    radius_kpc = np.asarray(data_format.data_format(galaxy, "Radius"), dtype=float)
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

    halo_velocity_squared = (
        constants.G * enclosed_halo_mass / radius_m / 1e6
    )

    baryonic = baryonic_velocity(galaxy)
    baryonic_squared = np.square(baryonic)
    total_velocity_squared = np.maximum(
        baryonic_squared + halo_velocity_squared,
        0.0,
    )
    return np.sqrt(total_velocity_squared)


def mond_velocity(galaxy):
    radius_kpc = np.asarray(data_format.data_format(galaxy, "Radius"), dtype=float)
    radius_m = radius_kpc * constants.kpc

    gas_velocity = np.asarray(data_format.data_format(galaxy, "Vgas"), dtype=float)
    disk_velocity = np.asarray(data_format.data_format(galaxy, "Vdisk"), dtype=float)
    bulge_velocity = np.asarray(data_format.data_format(galaxy, "Vbul"), dtype=float)

    baryonic_velocity_squared = (
        gas_velocity**2
        + constants.stellar_mass_to_light_disk * disk_velocity**2
        + constants.stellar_mass_to_light_bulge * bulge_velocity**2
    )

    newtonian_acceleration = (
        np.maximum(baryonic_velocity_squared, 0.0) * 1e6 / radius_m
    )

    mond_acceleration = 0.5 * (
        newtonian_acceleration
        + np.sqrt(
            newtonian_acceleration**2
            + 4 * newtonian_acceleration * constants.mond_acceleration_scale
        )
    )

    return np.sqrt(mond_acceleration * radius_m) / 1000


def plot_residuals(galaxy="NGC3198"):
    radius = np.asarray(data_format.data_format(galaxy, "Radius"), dtype=float)
    observed = np.asarray(data_format.data_format(galaxy, "Vobs"), dtype=float)
    errors = np.asarray(data_format.data_format(galaxy, "errV"), dtype=float)

    model_velocities = {
        "Baryons only": baryonic_velocity(galaxy),
        "Baryons + NFW": nfw_total_velocity(galaxy),
        "MOND": mond_velocity(galaxy),
    }

    fig, ax = plt.subplots(figsize=(10, 5))

    for label, predicted in model_velocities.items():
        predicted = np.asarray(predicted, dtype=float)

        if not (
            len(radius) == len(observed) == len(errors) == len(predicted)
        ):
            raise ValueError(f"Array lengths do not match for {label}")

        residual = observed - predicted
        ax.plot(
            radius,
            residual,
            marker="o",
            markersize=3,
            linewidth=1.2,
            label=label,
        )

    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Radius (kpc)")
    ax.set_ylabel(r"$V_{\mathrm{obs}}-V_{\mathrm{model}}$ (km/s)")
    ax.set_title(f"Rotation-curve residuals: {galaxy}")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_residuals()