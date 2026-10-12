import numpy as np

import data_format
from residuals import baryonic_velocity, mond_velocity, nfw_total_velocity


def calculate_chi_squared(galaxy="NGC3198"):
    observed = np.asarray(data_format.data_format(galaxy, "Vobs"), dtype=float)
    errors = np.asarray(data_format.data_format(galaxy, "errV"), dtype=float)

    if len(observed) == 0:
        raise ValueError(f"No data found for galaxy {galaxy}")
    if len(observed) != len(errors):
        raise ValueError("Observed velocities and uncertainties have different lengths")
    if np.any(errors <= 0):
        raise ValueError("All velocity uncertainties must be positive")

    model_velocities = {
        "Baryons only": baryonic_velocity(galaxy),
        "Baryons + NFW": nfw_total_velocity(galaxy),
        "MOND": mond_velocity(galaxy),
    }

    degrees_of_freedom = len(observed)
    results = {}

    for label, predicted in model_velocities.items():
        predicted = np.asarray(predicted, dtype=float)
        if len(predicted) != len(observed):
            raise ValueError(f"Array lengths do not match for {label}")

        chi_squared = float(np.sum(((observed - predicted) / errors) ** 2))
        reduced_chi_squared = chi_squared / degrees_of_freedom
        results[label] = (chi_squared, reduced_chi_squared)

    return results


if __name__ == "__main__":
    results = calculate_chi_squared()
    print(f"{'Model':<18} {'Chi-squared':>14} {'Reduced chi-squared':>21}")
    for label, (chi_squared, reduced_chi_squared) in results.items():
        print(f"{label:<18} {chi_squared:>14,.0f} {reduced_chi_squared:>21,.1f}")
