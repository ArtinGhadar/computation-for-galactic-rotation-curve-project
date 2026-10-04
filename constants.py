import math


G = 6.67430e-11
Msun = 1.98847e30
pc = 3.085677581e16
kpc = 1000 * pc
cm = 1e-2
GeV_c2 = 1.782661921e-27
rho_0 = 0.3 * GeV_c2 / cm**3

stellar_mass_to_light_disk = 0.5
stellar_mass_to_light_bulge = 0.7
mond_acceleration_scale = 1.2e-10
pseudo_iso_central_density = 0.1 * Msun / pc**3
pseudo_iso_core_radius = 5 * kpc

halo_mass = 1e12 * Msun
halo_radius = 200 * kpc
halo_concentration = 10
halo_scale_radius = halo_radius / halo_concentration
halo_scale_density = halo_mass / (
	4 * 3.141592653589793 * halo_scale_radius**3
	* (math.log(1 + halo_concentration)
	   - halo_concentration / (1 + halo_concentration))
)