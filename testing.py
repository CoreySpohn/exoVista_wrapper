"""
The goal of this file is to be the primary testing ground. Basic flow should be
1. Get relevant data
2. Create universe based on the saved data
3. Make plots of what is desired
"""
import pickle
from pathlib import Path

import astropy.constants as constants
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import rvsearch
from astropy.time import Time
from rvsearch import driver, inject, plots, search, utils
from tqdm import tqdm

import exoVista_wrapper.functions as ewf
from exoVista_wrapper.universe import System, Universe

# ewf.get_data(universes=[1])

# u1 = Universe("data/1/", cache=True)

system = System("./data/1/17-HIP_169-TYC_-mv_9.24-L_0.11-d_15.20-Teff_3953.67.fits")
for planet in system.planets:
    print(f"{planet.L_type} {planet.Rp_type}")
planet = system.planets[5]
time_ind = 3
# time = Time(planet.t[time_ind].value, format="decimalyear")
time = planet._t[time_ind]
mean_anom = planet.mean_anom(time)
x, y, z, M, vx, vy, vz, vs = [], [], [], [], [], [], [], []
# for t in planet.t:
# M.append(planet.mean_anom(t).to(u.rad).value)
# _x, _y, _z = planet.calc_position_vectors(t)
# x.append(_x.to(u.AU).value)
# y.append(_y.to(u.AU).value)
# z.append(_z.to(u.AU).value)
# _x, _y, _z = planet.calc_vectors(t, return_v=True, return_r=False)
# vx.append(_x.to(u.AU / u.yr).value)
# vy.append(_y.to(u.AU / u.yr).value)
# vz.append(_z.to(u.AU / u.yr).value)
# vs.append(planet.calc_vs(t).to(u.AU / u.yr).value[0])

# system_rv = system.get_rv(times)
# times = planet._t.decompose().value
times = (np.linspace(0, 10, int(365 * 10 / 2)) * u.yr).decompose().value
# system.propagate_system(times)

rv_error = 0.1 * u.m / u.s
rv_df = system.simulate_rv_observations(times, rv_error)
data_path = Path("./data/")
rv_curve_path = Path(data_path, "rv_curve").with_suffix(".csv")

IWA = 0.058 * u.arcsec
OWA = 6 * u.arcsec

min_period = (
    2
    * np.pi
    * np.sqrt((system.star.dist * np.tan(IWA)) ** 3 / (constants.G * system.star.mass))
).to(u.d)
max_period = (
    2
    * np.pi
    * np.sqrt((system.star.dist * np.tan(OWA)) ** 3 / (constants.G * system.star.mass))
).to(u.d)

searcher = search.Search(
    rv_df,
    min_per=min_period.value,
    oversampling=10,
    # max_per=min(max_period.value, 10000),
    workers=2,
    mcmc=True,
    verbose=True,
    max_planets=4,
    # mstar=(system.star.mass.to(u.M_sun).value, 0),
)
breakpoint()
searcher.run_search()

# plt.scatter(times, system.star._vz.decompose().value, label="exoVista")
# plt.scatter(times, system.star.vectors.vz, label="rebound")

# cmap = plt.get_cmap("viridis")
# cmap2 = plt.get_cmap("plasma")
# for i, time in enumerate(tqdm(times, desc="Plotting")):
#     fig, ax = plt.subplots()
#     ax = ewf.plot_planet_positions(ax, system, cmap, i)
# fig.savefig(f".tmp/{i:05}.png")
# plt.close()

# cvals = np.linspace(0, 1, len(system.planets))
# for j, planet in enumerate(system.planets):
#     # color = cmap2(cvals[j])
#     ms = ewf.planet_marker_size(
#         planet._z[j].decompose(),
#         planet._z.decompose(),
#         base_size=5 + planet.radius.to(u.R_earth).value,
#         factor=0.1,
#     )
#     ax.scatter(
#         planet._x[i].to(u.AU),
#         planet._y[i].to(u.AU),
#         label=f"Planet {j}",
#         color=planet.subtype_color,
#         s=ms,
#     )
#     ax.set_xlim([-8, 8])
#     ax.set_ylim([-8, 8])

rv_error = 0.1 * u.m / u.s
rv_df = system.simulate_rv_observations(times, rv_error)
data_path = Path("./data/")
rv_curve_path = Path(data_path, "rv_curve").with_suffix(".csv")
# params = ewf.gen_radvel_params(planet, rv_error)
# post = ewf.gen_posterior(rv_df, params, rv_error)
# chains = ewf.gen_chains(post)

# fig, [axs, posax] = plt.subplots(ncols=2, nrows=2)
# axs[0].plot(planet.t, system_rv.to(u.AU / u.yr), color="b")
# axs[0].plot(planet.t, system.star.vz, color="k")
# axs[0].set_ylabel("RV")
# axs[1].plot(planet.t, vy, color="b")
# axs[1].plot(planet.t, planet.vy, color="k")
# axs[1].set_ylabel("y")
# posax[0].plot(planet.t, vz, color="b")
# posax[0].plot(planet.t, planet.vz, color="k")
# posax[0].set_ylabel("z")
# posax[1].plot(planet.t, z, color="b")
# posax[1].plot(planet.t, planet.z, color="k")
# posax[1].set_ylabel("z")
breakpoint()
