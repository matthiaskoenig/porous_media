sim_flux_path: Path = DATA_DIR / "simliva" / "iri_flux_study_0" / "vtk"
sim_277k_path: Path = DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "vtk"
sim_310k_path: Path = DATA_DIR / "simliva" / "006_T_310_15K_P0__0Pa_t_24h" / "vtk"

# interpolate_xdmf(
#     xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
#     xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
#     times_interpolate=np.linspace(0, 600 * 60, num=11)  # [s] (11 points in 600 min) # static image
# )
#
# interpolate_xdmf(
#     xdmf_in=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results.xdmf",
#     xdmf_out=DATA_DIR / "simliva" / "005_T_277_15K_P0__0Pa_t_24h" / "results_interpolated_11.xdmf",
#     times_interpolate=np.linspace(0, 600 * 60, num=11)
# )
#
# interpolate_xdmf(
#     xdmf_in=sim_flux_path.parent / "results.xdmf",
#     xdmf_out=sim_flux_path.parent / "results_interpolated.xdmf",
#     times_interpolate=np.linspace(0, 10000, num=200)
# )

nums = [10, 200]
for num in nums:
    interpolate_xdmf(
        xdmf_in=sim_flux_path.parent / "results.xdmf",
        xdmf_out=sim_flux_path.parent / f"results_interpolated_{num}.xdmf",
        times_interpolate=np.linspace(0, 10000, num=num),
    )
