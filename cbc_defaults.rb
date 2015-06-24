######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2012-12-18 15:10:23 +0000                                           #
#                                                                    #
######################################################################

@defaults_file_description = "This is the defaults file for the standard Cyclone Base Case (Dimits et. al.)"


######################################
# Defaults for namelist parameters
#######################################

@beta = 0.0    # Ratio of particle to magnetic pressure (reference Beta, not total beta):  beta=n_0 T_0 /( B^2 / (8 pi))
@zeff = 1.0    # Effective ionic charge.


######################################
# Defaults for namelist kt_grids_knobs
#######################################

@grid_option = "box"    # The general layout of the perpendicular grid.


######################################
# Defaults for namelist kt_grids_box_parameters
#######################################

@nx = 128    # The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.
@ny = 96    # The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.
@jtwist = 4    # L_x = L_y  jtwist / (2 pi shat)
@y0 = 10.0    # The length of the box in the y direction (measured in the Larmour radius of species 1)
@x0 = 10.0    # The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)


######################################
# Defaults for namelist kt_grids_range_parameters
#######################################

@theta0_min = 0.0	# Set this instead of akx. Lower limit of kx range (usually 0 since usually fastest growing).
@theta0_max = 0.0	# Upper limit of kx range.
@aky_min = 0.0		# Lower limit of ky range
@aky_max = 2.0		# Upper limit of ky range


######################################
# Defaults for namelist kt_grids_single_parameters
#######################################

@theta0 = 0.0	# kx value (usually 0 since usually fastest growing)
@aky = 0.05	# ky value for reference species

######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 20   # Number of points along field line (theta) per 2 pi segment
@nperiod = 1    # Number of 2 pi segments along equilibrium magnetic field.
@eps = 0.18    # eps=r/R
@epsl = 2.0    # epsl=2 a/R
@shat = 0.8    # 
@pk = 1.44    # pk = 2 a / q R
@shift = 0.0    # shift = -R q**2 dbeta/drho (>0)


######################################
# Defaults for namelist theta_grid_knobs
#######################################

@equilibrium_option = "s-alpha"    # Controls which geometric assumptions are used in the run.


######################################
# Defaults for namelist theta_grid_salpha_knobs
#######################################

@model_option = "default"    # 


######################################
# Defaults for namelist le_grids_knobs
#######################################

@ngauss = 8    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 16    # Total number of energy grid points


######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@omprimfac = 1.0    # 
@boundary_option = "linked"    # Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi).
@adiabatic_option = "iphi00=2"    # The form of the adiabatic response (if a species is being modeled as adiabatic).
@g_exb = 0.0    # 
@nonad_zero = ".true."


######################################
# Defaults for namelist fields_knobs
#######################################

@field_option = "local"    # Controls which time-advance algorithm is used for the linear terms.


######################################
# Defaults for namelist knobs
#######################################

@wstar_units = ".false."    # For linear runs only. Evolves each k_y with a different timestep.
@fphi = 1.0    # Multiplies Phi (electrostatic potential).
@fapar = 0.0    # Multiplies A_par. Use 1 for finite beta (electromagnetic), 0 otherwise (electrostatic)
@faperp = 0.0    # Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead
@delt = 0.01    # Time step
@nstep = 100000    # Maximum number of timesteps


######################################
# Defaults for namelist reinit_knobs
#######################################

@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted 
@delt_minimum = 1.0e-06    # The minimum time step is delt_minimum.


######################################
# Defaults for namelist layouts_knobs
#######################################

@layout = "xyles"    # 'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.


######################################
# Defaults for namelist collisions_knobs
#######################################

@collision_model = "default"    # Collision model used in the simulation. Options: 'default', 'none', 'lorentz', 'ediffuse'


######################################
# Defaults for namelist nonlinear_terms_knobs
#######################################

@nonlinear_mode = "on"    # Include nonlinear terms? ('on','off')
@flow_mode = "off"    # 
@cfl = 0.5    # The maximum delt < cfl * min(Delta_perp/v_perp)


######################################
# Defaults for namelist species_knobs
#######################################

@nspec = 1    # Number of kinetic species evolved.


######################################
# Defaults for namelist species_parameters_1
#######################################

@z_1 = 1.0    # Charge
@mass_1 = 1.0    # Mass
@dens_1 = 1.0    # Density	
@temp_1 = 1.0    # Temperature
@tprim_1 = 6.9    # -1/T (dT/drho)
@fprim_1 = 2.2    # -1/n (dn/drho)
@uprim_1 = 0.0    # ?
@vnewk_1 = 0.01    # collisionality parameter
@type_1 = "ion"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist dist_fn_species_knobs_1
#######################################

@fexpr_1 = 0.45    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_1 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05


######################################
# Defaults for namelist species_parameters_2
#######################################

@z_2 = -1.0    # Charge
@mass_2 = 1.0/(2.0*1836.0)   # Mass
@dens_2 = 1.0    # Density	
@temp_2 = 1.0    # Temperature
@tprim_2 = 6.9    # -1/T (dT/drho)
@fprim_2 = 2.2    # -1/n (dn/drho)
@uprim_2 = 0.0    # ?
@vnewk_2 = 0.6    # collisionality parameter
@type_2 = "electron"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist dist_fn_species_knobs_2
#######################################

@fexpr_2 = 0.45    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_2 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05

######################################
# Defaults for namelist init_g_knobs
#######################################

@chop_side = ".false."    # Rarely needed. Forces asymmetry into initial condition.
@phiinit = 0.001    # Average amplitude of initial perturbation of each Fourier mode.
@ginit_option = "noise"    # Sets the way that the distribution function is initialized.


######################################
# Defaults for namelist gs2_diagnostics_knobs
#######################################

@print_flux_line = ".false."    # Instantaneous fluxes output to screen
@write_nl_flux = ".true."    # Write nonlinear fluxes as a function of time.
@print_line = ".false."    # Estimated frequencies and growth rates to the screen/stdout
@write_line = ".false."    # If (write_ascii = T) write estimated frequencies and growth rates to the output file
@write_avg_moments = ".F."    # Write flux surface averaged low-order moments of g to runname.out.nc and runname.moments (if write_ascii = T)
@write_omega = ".false."    # If (write_ascii = T) instantaneous omega to output file. Very heavy output
@write_omavg = ".false."    # If (write_ascii = T) time-averaged growth rate and frequency to the output file.
@write_final_fields = ".true."    # If (write_ascii = T) Phi(theta) written to '.fields'
@write_final_moments = ".true."    # write final n, T
@nsave = 500    # Write restart files every nsave timesteps
@nwrite = 100    # Output diagnostic data every nwrite
@navg = 10    # Any time averages performed over navg
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
@omegatinst = 500.0    # Recommended value: 500.
@save_for_restart = ".true."    # Write restart files.
 @write_verr = ".true." # Write velocity space diagnostics to '.lpc' and '.verr' files
# write_g not specified --- Write the distribution function to the '.dist' (NetCDF?)
# write_gyx not specified --- Write dist fn at a given physical spacial point to a file
@write_hrate = ".false." # Write heating rate, collisonal entropy generation etc to '.heat'
# write_final_epar not specified --- If (write_ascii = T) E_parallel(theta) written to runname.eigenfunc
# write_lorentzian not specified --- Frequency Sweep Data
 @write_eigenfunc = ".true." # If (write_ascii = T) Normalized phi written to runname.eigenfunc
# write_parity not specified --- Writes parities in dist fn and particle fluxes
# write_flux_line not specified --- 
# write_ascii not specified --- 
# write_kpar not specified --- 
# write_gs not specified --- 
# write_gg not specified --- 
# write_lpoly not specified --- 
# write_fields not specified --- 
# write_final_antot not specified --- 
# write_cerr not specified --- 
# write_max_verr not specified --- 
# nmovie not specified --- 
# igomega not specified --- 
# exit_when_converged not specified --- 
# write_full_moments_notgc not specified --- 
# write_cross_phase not specified --- 
# dump_check1 not specified --- 
# dump_check2 not specified --- 
# dump_fields_periodically not specified --- 
# make_movie not specified --- 
 @write_phi_over_time = ".false." # Write entire Phi field to NetCDF file every nwrite.
# write_apar_over_time not specified --- Write entire A_parallel field to NetCDF file every nwrite.
# write_bpar_over_time not specified --- Write entire B_parallel field to NetCDF file every nwrite.
# write_symmetry not specified --- Test the symmetry properties of the GK eqn.
# save_distfn not specified --- Save dist_fn with lots of detail.
# write_correlation_extend not specified --- Extend domain of correlation function calculation.
# nwrite_mult not specified --- Large datasets written every nwrite_mult * nwrite timesteps.
# write_correlation not specified --- Write parallel correlation.
# write_moments not specified --- 
# write_final_db not specified --- Write final delta B.
@write_full_moments_notgc = ".true."
######################################
# Defaults for namelist diagnostics_config
#######################################

@nwrite_new = 100    # Diagnostic quantities are written every nwrite timesteps.
@write_omega = ".false."    # Write growth rates and frequencies to the netcdf file
@navg = 10    # Any time averages performed over navg
@omegatinst = 500.0    # Growth rates > omegatinst assumed numerical instability.
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol

#################
# Optimizations #
#################

@opt_init_bc = ".true."    # 
@opt_source = ".true."    # If true then use an optimised linear source calculation which uses pre-calculated coefficients.
@field_subgath = ".true."    # Set to TRUE to use allgatherv to fetch part of the field update calculated on other procs. FALSE uses a sum_allreduce instead.
@dump_response = ".true."    # 
@response_dir = "response"    # 
@do_smart_update = ".true."    # 
@field_local_allreduce = ".true."    # 
@field_local_allreduce_sub = ".true."    # 
@unbalanced_xxf = ".true."    # This allows GS2 to set up an unbalanced xxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_xxf = 0.5    # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@unbalanced_yxf = ".true."    # This allows GS2 to set up an unbalanced yxxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_yxf = 0.5    # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@opt_redist_nbk = ".true."    # This enables the use of non-blocking communication in redistribute routines.
@opt_redist_init = ".true."    # This enables optimized initialization routines for creating redistribution objects.
@intmom_sub = ".true."    # This enables use of sub-communicators to do reduction associated with calculation of moments of distribution function. Most advantageous for collisional runs without LE layouts.
@intspec_sub = ".true."    # This enables use of sub-communicators to do reduction associated with calculation of species integrated moments of distribution function.
@opt_redist_persist = ".true."    # Set to true to use persistent (non-blocking) comms in the redistribute routines.
@opt_redist_persist_overlap = ".true."    # Set to true to try to overlap the mpi and local parts of the gather/scatter routines.
@use_le_layout = ".false."       # Use new LE layout for collisions

