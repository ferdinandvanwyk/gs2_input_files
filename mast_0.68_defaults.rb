######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2015-12-05 16:28:13 +0000                                           #
#                                                                    #
######################################################################

@defaults_file_description = ""


######################################
# Defaults for namelist parameters
#######################################

@beta = 0.005578844575006868    # Ratio of particle to magnetic pressure (reference Beta, not total beta):  beta=n_0 T_0 /( B^2 / (8 pi))
@zeff = 1.5651515769703725    # Effective ionic charge.


######################################
# Defaults for namelist kt_grids_knobs
#######################################

@grid_option = "range"    # The general layout of the perpendicular grid.


######################################
# Defaults for namelist kt_grids_box_parameters
#######################################

@naky = 16    # The actual number of ky modes (do not use for nonlinear runs, use ny)
@ntheta0 = 1    # 


######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 24    # Number of points along field line (theta) per 2 pi segment
@nperiod = 2    # Number of 2 pi segments along equilibrium magnetic field.
@shat = 3.7415049076080322    # 
@rhoc = 0.7811812162399292    # Flux surface label. Usually rho = diameter/diameter of LCFS
@qinp = 2.1158999863835426    # Sets value of the safety factor when using local toroidal equilibrium model.
@akappa = 1.4524037731998753    # Sets local elongation when local toroidal equilibrium is specified.
@akappri = 0.42305052280426025    # akappri = dkappa/drho
@tri = 0.20206886094940368    # tri = arcsin[(R(max(Z))-R_major)/r_mid]
@tripri = 0.47173187136650085    # tripri =  dtri/drho
@shift = -0.29601976607114266    # Sets Shafranov shift. See online help for definition.
@rmaj = 1.4937691051507538    # Major radius/a (Position of magnetic axis)
@r_geo = 1.4091217861022285    # Major radius/a (centerpoint of LCFS)


######################################
# Defaults for namelist theta_grid_knobs
#######################################

@equilibrium_option = "eik"    # Controls which geometric assumptions are used in the run.


######################################
# Defaults for namelist theta_grid_eik_knobs
#######################################

@itor = 1    # 
@iflux = 0    # 
@irho = 2    # Chooses definition of flux surface coordinate.
@local_eq = ".true."    # 
@bishop = 4    # 
@s_hat_input = 3.7415049076080322    # 
@beta_prime_input = -0.1962238997220993    # The gradient of the pressure.
@delrho = 0.001    # 
@equal_arc = ".false."    # 


######################################
# Defaults for namelist le_grids_knobs
#######################################

@ngauss = 16    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 16    # Total number of energy grid points


######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@boundary_option = "default"    # Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi).
@g_exb = 0.14179531982979182    # 
@apfac = 1.0    # 
@driftknob = 1.0    # 


######################################
# Defaults for namelist fields_knobs
#######################################

@field_option = "local"    # Controls which time-advance algorithm is used for the linear terms.
@field_subgath = ".true."    # Set to TRUE to use allgatherv to fetch part of the field update calculated on other procs. FALSE uses a sum_allreduce instead.
@do_smart_update = ".true."    # 
@field_local_allreduce = ".true."    # 
@field_local_allreduce_sub = ".true."    # 


######################################
# Defaults for namelist knobs
#######################################

@fphi = 1.0    # Multiplies Phi (electrostatic potential).
@fapar = 0.0    # Multiplies A_par. Use 1 for finite beta (electromagnetic), 0 otherwise (electrostatic)
@faperp = 0.0    # Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead
@delt = 0.01    # Time step
@nstep = 10000    # Maximum number of timesteps
@avail_cpu_time = 600    # Specify the available wall clock time in seconds. GS2 will exit before this time.


######################################
# Defaults for namelist reinit_knobs
#######################################

@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted
@abort_rapid_time_step_change = ".false."    # If true (default), exit if time step changes rapidly.


######################################
# Defaults for namelist layouts_knobs
#######################################

@layout = "lexys"    # 'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.
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


######################################
# Defaults for namelist collisions_knobs
#######################################

@collision_model = "default"    # Collision model used in the simulation. Options: 'default', 'none', 'lorentz', 'ediffuse'
@use_le_layout = ".false."    # Use more efficient layouts for le grid


######################################
# Defaults for namelist hyper_knobs
#######################################

@const_amp = ".true."    # Detrmines whether damping rate depends on amplitude variations. Recommend FALSE for nonlinear, TRUE for linear.


######################################
# Defaults for namelist nonlinear_terms_knobs
#######################################

@nonlinear_mode = "off"    # Include nonlinear terms? ('on','off')
@flow_mode = "off"    # 
@cfl = 0.5    # The maximum delt < cfl * min(Delta_perp/v_perp)


######################################
# Defaults for namelist species_knobs
#######################################

@nspec = 2    # Number of kinetic species evolved.


######################################
# Defaults for namelist species_parameters_1
#######################################

@z_1 = 1    # Charge
@mass_1 = 1.0    # Mass
@dens_1 = 1.0    # Density	
@temp_1 = 1.0    # Temperature
@tprim_1 = 5.3815484046936035    # -1/T (dT/drho)
@fprim_1 = 2.4859049320220947    # -1/n (dn/drho)
@uprim_1 = 0.0    # ?
@vnewk_1 = 0.011119120051652425    # collisionality parameter
@type_1 = "ion"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist species_parameters_2
#######################################

@z_2 = -1    # Charge
@mass_2 = 0.0002723311546840959    # Mass
@dens_2 = 1.1368033814881884    # Density	
@temp_2 = 1.095662143687228    # Temperature
@tprim_2 = 5.594787120819092    # -1/T (dT/drho)
@fprim_2 = 2.2827329635620117    # -1/n (dn/drho)
@uprim_2 = 0.0    # ?
@vnewk_2 = 0.5880413929789964    # collisionality parameter
@type_2 = "electron"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist dist_fn_species_knobs_1
#######################################

@fexpr_1 = 0.48    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_1 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05


######################################
# Defaults for namelist dist_fn_species_knobs_2
#######################################

@fexpr_2 = 0.48    # Temporal implicitness parameter. Recommended value: 0.48
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

@write_nl_flux = ".true."    # Write nonlinear fluxes as a function of time.
@print_line = ".false."    # Estimated frequencies and growth rates to the screen/stdout
@write_verr = ".true."    # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_line = ".true."    # If (write_ascii = T) write estimated frequencies and growth rates to the output file
@write_omega = ".true."    # If (write_ascii = T) instantaneous omega to output file. Very heavy output
@write_omavg = ".true."    # If (write_ascii = T) time-averaged growth rate and frequency to the output file.
@write_eigenfunc = ".true."    # If (write_ascii = T) Normalized phi written to runname.eigenfunc
@write_final_fields = ".true."    # If (write_ascii = T) Phi(theta) written to '.fields'
@write_final_moments = ".true."    # write final n, T
@nsave = 500    # Write restart files every nsave timesteps
@nwrite = 10    # Output diagnostic data every nwrite
@navg = 10    # Any time averages performed over navg
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
@omegatinst = 500.0    # Recommended value: 500.
@save_for_restart = ".false."    # Write restart files.
@write_full_moments_notgc = ".true."    # 
@write_cross_phase = ".true."    # 
@dump_check1 = ".false."    # 
@write_phi_over_time = ".true."    # Write entire Phi field to NetCDF file every nwrite.


######################################
# Defaults for namelist kt_grids_range_parameters
#######################################

@naky = 32    # The number of 'actual' ky modes.
@ntheta0 = 1    # 
@aky_min = 0.0    # 
@aky_max = 3.0    # 
@theta0_min = 0.0    # 
@theta0_max = 0.0    # 


######################################
# Defaults for namelist diagnostics_config
#######################################

@write_phi_over_time = ".false."    # Write entire phi field to NetCDF file every nwrite.
@write_omega = ".true."    # Write growth rates and frequencies to the netcdf file
@navg = 10    # Any time averages performed over navg
@omegatinst = 500.0    # Growth rates > omegatinst assumed numerical instability.
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
@print_line = ".false."    # Estimated frequencies and growth rates to the screen/stdout
@write_line = ".true."    # Write estimated frequencies and growth rates to the output file
@write_full_moments_notgc = ".true."    # 
@write_verr = ".true."    # Write velocity space diagnostics
@write_cross_phase = ".true."    # Write cross phase between electron temperature and density.
@write_eigenfunc = ".true."    # If (write_ascii = T) Normalized phi written to runname.eigenfunc
@write_final_fields = ".true."    # If (write_ascii = T) Phi(theta) written to '.fields'
@write_final_moments = ".true."    # write final n, T
