######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2014-10-16 12:20:24 +0100                                           #
#                                                                    #
######################################################################

@defaults_file_description = ""


######################################
# Defaults for namelist parameters
#######################################

@beta = 0.0047310768    # Ratio of particle to magnetic pressure (reference Beta, not total beta):  beta=n_0 T_0 /( B^2 / (8 pi))
@zeff = 1.5899834    # Effective ionic charge.


######################################
# Defaults for namelist kt_grids_knobs
#######################################

@grid_option = "single"    # The general layout of the perpendicular grid.


######################################
# Defaults for namelist kt_grids_box_parameters
#######################################

@nx = 256    # The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.
@ny = 96    # The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.
@jtwist = 313    # L_x = L_y  jtwist / (2 pi shat)
@y0 = 0.4    # The length of the box in the y direction (measured in the Larmour radius of species 1)
@x0 = 10.0    # The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)
@naky = 2    # The actual number of ky modes (do not use for nonlinear runs, use ny)
@ntheta0 = 1024    # 


######################################
# Defaults for namelist kt_grids_single_parameters
#######################################

@aky = 0.5    # The actual value of ky rho
@theta0 = 0.0    # 
@akx = 0.14644527086125353    # 


######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 20    # Number of points along field line (theta) per 2 pi segment
@nperiod = 3    # Number of 2 pi segments along equilibrium magnetic field.
@shat = 3.9955695    # 
@rhoc = 0.7966436    # Flux surface label. Usually rho = diameter/diameter of LCFS
@qinp = 2.31493    # Sets value of the safety factor when using local toroidal equilibrium model.
@akappa = 1.45734    # Sets local elongation when local toroidal equilibrium is specified.
@akappri = 0.44877291    # akappri = dkappa/drho
@tri = 0.205939    # tri = arcsin[(R(max(Z))-R_major)/r_mid]
@tripri = 0.46296135    # tripri =  dtri/drho
@shift = -0.30708938    # shift = -R q**2 dbeta/drho (>0)
@rmaj = 1.4891066    # Major radius/a (Position of magnetic axis)
@r_geo = 1.4091217    # Major radius/a (centerpoint of LCFS)


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
@s_hat_input = 3.9955695    # 
@beta_prime_input = -0.1212404    # The gradient of the pressure.
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
@g_exb = 0.15913171    # 
@apfac = 1.0    # 
@driftknob = 1.0    # 
@g_exb_start_timestep = 1000    # Flow shear is switched on at this time step.
@g_exb_error_limit = 0.1    # With flow shear in single mode, constrains relative error in phi^2.


######################################
# Defaults for namelist fields_knobs
#######################################

@field_option = "implicit"    # Controls which time-advance algorithm is used for the linear terms.
@field_subgath = ".false."    # Set to TRUE to use allgatherv to fetch part of the field update calculated on other procs. FALSE uses a sum_allreduce instead.
@do_smart_update = ".false."    # 
@field_local_allreduce = ".false."    # 
@field_local_allreduce_sub = ".false."    # 


######################################
# Defaults for namelist knobs
#######################################

@fphi = 1.0    # Multiplies Phi (electrostatic potential).
@fapar = 0.0    # Multiplies A_par. Use 1 for finite beta (electromagnetic), 0 otherwise (electrostatic)
@faperp = 0.0    # Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead
@delt = 0.01    # Time step
@nstep = 5000    # Maximum number of timesteps
@avail_cpu_time = 1500    # Specify the available wall clock time in seconds. GS2 will exit before this time.


######################################
# Defaults for namelist reinit_knobs
#######################################

@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted
@delt_minimum = 1.0e-06    # The minimum time step is delt_minimum.
@abort_rapid_time_step_change = ".false."    # If true (default), exit if time step changes rapidly.


######################################
# Defaults for namelist layouts_knobs
#######################################

@layout = "lxyes"    # 'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.
@unbalanced_xxf = ".false."    # This allows GS2 to set up an unbalanced xxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_xxf = 0.5    # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@unbalanced_yxf = ".false."    # This allows GS2 to set up an unbalanced yxxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_yxf = 0.5    # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@opt_redist_nbk = ".false."    # This enables the use of non-blocking communication in redistribute routines.
@opt_redist_init = ".false."    # This enables optimized initialization routines for creating redistribution objects.
@intmom_sub = ".false."    # This enables use of sub-communicators to do reduction associated with calculation of moments of distribution function. Most advantageous for collisional runs without LE layouts.
@intspec_sub = ".false."    # This enables use of sub-communicators to do reduction associated with calculation of species integrated moments of distribution function.
@opt_redist_persist = ".false."    # Set to true to use persistent (non-blocking) comms in the redistribute routines.
@opt_redist_persist_overlap = ".false."    # Set to true to try to overlap the mpi and local parts of the gather/scatter routines.


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

@z_1 = 1.0    # Charge
@mass_1 = 1.0    # Mass
@dens_1 = 1.0    # Density	
@temp_1 = 1.0    # Temperature
@tprim_1 = 5.535046    # -1/T (dT/drho)
@fprim_1 = 2.64278    # -1/n (dn/drho)
@uprim_1 = 0.0    # ?
@vnewk_1 = 0.02098588    # collisionality parameter
@type_1 = "ion"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist species_parameters_2
#######################################

@z_2 = -1.0    # Charge
@mass_2 = 0.0002723311546840959    # Mass
@dens_2 = 1.0    # Density	
@temp_2 = 1.091722    # Temperature
@tprim_2 = 5.773614    # -1/T (dT/drho)
@fprim_2 = 2.64278    # -1/n (dn/drho)
@uprim_2 = 0.0    # ?
@vnewk_2 = 0.5900574    # collisionality parameter
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
@restart_file = "gs2.nc"    # Base of filenames with restart data.
@ginit_option = "noise"    # Sets the way that the distribution function is initialized.
@restart_dir = "v/id_25/nc"    # 


######################################
# Defaults for namelist gs2_diagnostics_knobs
#######################################

@write_nl_flux = ".true."    # Write nonlinear fluxes as a function of time.
@print_line = ".false."    # Estimated frequencies and growth rates to the screen/stdout
@write_verr = ".true."    # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_line = ".false."    # If (write_ascii = T) write estimated frequencies and growth rates to the output file
@write_omega = ".false."    # If (write_ascii = T) instantaneous omega to output file. Very heavy output
@write_omavg = ".false."    # If (write_ascii = T) time-averaged growth rate and frequency to the output file.
@write_eigenfunc = ".true."    # If (write_ascii = T) Normalized phi written to runname.eigenfunc
@write_final_fields = ".true."    # If (write_ascii = T) Phi(theta) written to '.fields'
@write_final_moments = ".true."    # write final n, T
@nsave = 500    # Write restart files every nsave timesteps
@nwrite = 10    # Output diagnostic data every nwrite
@navg = 10    # Any time averages performed over navg
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
@omegatinst = 500.0    # Recommended value: 500.
@save_for_restart = ".false."    # Write restart files.
@write_cross_phase = ".true."    # 
@dump_check1 = ".false."    # 
@write_phi_over_time = ".false."    # Write entire Phi field to NetCDF file every nwrite.
@write_full_moments_notgc = ".true."


######################################
# Defaults for namelist diagnostics_config
#######################################

@write_phi_over_time = ".false."    # Write entire phi field to NetCDF file every nwrite.
@write_omega = ".false."    # Write growth rates and frequencies to the netcdf file
@navg = 10    # Any time averages performed over navg
@omegatinst = 500.0    # Growth rates > omegatinst assumed numerical instability.
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
