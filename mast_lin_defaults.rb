######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2012-12-11 16:21:07 +0000                                           #
#                                                                    #
######################################################################

@defaults_file_description = "MAST shot 27268/27274 (only difference is position of BES). This defaults file is for the BES further out => more turbulence."

######################################
# Defaults for namelist parameters
#######################################

@beta = 0.0047310768    # Ratio of particle to magnetic pressure (reference Beta, not total beta):  beta=n_0 T_0 /( B^2 / (8 pi))
@zeff = 1.5899834    # Effective ionic charge.

######################################
# Defaults for namelist kt_grids_knobs
#######################################

@grid_option = "range"    # The general layout of the perpendicular grid.

######################################
# Defaults for namelist kt_grids_box_parameters
#######################################

@nx = 256    # The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.
@ny = 96    # The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.
@jtwist = 160    # L_x = L_y  jtwist / (2 pi shat)
@y0 = 10.0    # The length of the box in the y direction (measured in the Larmour radius of species 1)
@x0 = 10.0    # The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)

######################################
# Defaults for namelist kt_grids_range_parameters
#######################################

@theta0_min = 0.0      # Set this instead of akx. Lower limit of kx range (usually 0 since usually fastest growing).
@theta0_max = 0.0      # Upper limit of kx range.
@aky_min = 0.0          # Lower limit of ky range
@aky_max = 2.0          # Upper limit of ky range

######################################
# Defaults for namelist kt_grids_single_parameters
#######################################

@theta0 = 0.0   # kx value (usually 0 since usually fastest growing)
@aky = 0.05     # ky value for reference species

######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 24    # Number of points along field line (theta) per 2 pi segment
@nperiod = 2    # Number of 2 pi segments along equilibrium magnetic field.
@shat =  3.9955695   # 
@rhoc = 0.79664360    # Flux surface label. Usually rho = diameter/diameter of LCFS
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
@equal_arc = ".false."    # 
@delrho = 1e-3

######################################
# Defaults for namelist le_grids_knobs
#######################################

@ngauss = 16    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 16    # Total number of energy grid points

######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@apfac = 1.0    # 
@driftknob = 1.0    # 
@g_exb = 0.15913171	# Velocity shear	
@g_exb_start_timestep = 1000    # Flow shear is switched on at this time step.
@boundary_option = "default"

######################################
# Defaults for namelist fields_knobs
#######################################

@field_option = "local"    # Controls which time-advance algorithm is used for the linear terms.
@field_subgath = ".true."     #Set to true to use allgatherv to fetch parts of the field update vector calculated on other procs. When false uses a sum_allreduce instead. This doesn't rely on sub-communicators so should work for any layout and processor count.
@do_smart_update = ".true." 
@field_local_allreduce = ".true." 
@field_local_allreduce_sub = ".true." 

######################################
# Defaults for namelist knobs
#######################################

@fphi = 1.0    # Multiplies Phi (electrostatic potential).
@fapar = 0.0    # Multiplies A_par. Use 1 for finite beta (electromagnetic), 0 otherwise (electrostatic)
@faperp = 0.0    # Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead
@delt = 0.01    # Time step
@delt_minimum = 1.0e-6	#Minimum time step allowed
@nstep = 3000    # Maximum number of timesteps
@abort_rapid_time_step_change = ".false." #Set whether error occurs for time step changing too rapidly. False better since will quit only of too small
@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted 

######################################
# Defaults for namelist layouts_knobs
#######################################

@layout = "lexys"    # 'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.
@unbalanced_xxf = ".true." # This allows GS2 to set up an unbalanced xxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_xxf = 0.5 # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@unbalanced_yxf = ".true." # This allows GS2 to set up an unbalanced yxxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_yxf = 0.5 # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@opt_redist_nbk = ".true." # This enables the use of non-blocking communication in redistribute routines.
@opt_redist_init = ".true." # This enables optimized initialization routines for creating redistribution objects.
@opt_redist_persist = ".true." 
@opt_redist_persist_overlap = ".true." 
@intmom_sub = ".true." # This enables use of sub-communicators to do reduction associated with calculation of moments of distribution function. Most advantageous for collisional runs without LE layouts.
@intspec_sub = ".true." # This enables use of sub-communicators to do reduction associated with calculation of species integrated moments of distribution function.

######################################
# Defaults for namelist collisions_knobs
#######################################

@collision_model = "default"    # Collision model used in the simulation. Options: 'default', 'none', 'lorentz', 'ediffuse'
@use_le_layout = ".false."       # Use new LE layout for collisions

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
@fprim_1 = 2.642780    # -1/n (dn/drho)
@uprim_1 = 0.0    # ?
@vnewk_1 = 2.098588e-02    # collisionality parameter
@type_1 = "ion"    # Type of species, e.g. 'ion', 'electron', 'beam'

######################################
# Defaults for namelist dist_fn_species_knobs_1
#######################################

@fexpr_1 = 0.48    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_1 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05

######################################
# Defaults for namelist species_parameters_2
#######################################

@z_2 = -1.0    # Charge
@mass_2 = 1.0/(2.0*1836.0)   # Mass
@dens_2 = 1.0    # Density	
@temp_2 = 1.091722    # Temperature
@tprim_2 = 5.773614    # -1/T (dT/drho)
@fprim_2 = 2.642780    # -1/n (dn/drho)
@uprim_2 = 0.0    # ?
@vnewk_2 = 5.900574e-01    # collisionality parameter
@type_2 = "electron"    # Type of species, e.g. 'ion', 'electron', 'beam'

######################################
# Defaults for namelist dist_fn_species_knobs_2
#######################################

@fexpr_2 = 0.48    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_2 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05

######################################
# Defaults for namelist init_g_knobs
#######################################

@chop_side = ".false."    # Rarely needed. Forces asymmetry into initial condition.
@phiinit = 1.0e-03    # Average amplitude of initial perturbation of each Fourier mode.
@ginit_option = "noise"    # Sets the way that the distribution function is initialized.
@restart_dir = "nc"    # 

######################################
# Defaults for namelist gs2_diagnostics_knobs
#######################################

@print_line = ".false."    # Estimated frequencies and growth rates to the screen/stdout
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
@write_phi_over_time = ".true."
@dump_check1 = ".false."    # 
@write_nl_flux = ".true."	# Write nonlinear fluxes
@write_verr = ".true." # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_cross_phase = ".true."
@write_full_moments_notgc = ".true."

######################################
# Defaults for namelist hyper_knobs
#######################################
@const_amp = ".true."

