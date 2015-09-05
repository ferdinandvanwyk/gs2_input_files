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

@nx = 256    # The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.
@ny = 96    # The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.
@jtwist = 600    # L_x = L_y  jtwist / (2 pi shat)
@y0 = 2    # The length of the box in the y direction (measured in the Larmour radius of species 1)
@x0 = 10.0    # The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)
@naky = 2    # The actual number of ky modes (do not use for nonlinear runs, use ny)
@ntheta0 = 16    # 

######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 24   # Number of points along field line (theta) per 2 pi segment
@nperiod = 2    # Number of 2 pi segments along equilibrium magnetic field.
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

@ngauss = 10    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 32    # Total number of energy grid points


######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@apfac = 1.0    # 
@driftknob = 1.0    # 
@g_exb = 0.0	# Velocity shear	
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
@nstep = 5000    # Maximum number of timesteps
@abort_rapid_time_step_change = ".false." #Set whether error occurs for time step changing too rapidly. False better since will quit only of too small
@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted 


######################################
# Defaults for namelist reinit_knobs
#######################################

@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted 
@delt_minimum = 1.0e-06    # The minimum time step is delt_minimum.


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

@nspec = 1    # Number of kinetic species evolved.


######################################
# Defaults for namelist species_parameters_1
#######################################

@z_1 = 1.0    # Charge
@mass_1 = 1.0    # Mass
@dens_1 = 1.0    # Density	
@temp_1 = 1.0    # Temperature
@tprim_1 = 6.92    # -1/T (dT/drho)
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
@tprim_2 = 6.92    # -1/T (dT/drho)
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
@write_phi_over_time = ".false."
@dump_check1 = ".false."    # 
@write_nl_flux = ".true."	# Write nonlinear fluxes
@write_verr = ".true." # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_cross_phase = ".true."
@write_full_moments_notgc = ".true."

######################################
# Defaults for namelist hyper_knobs
#######################################
@const_amp = ".true."

