######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2014-11-22 22:15:55 +0000                                           #
#                                                                    #
######################################################################

@defaults_file_description = ""


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

@nx = 8    # The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.
@ny = 8
@jtwist = 1    # L_x = L_y  jtwist / (2 pi shat)
@y0 = 2.0    # The length of the box in the y direction (measured in the Larmour radius of species 1)
@x0 = 10.0    # The length of the box in the x direction (measured in the Larmour radius of species 1) if shat is 0 (ie 1e-6)


######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 8    # Number of points along field line (theta) per 2 pi segment
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

@ngauss = 5    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 8    # Total number of energy grid points


######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@omprimfac = 1.0    # 
@boundary_option = "linked"    # Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi).
@adiabatic_option = "iphi00=2"    # The form of the adiabatic response (if a species is being modeled as adiabatic).
@g_exb = 0.0    # 
@nonad_zero = ".true."    # If true switches on new parallel boundary condition where h=0 at incoming boundary instead of g=0.


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
@delt = 0.5    # Time step
@nstep = 50    # Maximum number of timesteps
@avail_cpu_time = 1200    # Specify the available wall clock time in seconds. GS2 will exit before this time.


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
@write_verr = ".true."    # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_line = ".false."    # If (write_ascii = T) write estimated frequencies and growth rates to the output file
@write_hrate = ".false."    # Write heating rate, collisonal entropy generation etc to '.heat'
@write_avg_moments = ".F."    # Write flux surface averaged low-order moments of g to runname.out.nc and runname.moments (if write_ascii = T)
@write_omega = ".false."    # If (write_ascii = T) instantaneous omega to output file. Very heavy output
@write_omavg = ".false."    # If (write_ascii = T) time-averaged growth rate and frequency to the output file.
@write_eigenfunc = ".true."    # If (write_ascii = T) Normalized phi written to runname.eigenfunc
@write_final_fields = ".true."    # If (write_ascii = T) Phi(theta) written to '.fields'
@write_final_moments = ".true."    # write final n, T
@nsave = 500    # Write restart files every nsave timesteps
@nwrite = 1    # Output diagnostic data every nwrite
@navg = 10    # Any time averages performed over navg
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
@omegatinst = 500.0    # Recommended value: 500.
@save_for_restart = ".false."    # Write restart files.
@write_phi_over_time = ".true."    # Write entire Phi field to NetCDF file every nwrite.
@write_ntot_over_time = ".true."    # Write entire Phi field to NetCDF file every nwrite.
@write_moments = ".true."
@write_line = ".true."
@write_fluxes_by_mode = ".true."

