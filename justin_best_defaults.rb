######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2014-07-26 11:26:35 +0100                                           #
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

@nx = 192    # The number of kx modes: the number of kx modes actually simulated (ntheta0) is equal to 2*(nx - 1)/3 + 1, due to the need to prevent aliasing.
@ny = 64    # The number of ky modes: the number of ky modes actually simulated (naky) is equal to (ny - 1)/3 + 1, due to the need to prevent aliasing.
@jtwist = 5    # L_x = L_y  jtwist / (2 pi shat)
@y0 = 10.0    # The length of the box in the y direction (measured in the Larmour radius of species 1)


######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 32    # Number of points along field line (theta) per 2 pi segment
@nperiod = 1    # Number of 2 pi segments along equilibrium magnetic field.
@shat = 0.8    # 
@rhoc = 0.54    # Flux surface label. Usually rho = diameter/diameter of LCFS
@qinp = 1.4    # Sets value of the safety factor when using local toroidal equilibrium model.
@akappa = 2.0    # Sets local elongation when local toroidal equilibrium is specified.
@akappri = 0.0    # akappri = dkappa/drho
@tri = 0.0    # tri = arcsin[(R(max(Z))-R_major)/r_mid]
@tripri = 0.0    # tripri =  dtri/drho
@shift = -0.2    # shift = -R q**2 dbeta/drho (>0)
@rmaj = 3.0    # Major radius/a (Position of magnetic axis)
@r_geo = 3.0    # Major radius/a (centerpoint of LCFS)


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
@local_eq = "T"    # 
@eqfile = "dskeq.cdf"    # 
@bishop = 4    # 
@s_hat_input = 0.8    # 
@beta_prime_input = 0.0    # The gradient of the pressure.
@delrho = 0.001    # 
@isym = 0    # 
@ppl_eq = "F"    # 
@gen_eq = "F"    # 
@efit_eq = "F"    # 
@equal_arc = "F"    # 
@writelots = "F"    # 


######################################
# Defaults for namelist le_grids_knobs
#######################################

@vcut = 2.5    # No. of standard deviations from the standard Maxwellian beyond which the distribution function will be set to 0
@ngauss = 5    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 12    # Total number of energy grid points


######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@boundary_option = "linked"    # Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi).
@adiabatic_option = "iphi00=2"    # The form of the adiabatic response (if a species is being modeled as adiabatic).
@g_exb = 0.0    # 
@mach = 0.0    # 


######################################
# Defaults for namelist fields_knobs
#######################################

@field_option = "local"    # Controls which time-advance algorithm is used for the linear terms.
@field_subgath = "T"    # Set to TRUE to use allgatherv to fetch part of the field update calculated on other procs. FALSE uses a sum_allreduce instead.


######################################
# Defaults for namelist knobs
#######################################

@fphi = 1.0    # Multiplies Phi (electrostatic potential).
@fapar = 0.0    # Multiplies A_par. Use 1 for finite beta (electromagnetic), 0 otherwise (electrostatic)
@faperp = 0.0    # Multiplies A_perp. Use 1 for high beta, 0 otherwise. Deprecated: use fbpar instead
@delt = 0.025    # Time step
@nstep = 100000    # Maximum number of timesteps


######################################
# Defaults for namelist reinit_knobs
#######################################

@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted
@delt_minimum = 0.0001    # The minimum time step is delt_minimum.


######################################
# Defaults for namelist layouts_knobs
#######################################

@layout = "xyles"    # 'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.
@local_field_solve = "F"    # Strongly affects initialization time on some parallel computers.
@unbalanced_xxf = "T"    # This allows GS2 to set up an unbalanced xxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_xxf = 0.5    # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@unbalanced_yxf = "T"    # This allows GS2 to set up an unbalanced yxxf processor grid (e.g. leaving some tasks with no work) in order to balance the work load on each.
@max_unbalanced_yxf = 0.5    # This sets the maximum level of difference between the largest and smallest block sizes. Must be between 0 and 1
@opt_redist_nbk = "T"    # This enables the use of non-blocking communication in redistribute routines.
@opt_redist_init = "T"    # This enables optimized initialization routines for creating redistribution objects.
@intmom_sub = "T"    # This enables use of sub-communicators to do reduction associated with calculation of moments of distribution function. Most advantageous for collisional runs without LE layouts.
@intspec_sub = "T"    # This enables use of sub-communicators to do reduction associated with calculation of species integrated moments of distribution function.


######################################
# Defaults for namelist collisions_knobs
#######################################

@collision_model = "none"    # Collision model used in the simulation. Options: 'default', 'none', 'lorentz', 'ediffuse'


######################################
# Defaults for namelist hyper_knobs
#######################################

@hyper_option = "visc_only"    # 
@isotropic_shear = ".false."    # 
@const_amp = ".false."    # Detrmines whether damping rate depends on amplitude variations. Recommend FALSE for nonlinear, TRUE for linear.
@d_hypervisc = 0.1    # Sets hyperviscosity parameter multiplying damping term. See Belli (2006) thesis for more information.


######################################
# Defaults for namelist nonlinear_terms_knobs
#######################################

@nonlinear_mode = "on"    # Include nonlinear terms? ('on','off')
@cfl = 0.25    # The maximum delt < cfl * min(Delta_perp/v_perp)


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
@tprim_1 = 2.3    # -1/T (dT/drho)
@fprim_1 = 0.733    # -1/n (dn/drho)
@uprim_1 = 0.0    # ?
@vnewk_1 = 0.0    # collisionality parameter
@type_1 = "ion"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist species_parameters_2
#######################################

@z_2 = -1.0    # Charge
@mass_2 = 0.00027    # Mass
@dens_2 = 1.0    # Density	
@temp_2 = 1.0    # Temperature
@tprim_2 = 2.3    # -1/T (dT/drho)
@fprim_2 = 0.733    # -1/n (dn/drho)
@uprim_2 = 0.0    # ?
@vnewk_2 = 0.0    # collisionality parameter
@type_2 = "electron"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist dist_fn_species_knobs_1
#######################################

@fexpr_1 = 0.45    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_1 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05


######################################
# Defaults for namelist dist_fn_species_knobs_2
#######################################

@fexpr_2 = 0.45    # Temporal implicitness parameter. Recommended value: 0.48
@bakdif_2 = 0.05    # Spatial implicitness parameter. Recommended value: 0.05


######################################
# Defaults for namelist init_g_knobs
#######################################

@clean_init = "T"    # phi = 0 at either end of domain.
@chop_side = "F"    # Rarely needed. Forces asymmetry into initial condition.
@phiinit = 0.001    # Average amplitude of initial perturbation of each Fourier mode.
@ginit_option = "noise"    # Sets the way that the distribution function is initialized.


######################################
# Defaults for namelist gs2_diagnostics_knobs
#######################################

@print_flux_line = "T"    # Instantaneous fluxes output to screen
@write_nl_flux = "T"    # Write nonlinear fluxes as a function of time.
@print_line = "F"    # Estimated frequencies and growth rates to the screen/stdout
@write_verr = "F"    # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_g = "F"    # Write the distribution function to the '.dist' (NetCDF?)
@write_line = "F"    # If (write_ascii = T) write estimated frequencies and growth rates to the output file
@write_omega = "F"    # If (write_ascii = T) instantaneous omega to output file. Very heavy output
@write_eigenfunc = "T"    # If (write_ascii = T) Normalized phi written to runname.eigenfunc
@write_final_fields = "F"    # If (write_ascii = T) Phi(theta) written to '.fields'
@nwrite = 100    # Output diagnostic data every nwrite
@navg = 50    # Any time averages performed over navg
@omegatol = -0.001    # The convergence has to be better than one part in 1/omegatol
@omegatinst = 500.0    # Recommended value: 500.
@save_for_restart = ".true."    # Write restart files.
@write_ascii = "F"    # 
@write_fields = "F"    # 
@write_symmetry = "T"    # Test the symmetry properties of the GK eqn.
