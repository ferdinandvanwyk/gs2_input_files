######################################################################
#   Automatically generated defaults file for GS2 CodeRunner module  #
#                                                                    #
# This defaults file specifies a set of defaults for GS2 which are   #
# used by CodeRunner to set up and run GS2 simulations.              #
#                                                                    #
# Created 2015-02-23 17:40:48 +0000                                           #
#                                                                    #
######################################################################

@defaults_file_description = ""


######################################
# Defaults for namelist parameters
#######################################

@beta = 0.0    # Ratio of particle to magnetic pressure (reference Beta, not total beta):  beta=n_0 T_0 /( B^2 / (8 pi))
@zeff = 1.7    # Effective ionic charge.


######################################
# Defaults for namelist kt_grids_knobs
#######################################

@grid_option = "single"    # The general layout of the perpendicular grid.


######################################
# Defaults for namelist theta_grid_parameters
#######################################

@ntheta = 24    # Number of points along field line (theta) per 2 pi segment
@nperiod = 2    # Number of 2 pi segments along equilibrium magnetic field.
@shat = 0.29167    # 
@rhoc = 0.378    # Flux surface label. Usually rho = diameter/diameter of LCFS
@qinp = 1.1244    # Sets value of the safety factor when using local toroidal equilibrium model.
@akappa = 1.3063    # Sets local elongation when local toroidal equilibrium is specified.
@akappri = 0.04    # akappri = dkappa/drho
@tri = 0.0171    # tri = arcsin[(R(max(Z))-R_major)/r_mid]
@tripri = 0.029065700682713667    # tripri =  dtri/drho
@shift = -0.14681    # Sets Shafranov shift. See online help for definition.
@rmaj = 3.1869    # Major radius/a (Position of magnetic axis)
@r_geo = 3.1869    # Major radius/a (centerpoint of LCFS)


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
@eqfile = "../../JET75225_EQDSK_CHEASE_RES_400"    # 
@bishop = 4    # 
@s_hat_input = 0.29167    # 
@beta_prime_input = -0.13479869432580815    # The gradient of the pressure.
@delrho = 0.001    # 
@isym = 0    # 
@ppl_eq = "F"    # 
@gen_eq = "F"    # 
@efit_eq = ".false."    # 
@equal_arc = "F"    # 
@writelots = "F"    # 
@chs_eq = ".false."    # Use equilbrium data from the CHEASE file ogyropsi.dat


######################################
# Defaults for namelist le_grids_knobs
#######################################

@vcut = 2.5    # No. of standard deviations from the standard Maxwellian beyond which the distribution function will be set to 0
@ngauss = 10    # Number of untrapped pitch-angles moving in one direction along field line.
@negrid = 32    # Total number of energy grid points


######################################
# Defaults for namelist dist_fn_knobs
#######################################

@gridfac = 1.0    # Affects boundary condition at end of theta grid.
@boundary_option = "default"    # Sets the boundary condition along the field line (i.e. the boundary conditions at theta = +- pi).
@adiabatic_option = "iphi00=2"    # The form of the adiabatic response (if a species is being modeled as adiabatic).


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
@fbpar = 0.0    # Multiplies B_parallel. Use 1 for high beta, 0 otherwise.
@delt = 0.01    # Time step
@nstep = 50000    # Maximum number of timesteps
@avail_cpu_time = 10800    # Specify the available wall clock time in seconds. GS2 will exit before this time.
@trinity_linear_fluxes = ".false."    # If true and running linearly, return linear diffusive flux estimates to Trinity.


######################################
# Defaults for namelist reinit_knobs
#######################################

@delt_adj = 2.0    # When the time step needs to be changed, it is adjusted
@delt_minimum = 1.0e-06    # The minimum time step is delt_minimum.
@delt_cushion = 10.0    # 


######################################
# Defaults for namelist layouts_knobs
#######################################

@layout = "lexys"    # 'yxles', 'lxyes', 'lyxes', 'lexys' Determines the way the grids are laid out in memory.
@local_field_solve = "F"    # Strongly affects initialization time on some parallel computers.


######################################
# Defaults for namelist collisions_knobs
#######################################

@collision_model = "none"    # Collision model used in the simulation. Options: 'default', 'none', 'lorentz', 'ediffuse'


######################################
# Defaults for namelist hyper_knobs
#######################################

@hyper_option = "none"    # 
@isotropic_shear = ".false."    # 
@const_amp = ".true."    # Detrmines whether damping rate depends on amplitude variations. Recommend FALSE for nonlinear, TRUE for linear.
@d_hypervisc = 0.0    # Sets hyperviscosity parameter multiplying damping term. See Belli (2006) thesis for more information.


######################################
# Defaults for namelist nonlinear_terms_knobs
#######################################

@nonlinear_mode = "off"    # Include nonlinear terms? ('on','off')
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
@temp_1 = 1.19    # Temperature
@tprim_1 = 2.56    # -1/T (dT/drho)
@fprim_1 = 1.02    # -1/n (dn/drho)
@uprim_1 = 0.0    # ?
@vnewk_1 = 0.0001077420372277686    # collisionality parameter
@type_1 = "ion"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist species_parameters_2
#######################################

@z_2 = -1.0    # Charge
@mass_2 = 0.000272    # Mass
@dens_2 = 1.0    # Density	
@temp_2 = 1.0    # Temperature
@tprim_2 = 1.47    # -1/T (dT/drho)
@fprim_2 = 1.02    # -1/n (dn/drho)
@uprim_2 = 0.0    # ?
@vnewk_2 = 0.009734429700734465    # collisionality parameter
@type_2 = "electron"    # Type of species, e.g. 'ion', 'electron', 'beam'


######################################
# Defaults for namelist dist_fn_species_knobs_1
#######################################

@fexpr_1 = 0.48    # Temporal implicitness parameter. Recommended value: 0.48
@fexpi_1 = 0.0    # 
@bakdif_1 = 0.02    # Spatial implicitness parameter. Recommended value: 0.05


######################################
# Defaults for namelist dist_fn_species_knobs_2
#######################################

@fexpr_2 = 0.48    # Temporal implicitness parameter. Recommended value: 0.48
@fexpi_2 = 0.0    # 
@bakdif_2 = 0.02    # Spatial implicitness parameter. Recommended value: 0.05


######################################
# Defaults for namelist init_g_knobs
#######################################

@chop_side = ".false."    # Rarely needed. Forces asymmetry into initial condition.
@phiinit = 0.001    # Average amplitude of initial perturbation of each Fourier mode.
@restart_file = "v_resubmit_id_32_nwrite_200_tprim_1_2.56_tprim_2_1.47_fprim_1_1.02_fprim_2_1.02_id_38.nc"    # Base of filenames with restart data.
@ginit_option = "noise"    # Sets the way that the distribution function is initialized.
@zf_init = 0.0    # 
@restart_dir = "nc"    # 
@phifrac = 0.1    # 


######################################
# Defaults for namelist gs2_diagnostics_knobs
#######################################

@print_flux_line = ".false."    # Instantaneous fluxes output to screen
@write_nl_flux = ".true."    # Write nonlinear fluxes as a function of time.
@print_line = "F"    # Estimated frequencies and growth rates to the screen/stdout
@write_verr = "T"    # Write velocity space diagnostics to '.lpc' and '.verr' files
@write_line = "F"    # If (write_ascii = T) write estimated frequencies and growth rates to the output file
@write_gyx = "F"    # Write dist fn at a given physical spacial point to a file
@write_hrate = "F"    # Write heating rate, collisonal entropy generation etc to '.heat'
@write_avg_moments = ".true."    # Write flux surface averaged low-order moments of g to runname.out.nc and runname.moments (if write_ascii = T)
@write_omega = ".true."    # If (write_ascii = T) instantaneous omega to output file. Very heavy output
@write_final_fields = ".true."    # If (write_ascii = T) Phi(theta) written to '.fields'
@write_final_moments = "T"    # write final n, T
@nwrite = 200    # Output diagnostic data every nwrite
@navg = 1000    # Any time averages performed over navg
@omegatol = 0.0001    # The convergence has to be better than one part in 1/omegatol
@omegatinst = 500.0    # Recommended value: 500.
@save_for_restart = ".false."    # Write restart files.
@write_correlation = "T"    # Write parallel correlation.


######################################
# Defaults for namelist kt_grids_range_parameters
#######################################

@theta0 = 0
@aky = 0.778    # 


######################################
# Defaults for namelist diagnostics_config
#######################################

@write_omega = ".true."    # Write growth rates and frequencies to the netcdf file
@navg = 1000    # Any time averages performed over navg
@omegatinst = 500.0    # Growth rates > omegatinst assumed numerical instability.
@omegatol = 0.0001    # The convergence has to be better than one part in 1/omegatol
