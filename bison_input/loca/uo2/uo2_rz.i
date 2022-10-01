################################################################################
#
# Description: LOCA UO2 RZ case
#
#
# External files:
#                axial peaking factor file
#
################################################################################

peak_temp = 1200     # K
reflood_time = 5     # s
flow_velocity = 0.05 # m/s
power = 30000 # W/m

[GlobalParams]
  order = SECOND
  family = LAGRANGE
  energy_per_fission = 3.2e-11
  displacements = 'disp_x disp_y'
  volumetric_locking_correction = false
[]

[Problem]
  coord_type = RZ
  type = ReferenceResidualProblem
  group_variables = 'disp_x disp_y'
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  [smeared_pellet_mesh]
    type = SmearedPelletMeshGenerator
    clad_mesh_density = customize
    clad_thickness = 6.1e-4
    pellet_mesh_density = customize
    ny_p = 100
    nx_c = 4
    nx_p = 12
    pellet_outer_radius = .00413
    ny_cu = 3
    ny_c = 100
    clad_bot_gap_height = 2.54e-3
    pellet_quantity = 1
    pellet_height = 3.66
    ny_cl = 3
    clad_top_gap_height = 0.18613
    clad_gap_width = 7.5e-5
    elem_type = QUAD8
  []
  patch_size = 20
  patch_update_strategy = auto
  partitioner = centroid
  centroid_partitioner_direction = y
[]

[DefaultElementQuality]
  aspect_ratio_upper_bound = 253
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temp]
    initial_condition = 300
  []
[]

[AuxVariables] 
  [fast_neutron_flux]
    block = clad
  []
  [fast_neutron_fluence]
    block = clad
  []
  [grain_radius]
    block = pellet
    initial_condition = 7.8e-6 # 2D grain radius
  []
  [effective_creep_strain]
    block = clad
    order = CONSTANT
    family = MONOMIAL
  []
  [fract_beta_phase] # Fraction of beta phase in Zry
    order  = CONSTANT
    family = MONOMIAL
  []
  [creep_rate]
    order = CONSTANT
    family = MONOMIAL
  []
  [creep_rate_aux]
    order = CONSTANT
    family = MONOMIAL
  []
  [bursted]
    order = CONSTANT
    family = MONOMIAL
  []
  [creep_strain_mag]
    order = CONSTANT
    family = MONOMIAL
  []
  [gap_cond]
    order = CONSTANT
    family = MONOMIAL
  []
  [coolant_htc]
    order = CONSTANT
    family = MONOMIAL
  []
  [coolant_temp]
    order = CONSTANT
    family = MONOMIAL
  []
  [hmode]
    order = CONSTANT
    family = MONOMIAL
  []
  [htype]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Functions]
  [h_mode]
    type = PiecewiseConstant
    x = '-10  0  5  200'
    y = '9    8 10  10'
  []
  [htc]
    type = PiecewiseLinear
    x = '-10   -1.0  0.0    5   200'
    y = '30000 30000 30000 0.01  0.01'
  []
  [coolant_temp]
    type = PiecewiseLinear
     x = '-10     0     5     200'
     y = '570    570    373   373'
  []
  [coolant_pressure]
    type = PiecewiseLinear
    x = '-10       0      5  200'
    y = '15.5e6 15.5e6  0.2e6 0.2e6'  
  []   
#  [axial_peaking_factors]
#    type = ConstantFunction
#    value = 1
#  []

  [axial_peaking_factors]
    type = PiecewiseBilinear
    data_file = axial_power_profile.csv
    scale_factor = 1
    axis = 1
  []
  [power_history]    
    type = ParsedFunction
    value = '${power} * if( (t > 0.0), 1.2*0.066/(t-0.0)^0.2, 1)'
  [../]
  [pressure_ramp]             
    type = PiecewiseLinear
    data_file = 'coolant_pressure.csv'
    format = columns
  []
[]



[Modules/TensorMechanics/Master]
  [pellets]
    block = pellet
    strain = FINITE
    incremental = true
    eigenstrain_names = 'fuel_thermal_strain fuel_volumetric_strain'
    cylindrical_axis_point1 = '0 0 0'
    cylindrical_axis_point2 = '0 1 0'
    generate_output = 'vonmises_stress hydrostatic_stress stress_xx stress_yy
      stress_zz elastic_strain_yy strain_xx strain_yy strain_zz hoop_stress'
    extra_vector_tags = 'ref'
  []
  [clad]
    block = clad
    strain = FINITE
    incremental = true
    eigenstrain_names = 'clad_thermal_eigenstrain'
    cylindrical_axis_point1 = '0 0 0'
    cylindrical_axis_point2 = '0 1 0'
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz
      creep_strain_xx creep_strain_yy creep_strain_xy creep_strain_zz
      elastic_strain_xx elastic_strain_yy elastic_strain_zz strain_xx strain_yy
      strain_zz hoop_stress' #plastic_strain_xx plastic_strain_yy plastic_strain_zz
    extra_vector_tags = 'ref'
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = temp
    extra_vector_tags = 'ref'
  []
  [heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
    extra_vector_tags = 'ref'
  []
  [heat_source]
     type = NeutronHeatSource
     variable = temp
     block = pellet
     fission_rate = fission_rate
     extra_vector_tags = 'ref'
  []
[]

[AuxKernels]
  [fast_neutron_flux]
    type = FastNeutronFluxAux
    variable = fast_neutron_flux
    block = clad
    axial_power_profile = axial_peaking_factors
    factor = 0.16e15 #n/m2-s
    execute_on = timestep_begin
  []
  [grain_radius]
    type = GrainRadiusAux
    block = pellet
    variable = grain_radius
    temperature = temp
    execute_on = linear
  []
  [fast_neutron_fluence]
    type = FastNeutronFluenceAux
    block = clad
    variable = fast_neutron_fluence
    fast_neutron_flux = fast_neutron_flux
    execute_on = timestep_begin
  []
  [creep_strain_mag]
    type = MaterialRealAux
    property = effective_creep_strain
    variable = creep_strain_mag
    block = clad
    execute_on = timestep_end
  []
  [effective_creep_strain]
    type = MaterialRealAux
    property = effective_creep_strain
    variable = effective_creep_strain
    block = clad
    execute_on = timestep_end
  []
  [conductance]
    type = MaterialRealAux
    property = gap_conductance
    variable = gap_cond
    boundary = 10
  []
  [coolant_htc]
    type = MaterialRealAux
    property = coolant_channel_htc
    variable = coolant_htc
    boundary = 2
  []
  [coolant_temp]
    type = MaterialRealAux
    property = coolant_temperature
    variable = coolant_temp
    boundary = 2
  []
  [hmode]
    type = MaterialRealAux
    property = coolant_channel_hmode
    variable = hmode
    boundary = 2
  []
  [htype]
    type = MaterialRealAux
    property = coolant_channel_htype
    variable = htype
    boundary = 2
  []
  [fract_bphase]
    type = MaterialRealAux
    variable = fract_beta_phase
    property = fract_beta_phase
    block = clad
  []
  [creep_rate]
    type = MaterialRealAux
    variable = creep_rate
    property = creep_rate
    block = clad
    execute_on = timestep_end
  []
  [creep_rate_aux]
    type = MaterialRealAux
    variable = creep_rate_aux
    property = creep_rate
    block = clad
    execute_on = timestep_end
  []
  [bursted]
    type = MaterialRealAux
    variable = bursted
    property = failed
    boundary = 2
    execute_on = timestep_end
  []
[]

# TODO: Have StandardLWRFuelRodOutputs create this when the feature in issue #1054 is
#       developed.
#      We are using 'plenum_temp' rather than 'plenum_temperature', which is generated
#      automatically by StandardLWRFuelRodOutputs, but computed in a different way.
[PlenumTemperature]
  [plenum_temp]
    boundary = 5
    inner_surfaces = '5'
    outer_surfaces = '10'
    temperature = temp
  []
[]

[Burnup]
  [burnup]
    block = pellet
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    num_radial = 81
    num_axial = 11
    a_lower = 0.00478
    a_upper = 3.66478
    fuel_inner_radius = 0.0
    fuel_outer_radius = 0.00413 # m
    fuel_volume_ratio = 1.0
    i_enrich = '0.0293 .9707 0 0 0 0' # 3.67% enriched U-235 #TODO: Looks like it's set for 2.93%!
    RPF = RPF
    density = 10431 #95 %TD  Assume TD = 10980 kg/cm3
  []
[]

[Contact]
  [pellet_clad_mechanical]
    primary = 5
    secondary = 10
    penalty = 1e11
    normalize_penalty = true
    model = frictionless
#    model = coulomb
    formulation = penalty
#    friction_coefficient = 1.0
    tangential_tolerance = 1e-3
    normal_smoothing_distance = 0.1
  []
[]

[ThermalContact]
  [thermal_contact]
    type = GapHeatTransferLWR
    variable = temp
    primary = 5
    secondary = 10
    initial_moles = initial_moles
    gas_released = fission_gas_released
    jump_distance_model = KENNARD
    plenum_pressure = plenum_pressure
    contact_pressure = contact_pressure
    roughness_fuel = 2e-6
    roughness_clad = 1e-6
    roughness_coef = 3.2
    normal_smoothing_distance = 0.1
    quadrature = true
  []
[]

[BCs]
  [no_x_all]
    type = DirichletBC
    variable = disp_x
    boundary = 12
    value = 0.0
  []
  [no_y_clad_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = '1'
    value = 0.0
  []
  [no_y_fuel_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = '1020'
    value = 0.0
  []
  [Pressure]
    [coolantPressure]
      boundary = '1 2 3'
      factor = 1.0 # Pa
      function = pressure_ramp
    []
  []
  [PlenumPressure]
    [plenumPressure]
      boundary = 9  # clad interior + fuel exterior
      initial_pressure = 2.0e6 # Pa
      startup_time = 0
      R = 8.3143
      output_initial_moles = initial_moles
      temperature = plenum_temp
      volume = plenum_volume
      material_input = fission_gas_released
      output = plenum_pressure
    []
  []
[]

[CoolantChannel]
# convective boundary condition at clad outer surface
  [./clad_outer_surface]
    boundary          =  2
    variable          =  temp
    inlet_temperature =  coolant_temp      # K
    inlet_pressure    =  coolant_pressure  # Pa
    inlet_massflux = 3600.0   # kg/m^2-sec
    rod_diameter = 0.952e-2  # m
    rod_pitch    = 1.26e-2    # m    
    heat_transfer_mode = h_mode
    heat_transfer_coefficient = htc # W/m^2-K
    htc_correlation_type = 1
    flooding_time     =  ${reflood_time}
    flooding_rate     = ${flow_velocity}
    initial_temperature = ${peak_temp} # K
    # axial_offset = 0.0
    axial_scale_factor = 1.83
    initial_power     = 4.0    # kW/m
    blockage_ratio    = 0.0    #
    fuel_stack_length = 3.66   # m
#   input_Tmin        = 1000
#   input_rewetting_htc = 1.0e5
    reflooding_model  = 1
    outputs = all
    output_properties = 'coolant_channel_htc output_heat_flux critical_heat_flux dnbr'
    compute_enthalpy  = false
  [../]
[]

[Materials]
  [fuel_thermal]  # temperature and burnup dependent thermal properties of UO2
    type = ThermalFuel
    block = pellet
    thermal_conductivity_model = NFIR
    temperature = temp
    burnup = burnup
  []
  [fuel_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    block = pellet
    youngs_modulus = 2.0e11
    poissons_ratio = 0.345
  []
  [fuel_elastic_stress]
    type = ComputeFiniteStrainElasticStress
    block = pellet
  []
  [fuel_thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = pellet
    thermal_expansion_coeff = 10.0e-6
    temperature = temp
    stress_free_temperature = 300
    eigenstrain_name = fuel_thermal_strain
  []
  [fuel_volumetric_swelling]
    type = UO2VolumetricSwellingEigenstrain
    gas_swelling_model_type = SIFGRS
    block = pellet
    temperature = temp
    burnup = burnup
    initial_fuel_density = 10431.0 #95 %TD  Assume TD = 10980 kg/cm3
    eigenstrain_name = fuel_volumetric_strain
  []
  [fission_gas_release]
    type = Sifgrs
    block = pellet
    temperature = temp
    fission_rate = fission_rate        # coupling to fission_rate aux variable
#    initial_grain_radius = 6.552e-6    # 2D grain radius 4.2e-6
    grain_radius = grain_radius
    gbs_model = true
    burnup = burnup
#    compute_swelling = true
    transient_option = 2
  []
  [fuel_density]
    type = Density
    block = pellet
    density = 10431 #95 %TD  Assume TD = 10980 kg/cm3
  []
  [clad_thermal]
    type = HeatConductionMaterial
    block = clad
    thermal_conductivity = 16.0
    specific_heat = 330.0
  []
  [clad_elasticity_tensor]
    type = ZryElasticityTensor
    block = clad
    temperature = temp
  []
  [clad_stress]
    type = ComputeMultipleInelasticStress
    tangent_operator = elastic
    inelastic_models = 'clad_zrycreep'
    block = clad
  []
  [clad_zrycreep]
    type = ZryCreepLOCAErbacherLimbackHoppeUpdate
    block = clad
    temperature = temp
    fast_neutron_flux = fast_neutron_flux
    fast_neutron_fluence = fast_neutron_fluence
    model_irradiation_creep = false
    model_primary_creep = false
    model_thermal_creep = true
    temperature_standard_thermal_creep_end = 700.0
    temperature_loca_creep_begin = 900.0
    max_inelastic_increment = 1e-4
  []
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = clad
    temperature = temp
    thermal_expansion_coeff = 5.0e-6
    stress_free_temperature = 300
    eigenstrain_name = clad_thermal_eigenstrain
  []
  [phase]
    type = ZrPhase
    block = clad
    temperature = temp
    numerical_method = 2
  []
  [failure_criterion]
    type = ZryCladdingFailure
    boundary = '2'
    failure_criterion = combined_overstress_and_plastic_instability
    hoop_stress = hoop_stress
    effective_strain_rate_creep = creep_rate
    temperature = temp
    fraction_beta_phase = fract_beta_phase
    outputs = all
    output_properties='failed burst_stress'
  []
  [clad_density]
    type = Density
    block = clad
    density = 6551.0
  []
[]

[Dampers]
  [limitT]
     type = MaxIncrement
     variable = temp
     max_increment = 50
  []
[]

[Executioner]
  type = Transient

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       superlu_dist'

  line_search = 'none'
  verbose = true

  # controls for linear iterations
  l_max_its = 100
  l_tol = 8e-3

  # controls for nonlinear iterations
  nl_max_its = 50
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-10

  dt         = 0.1
  start_time = -50.0
  end_time   =  600

  dtmax = 10
  dtmin = 0.0001
  
  #[./TimeStepper]
  #  type = IterationAdaptiveDT
  #  dt = 0.1
  #  optimal_iterations = 200
  #  linear_iteration_ratio = 100
  #  time_t  = '-10 0    5   10  100  300'
  #  time_dt = '1.0 1.0 0.1 0.1  1.0  1.0'

   # #timestep_limiting_function = time_step_function
   # #max_function_change = 1e20
   # #force_step_every_function_point = true
  # [../]  
  
  [TimeStepper]
    type = PostprocessorDT
    postprocessor = material_timestep
    dt = 0.01
  []
  
  [Quadrature]
    order = FIFTH
    side_order = SEVENTH
  []

[]

[Postprocessors]
  [ave_temp_interior]            # average temperature of the cladding interior and all pellet exteriors
     type = SideAverageValue
     boundary = 9
     variable = temp
     execute_on = 'initial linear'
   []
  [avg_clad_temp]               # average temperature of cladding interior
    type = SideAverageValue
    boundary = 7
    variable = temp
    execute_on = 'initial timestep_end'
  []
   [fis_gas_released]
     type = ElementIntegralFisGasReleasedSifgrs
     block = pellet
   []
  [fis_gas_grain]
    type = ElementIntegralFisGasGrainSifgrs
    block = pellet
    outputs = exodus
    execute_on=linear
  []
  [fis_gas_boundary]
    type = ElementIntegralFisGasBoundarySifgrs
    block = pellet
    outputs = exodus
    execute_on=linear
  []
  [max_betaph_fract]
    type = ElementExtremeValue
    value_type = max
    variable = fract_beta_phase
  []
  [flux_from_clad]           # area integrated heat flux from the cladding
    type = SideDiffusiveFluxIntegral
    variable = temp
    boundary = 5
    diffusivity = thermal_conductivity
    execute_on=timestep_end
  []
  [flux_from_fuel]          # area integrated heat flux from the fuel
    type = SideDiffusiveFluxIntegral
    variable = temp
    boundary = 10
    diffusivity = thermal_conductivity
    execute_on=timestep_end
  []
  [average_fission_rate]
    type = ElementAverageValue
    block = pellet
    variable = fission_rate
    execute_on=timestep_end
  []
  [rod_ave_lin_pow]
    type = ElementIntegralPower
    block = pellet
    fission_rate = fission_rate
    variable = temp
    execute_on=timestep_end
  []
  [rod_input_power]
    type = FunctionValuePostprocessor
    function = power_history
    scale_factor = 3.66                   # rod height
    execute_on=timestep_end
  []
  [material_timestep]
    type = MaterialTimeStepPostprocessor
    block = clad
  []
  [max_creep_rate]
    type = ElementExtremeValue
    block = clad
    value_type = max
    variable = creep_rate_aux
  []
  [bursted]
    type = ElementExtremeValue
    block = clad
    value_type = max
    variable = bursted
  []
[]

[UserObjects]
  [terminator]
    type = Terminator
    expression = 'bursted > 0'
  []
[]

[StandardLWRFuelRodOutputs]
  fuel_pellet_blocks = 3
  temperature = temp
[]

[PerformanceMetricOutputs]
[]

[Outputs]
  exodus = true
  csv = true
  color = false
  perf_graph = true
  [console]
    type = Console
    output_linear = true
    max_rows = 40
  []
[]

[Debug]
  show_var_residual = 'disp_x disp_y temp'
  show_var_residual_norms = true
[]
