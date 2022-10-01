################################################################################
#
# Description: PWR case
#
#
# External files:
#                power history file BEN013_power.csv
#                axial peaking factor file BEN013_axial_peaking.csv
#
################################################################################
[GlobalParams]
  density = 10412.0  # 95% TD  
  displacements = 'disp_x disp_y'
  order = SECOND
  energy_per_fission = 3.2e-11
[]

[Problem]
  coord_type = RZ
  type = AugmentedLagrangianContactProblem
  reference_vector = 'ref'
  extra_tag_vectors = 'ref'
  maximum_lagrangian_update_iterations = 200
[]

[Mesh]
  [smeared_pellet_mesh]
    type = SmearedPelletMeshGenerator
    clad_mesh_density = customize
    clad_thickness = 5.7e-4
    pellet_mesh_density = customize
    ny_p = 200
    nx_c = 4
    nx_p = 12
    pellet_outer_radius = .0041
    ny_cu = 3
    ny_c = 200
    clad_bot_gap_height = 2.54e-3
    pellet_quantity = 1
    pellet_height = 3.66
    ny_cl = 3
    clad_top_gap_height = 0.28581
    clad_gap_width = 8.7e-5
    elem_type = QUAD8
  []
  patch_size = 20
  patch_update_strategy = auto
  partitioner = centroid
  centroid_partitioner_direction = y
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temp]
    initial_condition = 293
  []
[]

[AuxVariables]
  [fast_neutron_flux]
    block = 1
  []
  [fast_neutron_fluence]
    block = 1
  []
  [grain_radius]
    block = 3
    initial_condition = 5.0e-6
  []
  [effective_creep_strain]
    block = 1
    order = CONSTANT
    family = MONOMIAL
  []
  [gap_cond]
    order = CONSTANT
    family = MONOMIAL
  []
  [oxide_thickness]
    order  = CONSTANT
    family = MONOMIAL
  []
  [rod_burnup]
    order = CONSTANT
    family = MONOMIAL
  []
  [buavg]
    order = CONSTANT
    family = MONOMIAL
  []  
[]

[Functions]
  [power_history]
    type = PiecewiseLinear
#   data_file = BEN013_power.csv
    data_file = power_history.csv
#   direction = right
    format = columns
  []
  [axial_peaking_factors]
    type = PiecewiseBilinear
    # data_file = BEN013_axial_peaking.csv
    data_file = axial_power_profile.csv
    scale_factor = 1
    axis = 1
  []
  [pressure_ramp]
    type = PiecewiseLinear
    x = '-100 0 177922434 177922794'
    y = '0.0065315 1 1 0.0065315'
  []
  [temp_ramp]
    type = PiecewiseLinear
    x = '-100 0 177922434 177922794'
    y = '293 557.15 557.15 293'
  [] 
[]

[Modules/TensorMechanics/Master]
  [pellets]
    block = 3
    strain = FINITE
    eigenstrain_names = 'fuel_relocation_strain fuel_thermal_strain fuel_volumetric_strain'
    generate_output = 'vonmises_stress hydrostatic_stress stress_xx stress_yy stress_zz strain_xx strain_yy strain_zz'
    extra_vector_tags = 'ref'
  []
  [clad]
    block = 1
    strain = FINITE
    eigenstrain_names = 'clad_thermal_eigenstrain clad_irradiation_strain'
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz creep_strain_xx creep_strain_yy creep_strain_xy creep_strain_zz strain_xx strain_yy strain_zz'
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
    block = 3
    fission_rate = fission_rate
    extra_vector_tags = 'ref'
  []
[]

[AuxKernels]
  [./buavg]
    type = SpatialUserObjectAux
    block = 3
    variable = buavg
    execute_on = timestep_end
    user_object = rod_burnup
  [../]
  [fast_neutron_flux]
    type = FastNeutronFluxAux
    variable = fast_neutron_flux
    block = 1
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    factor = 3e13
    execute_on = timestep_begin
  []
  [fast_neutron_fluence]
    type = FastNeutronFluenceAux
    variable = fast_neutron_fluence
    block = 1
    fast_neutron_flux = fast_neutron_flux
    execute_on = timestep_begin
  []   
  [grain_radius]
     type = GrainRadiusAux
     block = 3
     variable = grain_radius
     temperature = temp
     execute_on = linear
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
  [oxide]
    type = MaterialRealAux
    property = oxide_scale_thickness
#      temperature = temp
#      fast_neutron_flux = fast_neutron_flux
    variable = oxide_thickness
    boundary = 2
#      use_coolant_channel = true   # true when oxide_thickness is coupled with coolant channel model
#      oxide_scale_factor = 1.0     # a scale factor to increase oxidation rate
#      model_option = 1
#      lithium_concentration = 1.5 # average Li concentration
#      tin_content = 1.45  # %
#      execute_on  = timestep_end
  []
[]

[Burnup]
  [burnup]
    block = 3
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    num_radial = 81
    num_axial = 11
    a_lower = 0.00478
    a_upper = 3.66
    fuel_inner_radius = 0.0
    fuel_outer_radius = 0.0041
    fuel_volume_ratio = 1
    i_enrich = '0.045 .955 0 0 0 0' # 
    RPF = RPF
  []
[]

[Contact]
  [pellet_clad_mechanical]
    primary = 5
    secondary = 10
    penalty = 1e9
    model = coulomb
    formulation = augmented_lagrange
    friction_coefficient = 0.4
    tangential_tolerance = 1e-3
    normal_smoothing_distance = 0.1
    al_penetration_tolerance = 1e-6
    al_incremental_slip_tolerance = 1e-6
    al_frictional_force_tolerance = 5e-2
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
      factor = 15.51320391e6
      function = pressure_ramp
    []
  []
  [PlenumPressure]
    [plenumPressure]
      boundary = 9
      initial_pressure = 2.72342913e6
      startup_time = 0
      R = 8.3143
      output_initial_moles = initial_moles
      temperature = plenum_temperature
      volume = plenum_volume
      material_input = fission_gas_released
      output = plenum_pressure
      displacements = 'disp_x disp_y'
    []
  []
[]

[CoolantChannel]
  [convective_clad_surface] # apply convective boundary to clad outer surface
    variable = temp
    boundary = '2'
    inlet_temperature = 558   # K
    inlet_pressure = 15.5e6   # Pa
    inlet_massflux = 3600.0   # kg/m^2-sec
    rod_diameter = 0.9515e-2  # m
    rod_pitch = 1.26e-2  # m
    linear_heat_rate  = power_history
    axial_power_profile = axial_peaking_factors
    outputs = all
    output_properties = 'coolant_ehtnalpy coolant_channel_htc  coolant_temperature output_heat_flux critical_heat_flux dnbr'
    oxide_thickness = oxide_thickness 
    axial_offset = 0.0
#   check_boundary_restricted = false
 []
[]

[Materials]
  [fuel_density]
    type = Density
    block = 3
  []
  [fuel_thermal]
    type = ThermalFuel
    block = 3
    thermal_conductivity_model = NFIR
    temperature = temp
    burnup = burnup
  []
  [fuel_elasticity_tensor]
    type = UO2ElasticityTensor
    block = 3
    temperature = temp
  []
  [fuel_elastic_stress]
    type = ComputeFiniteStrainElasticStress
    block = 3
  []
  [fuel_thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = 3
    thermal_expansion_coeff = 10.0e-6
    temperature = temp
    stress_free_temperature = 293.0
    eigenstrain_name = fuel_thermal_strain
  []
  [fuel_relocation]
    type = UO2RelocationEigenstrain
    block = 3
    burnup_function = burnup
    diameter = 0.0095631 #Fuel pellet diameter in m
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    gap = 174.0e-6 #diametral gap in m
    relocation_activation1 = 5000
    burnup_relocation_stop = 0.029
    eigenstrain_name = fuel_relocation_strain
  []
  [fuel_volumetric_swelling]
    type = UO2VolumetricSwellingEigenstrain
    block = 3
    temperature = temp
    burnup = burnup
    initial_fuel_density = 10412.0
    total_densification = 0.01
    initial_porosity = 0.05
    eigenstrain_name = fuel_volumetric_strain
  []
  [ZryOxidation]
    type = ZryOxidation
    boundary = 2   
    clad_inner_radius = 4.187e-3
    clad_outer_radius = 4.757e-3    
    use_coolant_channel = true
    temperature = temp
    normal_operating_temperature_model = pnnl_m5
    fast_neutron_flux = fast_neutron_flux
  []
  [clad_thermal]
    type = HeatConductionMaterial
    block = 1
    thermal_conductivity = 16.0
    specific_heat = 330.0
  []
  [clad_elasticity_tensor]
    type = ZryElasticityTensor
    block = clad
  []
  [clad_stress]
    type = ComputeMultipleInelasticStress
    tangent_operator = elastic
    inelastic_models = 'clad_zrycreep'
    block = clad
  []
  [clad_zrycreep]
    type = ZryCreepLimbackHoppeUpdate
    block = clad
    temperature = temp
    fast_neutron_flux = fast_neutron_flux
    fast_neutron_fluence = fast_neutron_fluence
    model_irradiation_creep = true
    model_primary_creep = true
    model_thermal_creep = true
  []
  [thermal_expansion]
    type = ZryThermalExpansionMATPROEigenstrain
    block = clad
    temperature = temp
    stress_free_temperature = 293.0
    eigenstrain_name = clad_thermal_eigenstrain
  []
  [irradiation_swelling]
    type = ZryIrradiationGrowthEigenstrain
    block = clad
    fast_neutron_fluence = fast_neutron_fluence
    zircaloy_material_type = stress_relief_annealed
    eigenstrain_name = clad_irradiation_strain
  []
  [clad_density]
    type = Density
    block = 1
    density = 6551.0
  []
  [fission_gas_release]
    type = Sifgrs
    block = 3
    temperature = temp
    fission_rate = fission_rate
    grain_radius = grain_radius
    gbs_model = true
    burnup = burnup
    transient_option = 2
  []
[]

[UserObjects]
  [./rod_burnup]
    type = LayeredAverage
    block = 3
    variable = burnup
    direction = y
    num_layers = 48
  [../]
[]


[Dampers]
  [limitT]
    type = MaxIncrement
    variable = temp
    max_increment = 50
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       superlu_dist'

  line_search = 'none'
  verbose = true

  l_max_its = 100
  l_tol = 8e-3

  nl_max_its = 100
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-8

  start_time = -100
  end_time = 177922794

  dtmax = 1e6
  dtmin = 1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e2
    optimal_iterations = 200
    linear_iteration_ratio = 100
    timestep_limiting_function = power_history
    max_function_change = 3e20
    force_step_every_function_point = true
  []

  [Quadrature]
    order = FIFTH
    side_order = SEVENTH
  []
[]

[Postprocessors]
  [clad_inner_vol]
   type = InternalVolume
    boundary = 7
  []
  [fis_gas_grain]
    type = ElementIntegralFisGasGrainSifgrs
    block = 3
    outputs = exodus
  []
  [fis_gas_boundary]
    type = ElementIntegralFisGasBoundarySifgrs
    block = 3
    outputs = exodus
  []
  [flux_from_clad]
    type = SideDiffusiveFluxIntegral
    variable = temp
    boundary = 5
    diffusivity = thermal_conductivity
  []
  [flux_from_fuel]
    type = SideDiffusiveFluxIntegral
    variable = temp
    boundary = 10
    diffusivity = thermal_conductivity
  []
  [average_fission_rate]
    type = ElementAverageValue
    block = 3
    variable = fission_rate
  []
  [rod_ave_lin_pow]
    type = ElementIntegralPower
    block = 3
    fission_rate = fission_rate
    variable = temp
  []
  [disp_y_3023]
    type = NodalVariableValue
    nodeid = 3022
    variable = disp_y
  []
[]

[StandardLWRFuelRodOutputs]
  temperature = temp
  fuel_pellet_blocks = 3
[]

[PerformanceMetricOutputs]
[]


[VectorPostprocessors]
  [axial_temp]
    type = LineValueSampler
    end_point = '0.00 3.66478 0.0'
    start_point = '0.0 0.00478 0.0'
    sort_by = y
    variable = temp
    execute_on = timestep_end
    num_points = 42
    outputs = 'outfile_axial_temp'
  []
  [axial_oxidation]
    type = LineValueSampler
    end_point = '4.757e-3 3.66478   0.0'
    start_point = '4.757e-3 0.00478  0.0'
    sort_by = y
    variable = oxide_thickness
    execute_on = timestep_end
    num_points = 42
    outputs = 'outfile_axial_oxidation'
  []
[]

[Outputs]
  exodus = true
  csv = true
  color = false
  print_linear_residuals = true
  perf_graph = true
  [console]
    type = Console
    max_rows = 40
  []
  [chkfile]
    type = CSV
    show = 'average_centerline_fuel_temperature fission_gas_released_percentage maximum_clad_elongation maximum_fuel_elongation'
#    execute_on = 'FINAL'
    sync_times = '3600 7200 10800 14400 177922434 177922794'
    sync_only = true
  []
  [outfile_axial_temp]
    type = CSV
    execute_on = 'TIMESTEP_END'
  []
  [outfile_axial_oxidation]
    type = CSV
    execute_on = 'FINAL'
  []
[]

[Debug]
  show_var_residual = 'disp_x disp_y temp'
  show_var_residual_norms = true
[]
