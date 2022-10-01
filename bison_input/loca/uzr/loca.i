#====================================================================================================
# A generalized plane strain input for helical metallic-fuel 
# 
# Fuel design parameters:
#
#   Total Fuel Area =  71.10 mm^2
#   Fuel Area =  50.10 mm^2
#   Fuel Area (excluding displacer) =  47.67 mm^2 
#
#   Equivalent inner radius =  0.880 mm
#   Equivalent outer fuel radius  =  3.993 mm
#   Equivalent outer radius  =  4.757 mm
#   Equivalent outer diameter  =  9.515 mm
#
#   Displacer width: 1.56 mm
#   Rod-to-rod pitch: 12.6 mm
#
# Initial Temperature: 300 K
#
# Element type: QUAD-8
# Number of elements: 576
# Number of nodes: 1833
# Number of DOFs: 5501
#
#
# B.C. for LOCA 
#
#====================================================================================================
peak_temp = 800      # K
reflood_time = 5     # s
flow_velocity = 0.05 # m/s
power = 30400 # W/m

[GlobalParams]
  temperature = temp
  displacements = 'disp_x disp_y'
  order = SECOND
  family = LAGRANGE
  energy_per_fission = 3.2e-11  # J/fission
[]

[Mesh]
  file = 2d_r2.e
[]

[Problem]
  type = ReferenceResidualProblem
  reference_vector = 'ref'
  extra_tag_vectors = 'ref'
  group_variables = 'disp_x disp_y'
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temp]
    initial_condition = 300.0     # set initial temp to ambient
  []
  [./scalar_strain_zz]
    order = SECOND
    family = SCALAR
  [../]
[]

[AuxVariables]
  [fract_beta_phase] # Fraction of beta phase in Zry
    order  = CONSTANT
    family = MONOMIAL
  []
  [fast_neutron_flux]
    block = '2 3'
  []
  [fast_neutron_fluence]
    block = '2 3'
  []  
  [oxide_thickness]
    order = CONSTANT
    family = MONOMIAL
  []
  [./saved_x]
    order = SECOND
    family = LAGRANGE
  [../]
  [./saved_y]
    order = SECOND
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [fract_bphase]
    type = MaterialRealAux
    variable = fract_beta_phase
    property = fract_beta_phase
    block = '2 3'
  []
  [oxide]
    type = MaterialRealAux
    property = oxide_scale_thickness
    variable = oxide_thickness
    boundary = 3
  []
  [fast_neutron_flux]
    type = FastNeutronFluxAux
    variable = fast_neutron_flux
    block = '2 3'
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    factor = 3e13
    execute_on = timestep_begin
  []
  [fast_neutron_fluence]
    type = FastNeutronFluenceAux
    variable = fast_neutron_fluence
    block = '2 3'
    fast_neutron_flux = fast_neutron_flux
    execute_on = timestep_begin
  [] 
  [./saved_x]
    type = TagVectorAux
    variable = 'saved_x'
    vector_tag = 'ref'
    v = 'disp_x'
    execute_on = timestep_end
  [../]
  [./saved_y]
    type = TagVectorAux
    variable = 'saved_y'
    vector_tag = 'ref'
    execute_on = timestep_end
    v = 'disp_y'
  [../]
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
  [axial_peaking_factors]
    type = PiecewiseBilinear
    data_file = axial_power_profile.csv
    scale_factor = 1
    axis = 2
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
  [fuel]
    block = '1'
    strain = FINITE
    eigenstrain_names = 'fuel_thermal_eigenstrain swelling'    
	#decomposition_method = EigenSolution	
	planar_formulation = GENERALIZED_PLANE_STRAIN
    scalar_out_of_plane_strain = scalar_strain_zz
	generate_output = 'vonmises_stress stress_xx stress_yy stress_zz hydrostatic_stress'    
    extra_vector_tags = 'ref'
  []
  [clad]
    block = '2 3'
    strain = FINITE
    planar_formulation = GENERALIZED_PLANE_STRAIN
    eigenstrain_names = 'clad_thermal_eigenstrain'
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz hydrostatic_stress'
    # decomposition_method = EigenSolution
	scalar_out_of_plane_strain = scalar_strain_zz
	extra_vector_tags = 'ref'
  []
[]

[Kernels]
  [heat]         # gradient term in heat conduction equation
    type = HeatConduction
    variable = temp
    extra_vector_tags = 'ref'
  []
  [heat_ie]       # time term in heat conduction equation
    type = HeatConductionTimeDerivative
    variable = temp
    extra_vector_tags = 'ref'
  []
  [heat_source]
    type=NeutronHeatSource
    variable = temp
    block = '1'
    extra_vector_tags = 'ref'
    fission_rate = fission_rate
  []
[]

[BCs]
# Define boundary conditions
  [no_y_all] # pin pellets and clad along axis of symmetry (y)
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  []
  [no_x_all] # pin pellets and clad along axis of symmetry (x)
    type = DirichletBC
    variable = disp_x
    boundary = 2
    value = 0.0
  []
  [Pressure] #  apply coolant pressure on clad outer walls
    [coolantPressure]
      boundary = '3'
      factor = 1.0
      function = pressure_ramp 
    []
  []  
[] 

[CoolantChannel]
# convective boundary condition at clad outer surface
  [./clad_outer_surface]
    boundary          =  3
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
    axial_offset = 1.83
    initial_power     = 4.0    # kW/m
    blockage_ratio    = 0.0    #
    fuel_stack_length = 3.66   # m
#   input_Tmin        = 1000
#   input_rewetting_htc = 1.0e5
    reflooding_model  = 1
    outputs = all
    output_properties = 'coolant_channel_htc output_heat_flux critical_heat_flux dnbr'
    compute_enthalpy  = false
    flow_direction = z
  [../]
[]

[Materials] 
   [burnup]
     type = UPuZrBurnup
     initial_X_Zr = 0.67
     initial_X_Pu = 0.0
     density = 9870
     block = '1'
     outputs = all
     output_properties = burnup
   []  
   [fission_rate]
    type = UZrFissionRate
    rod_linear_power = power_history
    axial_power_profile = axial_peaking_factors
    pellet_radius       = 3.993e-3
    pellet_inner_radius = 0.88e-3
    axial_offset = 1.83
    outputs = all
    output_properties = fission_rate
    block = '1'
  []  
  [thermalUZr]
    enable = true
    type = UZrThermal
    temperature = temp    
    block = '1'
    porosity = porosity
    spheat_model = utk
    thcond_model = utk
  [] 
  [fuel_density]
    type = Density
    block = '1'
    density = 9870.0
  []  
  [swelling]
     type = UZrVolumetricSwellingEigenstrain
     # burnup = 0
     hydrostatic_stress = hydrostatic_stress
     eigenstrain_name = 'swelling'
     output_properties = 'porosity gaseous_porosity'
     gas_swelling_scale_factor = 0.0
     block = '1'
     outputs = all
  []  
  [fuel_elasticity_tensor]
    type = UZrElasticityTensor
    block = '1'
    temperature = temp
  []
  [fuel_stress]
    type = ComputeMultipleInelasticStress
    block = '1'
    inelastic_models = 'fuel_creep'
  []
  [fuel_creep]
    type = UZrCreepUpdate
    temperature = temp
    porosity = porosity
    fission_rate = fission_rate
    max_inelastic_increment = 1e-2
    creep_model = bmi
    block = '1'
  []
  [fuel_plasticity]
    enable = false
    type = UZrPlasticityUpdate
    temperature = temp
    initial_fast_fluence = 0.0    
    strain_rate = 1.0e-4    
    block = '1'
  []
  [fuel_thermal_expansion]
     type = UZrThermalExpansionEigenstrain
     block = '1'
     temperature = temp
     stress_free_temperature = 300.0
     eigenstrain_name = 'fuel_thermal_eigenstrain'
  [] 
  [clad_thermal]
    type = ThermalZry
    block = '2 3'
  []
  [clad_elasticity_tensor]
    type = ZryElasticityTensor
    block = '2 3'
  []
  [clad_creep_model]
    type = ZryCreepLOCAErbacherLimbackHoppeUpdate
    block = '2 3'
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
  [clad_stress]
    type = ComputeMultipleInelasticStress
    block = '2 3'
    tangent_operator = elastic
    inelastic_models = 'clad_creep_model'
  []
  [clad_thermal_expansion]
     type = ZryThermalExpansionMATPROEigenstrain
     block = '2 3'
     temperature = temp
     stress_free_temperature = 300.0
     eigenstrain_name = 'clad_thermal_eigenstrain'
  []  
  [clad_oxidation]
     type = ZryOxidation
     clad_inner_radius = 3.993e-3
     clad_outer_radius = 4.757e-3
     temperature = temp
     fast_neutron_flux = fast_neutron_flux
     use_coolant_channel = true
     normal_operating_temperature_model = pnnl_m5
     boundary = '3'   
  []
  [clad_density]
    type = Density
    block = '2 3'
    density = 6551.0
  [] 
  [phase]
    type = ZrPhase
    block = '2 3'
    temperature = temp
    numerical_method = 2
  []
[]

[Dampers]
  [limitT]
    type = MaxIncrement
    max_increment = 100.0
    variable = temp
  []
  [limitX]
    type = MaxIncrement
    max_increment = 1e-5
    variable = disp_x
  []
[]
[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist'
  line_search = 'none'

  l_max_its = 100
  l_tol = 8e-3

  nl_max_its = 20
  nl_rel_tol = 1e-2
  nl_abs_tol = 1e-8

  dt         = 0.1
  start_time = -50.0
  end_time   =  600

  dtmax = 10
  dtmin = 0.0001
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 200
    linear_iteration_ratio = 100
    time_t  = '-10 0    5   10  100  300'
    time_dt = '1.0 1.0 0.1 0.1  1.0  1.0'

    #timestep_limiting_function = time_step_function
    #max_function_change = 1e20
    #force_step_every_function_point = true
  [../]  
  [Quadrature]
    order = THIRD
  [] 
[]
[Outputs]
  perf_graph = true
  exodus = true
  csv = true
  [console]
    type = Console
    max_rows = 25
  []
[]

[Postprocessors]
  [flux_from_clad]                     # area integrated heat flux from the cladding
    type = SideDiffusiveFluxIntegral
    variable = temp
    boundary = 3
    diffusivity = thermal_conductivity
    execute_on = timestep_end
  [] 
  [_dt]                     # time step
    type = TimestepSize
    execute_on = timestep_end
  []
  [num_lin_it]
    enable = false
    type = NumLinearIterations
  []
  [num_nonlin_it]
    enable = false
    type = NumNonlinearIterations
  []
  [tot_lin_it]
    enable = false
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  []
  [tot_nonlin_it]
    enable = false
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  []
  [alive_time]
    enable = false
    type = PerfGraphData
    section_name = Root
    data_type = TOTAL
  []
  [rod_total_power]
     type = ElementIntegralPower
     variable = temp
     fission_rate = fission_rate
     block = '1'
     execute_on = timestep_end
  []  
  [porosity]
    type = ElementAverageValue
    variable = porosity
    block = '1'
  []
  [peak_oxide_thickness]
    type = ElementExtremeValue
    variable = oxide_thickness
    block = '2'
  []
  [peak_heat_flux]
    type = ElementExtremeValue
    variable = output_heat_flux
    block = '2'
  []
  [peak_temperature]
    type = NodalMaxValue
    variable = temp
    block = '1'
  []
[]
[Postprocessors]
  [node_temp_1813]
    type = NodalVariableValue
    variable =  temp
    nodeid   =  1813
  []
  [node_temp_1182]
    type = NodalVariableValue
    variable =  temp
    nodeid   =  1182
  []
  [node_temp_1255]
    type = NodalVariableValue
    variable =  temp
    nodeid   =  1255
  []
  [node_temp_63]
    type = NodalVariableValue
    variable =  temp
    nodeid   =  63
  []
  [node_temp_0]
    type = NodalVariableValue
    variable =  temp
    nodeid   =  0
  []
[]
[Postprocessors]
  [element_oxide_thickness_353]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  353
  []
  [element_oxide_thickness_357]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  357
  []
  [element_oxide_thickness_361]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  361
  []
  [element_oxide_thickness_365]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  365
  []
  [element_oxide_thickness_369]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  369
  []
  [element_oxide_thickness_373]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  373
  []
  [element_oxide_thickness_377]
    type = ElementalVariableValue
    variable =  oxide_thickness
    elementid   =  377
  []
[]
[Postprocessors]
  [element_vonmises_stress_353]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  353
  []
  [element_vonmises_stress_357]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  357
  []
  [element_vonmises_stress_361]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  361
  []
  [element_vonmises_stress_365]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  365
  []
  [element_vonmises_stress_369]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  369
  []
  [element_vonmises_stress_373]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  373
  []
  [element_vonmises_stress_377]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  377
  []
[]

[Postprocessors]
  [element_vonmises_stress_0]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  0
  []
  [element_vonmises_stress_8]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  8
  []
  [element_vonmises_stress_12]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  12
  []
  [element_vonmises_stress_15]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  15
  []
  [element_vonmises_stress_23]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  23
  []
  [element_vonmises_stress_565]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  565
  []
  [element_vonmises_stress_566]
    type = ElementalVariableValue
    variable =  vonmises_stress
    elementid   =  566
  []
[]

[Postprocessors]
  [element_hydrostatic_stress_0]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  0
  []
  [element_hydrostatic_stress_8]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  8
  []
  [element_hydrostatic_stress_12]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  12
  []
  [element_hydrostatic_stress_15]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  15
  []
  [element_hydrostatic_stress_23]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  23
  []
  [element_hydrostatic_stress_565]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  565
  []
  [element_hydrostatic_stress_566]
    type = ElementalVariableValue
    variable =  hydrostatic_stress
    elementid   =  566
  []
[]
[Postprocessors]  
  [element_stress_zz_353]
    type = ElementalVariableValue
    variable =  stress_zz
    elementid   =  353
  []  
  [element_output_heat_flux_353]
    type = ElementalVariableValue
    variable =  output_heat_flux
    elementid   =  353
  []
  [element_critical_heat_flux_353]
    type = ElementalVariableValue
    variable =  critical_heat_flux
    elementid   =  353
  []
  [element_dnbr_353]
    type = ElementalVariableValue
    variable =  dnbr
    elementid   =  353
  []
[]
[Postprocessors]  
  [element_stress_xx_377]
    type = ElementalVariableValue
    variable =  stress_xx
    elementid   =  377
  []
  [element_output_heat_flux_377]
    type = ElementalVariableValue
    variable =  output_heat_flux
    elementid   =  377
  []
  [element_critical_heat_flux_377]
    type = ElementalVariableValue
    variable =  critical_heat_flux
    elementid   =  377
  []
  [element_dnbr_377]
    type = ElementalVariableValue
    variable =  dnbr
    elementid   =  377
  []
[]
