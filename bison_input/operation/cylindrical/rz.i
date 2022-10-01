#====================================================================================================
# A multi-slice model for helical metallic-fuel (Generalized Plane Strain Model)
# Fuel design parameters:
#
#   Equivalent inner radius =  0.880 mm
#   Total Fuel Area =  71.10 mm^2
#   Fuel Area =  50.10 mm^2
#   Fuel Area (excluding displacer) =  47.67 mm^2 
#   Equivalent outer fuel radius  =  3.993 mm
#   Equivalent outer radius  =  4.757 mm
#   Equivalent outer diameter  =  9.515 mm
#
#   Displacer width:   1.56 mm
#   Rod-to-rod pitch: 12.6 mm
#
#   Rod length:        3.66 m 
#
# Initial Temperature: 300 K
#
#  Element type: QUAD-8
#
#
#====================================================================================================
axial_offset = 0.0
axial_scale_factor = 1.0

end_time = 1.8e8
#end_time = 86400
num_steps = 3000

solid_swelling_scale_factor = 1.0
gas_swelling_scale_factor = 0.0

[GlobalParams]
  temperature = temp
  displacements = 'disp_x disp_y'
  order = SECOND
  family = LAGRANGE
  energy_per_fission = 3.2e-11  # J/fission
[]

[Problem]
  coord_type = RZ
  type = ReferenceResidualProblem
  reference_vector = 'ref'
  extra_tag_vectors = 'ref'
  group_variables = 'disp_x disp_y'
[]

[Mesh]
  [left]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 12
    ny = 200
    xmax = 0.00393
    ymax = 3.66
    elem_type = QUAD8
  []
  [right]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 200
    xmin = 0.00393
    xmax = 0.004757
    ymax = 3.66
    elem_type = QUAD8
  []
  [smg]
    type = StitchedMeshGenerator
    inputs = 'left right'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'right left'
    #parallel_type = 'replicated'
  []
  [fuel]
    type = SubdomainBoundingBoxGenerator
    input = smg
    bottom_left = '0 0 0'
    block_id = 1
    top_right = '0.00393 3.66 0'
  []
  [clad]
    type = SubdomainBoundingBoxGenerator
    input = fuel
    bottom_left = '0.00393 0 0'
    block_id = 2
    top_right = '0.004757 3.66 0'
  []
[]

[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [temp]
    initial_condition = 300.0     # set initial temp to ambient
  []  
[]

[AuxVariables]
  [fast_neutron_flux]
    block = '2'
  []
  [fast_neutron_fluence]
    block = '2'
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
  [oxide]
    type = MaterialRealAux
    property = oxide_scale_thickness
    variable = oxide_thickness
    boundary = 'right'
  []
  [fast_neutron_flux]
    type = FastNeutronFluxAux
    variable = fast_neutron_flux
    block = '2'
    rod_ave_lin_pow = power_history
    axial_power_profile = axial_peaking_factors
    factor = 3e13
    execute_on = timestep_begin
  []
  [fast_neutron_fluence]
    type = FastNeutronFluenceAux
    variable = fast_neutron_fluence
    block = '2'
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
  [f1]
    type = ConstantFunction
    value = 1
  []
  [const_power_history]
    type = PiecewiseLinear
    x = '0 10800.0'
    y = '0 20000.0'
  []
  [const_peaking_factors]
    type = ConstantFunction
    value = 1
  []  
  [power_history]
    type = PiecewiseLinear
    data_file = power_history.csv
    format = columns
  []
  [axial_peaking_factors]
    type = PiecewiseBilinear
    data_file = axial_power_profile.csv
    axis = 1
  []  
  [pressure_ramp]              # reads and interpolates input data defining amplitude curve for fill gas pressure
    type = PiecewiseLinear
    x = '0 10800'
    y = '0 1'
  []
[]

[Modules/TensorMechanics/Master]
  [fuel]
    block = '1'
    strain = FINITE
    eigenstrain_names = 'fuel_thermal_eigenstrain swelling'
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz hydrostatic_stress'
    # decomposition_method = EigenSolution
    extra_vector_tags = 'ref'
  []
  [clad]
    block = '2'
    strain = FINITE
    eigenstrain_names = 'clad_thermal_eigenstrain'
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz hydrostatic_stress'
    # decomposition_method = EigenSolution
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
    fission_rate = fission_rate
    extra_vector_tags = 'ref'    
  []
[]


[BCs]
# Define boundary conditions
  [no_y_all] # x-axis
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  []
  [no_x_all] # y-axis
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  []  
  [Pressure] #  apply coolant pressure on clad outer walls
    [coolantPressure]
      boundary = 'right'
      factor = 15.5e6
      function = pressure_ramp   # use the pressure_ramp function defined above
    []
  []  
[]
[CoolantChannel]
  [convective_clad_surface] # apply convective boundary to clad outer surface
    variable = temp
    boundary = 'right'
    inlet_temperature = 558   # K
    inlet_pressure = 15.5e6   # Pa
    inlet_massflux = 3600.0   # kg/m^2-sec
    rod_diameter = 0.9515e-2  # m
    rod_pitch = 1.26e-2  # m
    linear_heat_rate  = power_history
    axial_power_profile = axial_peaking_factors
    outputs = all
    output_properties = 'coolant_ehtnalpy coolant_channel_htc  coolant_channel_hmode  coolant_channel_htype    coolant_temperature output_heat_flux critical_heat_flux dnbr'
    oxide_thickness = oxide_thickness 
#   flow_direction =  z
    axial_scale_factor = ${axial_scale_factor}
    axial_offset = ${axial_offset}
 []
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
    pellet_radius = 3.993e-3
    outputs = all
    output_properties = fission_rate
    axial_scale_factor = ${axial_scale_factor}
    axial_offset = ${axial_offset}
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
     solid_swelling_scale_factor = ${solid_swelling_scale_factor}
     gas_swelling_scale_factor = ${gas_swelling_scale_factor}
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
    creep_model = bmi
    max_inelastic_increment = 1e-2
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
    block = '2'
    temperature = temp
  []
  [clad_elasticity_tensor]
    type = ZryElasticityTensor
    block = '2'
  []
  [clad_creep_model]
    type = ZryCreepLimbackHoppeUpdate
    block = '2'
    fast_neutron_flux = fast_neutron_flux
    temperature = temp
    zircaloy_material_type = stress_relief_annealed
    model_irradiation_creep = true
    model_primary_creep = true
    model_thermal_creep = true
  []
  [clad_stress]
    type = ComputeMultipleInelasticStress
    block = '2'
    tangent_operator = elastic
    inelastic_models = 'clad_creep_model'
  []
  [clad_thermal_expansion]
     type = ZryThermalExpansionMATPROEigenstrain
     block = '2'
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
     boundary = 'right'   
  []
  [clad_density]
    type = Density
    block = '2'
    density = 6551.0
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

  start_time = 0.0
  end_time = ${end_time}
  num_steps = ${num_steps}
  
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 2.0e2
    time_t  = '1e4 1e5 1e6'
    time_dt = '1e3 1e4 1e5'
  []
  dtmax = 2e6
  dtmin = 1
  
#  optimal_iterations = 6
#  iteration_window = 2
#  linear_iteration_ratio = 100

  [Quadrature]
    order = THIRD
  []
[]

[Postprocessors]
  # Define postprocessors (some are required as specified above; others are optional; many others are available)  
  [flux_from_clad]           # area integrated heat flux from the cladding
    type = SideDiffusiveFluxIntegral
    variable = temp
    boundary = 'right'
    diffusivity = thermal_conductivity
    execute_on = timestep_end
  []   
  [_dt]                     # time step
    type = TimestepSize
    execute_on = timestep_end
  []
  [num_lin_it]
    type = NumLinearIterations
  []
  [num_nonlin_it]
    type = NumNonlinearIterations
  []
  [tot_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  []
  [tot_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
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
  [peak_heatflux]
    type = ElementExtremeValue
    variable = output_heat_flux
    block = '2'
  []
  [peak_oxide_thickness]
    type = ElementExtremeValue
    variable = oxide_thickness
    block = '2'
  []
  [peak_fuel_temperature]
    type = NodalMaxValue
    variable = temp
    block = '1'
  []
  [total_volume]
    type = FunctionElementIntegral
    use_displaced_mesh = true
    function = f1
  []
  [fuel_volume]
    type = FunctionElementIntegral
    use_displaced_mesh = true
    function = f1
    block = '1'
  []
  [clad_volume]
    type = FunctionElementIntegral
    use_displaced_mesh = true
    function = f1
    block = '2'
  []
[]

#[VectorPostprocessors]
#  [axial_temp]
#    type = LineValueSampler
#    end_point = '0.0 0.0 2.0'
#    start_point = '0.0 0.0 0.0'
#    sort_by = z
#    variable = temp
#    execute_on = timestep_end
#    num_points = 21
#    outputs = 'outfile_axial_temp'
#  []
#[]

[Outputs]
  perf_graph = true
  exodus = true
  csv = true
  [console]
    type = Console
    max_rows = 25
  []
 # [outfile_axial_temp]
 #   type = CSV
 #   execute_on = 'FINAL'
 # []
[]
