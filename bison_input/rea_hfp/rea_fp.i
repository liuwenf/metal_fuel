#==========================================================================================
#  3-D case
#
#   Equivalent inner radius =  0.880 mm
#   Total Fuel Area =  71.10 mm^2
#   Fuel Area =  50.10 mm^2
#   Fuel Area (excluding displacer) =  47.67 mm^2 
#   Equivalent outer fuel radius  =  3.993 mm
#   Equivalent outer radius  =  4.757 mm
#   Equivalent outer diameter  =  9.515 mm
#
# Boundaries:
#
#  sideset 1, x-axis at the bottom 
#          2, y-axis at the bottom
#          3, outer surrounding surface 
#          4, bottom surface
#          5, centerline
#          6, top surface (top)
#          7, x-axis at the top surface
#          8, y-axis at the top surface
#
#==========================================================================================

axial_offset = 2.379
fuel_inner_radius = 0.88e-3
fuel_outer_radius = 3.993e-3
rod_diameter = 9.515e-3
clad_outer_radius = 4.757e-3
fast_flux_factor = 1.1
peak_power_factor = 2.1
start_time = -20.0
end_time = 20.0

top_plus = 0.2001
top_minus = 0.1999

[GlobalParams]
  temperature = temp
  displacements = 'disp_x disp_y disp_z'
  order = FIRST
  family = LAGRANGE
  energy_per_fission = 3.2e-11  # J/fission
[]

[Mesh]
 [file]
   file = 3d_cm_n1_s4.e 
   type = FileMeshGenerator   
  []
  [./x_axis]
   type = BoundingBoxNodeSetGenerator
   input = file
   new_boundary = 1
   top_right = '0.0064 0.00001 0.00001'
   bottom_left = '-0.0064 -0.00001 -0.00001' 
 [../]
 [./y_axis]
   type = BoundingBoxNodeSetGenerator
   input = x_axis
   new_boundary = 2
   bottom_left  = '-0.00001 -0.0064 -0.00001'
   top_right = '0.00001 0.0064 0.00001' 
 [../]
  [./center]
   type = BoundingBoxNodeSetGenerator
   input = y_axis
   new_boundary = 5
   top_right = '-0.00001 -0.00001 -0.0001'
   bottom_left = '0.00001 0.00001 4.0' 
 [../]
  [./bottom]
   type = BoundingBoxNodeSetGenerator
   input = center
   new_boundary = 4
   top_right = '-0.0064 -0.0064 -0.00001'
   bottom_left = '0.0064 0.0064 0.00001' 
 [../]
 [./top]
   type = BoundingBoxNodeSetGenerator
   input = bottom
   new_boundary = 6
   top_right = '-0.0064 -0.0064 ${top_minus}'
   bottom_left = '0.0064 0.0064 ${top_plus}' 
 [../]
 [./top_x]
   type = BoundingBoxNodeSetGenerator
   input = top
   new_boundary = 7
   top_right = '0.0064 0.00001 ${top_plus}'
   bottom_left = '-0.0064 -0.00001 ${top_minus}' 
 [../]
 [./top_y]
   type = BoundingBoxNodeSetGenerator
   input = top_x
   new_boundary = 8
   top_right = '0.00001 0.0064 ${top_plus}' 
   bottom_left = '-0.00001 -0.0064 ${top_minus}'
 [../] 
 [./top_nodes]
   type = BoundingBoxNodeSetGenerator
   input = top_y
   new_boundary = top_nodes
   top_right = '0.0002 0.0064 4'
   bottom_left = '-0.0002 0.0062 0' 
 [../]
 [./bottom_nodes]
   type = BoundingBoxNodeSetGenerator
   input = top_nodes
   new_boundary = bottom_nodes
   top_right = '0.0002 -0.0064 4'
   bottom_left = '-0.0002 -0.0062 0' 
 [../]
 [./right_nodes]
   type = BoundingBoxNodeSetGenerator
   input = bottom_nodes
   new_boundary = right_nodes
   top_right = '0.0064 0.0002 4'
   bottom_left = '0.0062 -0.0002 0' 
  [../]
  [./left_nodes]
   type = BoundingBoxNodeSetGenerator
   input = right_nodes
   new_boundary = left_nodes
   top_right = '-0.0064 0.0002 4'
   bottom_left = '-0.0062 -0.0002 0' 
 [../]
 uniform_refine = 0
[]
[Problem]
  type = ReferenceResidualProblem
  reference_vector = 'ref'
  extra_tag_vectors = 'ref'
[]
[Variables]
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
  [temp]
    initial_condition = 300.0     # set initial temp to ambient
  []
[]

[AuxVariables]
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
[]

[AuxKernels]  
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
    axial_power_profile = const_peaking_factors
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
[]

[Functions]
  [f1]
    type = ConstantFunction
    value = 1
  []
  [const_peaking_factors]
    type = ConstantFunction
    value = ${fast_flux_factor}
  []
  [power_history]
    type = PiecewiseLinear
    data_file = power_history.csv
    format = columns
    scale_factor = ${peak_power_factor}
  []
  [axial_peaking_factors]
    type = PiecewiseBilinear
    data_file = axial_power_profile.csv    
    axis = 2
  []
  [coolant_temp]
    type = PiecewiseBilinear
    data_file = rea_hfp_coolant_temp.csv    
    axis = 2
  [] 
  [pressure_ramp] 
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
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz hydrostatic_stress strain_xx strain_yy strain_zz'
    decomposition_method = EigenSolution
    extra_vector_tags = 'ref'
  []
  [clad]
    block = '2 3'
    strain = FINITE
    eigenstrain_names = 'clad_thermal_eigenstrain'
    generate_output = 'vonmises_stress stress_xx stress_yy stress_zz hydrostatic_stress strain_xx strain_yy strain_zz'
    decomposition_method = EigenSolution
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
  [no_y_all] # pin pellets and clad along axis of symmetry (y)
    type = DirichletBC
    variable = disp_y
    boundary = '1 7'
    value = 0.0
  []
  [no_x_all] # pin pellets and clad along axis of symmetry (x)
    type = DirichletBC
    variable = disp_x
    boundary = '2 8'
    value = 0.0
  []
  [pin_top_bottom]  
    enable = false
    type = DirichletBC
    variable = disp_y
    boundary = 'top_nodes bottom_nodes'
    value = 0.0
  []
  [pin_left_right]
    enable = false
    type = DirichletBC
    variable = disp_x
    boundary = 'left_nodes right_nodes'
    value = 0.0
  []
  [no_x_centerline] # pin pellets and clad along axis of symmetry (x)
    type = DirichletBC
    variable = disp_x
    boundary = '5'
    value = 0.0
  []
  [no_y_centerline] # pin pellets and clad along axis of symmetry (x)
    type = DirichletBC
    variable = disp_y
    boundary = '5'
    value = 0.0
  []
  [no_z_all] # pin pellets and clad along axis of symmetry (x)
    type = DirichletBC
    variable = disp_z
    boundary = 4
    value = 0.0
  []
  [Pressure] #  apply coolant pressure on clad outer walls
    [coolantPressure]
      boundary = '3'
      factor = 15.5e6
      function = pressure_ramp   # use the pressure_ramp function defined above
    []
  []  
[]

[Constraints]
  [disp_z]
    type = EqualValueBoundaryConstraint
    variable = disp_z
    secondary = '6'
    penalty = 1e7
  []
[]

[CoolantChannel]
  [convective_clad_surface] # apply convective boundary to clad outer surface
    variable = temp
    boundary = '3'
    inlet_temperature = coolant_temp # K
    inlet_pressure = 15.5e6   # Pa
    inlet_massflux = 3600.0   # kg/m^2-sec
    rod_diameter = ${rod_diameter}  # m
    rod_pitch = 1.26e-2       # m
    linear_heat_rate  = power_history
    axial_power_profile = axial_peaking_factors
    axial_offset = ${axial_offset}    
    flow_direction = z
    outputs = all
    output_properties = 'coolant_ehtnalpy coolant_channel_htc  coolant_channel_hmode  coolant_channel_htype  coolant_temperature output_heat_flux critical_heat_flux dnbr'
    compute_enthalpy = false
    model_post_chf = true
    chf_scalef        = 1.0
    htc_scalef        = 1.0
    htc_correlation_type = 8
    chf_correlation_type = 4
  # mode_to_apply_scaling = 6    
  # model_radiation = true
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
    pellet_radius       = ${fuel_outer_radius}
    pellet_inner_radius = ${fuel_inner_radius}
    axial_offset = ${axial_offset}
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
    type = ZryCreepHayesHoppeUpdate
    block = '2 3'
    fast_neutron_flux = fast_neutron_flux
    temperature = temp
    zircaloy_material_type = stress_relief_annealed
    model_irradiation_creep = true
    model_thermal_creep = true
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
     clad_inner_radius = ${fuel_outer_radius} 
     clad_outer_radius = ${clad_outer_radius}
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


  dtmin = 1.0e-7
  dtmax = 1.0
  start_time = ${start_time}
  end_time   = ${end_time}

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    optimal_iterations = 200
    linear_iteration_ratio = 100
    timestep_limiting_function = power_history
    max_function_change = 10000
    force_step_every_function_point = true
  [../]
  [Quadrature]
    order = THIRD
  []
[]

[Postprocessors]
  # Define postprocessors (some are required as specified above; others are optional; many others are available)
  
  [flux_from_clad]           # area integrated heat flux from the cladding
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
  [alive_time]
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
  [peak_fuel_temperature]
    type = NodalMaxValue
    variable = temp
    block = '1 3'
  []
  [peak_clad_temperature]
    type = NodalMaxValue
    variable = temp
    block = '2'
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
