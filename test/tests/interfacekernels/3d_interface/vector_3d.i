[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  xmax = 2
  ny = 2
  ymax = 2
  nz = 2
  zmax = 2
  elem_type = HEX20
[]

[MeshModifiers]
  [./subdomain1]
    type = SubdomainBoundingBox
    bottom_left = '0 0 0'
    top_right = '1 1 1'
    block_id = 1
  [../]
  [./break_boundary]
    type = BreakBoundaryOnSubdomain
    depends_on = subdomain1
  [../]
  [./interface]
    type = SideSetsBetweenSubdomains
    depends_on = break_boundary
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = NEDELEC_ONE
    block = 0
  [../]

  [./v]
    order = FIRST
    family = NEDELEC_ONE
    block = 1
  [../]
[]

[Kernels]
  [./curl_u_plus_u]
    type = VectorFEWave
    variable = u
    block = 0
  [../]
  [./curl_v_plus_v]
    type = VectorFEWave
    variable = v
    block = 1
  [../]
  [./u_source]
    type = VectorBodyForce
    variable = u
    function_x = 1
    function_y = 1
    function_z = 1
  [../]
[]

[InterfaceKernels]
  [./parallel]
    type = VectorPenaltyInterfaceDiffusion
    variable = u
    neighbor_var = v
    boundary = master0_interface
    penalty = 1e6
  [../]
[]

[BCs]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
[]
