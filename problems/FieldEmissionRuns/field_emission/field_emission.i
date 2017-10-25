gap = 10E-6 #m
appliedVoltage = 0.06 #kV
steadyStateTime = 1 # s

resistance = 1.09E0 #Ohms
area = 1E-4 # Formerly 3.14e-6

position_units = 1 #m
time_units = 1E-9 #s

NX = 1500
dom0Size = ${/ ${gap} ${position_units}}

ssT = ${/ ${steadyStateTime} ${time_units}}

[GlobalParams]
        offset = 20
		time_units = ${time_units}
		position_units = ${position_units}
        potential_units = kV
        use_moles = true
[]

[Mesh]
		type = GeneratedMesh	# Can generate simple lines, rectangles and rectangular prisms
		dim = 1								# Dimension of the mesh
		nx = ${NX}							# Number of elements in the x direction
		xmax = ${dom0Size}				# Length of test chamber
[]


[Problem]
        type = FEProblem
        # kernel_coverage_check = false
[]

[Preconditioning]
        [./smp]
                type = SMP
                full = true
        [../]
[]

[Executioner]
	type = Transient
	end_time = 1E6
		
	trans_ss_check = 1
	ss_check_tol = 1E-15
	ss_tmin = ${ssT}
		
		petsc_options = '-snes_ksp_ew -superlu_dist' # -snes_converged_reason -snes_linesearch_monitor'
		solve_type = NEWTON
		
		petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
		petsc_options_value = 'lu superlu_dist'
		
        nl_rel_tol = 1E-11
        nl_abs_tol = 1E-11
        dtmin = ${/ 1e-15 ${time_units}}
        #dtmax = ${/ 1e-9 ${time_units}}
        [./TimeStepper]
                type = IterationAdaptiveDT
                cutback_factor = 0.4
                dt = ${/ 1e-12 ${time_units}}
                growth_factor = 1.2
                optimal_iterations = 15
        [../]
[]

[Outputs]
        print_perf_log = true
        print_linear_residuals = false
        [./out]
                type = Exodus
                execute_on = 'final'
        [../]
[]

[Debug]
        show_var_residual_norms = false # true
[]

[Postprocessors]
#        [./current_cathode]
#                type = SideAverageValue
#                execute_on = timestep_end
#				boundary = left
#                variable = tot_gas_current
#        [../]
[]

[UserObjects]
 	[./current_density_user_object]
 		type = CurrentDensityShapeSideUserObject
 		boundary = left
 		potential = potential
 		em = em
 		ip = Arp
 		mean_en = mean_en
 		execute_on = 'linear nonlinear'
 	[../]
        [./data_provider]
                type = ProvideMobility
                electrode_area = ${area}
				ballast_resist = ${resistance}
                e = 1.6e-19
        [../]
[]

[Kernels]
        [./em_time_deriv]
                type = ElectronTimeDerivative
                variable = em
                block = 0
        [../]
        [./em_advection]
                type = EFieldAdvectionElectrons
                variable = em
                potential = potential
                mean_en = mean_en
                block = 0
        [../]
        [./em_diffusion]
                type = CoeffDiffusionElectrons
                variable = em
                mean_en = mean_en
                block = 0
        [../]
        [./em_ionization]
                type = ElectronsFromIonization
                em = em
                variable = em
                potential = potential
                mean_en = mean_en
                block = 0
        [../]
        [./em_log_stabilization]
                type = LogStabilizationMoles
                variable = em
                block = 0
        [../]
        # [./em_advection_stabilization]
        #               type = EFieldArtDiff
        #               variable = em
        #               potential = potential
        #               block = 0
        # [../]

        [./potential_diffusion_dom1]
                type = CoeffDiffusionLin
                variable = potential
                block = 0
        [../]

        [./Arp_charge_source]
                type = ChargeSourceMoles_KV
                variable = potential
                charged = Arp
                block = 0
        [../]
        [./em_charge_source]
                type = ChargeSourceMoles_KV
                variable = potential
                charged = em
                block = 0
        [../]

        [./Arp_time_deriv]
                type = ElectronTimeDerivative
                variable = Arp
                block = 0
        [../]
        [./Arp_advection]
                type = EFieldAdvection
                variable = Arp
                potential = potential
                block = 0
        [../]
        [./Arp_diffusion]
                type = CoeffDiffusion
                variable = Arp
                block = 0
        [../]
        [./Arp_ionization]
                type = IonsFromIonization
                variable = Arp
                potential = potential
                em = em
                mean_en = mean_en
                block = 0
        [../]
        [./Arp_log_stabilization]
                type = LogStabilizationMoles
                variable = Arp
                block = 0
        [../]
        # [./Arp_advection_stabilization]
        #               type = EFieldArtDiff
        #               variable = Arp
        #               potential = potential
        #               block = 0
        # [../]


        [./mean_en_time_deriv]
                type = ElectronTimeDerivative
                variable = mean_en
                block = 0
        [../]
        [./mean_en_advection]
                type = EFieldAdvectionEnergy
                variable = mean_en
                potential = potential
                em = em
                block = 0
        [../]
        [./mean_en_diffusion]
                type = CoeffDiffusionEnergy
                variable = mean_en
                em = em
                block = 0
        [../]
        [./mean_en_joule_heating]
                type = JouleHeating
                variable = mean_en
                potential = potential
                em = em
                block = 0
        [../]
        [./mean_en_ionization]
                type = ElectronEnergyLossFromIonization
                variable = mean_en
                potential = potential
                em = em
                block = 0
        [../]
        [./mean_en_elastic]
                type = ElectronEnergyLossFromElastic
                variable = mean_en
                potential = potential
                em = em
                block = 0
        [../]
        [./mean_en_excitation]
                type = ElectronEnergyLossFromExcitation
                variable = mean_en
                potential = potential
                em = em
                block = 0
        [../]
        [./mean_en_log_stabilization]
                type = LogStabilizationMoles
                variable = mean_en
                block = 0
                offset = 15
        [../]
        # [./mean_en_advection_stabilization]
        #               type = EFieldArtDiff
        #               variable = mean_en
        #               potential = potential
        #               block = 0
        # [../]
[]

[Variables]
	[./potential]
	[../]
	[./em]
		block = 0
		initial_condition = -20 
		#15
	[../]
	[./Arp]
		block = 0
		initial_condition = -20
		#15
	[../]
	[./mean_en]
		block = 0
		initial_condition = -20
		#14
	[../]

[]

[AuxVariables]
        [./e_temp]
                block = 0
                order = CONSTANT
                family = MONOMIAL
        [../]
        [./x]
                order = CONSTANT
                family = MONOMIAL
        [../]
        [./x_node]
        [../]
        [./rho]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./em_lin]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./Arp_lin]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./Efield]
                order = CONSTANT
                family = MONOMIAL
        [../]
        [./Current_em]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./Current_Arp]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./tot_gas_current]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./EFieldAdvAux_em]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./DiffusiveFlux_em]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./PowerDep_em]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./PowerDep_Arp]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./ProcRate_el]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./ProcRate_ex]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
        [./ProcRate_iz]
                order = CONSTANT
                family = MONOMIAL
                block = 0
        [../]
[]

[AuxKernels]
        [./PowerDep_em]
                type = PowerDep
                density_log = em
                potential = potential
                art_diff = false
                potential_units = kV
                variable = PowerDep_em
                block = 0
        [../]
        [./PowerDep_Arp]
                type = PowerDep
                density_log = Arp
                potential = potential
                art_diff = false
                potential_units = kV
                variable = PowerDep_Arp
                block = 0
        [../]
        [./ProcRate_el]
                type = ProcRate
                em = em
                potential = potential
                proc = el
                variable = ProcRate_el
                block = 0
        [../]
        [./ProcRate_ex]
                type = ProcRate
                em = em
                potential = potential
                proc = ex
                variable = ProcRate_ex
                block = 0
        [../]
        [./ProcRate_iz]
                type = ProcRate
                em = em
                potential = potential
                proc = iz
                variable = ProcRate_iz
                block = 0
        [../]
        [./e_temp]
                type = ElectronTemperature
                variable = e_temp
                electron_density = em
                mean_en = mean_en
                block = 0
        [../]
        [./x_g]
                type = Position
                variable = x
                block = 0
        [../]
        [./x_ng]
                type = Position
                variable = x_node
                block = 0
        [../]
        [./rho]
                type = ParsedAux
                variable = rho
                args = 'em_lin Arp_lin'
                function = 'Arp_lin - em_lin'
                execute_on = 'timestep_end'
                block = 0
        [../]
        [./tot_gas_current]
                type = ParsedAux
                variable = tot_gas_current
                args = 'Current_em Current_Arp'
                function = 'Current_em + Current_Arp'
                execute_on = 'timestep_end'
                block = 0
        [../]
        [./em_lin]
                type = Density
#               convert_moles = true
                variable = em_lin
                density_log = em
                block = 0
        [../]
        [./Arp_lin]
                type = Density
#               convert_moles = true
                variable = Arp_lin
                density_log = Arp
                block = 0
        [../]
        [./Efield_g]
                type = Efield
                component = 0
                potential = potential
                variable = Efield
                block = 0
        [../]
        [./Current_em]
                type = Current
                potential = potential
                density_log = em
                variable = Current_em
                art_diff = false
                block = 0
        [../]
        [./Current_Arp]
                type = Current
                potential = potential
                density_log = Arp
                variable = Current_Arp
                art_diff = false
                block = 0
        [../]
        [./EFieldAdvAux_em]
                type = EFieldAdvAux
                potential = potential
                density_log = em
                variable = EFieldAdvAux_em
                block = 0
        [../]
        [./DiffusiveFlux_em]
                type = DiffusiveFlux
                density_log = em
                variable = DiffusiveFlux_em
                block = 0
        [../]
[]

[BCs]

	[./potential_dirichlet_left]
		type = DirichletBC
		variable = potential
		boundary = left
		value = -${appliedVoltage}
	[../]

#	[./potential_neumann_left]
#		type = NeumannCircuitPotentialNew
#		variable = potential
#		boundary = left
#		ip = Arp
#		data_provider = data_provider
#		em = em
#		mean_en = mean_en
#		current = current_cathode
#		surface_potential = potential_bc_func
#		surface = cathode
#		resist = ${resistance}
#		position_units = ${position_units}
#		use_moles = false
#		A = ${area}
#		r = 0
#	[../]
# 	[./potential_left]
# 		boundary = left
# 		type = PenaltyCircuitPotential
# 		surface_potential = -${appliedVoltage}
# 		penalty = 1
# 		variable = potential
# 		current = current_density_user_object
# 		surface = 'cathode'
# 		data_provider = data_provider
# 		em = em
# 		ip = Arp
# 		mean_en = mean_en
# 		area = ${area}
# 		resistance = ${resistance}
# 	[../]
	[./potential_dirichlet_right]
		type = DirichletBC
		variable = potential
		boundary = right
		value = 0
	[../]


        [./em_physical_left]
                type = HagelaarElectronBC
                variable = em
                boundary = 'left'
                potential = potential
                mean_en = mean_en
                r = 0
        [../]
       [./sec_electrons_left]
               type = SecondaryElectronBC
               variable = em
               boundary = 'left'
               potential = potential
               ip = Arp
               mean_en = mean_en
               r = 1
       [../]
#        [./FieldEmission_left]
#                type = FieldEmissionBC
#                variable = em
#                boundary = 'left'
#                potential = potential
#                ip = Arp
#                mean_en = mean_en
#                r = 1
#       [../]
		
        [./em_physical_right]
                type = HagelaarElectronBC
                variable = em
                boundary = right
                potential = potential
                mean_en = mean_en
                r = 0.99
        [../]


        [./Arp_physical_left_diffusion]
                type = HagelaarIonDiffusionBC
                variable = Arp
                boundary = 'left'
                r = 0
        [../]
        [./Arp_physical_left_advection]
                type = HagelaarIonAdvectionBC
                variable = Arp
                boundary = 'left'
                potential = potential
                r = 0
        [../]
        [./Arp_physical_right_diffusion]
                type = HagelaarIonDiffusionBC
                variable = Arp
                boundary = right
                r = 0
        [../]
        [./Arp_physical_right_advection]
                type = HagelaarIonAdvectionBC
                variable = Arp
                boundary = right
                potential = potential
                r = 0
        [../]
		
        [./mean_en_physical_left]
                type = HagelaarEnergyBC
                variable = mean_en
                boundary = 'left'
                potential = potential
                em = em
                ip = Arp
                r = 0
        [../]
#		[./mean_en_sec_left]
#                type = SecondaryEnergyBC
#               variable = mean_en
#                boundary = 'left'
#                potential = potential
#                em = em
#                ip = Arp
#                r = 1
#        [../]
#		[./mean_en_field_left]
#                type = FieldEnergyBC
#                variable = mean_en
#                boundary = 'left'
#                potential = potential
#                em = em
#                ip = Arp
#                r = 1
#       [../]

        [./mean_en_physical_right]
                type = HagelaarEnergyBC
                variable = mean_en
                boundary = right
                potential = potential
                em = em
                ip = Arp
                r = 0.99
        [../]

[]

[ICs]
        [./potential_ic]
                type = FunctionIC
                variable = potential
                function = potential_ic_func
        [../]
		
        # [./em_ic]
        #               type = RandomIC
        #               variable = em
        #               block = 0
        # [../]
        # [./Arp_ic]
        #               type = RandomIC
        #               variable = Arp
        #               block = 0
        # [../]
        # [./mean_en_ic]
        #               type = RandomIC
        #               variable = mean_en
        #               block = 0
        # [../]
        # [./potential_ic]
        #               type = RandomIC
        #               variable = potential
        # [../]
[]

[Functions]
        [./potential_bc_func]
		type = ParsedFunction
		value = -${appliedVoltage}
        [../]
        [./potential_ic_func]
                type = ParsedFunction
                value = '-${appliedVoltage} * (${dom0Size} - x) / ${dom0Size}'
        [../]
        [./cathode_temperature]
                type = ParsedFunction
                value = 300
        [../]
[]

[Materials]
        [./gas_block]
                type = Gas
                interp_trans_coeffs = true
                interp_elastic_coeff = true
                ramp_trans_coeffs = false
                em = em
                potential = potential
                ip = Arp
                mean_en = mean_en
                user_se_coeff = 0.1
                user_work_function = 4.5 # eV
                user_field_enhancement = 55 #8.4
                property_tables_file = td_argon_mean_en.txt
                block = 0
        [../]
[]
