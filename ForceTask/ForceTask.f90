program ForceTask
	use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
	use ConstantSeedInitialize
    use SynapsesFactoryModule
    use jointAnkleForceTaskClass
    use AfferentPoolClass
    implicit none 
	
	type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength, L_soma_SOL, L_terminal_LG, L_terminal_MG, L_terminal_SOL, L_afferent_SOL
    integer :: i, j, poolIndex
	real(wp), dimension(:), allocatable :: t, MNv_mV_SOL,MNv_mV_MG,MNv_mV_LG, torque, forceLG, forceSOL, forceMG, emgSOL, emgMG, emgLG
    real(wp), dimension(:), allocatable :: tspks_LG_terminal, tspks_MG_terminal, tspks_SOL_soma, tspks_SOL_terminal, condvel
	real(wp), dimension(:), allocatable :: indiceMN_LG_terminal, indiceMN_MG_terminal, indiceMN_SOL_soma, indiceMN_SOL_terminal
	real(wp), dimension(:), allocatable :: IaFRSOL, IaFRMG, IaFRLG, indiceIA_SOL, tspks_SOL_afferentIA, velocitySOL, length
	real(wp) :: tic, toc, FRbasal, t1
    type(gpf) :: gp
	real(wp) :: FR
    integer :: trial, trials, k, indice
    integer :: scenario, cenvel, Sconf, velconf
	integer :: GammaOrder
	character(len = 80) :: pool, muscle, fileNumber, group, ForceMVC, paramTag
	character(len = 80) :: conffile = 'confForceTask'
	character(len = 80) :: configurationFile, iterS, itervel
	type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools    
    type(jointAnkleForceTask) :: ankle
    real(wp) :: angle, STD, M, CV, velocity_m_s
	
	
	trial = 1
	trials = 10
		
	do while (trial <= trials)
	do Sconf = 1,3
		
		scenario = Sconf - 1
		write(iterS, '(I3)')scenario
			
		call init_random_seed()
	
			configurationFile = trim(conffile) //'_S'// trim(adjustl(iterS))//'.rmto'
			print*, configurationFile
			conf = Configuration(configurationFile)
	
	
			allocate(afferentPools(3))
			pool = 'Ia'
			muscle = 'SOL'
			afferentPools(1) = AfferentPool(conf, pool, muscle)

			pool = 'Ia'
			muscle = 'MG'
			afferentPools(2) = AfferentPool(conf, pool, muscle)

			pool = 'Ia'
			muscle = 'LG'
			afferentPools(3) = AfferentPool(conf, pool, muscle)
		
	
			allocate(neuralTractPools(1))
			pool = 'CMExt'
			neuralTractPools(1) = NeuralTract(conf, pool)
			
			allocate(motorUnitPools(3))
			pool = 'SOL'
			motorUnitPools(1) = MotorUnitPool(conf, pool) 
			
			pool = 'MG'
			motorUnitPools(2) = MotorUnitPool(conf, pool)
			
			pool = 'LG'
			motorUnitPools(3) = MotorUnitPool(conf, pool)
			
			ankle = jointAnkleForceTask(conf, motorUnitPools)
			
			allocate(interneuronPools(0))
				
			synapticNoisePools = synapseFactory(conf, neuralTractPools, motorUnitPools, interneuronPools, afferentPools)
			
			tf = conf%simDuration_ms
			dt = conf%timeStep_ms
			timeLength = nint(tf/dt)
			
			allocate(t(timeLength))
			allocate(MNv_mV_SOL(timeLength))
			allocate(MNv_mV_MG(timeLength))
			allocate(MNv_mV_LG(timeLength))
			allocate(torque(timeLength))
			allocate(emgSOL(timeLength))
			allocate(forceSOL(timeLength))
			allocate(emgMG(timeLength))
			allocate(forceMG(timeLength))
			allocate(emgLG(timeLength))
			allocate(forceLG(timeLength))
			allocate(IaFRSOL(timeLength))
			allocate(IaFRMG(timeLength))
			allocate(IaFRLG(timeLength))
			allocate(velocitySOL(timeLength))
			allocate(length(timeLength))
			
			t = [(dt*(i-1), i=1, timeLength)]
	
			t1 = 2000.0;
			FRbasal = 40.0	

				print*, 'Trial = ', trial 
				do j = 1, size(afferentPools)
					call afferentPools(j)%reset()
				end do
				
				do j = 1, size(neuralTractPools)
					call neuralTractPools(j)%reset()
				end do

				do j = 1, size(motorUnitPools)
					call motorUnitPools(j)%reset()
				end do
				
				do j = 1, size(interneuronPools)
					call interneuronPools(j)%reset()
				end do
			
				call ankle%reset()
				
				call cpu_time(tic)
				do i = 1, size(t)
			
					angle = 0
					call ankle%atualizeAnkle(t(i), angle) 
				
					if (i <= t1/dt) then
						FR = FRbasal
						GammaOrder = 10
					else if (i > t1/dt) then
						FR = 250.0
						GammaOrder = 1
					end if
				
				
					do j = 1, size(neuralTractPools)
						call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
					end do
					
					do j = 1, 3
						call motorUnitPools(j)%atualizeMotorUnitPool(t(i), 25.0_wp, 25.0_wp)
						MNv_mV_SOL(i) = motorUnitPools(1)%v_mV(30)
						MNv_mV_LG(i) = motorUnitPools(2)%v_mV(30)
						MNv_mV_MG(i) = motorUnitPools(3)%v_mV(30)
					end do
					
					do j = 1, 3
						call afferentPools(j)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IaFR_Hz)
					end do
					IaFRSOL(i) = motorUnitPools(1)%spindle%IaFR_Hz
			 		IaFRLG(i) = motorUnitPools(2)%spindle%IaFR_Hz
					IaFRMG(i) = motorUnitPools(3)%spindle%IaFR_Hz
					call ankle%computeTorque(t(i))

					
				end do 
				call cpu_time(toc)
				
				print '(F15.6, A)', toc - tic, ' seconds'
				print*, 'GammaOrder = ', GammaOrder
				print*, 'FR = ', FR
			
				call motorUnitPools(1)%listSpikes()
				call motorUnitPools(2)%listSpikes()
				call motorUnitPools(3)%listSpikes()
				call neuralTractPools(1)%listSpikes()
				call motorUnitPools(1)%getMotorUnitPoolEMG()
				call motorUnitPools(2)%getMotorUnitPoolEMG()
				call motorUnitPools(3)%getMotorUnitPoolEMG()
						
				do j = 1, size(afferentPools)
					call afferentPools(j)%listSpikes()
				end do
				
				do i = 1,timeLength
					torque(i) = ankle%ankleTorque_Nm(i)
					forceSOL(i) = motorUnitPools(1)%HillMuscle%force(i)
					forceMG(i) = motorUnitPools(2)%HillMuscle%force(i)
					forceLG(i) = motorUnitPools(3)%HillMuscle%force(i)
					emgSOL(i) = motorUnitPools(1)%emg(i)
					emgMG(i) = motorUnitPools(2)%emg(i)
					emgLG(i) = motorUnitPools(3)%emg(i)
					velocitySOL(i) = motorUnitPools(1)%HillMuscle%velocity_m_ms(i)
					length(i) = motorUnitPools(1)%HillMuscle%length_m(i)

					!print *, torque(i)
				end do
				
				L_soma_SOL = size(motorUnitPools(1)%poolSomaSpikes(:,1))
				L_terminal_SOL = size(motorUnitPools(1)%poolTerminalSpikes(:,1))
				L_terminal_MG = size(motorUnitPools(2)%poolTerminalSpikes(:,1))
				L_terminal_LG = size(motorUnitPools(3)%poolTerminalSpikes(:,1))
				L_afferent_SOL = size(afferentPools(1)%poolTerminalSpikes(:,1))
				
				allocate(tspks_SOL_soma(L_soma_SOL))
				allocate(tspks_SOL_terminal(L_terminal_SOL))
				allocate(tspks_MG_terminal(L_terminal_MG))
				allocate(tspks_LG_terminal(L_terminal_LG))
				allocate(tspks_SOL_afferentIA(L_afferent_SOL))
				allocate(indiceIA_SOL(L_afferent_SOL))
				allocate(indiceMN_SOL_soma(L_soma_SOL))
				allocate(indiceMN_SOL_terminal(L_terminal_SOL))
				allocate(indiceMN_MG_terminal(L_terminal_MG))
				allocate(indiceMN_LG_terminal(L_terminal_LG))
				
				
				do j = 1, L_soma_SOL
					tspks_SOL_soma(j) = motorUnitPools(1)%poolSomaSpikes(j,1)
					indiceMN_SOL_soma(j) = motorUnitPools(1)%poolSomaSpikes(j,2)
				end do
				
				do j = 1, L_terminal_SOL
					tspks_SOL_terminal(j) = motorUnitPools(1)%poolTerminalSpikes(j,1)
					indiceMN_SOL_terminal(j) = motorUnitPools(1)%poolTerminalSpikes(j,2)
				end do
				do j = 1, L_terminal_MG
					tspks_MG_terminal(j) = motorUnitPools(2)%poolTerminalSpikes(j,1)
					indiceMN_MG_terminal(j) = motorUnitPools(2)%poolTerminalSpikes(j,2)
				end do
				
				do j = 1, L_terminal_LG
					tspks_LG_terminal(j) = motorUnitPools(3)%poolTerminalSpikes(j,1)
					indiceMN_LG_terminal(j) = motorUnitPools(3)%poolTerminalSpikes(j,2)
				end do
				
				do j = 1, L_afferent_SOL
					tspks_SOL_afferentIA(j) = afferentPools(1)%poolTerminalSpikes(j,1)
					indiceIA_SOL(j) = afferentPools(1)%poolTerminalSpikes(j,2)
				end do
					
				write(fileNumber, '(I5.1)')trial
				filename =  'Torque_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(torque(i))
				end do
				
				filename =  'ForceSOL_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceSOL(i))
				end do
				
				filename =  'ForceMG_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceMG(i))
				end do
				
				filename =  'ForceLG_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceLG(i))
				end do
				
				filename =  'EMGSOL_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgSOL(i))
				end do
				
				filename =  'EMGMG_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgMG(i))
				end do
				
				filename =  'EMGLG_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgLG(i))
				end do
				
				filename =  'spksterminalSOL_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_SOL
					write(10,*) tspks_SOL_terminal(i), indiceMN_SOL_terminal(i)
				end do
				
				filename =  'spkssomaSOL_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_soma_SOL
					write(10,*) tspks_SOL_soma(i),indiceMN_SOL_soma(i)
				end do
				
				filename =  'spksterminalMG_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_MG
					write(10,*) tspks_MG_terminal(i), indiceMN_MG_terminal(i)
				end do
				
				filename =  'spksterminalLG_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_LG
					write(10,*) tspks_LG_terminal(i), indiceMN_LG_terminal(i)
				end do
				
				filename = 'IAFRSOL_FT_10MVC_trial'// trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRSOL(i))
				end do
				
				filename = 'IAFRMG_FT_10MVC_trial'// trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRMG(i))
				end do
				
				filename = 'IAFRLG_FT_10MVC_trial'// trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRLG(i))
				end do
				
				filename = 'SPKSIA_FT_10MVC_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_afferent_SOL
					write(10,*) tspks_SOL_afferentIA(i), indiceIA_SOL(i)
				end do
				
				filename = 'lengthSOL_FT_10MVC_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(length(i))
				end do
				
				filename = 'VELSOL_FT_10MVC_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(velocitySOL(i))
				end do
		
		
				allocate(condvel(size(motorUnitPools(1)%poolTerminalSpikes(:,2))))
				do i = 1, size(motorUnitPools(1)%poolTerminalSpikes(:,2))
					indice = int(motorUnitPools(1)%poolTerminalSpikes(i,2))
					paramTag = 'axonDelayCondVel'
					pool = 'SOL'
					velocity_m_s = conf%getparameterGauss(paramTag, pool, indice)
					condvel(i) = velocity_m_s
				end do 
			
			
				filename =  'AxonCondVel_SOL_FT_10MVC_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_SOL
					write(10,*) indiceMN_SOL_terminal(i),condvel(i)
				end do
				
				deallocate(tspks_SOL_soma)
				deallocate(tspks_SOL_terminal)
				deallocate(tspks_MG_terminal)
				deallocate(tspks_LG_terminal)
				deallocate(indiceMN_SOL_soma)
				deallocate(indiceMN_SOL_terminal)
				deallocate(indiceMN_MG_terminal)
				deallocate(indiceMN_LG_terminal)
				deallocate(condvel)
				deallocate(tspks_SOL_afferentIA)
				deallocate(indiceIA_SOL)
				deallocate(neuralTractPools)
				deallocate(motorUnitPools)  
				deallocate(interneuronPools)
				deallocate(afferentPools)	
				deallocate(t)
				deallocate(MNv_mV_SOL)
				deallocate(MNv_mV_MG)
				deallocate(MNv_mV_LG)
				deallocate(torque)
				deallocate(emgSOL)
				deallocate(forceSOL)
				deallocate(emgMG)
				deallocate(forceMG)
				deallocate(emgLG)
				deallocate(forceLG)
				deallocate(IaFRSOL)
				deallocate(IaFRMG)
				deallocate(IaFRLG)
				deallocate(velocitySOL)
				deallocate(length)
				
				trial = trial + 1
			end do
								
	end do 
end program ForceTask