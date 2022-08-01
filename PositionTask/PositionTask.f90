program PositionTask	
    use MotorUnitPoolClass
    use NeuralTractClass
    use InterneuronPoolClass
    use SynapticNoiseClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
	use ConstantSeedInitialize
	use SynapsesFactoryModule
	use MuscleHillClass
	use jointAnklePositionTaskClass
    use AfferentPoolClass
	implicit none
	
	type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength, L_soma_SOL, L_soma_MG, L_soma_LG, L_soma_TA, L_terminal_LG, L_terminal_MG, L_terminal_SOL, L_terminal_TA, L_afferent_SOL
    integer :: i, j, poolIndex
	real(wp), dimension(:), allocatable :: t, MNv_mV, angle, torque, muscletorque, ankletorque, passivetorque, gravitationaltorque, load,forceSOL, forceMG,forceLG, forceTA
	real(wp), dimension(:), allocatable :: emgSOL, emgMG, emgLG, emgTA
	real(wp), dimension(:), allocatable :: IaFRSOL, indiceIA_SOL, tspks_SOL_afferentIA, velocitySOL, length
	real(wp), dimension(:), allocatable :: IaFRMG, indiceIA_MG, tspks_MG_afferentIA
	real(wp), dimension(:), allocatable :: IaFRLG, indiceIA_LG, tspks_LG_afferentIA
	real(wp), dimension(:), allocatable :: IaFRTA, indiceIA_TA, tspks_TA_afferentIA
	real(wp), dimension(:), allocatable :: condvel
	real(wp), dimension(:), allocatable :: tspks_SOL_soma,tspks_MG_soma,tspks_LG_soma,tspks_TA_soma,indiceMN_SOL_SOMA, indiceMN_MG_SOMA, indiceMN_LG_SOMA, indiceMN_TA_SOMA
	real(wp), dimension(:), allocatable :: indiceMN_SOL_terminal, tspks_SOL_terminal, indiceMN_MG_terminal, tspks_MG_terminal, indiceMN_LG_terminal, tspks_LG_terminal,indiceMN_TA_terminal, tspks_TA_terminal
	real(wp) :: tic, toc
    type(gpf) :: gp 
    real(wp) :: FR, FRbasal, t1, velocity_m_s
    integer :: GammaOrder, trial, trials, k, indice
	integer :: cenMU, cenvel, MUconf, velconf
	character(len = 80) :: pool, muscle, group, fileNumber,paramTag
	character(len = 80) :: conffile = 'confPositionTask'
	character(len = 80) :: configurationFile, iterMU, itervel
	type(MotorUnitPool), dimension(:), allocatable, target :: motorUnitPools
    type(NeuralTract), dimension(:), allocatable :: neuralTractPools    
    type(InterneuronPool), dimension(:), allocatable, target :: interneuronPools    
    type(SynapticNoise), dimension(:), allocatable:: synapticNoisePools   
    type(AfferentPool), dimension(:), allocatable:: afferentPools 
	type(jointAnklePositionTask) :: ankle

		
	trial = 1
	trials = 10
	do while (trial <= trials )	
	do Sconf = 1,3
		
		scenario = Sconf - 1
		write(iterS, '(I3)')scenario
			
		call init_random_seed()
		configurationFile = trim(conffile) //'_S'// trim(adjustl(iterS))//'.rmto'
		print*, configurationFile
		conf = Configuration(configurationFile)
	
			allocate(afferentPools(4))
			pool = 'Ia'
			muscle = 'SOL'
			afferentPools(1) = AfferentPool(conf, pool, muscle)
		
			pool = 'Ia'
			muscle = 'MG'
			afferentPools(2) = AfferentPool(conf, pool, muscle)

			pool = 'Ia'
			muscle = 'LG'
			afferentPools(3) = AfferentPool(conf, pool, muscle)
			
			pool = 'Ia'
			muscle = 'TA'
			afferentPools(4) = AfferentPool(conf, pool, muscle)
		
				
			allocate(neuralTractPools(1))
			pool = 'CMExt'
			neuralTractPools(1) = NeuralTract(conf, pool)
			
			allocate(motorUnitPools(4))
			pool = 'SOL'
			motorUnitPools(1) = MotorUnitPool(conf, pool) 	
			
			pool = 'MG'
			motorUnitPools(2) = MotorUnitPool(conf, pool) 	
			
			pool = 'LG'
			motorUnitPools(3) = MotorUnitPool(conf, pool) 	
			
			pool = 'TA'
			motorUnitPools(4) = MotorUnitPool(conf, pool)
			
			ankle = jointAnklePositionTask(conf, motorUnitPools)
			
			allocate(interneuronPools(0))
			
			synapticNoisePools = synapseFactory(conf, neuralTractPools, motorUnitPools, interneuronPools, afferentPools)
			
			tf = conf%simDuration_ms
			dt = conf%timeStep_ms
			timeLength = int(tf/dt)
			
			
			allocate(t(timeLength))
			allocate(MNv_mV(timeLength))
			allocate(angle(timeLength))
			allocate(ankletorque(timeLength))
			allocate(muscletorque(timeLength))
			allocate(torque(timeLength))
			allocate(passivetorque(timeLength))
			allocate(gravitationaltorque(timeLength))
			allocate(load(timeLength))
			allocate(emgSOL(timeLength))
			allocate(forceSOL(timeLength))
			allocate(emgMG(timeLength))
			allocate(forceMG(timeLength))
			allocate(emgLG(timeLength))
			allocate(forceLG(timeLength))
			allocate(emgTA(timeLength))
			allocate(forceTA(timeLength))
			allocate(IaFRSOL(timeLength))
			allocate(IaFRMG(timeLength))
			allocate(IaFRLG(timeLength))
			allocate(IaFRTA(timeLength))
			allocate(velocitySOL(timeLength))
			allocate(length(timeLength))
			
			t = [(dt*(i-1), i=1, timeLength)]
			print *, timeLength, size(t)
			

			FRbasal = 40.0
			t1 = 2000
			

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
					if (i <= t1/dt) then
								FR = FRbasal
								GammaOrder = 10
							else if (i > t1/dt) then
								FR = 1000/13.8
								GammaOrder = 7
								
							end if
					do j = 1, size(neuralTractPools)
						call neuralTractPools(j)%atualizePool(t(i), FR, GammaOrder)
					end do
					
					do j = 1,3
						call motorUnitPools(j)%atualizeMotorUnitPool(t(i),50.0_wp, 48.0_wp)
						MNv_mV(i) = motorUnitPools(j)%v_mV(2)
					end do
					call motorUnitPools(4)%atualizeMotorUnitPool(t(i), 50.0_wp, 50.0_wp)
				
					do j = 1, 4
				
						call afferentPools(j)%atualizeAfferentPool(t(i), motorUnitPools(j)%spindle%IaFR_Hz)
					end do
				
					call ankle%atualizeAnkle(t(i))
					
					torque(i) = ankle%ankleTorque_Nm(i)
					muscletorque(i) = ankle%muscleTorque_Nm(i)
					torque(i) = ankle%Torque_Nm(i)
					IaFRSOL(i) = motorUnitPools(1)%spindle%IaFR_Hz
					IaFRMG(i) = motorUnitPools(2)%spindle%IaFR_Hz
					IaFRLG(i) = motorUnitPools(3)%spindle%IaFR_Hz
					IaFRTA(i) = motorUnitPools(4)%spindle%IaFR_Hz			
					
					
				end do
				call cpu_time(toc)
				
				print '(F15.6, A)', toc - tic, ' seconds'
				print*, 'GammaOrder = ', GammaOrder
				print*, 'FR = ', FR
				
				call neuralTractPools(1)%listSpikes()
				do j = 1, size(motorUnitPools)
					call motorUnitPools(j)%listSpikes()
					call motorUnitPools(j)%getMotorUnitPoolEMG()
				end do
				
				do j = 1, size(afferentPools)
					call afferentPools(j)%listSpikes()
				end do
		
				do j = 1, size(interneuronPools)
					call interneuronPools(j)%listSpikes()
				end do
				
				do i = 1,timeLength
					ankletorque(i) = ankle%ankleTorque_Nm(i)
					muscletorque(i) = ankle%muscleTorque_Nm(i)
					torque(i) = ankle%Torque_Nm(i)
					passivetorque(i) = ankle%passiveTorque_Nm(i)
					gravitationaltorque(i) = ankle%gravitationalTorque_Nm(i)
					load(i) = ankle%load(i)
					
					angle(i) = ankle%ankleAngle_rad(i)*180.0/pi
					forceSOL(i) = motorUnitPools(1)%HillMuscle%force(i)
					velocitySOL(i) = motorUnitPools(1)%HillMuscle%velocity_m_ms(i)
					forceMG(i) = motorUnitPools(2)%HillMuscle%force(i)
					forceLG(i) = motorUnitPools(3)%HillMuscle%force(i)
					forceTA(i) = motorUnitPools(4)%HillMuscle%force(i)
					emgSOL(i) = motorUnitPools(1)%emg(i)
					length(i) = motorUnitPools(1)%HillMuscle%length_m(i)
					emgMG(i) = motorUnitPools(2)%emg(i)
					emgLG(i) = motorUnitPools(3)%emg(i)
					emgTA(i) = motorUnitPools(4)%emg(i)
					
				end do
				
				L_soma_SOL = size(motorUnitPools(1)%poolSomaSpikes(:,1))				
				L_soma_MG = size(motorUnitPools(2)%poolSomaSpikes(:,1))
				L_soma_LG = size(motorUnitPools(3)%poolSomaSpikes(:,1))
				L_soma_TA = size(motorUnitPools(4)%poolSomaSpikes(:,1))
				L_terminal_SOL = size(motorUnitPools(1)%poolTerminalSpikes(:,1))
				L_terminal_MG = size(motorUnitPools(2)%poolTerminalSpikes(:,1))
				L_terminal_LG = size(motorUnitPools(3)%poolTerminalSpikes(:,1))
				L_terminal_TA = size(motorUnitPools(4)%poolTerminalSpikes(:,1))
				L_afferent_SOL = size(afferentPools(1)%poolTerminalSpikes(:,1))
				
				allocate(tspks_SOL_soma(L_soma_SOL))
				allocate(tspks_MG_soma(L_soma_MG))
				allocate(tspks_LG_soma(L_soma_LG))
				allocate(tspks_TA_soma(L_soma_TA))			
				allocate(tspks_SOL_terminal(L_terminal_SOL))
				allocate(tspks_MG_terminal(L_terminal_MG))
				allocate(tspks_LG_terminal(L_terminal_LG))
				allocate(tspks_TA_terminal(L_terminal_TA))			
				allocate(indiceMN_SOL_soma(L_soma_SOL))
				allocate(indiceMN_MG_soma(L_soma_MG))
				allocate(indiceMN_LG_soma(L_soma_LG))
				allocate(indiceMN_TA_soma(L_soma_TA))			
				allocate(indiceMN_SOL_terminal(L_terminal_SOL))
				allocate(indiceMN_MG_terminal(L_terminal_MG))
				allocate(indiceMN_LG_terminal(L_terminal_LG))
				allocate(indiceMN_TA_terminal(L_terminal_TA))
				allocate(tspks_SOL_afferentIA(L_afferent_SOL))
				allocate(indiceIA_SOL(L_afferent_SOL))
				
				do j = 1, L_soma_SOL
					tspks_SOL_soma(j) = motorUnitPools(1)%poolSomaSpikes(j,1)
					indiceMN_SOL_soma(j) = motorUnitPools(1)%poolSomaSpikes(j,2)
				end do
				
				do j = 1, L_soma_MG
					tspks_MG_soma(j) = motorUnitPools(2)%poolSomaSpikes(j,1)
					indiceMN_MG_soma(j) = motorUnitPools(2)%poolSomaSpikes(j,2)
				end do
								
				do j = 1, L_soma_LG
					tspks_LG_soma(j) = motorUnitPools(3)%poolSomaSpikes(j,1)
					indiceMN_LG_soma(j) = motorUnitPools(3)%poolSomaSpikes(j,2)
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
				
				do j = 1, L_terminal_TA
					tspks_TA_terminal(j) = motorUnitPools(4)%poolTerminalSpikes(j,1)
					indiceMN_TA_terminal(j) = motorUnitPools(4)%poolTerminalSpikes(j,2)
				end do
				
				do j = 1, L_afferent_SOL
					tspks_SOL_afferentIA(j) = afferentPools(1)%poolTerminalSpikes(j,1)
					indiceIA_SOL(j) = afferentPools(1)%poolTerminalSpikes(j,2)
				end do
		
				write(fileNumber, '(I2.1)')trial
				filename = 'AnkleTorquePT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(ankletorque(i))
				end do
					
				filename = 'MuscleTorquePT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(muscletorque(i))
				end do
								
				filename = 'TorquePT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(torque(i))
				end do
				
				filename = 'pasTorquePT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(passivetorque(i))
				end do
				
				filename = 'gravTorquePT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(gravitationaltorque(i))
				end do
				
				filename = 'load_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(load(i))
				end do
				
				filename = 'Angle_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(angle(i))
				end do
				
				filename = 'ForceSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceSOL(i))
				end do
				
				filename = 'VELSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(velocitySOL(i))
				end do
				
				filename = 'lengthSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(length(i))
				end do
				
				filename = 'EMGSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgSOL(i))
				end do
						
				filename =  'spkssomaSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_soma_SOL
					write(10,*) tspks_SOL_soma(i),indiceMN_SOL_soma(i)
				end do
				
				filename = 'SPKSTerminal_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_SOL
					write(10,*) tspks_SOL_terminal(i), indiceMN_SOL_terminal(i)
				end do
				
				filename = 'SPKSIA_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_afferent_SOL
					write(10,*) tspks_SOL_afferentIA(i), indiceIA_SOL(i)
				end do
				
				filename = 'IAFRSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRSOL(i))
				end do

				filename = 'ForceMGPT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceMG(i))
				end do
				
				filename = 'VELSOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(velocitySOL(i))
				end do
				
				filename = 'IAFRMG_PT_10MVC_AxonDelayOK_trial'// trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRMG(i))
				end do
							
				filename =  'spkssomaMG_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_soma_MG
					write(10,*) tspks_MG_soma(i),indiceMN_MG_soma(i)
				end do
							
				filename =  'spksterminalMG_PT_10MVC_AxonDelayOK_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_MG
					write(10,*) tspks_MG_terminal(i), indiceMN_MG_terminal(i)
				end do
				
				filename = 'EMGMGPT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgMG(i))
				end do

				filename = 'ForceLG_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceLG(i))
				end do
						
				filename = 'IAFRLG_PT_10MVC_AxonDelayOK_trial'// trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRLG(i))
				end do
				
				filename =  'spkssomaLG_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_soma_LG
					write(10,*) tspks_LG_soma(i),indiceMN_LG_soma(i)
				end do
				
				filename =  'spksterminalLG_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_LG
					write(10,*) tspks_LG_terminal(i), indiceMN_LG_terminal(i)
				end do
				
				filename = 'EMGLG_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgLG(i))
				end do
				
				filename = 'ForceTA_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(forceTA(i))
						
				filename = 'IAFRTA_PT_10MVC_randomseed_trial'// trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(IaFRTA(i))
				end do
	
				
				filename =  'spkssomaTA_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_soma_TA
					write(10,*) tspks_TA_soma(i),indiceMN_TA_soma(i)
				end do
				
				filename =  'spksterminalTA_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_TA
					write(10,*) tspks_TA_terminal(i), indiceMN_TA_terminal(i)
				end do
				
				filename = 'EMGTA_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(t)
					write(10,*) t(i),(emgTA(i))
				end do
							
				filename = 'SPKSCMext_PT_10MVC_trial'//trim(adjustl(fileNumber))//'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, size(neuralTractPools(1)%poolTerminalSpikes(:,1))
					write(10,*)neuralTractPools(1)%poolTerminalSpikes(i,1), neuralTractPools(1)%poolTerminalSpikes(i,2)
				end do
				
				allocate(condvel(size(motorUnitPools(1)%poolTerminalSpikes(:,2))))

				do i = 1, size(motorUnitPools(1)%poolTerminalSpikes(:,2))
					indice = int(motorUnitPools(1)%poolTerminalSpikes(i,2))
					paramTag = 'axonDelayCondVel'
					pool = 'SOL'
					velocity_m_s = conf%getparameterGauss(paramTag, pool, indice)
					condvel(i) = velocity_m_s
				end do 
					
				
				
				filename =  'AxonCondVel_SOL_PT_10MVC_randomseed_trial'//trim(adjustl(fileNumber)) //'_S'// trim(adjustl(iterS))//'.txt'
				open(UNIT = 10, FILE = filename, status = 'replace')
				do i = 1, L_terminal_SOL
					write(10,*) indiceMN_SOL_terminal(i),condvel(i)
				end do
									
				deallocate(tspks_SOL_soma)
				deallocate(tspks_MG_soma)
				deallocate(tspks_LG_soma)
				deallocate(tspks_TA_soma)
				deallocate(tspks_SOL_terminal)
				deallocate(tspks_MG_terminal)
				deallocate(tspks_LG_terminal)
				deallocate(tspks_TA_terminal)
				deallocate(indiceMN_SOL_soma)
				deallocate(indiceMN_MG_soma)
				deallocate(indiceMN_LG_soma)
				deallocate(indiceMN_TA_soma)
				deallocate(indiceMN_SOL_terminal)
				deallocate(indiceMN_MG_terminal)
				deallocate(indiceMN_LG_terminal)
				deallocate(indiceMN_TA_terminal)
				deallocate(condvel)
				deallocate(tspks_SOL_afferentIA)
				deallocate(indiceIA_SOL)
				deallocate(neuralTractPools)
				deallocate(motorUnitPools)  
				deallocate(interneuronPools)
				deallocate(afferentPools)	
				deallocate(t)
				deallocate(MNv_mV)
				deallocate(angle)
				deallocate(ankletorque)
				deallocate(muscletorque)
				deallocate(torque)
				deallocate(passivetorque)
				deallocate(gravitationaltorque)
				deallocate(load)
				deallocate(emgSOL)
				deallocate(forceSOL)
				deallocate(emgMG)
				deallocate(forceMG)
				deallocate(forceLG)
				deallocate(emgLG)
				deallocate(forceTA)	
				deallocate(emgTA)
				deallocate(IaFRSOL)
				deallocate(IaFRMG)
				deallocate(IaFRLG)
				deallocate(IaFRTA)
				deallocate(velocitySOL)
				deallocate(length)
						
				trial = trial + 1
				
			end do

	
		end do 
	end do 
	
	end program
	