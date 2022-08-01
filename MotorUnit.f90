! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2020  Renato Naville Watanabe
!						  Marina Cardoso de Oliveira

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

!     Contact: renato.watanabe@ufabc.edu.br
! '''

module MotorUnitClass
    ! '''
    ! Class that implements a motor unit model. Encompasses a motoneuron
    ! and a muscle unit.
    ! '''
    use CompartmentClass
    use ConfigurationClass
    use AxonDelayClass
    use DynamicalArrays
    use CharacterMatrixClass
    use SynapsePointerClass
    use randomGen
    implicit none
    private
    public :: MotorUnit

    type MotorUnit
        type(Configuration), pointer :: conf
        character(len = 2) ::  neuronKind
        character(len = 6) :: pool  
        character(len = 80) :: MUSpatialDistribution
        real(wp) :: tSomaSpike
        integer :: index
        type(Compartment), dimension(:), allocatable :: Compartments ! ## Dictionary of Compartment of the Motor Unit.
        real(wp) :: threshold_mV, position_mm, distance_mm, attenuationToSkin
        real(wp) :: timeWidening, ampEMG_mV, timeCteEMG_ms
        real(wp) :: muSectionPosition_mm(2)
        integer :: hrType, compNumber, somaIndex, lastCompIndex
        real(wp), dimension(:), allocatable :: v_mV, tSpikes, capacitanceInv
        real(wp), dimension(:), allocatable :: iIonic, iInjected, EqCurrent_nA
        real(wp), dimension(:,:), allocatable :: G 
        real(wp) :: MNRefPer_ms, stimulusPositiontoTerminal, nerveLength   
        character(len = 3) :: nerve
        type(AxonDelay) :: Delay
        character(len = 5) :: stimulusCompartment
        real(wp) :: stimulusMeanFrequency_Hz, stimulusPulseDuration_ms, stimulusIntensity_mA
        real(wp) :: stimulusStart_ms,stimulusStop_ms, stimulusModulationStart_ms, stimulusModulationStop_ms
        real(wp), dimension(:), allocatable :: nerveStimulus_mA
        real(wp), dimension(:), allocatable :: somaSpikeTrain, lastCompSpikeTrain, terminalSpikeTrain
        real(wp) :: TwitchTc_ms, TwitchAmp_N, bSat, twTet
		real(wp) :: TwitchTcRaikova_ms, bSatRaikova, TwitchHrRaikova_ms, twTetRaikova
        type(CharacterMatrix) :: SynapsesOut
        integer, dimension(:), allocatable :: indicesOfSynapsesOnTarget
        type(SynapsePointer), dimension(:), allocatable :: transmitSpikesThroughSynapses
        real(wp) :: axonStimModulation

        contains
            procedure :: atualizeMotorUnit
            procedure :: atualizeCompartments
            procedure :: addCompartmentSpike
            procedure :: atualizeDelay
            procedure :: reset
            procedure :: getEMG
            procedure :: transmitSpikes
            procedure :: createStimulus

    end type MotorUnit

    interface MotorUnit
        module procedure init_MotorUnit
    end interface MotorUnit
    
    contains

    type(MotorUnit) function init_MotorUnit(conf, pool, index, neuronKind, muscleThickness, skinThickness)
        ! '''
        ! Constructor

        ! - Inputs:
        !     + **conf**: Configuration object with the simulation parameters.

        !     + **pool**: string with Motor unit pool to which the motor
        !     unit belongs.

        !     + **index**: integer corresponding to the motor unit order in
        !     the pool, according to the Henneman's principle (size principle).

        !     + **neuronKind**: string with the type of the motor unit. It can
        !     be *S* (slow), *FR* (fast and resistant), and
        !     *FF* (fast and fatigable).
        !     
        !     + **muscleThickness** :float, muscle Thickness, in mm
        !     
        !     + **skinThickness** : float, skin Thickness, in mm   
        ! '''
        class(Configuration), intent(in), target :: conf    
        character(*), intent(in) :: pool
        integer, intent(in) :: index
        character(*), intent(in) ::  neuronKind
        real(wp), intent(in) :: muscleThickness, skinThickness
        character(len=80) :: paramTag, paramChar
        real(wp) :: paramReal
        integer :: NumberOfAxonNodes
        real(wp) :: radius, angle, x, y
        real(wp) :: randomNumber
        character(len = 9):: compKind
        real(wp), dimension(:), allocatable :: gCoupling_muS
        integer :: i, j
        real(wp) :: rAxis1, rAxis2, cytR
        real(wp), dimension(:), allocatable :: gLeak, capacitance_nF, EqPot, IPump, compLength
        real(wp), dimension(:,:), allocatable :: GC, GL
        real(wp) :: dynamicNerveLength, delayLength
        real(wp) :: stimulusInDynamiAxon
        integer :: simDurationSteps, numberOfSteps, startStep, stimPulseDurationSteps
        real(wp) :: stimulusFrequency_Hz, stimulusPeriod_ms
        character(len = 80) :: activationModel

        init_MotorUnit%index = index         
        
        ! ## Configuration object with the simulation parameters.
        init_MotorUnit%conf => conf

        ! ## String with the type of the motor unit. It can be
        ! ## *S* (slow), *FR* (fast and resistant) and
        ! ## *FF** (fast and fatigable).
        init_MotorUnit%neuronKind = neuronKind

        init_MotorUnit%pool = pool
        
        ! # Neural compartments
        ! ## The instant of the last spike of the Motor unit
        ! ## at the Soma compartment.
        init_MotorUnit%tSomaSpike = -1e10

        
        paramTag = 'NumberAxonNodes'
        paramChar =  init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)paramReal
        NumberOfAxonNodes = int(paramReal)
        
        ! ## Integer corresponding to the motor unit order in the pool, according to the Henneman's principle (size principle).
                   
        ! ## Anatomical position of the neuron in the spinal cord, in mm.
        paramTag = 'position'
        paramChar =  init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%position_mm 
          
        ! ## Value of the membrane potential, in mV, that is considered a spike.
        paramTag = 'threshold'
        paramChar =  init_MotorUnit%conf%parameterSet(paramTag, init_MotorUnit%pool, init_MotorUnit%index)
        read(paramChar, *)init_MotorUnit%threshold_mV 
        
        ! # EMG data
        paramTag = 'MUSpatialDistribution'
        init_MotorUnit%MUSpatialDistribution = init_MotorUnit%conf%parameterSet(paramTag,pool, index)       
        if (trim(init_MotorUnit%MUSpatialDistribution).eq.'random') then
            call random_number(randomNumber)
            radius = sqrt(((muscleThickness/2.0)**2) * randomNumber)
            call random_number(randomNumber)
            angle = 2.0 * pi * randomNumber
        end if
        
        x = radius * sin(angle)
        y = radius * cos(angle)
        ! ## Anatomical coordinate of the muscle unit in a muscle section, in (mm,mm).
        init_MotorUnit%muSectionPosition_mm(1) = x
        init_MotorUnit%muSectionPosition_mm(2) = y

        ! ## Distance of the MU to the EMG elctrode, in mm.
        init_MotorUnit%distance_mm = sqrt((x + muscleThickness/2.0 + skinThickness)**2 + y**2)

        ! ## Attenuation of the MUAP amplitude, as measured in the electrode.
        init_MotorUnit%attenuationToSkin = exp(-init_MotorUnit%distance_mm / init_MotorUnit%conf%EMGAttenuation_mm1)

        ! ## Widening of the MUAP duration, as measured in the electrode.
        init_MotorUnit%timeWidening = 1 + init_MotorUnit%conf%EMGWidening_mm1 * init_MotorUnit%distance_mm
        
        ! ## Type of the Hermitez-Rodiguez curve. It can be 1 or 2.
        call random_number(randomNumber)
        init_MotorUnit%hrType = int(randomNumber+1.5)
        
        ! ## MUAP amplitude in mV.
        paramTag = 'EMGAmplitude'
        paramChar =  init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%ampEMG_mV
		!print*,'Parametro amplitude = ',init_MotorUnit%ampEMG_mV
        init_MotorUnit%ampEMG_mV = (1+0.05*randn())*init_MotorUnit%ampEMG_mV
        init_MotorUnit%ampEMG_mV = init_MotorUnit%ampEMG_mV * init_MotorUnit%attenuationToSkin

        ! ## MUAP time constant, in ms.
        paramTag = 'EMGDuration'
        paramChar =  init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%timeCteEMG_ms
		
        init_MotorUnit%timeCteEMG_ms = (1+0.05*randn())*init_MotorUnit%timeCteEMG_ms
		init_MotorUnit%timeCteEMG_ms = init_MotorUnit%timeCteEMG_ms * init_MotorUnit%timeWidening
       
        ! ## Number of compartments.
        init_MotorUnit%compNumber = 2+NumberOfAxonNodes
        
        
        if (allocated(init_MotorUnit%Compartments)) then
            deallocate(init_MotorUnit%Compartments)    
        end if
        allocate(init_MotorUnit%Compartments(init_MotorUnit%compNumber))
        
        compKind = 'dendrite'
        init_MotorUnit%Compartments(1) = Compartment(compKind, init_MotorUnit%conf, &
                                                     init_MotorUnit%pool, init_MotorUnit%index, &
                                                     init_MotorUnit%neuronKind)

        compKind = 'soma'
        init_MotorUnit%Compartments(2) = Compartment(compKind, init_MotorUnit%conf, &
                                                    init_MotorUnit%pool, init_MotorUnit%index, &
                                                    init_MotorUnit%neuronKind)                                                     
        ! ## index of the soma compartment.
        init_MotorUnit%somaIndex = 2
        do i = 3, NumberOfAxonNodes+2, 2

            compKind = 'internode'
            init_MotorUnit%Compartments(i) = Compartment(compKind, init_MotorUnit%conf, &
                                                         init_MotorUnit%pool, init_MotorUnit%index, &
                                                         init_MotorUnit%neuronKind)
			print *,'i=',i,'[MotorUnit.f90]#tst#', compkind	
			

			compKind = 'node'
            init_MotorUnit%Compartments(i+1) = Compartment(compKind, init_MotorUnit%conf, &
                                                           init_MotorUnit%pool, init_MotorUnit%index, &
                                                           init_MotorUnit%neuronKind)
			print *,'i=',i+1,'[MotorUnit.f90]#tst#', compkind												   
			
        end do

        
        ! ## Vector with membrane potential,in mV, of all compartments. 
        allocate(init_MotorUnit%v_mV(init_MotorUnit%compNumber))
        do i = 1, init_MotorUnit%compNumber
            init_MotorUnit%v_mV(i) = 0.0    
        end do
        ! ## Vector with the last instant of spike of all compartments. 
        allocate(init_MotorUnit%tSpikes(init_MotorUnit%compNumber))
        do i = 1, init_MotorUnit%compNumber
            init_MotorUnit%tSpikes(i) = -1e10   
        end do
        
        allocate(gCoupling_muS(init_MotorUnit%compNumber))        
        
        paramTag = 'cytR'
        paramChar =  init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)cytR

        do i = 1, (init_MotorUnit%compNumber-1)
            ! '''
            ! Calculates the coupling conductance between two compartments.

            ! - Inputs: 
            !     + **cytR**: Cytoplasmatic resistivity in \f$\Omega\f$.cm.

            !     + **lComp1, lComp2**: length of the compartments in \f$\mu\f$m.

            !     + **dComp1, dComp2**: diameter of the compartments in \f$\mu\f$m.

            ! - Output:
            !     + coupling conductance in \f$\mu\f$S.

            ! The coupling conductance between compartment 1 and 2 is
            ! computed by the following equation:

            ! \f{equation}{
            !     g_c = \frac{2.10^2}{\frac{R_{cyt}l_1}{\pi r_1^2}+\frac{R_{cyt}l_2}{\pi r_2^2}}
            ! \f}
            ! where \f$g_c\f$ is the coupling conductance [\f$\mu\f$S], \f$R_{cyt}\f$ is the
            ! cytoplasmatic resistivity [\f$\Omega\f$.cm], \f$l_1\f$ and \f$l_2\f$
            ! are the lengths [\f$\mu\f$m] of compartments 1 and 2, respectively and
            ! \f$r_1\f$ and \f$r_2\f$ are the radius [\f$\mu\f$m] of compartments 1 and
            ! 2, respectively.
            ! '''
            rAxis2=(cytR * init_MotorUnit%Compartments(i+1)%length_mum)/(pi*(init_MotorUnit%Compartments(i+1)%diameter_mum/2.0)**2)
            rAxis1=(cytR * init_MotorUnit%Compartments(i)%length_mum)/(pi*(init_MotorUnit%Compartments(i)%diameter_mum/2.0)**2)
                   
            gCoupling_muS(i) = 200.0 / (rAxis1 + rAxis2)
        end do
        gCoupling_muS(init_MotorUnit%compNumber) = 0.0
        
        
        allocate(gLeak(init_MotorUnit%compNumber))
        allocate(capacitance_nF(init_MotorUnit%compNumber))
        allocate(EqPot(init_MotorUnit%compNumber))
        allocate(IPump(init_MotorUnit%compNumber))
        allocate(compLength(init_MotorUnit%compNumber))
        
        do i = 1, init_MotorUnit%compNumber
            capacitance_nF(i) = init_MotorUnit%Compartments(i)%capacitance_nF
            gLeak(i) = init_MotorUnit%Compartments(i)%gLeak_muS
            EqPot(i) = init_MotorUnit%Compartments(i)%EqPot_mV
            IPump(i) = init_MotorUnit%Compartments(i)%IPump_nA
            compLength(i) = init_MotorUnit%Compartments(i)%length_mum
            init_MotorUnit%v_mV(i) = init_MotorUnit%Compartments(i)%EqPot_mV
        end do
        
        ! ## Vector with  the inverse of the capacitance of all compartments. 
        init_MotorUnit%capacitanceInv = 1.0 / capacitance_nF

        
        ! ## Vector with current, in nA,  of each compartment coming from other elements of the model. For example 
        ! ## from ionic channels and synapses.       
        allocate(init_MotorUnit%iIonic(init_MotorUnit%compNumber))
        ! ## Vector with the current, in nA, injected in each compartment.
        allocate(init_MotorUnit%iInjected(init_MotorUnit%compNumber))
        
        init_MotorUnit%iIonic(:) = 0.0
        init_MotorUnit%iInjected(:) = 0.0
        
        

        allocate(GC(init_MotorUnit%compNumber,init_MotorUnit%compNumber))
        
        GC(:,:) = 0.0
    
        do i = 1, init_MotorUnit%compNumber
            ! '''
            ! Computes the Coupling Matrix to be used in the dVdt function of the N compartments of the motor unit. 
            ! The Matrix uses the values obtained with the function calcGcoupling.
        
            ! - Inputs: 
            !     + **gc**: the vector with N elements, with the coupling conductance of each compartment of the Motor Unit.

            ! - Output:
            !     + the GC matrix


            ! \f{equation}{
            ! GC = \left[\begin{array}{cccccccc}
            ! -g_c[0]&g_c[0]&0&...&...&0&0&0\\
            ! g_c[0]&-g_c[0]-g_c[1]&g_c[1]&0&...&...&0&0\\
            ! \vdots&&\ddots&&...&&0&0 \\
            ! 0&...&g_c[i-1]&-g_c[i-1]-g_c[i]&g_c[i]&0&...&0\\
            ! 0&0&0&...&...&&&0\\
            ! 0&&...&&g_c[N-2]&-g_c[N-2]-g_c[N-1]&g_c[N-1]&0\\
            ! 0&...&0&&&0&g_c[N-1]&-g_c[N-1]\end{array}\right] 
            ! \f} 
            ! '''
            if (i.eq.1) then
                GC(i,i:i+1) = [-gCoupling_muS(i), gCoupling_muS(i)]
            else if (i.eq.init_MotorUnit%compNumber) then
               GC(i,i-1:i) = [gCoupling_muS(i-1), -gCoupling_muS(i-1)]
            else
               GC(i,i-1:i+1) = [gCoupling_muS(i-1),-gCoupling_muS(i-1)-gCoupling_muS(i), gCoupling_muS(i)]
            end if
        end do        
        
        allocate(GL(init_MotorUnit%compNumber, init_MotorUnit%compNumber))
        GL(:,:) = 0.0
        do i = 1, init_MotorUnit%compNumber
            GL(i,i) = -gLeak(i)
        end do
        
        ! ## Matrix of the conductance of the motoneuron. Multiplied by the vector self.v_mV,
        ! ## results in the passive currents of each compartment.
        init_MotorUnit%G = GC + GL

        allocate(init_MotorUnit%EqCurrent_nA(init_MotorUnit%compNumber))
        init_MotorUnit%EqCurrent_nA = matmul(-GL, EqPot) + IPump 

        
        

        ! ## index of the last compartment.
        init_MotorUnit%lastCompIndex = 1*init_MotorUnit%compNumber
        
        ! ## Refractory period, in ms, of the motoneuron.
        paramTag = 'MNSomaRefPer'
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%MNRefPer_ms
        
        ! # delay
        ! ## String with type of the nerve. It can be PTN (posterior tibial nerve) or CPN
        ! ## (common peroneal nerve).
        if (trim(init_MotorUnit%pool).eq.'SOL'.or.trim(init_MotorUnit%pool).eq.'MG'.or.trim(init_MotorUnit%pool).eq.'LG') then
            init_MotorUnit%nerve = 'PTN'
        else if (trim(init_MotorUnit%pool).eq.'TA') then
            init_MotorUnit%nerve = 'CPN'
        end if

        ! ## Distance, in m, of the stimulus position to the terminal.
        paramTag =  'stimDistToTerm_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusPositiontoTerminal
        ! ## AxonDelay object of the motor unit.
        if (NumberOfAxonNodes.eq.0) then
            dynamicNerveLength = 0
        else
            dynamicNerveLength = sum(compLength(3:init_MotorUnit%compNumber)) * 1E-6
        end if

        paramTag = 'nerveLength_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%nerveLength

        delayLength =  init_MotorUnit%nerveLength - dynamicNerveLength

        if (init_MotorUnit%stimulusPositiontoTerminal < delayLength) then
            init_MotorUnit%Delay = AxonDelay(init_MotorUnit%conf, init_MotorUnit%nerve, &
                                             init_MotorUnit%pool, delayLength, & 
                                             init_MotorUnit%stimulusPositiontoTerminal, &
                                             init_MotorUnit%index)
            init_MotorUnit%stimulusCompartment = 'delay'
			
        else
            stimulusInDynamiAxon = -1.0
            init_MotorUnit%Delay = AxonDelay(init_MotorUnit%conf, init_MotorUnit%nerve, &
                                             init_MotorUnit%pool, delayLength, &
                                             stimulusInDynamiAxon, init_MotorUnit%index)
            init_MotorUnit%stimulusCompartment = 'dynam'    
        end if
        ! # Nerve stimulus function
		
        
		paramTag = 'stimFrequency_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusMeanFrequency_Hz
        
        paramTag = 'stimPulseDuration_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusPulseDuration_ms
        
        paramTag = 'stimIntensity_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusIntensity_mA
        
        paramTag = 'stimStart_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusStart_ms
        
        paramTag = 'stimStop_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusStop_ms
        
        paramTag = 'stimModulationStart_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusModulationStart_ms
        
        paramTag = 'stimModulationStop_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%stimulusModulationStop_ms
        if (init_MotorUnit%stimulusStop_ms >= init_MotorUnit%conf%simDuration_ms) then
            init_MotorUnit%stimulusStop_ms = init_MotorUnit%conf%simDuration_ms - 1
        end if
        
        paramTag = 'stimModulation_' // trim(init_MotorUnit%nerve)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_MotorUnit%axonStimModulation
        
        
        
        startStep = nint(init_MotorUnit%stimulusStart_ms / init_MotorUnit%conf%timeStep_ms)
        
        ! ## Vector with the nerve stimulus, in mA.
        simDurationSteps = nint(init_MotorUnit%conf%simDuration_ms/init_MotorUnit%conf%timeStep_ms)

        allocate(init_MotorUnit%nerveStimulus_mA(simDurationSteps))
        init_MotorUnit%nerveStimulus_mA(:) = 0.0
        
        stimPulseDurationSteps = nint(init_MotorUnit%stimulusPulseDuration_ms/init_MotorUnit%conf%timeStep_ms)
        do i = 1, simDurationSteps
            if ((i * init_MotorUnit%conf%timeStep_ms >= init_MotorUnit%stimulusStart_ms).and.&
                (i * init_MotorUnit%conf%timeStep_ms <= init_MotorUnit%stimulusStop_ms)) then
                    if ((i * init_MotorUnit%conf%timeStep_ms > init_MotorUnit%stimulusModulationStart_ms).and.&
                        (i * init_MotorUnit%conf%timeStep_ms < init_MotorUnit%stimulusModulationStop_ms)) then 
                        stimulusFrequency_Hz = init_MotorUnit%stimulusMeanFrequency_Hz + init_MotorUnit%axonStimModulation
                    else
                        stimulusFrequency_Hz = init_MotorUnit%stimulusMeanFrequency_Hz
                    end if
                    if (stimulusFrequency_Hz > 0.0) then
                        stimulusPeriod_ms = 1000.0 / stimulusFrequency_Hz
                        numberOfSteps = nint(stimulusPeriod_ms / init_MotorUnit%conf%timeStep_ms)
                        if (mod(i - startStep,numberOfSteps) == 0) then 
                            init_MotorUnit%nerveStimulus_mA(i:i+stimPulseDurationSteps) = init_MotorUnit%stimulusIntensity_mA
                        end if
                    end if
            end if
        end do
        
        
        paramTag = 'activationModel' ! 	#SOCDS
        activationModel = init_MotorUnit%conf%parameterSet(paramTag, pool,index)
        ! ## Contraction time of the twitch muscle unit, in ms.
		if (trim(activationModel).eq.'SOCDS') then
			paramTag = 'twitchTimePeak'
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
			read(paramChar, *)init_MotorUnit%TwitchTc_ms !L! twitchTimePeak
		
		! ## Amplitude of the muscle unit twitch, in N.
			paramTag = 'twitchPeak'
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
			read(paramChar, *)init_MotorUnit%TwitchAmp_N
		

		! ## Parameter of the saturation.
			paramTag = 'bSat'// trim(activationModel)
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
			read(paramChar, *)init_MotorUnit%bSat
        
		! ## Twitch- tetanus relationship
        paramTag = 'twTet' // trim(activationModel)
        paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)       
			read(paramChar,*)init_MotorUnit%twTet
        
		else
			paramTag = 'twitchTPeakRaikova' !alt!
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index) !alt!
			read(paramChar, *)init_MotorUnit%TwitchTcRaikova_ms !alt! twitchTimePeak
			
			paramTag = 'twitchTHalfRaikova' !alt!
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index) !alt!
			read(paramChar, *)init_MotorUnit%TwitchHrRaikova_ms !alt! half relaxation time
			
			paramTag = 'bSatRaikova'!// trim(activationModel)
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
			read(paramChar, *)init_MotorUnit%bSatRaikova
			
			! ## Twitch- tetanus relationship
			paramTag = 'twTet' // trim(activationModel)
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)       
			read(paramChar,*)init_MotorUnit%twTetRaikova
			
			paramTag = 'twitchPeak'
			paramChar = init_MotorUnit%conf%parameterSet(paramTag, pool, index)
			read(paramChar, *)init_MotorUnit%TwitchAmp_N
			! !!! 
		end if
		
        ! ## EMG data
        
        !TODO:
        ! ## Build synapses       
        init_MotorUnit%SynapsesOut = CharacterMatrix()
            

    end function

    subroutine atualizeMotorUnit(self, t, v_mV)
        ! '''
        ! Atualize the dynamical and nondynamical (delay) parts of the motor unit.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(MotorUnit), intent(inout) :: self 
        real(wp), intent(in) :: t
        real(wp), dimension(self%compNumber) :: v_mV

        call self%atualizeCompartments(t, v_mV)
        call self%atualizeDelay(t)
    end subroutine

    subroutine atualizeCompartments(self, t, v_mV)
        ! '''
        ! Atualize all neural compartments.

        ! - Inputs:
        !     + **t**: current instant, in ms.

        ! '''        
        class(MotorUnit), intent(inout) :: self 
        real(wp), intent(in) :: t
        real(wp), dimension(self%compNumber) :: v_mV
        integer :: i
        
        self%v_mV(:) = v_mV
		 
        do i = self%somaIndex, self%compNumber
		if ((self%v_mV(i) > self%threshold_mV).and.(t-self%tSpikes(i) > self%MNRefPer_ms)) then 
                call self%addCompartmentSpike(t, i)
		end if
        end do  
        
    end subroutine

    subroutine addCompartmentSpike(self, t, comp)
        ! '''
        ! When the soma potential is above the threshold a spike is added to the soma.

        ! - Inputs:
        !     + **t**: current instant, in ms.

        !     + **comp**: integer with the compartment index.
        ! '''
        class(MotorUnit), intent(inout) :: self 
        real(wp), intent(in) :: t
        integer, intent(in) :: comp
        integer :: i, j

        self%tSpikes(comp) = t
		
        if (comp.eq.self%somaIndex) then
		
            call AddToList(self%somaSpikeTrain, t)
			
            call self%transmitSpikes(t)
        end if

        if (comp.eq.self%lastCompIndex) then
            call AddToList(self%lastCompSpikeTrain, t)
            call self%Delay%addSpinalSpike(t)
        end if
        
        do i = 1, self%Compartments(comp)%numberChannels
            do j = 1, self%Compartments(comp)%Channels(i)%lenStates
                call self%Compartments(comp)%Channels(i)%condState(j)%changeState(t)
            end do
        end do
    end subroutine
              
    subroutine atualizeDelay(self, t)
        ! '''
        ! Atualize the terminal spike train, by considering the Delay of the nerve.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(MotorUnit), intent(inout) :: self 
        real(wp), intent(in) :: t
        integer :: i, j
		
        do while (self%Delay%terminalSpikeTrain_fila%n > 0 .and. abs(t - self%Delay%terminalSpikeTrain_fila%top()) < 1e-3) 
			call self%Delay%terminalSpikeTrain_fila%pop()
			
		   call AddToList(self%terminalSpikeTrain, t)
			
        end do
        
		do while (self%Delay%terminalSpikeTrain_fila%n > 0 .and. abs(self%Delay%terminalSpikeTrain_fila%top()) < t) 
			call self%Delay%terminalSpikeTrain_fila%pop()
		end do
        
        ! # Check whether there is antidromic impulse reaching soma or RC
        if (allocated(self%Delay%antidromicSpikeTrain)) then
            if (self%Delay%indexAntidromicSpike <= size(self%Delay%antidromicSpikeTrain)) then
                if (abs(t - self%Delay%antidromicSpikeTrain(self%Delay%indexAntidromicSpike)) < 1e-2) then 

                    ! # Considers only MN-RC connections
                    call self%transmitSpikes(t)
                    
                    ! # Refractory period of MN soma

                    if (t - self%tSpikes(self%somaIndex) > self%MNRefPer_ms) then

					   self%tSpikes(self%somaIndex) = t
                        call AddToList(self%somaSpikeTrain, t)
                        self%Delay%indexAntidromicSpike = self%Delay%indexAntidromicSpike  + 1
                        do i = 1, self%Compartments(self%somaIndex)%numberChannels
                            do j = 1, self%Compartments(self%somaIndex)%Channels(i)%lenStates
                                call self%Compartments(self%somaIndex)%Channels(i)%condState(j)%changeState(t)
                            end do
                        end do
                    end if
                end if
            end if        
        end if
        i = nint(t/self%conf%timeStep_ms) + 1
        
        if (trim(self%stimulusCompartment).eq.'delay') then
            call self%Delay%atualizeStimulus(t, self%nerveStimulus_mA(i))
        end if

    end subroutine

    
   subroutine transmitSpikes(self, t)
            !     '''
            !     - Inputs:
            !         + **t**: current instant, in ms.
            !     '''
            class(MotorUnit), intent(inout) :: self
            real(wp), intent(in) :: t     
            integer :: i
            
            if (allocated(self%indicesOfSynapsesOnTarget)) then
                do i = 1, size(self%indicesOfSynapsesOnTarget)
                    call self%transmitSpikesThroughSynapses(i)%synapse%receiveSpike(t, self%indicesOfSynapsesOnTarget(i))
                end do
            end if
            
    end subroutine

    real(wp) function getEMG(self, t) result(emg)
        ! '''
        ! Computes the EMG signal of the motor unit
        ! '''
        class(MotorUnit), intent(inout) :: self     
        real(wp), intent(in) :: t
        real(wp), dimension(:), allocatable :: numberOfSpikesUntilt
        real(wp) :: ta, newSpike
        integer :: i
        integer :: numberOfSpikes, sizeOfNumberOfSpikesUntilt

        sizeOfNumberOfSpikesUntilt = 0
        emg = 0.0
        newSpike = 0.0
        ta = 0.0        
        if (allocated(numberOfSpikesUntilt)) deallocate(numberOfSpikesUntilt)
        
        if (allocated(self%terminalSpikeTrain)) then 
            numberOfSpikes = size(self%terminalSpikeTrain)
            if (numberOfSpikes == 0) then
                emg = 0.0
            else
                do i = 1, numberOfSpikes
                    newSpike = self%terminalSpikeTrain(i)
                    if (newSpike < t .and. newSpike >= t-9*self%timeCteEMG_ms) then                                             
                        call AddToList(numberOfSpikesUntilt, newSpike)
                    end if
                end do
            end if
        end if
	
        
        if (allocated(numberOfSpikesUntilt)) then
            sizeOfNumberOfSpikesUntilt = size(numberOfSpikesUntilt)
            do i = 1, sizeOfNumberOfSpikesUntilt
                ta = t - numberOfSpikesUntilt(i) - 3.0 * self%timeCteEMG_ms
                if (ta <= 6 * self%timeCteEMG_ms) then                
                    if (self%hrType == 1) then
                        emg = emg + 1.19 * self%ampEMG_mV * ta * &
                                    exp(-(ta/self%timeCteEMG_ms)**2) / self%timeCteEMG_ms
						
                    else
                        emg = emg + 0.69 * self%ampEMG_mV * &
                                    (1 - 2*((ta / self%timeCteEMG_ms)**2)) *&
                                     exp(-(ta/self%timeCteEMG_ms)**2)
						
                    end if
                end if
            end do
        end if

        if (allocated(numberOfSpikesUntilt)) deallocate(numberOfSpikesUntilt)        
    end function

    subroutine reset(self)
        class(MotorUnit), intent(inout) :: self     
        integer :: i
        
        self%tSomaSpike = -1e10
        do i = 1, self%compNumber
            self%v_mV(i) = self%Compartments(i)%EqPot_mV
            call self%Compartments(i)%reset()
        end do
        
        call self%Delay%reset()
        
        self%tSpikes(:) = -1e10
        self%iIonic(:) = 0.0
        self%iInjected(:) = 0.0

        if (allocated(self%somaSpikeTrain)) deallocate(self%somaSpikeTrain)
        ! ## Vector with the instants of spikes at the last compartment.
        if (allocated(self%lastCompSpikeTrain)) deallocate(self%lastCompSpikeTrain)
        ! ## Vector with the instants of spikes at the terminal.
        if (allocated(self%terminalSpikeTrain)) deallocate(self%terminalSpikeTrain)

        !print *, self%pool, self%index
    end subroutine

    subroutine createStimulus(self)
        ! '''

        ! '''
        class(MotorUnit), intent(inout) :: self
        character(len = 80) :: paramtag, paramChar
        integer :: startStep, simDurationSteps, stimPulseDurationSteps, numberOfSteps, i
        real(wp) :: stimulusFrequency_Hz, stimulusPeriod_ms
        

        paramTag = 'stimFrequency_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusMeanFrequency_Hz
        
        paramTag = 'stimPulseDuration_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusPulseDuration_ms
        
        paramTag = 'stimIntensity_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusIntensity_mA
        
        paramTag = 'stimStart_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusStart_ms
        
        paramTag = 'stimStop_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusStop_ms
        
        paramTag = 'stimModulationStart_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusModulationStart_ms
        
        paramTag = 'stimModulationStop_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%stimulusModulationStop_ms
        if (self%stimulusStop_ms >= self%conf%simDuration_ms) then
            self%stimulusStop_ms = self%conf%simDuration_ms - 1
        end if
        
        paramTag = 'stimModulation_' // trim(self%nerve)
        paramChar = self%conf%parameterSet(paramTag, self%pool, self%index)
        read(paramChar, *)self%axonStimModulation        
        
        startStep = nint(self%stimulusStart_ms / self%conf%timeStep_ms)
        
        ! ## Vector with the nerve stimulus, in mA.
        simDurationSteps = nint(self%conf%simDuration_ms/self%conf%timeStep_ms)

        self%nerveStimulus_mA(:) = 0.0
        
        stimPulseDurationSteps = nint(self%stimulusPulseDuration_ms/self%conf%timeStep_ms)
        do i = 1, simDurationSteps
            if ((i * self%conf%timeStep_ms >= self%stimulusStart_ms).and.&
                (i * self%conf%timeStep_ms <= self%stimulusStop_ms)) then
                    if ((i * self%conf%timeStep_ms > self%stimulusModulationStart_ms).and.&
                        (i * self%conf%timeStep_ms < self%stimulusModulationStop_ms)) then 
                        stimulusFrequency_Hz = self%stimulusMeanFrequency_Hz + self%axonStimModulation
                    else
                        stimulusFrequency_Hz = self%stimulusMeanFrequency_Hz
                    end if
                    if (stimulusFrequency_Hz > 0.0) then
                        stimulusPeriod_ms = 1000.0 / stimulusFrequency_Hz
                        numberOfSteps = nint(stimulusPeriod_ms / self%conf%timeStep_ms)
                        if (mod(i - startStep,numberOfSteps) == 0) then 
                            self%nerveStimulus_mA(i:i+stimPulseDurationSteps) = self%stimulusIntensity_mA
                        end if
                    end if
            end if
        end do
    end subroutine
    
end module MotorUnitClass
    

    
    


