! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2019  Renato Naville Watanabe

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

module AfferentUnitClass
    ! '''
    ! Class that implements a motor unit model. Encompasses a motoneuron
    ! and a muscle unit.
    ! '''
    use CompartmentClass
    use AxonDelayClass
    use PointProcessGeneratorClass
    use ConfigurationClass
    use CharacterMatrixClass
    use SynapsePointerClass
    use DynamicalArrays
    use randomGen
    implicit none
    private
    public :: AfferentUnit

    type AfferentUnit
        type(Configuration), pointer :: conf
        character(len = 6) :: pool, muscle, neuronKind
        integer :: index, compNumber, lastCompIndex
        real(wp) :: timeStep_ms, timeStepByTwo_ms, timeStepBySix_ms
        type(Compartment), dimension(:), allocatable :: Compartments
        real(wp) :: threshold_mV, AFRefPer_ms
        real(wp), dimension(:), allocatable :: v_mV, tSpikes, capacitanceInv
        real(wp), dimension(:), allocatable :: iIonic, iInjected, EqCurrent_nA
        real(wp), dimension(:,:), allocatable :: G
        character(len = 3) :: nerve
        real(wp) :: stimulusPositiontoTerminal, frequencyThreshold_Hz
        type(AxonDelay) :: Delay
        integer :: stimulusCompartment
        real(wp), dimension(:), allocatable :: nerveStimulus_mA, terminalSpikeTrain, lastCompSpikeTrain
        integer :: GammaOrder
        type(PointProcessGenerator) :: spikesGenerator
        type(CharacterMatrix) :: SynapsesOut
        integer, dimension(:), allocatable :: indicesOfSynapsesOnTarget
        type(SynapsePointer), dimension(:), allocatable :: transmitSpikesThroughSynapses
        real(wp) :: axonStimModulation, nerveLength
        real(wp) :: stimulusMeanFrequency_Hz, stimulusPulseDuration_ms, stimulusIntensity_mA
        real(wp) :: stimulusStop_ms, stimulusStart_ms, stimulusModulationStart_ms
        real(wp) :: stimulusModulationStop_ms
        real(wp) :: position_mm, initialFiringRate

        contains
            procedure :: atualizeAfferentUnit
            procedure :: dVdt
            procedure :: transmitSpikes
            procedure :: atualizeDelay
            procedure :: addCompartmentSpike
            procedure :: createStimulus
            procedure :: atualizeCompartments
            procedure :: reset

    end type AfferentUnit

    interface AfferentUnit
        module procedure init_AfferentUnit
    end interface AfferentUnit

    contains

        type(AfferentUnit) function init_AfferentUnit(conf, pool, muscle, index)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pool**: string with Motor unit pool to which the motor
            !     unit belongs.

            !     + **muscle**: 

            !     + **index**: integer corresponding to the motor unit order in
            !     the pool, according to the Henneman's principle (size principle).
            ! '''
            class(Configuration), intent(in), target :: conf
            character(len = 6), intent(in) :: pool
            character(len = 6), intent(in) :: muscle
            integer, intent(in) :: index
            character(len = 80) :: paramTag, paramChar, paramTag2
            real(wp) :: paramReal
            integer :: i
            character(len = 9):: compKind
            real(wp), dimension(:), allocatable :: gCoupling_muS
            real(wp) :: rAxis1, rAxis2, cytR
            real(wp), dimension(:), allocatable :: gLeak, capacitance_nF, EqPot, IPump, compLength
            real(wp), dimension(:,:), allocatable :: GC, GL
            real(wp) :: dynamicNerveLength, distanceToTerminal, delayLength
            integer :: simDuration, NumberOfAxonNodes
            real(wp) :: randNumber
            integer :: N
            
            init_AfferentUnit%position_mm = 0.0
            ! ## Configuration object with the simulation parameters.
            init_AfferentUnit%conf => conf

            init_AfferentUnit%timeStep_ms = init_AfferentUnit%conf%timeStep_ms
            init_AfferentUnit%timeStepByTwo_ms = init_AfferentUnit%conf%timeStepByTwo_ms
            init_AfferentUnit%timeStepBySix_ms = init_AfferentUnit%conf%timeStepBySix_ms
            
            init_AfferentUnit%neuronKind = muscle

            init_AfferentUnit%muscle = muscle
            ! # Neural compartments

            init_AfferentUnit%pool = pool
            paramTag2 = trim(pool) // '-' // trim(muscle)
            paramTag = 'NumberAxonNodes'
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)paramReal
            NumberOfAxonNodes = nint(paramReal)


            ! ## Integer corresponding to the motor unit order in the pool, according to the Henneman's principle (size principle).
            init_AfferentUnit%index = index
            ! ## Dictionary of Compartment of the Motor Unit.
            allocate(init_AfferentUnit%Compartments(NumberOfAxonNodes*2))

            do i = 1, NumberOfAxonNodes, 2
                compKind = 'internode'
                init_AfferentUnit%Compartments(i) = Compartment(compKind, init_AfferentUnit%conf, &
                                                            init_AfferentUnit%pool, init_AfferentUnit%index, &
                                                            init_AfferentUnit%neuronKind)
                compKind = 'node'
                init_AfferentUnit%Compartments(i+1) = Compartment(compKind, init_AfferentUnit%conf, &
                                                            init_AfferentUnit%pool, init_AfferentUnit%index, &
                                                            init_AfferentUnit%neuronKind)
            end do
            ! ## Number of compartments.
            init_AfferentUnit%compNumber = size(init_AfferentUnit%Compartments)
            ! ## Value of the membrane potential, in mV, that is considered a spike.
            if (init_AfferentUnit%compNumber>0) then
                paramTag = 'threshold'
                paramTag2 = trim(pool) // '-' // trim(muscle)
                paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
                read(paramChar, *)init_AfferentUnit%threshold_mV
            else
                init_AfferentUnit%threshold_mV = 0.0
            end if
            ! ## Vector with membrane potential,in mV, of all compartments. 
            
            allocate(init_AfferentUnit%v_mV(init_AfferentUnit%compNumber))
            ! ## Vector with the last instant of spike of all compartments. 
            allocate(init_AfferentUnit%tSpikes(init_AfferentUnit%compNumber))
            init_AfferentUnit%tSpikes(:) = 0.0


            allocate(gCoupling_muS(init_AfferentUnit%compNumber))
            
            paramTag = 'cytR'
            paramChar =  init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)cytR

            do i = 1, init_AfferentUnit%compNumber - 1
                rAxis1 = (cytR * init_AfferentUnit%Compartments(i)%length_mum) / &
                (pi * (init_AfferentUnit%Compartments(i)%diameter_mum/2)**2)
                rAxis2 = (cytR * init_AfferentUnit%Compartments(i+1)%length_mum) / &
                (pi * (init_AfferentUnit%Compartments(i+1)%diameter_mum/2)**2)    
    
                gCoupling_muS(i) = 200 / (rAxis1 + rAxis2)
            end do
            if (init_AfferentUnit%compNumber>0) then
                gCoupling_muS(init_AfferentUnit%compNumber) = 0.0
            end if
            
            allocate(gLeak(init_AfferentUnit%compNumber))
            allocate(capacitance_nF(init_AfferentUnit%compNumber))
            allocate(EqPot(init_AfferentUnit%compNumber))
            allocate(IPump(init_AfferentUnit%compNumber))
            allocate(compLength(init_AfferentUnit%compNumber))
            
            do i = 1, init_AfferentUnit%compNumber                          
                capacitance_nF(i) = init_AfferentUnit%Compartments(i)%capacitance_nF
                gLeak(i) = init_AfferentUnit%Compartments(i)%gLeak_muS
                EqPot(i) = init_AfferentUnit%Compartments(i)%EqPot_mV
                IPump(i) = init_AfferentUnit%Compartments(i)%IPump_nA
                compLength(i) = init_AfferentUnit%Compartments(i)%length_mum
                init_AfferentUnit%v_mV(i) = init_AfferentUnit%Compartments(i)%EqPot_mV
            end do
            
            ! ## Vector with  the inverse of the capacitance of all compartments. 
            
            if (size(capacitance_nF)>0) then
                init_AfferentUnit%capacitanceInv = 1.0 / capacitance_nF
            else
                allocate(init_AfferentUnit%capacitanceInv(init_AfferentUnit%compNumber))
            end if
            
            ! ## Vector with current, in nA,  of each compartment coming from other elements of the model. For example 
            ! ## from ionic channels and synapses.       
            
            allocate(init_AfferentUnit%iIonic(init_AfferentUnit%compNumber))
            
            ! ## Vector with the current, in nA, injected in each compartment.
            allocate(init_AfferentUnit%iInjected(init_AfferentUnit%compNumber))
            allocate(GC(init_AfferentUnit%compNumber,init_AfferentUnit%compNumber))
            allocate(GL(init_AfferentUnit%compNumber,init_AfferentUnit%compNumber))
            
            if (init_AfferentUnit%compNumber>0) then
                init_AfferentUnit%iInjected(:) = 0.0
                init_AfferentUnit%iIonic(:) = 0.0
                GC(:,:) = 0.0
                GL(:,:) = 0.0
            end if

            do i = 1, init_AfferentUnit%compNumber
                if (i == 1) then
                    GC(i,i:i+1) = [-gCoupling_muS(i), gCoupling_muS(i)] 
                else if (i == init_AfferentUnit%compNumber - 1) then
                    GC(i,i-1:i) = [gCoupling_muS(i-1), -gCoupling_muS(i-1)]  
                else
                    GC(i,i-1:i+1) = [gCoupling_muS(i-1), -gCoupling_muS(i-1)-gCoupling_muS(i), &
                                     gCoupling_muS(i)]                    
                end if
                GL(i,i) = gLeak(i)
            end do
            
            ! ## Matrix of the conductance of the motoneuron. Multiplied by the vector self.v_mV,
            ! ## results in the passive currents of each compartment.
            if (init_AfferentUnit%compNumber > 0) then
                init_AfferentUnit%G = GC + GL
                init_AfferentUnit%EqCurrent_nA = matmul(-GL, EqPot) + IPump 
                ! ## index of the last compartment.
            else
                allocate(init_AfferentUnit%G(0,0))
                allocate(init_AfferentUnit%EqCurrent_nA(0))
            end if    
            init_AfferentUnit%lastCompIndex = init_AfferentUnit%compNumber
            
            ! ## Refractory period, in ms, of the motoneuron.
            paramTag = 'AFRefPer'
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%AFRefPer_ms
            
            ! # delay
            ! ## String with type of the nerve. It can be PTN (posterior tibial nerve) or CPN
            ! ## (common peroneal nerve).
            if (init_AfferentUnit%muscle == 'SOL' .or. init_AfferentUnit%muscle == 'MG' .or. &
                init_AfferentUnit%muscle == 'LG') then
                    init_AfferentUnit%nerve = 'PTN'
            else if (init_AfferentUnit%muscle == 'TA') then
                init_AfferentUnit%nerve = 'CPN'
            end if

            ! ## AxonDelay object of the motor unit.
            if (NumberOfAxonNodes == 0) then
                dynamicNerveLength = 0.0
            else
                dynamicNerveLength = sum(compLength) * 1e-6
            end if
            
            paramTag = 'nerveLength_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%nerveLength

            ! ## Distance, in m, of the stimulus position to the terminal. 
            
            paramTag = 'stimDistToTerm_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)distanceToTerminal 
            
            init_AfferentUnit%stimulusPositiontoTerminal = init_AfferentUnit%nerveLength - distanceToTerminal  

            ! ##Frequency threshold of the afferent to th proprioceptor input
            paramTag = 'frequencyThreshold'
            paramTag2 = trim(pool) // '-' // trim(muscle)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%frequencyThreshold_Hz
            
            delayLength =  init_AfferentUnit%nerveLength - dynamicNerveLength

            if (init_AfferentUnit%stimulusPositiontoTerminal < delayLength) then
                paramTag = trim(pool) // '-' // trim(init_AfferentUnit%muscle)
                init_AfferentUnit%Delay = AxonDelay(conf, init_AfferentUnit%nerve, paramtag, delayLength, &
                init_AfferentUnit%stimulusPositiontoTerminal, index)
                init_AfferentUnit%stimulusCompartment = -1
            else
                paramTag = trim(pool) // '-' // trim(init_AfferentUnit%muscle)
                init_AfferentUnit%Delay = AxonDelay(conf, init_AfferentUnit%nerve, paramTag, delayLength, -1.0_wp, index)
                init_AfferentUnit%stimulusCompartment = init_AfferentUnit%compNumber
            end if
            ! # Nerve stimulus function    
            paramTag = 'stimFrequency_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusMeanFrequency_Hz
            
            paramTag = 'stimPulseDuration_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusPulseDuration_ms

            paramTag = 'stimIntensity_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusIntensity_mA

            paramTag = 'stimStart_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusStart_ms
            
            paramtag = 'stimStop_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusStop_ms
            
            paramTag = 'stimModulationStart_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusModulationStart_ms

            paramTag = 'stimModulationStop_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%stimulusModulationStop_ms

            paramTag = 'stimModulation_' // trim(init_AfferentUnit%nerve)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, index)
            read(paramChar, *)init_AfferentUnit%axonStimModulation
            ! ## Vector with the nerve stimulus, in mA.
            simDuration = nint(init_AfferentUnit%conf%simDuration_ms/init_AfferentUnit%conf%timeStep_ms)
            allocate(init_AfferentUnit%nerveStimulus_mA(simDuration))
            
            call init_AfferentUnit%createStimulus()
            
            ! ## Vector with the instants of spikes at the last compartment.
            if (allocated(init_AfferentUnit%lastCompSpikeTrain)) then
                deallocate(init_AfferentUnit%lastCompSpikeTrain)
            end if
            ! ## Vector with the instants of spikes at the terminal.
            
            if (allocated(init_AfferentUnit%terminalSpikeTrain)) then
                deallocate(init_AfferentUnit%terminalSpikeTrain)
            end if

            ! Initial firing rate
            
            randNumber = randn()
            init_AfferentUnit%initialFiringRate = 5.0 + 0.5*randNumber
            if (init_AfferentUnit%initialFiringRate < 0.0) then
                init_AfferentUnit%initialFiringRate = 0.0
            end if
            
            paramTag = 'GammaOrder_' // trim(init_AfferentUnit%pool) // '-' // trim(init_AfferentUnit%muscle)
            paramChar = init_AfferentUnit%conf%parameterSet(paramTag, paramTag2, 0)
            read(paramChar, *)paramReal
            
            init_AfferentUnit%GammaOrder = nint(paramReal)
            ! ## A PointProcessGenerator object, corresponding the generator of
            ! ## spikes of the neural tract unit.   
            init_AfferentUnit%spikesGenerator = PointProcessGenerator(index) 
            
            
            ! ## Build synapses       
            
            init_AfferentUnit%SynapsesOut = CharacterMatrix()

            if (allocated(gCoupling_muS)) deallocate(gCoupling_muS) 
            if (allocated(gLeak)) deallocate(gLeak)
            if (allocated(capacitance_nF)) deallocate(capacitance_nF)
            if (allocated(EqPot)) deallocate(EqPot)
            if (allocated(IPump)) deallocate(IPump)            
            if (allocated(compLength)) deallocate(compLength)           
            if (allocated(GC)) deallocate(GC)           
            if (allocated(GL)) deallocate(GL)           

            

        end function 

        subroutine atualizeAfferentUnit(self, t, proprioceptorFR) 
            ! '''
            ! Atualize the dynamical and nondynamical (delay) parts of the motor unit.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **proprioceptorFR**: proprioceptor firing rate, in Hz.
            ! ''' 
            class(AfferentUnit), intent(inout) :: self
            real(wp), intent(in) :: t, proprioceptorFR

            call self%spikesGenerator%atualizeGenerator(t, proprioceptorFR, self%GammaOrder)
           ! DISPTST1(5, proprioceptorFR)
            if (allocated(self%spikesGenerator%points)) then 
                if (abs(t - self%spikesGenerator%points(size(self%spikesGenerator%points))) < 1e-3) then
                    call self%Delay%addSpinalSpike(t)
                end if
            end if

            if (self%compNumber > 0) then
                call self%atualizeCompartments(t)
            end if

            call self%atualizeDelay(t)
        end subroutine

        subroutine atualizeCompartments(self, t)
            ! '''
            ! Atualize all neural compartments.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            ! '''
            class(AfferentUnit), intent(inout) :: self
            real(wp), intent(in) :: t
            real(wp), dimension(self%compNumber) :: k1, k2, k3, k4
            real(wp) :: vmax, vmin
            integer :: i

            vmin = -30.0
            vmax = 120.0      

            k1 = self%dVdt(t, self%v_mV)        
            k2 = self%dVdt(t + self%conf%timeStepByTwo_ms, self%v_mV + self%conf%timeStepByTwo_ms * k1)
            k3 = self%dVdt(t + self%conf%timeStepByTwo_ms, self%v_mV + self%conf%timeStepByTwo_ms * k2)
            k4 = self%dVdt(t + self%conf%timeStep_ms, self%v_mV + self%conf%timeStep_ms * k3)
            
            self%v_mV = self%v_mV + self%conf%timeStepBySix_ms * (k1+ 2*k2 + 2*k3 + k4)
                    
            self%v_mV = merge(self%v_mV, vmax, self%v_mV.lt.vmax)
            self%v_mV = merge(self%v_mV, vmin, self%v_mV.gt.vmin)

            do i = 1, self%compNumber
                if ((self%v_mV(i) > self%threshold_mV).and.(t-self%tSpikes(i) > self%AFRefPer_ms)) then 
                    call self%addCompartmentSpike(t, i)
                end if
            end do  
        end subroutine
        
        
        function dVdt(self, t, V)
            ! '''
            ! Compute the potential derivative of all compartments of the motor unit.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **V**: Vector with the current potential value of all neural
            !     compartments of the motor unit.
            
            ! \f{equation}{
            !     \frac{dV}{dt} = (I_{active} + GV+ I_{inj} + I_{eq})C_inv   
            ! }
            ! where all the variables are vectors with the number of elements equal
            ! to the number of compartments and \f$G\f$ is the conductance matrix built
            ! in the compGCouplingMatrix function.
            ! '''
            class(AfferentUnit), intent(inout) :: self
            real(wp), intent(in) :: t
            real(wp), intent(in), dimension(self%compNumber) :: V
            real(wp), dimension(self%compNumber) :: dVdt
            integer :: i

            do i = 1, self%compNumber
                self%iIonic(i) = self%Compartments(i)%computeCurrent(t, V(i))
            end do
                
            dVdt =  (self%iIonic + matmul(self%G, V) + self%iInjected + self%EqCurrent_nA) * self%capacitanceInv
        end function
        
        
        subroutine addCompartmentSpike(self, t, comp)
            ! '''
            ! When the soma potential is above the threshold a spike is added tom the soma.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **comp**: integer with the compartment index.
            ! '''
            class(AfferentUnit), intent(inout) :: self
            real(wp), intent(in) :: t
            integer, intent(in) :: comp
            integer :: i, j    

            self%tSpikes(comp) = t
            
            if (comp == self%lastCompIndex) then
                call AddToList(self%lastCompSpikeTrain, t)
                call self%Delay%addSpinalSpike(t)
            end if
            
            
            do i = 1, size(self%Compartments(comp)%Channels) 
                do j = 1, size(self%Compartments(comp)%Channels(i)%condState)
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
            class(AfferentUnit), intent(inout) :: self
            real(wp), intent(in) :: t
            integer :: simDuration
            
            !DISPTST2(5,size(self%Delay%orthodromicSpikeTrainStimTerminal), self%Delay%indexOrthodromicSpikeStimTerminal )
			if (allocated(self%Delay%orthodromicSpikeTrainStimTerminal)) then
                if(size(self%Delay%orthodromicSpikeTrainStimTerminal) >= self%Delay%indexOrthodromicSpikeStimTerminal) then !tst!
				if (abs(t - self%Delay%orthodromicSpikeTrainStimTerminal(self%Delay%indexOrthodromicSpikeStimTerminal)) < 1e-3) then
                    call AddToList(self%terminalSpikeTrain,t)
					!DISPTST2(5,size(self%Delay%orthodromicSpikeTrainStimTerminal), self%Delay%indexOrthodromicSpikeStimTerminal )
                    self%Delay%indexOrthodromicSpikeStimTerminal = self%Delay%indexOrthodromicSpikeStimTerminal + 1 
                    call self%transmitSpikes(t)
				end if 		!tst!
                end if
            end if            

            if (self%stimulusCompartment == -1) then
                simDuration = nint(t/self%conf%timeStep_ms) + 1
                call self%Delay%atualizeStimulus(t, self%nerveStimulus_mA(simDuration))
            end if
        end subroutine

        subroutine transmitSpikes(self, t)
            ! '''
            ! - Inputs:
            !     + **t**: current instant, in ms.
            ! '''
            class(AfferentUnit), intent(inout) :: self
            real(wp), intent(in) :: t
            integer :: i
            
            if (allocated(self%indicesOfSynapsesOnTarget)) then
                do i = 1, size(self%indicesOfSynapsesOnTarget)
                    call self%transmitSpikesThroughSynapses(i)%synapse%receiveSpike(t, self%indicesOfSynapsesOnTarget(i))
                end do
            end if
        end subroutine

        subroutine createStimulus(self)
            ! '''
            ! '''
            class(AfferentUnit), intent(inout) :: self
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

        subroutine reset(self)
            ! '''

            ! '''
            class(AfferentUnit), intent(inout) :: self
            integer :: i
            
            do i = 1, size(self%Compartments)
                self%v_mV(i) = self%Compartments(i)%EqPot_mV
            end do
            
            call self%Delay%reset()
            
            if (self%compNumber > 0) then
                self%tSpikes(:) = 0.0
            end if
            
            if (allocated(self%lastCompSpikeTrain)) deallocate(self%lastCompSpikeTrain)
            ! ## Vector with the instants of spikes at the terminal.
            if (allocated(self%terminalSpikeTrain)) deallocate(self%terminalSpikeTrain)
            
            call self%spikesGenerator%reset()
            
            
        end subroutine


        
end module AfferentUnitClass
