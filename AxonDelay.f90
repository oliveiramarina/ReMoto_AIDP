! # '''
! #     Neuromuscular simulator in Fortran.
! #     Copyright (C) 2021  Renato Naville Watanabe
! #							Marina Cardoso de Oliveira

! #     This program is free software: you can redistribute it and/or modify
! #     it under the terms of the GNU General Public License as published by
! #     the Free Software Foundation, either version 3 of the License, or
! #     any later version.

! #     This program is distributed in the hope that it will be useful,
! #     but WITHOUT ANY WARRANTY; without even the implied warranty of
! #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! #     GNU General Public License for more details.

! #     You should have received a copy of the GNU General Public License
! #     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! #     Contact: renato.watanabe@ufabc.edu.br
! # '''
module AxonDelayClass
    ! '''
    ! Class that implements a delay correspondent to the nerve. This class corresponds to the part of the axon that is
    ! modeled with no dynamics. Ideally this  class would not exist and all the axon would be modelled in the motor unit
    ! or sensory class with the proper dynamics. 
    ! '''
    use ConfigurationClass
    use DynamicalArrays
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    public :: AxonDelay

    type AxonDelay
        type(Configuration), pointer :: conf
        integer :: index
        real :: length_m, velocity_m_s, latencyStimulusSpinal_ms
        real(wp) :: latencySpinalTerminal_ms, latencyStimulusTerminal_ms
        real(wp) :: axonSpikeTrain
		type (fila) :: terminalSpikeTrain_fila
		real(wp) :: terminalSpikeTrain_Temp1, terminalSpikeTrain_Temp2
		real(wp), dimension(:), allocatable :: orthodromicSpikeTrain, antidromicSpikeTrain
        real(wp), dimension(:), allocatable :: orthodromicSpikeTrainStimTerminal
        integer :: indexOrthodromicSpike, indexAntidromicSpike, indexOrthodromicSpikeStimTerminal
        real(wp) :: electricCharge_muC, threshold_muC, refractoryPeriod_ms, leakageTimeConstant_ms
        character(len = 6) :: pool

        contains
            procedure :: atualizeStimulus
            procedure :: reset
            procedure :: addTerminalSpike
            procedure :: addAntidromicSpike
            procedure :: addSpinalSpike
    end type AxonDelay

    interface AxonDelay
        module procedure init_AxonDelay
    end interface AxonDelay

    contains

    type(AxonDelay) function init_AxonDelay(conf, nerve, pool, length, stimulusPositiontoTerminal, index)
        ! '''
        ! Constructor


        ! - Inputs:
        !     + **conf**: Configuration object with the simulation parameters.

        !     + **nerve**: string with type of the nerve. It can be *PTN* 
        !     (posterior tibial nerve) or *CPN* (common peroneal nerve).

        !     + **pool**: string with Motor unit pool to which the motor unit belongs.

        !     + **length**:  float, length of the part of the nerve that is
        !     modelled as a delay, in m.

        !     + **stimulusPositiontoTerminal**: float, distance, in m, of the stimulus
        !     position to the terminal, in m. If -1, indicates it is not in the Delay.

        !     + **index**: integer corresponding to the motor unit order in the pool, according to 
        !     the Henneman's principle (size principle).
        ! '''
        class(Configuration), intent(in), target :: conf        
        character(len = 3), intent(in) :: nerve 
        character(len = 6), intent(in) :: pool
        real(wp), intent(in) :: length, stimulusPositiontoTerminal
        integer, intent(in) :: index
        character(len=80) :: paramTag, paramChar

        init_AxonDelay%conf => conf
        init_AxonDelay%pool = pool
        ! ## Integer corresponding to the motor unit order in the pool, according to 
        ! ## the Henneman's principle (size principle).
        init_AxonDelay%index = index

        ! ## Length, in m, of the part of the nerve that is modelled as a delay.
        init_AxonDelay%length_m = length
        ! ## Velocity of conduction, in m/s, of the part of the nerve that is not modelled as a delay.     
        paramTag = 'axonDelayCondVel'
		init_AxonDelay%velocity_m_s = init_AxonDelay%conf%getparameterGauss(paramTag, pool, index)
						 
        ! ## time, in ms, that the signal takes to travel between the stimulus and the spinal cord.        
        init_AxonDelay%latencyStimulusSpinal_ms = nint((init_AxonDelay%length_m - &
            stimulusPositiontoTerminal)/init_AxonDelay%velocity_m_s*1000.0/init_AxonDelay%conf%timeStep_ms) *&
            init_AxonDelay%conf%timeStep_ms
        ! ## time, in ms, that the signal takes to travel between the spinal cord and the terminal.
        init_AxonDelay%latencySpinalTerminal_ms = nint((init_AxonDelay%length_m)/init_AxonDelay%velocity_m_s/init_AxonDelay%conf%timeStep_ms*1000.0) * init_AxonDelay%conf%timeStep_ms
        ! ## time, in ms, that the signal takes to travel between the stimulus and the terminal.
        init_AxonDelay%latencyStimulusTerminal_ms = nint((stimulusPositiontoTerminal)/&
        init_AxonDelay%velocity_m_s/init_AxonDelay%conf%timeStep_ms*1000.0) * init_AxonDelay%conf%timeStep_ms
			
        ! ## Float with instant, in ms, of the last spike in the terminal. 
        init_AxonDelay%axonSpikeTrain = -1e6
		init_AxonDelay%terminalSpikeTrain_Temp1 = -1e6
		init_AxonDelay%terminalSpikeTrain_Temp2 = -1e6
        
        init_AxonDelay%indexOrthodromicSpike = 1
        init_AxonDelay%indexOrthodromicSpikeStimTerminal = 1
        init_AxonDelay%indexAntidromicSpike = 1

        init_AxonDelay%electricCharge_muC = 0
        
        
        paramTag = 'axonDelayThreshold'
        paramChar = init_AxonDelay%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_AxonDelay%threshold_muC
        
        paramTag = 'axonDelayRefPeriod_' // nerve
        paramChar = init_AxonDelay%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_AxonDelay%refractoryPeriod_ms
        
        paramTag = 'axonDelayLeakTimeConstant'
        paramChar = init_AxonDelay%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_AxonDelay%leakageTimeConstant_ms
    end function

    subroutine addTerminalSpike(self, t, latency)
        ! '''
        ! Indicates to the AxonDelay object that a spike has occurred in the Terminal.

        ! - Inputs:
        !     + **t**: current instant, in ms.

        !     + **latency**: time elapsed until the spike take effect, in ms.
        ! '''
        class(AxonDelay), intent(inout) :: self
        real(wp), intent(in) :: t, latency

			call self%terminalSpikeTrain_fila%enqueue(t + latency)
			call AddToList(self%orthodromicSpikeTrainStimTerminal, t + latency)
				
    end subroutine

    subroutine addSpinalSpike(self, t)
        ! '''
        ! Indicates to the AxonDelay object that a spike has occurred in the last
        ! dynamical compartment of the motor unit.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(AxonDelay), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp) :: timeInSpinalCord

        timeInSpinalCord = t + self%latencyStimulusSpinal_ms
        call AddToList(self%orthodromicSpikeTrain, timeInSpinalCord)
    end subroutine

    subroutine addAntidromicSpike(self, t)
        ! '''
        ! Indicates to the AxonDelay object that a spike has occurred in the last
        ! dynamical compartment of the motor unit.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''

        class(AxonDelay), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp) :: timeInSpinalCord

        timeInSpinalCord = t + self%latencyStimulusSpinal_ms

        call AddToList(self%antidromicSpikeTrain, timeInSpinalCord)
        
    end subroutine

    subroutine atualizeStimulus(self, t, stimulus)
        ! '''
        ! Inputs:
        ! 
        ! + **t**: float, instant, in ms
        ! 
        ! + **stimulus**:float , stimulus inthe nerve, in mA
        ! 
        ! '''
        class(AxonDelay), intent(inout) :: self
        real(wp), intent(in) :: t, stimulus

        
        if ((t - self%axonSpikeTrain) > self%refractoryPeriod_ms) then
            self%electricCharge_muC = (stimulus * self%conf%timeStep_ms +&
                                   self%electricCharge_muC * &
                                   exp(-self%conf%timeStep_ms /self%leakageTimeConstant_ms))
        
            if (self%electricCharge_muC >= self%threshold_muC) then
                self%electricCharge_muC = 0
                call self%addTerminalSpike(t, self%latencyStimulusTerminal_ms)
                call self%addAntidromicSpike(t)
                self%axonSpikeTrain = t
            end if
            if (allocated(self%orthodromicSpikeTrain)) then
                if (self%indexOrthodromicSpike <= size(self%orthodromicSpikeTrain)) then
                    if (t > self%orthodromicSpikeTrain(self%indexOrthodromicSpike)) then
                        if (allocated(self%antidromicSpikeTrain)) then 
                            if (self%indexAntidromicSpike <= size(self%antidromicSpikeTrain)) then 
                                if (abs(self%orthodromicSpikeTrain(self%indexOrthodromicSpike) - &
                                        self%antidromicSpikeTrain(self%indexAntidromicSpike)) < self%latencyStimulusSpinal_ms) then
                                    self%indexOrthodromicSpike = self%indexOrthodromicSpike + 1
                                    self%indexAntidromicSpike = self%indexAntidromicSpike + 1
                                else
                                    self%electricCharge_muC = 0
                                    call self%addTerminalSpike(t, self%latencyStimulusTerminal_ms)
                                    self%axonSpikeTrain = t
                                    self%indexOrthodromicSpike = self%indexOrthodromicSpike + 1
								end if
                            end if
                        else
                            self%electricCharge_muC = 0
                            call self%addTerminalSpike(t, self%latencyStimulusTerminal_ms)
                            self%axonSpikeTrain = t
                            self%indexOrthodromicSpike = self%indexOrthodromicSpike + 1
						end if                        
                    end if
                end if
            end if
        end if
    end subroutine

    subroutine reset(self)
        ! '''

        ! '''
        class(AxonDelay), intent(inout) :: self
		real(wp):: temp
        
        self%electricCharge_muC = 0
		do while (self%terminalSpikeTrain_fila%n > 0) 
			
			call self%terminalSpikeTrain_fila%pop()
		
		end do
			
        if (allocated(self%terminalSpikeTrain_fila%prioritytime)) deallocate(self%terminalSpikeTrain_fila%prioritytime) 
			
		self%axonSpikeTrain = -1e6

        if (allocated(self%orthodromicSpikeTrain)) deallocate(self%orthodromicSpikeTrain)
        if (allocated(self%antidromicSpikeTrain)) deallocate(self%antidromicSpikeTrain)
        if (allocated(self%orthodromicSpikeTrainStimTerminal)) deallocate(self%orthodromicSpikeTrainStimTerminal)

        self%indexOrthodromicSpike = 1
        self%indexOrthodromicSpikeStimTerminal = 1
        self%indexAntidromicSpike = 1
    end subroutine       
end module AxonDelayClass


        
