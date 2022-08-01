! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2018  Renato Naville Watanabe

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

module PulseConductanceStateClass
    ! '''
    ! Implements the Destexhe pulse approximation of the solution of 
    ! the states of the Hodgkin-Huxley neuron model.
    ! '''
    use ConfigurationClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    public :: PulseConductanceState


    type PulseConductanceState
        type(Configuration), pointer :: conf
        character(len = 2) :: stateKind
        real(wp) :: value
        logical :: state
        real(wp) :: beta_ms1, alpha_ms1, PulseDur_ms
        real(wp) :: AlphaExp, BetaExp
        real(wp) :: endOfPulse_ms
        character(len=12) :: actType


        contains
            procedure :: changeState
            procedure :: computeStateValue
            procedure :: reset
    end type PulseConductanceState

    interface PulseConductanceState
        module procedure init_PulseConductanceState
    end interface PulseConductanceState

    contains

        type(PulseConductanceState) function init_PulseConductanceState(stateKind, conf, pool, neuronKind, compKind, index)
            ! '''
            ! Initializes the pulse conductance state.

            ! Variables:
            !     + **stateKind**: string with type of the state (m, h, n, q).
                
            !     + **conf**:  an instance of the Configuration class with the functions to correctly parameterize the model. See the Configuration class.
                
            !     + **pool**: string with the pool that this state belongs.
                
            !     + **neuronKind**: string with the type of the motor unit. It used for 
            !     motoneurons. It can be *S* (slow), *FR* (fast and resistant), and *FF* 
            !     (fast and fatigable). 

            !     + **compKind**: The kind of compartment that the Channel belongs. 
            !     For now, it can be *soma*, *dendrite*, *node* or *internode*.

            !     + **index**: the index of the unit that this state belongs.                    
            ! ''' 
            character(*), intent(in) :: stateKind, neuronKind
            class(Configuration), intent(in), target :: conf
            character(*), intent(in) :: pool
            character(*), intent(in) :: compKind
            integer, intent(in) :: index
            character(len=80) :: paramTag, paramChar

            init_PulseConductanceState%stateKind = stateKind
            init_PulseConductanceState%value = 0.0
            
            init_PulseConductanceState%state = .false.
            
            init_PulseConductanceState%conf => conf

            paramTag = 'beta_' // trim(stateKind) // ':' // trim(pool) // '-' // trim(neuronKind) // '@' // trim(compKind)
            paramChar = init_PulseConductanceState%conf%parameterSet(paramTag, pool, index)
            read(paramChar, *)init_PulseConductanceState%beta_ms1
            
            paramTag = 'alpha_' // trim(stateKind) // ':' // trim(pool) // '-' // trim(neuronKind) // '@' // trim(compKind)
            paramChar = init_PulseConductanceState%conf%parameterSet(paramTag, pool, index)
            read(paramChar, *)init_PulseConductanceState%alpha_ms1
            
            paramTag = 'PulseDur_' // trim(stateKind)
            paramChar = init_PulseConductanceState%conf%parameterSet(paramTag, pool, index)
            read(paramChar, *)init_PulseConductanceState%PulseDur_ms    
            

            init_PulseConductanceState%AlphaExp = exp(-init_PulseConductanceState%alpha_ms1 * conf%timeStep_ms)
            init_PulseConductanceState%BetaExp = exp(-init_PulseConductanceState%beta_ms1 * conf%timeStep_ms)

            init_PulseConductanceState%endOfPulse_ms = init_PulseConductanceState%PulseDur_ms

            if (trim(init_PulseConductanceState%stateKind).eq.'m') init_PulseConductanceState%actType = 'activation'
            if (trim(init_PulseConductanceState%stateKind).eq.'h') init_PulseConductanceState%actType = 'inactivation'
            if (trim(init_PulseConductanceState%stateKind).eq.'n') init_PulseConductanceState%actType = 'activation'
            if (trim(init_PulseConductanceState%stateKind).eq.'q') init_PulseConductanceState%actType = 'activation'
            if (trim(init_PulseConductanceState%stateKind).eq.'mp') init_PulseConductanceState%actType = 'activation'
            if (trim(init_PulseConductanceState%stateKind).eq.'s') init_PulseConductanceState%actType = 'activation'
            if (trim(init_PulseConductanceState%stateKind).eq.'qh') init_PulseConductanceState%actType = 'inactivation'
        end function init_PulseConductanceState

        subroutine changeState(self, t)
            ! '''
            ! Void function that modify the current situation (true/false)
            ! of the state.

            ! - Inputs:
            !     + **t**: current instant, in ms.
            ! '''
            class(PulseConductanceState), intent(inout) :: self
            real(wp) :: t

            self%state = .not.self%state
            self%endOfPulse_ms = self%PulseDur_ms + t

        end subroutine

        subroutine computeStateValue(self, t)
            ! '''
            ! Compute the state value by using the approximation of Destexhe (1997) to
            ! compute the Hodgkin-Huxley states.

            ! - Input:
            !     + **t**: current instant, in ms.

            ! The value of the state \f$v\f$ is computed according to the following
            ! equation before and after the pulse if theactType is activation:

            ! \f{equation}{
            !     v(t) = v_0\exp[-\beta(t-t_0)]
            ! \f} 

            ! and according to the following equation during the pulse:

            ! \f{equation}{
            !     v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
            ! \f} 
            ! where \f$t_0\f$ is the time at which the pulse changed
            ! the value (on to off or off to on) and \f$v_0\f$ is value
            ! of the state at that time.
            !The value of the state \f$v\f$ is computed according to the following
            ! equation before and after the pulse if the actType is inactivation:

            ! \f{equation}{
            !     v(t) = v_0\exp[-\beta(t-t_0)]
            ! \f} 

            ! and according to the following equation during the pulse:

            ! \f{equation}{
            !     v(t) = 1 + (v_0 - 1)\exp[-\alpha(t-t_0)]
            ! \f} 
            ! where \f$t_0\f$ is the time at which the pulse changed
            ! the value (on to off or off to on) and \f$v_0\f$ is value
            ! of the state at that time.
            ! '''
            
            class(PulseConductanceState), intent(inout) :: self
            real(wp) :: t

            if (trim(self%actType).eq.'activation') then
                if (.not.self%state) then
                    self%value = self%value*self%BetaExp
                else if (t > self%endOfPulse_ms) then
                    call self%changeState(t)
                    self%value = self%value*self%BetaExp
                else
                    self%value = (self%value - 1) * self%AlphaExp + 1
                end if
            else
                if (.not.self%state) then
                    self%value = (self%value - 1) * self%AlphaExp + 1
                else if (t > self%endOfPulse_ms) then
                        call self%changeState(t)
                        self%value = (self%value - 1) * self%AlphaExp + 1
                else
                     self%value = self%value*self%BetaExp 
                end if 
            end if
            
        end subroutine


        subroutine reset(self)
            ! '''
            !     reset the pulse state
            ! '''

            class(PulseConductanceState), intent(inout) :: self    
            
            self%value = 0.0
            self%endOfPulse_ms = self%PulseDur_ms
            self%state = .false.
        end subroutine
end module PulseConductanceStateClass