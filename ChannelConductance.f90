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

module ChannelConductanceClass
    ! '''
    ! Class that implements a model of the ionic Channels in a compartment.
    ! '''
    use ConfigurationClass
    use PulseConductanceStateClass
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    public :: ChannelConductance

    type ChannelConductance
        type(Configuration), pointer :: conf
        character(len = 2) :: neuronKind
        character(len = 6) :: channelKind
        real(wp) :: EqPot_mV, gmax_muS
        character(len = 5) :: stateType
        type(PulseConductanceState), dimension(:), allocatable :: condState 
        integer :: lenStates

        contains
            procedure :: computeCurrent
            procedure :: reset

    end type ChannelConductance

    interface ChannelConductance
        module procedure init_ChannelConductance
    end interface ChannelConductance

    contains

    type(ChannelConductance) function init_ChannelConductance(channelKind, conf, compArea, pool, neuronKind, compKind, index)
        ! '''
        ! Constructor
        
        ! Builds an ionic channel conductance.

        ! -Inputs: 
        !     + **channelKind**: string with the type of the ionic channel. For now it 
        !     can be *Na* (Sodium), *Ks* (slow Potassium), *Kf* (fast Potassium) or 
        !     *Ca* (Calcium).

        !     + **conf**: instance of the Configuration class (see Configuration file).

        !     + **compArea**: float with the area of the compartment that the Channel belongs, in \f$\text{cm}^2\f$.

        !     + **pool**: the pool that this state belongs.

        !     + **neuronKind**: string with the type of the motor unit. It used for 
        !     motoneurons. It can be *S* (slow), *FR* (fast and resistant), and *FF* 
        !     (fast and fatigable).

        !     + **compKind**: The kind of compartment that the Channel belongs. 
        !     For now, it can be *soma*, *dendrite*, *node* or *internode*.

        !     + **index**: the index of the unit that this state belongs.          
        ! '''
        character(*), intent(in) ::  neuronKind
        class(Configuration), intent(in), target :: conf
        real(wp) :: compArea
        character(*), intent(in) :: channelKind, pool
        character(*), intent(in) :: compKind
        integer, intent(in) :: index
        character(len=80) :: paramTag, paramChar
        real(wp) :: channelDensity
        character(len = 2) :: stateType
        
        init_ChannelConductance%conf => conf
        ! ## string with the type of the ionic channel. For now it 
        ! ## can be *Na* (Sodium), *Ks* (slow Potassium), *Kf* (fast Potassium) or 
        ! ## *Ca* (Calcium).
        init_ChannelConductance%channelKind = channelKind
        ! ## Equilibrium Potential of the ionic channel, mV.
        paramTag = 'EqPot_' // trim(channelKind) // '@' // trim(compKind)
        paramChar = init_ChannelConductance%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_ChannelConductance%EqPot_mV
        ! ## Maximal conductance, in \f$\mu\f$S, of the ionic channel. 
        paramTag = 'gmax_' // trim(channelKind) // ':' // trim(pool) // '-' // trim(neuronKind) // '@' // trim(compKind)
        paramChar = init_ChannelConductance%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)channelDensity 
        init_ChannelConductance%gmax_muS = compArea * channelDensity
        !  ## String with type of dynamics of the states. For now it accepts the string pulse.
        paramTag = 'StateType'
        init_ChannelConductance%stateType = init_ChannelConductance%conf%parameterSet(paramTag, pool, index)
        
        if (trim(init_ChannelConductance%stateType) == 'pulse') then
            if (trim(init_ChannelConductance%channelKind).eq.'Kf') then
                allocate(init_ChannelConductance%condState(1))
                stateType = 'n'
                init_ChannelConductance%condState(1) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
            end if
            if (trim(init_ChannelConductance%channelKind).eq.'Ks') then
                allocate(init_ChannelConductance%condState(1))
                statetype = 'q'
                init_ChannelConductance%condState(1) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
            end if
            if (trim(init_ChannelConductance%channelKind).eq.'Na') then 
                allocate(init_ChannelConductance%condState(2))
                stateType = 'm'
                init_ChannelConductance%condState(1) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
                stateType = 'h'
                init_ChannelConductance%condState(2) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
            end if
            !TODO:
            ! if (init_ChannelConductance%kind == 'Ca')
            !     pass  # to be implemented
            if (init_ChannelConductance%channelKind.eq.'Nap') then
                allocate(init_ChannelConductance%condState(1))
                stateType = 'mp'
                init_ChannelConductance%condState(1) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
            end if
            if (init_ChannelConductance%channelKind.eq.'KsAxon') then
                allocate(init_ChannelConductance%condState(1))
                stateType = 's'
                init_ChannelConductance%condState(1) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
            end if
            if (init_ChannelConductance%channelKind.eq.'H') then
                allocate(init_ChannelConductance%condState(1))
                stateType = 'qh'
                init_ChannelConductance%condState(1) = PulseConductanceState(stateType, conf, pool, neuronKind, compKind, index)
            end if
        end if 
        
        ! ## Integer with the number of states in the ionic channel.    
        init_ChannelConductance%lenStates = size(init_ChannelConductance%condState)          
    end function

    real(wp) function computeCurrent(self, t, V_mV) result(current) 
        ! '''
        ! Computes the current generated by the ionic Channel
        
        ! - Inputs:
        !     + **t**: instant in ms.
        !     + **V_mV**: membrane potential of the compartment in mV.
        
        ! - Outputs:
        !     + Ionic current, in nA
        ! '''        
        class(ChannelConductance), intent(inout) :: self
        real(wp), intent(in) :: t, V_mV
        integer :: i    
        


        do i = 1, self%lenStates
            call self%condState(i)%computeStateValue(t) 
        end do
        
        if (trim(self%channelKind)=='Kf') then
            !     It is computed as:

            ! \f{equation}{
            !     g = g_{max}n^4(E_0-V)
            ! \f}
            current = self%gmax_muS*(self%condState(1)%value**4) * (self%EqPot_mV - V_mV)            
        end if
        if (trim(self%channelKind)=='Ks') then
            !It is computed as:

            ! \f{equation}{
            !     g = g_{max}q^2(E_0-V)
            ! \f}
            current = self%gmax_muS*(self%condState(1)%value**2) * (self%EqPot_mV - V_mV)
        end if
        if (trim(self%channelKind)=='Na') then
            !     It is computed as:

            ! \f{equation}{
            !     g = g_{max}m^3h(E_0-V)
            ! \f} 
            current = self%gmax_muS*(self%condState(1)%value**3)*self%condState(2)%value*(self%EqPot_mV - V_mV)        
        end if
        !TODO:
        ! if (init_ChannelConductance%kind == 'Ca')
        
        !end if
        if (self%channelKind=='Nap') then
            !     It is computed as:

            ! \f{equation}{
            !     g = g_{max}m_p^3(E_0-V)
            ! \f}
            current = self%gmax_muS*(self%condState(1)%value**3)*(self%EqPot_mV - V_mV)    
        end if
        if (self%channelKind=='KsAxon') then
            !It is computed as:

            ! \f{equation}{
            !     g = g_{max}s(E_0-V)
            ! \f}
            current = self%gmax_muS * self%condState(1)%value * (self%EqPot_mV - V_mV)    
        end if
        if (self%channelKind=='H') then
            ! It is computed as:

            ! \f{equation}{
            !     g = g_{max}q_h(E_0-V)
            ! \f}
            current = self%gmax_muS*self%condState(1)%value*(self%EqPot_mV - V_mV)    
        end if
    end function 

    subroutine reset(self)
        ! '''

        ! '''
        class(ChannelConductance), intent(inout) :: self
        integer :: i

        do i = 1, self%lenStates 
            call self%condState(i)%reset()
        end do
    end subroutine


end module ChannelConductanceClass





    
        
    

    

   
    
            
    
    
         
    
    
         
        
