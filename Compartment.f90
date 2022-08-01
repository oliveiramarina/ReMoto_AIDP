! # '''
! #     Neuromuscular simulator in Fortran.
! #     Copyright (C) 2018  Renato Naville Watanabe

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

module CompartmentClass
    ! '''
    ! Class that implements a neural compartment. For now it is implemented
    ! *dendrite* and *soma*.
    ! '''
    use ConfigurationClass
    use ChannelConductanceClass    
    use SynapseClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4*atan(1.0_wp)    
    public :: Compartment

    type Compartment
        character(len = 9) :: compKind
        type(ChannelConductance), dimension(:), allocatable :: Channels ! ## List of ChannelConductance objects in the Compartment.
        character(len = 2)::  neuronKind
        integer :: index, numberChannels
        type(Configuration), pointer :: conf
        real(wp) :: length_mum, diameter_mum, capacitance_nF
        real(wp) :: EqPot_mV, IPump_nA, gLeak_muS
        type(Synapse), dimension(2) :: SynapsesIn
        
        contains

            procedure :: reset
            procedure :: computeCurrent

    end type Compartment

    interface Compartment
        module procedure init_Compartment
    end interface Compartment

    contains

    type(Compartment) function init_Compartment(compKind, conf, pool, index, neuronKind)
        ! '''
        ! Constructor

        ! - Inputs:
        !     + **kind**: The kind of compartment. For now, it can be *soma*, *dendrite*, 
        !     *node* or *internode*.

        !     + **conf**: Configuration object with the simulation parameters.

        !     + **pool**: string with Motor unit pool to which the motor unit belongs.

        !     + **index**: integer corresponding to the motor unit order in the pool, according to 
        !     the Henneman's principle (size principle).

        !     + **neuronKind**: string with the type of the motor unit. It can be *S* (slow), *FR* (fast and resistant), 
        !     and *FF* (fast and fatigable).
        ! '''
        character(*), intent(in) :: compKind
        class(Configuration), intent(in), target :: conf
        character(*), intent(in) :: pool
        integer, intent(in) :: index
        character(*), intent(in) ::  neuronKind
        character(len=80) :: paramTag, paramChar
        real(wp) :: area_cm2, specifRes_Ohmcm2, membraneCapacitance
        character(len = 6) :: channelKind
        character(len = 10) :: synapseKind
        
        
        init_Compartment%conf => conf
        ! ## String with the type of the motor unit. It can be *S* (slow), *FR* (fast and resistant), 
        ! ## and *FF* (fast and fatigable).
        init_Compartment%neuronKind = neuronKind
        ! ## List of summed synapses (see Lytton, 1996) that the Compartment receive from other neural components.
        synapseKind = 'excitatory'
        init_Compartment%SynapsesIn(1) = Synapse(conf, pool, index, compKind, synapseKind, neuronKind)
        synapseKind = 'inhibitory'
        init_Compartment%SynapsesIn(2) = Synapse(conf, pool, index, compKind, synapseKind, neuronKind)
        
        ! ## The kind of compartment. For now, it can be *soma* or *dendrite*, *node* or *internode*..
        init_Compartment%compKind = compKind
        
        ! ## Integer corresponding to the motor unit order in the pool, according to 
        ! ## the Henneman's principle (size principle).
        init_Compartment%index = index
        
        ! ## Length of the compartment, in \f$\mu\f$m.
        paramTag = 'l@' // trim(init_Compartment%compKind)
        paramChar =  init_Compartment%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Compartment%length_mum
        ! ## Diameter of the compartment, in \f$\mu\f$m.
        paramTag = 'd@' // trim(init_Compartment%compKind)
        paramChar =  init_Compartment%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Compartment%diameter_mum
        area_cm2 = init_Compartment%length_mum * pi * init_Compartment%diameter_mum * 1e-8
        paramTag = 'res@' // trim(init_Compartment%compKind)
        paramChar = init_Compartment%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)specifRes_Ohmcm2
        ! ## Capacitance of the compartment, in nF.
        paramTag = 'membCapac'
        paramChar = init_Compartment%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)membraneCapacitance
        init_Compartment%capacitance_nF = membraneCapacitance * area_cm2 * 1e3
        
        ! ## Equilibrium potential, in mV.
        paramTag = 'EqPot@' // trim(init_Compartment%compKind)
        paramChar = init_Compartment%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Compartment%EqPot_mV

        ! ## Pump current in the compartment, in nA.
        paramTag = 'IPump@' // trim(init_Compartment%compKind)
        paramChar = init_Compartment%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Compartment%IPump_nA

        ! ## Leak conductance of the compartment, in \f$\mu\f$S.
        init_Compartment%gLeak_muS = (1e6 * area_cm2) / specifRes_Ohmcm2

        init_Compartment%numberChannels = 0    
        if (trim(init_Compartment%compKind).eq.'soma') then
            allocate(init_Compartment%Channels(3))
            channelKind = 'Kf'
            init_Compartment%Channels(1) = ChannelConductance(channelKind, init_Compartment%conf,&
                                                              area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1                                                            
            channelKind = 'Ks'
            init_Compartment%Channels(2) = ChannelConductance(channelKind, init_Compartment%conf,&
                                                              area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
            channelKind = 'Na'
            init_Compartment%Channels(3) = ChannelConductance(channelKind, init_Compartment%conf,&
                                                              area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1                                                              
        !TODO: else if (kind == 'dendrite'):
        !    pass
		!Bostock 2012
        else if (trim(init_Compartment%compKind).eq.'node') then
            allocate(init_Compartment%Channels(4))
            channelKind = 'Na'
            init_Compartment%Channels(1) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
            channelKind = 'Nap'    
            init_Compartment%Channels(2) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
            channelKind = 'Kf'
            init_Compartment%Channels(3) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
            channelKind = 'KsAxon'
            init_Compartment%Channels(4) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
        else if (trim(init_Compartment%compKind).eq.'internode') then
            allocate(init_Compartment%Channels(3))
            channelKind = 'Kf'
            init_Compartment%Channels(1) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
            channelKind = 'KsAxon'
            init_Compartment%Channels(2) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
            channelKind = 'H'
            init_Compartment%Channels(3) = ChannelConductance(channelKind, conf, area_cm2, pool, neuronKind, compKind, index)
            init_Compartment%numberChannels = init_Compartment%numberChannels + 1
        end if

        ! ## Integer with the number of ionic channels.
        

    end function

    real(wp) function computeCurrent(self, t, V_mV) result(current)
        ! '''
        ! Computes the active currents of the compartment. Active currents are the currents from the ionic channels
        ! and from the synapses.

        ! - Inputs:
        !     + **t**: current instant, in ms.

        !     + **V_mV**: membrane potential, in mV.
        ! '''
        class(Compartment), intent(inout):: self
        real(wp), intent(in) :: t, V_mV
        integer :: i
        
        
        current = 0.0
        
        
        if (self%SynapsesIn(1)%numberOfIncomingSynapses > 0) current = current + self%SynapsesIn(1)%computeCurrent(t, V_mV)
        
        if (self%SynapsesIn(2)%numberOfIncomingSynapses > 0) current = current + self%SynapsesIn(2)%computeCurrent(t, V_mV)
        
        do i = 1, self%numberChannels 
            current = current + self%Channels(i)%computeCurrent(t, V_mV)
        end do
        
    end function

    subroutine reset(self)
        ! '''

        ! '''
        class(Compartment), intent(inout):: self
        integer :: i

        do i = 1, size(self%SynapsesIn)
            call self%SynapsesIn(i)%reset()
        end do
        do i = 1, self%numberChannels
            call self%Channels(i)%reset()
        end do
    end subroutine


end module CompartmentClass



        

        
    
    