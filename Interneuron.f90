! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2018  Renato Naville Watanabe
!                         Pablo Alejandro

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

module InterneuronClass
    ! '''
    ! Class that implements a motor unit model. Encompasses a motoneuron
    ! and a muscle unit.
    ! '''
    use CompartmentClass
    use ConfigurationClass
    use CharacterMatrixClass
    use SynapsePointerClass
    use DynamicalArrays
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: PI = 4 * atan(1.0_wp)    
    public :: Interneuron

    type Interneuron
        type(Configuration), pointer :: conf    
        character(len = 15) :: pool
        character(len = 2) :: neuronKind
        real(wp) :: tSomaSpike
        real(wp), dimension(:), allocatable :: somaSpikeTrain
        integer :: index, somaIndex, compNumber
        type(Compartment), dimension(:), allocatable :: Compartments 
        real(wp), dimension(:), allocatable :: v_mV,  capacitanceInv
        real(wp), dimension(:), allocatable :: iIonic, iInjected
        real(wp), dimension(:,:), allocatable :: G
        real(wp), dimension(:), allocatable :: EqCurrent_nA
        type(CharacterMatrix) :: SynapsesOut
        integer, dimension(:), allocatable :: indicesOfSynapsesOnTarget
        type(SynapsePointer), dimension(:), allocatable :: transmitSpikesThroughSynapses
        real(wp) :: threshold_mV, position_mm, RefPer_ms

        contains
            procedure :: atualizeInterneuron
            procedure :: transmitSpikes
            procedure :: atualizeCompartments
            procedure :: addSomaSpike
            procedure :: reset

    end type Interneuron

    interface Interneuron
        module procedure init_Interneuron
    end interface Interneuron

    contains

    type(Interneuron) function init_Interneuron(conf, pool, index)
        ! '''
        ! Constructor

        ! - Inputs:
        !     + **conf**: Configuration object with the simulation parameters.

        !     + **pool**: string with Interneuron pool to which the motor
        !     unit belongs.  It can
        !     be *RC* (Renshaw cell), *IaIn* (Ia Interneuron), *IbIn* (Ib Interneuron) and 
        !     *gII*.

        !     + **index**: integer corresponding to the motor unit order in
        !     the pool, according to the Henneman's principle (size principle).
        ! '''

        ! ## Configuration object with the simulation parameters.

        class(Configuration), intent(in), target :: conf    
        character(len = 15), intent(in) :: pool
        integer, intent(in) :: index
        character(len=80) :: paramTag, paramChar
        integer :: i
        character(len = 9):: compKind
        real(wp), dimension(:), allocatable :: gLeak, capacitance_nF, EqPot
        real(wp), dimension(:,:), allocatable :: GL

        init_Interneuron%conf => conf

        init_Interneuron%pool = pool
        
        init_Interneuron%neuronKind = ' '
        ! # Neural compartments
        ! ## The instant of the last spike of the Motor unit
        ! ## at the Soma compartment.
        init_Interneuron%tSomaSpike = -1e6
        
        ! ## Vector with the instants of spikes at the soma.
        if (allocated(init_Interneuron%somaSpikeTrain)) then 
            deallocate(init_Interneuron%somaSpikeTrain)
        end if
        ! ## Integer corresponding to the Interneuron order in the pool.
        init_Interneuron%index = index
        ! ## Vector of Compartment of the Motor Unit.
        if (allocated(init_Interneuron%Compartments)) then
            deallocate(init_Interneuron%Compartments)
        end if

        ! ## Value of the membrane potential, in mV, that is considered a spike.
        paramTag = 'threshold'
        paramChar = init_Interneuron%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Interneuron%threshold_mV

        ! ## Anatomical position of the neuron, in mm.
        paramTag = 'position'
        paramChar = init_Interneuron%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Interneuron%position_mm
        
        allocate(init_Interneuron%Compartments(1))
        compKind = 'soma'
        init_Interneuron%Compartments(1)= Compartment(compKind, &
                            init_Interneuron%conf, init_Interneuron%pool,&
                            init_Interneuron%index, init_Interneuron%neuronKind)
        init_Interneuron%somaIndex = 1
        ! ## Number of compartments.
        init_Interneuron%compNumber = size(init_Interneuron%Compartments)
        ! ## Vector with membrane potential,in mV, of all compartments.
        allocate(init_Interneuron%v_mV(init_Interneuron%compNumber))
        

        allocate(gLeak(init_Interneuron%compNumber))
        allocate(capacitance_nF(init_Interneuron%compNumber))
        allocate(EqPot(init_Interneuron%compNumber))
        

        
        do i = 1, init_Interneuron%compNumber
            capacitance_nF(i) = init_Interneuron%Compartments(i)%capacitance_nF
            gLeak(i) = init_Interneuron%Compartments(i)%gLeak_muS
            EqPot(i) = init_Interneuron%Compartments(i)%EqPot_mV
            init_Interneuron%v_mV(i) = EqPot(i)
        end do

        ! ## Vector with  the inverse of the capacitance of all compartments.
        init_Interneuron%capacitanceInv = 1.0 / capacitance_nF

        ! ## Vector with current, in nA,  of each compartment coming from other elements of the model. For example
        ! ## from ionic channels and synapses.
        allocate(init_Interneuron%iIonic(init_Interneuron%compNumber))
        init_Interneuron%iIonic(:) = 0.0
        ! ## Vector with the current, in nA, injected in each compartment.
        allocate(init_Interneuron%iInjected(init_Interneuron%compNumber))
        init_Interneuron%iInjected(:) = 0.0
        
        allocate(GL(init_Interneuron%compNumber, init_Interneuron%compNumber))
        GL(:,:) = 0.0

        do i = 1, init_Interneuron%compNumber
            GL(i,i) = -gLeak(i)
        end do
        
        

        ! ## Matrix of the conductance of the motoneuron. Multiplied by the vector self.v_mV,
        ! ## results in the passive currents of each compartment.
        init_Interneuron%G = GL

        allocate(init_Interneuron%EqCurrent_nA(init_Interneuron%compNumber))
        init_Interneuron%EqCurrent_nA = matmul(-GL, EqPot)

        ! ## index of the soma compartment.
        
        ! ## Refractory period, in ms, of the motoneuron.
        paramTag = trim(init_Interneuron%pool) // 'SomaRefPer'
        paramChar = init_Interneuron%conf%parameterSet(paramTag, pool, index)
        read(paramChar, *)init_Interneuron%RefPer_ms
    
        ! ## Build synapses       
         
        init_Interneuron%SynapsesOut = CharacterMatrix()

        if (allocated(gLeak)) deallocate(gLeak)
        if (allocated(capacitance_nF)) deallocate(capacitance_nF)
        if (allocated(EqPot)) deallocate(EqPot)
        if (allocated(GL)) deallocate(GL)
    end function

    subroutine atualizeInterneuron(self, t, v_mV)
        ! '''
        ! Atualize the dynamical and nondynamical (delay) parts of the motor unit.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(Interneuron), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in) :: v_mV(self%compNumber)

        call self%atualizeCompartments(t, v_mV)
    end subroutine

    subroutine atualizeCompartments(self, t, v_mV)
        ! '''
        ! Atualize all neural compartments.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(Interneuron), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in) :: v_mV(self%compNumber)
        integer :: i
        
        self%v_mV(:) = v_mV

        if ((self%v_mV(self%somaIndex) > self%threshold_mV).and.(t-self%tSomaSpike > self%RefPer_ms)) then 
            call self%addSomaSpike(t)
        end if         
    end subroutine

    subroutine addSomaSpike(self, t)
        ! '''
        ! When the soma potential is above the threshold a spike is added to the soma.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(Interneuron), intent(inout) :: self
        real(wp), intent(in) :: t
        integer :: i, j

        self%tSomaSpike = t
        call AddToList(self%somaSpikeTrain, t)
        call self%transmitSpikes(t)

        do i = 1, self%Compartments(self%somaIndex)%numberChannels
            do j = 1, self%Compartments(self%somaIndex)%Channels(i)%lenStates
                call self%Compartments(self%somaIndex)%Channels(i)%condState(j)%changeState(t)
            end do
        end do
    end subroutine

    subroutine transmitSpikes(self, t)
            !     '''
            !     - Inputs:
            !         + **t**: current instant, in ms.
            !     '''
            class(Interneuron), intent(inout) :: self
            real(wp), intent(in) :: t
            integer :: i
            
            if (allocated(self%indicesOfSynapsesOnTarget)) then
                do i = 1, size(self%indicesOfSynapsesOnTarget)
                    call self%transmitSpikesThroughSynapses(i)%synapse%receiveSpike(t, self%indicesOfSynapsesOnTarget(i))
                end do
            end if
            
    end subroutine

    subroutine reset(self)
        ! '''

        ! '''
        class(Interneuron), intent(inout) :: self
        integer :: i

        self%tSomaSpike = -1e6
        
        do i = 1, self%compNumber
            self%v_mV(i) = self%Compartments(i)%EqPot_mV
            call self%Compartments(i)%reset()
        end do

        self%iIonic(:) = 0.0
        self%iInjected(:) = 0.0

        if (allocated(self%somaSpikeTrain)) deallocate(self%somaSpikeTrain)       
        
        !print *, self%pool, self%index
        
    end subroutine

end module InterneuronClass


        
    
    

        
