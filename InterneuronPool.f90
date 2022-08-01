! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2021  Renato Naville Watanabe
!                         Pablo Alejandro
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
module InterneuronPoolClass
    use ConfigurationClass
!L!    use mkl_spblas
    use DynamicalArrays
    use InterneuronClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    public :: InterneuronPool

    type InterneuronPool
        character(len = 2) :: poolKind
        character(len = 15) :: pool
        type(Configuration), pointer :: conf
        integer :: Nnumber, totalNumberOfCompartments
        type(Interneuron), dimension(:), allocatable :: unit
        real(wp), dimension(:,:), allocatable :: poolSomaSpikes
        real(wp), dimension(:), allocatable :: v_mV, iInjected, capacitanceInv
        real(wp), dimension(:), allocatable :: iIonic, EqCurrent_nA
        real(wp), dimension(:,:), allocatable :: G
        integer :: spIndexing, spRows, spCols, spNumberOfElements, spOperation
        integer, dimension(:), allocatable :: spRowStart, spRowEnd, spColIdx
        real(wp), dimension(:), allocatable :: spValues

        contains
            procedure :: dVdt
            procedure :: atualizeInterneuronPool
            procedure :: listSpikes
            procedure :: reset

    end type InterneuronPool

    interface InterneuronPool
        module procedure init_InterneuronPool
    end interface InterneuronPool

    contains

        type(InterneuronPool) function init_InterNeuronPool(conf, pool, group)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pool**: string with Interneuron pool to which the motor unit belongs.
            ! '''
            class(Configuration), intent(in), target :: conf
            character(len = 15), intent(in) :: pool
            character(len = 4), intent(in) :: group
            character(len = 80) :: paramTag, paramChar
            real(wp) :: paramReal
            integer :: i, j, stat

            ! ## Indicates that is Motor Unit pool.
            init_InterNeuronPool%poolKind = 'IN'

            ! ## Configuration object with the simulation parameters.
            init_InterNeuronPool%conf => conf
            ! ## String with Motor unit pool to which the motor unit belongs.
            init_InterNeuronPool%pool = trim(pool) // '_' // trim(group)
            ! ## Number of Neurons.
            
            paramTag = 'Number_' // trim(init_InterNeuronPool%pool)
            paramChar = init_InterNeuronPool%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)paramReal
            init_InterNeuronPool%Nnumber = nint(paramReal)

            ! ## List of Interneuron objects.
            allocate(init_InterNeuronPool%unit(init_InterNeuronPool%Nnumber))

            do i = 1, init_InterNeuronPool%Nnumber
                init_InterNeuronPool%unit(i) = Interneuron(init_InterNeuronPool%conf, init_InterNeuronPool%pool, i)
            end do

            ! ## Vector with the instants of spikes in the soma compartment, in ms.
            if (allocated(init_InterNeuronPool%poolSomaSpikes)) then
                deallocate(init_InterNeuronPool%poolSomaSpikes)    
            end if
            ! ##
            
            init_InterNeuronPool%totalNumberOfCompartments = 0

            do i = 1, init_InterNeuronPool%Nnumber
                init_InterNeuronPool%totalNumberOfCompartments = init_InterNeuronPool%totalNumberOfCompartments &
                    + init_InterNeuronPool%unit(i)%compNumber
            end do

            allocate(init_InterNeuronPool%v_mV(init_InterNeuronPool%totalNumberOfCompartments))
                
            allocate(init_InterNeuronPool%G(init_InterNeuronPool%totalNumberOfCompartments,&
                                                init_InterNeuronPool%totalNumberOfCompartments)) 
            allocate(init_InterNeuronPool%iInjected(init_InterNeuronPool%totalNumberOfCompartments))
            allocate(init_InterNeuronPool%capacitanceInv(init_InterNeuronPool%totalNumberOfCompartments))
            allocate(init_InterNeuronPool%iIonic(init_InterNeuronPool%totalNumberOfCompartments))
            allocate(init_InterNeuronPool%EqCurrent_nA(init_InterNeuronPool%totalNumberOfCompartments))

            

            init_InterNeuronPool%iInjected(:) = 0.0
            init_InterNeuronPool%iIonic(:) = 0.0
            ! # Retrieving data from Interneuron class
            init_InterNeuronPool%G(:,:) = 0.0
            do i = 1, init_InterNeuronPool%Nnumber
                init_InterNeuronPool%v_mV((i-1)*init_InterNeuronPool%unit(i)%compNumber+1:&
                                                    i*init_InterNeuronPool%unit(i)%compNumber) =&
                                                     init_InterNeuronPool%unit(i)%v_mV
                ! # With only one compartment, it is a diagonal matrix
                init_InterNeuronPool%G(i:i,i:i) = init_InterNeuronPool%unit(i)%G
                init_InterNeuronPool%capacitanceInv((i-1)*init_InterNeuronPool%unit(i)%compNumber+1:&
                                                    i*init_InterNeuronPool%unit(i)%compNumber) &
                                                    = init_InterNeuronPool%unit(i)%capacitanceInv
                init_InterNeuronPool%EqCurrent_nA((i-1)*init_InterNeuronPool%unit(i)%compNumber+1:&
                                                    i*init_InterNeuronPool%unit(i)%compNumber) &
                                                    = init_InterNeuronPool%unit(i)%EqCurrent_nA
            end do

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            init_InterNeuronPool%spIndexing = 1
            init_InterNeuronPool%spRows = init_InterNeuronPool%totalNumberOfCompartments
            init_InterNeuronPool%spCols = init_InterNeuronPool%totalNumberOfCompartments
            init_InterNeuronPool%spNumberOfElements = 0
            allocate(init_InterNeuronPool%spRowStart(init_InterNeuronPool%spRows))
            allocate(init_InterNeuronPool%spRowEnd(init_InterNeuronPool%spRows))
            init_InterNeuronPool%spRowStart(:) = 0
            do i = 1, init_InterNeuronPool%totalNumberOfCompartments
                do j = 1, init_InterNeuronPool%totalNumberOfCompartments
                    if (abs(init_InterNeuronPool%G(i,j))>1e-10) then
                        init_InterNeuronPool%spNumberOfElements = init_InterNeuronPool%spNumberOfElements + 1
                        if (init_InterNeuronPool%spRowStart(i) == 0) then
                            init_InterNeuronPool%spRowStart(i) = init_InterNeuronPool%spNumberOfElements
                        end if                    
                        call AddToList(init_InterNeuronPool%spValues, init_InterNeuronPool%G(i,j))
                        call integerAddToList(init_InterNeuronPool%spColIdx, j)
                    end if                

                end do
                init_InterNeuronPool%spRowEnd(i) = init_InterNeuronPool%spNumberOfElements+1
            end do
            print '(A)', 'Interneuron Pool of ' // trim(pool) // ' ' // trim(group) // ' built'


        end function

        subroutine atualizeInterneuronPool(self, t)
            ! '''
            ! Update all parts of the Motor Unit pool. It consists
            ! to update all motor units, the activation signal and
            ! the muscle force.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            ! '''
            class(InterneuronPool), intent(inout) :: self
            real(wp), intent(in) :: t
            real(wp), dimension(self%totalNumberOfCompartments) :: k1, k2, k3, k4
            real(wp), dimension(self%totalNumberOfCompartments) :: newStateTemp
            real(wp), dimension(:), allocatable :: vUnit
            real(wp) :: newtTemp
            integer :: i
            real(wp) :: vmax, vmin
            
            vmin = -30.0
            vmax = 120.0      

            k1 = self%dVdt(t, self%v_mV)        
            newStateTemp = self%v_mV + self%conf%timeStepByTwo_ms * k1
            newtTemp = t + self%conf%timeStepByTwo_ms
            k2 = self%dVdt(newtTemp, newStateTemp)
            newStateTemp = self%v_mV + self%conf%timeStepByTwo_ms * k2
            k3 = self%dVdt(newtTemp, newStateTemp)
            newStateTemp = self%v_mV + self%conf%timeStep_ms * k3
            newtTemp = t + self%conf%timeStep_ms
            k4 = self%dVdt(newtTemp, newStateTemp)
            
            self%v_mV = self%v_mV + self%conf%timeStepBySix_ms * (k1 + 2.0*k2 + 2.0*k3 + k4)

            self%v_mV = merge(self%v_mV, vmax, self%v_mV < vmax)
            self%v_mV = merge(self%v_mV, vmin, self%v_mV > vmin)

            allocate(vUnit(self%unit(1)%compNumber))
            do i = 1, self%Nnumber
                vUnit = self%v_mV((i-1)*self%unit(i)%compNumber+1:i*self%unit(i)%compNumber)
                call self%unit(i)%atualizeInterneuron(t, vUnit)
            end do
        end subroutine
    
        function dVdt(self, t, V)
            class(InterneuronPool), intent(inout) :: self
            real(wp), intent(in) :: t
            real(wp), intent(in) :: V(self%totalNumberOfCompartments)
            real(wp), dimension(self%totalNumberOfCompartments) :: dVdt
            integer :: i, j, stat
            real(wp), dimension(self%totalNumberOfCompartments) :: matInt
            
            do i = 1, self%Nnumber
                do j = 1, self%unit(i)%compNumber
                    self%iIonic((i-1)*self%unit(i)%compNumber+j) = &
                            self%unit(i)%Compartments(j)%computeCurrent(t, V((i-1)*self%unit(i)%compNumber+j))
                end do
            end do

			matInt = matmul(self%G, V)
            
            dVdt =  (self%iIonic + matInt + self%iInjected + self%EqCurrent_nA) * &
                    self%capacitanceInv
        end function

    subroutine listSpikes(self)
        ! '''
        ! List the spikes that occurred in the soma and in
        ! the terminal of the different motor units.
        ! '''
        class(InterneuronPool) , intent(inout) :: self
        integer :: i
        integer, dimension(self%Nnumber) :: numberOfNewSpikesSoma
        integer :: numberOfSpikesSoma
        integer :: initInd, endInd

        do i = 1, self%Nnumber
            if (allocated(self%unit(i)%somaSpikeTrain)) then 
                numberOfNewSpikesSoma(i) = size(self%unit(i)%somaSpikeTrain)
            else 
                numberOfNewSpikesSoma(i) = 0
            end if
        end do

        numberOfSpikesSoma = sum(numberOfNewSpikesSoma)
        allocate(self%poolSomaSpikes(numberOfSpikesSoma,2))        

        initInd = 1
        do i = 1, self%Nnumber
            if (allocated(self%unit(i)%somaSpikeTrain)) then
                endInd = initInd + numberOfNewSpikesSoma(i) - 1
                self%poolSomaSpikes(initInd:endInd,1) = self%unit(i)%somaSpikeTrain
                self%poolSomaSpikes(initInd:endInd,2) = i
                initInd = endInd+1                
            end if
        end do 
    end subroutine

    subroutine reset(self)
        ! '''
        ! '''
        class(InterneuronPool) , intent(inout) :: self
        integer :: i

        if (allocated(self%poolSomaSpikes)) deallocate(self%poolSomaSpikes)

        do i = 1, self%Nnumber
            call self%unit(i)%reset()
        end do

        do i = 1, self%Nnumber
                self%v_mV((i-1)*self%unit(i)%compNumber+1:&
                         i*self%unit(i)%compNumber) = &
                         self%unit(i)%v_mV
        end do

        print *, 'Interneuron pool reseted'
        ! read(*,*)
    end subroutine

end module InterneuronPoolClass
