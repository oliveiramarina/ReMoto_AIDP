! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2021  Renato Naville Watanabe
!                         Pablo Alejandro
!						  Marina Oliveira
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

module MotorUnitPoolClass
    ! '''
    ! Class that implements a motor unit pool. Encompasses a set of motor
    ! units that controls a single  muscle.
    ! '''
    use ConfigurationClass
    use DynamicalArrays
    use MotorUnitClass
    use MuscularActivationClass
    use MuscleNoHillClass
    use MuscleHillClass
    use MuscleSpindleClass
    use GolgiTendonOrganClass
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)

	real(wp),private :: pic, poc		!tst!
    public :: MotorUnitPool

    type MotorUnitPool
        real(wp) :: t
        character(len = 2) :: poolKind
        type(Configuration), pointer :: conf
        character(len = 6) :: pool
        integer :: MUnumber, totalNumberOfCompartments
        real(wp) :: muscleThickness_mm
        type(MotorUnit), dimension(:), allocatable:: unit
        real(wp), dimension(:), allocatable :: v_mV
        real(wp), dimension(:,:), allocatable :: G
        real(wp), dimension(:), allocatable :: iInjected, capacitanceInv, iIonic, EqCurrent_nA
        real(wp), dimension(:,:), allocatable :: poolSomaSpikes! ## Vector with the instants of spikes in the soma compartment, in ms.
        real(wp), dimension(:,:), allocatable :: poolLastCompSpikes   ! ## Vector with the instants of spikes in the last dynamical compartment, in ms.
        real(wp), dimension(:,:), allocatable :: poolTerminalSpikes ! ## Vector with the instants of spikes in the terminal, in ms.
        real(wp), dimension(:), allocatable :: emg
        type(MuscularActivation) :: Activation
        type(MuscleNoHill) :: NoHillMuscle
        type(MuscleHill) :: HillMuscle
        type(MuscleSpindle) :: spindle
        type(GolgiTendonOrgan) :: GTO
        character(len=80) :: hillModel
        integer(c_int) :: spIndexing, spRows, spCols, spNumberOfElements, spOperation
        integer(c_int), dimension(:), allocatable :: spRowStart, spRowEnd, spColIdx
        real(wp), dimension(:), allocatable :: spValues


        contains
            procedure :: dVdt
            procedure :: atualizeMotorUnitPool
            procedure :: getMotorUnitPoolInstantEMG
            procedure :: getMotorUnitPoolEMG
            procedure :: listSpikes
            procedure :: reset

    end type MotorUnitPool
    interface MotorUnitPool
        module procedure init_MotorUnitPool
    end interface MotorUnitPool

    contains

    type (MotorUnitPool) function init_MotorUnitPool(conf, pool)
        ! '''
        ! Constructor
        ! - Inputs:
        !     + **conf**: Configuration object with the simulation parameters.
        !     + **pool**: string with Motor unit pool to which the motor unit belongs.
        ! '''
        class(Configuration), intent(in), target :: conf
        character(*), intent(in):: pool
        character(len = 80) :: paramTag, paramChar
        integer :: MUnumber_S, MUnumber_FR, MUnumber_FF
        real(wp) :: paramReal
        integer(c_int) :: i, j
        character(len = 2) ::  neuronKind
        integer :: simDurationSteps, lastIndex, stat

        init_MotorUnitPool%t = 0.0

        ! ## Indicates that is Motor Unit pool.
        init_MotorUnitPool%poolKind = 'MU'

        ! ## Configuration object with the simulation parameters.
        init_MotorUnitPool%conf => conf

        ! ## String with Motor unit pool to which the motor unit belongs.
        init_MotorUnitPool%pool = pool

        paramTag = 'MUnumber_' // trim(pool) // '-S'
        paramChar = init_MotorUnitPool%conf%parameterSet(paramTag, pool, 0)
        read(paramChar, *)paramReal
        MUnumber_S = int(paramReal)

        paramTag = 'MUnumber_' // trim(pool) // '-FR'
        paramChar = init_MotorUnitPool%conf%parameterSet(paramTag, pool, 0)
        read(paramChar, *)paramReal
        MUnumber_FR = int(paramReal)

        paramTag = 'MUnumber_' // trim(pool) // '-FF'
        paramChar = init_MotorUnitPool%conf%parameterSet(paramTag, pool, 0)
        read(paramChar, *)paramReal
        MUnumber_FF = int(paramReal)

        ! ## Number of motor units.
        init_MotorUnitPool%MUnumber = MUnumber_S + MUnumber_FR + MUnumber_FF
        ! ## Muscle thickness, in mm.
        paramTag = 'thickness:' // trim(pool)
        paramChar = init_MotorUnitPool%conf%parameterSet(paramTag, pool, 0)
        read(paramChar, *)init_MotorUnitPool%muscleThickness_mm

        ! ## Vector of MotorUnit objects.
        allocate(init_MotorUnitPool%unit(init_MotorUnitPool%MUnumber))

		call cpu_time(pic)    !tst!
		print *, init_MotorUnitPool%MUnumber, ' <-- init_MotorUnitPool%MUnumber MotorUnit'		!tst! 900
		print *, 'init_MotorUnitPool%MUnumber=',init_MotorUnitPool%MUnumber,';MUnumber_S=',MUnumber_S,';MUnumber_FR=',MUnumber_FR,';MUnumber_FF=',MUnumber_FF,'[MotorUnitPool.f90]#tst#'		!tst! 900
		
		
        do i = 1, init_MotorUnitPool%MUnumber

            if (i <= MUnumber_S) then
                neuronKind = 'S'
                init_MotorUnitPool%unit(i) = MotorUnit(init_MotorUnitPool%conf, init_MotorUnitPool%pool,&
                                                       i, neuronKind, init_MotorUnitPool%muscleThickness_mm,&
                                                       init_MotorUnitPool%conf%skinThickness_mm)
            else if (i <= MUnumber_S + MUnumber_FR) then
                neuronKind = 'FR'
                init_MotorUnitPool%unit(i) = MotorUnit(init_MotorUnitPool%conf, init_MotorUnitPool%pool,&
                                                       i, neuronKind, init_MotorUnitPool%muscleThickness_mm,&
                                                       init_MotorUnitPool%conf%skinThickness_mm)
            else
                neuronKind = 'FF'
                init_MotorUnitPool%unit(i) = MotorUnit(init_MotorUnitPool%conf, init_MotorUnitPool%pool,&
                                                       i, neuronKind, init_MotorUnitPool%muscleThickness_mm,&
                                                       init_MotorUnitPool%conf%skinThickness_mm)
            end if
        end do

        call cpu_time(poc)    !tst!
		print '(F15.6, A)', poc - pic, ' seconds (... init_MotorUnitPool%MUnumber MotorUnit)'		!tst!  306.968750 seconds para 900

        ! # This is used to get values from MotorUnit.py and make computations
        ! # in MotorUnitPool.py

        init_MotorUnitPool%totalNumberOfCompartments = 0

        do i = 1, init_MotorUnitPool%MUnumber
            init_MotorUnitPool%totalNumberOfCompartments = init_MotorUnitPool%totalNumberOfCompartments &
                                                        + init_MotorUnitPool%unit(i)%compNumber
			!print *,'i=',i,';init_MotorUnitPool%unit(i)%compNumber=',init_MotorUnitPool%unit(i)%compNumber,'[MotorUnitPool.f90]#tst#'		!tst!
        end do


        allocate(init_MotorUnitPool%v_mV(init_MotorUnitPool%totalNumberOfCompartments))
        init_MotorUnitPool%v_mV(:) = 0.0

        allocate(init_MotorUnitPool%G(init_MotorUnitPool%totalNumberOfCompartments, init_MotorUnitPool%totalNumberOfCompartments))
        allocate(init_MotorUnitPool%iInjected(init_MotorUnitPool%totalNumberOfCompartments))
        init_MotorUnitPool%iInjected(:) = 0.0
        allocate(init_MotorUnitPool%capacitanceInv(init_MotorUnitPool%totalNumberOfCompartments))
        allocate(init_MotorUnitPool%iIonic(init_MotorUnitPool%totalNumberOfCompartments))
        init_MotorUnitPool%iIonic(:) = 0.0
        allocate(init_MotorUnitPool%EqCurrent_nA(init_MotorUnitPool%totalNumberOfCompartments))

        init_MotorUnitPool%G(:,:) = 0.0

        ! # Retrieving data from Motorneuron class
        ! # Vectors or matrices from Motorneuron compartments are copied,
        ! # populating larger vectors or matrices that will be used for computations

		!call cpu_time(pic)    !tst!
		!print *, init_MotorUnitPool%MUnumber, ' <-- init_MotorUnitPool%MUnumber v_mV capacitanceInv' 		!tst!

        do i = 1, init_MotorUnitPool%MUnumber
            init_MotorUnitPool%v_mV((i-1)*init_MotorUnitPool%unit(i)%compNumber+1:&
                                    (i)*init_MotorUnitPool%unit(i)%compNumber)&
                                     = init_MotorUnitPool%unit(i)%v_mV
            ! # Consists of smaller matrices on its diagonal
            init_MotorUnitPool%G((i-1)*init_MotorUnitPool%unit(i)%compNumber+1:&
                                 (i)*init_MotorUnitPool%unit(i)%compNumber, &
                                 (i-1)*init_MotorUnitPool%unit(i)%compNumber+1:&
                                 i*init_MotorUnitPool%unit(i)%compNumber) &
                                 = init_MotorUnitPool%unit(i)%G
            init_MotorUnitPool%capacitanceInv((i-1)*init_MotorUnitPool%unit(i)%compNumber+1: &
                                              i*init_MotorUnitPool%unit(i)%compNumber) &
                                              = init_MotorUnitPool%unit(i)%capacitanceInv
            init_MotorUnitPool%EqCurrent_nA((i-1)*init_MotorUnitPool%unit(i)%compNumber+1: &
                                            i*init_MotorUnitPool%unit(i)%compNumber) &
                                            = init_MotorUnitPool%unit(i)%EqCurrent_nA
        end do

        !call cpu_time(poc)    !tst!
		!print '(F15.6, A)', poc - pic, ' seconds (... init_MotorUnitPool%MUnumber v_mV capacitanceInv)'		!tst!

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        init_MotorUnitPool%spIndexing = 1
        init_MotorUnitPool%spRows = init_MotorUnitPool%totalNumberOfCompartments
        init_MotorUnitPool%spCols = init_MotorUnitPool%totalNumberOfCompartments
        init_MotorUnitPool%spNumberOfElements = 0
        allocate(init_MotorUnitPool%spRowStart(init_MotorUnitPool%spRows))
        allocate(init_MotorUnitPool%spRowEnd(init_MotorUnitPool%spRows))
        init_MotorUnitPool%spRowStart(:) = 0

		!call cpu_time(pic)    !tst!
		!print *, init_MotorUnitPool%totalNumberOfCompartments, ' <-- init_MotorUnitPool%totalNumberOfCompartments' 		!tst!

        do i = 1, init_MotorUnitPool%totalNumberOfCompartments
            do j = 1, init_MotorUnitPool%totalNumberOfCompartments
                if (abs(init_MotorUnitPool%G(i,j))>1e-10) then
                    init_MotorUnitPool%spNumberOfElements = init_MotorUnitPool%spNumberOfElements + 1
                    if (init_MotorUnitPool%spRowStart(i) == 0) then
                        init_MotorUnitPool%spRowStart(i) = init_MotorUnitPool%spNumberOfElements
                    end if
                    call AddToList(init_MotorUnitPool%spValues, init_MotorUnitPool%G(i,j))			!L! obs: spValues !!!
                    call integerAddToList(init_MotorUnitPool%spColIdx, j)
                end if
            end do
            init_MotorUnitPool%spRowEnd(i) = init_MotorUnitPool%spNumberOfElements+1
        end do

        !call cpu_time(poc)    !tst!
		!print '(F15.6, A)', poc - pic, ' seconds (... totalNumberOfCompartments)'		!tst! 0.046875 seconds

!L!     stat = mkl_sparse_d_create_csr(init_MotorUnitPool%GSp, &
!L!                                    init_MotorUnitPool%spIndexing, &
!L!                                    init_MotorUnitPool%spRows, &
!L!                                    init_MotorUnitPool%spCols, &
!L!                                    init_MotorUnitPool%spRowStart, &
!L!                                    init_MotorUnitPool%spRowEnd, &
!L!                                    init_MotorUnitPool%spColIdx, &
!L!                                    init_MotorUnitPool%spValues)


        ! Values for the matrix-vector operation matInt = GV
!L!     init_MotorUnitPool%spDescr%type = SPARSE_MATRIX_TYPE_GENERAL
!L!     init_MotorUnitPool%spAlpha = 1.0
!L!     init_MotorUnitPool%spOperation = SPARSE_OPERATION_NON_TRANSPOSE
!L!     init_MotorUnitPool%spBeta = 0.0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! #activation signal
        init_MotorUnitPool%Activation = MuscularActivation(init_MotorUnitPool%conf,init_MotorUnitPool%pool,&
                                             init_MotorUnitPool%MUnumber, init_MotorUnitPool%unit)

        ! #Force
        ! ## String indicating whther a Hill  model is used or not. For now, it can be *No*.
        paramTag = 'hillModel'
        init_MotorUnitPool%hillModel = init_MotorUnitPool%conf%parameterSet(paramTag, pool, 0)
        if (trim(init_MotorUnitPool%hillModel).eq.'No') then
            init_MotorUnitPool%NoHillMuscle = MuscleNoHill(init_MotorUnitPool%conf, &
                                                           init_MotorUnitPool%pool, &
                                                           init_MotorUnitPool%MUnumber, &
                                                           MUnumber_S, init_MotorUnitPool%unit)
        else
            init_MotorUnitPool%HillMuscle = MuscleHill(init_MotorUnitPool%conf, &
                                                   init_MotorUnitPool%pool, &
                                                   init_MotorUnitPool%MUnumber, &
                                                   MUnumber_S, init_MotorUnitPool%unit)
        end if
        ! # EMG
        ! ## EMG along time, in mV.
        simDurationSteps = nint(init_MotorUnitPool%conf%simDuration_ms/init_MotorUnitPool%conf%timeStep_ms)

        allocate(init_MotorUnitPool%emg(simDurationSteps))
        init_MotorUnitPool%emg(:) = 0.0

        ! Spindle
        init_MotorUnitPool%spindle = MuscleSpindle(init_MotorUnitPool%conf, init_MotorUnitPool%pool)

        ! GTO
        init_MotorUnitPool%GTO = GolgiTendonOrgan(init_MotorUnitPool%conf, init_MotorUnitPool%pool)

        ! ##
        print '(A)', 'Motor Unit Pool ' // trim(pool) // ' built'



    end function

    function dVdt(self, t, V)
        class(MotorUnitPool), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in) :: V(:)
        real(wp), dimension(self%totalNumberOfCompartments) :: dVdt
        integer :: i, j, stat
        real(wp), dimension(:), allocatable :: matInt

	    real(wp) :: pic, poc		!tst!

        allocate(matInt(self%totalNumberOfCompartments))

        do i = 1, self%MUnumber
            do j = 1, self%unit(i)%compNumber
                self%iIonic((i-1)*self%unit(i)%compNumber+j) = &
                    self%unit(i)%Compartments(j)%computeCurrent(t, &
                    V((i-1)*self%unit(i)%compNumber+j))
            end do
        end do

        matInt = matmul(self%G, V)		!<--- !L!
        dVdt = (self%iIonic + matInt + self%iInjected &
                + self%EqCurrent_nA) * self%capacitanceInv
        if (allocated(matInt)) deallocate(matInt)

    end function

    subroutine atualizeMotorUnitPool(self, t, dynGamma, statGamma)
        ! '''
        ! Update all parts of the Motor Unit pool. It consists
        ! to update all motor units, the activation signal and
        ! the muscle force.
        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(MotorUnitPool) , intent(inout), target:: self
        real(wp), intent(in) :: t
        real(wp), intent(in) :: dynGamma, statGamma
        real(wp), dimension(:), allocatable :: k1, k2, k3, k4, newStateTemp
        real(wp), dimension(:), allocatable :: vUnit
        real(wp) :: newtTemp
        integer :: i
        real(wp) :: vmax, vmin
        real(wp) :: length, velocity, acceleration, tendonForce

        allocate(k1(self%totalNumberOfCompartments))
        allocate(k2(self%totalNumberOfCompartments))
        allocate(k3(self%totalNumberOfCompartments))
        allocate(k4(self%totalNumberOfCompartments))
        allocate(newStateTemp(self%totalNumberOfCompartments))

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

        self%v_mV = self%v_mV + self%conf%timeStepBySix_ms * (k1+ 2.0*k2 + 2.0*k3 + k4)

        self%v_mV = merge(self%v_mV, vmax, self%v_mV < vmax)
        self%v_mV = merge(self%v_mV, vmin, self%v_mV > vmin)

        allocate(vUnit(self%unit(1)%compNumber))
        do i = 1, self%MUnumber
            vUnit = self%v_mV((i-1)*self%unit(i)%compNumber+1:i*self%unit(i)%compNumber)
            call self%unit(i)%atualizeMotorUnit(t, vUnit)
        end do

        self%Activation%unit => self%unit(:)

        call self%Activation%atualizeActivationSignal(t)
        if (trim(self%hillModel).eq.'No') then
            call self%NoHillMuscle%atualizeForce(self%Activation%activation_Sat)
            length = self%NoHillMuscle%lengthNorm
            velocity = self%NoHillMuscle%velocityNorm
            acceleration = self%NoHillMuscle%accelerationNorm
            tendonForce = self%NoHillMuscle%force(nint(t/self%conf%timeStep_ms)+1)
        else
            call self%HillMuscle%atualizeForce(self%Activation%activation_Sat)
            length = self%HillMuscle%lengthNorm
            velocity = self%HillMuscle%velocityNorm
            acceleration = self%HillMuscle%accelerationNorm
            tendonForce = self%HillMuscle%tendonForce_N(nint(t/self%conf%timeStep_ms)+1)
        end if

        call self%spindle%atualizeMuscleSpindle(t, length, velocity, acceleration, dynGamma, statGamma)

        call self%GTO%atualizeGolgiTendonOrgan(t, tendonForce) !Feedback

        deallocate(k1)
        deallocate(k2)
        deallocate(k3)
        deallocate(k4)
        deallocate(vUnit)
    end subroutine

    subroutine listSpikes(self)
        ! '''
        ! List the spikes that occurred in the soma and in
        ! the terminal of the different motor units.
        ! '''
        class(MotorUnitPool) , intent(inout) :: self
        integer :: i
        integer, dimension(self%MUnumber) :: numberOfNewSpikesSoma, numberOfNewSpikesLastComp, numberOfNewSpikesTerminal
        integer :: numberOfSpikesSoma, numberOfSpikesLastComp, numberOfSpikesTerminal
        integer :: initInd, endInd

        do i = 1, self%MUnumber
            if (allocated(self%unit(i)%lastCompSpikeTrain)) then
                numberOfNewSpikesLastComp(i) = size(self%unit(i)%lastCompSpikeTrain)
            else
                numberOfNewSpikesLastComp(i) = 0
            end if
            if (allocated(self%unit(i)%somaSpikeTrain)) then
                numberOfNewSpikesSoma(i) = size(self%unit(i)%somaSpikeTrain)
            else
                numberOfNewSpikesSoma(i) = 0
            end if
            if (allocated(self%unit(i)%terminalSpikeTrain)) then
                numberOfNewSpikesTerminal(i) = size(self%unit(i)%terminalSpikeTrain)
            else
                numberOfNewSpikesTerminal(i) = 0
            end if
        end do

        numberOfSpikesSoma = sum(numberOfNewSpikesSoma)
        numberOfSpikesTerminal = sum(numberOfNewSpikesTerminal)
        numberOfSpikesLastComp = sum(numberOfNewSpikesLastComp)
		
		print *, 'numberOfSpikesTerminal = ', numberOfSpikesTerminal

        allocate(self%poolSomaSpikes(numberOfSpikesSoma,2))
        allocate(self%poolLastCompSpikes(numberOfSpikesLastComp,2))
        allocate(self%poolTerminalSpikes(numberOfSpikesTerminal,2))

        initInd = 1
        do i = 1, self%MUnumber
            if (allocated(self%unit(i)%terminalSpikeTrain)) then
                endInd = initInd + numberOfNewSpikesTerminal(i) - 1
                self%poolTerminalSpikes(initInd:endInd,1) = self%unit(i)%terminalSpikeTrain
                self%poolTerminalSpikes(initInd:endInd,2) = i
                initInd = endInd+1
            end if
        end do

        initInd = 1
        do i = 1, self%MUnumber
            if (allocated(self%unit(i)%somaSpikeTrain)) then
                endInd = initInd + numberOfNewSpikesSoma(i) - 1
                self%poolSomaSpikes(initInd:endInd,1) = self%unit(i)%somaSpikeTrain
                self%poolSomaSpikes(initInd:endInd,2) = i
                initInd = endInd + 1
            end if
        end do

        initInd = 1
        do i = 1, self%MUnumber
            if (allocated(self%unit(i)%lastCompSpikeTrain)) then
                endInd = initInd + numberOfNewSpikesLastComp(i) - 1
                self%poolLastCompSpikes(initInd:endInd,1) = self%unit(i)%lastCompSpikeTrain
                self%poolLastCompSpikes(initInd:endInd,2) = i
                initInd = endInd + 1
            end if
        end do
    end subroutine

    real(wp) function getMotorUnitPoolInstantEMG(self, t) result(emg)
        ! '''
        ! '''
        class(MotorUnitPool), intent(inout) :: self
        real(wp), intent(in) :: t
        integer :: j
        real(wp), dimension(self%MUnumber) :: newEMG


        do j = 1, self%MUnumber
            newEMG(j)  = self%unit(j)%getEMG(t)
        end do

        emg = sum(newEMG)
    end function

    subroutine getMotorUnitPoolEMG(self)
        ! '''
        ! '''
        class(MotorUnitPool), intent(inout) :: self
        integer :: i
        real(wp) :: instant
        integer :: simDuration

        simDuration = size(self%emg)

        do i = 1, simDuration
            instant =  (i-1) * self%conf%timeStep_ms
            self%emg(i) = self%getMotorUnitPoolInstantEMG(instant)
        end do
    end subroutine

    subroutine reset(self)
        ! '''
        ! '''
        class(MotorUnitPool), intent(inout) :: self
        integer :: i

        if (allocated(self%poolSomaSpikes)) deallocate(self%poolSomaSpikes)
        if (allocated(self%poolLastCompSpikes)) deallocate(self%poolLastCompSpikes)
        if (allocated(self%poolTerminalSpikes)) deallocate(self%poolTerminalSpikes)
        self%emg(:) = 0.0

        do i = 1, self%MUnumber
            call self%unit(i)%reset()
        end do

        do i = 1, self%MUnumber
            self%v_mV((i-1)*self%unit(i)%compNumber+1:&
                       i*self%unit(i)%compNumber) = self%unit(i)%v_mV
        end do

        call self%Activation%reset()
        if (trim(self%hillModel)=='No') then
            call self%NoHillMuscle%reset()
        else
            call self%HillMuscle%reset()
        end if
        call self%spindle%reset()
        print *, 'Motorneuron Pool reseted'
    end subroutine


end module MotorUnitPoolClass