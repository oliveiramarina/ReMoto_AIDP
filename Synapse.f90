! '''
!     Neuromuscular simulator in Python.
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

module SynapseClass
    ! '''
    ! Implements the synapse model from Destexhe (1994)
    ! using the computational method from Lytton (1996).
    ! '''
    use ConfigurationClass
    use CharacterArrayClass
    use QueueClass
    use DynamicalArrays
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: pi = 4 * atan(1.0_wp)
    public :: Synapse

    type Synapse
        character(len = 9) :: pool
        character(len = 10) :: synapseKind
        character(len = 2) :: neuronKind
        integer ::index
        type(Configuration), pointer :: conf
        real(wp) :: timeStep_ms
        real(wp), dimension(:), allocatable :: gmax_muS, delay_ms, variation, timeConstant_ms
        type(CharacterArray) :: dynamics
        real(wp) :: gMaxTot_muS, rInf, tauOn, tauOff
        integer :: numberOfIncomingSynapses
        real(wp) :: expFinish
        real(wp) :: Non, Ron, Roff, t0
        real(wp), dimension(:), allocatable :: tBeginOfPulse
        logical, dimension(:), allocatable :: conductanceState
        real(wp), dimension(:), allocatable :: tLastPulse, tEndOfPulse
        ! List with the fractions of postsynaptic receptors that are bound to neurotransmitters of the individual synapses.
        real(wp), dimension(:), allocatable :: ri
        ! ## List with the instants of spike arriving at each conductance, in ms.
        real(wp), dimension(:), allocatable :: ti, dynamicGmax
        ! List of individual conductance constribution to the global synaptic conductance
            ! ## (\f$S_{indCont} = \frac{g_{i_{max}}{G_{max}}\f$).
        real(wp), dimension(:), allocatable :: synContrib
        type(Queue) :: inQueue, outQueue
        real(wp) :: EqPot_mV, alpha_ms1, beta_ms1, Tmax_mM, tPeak_ms
        real(wp) :: lastUpdateInstant

        contains
            procedure :: computeConductance
            procedure :: computeCurrent
            procedure :: startConductance
            procedure :: stopConductance
            procedure :: addConductance
            procedure :: reset
            procedure :: receiveSpike

    end type Synapse

    interface Synapse
        module procedure init_Synapse
    end interface Synapse

    contains

        type(Synapse) function init_Synapse(conf, pool, index, compartment, synapseKind, neuronKind)
            ! '''
            ! Constructor

            ! - Input:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pool**: string with identification of the pool to which
            !     the synapse belongs.

            !     + **index**: integer identification of the unit in the pool.

            !     + **compartment**: integer identification of the compartment of the unit
            !     where the synapse is.

            !     + **kind**: string with the type of synapse. It can be *excitatory* or *inhibitory*.

            !     + **neuronKind**:
            ! '''
            class(Configuration), intent(in), target :: conf
            character(*), intent(in) :: pool
            integer, intent(in) :: index
            character(*), intent(in) :: compartment
            character(*), intent(in) :: synapseKind
            character(*), intent(in) ::  neuronKind
            character(len = 80) :: paramTag, paramChar


            init_Synapse%pool = pool
            init_Synapse%synapseKind = synapseKind
            init_Synapse%neuronKind = neuronKind
            init_Synapse%index = index
            init_Synapse%conf => conf

            init_Synapse%timeStep_ms = init_Synapse%conf%timeStep_ms

            paramTag = 'EqPotSyn:' // trim(pool) // '-' // trim(init_Synapse%neuronKind) // '|' // init_Synapse%synapseKind
            paramChar = init_Synapse%conf%parameterSet(paramTag, pool, index)
            read(paramChar,*)init_Synapse%EqPot_mV

            paramTag = 'alphaSyn:' // trim(pool) // '-'  // trim(init_Synapse%neuronKind) // '|' // trim(init_Synapse%synapseKind)
            paramChar = init_Synapse%conf%parameterSet(paramTag, pool, index)
            read(paramChar,*)init_Synapse%alpha_ms1

            paramTag = 'betaSyn:' // trim(pool) // '-'  // trim(init_Synapse%neuronKind) // '|' // trim(init_Synapse%synapseKind)
            paramChar = init_Synapse%conf%parameterSet(paramTag, pool, index)
            read(paramChar,*)init_Synapse%beta_ms1

            paramTag = 'TmaxSyn:' // trim(pool) // '-'  // trim(init_Synapse%neuronKind) // '|' // trim(init_Synapse%synapseKind)
            paramChar = init_Synapse%conf%parameterSet(paramTag, pool, index)
            read(paramChar,*)init_Synapse%Tmax_mM
            ! ## Pulse duration, in ms.
            paramTag = 'tPeakSyn:' // trim(pool) // '-' // trim(init_Synapse%neuronKind) // '|' // trim(init_Synapse%synapseKind)
            paramChar = init_Synapse%conf%parameterSet(paramTag, pool, index)
            read(paramChar,*)init_Synapse%tPeak_ms


            init_Synapse%dynamics = CharacterArray()

            ! ## The sum of individual conductances of all synapses in
            ! ## the compartment, in \f$\mu\f$S (\f$G_{max} = \limits\sum_{i=1}^Ng_i\f$).
            init_Synapse%gMaxTot_muS = 0.0
            init_Synapse%numberOfIncomingSynapses = 0

            ! ## The fraction of postsynaptic receptors
            ! ## that would be bound to neurotransmitters
            ! ## after an infinite amount of time with
            ! ## neurotransmitter being released.
            init_Synapse%rInf = (init_Synapse%alpha_ms1 * init_Synapse%Tmax_mM) /&
                                (init_Synapse%alpha_ms1 * init_Synapse%Tmax_mM + init_Synapse%beta_ms1)
            ! ## Time constant during a pulse, in ms.
            ! ## \f$\tau_{on}=\frac{1}{\alpha.T_{max} +\beta}\f$
            init_Synapse%tauOn = 1.0 / (init_Synapse%alpha_ms1 * init_Synapse%Tmax_mM + init_Synapse%beta_ms1)
            ! ## Time constant after a pulse, in ms.
            ! ## \f$\tau_{off}=\frac{1}{\beta}\f$
            init_Synapse%tauOff = 1.0 / init_Synapse%beta_ms1
            ! ## Is the value of the exponential at the
            ! ## end of the pulse. It is computed as
            ! ## \f$\exp(T_{dur}/\tau_{on})\f$.
            init_Synapse%expFinish = exp(-init_Synapse%tPeak_ms / init_Synapse%tauOn)


            ! ## Sum of the fractions of the individual conductances that are
            ! ## receiving neurotransmitter (during pulse) relative to
            ! ## the \f$G_{max}\f$. (\f$N_{on}=\limits\sum_{i=1}g_{i_{on}}/G_{max}).
            init_Synapse%Non = 0.0
            ! ## Sum of the fraction of postsynaptic receptors
            ! ## that are bound to neurotransmitters of all the individual synapses
            ! ## that have neurotransmitters being released (during the pulse).
            init_Synapse%Ron = 0.0
            ! ## Sum of the fraction of postsynaptic receptors
            ! ## that are bound to neurotransmitters of all the individual synapses
            ! ## that do not have neurotransmitters being released (before and after
            ! ## the pulse).
            init_Synapse%Roff = 0.0
            ! ## Instant that the last spike arrived to the compartment.
            init_Synapse%t0 = 0.0

            init_Synapse%lastUpdateInstant = 0.0

            init_Synapse%inQueue = Queue()
            init_Synapse%outQueue = Queue()

            if (allocated(init_synapse%gmax_muS)) deallocate(init_synapse%gmax_muS)
            if (allocated(init_synapse%delay_ms)) deallocate(init_synapse%delay_ms)
            if (allocated(init_synapse%variation)) deallocate(init_synapse%variation)
            if (allocated(init_synapse%timeConstant_ms)) deallocate(init_synapse%timeConstant_ms)

            if (allocated(init_synapse%tBeginOfPulse)) deallocate(init_synapse%tBeginOfPulse)
            if (allocated(init_synapse%tEndOfPulse)) deallocate(init_synapse%tEndOfPulse)
            if (allocated(init_synapse%tLastPulse)) deallocate(init_synapse%tLastPulse)
            if (allocated(init_synapse%conductanceState)) deallocate(init_synapse%conductanceState)
            if (allocated(init_synapse%ri)) deallocate(init_synapse%ri)
            if (allocated(init_synapse%ti)) deallocate(init_synapse%ti)
            if (allocated(init_synapse%dynamicGmax)) deallocate(init_synapse%dynamicGmax)
            if (allocated(init_synapse%synContrib)) deallocate(init_synapse%synContrib)

            !self.startDynamicFunction = []
            !self.stopDynamicFunction = []

            !init_Synapse%computeCurrent => init_Synapse%computeCurrentInic
        end function

        real(wp) function computeCurrent(self, t, V_mV) result(current)
            ! '''
            ! Computes the current on the compartment due to the synapse.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **V_mV**: membrane potential of the compartment that the
            !     synapse belongs, in mV.

            ! - Output:
            !     + The current on the compartment due to the synapse.
            ! '''
            class(Synapse), intent(inout) :: self
            real(wp), intent(in) :: t, V_mV

            ! if (.not.allocated(self%tEndOfPulse)) then
            !     allocate(self%tBeginOfPulse(size(self%gmax_muS)))
            !     self%tBeginOfPulse(:) = -1e4

            !     allocate(self%tEndOfPulse(size(self%gmax_muS)))
            !     self%tEndOfPulse(:) = -1e4

            !     allocate(self%tLastPulse(size(self%gmax_muS)))
            !     self%tLastPulse(:) = -1e4

            !     allocate(self%conductanceState(size(self%gmax_muS)))
            !     self%conductanceState(:) = .false.

            !     allocate(self%ri(size(self%gmax_muS)))
            !     self%ri(:) = 0.0

            !     allocate(self%ti(size(self%gmax_muS)))
            !     self%ti(:) = 0.0

            !     allocate(self%dynamicGmax(size(self%gmax_muS)))
            !     self%dynamicGmax(:) = self%gmax_muS

            !     allocate(self%synContrib(size(self%gmax_muS)))
            !     self%synContrib = self%gmax_muS / self%gMaxTot_muS
            !     !self%computeCurrent => self%computeCurrent2

            ! end if

            current = self%computeConductance(t) * (self%EqPot_mV - V_mV)
        end function

        ! real(wp) function computeCurrent2(self, t, V_mV) result(current)
        !     ! '''
        !     ! The same function of computeCurrent. It overrides this function for
        !     ! computational efficiency.

        !     ! - Inputs:
        !     !     + **t**: current instant, in ms.

        !     !     + **V_mV**: membrane potential of the compartment that the
        !     !     synapse belongs, in mV.
        !     ! '''
        !     class(Synapse), intent(inout) :: self
        !     real(wp), intent(in) :: t, V_mV

        !     current = self%computeConductance(t) * (self%EqPot_mV - V_mV)
        ! end function

        real(wp) function computeConductance(self, t) result(conductance)
            ! '''

            ! - Inputs:
            !     + **t**: current instant, in ms.
            ! '''
            class(Synapse), intent(inout) :: self
            real(wp), intent(in) :: t
            integer, dimension(:), allocatable :: idxBeginPulse, idxEndPulse
            logical :: continueFlag
            integer :: newPulse

            self%Ron = self%Ron * exp((self%lastUpdateInstant-t) / self%tauOn) + &
                       self%Non * self%rInf * (1.0 - exp((self%lastUpdateInstant-t) / self%tauOn))
            self%Roff = self%Roff*(exp((self%lastUpdateInstant-t) / self%tauOff))
            self%lastUpdateInstant = t

            if (allocated(idxBeginPulse)) deallocate(idxBeginPulse)
            continueFlag = .true.
            do while (self%inQueue%endQueue>0.and.continueFlag)
                if (t >= self%tBeginOfPulse(self%inQueue%item(1))) then
                    newPulse = self%inQueue%popleft()
                    call integerAddToList(idxBeginPulse, newPulse)
                else
                    continueFlag = .false.
                end if
            end do

            if (allocated(idxEndPulse)) deallocate(idxEndPulse)
            continueFlag = .true.
            do while (self%outQueue%endQueue>0.and.continueFlag)
                if  (t >= self%tEndOfPulse(self%outQueue%item(1))) then
                    newPulse = self%outQueue%popleft()
                    call integerAddToList(idxEndPulse, newPulse)
                else
                    continueFlag = .false.
                end if
            end do

            if (allocated(idxBeginPulse)) call self%startConductance(t, idxBeginPulse)

            if (allocated(idxEndPulse)) call self%stopConductance(t, idxEndPulse)

            conductance = self%gMaxTot_muS * (self%Ron + self%Roff)

            if (allocated(idxBeginPulse)) deallocate(idxBeginPulse)
            if (allocated(idxEndPulse)) deallocate(idxEndPulse)
        end function

        subroutine startConductance(self, t, idxBeginPulse)
            ! '''
            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **idxBeginPulse**: integer with the index of the conductance
            !         that the pulse begin at time **t**.
            ! '''
            class(Synapse), intent(inout) :: self
            real(wp), intent(in) :: t
            integer, intent(in) :: idxBeginPulse(:)
            integer, dimension(:), allocatable :: idxTurningOnCond
            integer, dimension(:), allocatable :: logCount
            real(wp) :: synGain
            integer :: i


            self%dynamicGmax(idxBeginPulse) =  self%gmax_muS(idxBeginPulse) +&
                                               exp((self%tLastPulse(idxBeginPulse) - t) / &
                                               self%timeConstant_ms(idxBeginPulse)) *&
                                               (self%dynamicGmax(idxBeginPulse) *&
                                                self%variation(idxBeginPulse) -&
                                                self%gmax_muS(idxBeginPulse))

            self%synContrib(idxBeginPulse) = self%dynamicGmax(idxBeginPulse) / self%gMaxTot_muS

            do i = 1, size(idxBeginPulse)
                if (self%conductanceState(idxBeginPulse(i))) call self%outQueue%remove(idxBeginPulse(i))
            end do

            call self%outQueue%extend(idxBeginPulse)

            if (allocated(idxTurningOnCond)) deallocate(idxTurningOnCond)

            allocate(logCount(size(idxBeginPulse)))
            logCount(:) = 0
            do i = 1, size(idxBeginPulse)
                if (.not.self%conductanceState(idxBeginPulse(i))) then
                    call integerAddToList(idxTurningOnCond, idxBeginPulse(i))
                end if
            end do


            if (allocated(idxTurningOnCond)) then
                self%conductanceState(idxTurningOnCond) = .true.
                self%ri(idxTurningOnCond) = self%ri(idxTurningOnCond)*&
                                exp((self%ti(idxTurningOnCond)+self%tPeak_ms - t) / self%tauOff)
                self%Non = self%Non + sum(self%synContrib(idxTurningOnCond))
                self%ti(idxTurningOnCond) = t
                synGain = dot_product(self%ri(idxTurningOnCond), self%synContrib(idxTurningOnCond))
                self%Ron = self%Ron + synGain
                self%Roff = self%Roff - synGain
            end if

            self%tEndOfPulse(idxBeginPulse) = t + self%tPeak_ms
            self%tLastPulse(idxBeginPulse) = self%tBeginOfPulse(idxBeginPulse)
            self%tBeginOfPulse(idxBeginPulse) = -1e4

            if (allocated(idxTurningOnCond)) deallocate(idxTurningOnCond)
            if (allocated(logCount)) deallocate(logCount)
        end subroutine

        subroutine stopConductance(self, t, idxEndPulse)
            ! '''
            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **idxEndPulse**: integer with the index of the conductance
            !         that the pulse end at time **t**.
            ! '''
            class(Synapse), intent(inout) :: self
            real(wp), intent(in) :: t
            integer, intent(in) :: idxEndPulse(:)
            real(wp) :: synLost

            self%ri(idxEndPulse) = self%rInf + (self%ri(idxEndPulse) - self%rInf) * self%expFinish
            synLost = dot_product(self%ri(idxEndPulse),self%synContrib(idxEndPulse))
            self%Ron = self%Ron - synLost
            self%Roff = self%Roff + synLost
            self%Non = self%Non - sum(self%synContrib(idxEndPulse))
            self%tEndOfPulse(idxEndPulse) = -1e4
            self%conductanceState(idxEndPulse) = .false.
        end subroutine

        subroutine receiveSpike(self, t, synapseNumber)
            ! '''
            ! - Inputs:
            !     + **t**:
            !     + **synapseNumber**:
            ! '''
            class(Synapse), intent(inout) :: self
            real(wp), intent(in) :: t
            integer, intent(in) :: synapseNumber

			!print *, '(receiveSpike) synapseNumber=', synapseNumber, allocated(self%tBeginOfPulse), size(self%tBeginOfPulse), '<== debug'		!debug!

            self%tBeginOfPulse(synapseNumber) = t + self%delay_ms(synapseNumber)
			call self%inQueue%append(synapseNumber)


        end subroutine

        subroutine addConductance(self, gmax, delay, dynamics, variation, timeConstant)
            ! '''
            ! Adds a synaptic conductance to the compartment. As the computation
            ! is performed once for each compartment at each time step, the data of
            ! each individual synapse is integrate in a big synapse.

            ! - Inputs:
            !     + **gmax**: the maximum conductance of the individual
            !     synase, in \f$\mu\f$S.

            !     + **delay**: transmission delay between the transmitter of the
            !     spike and the receiver compartment, in ms.

            !     + **dynamics**: type of the synapse dynamics. For now it
            !     can be *None*.

            ! '''

            class(Synapse), intent(inout) :: self
            real(wp), intent(in) :: gmax, delay
            character(len = 80), intent(in) ::dynamics
            real(wp), intent(in) :: variation, timeConstant

            self%gMaxTot_muS = self%gMaxTot_muS + gmax
            self%numberOfIncomingSynapses = self%numberOfIncomingSynapses + 1
            call AddToList(self%gmax_muS, gmax)
            call AddToList(self%delay_ms, delay)
            call self%dynamics%AddToList(trim(dynamics))
            if (trim(dynamics) == 'Depressing') then
                call AddToList(self%variation, 1.0 - variation)
            else
                call AddToList(self%variation, 1.0 + variation)
            end if
            call AddToList(self%timeConstant_ms, timeConstant)

			!+alt.3! below, to avoid "Segmentation fault - invalid memory reference" at "receiveSpike" subroutine
			call AddToList(self%tBeginOfPulse,    real(-1e4,wp)	         )	!+alt.3!	!self%tBeginOfPulse(:)    = -1e4
			call AddToList(self%tEndOfPulse,      real(-1e4,wp)	         )	!+alt.3!	!self%tEndOfPulse(:)      = -1e4
			call AddToList(self%tLastPulse,       real(-1e4,wp)	         )	!+alt.3!	!self%tLastPulse(:)       = -1e4
			call boolAddToList(self%conductanceState, .false.            )	!+alt.3!	!self%conductanceState(:) = .false.
			call AddToList(self%ri,               real(0.0,wp)	         )	!+alt.3!	!self%ri(:)               = 0.0
			call AddToList(self%ti,               real(0.0,wp)	         )	!+alt.3!	!self%ti(:)               = 0.0
			call AddToList(self%dynamicGmax,      gmax                   ) 	!+alt.3!	!self%dynamicGmax(:)      = self%gmax_muS
			call AddToList(self%synContrib,       gmax / self%gMaxTot_muS)	!+alt.3!	!self%synContrib          = self%gmax_muS / self%gMaxTot_muS

        end subroutine

        subroutine reset(self)
            ! '''

            ! '''
            class(Synapse), intent(inout) :: self


            if (.not.allocated(self%tEndOfPulse)) then
                allocate(self%tBeginOfPulse(self%numberOfIncomingSynapses))
                allocate(self%tEndOfPulse(self%numberOfIncomingSynapses))
                allocate(self%tLastPulse(self%numberOfIncomingSynapses))
                allocate(self%conductanceState(self%numberOfIncomingSynapses))
                allocate(self%ri(self%numberOfIncomingSynapses))
                allocate(self%ti(self%numberOfIncomingSynapses))
                allocate(self%dynamicGmax(self%numberOfIncomingSynapses))
                allocate(self%synContrib(self%numberOfIncomingSynapses))
            end if

            call self%inQueue%clear()
            call self%outQueue%clear()

            if (self%numberOfIncomingSynapses > 0) then
                self%tBeginOfPulse(:) = -1e4
                self%tEndOfPulse(:) = -1e4
                self%tLastPulse(:) = -1e4
                self%conductanceState(:) = .false.
                self%ri(:) = 0.0
                self%ti(:) = 0.0
                self%dynamicGmax(:) = self%gmax_muS
                self%synContrib = self%gmax_muS / self%gMaxTot_muS
            end if

            self%Non = 0.0
            self%Ron = 0.0
            self%Roff = 0.0
            self%t0 = 0.0
            self%lastUpdateInstant = 0.0

        end subroutine

end module SynapseClass

! #@jit
! def compRon(Non, rInf, Ron, t0, t, tauOn):
!     '''
!     Computes the fraction of postsynaptic receptors
!     that are bound to neurotransmitters of all the individual synapses
!     that have neurotransmitters being released (during the pulse).

!     - Inputs:
!         + **Non**: sum of the fractions of the individual conductances that are
!         receiving neurotransmitter (during pulse) relative to
!         the \f$G_{max}\f$ (\f$N_{on}=\limits\sum_{i=1}g_{i_{on}}/G_{max}\f$).

!         + **rInf**: the fraction of postsynaptic receptors that
!         would be bound to neurotransmitters after an infinite
!         amount of time with neurotransmitter being released.

!         + **Ron**: sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that have neurotransmitters being released (during the pulse).

!         + **t0**: instant that the last spike arrived to the compartment.

!         + **t**: current instant, in ms.

!         + **tauOn**: Time constant during a pulse, in ms.
!         \f$\tau_{on}=\frac{1}{\alpha.T_{max} +\beta}\f$.
!     - Outputs:
!         + The fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that have neurotransmitters being released

!     It is computed by the following equation:

!     \f{equation}{
!         R_{on_{newValue}} = N_{on}r_{\infty}\Bigg[1-\exp\left(-\frac{t-t_0}{\tau_{on}}\right)\Bigg] + R_{on_{oldValue}}\exp\left(-\frac{t-t_0}{\tau_{on}}\right)
!     \f}
!     '''
!     return Non * rInf + (Ron - Non * rInf) * np.exp((t0 - t) / tauOn)


! #@jit
! def compRoff(Roff, t0, t, tauOff):
!     '''
!     Computes the fraction of postsynaptic receptors
!     that are bound to neurotransmitters of all the individual synapses
!     that do not have neurotransmitters being released (before and after
!     the pulse).

!     - Inputs:
!         + **Roff**: sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that do not have neurotransmitters being released (before and after
!         the pulse).

!         + **t0**: instant that the last spike arrived to the compartment.

!         + **t**: current instant, in ms.

!         + **tauOff**: time constant after a pulse, in ms.

!     + Output:
!         + The fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that do not have neurotransmitters being released.

!     It is computed by the following formula:

!     \f{equation}{
!         R_{off_{newValue}} = R_{off_{oldValue}}\exp\left(-\frac{t - t0}{\tau_{off}} \right)
!     \f}

!     '''
!     return Roff * np.exp((t0 - t) / tauOff)

! #@jit
! def compRiStart(ri, t, ti, tPeak, tauOff):
!     '''
!     Computes the fraction of bound postsynaptic receptors
!     to neurotransmitters in individual synapses when the
!     neurotransmitter begin (begin of the pulse).

!     - Inputs:
!         + **ri**: the fraction of postsynaptic receptors that
!         were bound to neurotransmitters at the last state change.

!         + **t**: current instant, in ms.

!         + **ti**: The instant that the last pulse began.

!         + **tPeak**: The duration of the pulse.

!         + **tauOff**: Time constant after a pulse, in ms.

!     - Output:
!         + individual synapse state value.

!     It is computed by the following equation:

!     \f{equation}{
!         r_{i_{newValue}} = r_{i_{oldValue}} \exp\left(\frac{t_i+T_{dur}-t}{\tau_{off}}\right)
!     \f}
!     '''
!     return ri * np.exp((ti + tPeak - t) / tauOff)

! #@jit
! def compRiStop(rInf, ri, expFinish):
!     '''
!     Computes the fraction of bound postsynaptic receptors
!     to neurotransmitters in individual synapses when the
!     neurotransmitter release stops (the pulse ends).

!     - Inputs:
!         + **rInf**: the fraction of postsynaptic receptors that
!         would be bound to neurotransmitters after an infinite
!         amount of time with neurotransmitter being released.

!         + **ri**: the fraction of postsynaptic receptors
!         that were bound to neurotransmitters at the last
!         state change.

!         + **expFinish**: Is the value of the exponential at the
!         end of the pulse (\f$\exp(T_{dur}/\tau_{on})\f$). It is
!         is computed before for computational efficiency.

!     - Output:
!         + individual synapse state value.

!     It is computed by the following equation:

!     \f{equation}{
!         r_{i_{newValue}} = r_{\infty} + (r_{i_{oldValue}} - r_{\infty}) \exp\left(\frac{T_{dur}}{\tau_{on}}\right)
!     \f}
!     '''
!     return rInf + (ri - rInf) * expFinish

! #@jit
! def compRonStart(Ron, ri, synContrib):
!     '''
!     Incorporates a new conductance to the set of
!     conductances during a pulse.

!     - Inputs:
!         + **Ron**: sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that have neurotransmitters being released (during the pulse).

!         + **ri**: fraction of postsynaptic receptors that are
!         bound to neurotransmitters of the individual synapses.

!         + **synContrib**: individual conductance constribution
!         to the global synaptic conductance.

!     + Output:
!         + The new value of the sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that have neurotransmitters being released (during the pulse).

!     It is computed as:

!     \f{equation}{
!         R_{on_{newValue}} = R_{on_{oldValue}} + r_iS_{indCont}
!     \f}
!     '''
!     return Ron + np.sum(ri * synContrib)

! #@jit
! def compRoffStart(Roff, ri, synContrib):
!     '''
!     Incorporates a new conductance to the set of
!     conductances that are not during a pulse.

!     - Inputs:
!         + **Roff**: sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that do not have neurotransmitters being released (before and after
!         the pulse).

!         + **ri**: fraction of postsynaptic receptors that are
!         bound to neurotransmitters of the individual synapses.

!         + **synContrib**: individual conductance constribution
!         to the global synaptic conductance.

!     + Output:
!         + The new value of the sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that do not have neurotransmitters being released (before and after
!         the pulse).

!     It is computed as:

!     \f{equation}{
!         R_{off_{newValue}} = R_{off_{oldValue}} - r_iS_{indCont}
!     \f}
!     '''
!     return Roff - np.sum(ri * synContrib)

! #@jit
! def compRonStop(Ron, ri, synContrib):
!     '''
!     Removes a conductance from the set of
!     conductances during a pulse.

!     - Inputs:
!         + **Ron**: sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that have neurotransmitters being released (during the pulse).

!         + **ri**: fraction of postsynaptic receptors that are
!         bound to neurotransmitters of the individual synapses.

!         + **synContrib**: individual conductance constribution
!         to the global synaptic conductance.

!     + Output:
!         + The new value of the sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that have neurotransmitters being released (during the pulse).

!     It is computed as:

!     \f{equation}{
!         R_{on_{newValue}} = R_{on_{oldValue}} - r_iS_{indCont}
!     \f}
!     '''

!     return Ron - np.sum(ri * synContrib)

! #@jit
! def compRoffStop(Roff, ri, synContrib):
!     '''
!     Removes a conductance from the set of
!     conductances that are not during a pulse.

!     - Inputs:
!         + **Roff**: sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that do not have neurotransmitters being released (before and after
!         the pulse).

!         + **ri**: fraction of postsynaptic receptors that are
!         bound to neurotransmitters of the individual synapses.

!         + **synContrib**: individual conductance constribution
!         to the global synaptic conductance.

!     + Output:
!         + The new value of the sum of the fraction of postsynaptic receptors
!         that are bound to neurotransmitters of all the individual synapses
!         that do not have neurotransmitters being released (before and after
!         the pulse).

!     It is computed as:

!     \f{equation}{
!         R_{off_{newValue}} = R_{off_{oldValue}} + r_iS_{indCont}
!     \f}
!     '''
!     return Roff + np.sum(ri * synContrib)
























