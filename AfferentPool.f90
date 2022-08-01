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

module AfferentPoolClass
    ! '''
    ! Class that implements an afferent pool. Encompasses a set of axons.
    ! '''

    use ConfigurationClass
    use AfferentUnitClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    public :: AfferentPool

    type AfferentPool
        character(len = 2) :: poolKind
        type(Configuration), pointer :: conf
        character(len = 6) :: pool, muscle
        integer :: AFnumber
        type(AfferentUnit), dimension(:), allocatable :: unit
        real(wp), dimension(:,:), allocatable :: poolLastCompSpikes, poolTerminalSpikes


        contains
            procedure :: atualizeAfferentPool
            procedure :: listSpikes
            procedure :: reset
    end type AfferentPool

    interface AfferentPool
        module procedure init_AfferentPool
    end interface AfferentPool

    contains

        type(AfferentPool) function init_AfferentPool(conf, pool, muscle)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pool**: string with Motor unit pool to which the motor unit belongs.
            ! '''   
            class(Configuration), intent(in), target :: conf
            character(len = 6), intent(in) :: pool
            character(len = 6), intent(in) :: muscle
            character(len=80) :: paramTag, paramChar
            real(wp) :: paramReal
            integer :: i
			
            ! ## Indicates that is Motor Unit pool.
            init_AfferentPool%poolKind = 'AF'

            ! ## Configuration object with the simulation parameters.
            init_AfferentPool%conf => conf
            ! ## String with Motor unit pool to which the motor unit belongs.
            init_AfferentPool%pool = pool
            
            init_AfferentPool%muscle = muscle

            paramTag = 'Number_' // trim(init_AfferentPool%pool) // '-' // trim(init_AfferentPool%muscle)
			  
            paramChar = init_AfferentPool%conf%parameterSet(paramTag, init_AfferentPool%pool, 0)
            read(paramChar, *)paramReal
           
			init_AfferentPool%AFnumber = nint(paramReal)
				            
            ! ## Dictionary of Axon objects.
            if (allocated(init_AfferentPool%unit)) then
                deallocate(init_AfferentPool%unit)
            end if    
            allocate(init_AfferentPool%unit(init_AfferentPool%AFnumber))
            
            do i = 1, init_AfferentPool%AFnumber
                init_AfferentPool%unit(i) = AfferentUnit(init_AfferentPool%conf, &
                                                         init_AfferentPool%pool, &
                                                         init_AfferentPool%muscle, &
                                                         i)
            end do
            
            ! ## Vector with the instants of spikes in the last dynamical compartment, in ms.
            if (allocated(init_AfferentPool%poolLastCompSpikes)) then
                deallocate(init_AfferentPool%poolLastCompSpikes)                
            end if

            ! ## Vector with the instants of spikes in the terminal, in ms.
            if (allocated(init_AfferentPool%poolTerminalSpikes)) then
                deallocate(init_AfferentPool%poolTerminalSpikes)
            end if

            ! ##
            print '(A)', 'Afferent Pool ' // trim(pool) // ' of muscle ' // trim(muscle) // ' built'

        end function

        subroutine atualizeAfferentPool(self, t, proprioceptorFR)
            ! '''
            ! Update all parts of the Motor Unit pool. It consists
            ! to update all motor units, the activation signal and
            ! the muscle force.

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **proprioceptorFR**: proprioceptor firing rate, in Hz.
            ! '''
            class(AfferentPool), intent(inout) :: self
            real(wp), intent(in) :: t, proprioceptorFR
            integer :: i
            real(wp) :: frequency            

            do i = 1, self%AFnumber
                if (proprioceptorFR >= self%unit(i)%frequencyThreshold_Hz) then 
                    frequency = proprioceptorFR - self%unit(i)%frequencyThreshold_Hz + self%unit(i)%initialFiringRate
                else
                    frequency = 0.0_wp
                end if               
                
                call self%unit(i)%atualizeAfferentUnit(t, frequency*self%conf%timeStep_ms/1000.0)
            end do
        end subroutine

        subroutine listSpikes(self)
            ! '''
            ! List the spikes that occurred in the soma and in
            ! the terminal of the different motor units.
            ! '''
            class(AfferentPool), intent(inout) :: self
            integer :: i
            integer, dimension(self%AFnumber) :: numberOfNewSpikesLastComp, numberOfNewSpikesTerminal
            integer :: numberOfSpikesLastComp, numberOfSpikesTerminal
            integer :: initInd, endInd

            do i = 1, self%AFnumber
                if (allocated(self%unit(i)%lastCompSpikeTrain)) then 
                    numberOfNewSpikesLastComp(i) = size(self%unit(i)%lastCompSpikeTrain)
                else
                    numberOfNewSpikesLastComp(i) = 0
                end if
                
                if (allocated(self%unit(i)%terminalSpikeTrain)) then
                    numberOfNewSpikesTerminal(i) = size(self%unit(i)%terminalSpikeTrain)
                else
                    numberOfNewSpikesTerminal(i) = 0
                end if 
            end do

            
            numberOfSpikesTerminal = sum(numberOfNewSpikesTerminal)
            numberOfSpikesLastComp = sum(numberOfNewSpikesLastComp)

            
            
            
            if (allocated(self%poolLastCompSpikes)) deallocate(self%poolLastCompSpikes)
            if (allocated(self%poolTerminalSpikes)) deallocate(self%poolTerminalSpikes)
            allocate(self%poolLastCompSpikes(numberOfSpikesLastComp,2))
            allocate(self%poolTerminalSpikes(numberOfSpikesTerminal,2))                

            initInd = 1
            do i = 1, self%AFnumber
                if (allocated(self%unit(i)%terminalSpikeTrain)) then
                    endInd = initInd + numberOfNewSpikesTerminal(i) - 1
                    self%poolTerminalSpikes(initInd:endInd,1) = self%unit(i)%terminalSpikeTrain
                    self%poolTerminalSpikes(initInd:endInd,2) = i
                    initInd = endInd+1                        
                end if
            end do 

            initInd = 1
            do i = 1, self%AFnumber
                if (allocated(self%unit(i)%lastCompSpikeTrain)) then
                    endInd = initInd + numberOfNewSpikesLastComp(i) - 1
                    self%poolLastCompSpikes(initInd:endInd,1) = self%unit(i)%lastCompSpikeTrain
                    self%poolLastCompSpikes(initInd:endInd,2) = i
                    initInd = endInd+1                
                end if
            end do 
        end subroutine

        subroutine reset(self)
            ! '''

            ! '''
            class(AfferentPool), intent(inout) :: self
            integer :: i

            if (allocated(self%poolLastCompSpikes)) deallocate(self%poolLastCompSpikes)
            if (allocated(self%poolTerminalSpikes)) deallocate(self%poolTerminalSpikes)
            do i = 1, self%AFnumber
                call self%unit(i)%reset()
            end do

            print '(A)', 'Afferent Pool reseted'
        end subroutine



end module AfferentPoolClass


        

        

    