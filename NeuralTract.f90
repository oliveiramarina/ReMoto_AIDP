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

module NeuralTractClass
    ! '''
    ! Class that implements a a neural tract, composed by the descending
    ! commands from the motor cortex.
    ! '''
    use NeuralTractUnitClass
    use ConfigurationClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    public :: NeuralTract

    type NeuralTract
        type(NeuralTractUnit), dimension(:), allocatable :: unit
        type(Configuration), pointer :: conf
        character(len = 6) :: poolKind, pool
        integer :: Number
        real(wp) , dimension(:,:), allocatable :: poolTerminalSpikes
        
        contains
            procedure :: atualizePool
            procedure :: listSpikes
            procedure :: reset
    end type NeuralTract

    interface NeuralTract
        module procedure init_NeuralTract
    end interface

    contains

        type(NeuralTract) function init_NeuralTract(conf, pool)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pool**: string with the name of the Neural tract.
            ! '''
            
            character(len = 6), intent(in) :: pool
            class(Configuration), intent(in), target :: conf
            character(len=80) :: paramTag, paramChar
            real(wp) :: paramReal
            integer :: i
            
            
            init_NeuralTract%conf => conf
            ! ## Indicates that is a neural tract.
            init_NeuralTract%poolKind = 'NT'
            ! ## String with the name of the Neural tract.
            init_NeuralTract%pool = pool
            ! ## The number of neural tract units.
            paramTag = 'Number_' // pool
            paramChar = init_NeuralTract%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)paramReal
            init_NeuralTract%Number = int(paramReal)
            
            ! ## List of NeuralTRactUnit objects.
            allocate(init_NeuralTract%unit(init_NeuralTract%Number))
            
            do i = 1, init_NeuralTract%Number
                init_NeuralTract%unit(i) = NeuralTractUnit(pool, i)
            end do         
            
            print '(A)', 'Descending Command ' // pool // ' built'    
                
        end function init_NeuralTract

        subroutine atualizePool(self, t, FR, GammaOrder)
            ! '''
            ! Update all neural tract units from the neural tract.
            
            ! - Inputs:
            !     + **t**: current instant, in ms.
            !     + **FR**: firing rate, in Hz     
            !     + **GammaOrder** : order of the Gamma distribution, integer      
            ! '''
            class(NeuralTract), intent(inout) :: self
            real(wp), intent(in) :: t, FR
            integer, intent(in) :: GammaOrder
            integer :: i
            real(wp) :: FiringRate

            FiringRate = FR*self%conf%timeStep_ms/1000.0

            do i = 1, self%Number 
                call self%unit(i)%atualizeNeuralTractUnit(t, FiringRate , GammaOrder)
            end do
            
        end subroutine

        subroutine listSpikes(self)
            ! '''
            ! List the spikes that occurred in neural tract units.
            ! '''
            class(NeuralTract), intent(inout) :: self
            integer :: i
            integer :: numberOfSpikes, initInd, endInd
            integer, dimension(self%Number) :: numberOfNewSpikes
            
            
            do i = 1, self%Number
                numberOfNewSpikes(i) = size(self%unit(i)%spikesGenerator%points)
            end do
            numberOfSpikes = sum(numberOfNewSpikes)

                      
            allocate(self%poolTerminalSpikes(numberOfSpikes,2))
            initInd = 1
            do i = 1, self%Number
                endInd = initInd + size(self%unit(i)%spikesGenerator%points) - 1
                self%poolTerminalSpikes(initInd:endInd,1) = self%unit(i)%spikesGenerator%points
                self%poolTerminalSpikes(initInd:endInd,2) = i
                initInd = endInd+1                
            end do                                
        end subroutine
    
        subroutine reset(self)
            class(NeuralTract), intent(inout) :: self
            integer :: i
            
            do i = 1, self%Number
                call self%unit(i)%reset()
            end do
            

            if (allocated(self%poolTerminalSpikes)) deallocate(self%poolTerminalSpikes)

            print *, 'Neural tract reseted'
            
        end subroutine
end module NeuralTractClass
