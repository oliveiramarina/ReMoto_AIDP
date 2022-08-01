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

module SynapticNoiseClass
    use NeuralTractUnitClass
    use ConfigurationClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    public :: SynapticNoise

    type SynapticNoise
        type(Configuration), pointer :: conf
        character(len = 6) :: poolKind, pool
        integer :: Number, GammaOrder, timeIndex
        type(NeuralTractUnit), dimension(:), allocatable :: unit
        real(wp), dimension(:,:), allocatable :: poolTerminalSpikes
        real(wp) :: FR

        contains
            procedure :: atualizePool
            procedure :: listSpikes
            procedure :: reset

    end type SynapticNoise

    interface SynapticNoise
        module procedure init_SynapticNoise
    end interface SynapticNoise

    contains

        type(SynapticNoise) function init_SynapticNoise(conf, pool)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **conf**: Configuration object with the simulation parameters.

            !     + **pool**: string with the name of the pool.
            ! '''
            class(Configuration), intent(in), target :: conf
            character(len = 6), intent(in) :: pool
            character(len = 80) :: paramTag, paramChar
            real(wp) :: paramReal
            integer :: i
            
            init_SynapticNoise%conf => conf
            init_SynapticNoise%poolKind = 'SN'

            ! ## String with the name of the pool.
            init_SynapticNoise%pool = pool
            ! ## The number of neural tract units.
            paramTag = 'Number_' // trim(pool)
            paramChar = init_SynapticNoise%conf%parameterSet(paramTag, pool, 1)
            read(paramChar, *)paramReal
            init_SynapticNoise%Number = nint(paramReal)
            
            ! ## List of NeuralTractUnit objects.
            allocate(init_SynapticNoise%unit(init_SynapticNoise%Number))
            
            paramTag = 'NoiseGammaOrder_' // trim(pool)
            paramChar = init_SynapticNoise%conf%parameterSet(paramTag, pool, 1)
            read(paramChar, *)paramReal
            init_SynapticNoise%GammaOrder = nint(paramReal)

            do i = 1, init_SynapticNoise%Number
                init_SynapticNoise%unit(i) = NeuralTractUnit(pool, i)
            end do
            ! ## Vector with the instants of spikes in the terminal, in ms.
            
            if (allocated(init_SynapticNoise%poolTerminalSpikes)) then 
                deallocate(init_SynapticNoise%poolTerminalSpikes)
            end if
            
            
            ! ## The  mean firing rate of the neural tract units.
            paramTag = 'NoiseFunction_' // trim(pool)
            paramChar = init_SynapticNoise%conf%parameterSet(paramTag, pool, 1)
            read(paramChar, *)paramReal
            init_SynapticNoise%FR = paramReal * init_SynapticNoise%conf%timeStep_ms/1000.0

        
            ! ##
            init_SynapticNoise%timeIndex = 1
            ! ##
            print '(A)', 'Synaptic Noise on ' // trim(pool) // ' built'
        end function

        subroutine atualizePool(self, t)
            ! '''
            ! Update all neural tract units from the neural tract.
            
            ! - Inputs:
            !     + **t**: current instant, in ms.
            ! '''
            class(SynapticNoise), intent(inout) :: self
            real(wp), intent(in) :: t
            integer :: i    
                     

            do i = 1, self%Number
                call self%unit(i)%atualizeNeuralTractUnit(t, self%FR, self%GammaOrder)
            end do
            self%timeIndex = self%timeIndex + 1
        end subroutine

        subroutine listSpikes(self)
            ! '''
            ! List the spikes that occurred in neural tract units.
            ! '''
            class(SynapticNoise), intent(inout) :: self
            integer :: i
            integer :: numberOfSpikes, initInd, endInd
            integer, dimension(self%Number) :: numberOfNewSpikes
            
            
            do i = 1, self%Number
                if (allocated(self%unit(i)%spikesGenerator%points)) then 
                    numberOfNewSpikes(i) = size(self%unit(i)%spikesGenerator%points)
                else 
                    numberOfNewSpikes(i) = 0
                end if
            end do
            numberOfSpikes = sum(numberOfNewSpikes)
                      
            allocate(self%poolTerminalSpikes(numberOfSpikes,2))
            initInd = 1
            do i = 1, self%Number
                if (numberOfNewSpikes(i) > 0) then 
                    endInd = initInd + numberOfNewSpikes(i) - 1
                    self%poolTerminalSpikes(initInd:endInd,1) = self%unit(i)%spikesGenerator%points
                    self%poolTerminalSpikes(initInd:endInd,2) = i
                    initInd = endInd+1                
                end if
            end do        
        end subroutine

        subroutine reset(self)
            class(SynapticNoise), intent(inout) :: self
            integer :: i
            do i = 1, self%Number
                call self%unit(i)%reset()    
            end do
        end subroutine
        

end module SynapticNoiseClass


    
    