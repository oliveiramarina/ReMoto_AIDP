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


module PointProcessGeneratorClass
    use DynamicalArrays
    implicit none    
    private
    integer, parameter :: wp = kind( 1.0d0 )
    public :: PointProcessGenerator

    type PointProcessGenerator
        integer :: index
        real(wp) :: threshold, lastSpike
        real(wp), dimension(:), allocatable :: points

        contains
            procedure :: newThreshold
            procedure :: atualizeGenerator
            procedure :: reset
    end type PointProcessGenerator

    interface PointProcessGenerator
        module procedure init_PointProcessGenerator
    end interface 

    contains

        type(PointProcessGenerator) function init_PointProcessGenerator(index)
            ! '''
            ! Constructor

            ! - Inputs:
            !     + **index**: integer corresponding to the unit order in the pool.
            ! '''
            integer, intent(in) :: index

            
            
            init_PointProcessGenerator%index = index
            call init_PointProcessGenerator%newThreshold(1)
            
            init_PointProcessGenerator%lastSpike = -1e6

            
        end function init_PointProcessGenerator

        subroutine newThreshold(self, gammaOrder)
            ! '''
            ! Generates a number according to a Gamma Distribution with an integer order **GammaOrder**.

            ! - Inputs:
            !     + **GammaOrder**: integer order of the Gamma distribution.

            
            ! - Outputs:
            !     + The number generated from the Gamma distribution.

            ! The number is generated according to:

            ! \f{equation}{
            !     \Gamma = -\frac{1}{\lambda}\ln(\limits\prod_{i=1}^{\lambda} U(0,1))
            ! \f}
            ! where \f$\lambda\f$ is the order of the Gamma distribution and U(a,b) is
            ! a uniform distribution from a to b.

            ! '''
            class(PointProcessGenerator), intent(inout) :: self
            integer, intent(in) :: gammaOrder
            real(wp) :: randomVector(gammaOrder)
            
            call random_number(randomVector)
			!DISPTST1(5, randomVector)
            
            self%threshold = - 1.0/gammaOrder * log(product(randomVector))
        end subroutine newThreshold

        subroutine atualizeGenerator(self, t, firingRate, gammaOrder)
            ! '''

            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **firingRate**: instant firing rate, in spikes/s.
            ! '''
            class(PointProcessGenerator), intent(inout) :: self
            integer, intent(in) :: gammaOrder
            real(wp), intent(in) :: t, firingRate

            if (self%threshold <= 0 .and. (t - self%lastSpike)>0.3) then
                call AddToList(self%points, t)
                call self%newThreshold(gammaOrder)
                self%lastSpike = t
            end if
            self%threshold = self%threshold - firingRate
        end subroutine

        subroutine reset(self)
            class(PointProcessGenerator), intent(inout) :: self
            
            if (allocated(self%points)) deallocate(self%points)
            call self%newThreshold(1)
            self%lastSpike = -1e6
        end subroutine   


end module PointProcessGeneratorClass


     

