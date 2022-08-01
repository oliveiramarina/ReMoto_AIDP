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

module GolgiTendonOrganClass
    ! '''
    ! Class that implements a muscle spindle model. 
    ! '''
    use ConfigurationClass
    implicit none
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    private
    public :: GolgiTendonOrgan

    type GolgiTendonOrgan
        type(Configuration), pointer:: conf
        character(len = 6) :: muscle
        real(wp) :: IbFR_Hz
        real(wp) :: x1, x2, Gg

        contains
            procedure :: atualizeGolgiTendonOrgan
            procedure :: reset
            
    end type GolgiTendonOrgan

    interface GolgiTendonOrgan
        module procedure init_GolgiTendonOrgan
    end interface GolgiTendonOrgan

    contains

    type(GolgiTendonOrgan) function init_GolgiTendonOrgan(conf, muscle)
        ! '''
        ! Constructor

        ! - Inputs:
        !     + **conf**: Configuration object with the simulation parameters.

        !     + **muscle**: string with the muscle to which the muscle spindle belongs.              
        ! '''
        class(Configuration), intent(in), target :: conf
        character(len = 6), intent(in) :: muscle       
        character(len = 80) :: paramTag, paramChar
        real(wp) :: AFnumber

        ! ## Configuration object with the simulation parameters.
        init_GolgiTendonOrgan%conf => conf

        init_GolgiTendonOrgan%muscle = muscle

        init_GolgiTendonOrgan%x1 = 0.0

        init_GolgiTendonOrgan%x2 = 0.0
                
        init_GolgiTendonOrgan%IbFR_Hz = 0.0

        paramTag = 'Number_Ib-' // trim(muscle)
        paramChar = init_GolgiTendonOrgan%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)AFnumber
                
        init_GolgiTendonOrgan%Gg = 60.0 / AFnumber !Hz

        print '(A)', 'Golgi Tendon Organ from muscle ' // trim(init_GolgiTendonOrgan%muscle) // ' built.'

    end function

    subroutine atualizeGolgiTendonOrgan(self, t, tendonForce)
        ! '''
        ! Atualize the dynamical and nondynamical (delay) parts of the motor unit.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(GolgiTendonOrgan), intent(inout) :: self
        real(wp), intent(in) :: t, tendonForce
        real(wp) :: Gf, Rgn, RgfSlope, RgfDotSlope, Rgf
        
        Gf = 4.0 !unidade N

        Rgn = self%Gg*log(tendonForce/Gf + 1)  !Equação 19 Tese de Doutorado LAE

        RgfSlope = self%x2
        RgfDotSlope = -0.4 * self%x1 - 2.2 * self%x2 + Rgn
        
        self%x1 = self%x1 + self%conf%timeStep_ms * RgfSlope
        self%x2 = self%x2 + self%conf%timeStep_ms * RgfDotSlope

        Rgf = -11.2 * self%x1 - 46.4 * self%x2 + 68 * Rgn

        if (Rgf > 0) then
            self%IbFR_Hz = Rgf
        else 
            self%IbFR_Hz = 0.0
        end if
        
    end subroutine
    
    subroutine reset(self)
        ! '''
        ! '''s
        class(GolgiTendonOrgan), intent(inout) :: self

        self%x1 = 0.0

        self%x2 = 0.0
    
        self%IbFR_Hz = 0.0
        
    end subroutine

end module GolgiTendonOrganClass