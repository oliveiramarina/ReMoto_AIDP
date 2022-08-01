! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2020  Renato Naville Watanabe
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

module MuscleNoHillClass
    use ConfigurationClass
    use MotorUnitClass
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: pi = 4 * atan(1.0_wp) 
    public :: MuscleNoHill

    type MuscleNoHill
        type(Configuration), pointer :: conf
        character(len = 6) :: pool
        integer :: MUnumber, MUtypeInumber
        real(wp), dimension(:), allocatable :: twTet, twitchAmp_N, maximumActivationForce
        real(wp), dimension(:), allocatable :: force, momentArm_m
        integer :: timeIndex, simDurationSteps 
        real(wp) :: lengthNorm, velocityNorm, accelerationNorm  
        real(wp) :: m0, m1, m2, m3, m4
        real(wp) :: n0, n1, n2, n3, n4
		character(len = 80) :: activationModel

        contains
            procedure :: atualizeMusculoTendonLength
            procedure :: atualizeForce
            procedure :: atualizeMomentArm
            procedure :: reset

    end type MuscleNoHill

    interface MuscleNoHill
        module procedure init_MuscleNoHill
    end interface MuscleNoHill

    contains

        type(MuscleNoHill) function init_MuscleNoHill(conf, pool, MUnumber, MUtypeInumber, unit)
            class(Configuration), intent(in), target :: conf
            character(len = 6), intent(in) :: pool
            integer, intent(in) :: MUnumber, MUtypeInumber
            class(MotorUnit), dimension(MUnumber), intent(in) :: unit
            integer :: i
            character(len = 80) :: paramChar, paramTag

            init_MuscleNoHill%conf => conf
            init_MuscleNoHill%pool = pool
            init_MuscleNoHill%MUnumber = MUnumber
            init_MuscleNoHill%MUtypeInumber = MUtypeInumber

            ! ## Twitch- tetanus relationship (see atualizeForceNoHill function explanation)
            allocate(init_MuscleNoHill%twTet(init_MuscleNoHill%MUnumber))
            ! ## Amplitude of the muscle unit twitch, in N (see atualizeForceNoHill function explanation).
            allocate(init_MuscleNoHill%twitchAmp_N(init_MuscleNoHill%MUnumber))

			paramTag = 'activationModel'
            init_MuscleNoHill%activationModel = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)

            do i = 1, init_MuscleNoHill%MUnumber
                init_MuscleNoHill%twitchAmp_N(i) = unit(i)%TwitchAmp_N
				if (trim(init_MuscleNoHill%activationModel).eq.'SOCDS') then
					init_MuscleNoHill%twTet(i) = unit(i)%twTet  
				else
					init_MuscleNoHill%twTet(i) = unit(i)%twTetRaikova
				end if
            end do
            ! ## This is used for normalization purposes. It is the maximum force that
            ! ## the muscle reach when the Hill model is not used.
            allocate(init_MuscleNoHill%maximumActivationForce(init_MuscleNoHill%MUnumber))
            init_MuscleNoHill%maximumActivationForce = init_MuscleNoHill%twitchAmp_N * init_MuscleNoHill%twTet
            ! ## Muscle force along time, in N.
            init_MuscleNoHill%simDurationSteps = nint(init_MuscleNoHill%conf%simDuration_ms/init_MuscleNoHill%conf%timeStep_ms)    
			print *, init_MuscleNoHill%simDurationSteps, 'simDurationSteps'
            allocate(init_MuscleNoHill%force(init_MuscleNoHill%simDurationSteps))
            init_MuscleNoHill%force(:) = 0.0

            init_MuscleNoHill%timeIndex = 1


            ! ##
            init_MuscleNoHill%lengthNorm = 1.0
            init_MuscleNoHill%velocityNorm = 0.0
            init_MuscleNoHill%accelerationNorm = 0.0

            !##
            
            allocate(init_MuscleNoHill%momentArm_m(init_MuscleNoHill%simDurationSteps))

            ! ##
            paramTag = 'm0:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleNoHill%m0
            ! ##
            paramTag = 'm1:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar,*)init_MuscleNoHill%m1
            ! ##
            paramTag = 'm2:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleNoHill%m2
            ! ## 
            paramTag =  'm3:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleNoHill%m3
            ! ##  
            paramTag = 'm4:' // trim(pool)    
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleNoHill%m4
            ! ##
            paramTag = 'n0:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0) 
            read(paramChar, *)init_MuscleNoHill%n0
            ! ##
            paramTag = 'n1:' // trim(pool)  
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleNoHill%n1
            ! ## 
            paramTag = 'n2:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0) 
            read(paramChar, *)init_MuscleNoHill%n2
            ! ## 
            paramTag = 'n3:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0) 
            read(paramChar, *)init_MuscleNoHill%n3
            ! ##  
            paramTag = 'n4:' // trim(pool)
            paramChar = init_MuscleNoHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleNoHill%n4
        end function

        subroutine atualizeForce(self, activation_Sat)
            ! '''
            ! Compute the muscle force when no muscle dynamics (Hill model) is used. This
            ! operation is vectorized. Each element of the vectors correspond to one motor
            ! unit. For each motor unit, the force is computed by the following formula:

            ! \f{equation}{
            !     F_{MU} = a_{sat}A_{MU}R_{MU}
            ! }
            ! where \f$a_{sat}\f$ is the saturated activation signal, \f$A_{MU}\f$ is the
            ! motor unit twitch amplitude, and  \f$R_{MU}\f$ is the relation between 
            ! the twitch amplitude and the tetanus of the motor unit. 

            ! Then the muscle force is obtained from: 

            ! \f{equation}{
            !     F = \limits\sum_{i=1}^N_{MU}F_{i}
            ! }
            ! where \f$N_{MU}\f$ is the number of motor units in the pool.
            ! '''
            class(MuscleNoHill), intent(inout) :: self
            real(wp), dimension(self%MUnumber), intent(in) :: activation_Sat


            self%force(self%timeIndex) = dot_product(activation_Sat, self%maximumActivationForce)
            self%timeIndex = self%timeIndex + 1
        end subroutine

        subroutine atualizeMusculoTendonLength(self, ankleAngle)
            ! ''' 
            ! '''
            class(MuscleNoHill), intent(inout) :: self
            real(wp), intent(in) :: ankleAngle
        end subroutine

        subroutine atualizeMomentArm(self, ankleAngle)
            ! '''
            ! '''
            class(MuscleNoHill), intent(inout) :: self
            real(wp), intent(in) :: ankleAngle
			if (size(self%momentArm_m) >= self%timeIndex) then 
				self%momentArm_m(self%timeIndex) = self%n0 + self%n1*ankleAngle + self%n2*(ankleAngle**2) + self%n3*(ankleAngle**3) + self%n4*(ankleAngle**4)
			end if
        end subroutine


        subroutine reset(self)
            ! '''

            ! '''
            class(MuscleNoHill), intent(inout) :: self

            self%force(:) = 0.0
            self%timeIndex = 1
        end subroutine      


end module MuscleNoHillClass



    
    

   
    
        