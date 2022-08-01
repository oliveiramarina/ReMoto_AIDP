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

!     Contact: renato.watanabe@ufabc.edu.br.
! '''
!Versão corrigida
module MuscleHillClass
    use ConfigurationClass
    use MotorUnitClass
    implicit none
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: PI = 4 * atan(1.0_wp)
    private
    public :: MuscleHill


    type MuscleHill
        type(Configuration), pointer :: conf
        character(len = 6) :: pool
        integer :: MUnumber, MUtypeInumber, timeIndex
        real(wp), dimension(:), allocatable :: twTet, twitchAmp_N
        real(wp) :: maximumActivationForce, optimalLength_m, pennationAngleAtOptimalLengthSin
        real(wp) :: maximumForce_N, elasticity, strain, viscosity, mass, tendonElasticity
        real(wp) :: tendonLinearOnsetLength, tendonCurvatureConstant, optimalTendonLength
        real(wp), dimension(:), allocatable :: force, tendonForce_N, contractileForce_N
        real(wp), dimension(:), allocatable :: elasticForce_N, viscousForce_N, length_m
        real(wp), dimension(:), allocatable :: velocity_m_ms, acceleration_m_ms2, tendonLength_m
        real(wp), dimension(:), allocatable :: pennationAngle_rad, activationTypeI, activationTypeII
        real(wp), dimension(:), allocatable :: musculoTendonLength_m, momentArm_m
        character(len = 80) :: paramChar, paramTag, activationModel
        real(wp) :: lengthNorm, velocityNorm, accelerationNorm
        real(wp) :: tendonLengthNorm, forceNorm, tendonForceNorm
        real(wp) :: b_TypeI, b_TypeII, p_TypeI, p_TypeII, w_TypeI, w_TypeII
        real(wp) :: d_TypeI, d_TypeII, a0_TypeI, a0_TypeII, a1_TypeI, a1_TypeII
        real(wp) :: a2_TypeI, a2_TypeII, c0_TypeI, c0_TypeII, c1_TypeI, c1_TypeII
        real(wp) :: Vmax_TypeI, Vmax_TypeII
        real(wp) :: m0, m1, m2, m3, m4, n0, n1, n2, n3, n4


        contains
            procedure :: atualizeForce
            procedure :: computeForceLengthTypeI
            procedure :: atualizeActivation
            procedure :: computePennationAngle
            procedure :: computeForceLengthTypeII
            procedure :: computeForceVelocityTypeI
            procedure :: computeForceVelocityTypeII
            procedure :: computeAcceleration
            procedure :: dLdt
            procedure :: atualizeMuscleForce
            procedure :: atualizeTendonForce
            procedure :: computeElasticElementForce
            procedure :: computeViscousElementForce
            procedure :: computeTypeIActiveForce
            procedure :: computeTypeIIActiveForce
            procedure :: atualizeLenghtsAndVelocity
            procedure :: atualizeMusculoTendonLength
            procedure :: atualizeMomentArm
            procedure :: reset

    end type MuscleHill

    interface MuscleHill
        module procedure init_MuscleHill
    end interface MuscleHill

    contains

        type(MuscleHill) function init_MuscleHill(conf, pool, MUnumber, MUtypeInumber, unit)
            class(Configuration), intent(in), target :: conf
            character(len = 6), intent(in) :: pool
            integer, intent(in) :: MUnumber, MUtypeInumber
            class(MotorUnit), dimension(MUnumber), intent(in) :: unit
            integer :: i, timeLength
            real(wp) :: optimalPennationAngle
            character(len = 80) ::paramChar, paramTag


            ! ##
            init_MuscleHill%conf => conf
            ! ##
            init_MuscleHill%pool = pool
            ! ##
            init_MuscleHill%MUnumber = MUnumber
            ! ##
            init_MuscleHill%MUtypeInumber = MUtypeInumber
            ! ##
            init_MuscleHill%timeIndex = 1

            ! ## Twitch-tetanus relationship (see atualizeForce function explanation)
            allocate(init_MuscleHill%twTet(init_MuscleHill%MUnumber))
            init_MuscleHill%twTet(:) = 0.0

            ! ## Amplitude of the muscle unit twitch, in N (see atualizeForce function explanation).
            allocate(init_MuscleHill%twitchAmp_N(init_MuscleHill%MUnumber))
            init_MuscleHill%twitchAmp_N(:) = 0.0

            !!!Modificado - Marina - Modelo Raikova
			paramTag = 'activationModel'
            init_MuscleHill%activationModel = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            do i = 1, init_MuscleHill%MUnumber
                init_MuscleHill%twitchAmp_N(i) = unit(i)%TwitchAmp_N
				if (trim(init_MuscleHill%activationModel).eq.'SOCDS') then
					init_MuscleHill%twTet(i) = unit(i)%twTet
					else
					init_MuscleHill%twTet(i) = unit(i)%twTetRaikova
				end if
            end do

            ! ## This is used for normalization purposes. It is the maximum force that
            ! ## the muscle reach when the Hill model is not used.
            init_MuscleHill%maximumActivationForce = dot_product(init_MuscleHill%twitchAmp_N, &
                                                                 init_MuscleHill%twTet)
            timeLength = int(conf%simDuration_ms/conf%timeStep_ms)
            ! ## Muscle force along time, in N.
            allocate(init_MuscleHill%force(timeLength))
            init_MuscleHill%force(:) = 0.0
            ! ##
            allocate(init_MuscleHill%tendonForce_N(timeLength))
            init_MuscleHill%tendonForce_N(:) = 0.0
            ! ##
            allocate(init_MuscleHill%contractileForce_N(timeLength))
            init_MuscleHill%contractileForce_N(:) = 0.0
            ! ##
            allocate(init_MuscleHill%elasticForce_N(timeLength))
            init_MuscleHill%elasticForce_N(:) = 0.0
            ! ##
            allocate(init_MuscleHill%viscousForce_N(timeLength))
            init_MuscleHill%viscousForce_N(:) = 0.0
            ! ##
            allocate(init_MuscleHill%length_m(timeLength))
            init_MuscleHill%length_m(:) = 0.0
            ! ##
            allocate(init_MuscleHill%velocity_m_ms(timeLength))
            init_MuscleHill%velocity_m_ms(:) = 0.0
            ! ##
            allocate(init_MuscleHill%acceleration_m_ms2(timeLength))
            init_MuscleHill%acceleration_m_ms2(:) = 0.0
            ! ##
            allocate(init_MuscleHill%tendonLength_m(timeLength))
            init_MuscleHill%tendonLength_m(:) = 0.0
            ! ##
            allocate(init_MuscleHill%pennationAngle_rad(timeLength))
            init_MuscleHill%pennationAngle_rad(:) = 0.0
            ! ##
            allocate(init_MuscleHill%activationTypeI(timeLength))
            init_MuscleHill%activationTypeI(:) = 0.0
            ! ##
            allocate(init_MuscleHill%activationTypeII(timeLength))
            init_MuscleHill%activationTypeII(:) = 0.0
            ! ##
            allocate(init_MuscleHill%musculoTendonLength_m(timeLength))
            init_MuscleHill%musculoTendonLength_m(:) = 0.0
            ! ##
            allocate(init_MuscleHill%momentArm_m(timeLength))
            init_MuscleHill%momentArm_m(:) = 0.0
            ! ##
            paramTag = 'optimalMuscleLength:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%optimalLength_m
            ! ##
            paramtag = 'optimalPennationAngle:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)optimalPennationAngle
            init_MuscleHill%pennationAngleAtOptimalLengthSin = sin(optimalPennationAngle)
            ! ## Maximum force of the Hill model, in N.
            paramTag = 'Fmax:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%maximumForce_N
            ! ##
            paramTag = 'muscleElasticity:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%elasticity
            ! ##
            paramtag = 'muscleStrain:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%strain
            ! ##
            paramTag = 'muscleViscosity:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%viscosity
            ! ##
            paramTag = 'muscleMass:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%mass
            ! ##
            paramTag = 'initialMuscleLength:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%length_m(1)
            ! ##
            paramTag = 'tendonElasticity:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%tendonElasticity
            ! ##
            paramTag = 'tendonLinearOnsetLength:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%tendonLinearOnsetLength
            ! ##
            paramTag = 'tendonCurvatureConstant:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%tendonCurvatureConstant
            ! ##
            paramTag = 'optimalTendonLength:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%optimalTendonLength

            ! ##
            init_MuscleHill%lengthNorm = init_MuscleHill%length_m(1)/init_MuscleHill%optimalLength_m 
			! ##
            init_MuscleHill%velocityNorm = 0.0
            ! ##
            init_MuscleHill%accelerationNorm = 0.0

            init_MuscleHill%pennationAngle_rad(1) = init_MuscleHill%computePennationAngle()
            init_MuscleHill%tendonLength_m(1) = init_MuscleHill%musculoTendonLength_m(1) - init_MuscleHill%length_m(1) * cos(init_MuscleHill%pennationAngle_rad(1))


            init_MuscleHill%tendonLengthNorm = init_MuscleHill%tendonLength_m(1)/init_MuscleHill%optimalTendonLength
            ! ##
            init_MuscleHill%forceNorm = 0.0
            ! ##
            init_MuscleHill%tendonForceNorm = 0.0

            ! ##
            paramTag = 'bTypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%b_TypeI
            ! ##
            paramTag = 'bTypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%b_TypeII
            ! ##
            paramTag =  'pTypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%p_TypeI
            ! ##
            paramTag = 'pTypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%p_TypeII
            ! ##
            paramTag = 'wTypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%w_TypeI
            ! ##
            paramTag =  'wTypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%w_TypeII
            ! ##
            paramTag = 'dTypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%d_TypeI
            ! ##
            paramTag = 'dTypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%d_TypeII
            ! ##
            paramTag = 'a0TypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%a0_TypeI
            ! ##
            paramTag = 'a0TypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%a0_TypeII
            ! ##
            paramTag = 'a1TypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%a1_TypeI
            ! ##
            paramTag = 'a1TypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%a1_TypeII
            !
            paramTag = 'a2TypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%a2_TypeI
            !##
            paramTag = 'a2TypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%a2_TypeII
            !##
            paramTag = 'c0TypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%c0_TypeI
            !##
            paramTag = 'c0TypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%c0_TypeII
            ! ##
            paramTag = 'c1TypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%c1_TypeI
            !##
            paramTag = 'c1TypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%c1_TypeII
            !##
            paramTag = 'VmaxTypeI:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%Vmax_TypeI
            !##
            paramTag = 'VmaxTypeII:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%Vmax_TypeII
            !##
            paramTag = 'm0:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%m0
            !##
            paramTag = 'm1:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%m1
            ! ##
            paramTag = 'm2:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%m2
            ! ##
            paramTag = 'm3:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%m3
            !##
            paramTag = 'm4:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%m4
            !##
            paramTag = 'n0:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%n0
            !##
            paramTag = 'n1:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%n1
            !##
            paramTag = 'n2:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%n2
            !##
            paramTag = 'n3:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%n3
            ! ##
            paramTag = 'n4:' // trim(pool)
            paramChar = init_MuscleHill%conf%parameterSet(paramTag, pool, 0)
            read(paramChar, *)init_MuscleHill%n4

            print '(A, F15.2, A)', 'Hill muscle of the ' // trim(pool) // &
            ' muscle with maximum force of ',  init_MuscleHill%maximumForce_N,  ' N  built.'

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
            class(MuscleHill), intent(inout) :: self
            real(wp), dimension(self%MUnumber) :: activation_Sat

            call self%atualizeLenghtsAndVelocity()

            self%lengthNorm = self%length_m(self%timeIndex) / self%optimalLength_m
            self%velocityNorm = self%velocity_m_ms(self%timeIndex) / self%optimalLength_m
            self%accelerationNorm = self%acceleration_m_ms2(self%timeIndex) / self%optimalLength_m
            self%pennationAngle_rad(self%timeIndex) = self%computePennationAngle()
            self%tendonLength_m(self%timeIndex) = self%musculoTendonLength_m(self%timeIndex) - self%length_m(self%timeIndex) * cos(self%pennationAngle_rad(self%timeIndex))
            self%tendonLengthNorm = self%tendonLength_m(self%timeIndex) / self%optimalTendonLength

            call self%atualizeActivation(activation_Sat)

            call self%atualizeMuscleForce()
            call self%atualizeTendonForce()

            self%force(self%timeIndex) = self%forceNorm * self%maximumForce_N
            self%tendonForce_N(self%timeIndex) = self%tendonForceNorm * self%maximumForce_N
            self%timeIndex = self%timeIndex + 1
        end subroutine

        subroutine atualizeActivation(self, activation_Sat)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self
            real(wp), dimension(self%MUnumber) :: activation_Sat
            self%activationTypeI(self%timeIndex) = sum(activation_Sat(1:self%MUtypeInumber) * self%twitchAmp_N(1:self%MUtypeInumber) * self%twTet(1:self%MUtypeInumber)) / self%maximumActivationForce
			self%activationTypeII(self%timeIndex) = sum(activation_Sat(self%MUtypeInumber+1:) * self%twitchAmp_N(self%MUtypeInumber+1:) * self%twTet(self%MUtypeInumber+1:)) / self%maximumActivationForce
		end subroutine

        real(wp) function computePennationAngle(self) result(pennationAngle)
            ! '''
			! 	!M! Definição
			!	FONTE: WIKIPEDIA
			! 	Pennation angle is the angle between the longitudinal axis of the entire muscle and its fibers.
			!	The longitudinal axis is the force generating axis of the muscle and pennate fibers lie at an oblique angle.
			!	As tension increases in the muscle fibers, the pennation angle also increases.
            ! '''
            class(MuscleHill), intent(inout) :: self

            pennationAngle = asin(self%pennationAngleAtOptimalLengthSin / self%lengthNorm)
        end function

        real(wp) function computeForceLengthTypeI(self) result(fl1curve)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            fl1curve = exp(-(abs(self%lengthNorm ** self%b_TypeI - 1.0) / self%w_TypeI) ** self%p_TypeI)
			end function

        real(wp) function computeForceLengthTypeII(self) result(fl2curve)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            fl2curve = exp(-(abs(self%lengthNorm ** self%b_TypeII - 1.0) / self%w_TypeII) ** self%p_TypeII)
		end function

        real(wp) function computeForceVelocityTypeI(self) result(fv)  
            ! '''													  
            ! '''														
            class(MuscleHill), intent(inout) :: self

            if 	(1000*self%velocityNorm > 0.1) then	
				fv = (self%d_TypeI - (self%a0_TypeI+self%a1_TypeI * self%lengthNorm + self%a2_TypeI * self%lengthNorm ** 2)*(1000*self%velocityNorm))/(self%d_TypeI + (1000*self%velocityNorm))
			else	
				fv = (self%Vmax_TypeI - (1000*self%velocityNorm)) / (self%Vmax_TypeI + (1000*self%velocityNorm) * &	
                (self%c0_TypeI + self%c1_TypeI * self%lengthNorm))		
            end if
			
			

        end function

        real(wp) function computeForceVelocityTypeII(self) result(fv) 
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            if 	(1000*self%velocityNorm > 0.1) then	
				fv = (self%d_TypeII - (self%a0_TypeII + self%a1_TypeII * self%lengthNorm + &
                self%a2_TypeII * self%lengthNorm**2)*(1000*self%velocityNorm)) / (self%d_TypeII + (1000*self%velocityNorm))
            else
				fv = (self%Vmax_TypeII - (1000*self%velocityNorm)) / (self%Vmax_TypeII + (1000*self%velocityNorm) * &
                (self%c0_TypeII + self%c1_TypeII * self%lengthNorm))
            end if
        end function
!
        real(wp) function computeAcceleration(self) result(acceleration)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            if (self%timeIndex>1) then
				self%acceleration_m_ms2(self%timeIndex) =  (self%tendonForce_N(self%timeIndex-1) -& 
                    self%force(self%timeIndex-1) * cos(self%pennationAngle_rad(self%timeIndex-1))) /& 
                    (self%mass * cos(self%pennationAngle_rad(self%timeIndex-1)))/ 1.0e6

            else
                self%acceleration_m_ms2(self%timeIndex) = 0.0
            end if

            acceleration = self%acceleration_m_ms2(self%timeIndex)

        end function


        function dLdt(self) result(derivative)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self
            real(wp), dimension(2) :: derivative

            derivative =  [self%velocity_m_ms(self%timeIndex-1), self%computeAcceleration()]
	        end function


        subroutine atualizeMuscleForce(self)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            self%forceNorm  = (self%computeElasticElementForce() + self%computeViscousElementForce() + self%computeTypeIActiveForce() + self%computeTypeIIActiveForce())

        end subroutine

        subroutine atualizeTendonForce(self)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self
            self%tendonForceNorm = self%tendonCurvatureConstant * self%tendonElasticity * log(exp((self%tendonLengthNorm - self%tendonLinearOnsetLength) / self%tendonCurvatureConstant) + 1.0)
        end subroutine

        real(wp) function computeElasticElementForce(self) result(elasticForce) !M! Equação 4.6 Dissertação Renato Naville Watanabe
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            elasticForce = exp(self%elasticity / self%strain * (self%lengthNorm - self%strain - 1.0))
		   end function

        real(wp) function computeViscousElementForce(self) result(viscousForce) !M! Equação 4.8 Dissertação Renato Naville Watanabe
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            viscousForce = self%viscosity * self%velocityNorm
        end function

        real(wp) function computeTypeIActiveForce(self) result(activeForce1) !Fce
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            activeForce1 = (self%activationTypeI(self%timeIndex) * self%computeForceLengthTypeI() * self%computeForceVelocityTypeI())
			
		end function

        real(wp) function computeTypeIIActiveForce(self) result(activeForce2)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self

            activeForce2 = (self%activationTypeII(self%timeIndex) * self%computeForceLengthTypeII() * &
                        self%computeForceVelocityTypeII())
		end function

        subroutine atualizeLenghtsAndVelocity(self)
            ! '''
            ! '''
            class(MuscleHill), intent(inout) :: self
            real(wp), dimension(2) :: derivLdt



            if (self%timeIndex>1) then
                derivLdt = self%dLdt()
                self%length_m(self%timeIndex) = self%length_m(self%timeIndex-1) + self%conf%timeStep_ms * derivLdt(1)
                self%velocity_m_ms(self%timeIndex) = self%velocity_m_ms(self%timeIndex-1) + self%conf%timeStep_ms * derivLdt(2)
            end if
        end subroutine

        subroutine atualizeMusculoTendonLength(self, ankleAngle)
            ! '''
            ! Computes the musculotendon complex length.
            !
            ! Inputs:
            !
            ! * ankleAngle - the angle of the ankle, in degrees (real).
            ! '''
            class(MuscleHill), intent(inout) :: self
            real(wp), intent(in) :: ankleAngle
			!M! Apêndice A Dissertação Renato Naville Watanabe
			if (size(self%musculoTendonLength_m) >= self%timeIndex) then
				self%musculoTendonLength_m(self%timeIndex) = (self%m0 + self%m1 * ankleAngle + self%m2 * (ankleAngle ** 2) + &
															self%m3 * (ankleAngle ** 3) + self%m4 * (ankleAngle ** 4))
			end if

        end subroutine

        subroutine atualizeMomentArm(self, ankleAngle)
            ! '''
            ! Computes the moment-arm of the muscle.
            !
            ! Inputs:
            !
            ! * ankleAngle - the angle of the ankle, in degrees (real).
            ! '''
            class(MuscleHill), intent(inout) :: self
            real(wp), intent(in) :: ankleAngle
			if (size(self%momentArm_m) >= self%timeIndex) then
				self%momentArm_m(self%timeIndex) = (self%n0 + self%n1 * ankleAngle + self%n2 * (ankleAngle ** 2) + &
													self%n3 * (ankleAngle ** 3) + self%n4 * (ankleAngle ** 4))
			end if
        end subroutine

        subroutine reset(self)
            ! '''

            ! '''
            class(MuscleHill), intent(inout) :: self
            character(len = 80) :: paramTag, paramChar


            self%timeIndex = 1
            self%tendonForce_N(:) = 0.0
            self%contractileForce_N(:) = 0.0
            self%elasticForce_N(:) = 0.0
            self%viscousForce_N(:) = 0.0
            self%length_m(:) = 0.0
            self%velocity_m_ms(:) = 0.0
            self%tendonLength_m(:) = 0.0
            self%pennationAngle_rad(:) = 0.0
            self%activationTypeI(:) = 0.0
            self%activationTypeII(:) = 0.0
            self%musculoTendonLength_m(:) = 0.0
            self%momentArm_m(:) = 0.0

            self%velocityNorm = 0.0
            self%forceNorm = 0.0
            self%tendonForceNorm = 0.0

            paramTag = 'initialMuscleLength:' // trim(self%pool)
            paramChar = self%conf%parameterSet(paramTag, self%pool, 0)
            read(paramChar, *)self%length_m(1)

            self%lengthNorm = self%length_m(1)/self%optimalLength_m

            self%pennationAngle_rad(1) = self%computePennationAngle()

            self%tendonLength_m(1) = self%musculoTendonLength_m(1) - self%length_m(1) * cos(self%pennationAngle_rad(1))
            self%tendonLengthNorm = self%tendonLength_m(1)/self%optimalTendonLength

        end subroutine
!
end module MuscleHillClass


