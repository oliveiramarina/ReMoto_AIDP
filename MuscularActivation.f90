! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2020  Renato Naville Watanabe
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

module MuscularActivationClass
    use ConfigurationClass
    use DynamicalArrays
    use MotorUnitClass
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: PI = 4 * atan(1.0_wp)  	
    public :: MuscularActivation

    type MuscularActivation
        type(Configuration), pointer :: conf
        character(len = 6) :: pool
        integer :: MUnumber
        character(len = 80) :: activationModel
        real(wp), dimension(:,:), allocatable :: ActMatrix
        real(wp), dimension(:), allocatable :: an, activation_nonSat, k, p, m, bSat, bSatRaikova,activation_Sat
        real(wp), dimension(:), allocatable :: activationRaikova
		real(wp), dimension(:), allocatable :: diracDeltaValue
        type(MotorUnit), pointer:: unit(:)
       
        contains
            procedure :: atualizeActivationSignal
            procedure :: reset
    end type MuscularActivation
    
    interface MuscularActivation
        module procedure init_MuscularActivation
    end interface MuscularActivation

    contains

        type(MuscularActivation) function init_MuscularActivation(conf, pool, MUnumber, unit)
            class(Configuration), intent(in), target :: conf        
            character(len = 6), intent(in) :: pool
            integer, intent(in) :: MUnumber
            class(MotorUnit), dimension(MUnumber),  intent(in), target:: unit
            character(len = 80) :: paramTag
            integer :: i, j, stat

            init_MuscularActivation%conf => conf
            init_MuscularActivation%pool = pool
            init_MuscularActivation%MUnumber = MUnumber
            nullify(init_MuscularActivation%unit)
            init_MuscularActivation%unit => unit(:)
            
            !## Model of the activation signal. For now, it can be *SOCDS* (second order critically damped system) or Raikova.
            paramTag = 'activationModel'
            init_MuscularActivation%activationModel = init_MuscularActivation%conf%parameterSet(paramTag, pool, 0)

            if (trim(init_MuscularActivation%activationModel).eq.'SOCDS') then
                allocate(init_MuscularActivation%ActMatrix(init_MuscularActivation%MUnumber, 3*init_MuscularActivation%MUnumber))
                init_MuscularActivation%ActMatrix(:,:) = 0.0
				
                do i = 1, init_MuscularActivation%MUnumber
                    init_MuscularActivation%ActMatrix(i,3*(i-1)+1:3*(i-1)+3) = &
                                    [2*exp(-init_MuscularActivation%conf%timeStep_ms/init_MuscularActivation%unit(i)%TwitchTc_ms), &
										-exp(-2*init_MuscularActivation%conf%timeStep_ms/init_MuscularActivation%unit(i)%TwitchTc_ms), & 
										(init_MuscularActivation%conf%timeStep_ms**2)/init_MuscularActivation%unit(i)%TwitchTc_ms*exp(1.0-init_MuscularActivation%conf%timeStep_ms/init_MuscularActivation%unit(i)%TwitchTc_ms)]
			   end do
                
                allocate(init_MuscularActivation%an(3*init_MuscularActivation%MUnumber))
                init_MuscularActivation%an(:) = 0.0
				
				allocate(init_MuscularActivation%bSat(init_MuscularActivation%MUnumber))
				do i = 1, init_MuscularActivation%MUnumber     
                init_MuscularActivation%bSat(i) = init_MuscularActivation%unit(i)%bSat
				end do
				init_MuscularActivation%diracDeltaValue = - init_MuscularActivation%bSat / init_MuscularActivation%conf%timeStep_ms
				
			print *, 'SOCDS activation model built'  
            else 
				
				allocate(init_MuscularActivation%k(init_MuscularActivation%MUnumber))
				init_MuscularActivation%k(:) = 0.0
				allocate(init_MuscularActivation%p(init_MuscularActivation%MUnumber))
				init_MuscularActivation%p(:) = 0.0
				allocate(init_MuscularActivation%m(init_MuscularActivation%MUnumber))
				init_MuscularActivation%m(:) = 0.0
				allocate(init_MuscularActivation%bSatRaikova(init_MuscularActivation%MUnumber))
				init_MuscularActivation%bSatRaikova(:) = 0.0
				allocate(init_MuscularActivation%activationRaikova(init_MuscularActivation%MUnumber))
				init_MuscularActivation%activationRaikova(:) = 0.0

				
				do i = 1, init_MuscularActivation%MUnumber
					init_MuscularActivation%k(i) =  log(2.0)/(-init_MuscularActivation%unit(i)%TwitchTcRaikova_ms * log((init_MuscularActivation%unit(i)%TwitchHrRaikova_ms)/(init_MuscularActivation%unit(i)%TwitchTcRaikova_ms)) + init_MuscularActivation%unit(i)%TwitchHrRaikova_ms - init_MuscularActivation%unit(i)%TwitchTcRaikova_ms)
					
					init_MuscularActivation%p(i) = exp(-init_MuscularActivation%k(i) * init_MuscularActivation%unit(i)%TwitchTcRaikova_ms * log(init_MuscularActivation%unit(i)%TwitchTcRaikova_ms) - (-init_MuscularActivation%k(i) * init_MuscularActivation%unit(i)%TwitchTcRaikova_ms))					
									
					init_MuscularActivation%m(i) = init_MuscularActivation%k(i) * init_MuscularActivation%unit(i)%TwitchTcRaikova_ms
					
					init_MuscularActivation%bSatRaikova(i) = init_MuscularActivation%unit(i)%bSatRaikova
					
				end do				
				print *, 'Raikova activation model built'                        			
			end if
			
            ! ## The non-saturated activation signal of all motor units (see actMatrix explanation).
            allocate(init_MuscularActivation%activation_nonSat(init_MuscularActivation%MUnumber))    
            init_MuscularActivation%activation_nonSat(:) = 0.0
            ! ## The parameter \f$b\f$ (see twitchSaturation function explanation) of 
            ! ## each motor unit.
            
                       
            ! ## The non-saturated activation signal of all motor units (see actMatrix explanation).
            allocate(init_MuscularActivation%activation_Sat(init_MuscularActivation%MUnumber))
            init_MuscularActivation%activation_Sat(:) = 0.0
            ! ## Dirac's delta approximation amplitude value. Is the inverse
            ! ## of the simulation time step (\f$1/T\f$). 
			
        end function

        subroutine atualizeActivationSignal(self, t)
            ! '''
            ! Update the activation signal of the motor units.
            !
            ! - Inputs:
            !     + **t**: current instant, in ms.        
            ! '''
            class(MuscularActivation), intent(inout) :: self
            real(wp), intent(in) :: t
			real(wp), dimension(self%MUnumber) :: temp!
            integer :: i, j, sizeTrain, stat            
			
            if (trim(self%activationModel).eq.'SOCDS') then
				do i = 1, self%MUnumber
					self%an(3*(i-1)+2) = self%an(3*(i-1)+1)
					self%an(3*(i-1)+1) = self%activation_nonSat(i)
					if (allocated(self%unit(i)%terminalSpikeTrain)) then                    
						sizeTrain = size(self%unit(i)%terminalSpikeTrain)				
						if (abs(t - self%unit(i)%terminalSpikeTrain(sizeTrain)) < 1e-6) then
							self%an(3*(i-1)+3) = self%diracDeltaValue(i)
						else
							self%an(3*(i-1)+3) = 0.0
						end if   
					else
						self%an(3*(i-1)+3) = 0.0
					end if
					
				end do  
			
				self%activation_nonSat = matmul(self%ActMatrix, self%an)
			else
				do i = 1, self%MUnumber
					temp(i) = 0
					self%activationRaikova(i) = 0
					if (allocated(self%unit(i)%terminalSpikeTrain)) then                    
						sizeTrain = size(self%unit(i)%terminalSpikeTrain)  
						do j = 1, sizeTrain
						if(t >= 0 .and. (t-self%unit(i)%terminalSpikeTrain(j))< 10.0*self%unit(i)%TwitchHrRaikova_ms) then
								temp(i) = self%p(i) * (t-self%unit(i)%terminalSpikeTrain(j))**self%m(i) * exp(-self%k(i) * (t-self%unit(i)%terminalSpikeTrain(j)));
								self%activationRaikova(i) = temp(i) + self%activationRaikova(i)
								end if							
						end do
					end if
					
				end do
				self%activation_nonSat =  - self%bSatRaikova * self%activationRaikova
				end if
				
            self%activation_Sat(:) = 2.0 / (1.0 + exp(self%activation_nonSat)) - 1.0
	
        end subroutine

        subroutine reset(self)
            ! '''
			!
            ! '''
            class(MuscularActivation), intent(inout) :: self
            
           if (trim(self%activationModel).eq.'SOCDS') then
				self%an(:) = 0.0            
           else
				self%activationRaikova(:) = 0.0
			end if
			self%activation_nonSat(:) = 0.0
			self%activation_Sat(:) = 0.0
        end subroutine

end module MuscularActivationClass

    
    
