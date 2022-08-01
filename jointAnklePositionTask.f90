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

module jointAnklePositionTaskClass
    use ConfigurationClass
    use MotorUnitPoolClass
    use MusclePointerClass
    implicit none
    private
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp), parameter :: g = 9.80655
    public :: jointAnklePositionTask
    

    type jointAnklePositionTask
        type(Configuration), pointer :: conf
        type(MusclePointer), dimension(:), allocatable :: muscles
        real(wp), dimension(:), allocatable :: ankleAngle_rad, ankleTorque_Nm, ankleOmega_rad_s, muscleTorque_Nm, passiveTorque_Nm, gravitationalTorque_Nm, Torque_Nm, load
        real(wp) :: footInertia, footMass, ankleFootDistance, jointViscosity, jointElasticity, maximumVoluntaryTorque, percentageOfMVC
    
        contains
            procedure :: atualizeAnkle
            procedure :: computeTorque
            procedure :: reset

    end type jointAnklePositionTask

    interface jointAnklePositionTask   
        module procedure init_jointAnklePositionTask
    end interface jointAnklePositionTask

    contains

        type(jointAnklePositionTask) function   init_jointAnklePositionTask(conf, pools)
            class(Configuration), intent(in), target :: conf
            class(MotorUnitPool), intent(in) :: pools(:)
            integer :: numberOfPools, i, timeLength

            init_jointAnklePositionTask%conf => conf
            
            numberOfPools = size(pools)
            allocate(init_jointAnklePositionTask%muscles(numberOfPools))

            do i = 1, numberOfPools
                if (pools(i)%pool == 'SOL' .or. pools(i)%pool == 'MG' .or. pools(i)%pool == 'LG' .or. pools(i)%pool == 'TA') then
                    init_jointAnklePositionTask%muscles(i) = MusclePointer()
                    call init_jointAnklePositionTask%muscles(i)%assignMuscle(pools(i))
                end if
            end do
            ! ##
            timeLength = nint(conf%simDuration_ms/conf%timeStep_ms)
            allocate(init_jointAnklePositionTask%ankleAngle_rad(timeLength))
            init_jointAnklePositionTask%ankleAngle_rad(:) = 5.0*pi/180.0 
            allocate(init_jointAnklePositionTask%ankleTorque_Nm(timeLength))
            init_jointAnklePositionTask%ankleTorque_Nm(:) = 0.0
			allocate(init_jointAnklePositionTask%muscleTorque_Nm(timeLength))
            init_jointAnklePositionTask%muscleTorque_Nm(:) = 0.0
			allocate(init_jointAnklePositionTask%passiveTorque_Nm(timeLength))
            init_jointAnklePositionTask%passiveTorque_Nm(:) = 0.0
			allocate(init_jointAnklePositionTask%gravitationalTorque_Nm(timeLength))
            init_jointAnklePositionTask%gravitationalTorque_Nm(:) = 0.0
			allocate(init_jointAnklePositionTask%Torque_Nm(timeLength))
            init_jointAnklePositionTask%Torque_Nm(:) = 0.0
			allocate(init_jointAnklePositionTask%load(timeLength))
            init_jointAnklePositionTask%load(:) = 0.0
            allocate(init_jointAnklePositionTask%ankleOmega_rad_s(timeLength))
            init_jointAnklePositionTask%ankleOmega_rad_s(:) = 0.0
																									

			!REF constantes: C. Maurer and R. J. Peterka, “A new interpretation of spontaneous
								!sway measures based on a simple model of human postural control,”
								!Journal of Neurophysiology, vol. 93, no. 1, pp. 189–200, 2005.
			
            init_jointAnklePositionTask%footMass = 2.01 					!kg															
            init_jointAnklePositionTask%ankleFootDistance = 0.052  			![m]
			init_jointAnklePositionTask%jointViscosity = 1.1				![N.m.s.rad-1]
			init_jointAnklePositionTask%jointElasticity = 320 				![Nm.rad-1]			--  0.65*init_jointAnklePositionTask%footMass*g*init_jointAnklePositionTask%ankleFootDistance			
            init_jointAnklePositionTask%footInertia = (4.0/3.0)*(init_jointAnklePositionTask%footMass*init_jointAnklePositionTask%ankleFootDistance**2) 	![kg.m2] (foot inertia)
			init_jointAnklePositionTask%maximumVoluntaryTorque = 12.00
			init_jointAnklePositionTask%percentageOfMVC = conf%percentageOfMVC
			print*, conf%percentageOfMVC, init_jointAnklePositionTask%percentageOfMVC
            print '(A)', 'Position Task built'

        end function 

        subroutine atualizeAnkle(self, t)
            ! '''
            ! Update the ankle joint.
            ! Atualizes the musculotendon length and the moment-arm of each muscle.
            ! Updates the angle and angular velocity by numerically solving 
            ! the differential equations with the Euler method.
            ! 
            ! Inputs:
            !
            ! * t: real 
            ! 
            
            class(jointAnklePositionTask), intent(inout) :: self
            real(wp), intent(in) ::t
            integer :: i, timeIndex
            real(wp) :: dthetadt, domegadt, angle


            timeIndex = nint(t/self%conf%timeStep_ms)+1
            
            angle = self%ankleAngle_rad(timeIndex)*180/pi
            if (self%muscles(1)%muscle%hillModel == 'No') then
                do i = 1, size(self%muscles)
                    
                    call self%muscles(i)%muscle%NoHillMuscle%atualizeMusculoTendonLength(angle)
                    call self%muscles(i)%muscle%NoHillMuscle%atualizeMomentArm(angle)
                end do
            else 
                do i = 1, size(self%muscles)
                    call self%muscles(i)%muscle%HillMuscle%atualizeMusculoTendonLength(angle)
                    call self%muscles(i)%muscle%HillMuscle%atualizeMomentArm(angle)
                end do
            end if
            
            call self%computeTorque(t)

            if (t > 1000.0) then
				if (timeIndex < size(self%ankleOmega_rad_s)) then 
					dthetadt = self%ankleOmega_rad_s(timeIndex)		
					domegadt = (self%ankleTorque_Nm(timeIndex))/self%footInertia 
					self%ankleOmega_rad_s(timeIndex + 1) = self%ankleOmega_rad_s(timeIndex) + self%conf%timeStep_ms*domegadt/1000.0
					self%ankleAngle_rad(timeIndex + 1) = self%ankleAngle_rad(timeIndex) + self%conf%timeStep_ms*dthetadt/1000.0
				end if
			else    
				if (timeIndex < size(self%ankleOmega_rad_s)) then
					self%ankleOmega_rad_s(timeIndex + 1) = 0.0
					self%ankleAngle_rad(timeIndex + 1) = 5.0*pi/180.0
				end if
            end if 
            
        end subroutine


        subroutine computeTorque(self, t)
            ! '''
            ! '''
            class(jointAnklePositionTask), intent(inout) :: self
            real(wp), intent(in) ::t
            real(wp) :: muscularTorque, velocity, acceleration, passiveTorque, gravitationalTorque, viscosity, disturbance, phase
            integer :: i
            integer :: timeIndex

            timeIndex = nint(t/self%conf%timeStep_ms) + 1
                
            muscularTorque = 0.0

            if (self%muscles(1)%muscle%hillModel == 'No') then
                do i = 1, size(self%muscles)
                    muscularTorque = muscularTorque + self%muscles(i)%muscle%NoHillMuscle%force(timeIndex) * self%muscles(i)%muscle%NoHillMuscle%momentArm_m(timeIndex)
                end do
            else
                do i = 1, size(self%muscles)
                    muscularTorque = muscularTorque + self%muscles(i)%muscle%HillMuscle%force(timeIndex) * self%muscles(i)%muscle%HillMuscle%momentArm_m(timeIndex)
                end do
            end if
            
			call random_number(phase) 
			disturbance =  0.005*sin(2*pi*50.0/200.0*t/1000 + 2*pi*phase)
			   passiveTorque = (self%jointViscosity*self%ankleOmega_rad_s(timeIndex)) + self%jointElasticity *(self%ankleAngle_rad(timeIndex)-0*pi/180)
               gravitationalTorque = self%footMass*g*self%ankleFootDistance*sin(pi/4-self%ankleAngle_rad(timeIndex))
               viscosity = self%jointViscosity*self%ankleOmega_rad_s(timeIndex)
			   self%ankleTorque_Nm(timeIndex) = muscularTorque - passiveTorque  - gravitationalTorque + self%percentageOfMVC* self%maximumVoluntaryTorque *cos(self%ankleAngle_rad(timeIndex))! + disturbance !Alt - percentageMVC*MVC
			
			self%muscleTorque_Nm(timeIndex) = muscularTorque
			self%Torque_Nm(timeIndex) = muscularTorque - passiveTorque - gravitationalTorque
			self%passiveTorque_Nm(timeIndex) = passiveTorque 
			self%gravitationalTorque_Nm(timeIndex) = gravitationalTorque
			self%load(timeIndex) = self%percentageOfMVC* self%maximumVoluntaryTorque *cos(self%ankleAngle_rad(timeIndex))
			end subroutine

		
        subroutine reset(self)
            ! '''
            ! '''
            class(jointAnklePositionTask), intent(inout) :: self
            
            self%ankleAngle_rad(:) = 5.0*pi/180.0
            self%ankleTorque_Nm(:) = 0.0
            self%ankleOmega_rad_s(:) = 0.0
        end subroutine
		
		        
		
			
end module jointAnklePositionTaskClass        
    