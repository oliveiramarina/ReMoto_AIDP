! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2021  Renato Naville Watanabe
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

module jointAnkleForceTaskClass
    use ConfigurationClass
    use MotorUnitPoolClass
    use MusclePointerClass
    implicit none
    private
    integer, parameter :: wp = kind( 1.0d0 )
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    public :: jointAnkleForceTask

    type jointAnkleForceTask
        type(Configuration), pointer :: conf
        type(MusclePointer), dimension(:), allocatable :: muscles
        real(wp), dimension(:), allocatable :: ankleAngle_rad, ankleTorque_Nm
        contains
            procedure :: atualizeAngle
            procedure :: atualizeAnkle
            procedure :: computeTorque
            procedure :: reset
    end type jointAnkleForceTask

    interface jointAnkleForceTask   
        module procedure init_jointAnkleForceTask
    end interface jointAnkleForceTask

    contains

        type(jointAnkleForceTask) function   init_jointAnkleForceTask(conf, pools)
            class(Configuration), intent(in), target :: conf
            class(MotorUnitPool), intent(in) :: pools(:)
            integer :: numberOfPools, i, timeLength
            init_jointAnkleForceTask%conf => conf
            
            numberOfPools = size(pools)
            allocate(init_jointAnkleForceTask%muscles(numberOfPools))

            do i = 1, numberOfPools
                if (pools(i)%pool == 'SOL' .or. pools(i)%pool == 'MG' .or. pools(i)%pool == 'LG' .or. pools(i)%pool == 'TA') then
                    init_jointAnkleForceTask%muscles(i) = MusclePointer()
                    call init_jointAnkleForceTask%muscles(i)%assignMuscle(pools(i))
                end if
            end do
            ! ##
            timeLength = nint(conf%simDuration_ms/conf%timeStep_ms)
            allocate(init_jointAnkleForceTask%ankleAngle_rad(timeLength))
            init_jointAnkleForceTask%ankleAngle_rad(:) = 0.0
            allocate(init_jointAnkleForceTask%ankleTorque_Nm(timeLength))
            init_jointAnkleForceTask%ankleTorque_Nm(:) = 0.0

            
            print '(A)', 'Ankle joint for Force Task built'
			

        end function 

        subroutine atualizeAnkle(self, t, ankleAngle)
            ! '''
            ! Update the ankle joint.
            
            ! - Inputs:
            !     + **t**: current instant, in ms.

            !     + **ankleAngle**: ankle angle, in rad. 
            ! '''
            class(jointAnkleForceTask), intent(inout) :: self
            real(wp), intent(in) ::t, ankleAngle
            integer :: i

            call self%atualizeAngle(t, ankleAngle)

            if (self%muscles(1)%muscle%hillModel == 'No') then
                do i = 1, size(self%muscles)
                    call self%muscles(i)%muscle%NoHillMuscle%atualizeMusculoTendonLength(ankleAngle)
                    call self%muscles(i)%muscle%NoHillMuscle%atualizeMomentArm(ankleAngle)
                end do
            else 
                do i = 1, size(self%muscles)
                    call self%muscles(i)%muscle%HillMuscle%atualizeMusculoTendonLength(ankleAngle)
                    call self%muscles(i)%muscle%HillMuscle%atualizeMomentArm(ankleAngle)
                end do
            end if

            
        end subroutine

        subroutine atualizeAngle(self, t, ankleAngle)
            ! '''
            ! '''
            class(jointAnkleForceTask), intent(inout) :: self
            real(wp), intent(in) ::t, ankleAngle
            integer :: timeIndex

            timeIndex = nint(t/self%conf%timeStep_ms)+1
            
            self%ankleAngle_rad(timeIndex) = ankleAngle
        end subroutine

    subroutine computeTorque(self, t)
        ! '''
        ! '''
        class(jointAnkleForceTask), intent(inout) :: self
        real(wp), intent(in) ::t
        real(wp) :: torque, velocity, acceleration
        integer :: i
        integer :: timeIndex

        timeIndex = nint(t/self%conf%timeStep_ms)+1
            
        torque = 0.0

        if (self%muscles(1)%muscle%hillModel == 'No') then
            do i = 1, size(self%muscles) 
                torque = torque + self%muscles(i)%muscle%NoHillMuscle%force(timeIndex) * self%muscles(i)%muscle%NoHillMuscle%momentArm_m(timeIndex)
            end do
        else
            do i = 1, size(self%muscles)
                torque = torque + self%muscles(i)%muscle%HillMuscle%force(timeIndex) * self%muscles(i)%muscle%HillMuscle%momentArm_m(timeIndex)
            end do
        end if

        if (timeIndex > 1) then
            velocity = (self%ankleAngle_rad(timeIndex) - self%ankleAngle_rad(timeIndex - 1))/self%conf%timeStep_ms
        else
            velocity = 0.0
        end if

        if (timeIndex > 2) then
			!M! Como foi calculado no JAVA?
            acceleration = (self%ankleAngle_rad(timeIndex) - 2*self%ankleAngle_rad(timeIndex - 1) + self%ankleAngle_rad(timeIndex - 2))/(self%conf%timeStep_ms**2)
        else
            acceleration = 0.0 
        end if
        
        torque = torque - 1100*velocity  - 320*self%ankleAngle_rad(timeIndex) + 0.007*acceleration
        
        self%ankleTorque_Nm(timeIndex) = torque

    end subroutine

    subroutine reset(self)
        ! '''
        ! '''
        class(jointAnkleForceTask), intent(inout) :: self
        
        self%ankleAngle_rad(:) = 0.0
        self%ankleTorque_Nm(:) = 0.0
    end subroutine
end module jointAnkleForceTaskClass



    
    
        
    

        
        
    