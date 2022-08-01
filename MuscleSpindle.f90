! '''
!     Neuromuscular simulator in Fortran.
!     Copyright (C) 2018  Renato Naville Watanabe
!
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

module MuscleSpindleClass
    ! '''
    ! Class that implements a muscle spindle model. 
    ! '''
    use ConfigurationClass
    implicit none
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    private
    public :: MuscleSpindle

    type MuscleSpindle
        type(Configuration), pointer:: conf
        character(len = 6) :: muscle
        real(wp) :: beta0Bag1, beta0Bag2, beta0Chain
        real(wp) :: beta1Bag1, beta2Bag2, beta2Chain
        real(wp) :: GAMMA1Bag1, GAMMA2Bag2, GAMMA2Chain
        real(wp) :: freq_bag1_Hz, freq_bag2_Hz, freq_Chain_Hz
        real(wp) :: tauBag1_s, tauBag2_s
        real(wp) :: KsrBag1, KsrBag2, KsrChain
        real(wp) :: MBag1, MBag2, MChain
        real(wp) :: L0SrBag1, L0SrBag2, L0SrChain
        real(wp) :: L0PrBag1, L0PrBag2, L0PrChain
        real(wp) :: LNSrBag1, LNSrBag2, LNSrChain
        real(wp) :: LNPrBag2, LNPrChain
        real(wp) :: RBag1, RBag2, RChain
        real(wp) :: KPrBag1, KPrBag2, KPrChain
        real(wp) :: GPrimaryBag1, GPrimaryBag2, GPrimaryChain
        real(wp) :: GSecondaryBag1, GSecondaryBag2, GSecondaryChain
        real(wp) :: SOcclusionFactor
        real(wp) :: betaBag1, betaBag2, betaChain
        real(wp) :: GAMMABag1, GAMMABag2, GAMMAChain
        real(wp) :: primaryPotentialBag1, primaryPotentialBag2, primaryPotentialChain
        real(wp) :: secondaryPotentialBag1, secondaryPotentialBag2, secondaryPotentialChain
        real(wp), dimension(3) :: fusimotorActivation
        real(wp), dimension(6) :: fiberTension
        real(wp) :: IaFR_Hz, IIFR_Hz    
        real(wp) :: XBag2, XChain
        real(wp) :: LSecondaryBag2, LSecondaryChain

        contains
            procedure :: atualizeMuscleSpindle
            procedure :: computeFusimotorActivation
            procedure :: computePrimaryActivity
            procedure :: computeIa
            procedure :: dfdt
            procedure :: computeFiberTension
            procedure :: dTdt
            procedure :: computeSecondaryActivity
            procedure :: computeII
            procedure :: reset
            
    end type MuscleSpindle

    interface MuscleSpindle
        module procedure init_MuscleSpindle
    end interface MuscleSpindle

    contains

    type(MuscleSpindle) function init_MuscleSpindle(conf, muscle)
        ! '''
        ! Constructor

        ! - Inputs:
        !     + **conf**: Configuration object with the simulation parameters.

        !     + **muscle**: string with the muscle to which the muscle spindle belongs.              
        ! '''
        class(Configuration), intent(in), target :: conf
        character(len = 6), intent(in) :: muscle
        character(len = 80) :: paramChar, paramTag

        ! ## Configuration object with the simulation parameters.
        init_MuscleSpindle%conf => conf

        init_MuscleSpindle%muscle = muscle
        
        paramTag = 'beta0Bag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%beta0Bag1

        paramTag = 'beta0Bag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%beta0Bag2

        paramTag = 'beta0Chain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%beta0Chain

        paramTag = 'beta1Bag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%beta1Bag1

        paramTag = 'beta2Bag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%beta2Bag2

        paramTag = 'beta2Chain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%beta2Chain
        
        paramTag = 'GAMMA1Bag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GAMMA1Bag1
        
        paramTag = 'GAMMA2Bag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GAMMA2Bag2
        
        paramTag = 'GAMMA2Chain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GAMMA2Chain

        paramTag = 'freqBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%freq_bag1_Hz
        
        paramTag = 'freqBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%freq_bag2_Hz
        
        paramTag = 'freqChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%freq_Chain_Hz

        paramTag = 'tauBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%tauBag1_s
        
        paramTag = 'tauBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%tauBag2_s

        paramTag = 'KsrBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%KsrBag1
        
        paramTag = 'KsrBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%KsrBag2
        
        paramTag = 'KsrChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%KsrChain

        paramTag = 'MBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%MBag1
        
        paramTag = 'MBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%MBag2
        
        paramTag = 'MChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%MChain

        paramTag = 'L0SrBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%L0SrBag1
        
        paramTag = 'L0SrBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%L0SrBag2
        
        paramTag = 'L0SrChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%L0SrChain

        paramTag = 'L0PrBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%L0PrBag1
        
        paramTag = 'L0PrBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%L0PrBag2
        
        paramTag = 'L0PrChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%L0PrChain

        paramTag = 'LNSrBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LNSrBag1
        
        paramTag = 'LNSrBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LNSrBag2
        
        paramTag = 'LNSrChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LNSrChain

        paramTag = 'LNPrBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LNPrBag2
        
        paramTag = 'LNPrChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LNPrChain

        paramTag = 'RBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%RBag1
        
        paramTag = 'RBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%RBag2
        
        paramTag = 'RChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%RChain

        paramTag = 'KPrBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%KPrBag1
        
        paramTag = 'KPrBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%KPrBag2
        
        paramTag = 'KPrChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%KPrChain

        paramTag = 'GPrimaryBag1'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GPrimaryBag1
        
        paramTag = 'GPrimaryBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GPrimaryBag2
        
        paramtag = 'GPrimaryChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GPrimaryChain

        paramTag = 'GSecondaryBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GSecondaryBag2
        
        paramtag = 'GSecondaryChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%GSecondaryChain

        paramTag = 'SOcclusionFactor'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%SOcclusionFactor

        paramTag = 'XBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%XBag2

        paramTag = 'XChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%XChain

        paramTag = 'LSecondaryBag2'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LSecondaryBag2

        paramTag = 'LSecondaryChain'
        paramChar = init_MuscleSpindle%conf%parameterSet(paramTag, muscle, 0)
        read(paramChar, *)init_MuscleSpindle%LSecondaryChain

        init_MuscleSpindle%betaBag1 = 0.0
        init_MuscleSpindle%betaBag2 = 0.0
        init_MuscleSpindle%betaChain = 0.0

        init_MuscleSpindle%GAMMABag1 = 0.0
        init_MuscleSpindle%GAMMABag2 = 0.0
        init_MuscleSpindle%GAMMAChain = 0.0

        init_MuscleSpindle%primaryPotentialBag1 = 0.0
        init_MuscleSpindle%primaryPotentialBag2 = 0.0
        init_MuscleSpindle%primaryPotentialChain = 0.0

        init_MuscleSpindle%secondaryPotentialBag1 = 0.0
        init_MuscleSpindle%secondaryPotentialBag2 = 0.0
        init_MuscleSpindle%secondaryPotentialChain = 0.0

        

        ! ## Vector with the activation of each fusimotor fiber. The first
        ! # element is the frequency of Bag1, the second of Bag2 and the
        ! # third of the Chain.  
        init_MuscleSpindle%fusimotorActivation = [0.0, 0.0, 0.0]

        ! ## Vector with the tensions and tensions derivatives of each 
        ! # fusimotor fiber. The first two elements correspond to the Bag1,
        ! #  the third and fourth to Bag2 and the last ones to the Chain. 
        init_MuscleSpindle%fiberTension = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        init_MuscleSpindle%IaFR_Hz = 0.0
        init_MuscleSpindle%IIFR_Hz = 0.0

        print '(A)', 'Muscle spindle from muscle ' // trim(init_MuscleSpindle%muscle) // ' built.'


    end function

    subroutine atualizeMuscleSpindle(self, t, fascicleLength, fascicleVelocity, &
                              fascicleAcceleration, gammaMNDynamicFR, gammaMNStaticFR)
        ! '''
        ! Atualize the dynamical and nondynamical (delay) parts of the motor unit.

        ! - Inputs:
        !     + **t**: current instant, in ms.
        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t, fascicleLength, fascicleVelocity, fascicleAcceleration
        real(wp), intent(in) :: gammaMNDynamicFR, gammaMNStaticFR

        call self%computeFusimotorActivation(t, gammaMNDynamicFR, gammaMNStaticFR)
                
        self%betaBag1 = self%beta0Bag1 + self%beta1Bag1 * self%fusimotorActivation(1)
        self%betaBag2 = self%beta0Bag2 + self%beta2Bag2 * self%fusimotorActivation(2)
        self%betaChain = self%beta0Chain + self%beta2Chain * self%fusimotorActivation(3)
        
        self%GAMMABag1 = self%GAMMA1Bag1 * self%fusimotorActivation(1)
        self%GAMMABag2 = self%GAMMA2Bag2 * self%fusimotorActivation(2)
        self%GAMMAChain = self%GAMMA2Chain * self%fusimotorActivation(3)
        
        call self%computeFiberTension(t, fascicleLength, fascicleVelocity, fascicleAcceleration)
        
        self%IaFR_Hz = self%computeIa(t)
        self%IIFR_Hz = self%computeII(t, fascicleLength)        
    end subroutine

    subroutine computeFusimotorActivation(self, t, gammaMNDynamicFR, gammaMNStaticFR)
        ! '''

        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t, gammaMNDynamicFR, gammaMNStaticFR
        real(wp), dimension(2) :: df


        df = self%dfdt(t, gammaMNDynamicFR, gammaMNStaticFR)
        self%fusimotorActivation(1) = self%fusimotorActivation(1) + self%conf%timeStep_ms / 1000.0 * df(1) 
        self%fusimotorActivation(2) = self%fusimotorActivation(2) + self%conf%timeStep_ms / 1000.0 * df(2)
        self%fusimotorActivation(3) = gammaMNStaticFR**2/(gammaMNStaticFR**2 + self%freq_Chain_Hz**2)
    end subroutine

    function dfdt(self, t, gammaMNDynamicFR, gammaMNStaticFR) result(df)
        ! '''
        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t, gammaMNDynamicFR, gammaMNStaticFR
        real(wp), dimension(2) :: df
        
        

        df(1) = ((gammaMNDynamicFR**2)/(gammaMNDynamicFR**2 + self%freq_bag1_Hz**2) - self%fusimotorActivation(1)) / self%tauBag1_s
        df(2) = ((gammaMNStaticFR**2)/(gammaMNStaticFR**2 + self%freq_bag2_Hz**2) - self%fusimotorActivation(2)) / self%tauBag2_s

    end function

    subroutine computeFiberTension(self, t, fascicleLength, fascicleVelocity, fascicleAcceleration)
        ! '''
        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t, fascicleLength, fascicleVelocity, fascicleAcceleration
        real(wp), dimension(6) :: dT

        dT = self%dTdt(t, fascicleLength, fascicleVelocity, fascicleAcceleration)
        self%fiberTension = self%fiberTension  + self%conf%timeStep_ms / 1000.0 * dT
    end subroutine

    function dTdt(self, t, fascicleLength, fascicleVelocity, fascicleAcceleration) result(dT)
        ! '''

        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t, fascicleLength, fascicleVelocity, fascicleAcceleration
        real(wp), dimension(6) :: dT
        real(wp) :: C
        
        if ((fascicleVelocity - self%fiberTension(2)/self%KsrBag1)  > 0.0) then
            C = 1.0
        else
            C = 0.42
        end if

        dT(1) = self%fiberTension(2)
        dT(2) = self%KsrBag1 / self%MBag1 * (C * self%betaBag1 * sign(1.0_wp, fascicleVelocity - &
                                            self%fiberTension(2)/self%KsrBag1) &
                                             * (abs(fascicleVelocity - self%fiberTension(2)/self%KsrBag1)**0.3) * &
                                             (fascicleLength - self%L0SrBag1 - self%fiberTension(1)/self%KsrBag1 - self%RBag1)&
                                             + self%KPrBag1 * (fascicleLength - self%L0SrBag1 - self%fiberTension(1)/self%KsrBag1&
                                              - self%L0PrBag1) &
                                             + self%MBag1 * fascicleAcceleration + self%GAMMABag1 - self%fiberTension(1))
        
        if ((fascicleVelocity - self%fiberTension(4)/self%KsrBag2)  > 0.0) then
            C = 1.0
        else
            C = 0.42
        end if
        
        dT(3) = self%fiberTension(4)
        dT(4) = self%KsrBag2 / self%MBag2 * (C * self%betaBag2 &
                                            * sign(1.0_wp,fascicleVelocity - self%fiberTension(4)/self%KsrBag2) &
                                        * (abs(fascicleVelocity - self%fiberTension(4)/self%KsrBag2)**0.3) * & 
                                        (fascicleLength - self%L0SrBag2 - self%fiberTension(3)/self%KsrBag2 - self%RBag2) &
                                        + self%KPrBag2 * (fascicleLength - self%L0SrBag2 - self%fiberTension(3)/self%KsrBag2 -&
                                         self%L0PrBag2) &
                                        + self%MBag2 * fascicleAcceleration + self%GAMMABag2 - self%fiberTension(3))
        
        if ((fascicleVelocity - self%fiberTension(6)/self%KsrChain)  > 0.0) then
            C = 1.0
        else
            C = 0.42
        end if

        dT(5) = self%fiberTension(6)
        dT(6) = self%KsrChain / self%MChain * (C * &
                                               self%betaChain * sign(1.0_wp, fascicleVelocity -&
                                                self%fiberTension(6)/self%KsrChain) &
                                        * (abs(fascicleVelocity - self%fiberTension(6)/self%KsrChain)**0.3) * &
                                        (fascicleLength - self%L0SrChain - self%fiberTension(5)/self%KsrChain - self%RChain) &
                                        + self%KPrChain * (fascicleLength - self%L0SrChain - self%fiberTension(5)/self%KsrChain -&
                                         self%L0PrChain) &
                                        + self%MChain * fascicleAcceleration + self%GAMMAChain - self%fiberTension(5))
    end function

    real(wp) function computeIa(self, t) result(IaFR)
        ! '''

        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp) :: smaller, larger

        call self%computePrimaryActivity(t)

        if (self%primaryPotentialBag1 >= (self%primaryPotentialBag2 + self%primaryPotentialChain)) then
            larger = self%primaryPotentialBag1
            smaller = self%primaryPotentialBag2 + self%primaryPotentialChain
        else
            smaller = self%primaryPotentialBag1
            larger = self%primaryPotentialBag2 + self%primaryPotentialChain
        end if
        IaFR = self%SOcclusionFactor * smaller + larger
    end function

    real(wp) function computeII(self, t, fascicleLength) result(IIFR)
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in) :: fascicleLength

        call self%computeSecondaryActivity(t, fascicleLength)

        IIFR = self%secondaryPotentialBag2 + self%secondaryPotentialChain
    end function

    subroutine computePrimaryActivity(self, t) 
        ! '''

        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t
        
        

        self%primaryPotentialBag1 = self%GPrimaryBag1 * (self%fiberTension(1) / self%KsrBag1 - &
                                                         self%LNSrBag1 + self%L0SrBag1)

        self%primaryPotentialBag2 = self%GPrimaryBag2 * (self%fiberTension(3) / self%KsrBag2 - &
                                                         self%LNSrBag2 + self%L0SrBag2)
                                                    
        self%primaryPotentialChain = self%GPrimaryChain * (self%fiberTension(5) / self%KsrChain - &
                                                           self%LNSrChain + self%L0SrChain)
    end subroutine

    subroutine computeSecondaryActivity(self, t, fascicleLength) 
        ! '''

        ! '''
        class(MuscleSpindle), intent(inout) :: self
        real(wp), intent(in) :: t
        real(wp), intent(in) :: fascicleLength

        self%secondaryPotentialBag2 = self%GSecondaryBag2 * ((self%XBag2 * self%LSecondaryBag2 / self%L0SrBag2) * &
                                      (self%fiberTension(3)/self%KsrBag2 - (self%LNSrBag2 - self%L0SrBag2)) + &
                                      (1 - self%XBag2) * self%LSecondaryBag2 / self%L0PrBag2 * &
                                      (fascicleLength - self%fiberTension(3)/self%KsrBag2 - self%L0SrBag2 - self%LNPrBag2)) 
                                                    
        self%secondaryPotentialChain = self%GSecondaryChain * ((self%XChain * self%LSecondaryChain / self%L0SrChain) * &
                                      (self%fiberTension(5)/self%KsrChain - (self%LNSrChain - self%L0SrChain)) + &
                                      (1 - self%XChain) * self%LSecondaryChain / self%L0PrChain * &
                                      (fascicleLength - self%fiberTension(5)/self%KsrChain - self%L0SrChain - self%LNPrChain)) 
    end subroutine


    subroutine reset(self)
        ! '''
        ! '''
        class(MuscleSpindle), intent(inout) :: self

        self%betaBag1 = 0.0
        self%betaBag2 = 0.0
        self%betaChain = 0.0

        self%GAMMABag1 = 0.0
        self%GAMMABag2 = 0.0
        self%GAMMAChain = 0.0

        self%primaryPotentialBag1 = 0.0
        self%primaryPotentialBag2 = 0.0
        self%primaryPotentialChain = 0.0

        self%secondaryPotentialBag1 = 0.0
        self%secondaryPotentialBag2 = 0.0
        self%secondaryPotentialChain = 0.0
      
        self%fusimotorActivation(:) = 0.0

        self%fiberTension(:) = 0.0

        self%IaFR_Hz = 0.0
        self%IIFR_Hz = 0.0
    end subroutine



end module MuscleSpindleClass       
