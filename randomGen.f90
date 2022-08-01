module randomGen
    implicit none
    integer, parameter :: wp = kind(1.0d0)
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    

    contains 

        real(wp) function randn() result(x)
            ! Generates a random number from a normal distribution
            ! using Box-Müller method
            
            real(wp) :: u1, u2

            call random_number(u1)
            call random_number(u2)

            x = sqrt(-2.0*log(u1))*sin(2.0*pi*u2)
            
        end function randn        
        
end module randomGen