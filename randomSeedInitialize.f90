module randomSeedInitialize
    implicit none
    contains 

        subroutine init_random_seed()

            integer :: i, n, clock
            integer, dimension(:), allocatable :: seed

            call random_seed(size = n)
            
            allocate(seed(n))

            call system_clock(count=clock)

            
            seed = clock + 37 * [(i - 1, i = 1, n)]
            call random_seed(put = seed)
            
            deallocate(seed)
        end subroutine init_random_seed

        integer function init_random_seed_returning() result(clock)

            integer :: i, n
            integer, dimension(:), allocatable :: seed

            call random_seed(size = n)
            
            allocate(seed(n))

            call system_clock(count=clock)

            
            seed = clock + 37 * [(i - 1, i = 1, n)]
            call random_seed(put = seed)
            
            deallocate(seed)
        end function init_random_seed_returning

        subroutine init_seed(clock)
            integer, intent(in) ::  clock
            integer :: i, n
            integer, dimension(:), allocatable :: seed

            call random_seed(size = n)
            
            allocate(seed(n))

            seed = clock + 37 * [(i - 1, i = 1, n)]
            call random_seed(put = seed)
            
            deallocate(seed)
        end subroutine init_seed

        
end module randomSeedInitialize