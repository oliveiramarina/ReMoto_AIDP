  module DynamicalArrays
    implicit none
	integer, parameter :: wp1 = kind(1.0d0)

	type fila
		real(wp1), allocatable :: prioritytime(:)
		integer                 :: n = 0

		contains
		procedure :: top
		procedure :: pop
		procedure :: enqueue
		procedure :: shiftup
		
	end type
	
	contains
	
		subroutine shiftup(this, a)
			class (fila)           :: this
			integer                 :: a, parent, child
			  associate (x => this%prioritytime)
			  parent = a
			  do while(parent*2 <= this%n)
				child = parent*2
				if (child + 1 <= this%n) then
				  if (x(child+1) < x(child)) then
					child = child +1
				  end if
				end if
				if (x(parent) > x(child)) then
				  x([child, parent]) = x([parent, child])
				  parent = child
				else
				  exit
				end if
			  end do
			  end associate
		end subroutine
		
		
		function top(this) result (res)
		  class(fila) :: this
		  real(wp1)   :: res
		  if (this%n >= 1) then
			res = this%prioritytime(1)
		  else
			res = -1e6
		  end if
		end function

		subroutine pop(this)
		  class(fila) :: this
		  this%prioritytime(1) = this%prioritytime(this%n)
		  this%n = this%n - 1
		  call this%shiftup(1)
		end subroutine

		subroutine enqueue(this, time)
		  class(fila), intent(inout) :: this
		  real(wp1)                    :: time
		  real(wp1), allocatable     :: tmp(:)
		  integer                     :: i
		  this%n = this%n +1
		  if (.not.allocated(this%prioritytime)) allocate(this%prioritytime(1))
		  if (size(this%prioritytime)<this%n) then
			allocate(tmp(2*size(this%prioritytime)))
			tmp(1:this%n-1) = this%prioritytime
			call move_alloc(tmp, this%prioritytime)
		  end if
		  this%prioritytime(this%n) = time
		  i = this%n
		  do
			i = i / 2
			if (i==0) exit
			call this%shiftup(i)
		  end do
		end subroutine

      subroutine AddToList(list, element)

          implicit none

          integer :: i, isize
          double precision, intent(in) :: element
          double precision, dimension(:), allocatable, intent(inout) :: list
          double precision, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize
              clist(i) = list(i)
              end do
              clist(isize+1) = element

              deallocate(list)

              allocate(list(isize+1))

              do i=1,isize + 1
                list(i) = clist(i)
              end do
              deallocate(clist)

          else
              allocate(list(1))
              list(1) = element
          end if


      end subroutine AddToList

      subroutine integerAddToList(list, element)

          implicit none

          integer :: i, isize
          integer, intent(in) :: element
          integer, dimension(:), allocatable, intent(inout) :: list
          integer, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize
                clist(i) = list(i)
              end do

              clist(isize+1) = element

              deallocate(list)
              allocate(list(isize+1))

              do i=1,isize + 1
                list(i) = clist(i)
              end do
              deallocate(clist)
          else
              allocate(list(1))
              list(1) = element
          end if


      end subroutine integerAddToList

      subroutine boolAddToList(list, element)

          implicit none

          integer :: i, isize
          logical, intent(in) :: element
          logical, dimension(:), allocatable, intent(inout) :: list
          logical, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize
                clist(i) = list(i)
              end do

              clist(isize+1) = element

              deallocate(list)
              allocate(list(isize+1))

              do i=1,isize + 1
                list(i) = clist(i)
              end do
              deallocate(clist)
          else
              allocate(list(1))
              list(1) = element
          end if


      end subroutine boolAddToList


	

  end module DynamicalArrays