! '''
!     Neuromuscular simulator in Python.
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

module CharacterMatrixClass
    use CharacterArrayClass
    implicit none
    private
    public :: CharacterMatrix

    type CharacterMatrix
        type(CharacterArray), dimension(:), allocatable :: item

        contains
            procedure :: append
    end type CharacterMatrix
    
    interface CharacterMatrix
        module procedure :: init_CharacterMatrix
    end interface

    contains

        type(CharacterMatrix) function init_CharacterMatrix()
            if (allocated(init_CharacterMatrix%item)) deallocate(init_CharacterMatrix%item)
        end function

        subroutine append(self, newLine)
            class(CharacterMatrix), intent(inout) :: self
            class(CharacterArray), intent(in) :: newLine
            type(CharacterMatrix) :: clist
            integer :: isize, i

            clist = CharacterMatrix()
            
            
            if(allocated(self%item)) then
                isize = size(self%item)
                allocate(clist%item(isize+1))
                do i=1,isize          
                    clist%item(i) = self%item(i)
                end do
                clist%item(isize+1) = newLine

                deallocate(self%item)
                allocate(self%item(isize+1))
                
                do i=1, isize + 1         
                    self%item(i) = clist%item(i)
                end do

            else
                allocate(self%item(1))
                self%item(1) = newLine
            end if
        end subroutine

end module CharacterMatrixClass