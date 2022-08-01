! # '''
! #     Neuromuscular simulator in Fortran.
! #     Copyright (C) 2018  Renato Naville Watanabe

! #     This program is free software: you can redistribute it and/or modify
! #     it under the terms of the GNU General Public License as published by
! #     the Free Software Foundation, either version 3 of the License, or
! #     any later version.

! #     This program is distributed in the hope that it will be useful,
! #     but WITHOUT ANY WARRANTY; without even the implied warranty of
! #     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! #     GNU General Public License for more details.

! #     You should have received a copy of the GNU General Public License
! #     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! #     Contact: renato.watanabe@ufabc.edu.br
! # '''

module SynapsePointerClass

    use SynapseClass
    implicit none
    private
    public :: SynapsePointer

    type SynapsePointer
        type(Synapse), pointer :: synapse

        contains
            procedure :: assignSynapse

    end type SynapsePointer

    interface SynapsePointer
        module procedure init_SynapsePointer
    end interface SynapsePointer

    contains

        type(SynapsePointer) function init_SynapsePointer()
            if (associated(init_SynapsePointer%synapse)) nullify(init_SynapsePointer%synapse) 
        end function

        subroutine assignSynapse(self, newSynapse)
            class(SynapsePointer), intent(inout) :: self
            class(Synapse), intent(in), target :: newSynapse

            self%synapse => newSynapse
        end subroutine

end module SynapsePointerClass