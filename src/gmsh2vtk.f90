!--------------------------------------------------------------------
!
! All the subroutines that handle GMSH Mesh Files 
!
!--------------------------------------------------------------------

      subroutine conv_gmsh2vtk(msh)
      use varmod
      implicit none

      type(meshType), intent(out)  :: msh

      return
      end subroutine conv_gmsh2vtk