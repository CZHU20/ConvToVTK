!--------------------------------------------------------------------
!
! All the subroutines that handle GMSH Mesh Files.
! Compatible with MSH file format version 4.1.
! https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
!
!--------------------------------------------------------------------

      ! include "./gmshMod.f90"

      subroutine conv_gmsh2vtk(msh)
      use varmod
      implicit none

      type(meshType), intent(out)  :: msh

      real(kind=8) :: rtemp
      integer :: fid
      integer :: i, j, n, itemp, tag, iFa
      character(len=strL) :: rLine

      ! Physical names
      integer, allocatable :: physicalnames_tag(:)

      ! Nodes
      integer :: numNodeBlocks, numNodesInBlock
      integer,allocatable :: NodeTag(:)

      ! Elements
      integer :: numElementBlocks, numElementsInBlocks

      ! ! Entities
      ! integer :: numPoints, numCurves, numSurfaces, numVolumes

      fid = 100
      write(stdout,ftab1) "Loading file "//TRIM(msh%fname)
      open(fid, file=TRIM(msh%fname))

      ! dimension
      write(stdout,ftab1,advance='no') &
      "Dimension of the body mesh (2/3):  "
      read(*,*) nsd

      ! metadata
      call findKwrd(fid, "$MeshFormat")
      read(fid,*) rtemp, i, j
      if ( i .ne. 0) then
         write(stdout,ftab4) "ERROR: "//TRIM(msh%fname)// &
         " is not stored in plain text."
         stop
      end if

      ! PhysicalNames
      call findKwrd(fid, "$PhysicalNames")
      read(fid,*) n
      j = 0
      do i = 1, n
         read(fid, *) itemp
         if ( itemp .eq. (nsd-1) ) j = j+1 ! boundary
      end do
      msh%nFa = j
      allocate(physicalnames_tag(msh%nFa), msh%fa(msh%nFa))
      rewind(fid)
      call findKwrd(fid, "$PhysicalNames")
      read(fid,*) n
      j = 0
      do i = 1, n
         read(fid, *) itemp, tag, rLine
         if ( itemp .eq. (nsd-1) ) then
            j = j+1 ! boundary
            physicalnames_tag(j) = tag
            msh%fa(j)%fname = TRIM(rLine)
         end if
      end do

      ! ! Entities
      ! call findKwrd(fid, "$Entities")
      ! read(fid,*) numPoints, numCurves, numSurfaces, numVolumes

      ! print *, "numPoints   =", numPoints    
      ! print *, "numCurves   =", numCurves    
      ! print *, "numSurfaces =", numSurfaces  
      ! print *, "numVolumes  =", numVolumes   
      
      ! Nodes
      call findKwrd(fid, "$Nodes")
      read(fid,*) numNodeBlocks, msh%nNo
      allocate(NodeTag(msh%nNo), msh%x(3,msh%nNo))
      do i = 1, numNodeBlocks
         read(fid,*) itemp, itemp, itemp, numNodesInBlock

         do j = 1, numNodesInBlock
            read(fid,*) NodeTag(j)
         end do
         do j = 1, numNodesInBlock
            read(fid,*) msh%x(:,NodeTag(j))
         end do
      end do

      ! Elements
      call findKwrd(fid, "$Elements")
      read(fid,*) numElementBlocks, msh%nEl, itemp, itemp
      


      close(fid)

      ! debug
      stop

      return
      end subroutine conv_gmsh2vtk