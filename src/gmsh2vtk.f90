!--------------------------------------------------------------------
!
! All the subroutines that handle GMSH Mesh Files.
! Compatible with MSH file format version 4.1.
! https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
!
!--------------------------------------------------------------------

      include "./gmshMod.f90"

      subroutine conv_gmsh2vtk(msh)
      use varmod
      use gmshMod
      implicit none

      type(meshType), intent(out)  :: msh

      real(kind=8) :: rtemp
      integer :: fid
      integer :: i, j, n, itemp, tag, iFa
      character(len=strL) :: rLine

      fid = 100
      write(stdout,ftab1) "Loading file "//TRIM(msh%fname)
      open(fid, file=TRIM(msh%fname))

      ! dimension
      write(stdout,ftab1,advance='no') &
      "Dimension of the body mesh (2/3):  "
      read(*,*) nsd

      ! metadata
      call gmsh_readmeta(msh,fid)

      ! PhysicalNames
      call gmsh_readphysicalnames(fid)


      ! ! Entities
      ! call findKwrd(fid, "$Entities")
      ! read(fid,*) numPoints, numCurves, numSurfaces, numVolumes

      ! print *, "numPoints   =", numPoints    
      ! print *, "numCurves   =", numCurves    
      ! print *, "numSurfaces =", numSurfaces  
      ! print *, "numVolumes  =", numVolumes   
      
      ! Nodes
      call gmsh_readnodes(fid)

      ! Elements
      call gmsh_readelements(fid)




      close(fid)

      ! debug
      stop

      return
      end subroutine conv_gmsh2vtk

!***********************************************************************
!     Read gmsh meta data
      subroutine gmsh_readmeta(msh,fid)
      use varmod
      use gmshMod

      implicit none

      type(meshType), intent(out)  :: msh
      integer, intent(inout) :: fid

      real(kind=8) :: rtemp
      integer :: i, j

      ! metadata
      call findKwrd(fid, "$MeshFormat")
      read(fid,*) rtemp, i, j
      if ( i .ne. 0) then
         write(stdout,ftab4) "ERROR: "//TRIM(msh%fname)// &
         " is not stored in plain text."
         stop
      end if

      return
      end subroutine gmsh_readmeta
!***********************************************************************
!     Read PhysicalNames
      subroutine gmsh_readphysicalnames(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, n

      call findKwrd(fid, "$PhysicalNames")
      read(fid,*) n

      gmshPhysicalNames%num = n
      allocate(gmshPhysicalNames%dimension(n))
      allocate(gmshPhysicalNames%physicalTag(n))
      allocate(gmshPhysicalNames%name(n))
      do i = 1, gmshPhysicalNames%num
         read(fid,*) gmshPhysicalNames%dimension(i), gmshPhysicalNames%physicalTag(i), gmshPhysicalNames%name(i)
      end do

      return
      end subroutine gmsh_readphysicalnames
!***********************************************************************
!     Read Nodes
      subroutine gmsh_readnodes(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, j, itemp
      integer, allocatable :: NodeTag(:)
   
      call findKwrd(fid, "$Nodes")
      read(fid,*) gmshNodes%numNodeBlocks, gmshNodes%numNodes, &
                  gmshNodes%minNodeTag, gmshNodes%maxNodeTag

      allocate(gmshNodes%numNodesInBlock(gmshNodes%numNodeBlocks))
      allocate(NodeTag(gmshNodes%numNodes))
      allocate(gmshNodes%coord(3,gmshNodes%numNodes))
      
      ! read coordiantes block by block
      do i = 1, gmshNodes%numNodeBlocks
         read(fid,*) itemp, itemp, itemp, gmshNodes%numNodesInBlock(i)

         do j = 1, gmshNodes%numNodesInBlock(i)
            read(fid,*) NodeTag(j)
         end do
         do j = 1, gmshNodes%numNodesInBlock(i)
            read(fid,*) gmshNodes%coord(:,NodeTag(j))
         end do
      end do
   
      return
      end subroutine gmsh_readnodes
!***********************************************************************
!     Read Elements
      subroutine gmsh_readelements(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, j, itemp, ie, eNoN, maxeNoN

      call findKwrd(fid, "$Elements")
      read(fid,*) gmshElements%numElementBlocks, gmshElements%numElements, &
                  gmshElements%minElementTag, gmshElements%maxElementTag

      allocate(gmshElements%ElementDim(gmshElements%numElementBlocks))
      allocate(gmshElements%EntityTag(gmshElements%numElementBlocks))
      allocate(gmshElements%eNoN(gmshElements%numElementBlocks))
      allocate(gmshElements%numElementsInBlock(gmshElements%numElementBlocks))

      ! find maximum eNoN
      do i = 1, gmshElements%numElementBlocks
         read(fid,*) gmshElements%ElementDim(i), gmshElements%EntityTag(i), &
                     itemp, gmshElements%numElementsInBlock(i)
         select case (itemp)
         case(1)
            gmshElements%eNoN(i) = 2 ! 2-node line.
         case(2)
            gmshElements%eNoN(i) = 3 ! 3-node triangle.
         case(3)
            gmshElements%eNoN(i) = 4 ! 4-node quadrangle.
         case(4)
            gmshElements%eNoN(i) = 4 ! 4-node tetrahedron.
         case(5)
            gmshElements%eNoN(i) = 8 ! 8-node hexahedron.
         case(8)
            gmshElements%eNoN(i) = 3 ! 3-node second order line
         case(9)
            gmshElements%eNoN(i) = 6 ! 6-node second order triangle
         case(10)
            gmshElements%eNoN(i) = 9 ! 9-node second order quadrangle 
         case(11)
            gmshElements%eNoN(i) = 10 ! 10-node second order tetrahedron 
         case default 
            write(stdout,ftab4) "ERROR: elelment type not defined!"
         end select

         do j = 1, gmshElements%numElementsInBlock(i)
            read(fid,*)
         end do
      end do

      maxeNoN = MAXVAL(gmshElements%eNoN)
      allocate(gmshElements%conn(maxeNoN,gmshElements%numElements))
      ie   = 0
      gmshElements%conn = 0
      rewind(fid)
      call findKwrd(fid, "$Elements")
      read(fid,*)
      do i = 1, gmshElements%numElementBlocks
         read(fid,*) 
         eNoN = gmshElements%eNoN(i)
         do j = 1, gmshElements%numElementsInBlock(i)
            ie = ie+1
            read(fid,*) itemp, gmshElements%conn(1:eNoN,ie)
         end do
      end do

      return
      end subroutine gmsh_readelements