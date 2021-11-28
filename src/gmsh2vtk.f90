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

      ! Entities
      call gmsh_readentities(fid)
      
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
      rewind(fid)
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
      
      rewind(fid)
      call findKwrd(fid, "$PhysicalNames")
      read(fid,*) n

      gmshPhysicalNames%num = n
      allocate(gmshPhysicalNames%dimension(n))
      allocate(gmshPhysicalNames%physicalTag(n))
      allocate(gmshPhysicalNames%name(n))
      do i = 1, gmshPhysicalNames%num
         read(fid,*) gmshPhysicalNames%dimension(i), gmshPhysicalNames%physicalTag(i), gmshPhysicalNames%name(i)
      end do

      if (nsd .ne. MAXVAL(gmshPhysicalNames%dimension)) then
         write(stdout,"(14X,A,I5)") "nsd = ", nsd
         write(stdout,"(14X,A,I5)") "MAXVAL(gmshPhysicalNames%dimension) = ", MAXVAL(gmshPhysicalNames%dimension)
         write(stdout,ftab4) "Both body mesh and surface mesh need to be assigned to a physical group."
         stop
      end if

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
         
      rewind(fid)
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

      rewind(fid)
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
         ! case(5)
         !    gmshElements%eNoN(i) = 8 ! 8-node hexahedron.
         ! case(8)
         !    gmshElements%eNoN(i) = 3 ! 3-node second order line
         ! case(9)
         !    gmshElements%eNoN(i) = 6 ! 6-node second order triangle
         ! case(10)
         !    gmshElements%eNoN(i) = 9 ! 9-node second order quadrangle 
         ! case(11)
         !    gmshElements%eNoN(i) = 10 ! 10-node second order tetrahedron 
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
!***********************************************************************
!     Read Entities
      subroutine gmsh_readentities(fid)
      use varmod
      use gmshMod

      implicit none

      integer, intent(inout) :: fid

      integer :: i, maxNumTags
      integer :: numPhysicalTags, PhysicalTags, Tag
      real(kind=8) :: rtemp

      rewind(fid)
      call findKwrd(fid, "$Entities")
      read(fid,*) gmshEntities%numPoints, gmshEntities%numCurves, &
                  gmshEntities%numSurfaces, gmshEntities%numVolumes

      maxNumTags = 1 !gmshPhysicalNames%num

      ! Point entites
      if (gmshEntities%numPoints > 0) then
         allocate(gmshEntities%Point_numPhysicalTags(gmshEntities%numPoints))
         allocate(gmshEntities%Point_PhysicalTags(gmshEntities%numPoints,maxNumTags))
         gmshEntities%Point_numPhysicalTags = 0
         gmshEntities%Point_PhysicalTags = 0

         do i = 1, gmshEntities%numPoints
            read(fid,*) Tag, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each point entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Point_numPhysicalTags(Tag) = 1
               gmshEntities%Point_PhysicalTags(Tag,1) = PhysicalTags
            end if
         end do
      end if

      ! Curve entites
      if (gmshEntities%numCurves > 0) then
         allocate(gmshEntities%Curve_numPhysicalTags(gmshEntities%numCurves))
         allocate(gmshEntities%Curve_PhysicalTags(gmshEntities%numCurves,maxNumTags))
         gmshEntities%Curve_numPhysicalTags = 0
         gmshEntities%Curve_PhysicalTags = 0

         do i = 1, gmshEntities%numCurves
            read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each curve entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Curve_numPhysicalTags(Tag) = 1
               gmshEntities%Curve_PhysicalTags(Tag,1) = PhysicalTags
            end if
         end do
      end if

      ! Surface entites
      if (gmshEntities%numSurfaces > 0) then
         allocate(gmshEntities%Surface_numPhysicalTags(gmshEntities%numSurfaces))
         allocate(gmshEntities%Surface_PhysicalTags(gmshEntities%numSurfaces,maxNumTags))
         gmshEntities%Surface_numPhysicalTags = 0
         gmshEntities%Surface_PhysicalTags = 0

         do i = 1, gmshEntities%numSurfaces
            read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each surface entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Surface_numPhysicalTags(Tag) = 1
               gmshEntities%Surface_PhysicalTags(Tag,1) = PhysicalTags
            end if
         end do
      end if

      ! Volume entites
      if (gmshEntities%numVolumes > 0) then
         allocate(gmshEntities%Volume_numPhysicalTags(gmshEntities%numVolumes))
         allocate(gmshEntities%Volume_PhysicalTags(gmshEntities%numVolumes,maxNumTags))
         gmshEntities%Volume_numPhysicalTags = 0
         gmshEntities%Volume_PhysicalTags = 0

         do i = 1, gmshEntities%numVolumes
            read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags
            if (numPhysicalTags > 1) then
               write(stdout,ftab4) "ERROR: each volume entity can only has one physical name."
               stop
            elseif (numPhysicalTags == 1) then
               BACKSPACE(fid)
               read(fid,*) Tag, rtemp, rtemp, rtemp, rtemp, rtemp, rtemp, numPhysicalTags, PhysicalTags
               gmshEntities%Volume_numPhysicalTags(Tag) = 1
               gmshEntities%Volume_PhysicalTags(Tag,1) = PhysicalTags
            end if
         end do
      end if

      ! debug

      ! print *, "numPoints   =", gmshEntities%numPoints    
      ! print *, "numCurves   =", gmshEntities%numCurves    
      ! print *, "numSurfaces =", gmshEntities%numSurfaces  
      ! print *, "numVolumes  =", gmshEntities%numVolumes  

      ! print *, repeat('---',10)
      ! print *, gmshEntities%Point_numPhysicalTags
      ! print *, gmshEntities%Point_PhysicalTags
      ! print *, repeat('---',10)
      ! print *, gmshEntities%Curve_numPhysicalTags
      ! print *, gmshEntities%Curve_PhysicalTags
      ! print *, repeat('---',10)
      ! print *, gmshEntities%Surface_numPhysicalTags
      ! print *, gmshEntities%Surface_PhysicalTags
      ! print *, repeat('---',10)
      ! print *, gmshEntities%Volume_numPhysicalTags
      ! print *, gmshEntities%Volume_PhysicalTags

      return
      end subroutine gmsh_readentities