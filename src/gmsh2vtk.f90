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

      integer :: fid, i

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

      ! Assemble body mesh
      call gmsh_bodymesh(msh)

      ! Assemble boundary mesh
      call gmsh_boundarymesh(msh)

      call selectel(msh)

      ! debug
      do i = 1, msh%nNo
         write(11111,*) i, msh%x(1:nsd,i)
      end do

      do i = 1, msh%nEl
         write(22222,"(5I8)") i, msh%IEN(:,i)
      end do

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
         case(15)
            gmshElements%eNoN(i) = 1 ! 1-node point.
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
!***********************************************************************
!     Assemble body mesh
      subroutine gmsh_bodymesh(msh)
      use varmod
      use gmshMod

      implicit none

      Type(meshType), intent(inout) :: msh

      integer :: i, j, id, physicalTag, Tag
      character(len=20) :: name

      j = 0
      do i = 1, gmshPhysicalNames%num
         if (nsd .eq. gmshPhysicalNames%dimension(i)) then
            id = i
            physicalTag = gmshPhysicalNames%physicalTag(i)
            j = j+1
            name = gmshPhysicalNames%name(i)
         end if 
      end do
      if (j .gt. 1) then
         write(stdout,ftab4) "ERROR: can only handle one body mesh."
         stop
      end if

      ! determine the surface tag for 2D problem
      ! volume tag for 3D problem
      if (nsd .eq. 2) then
         do i = 1, gmshEntities%numSurfaces
            do j = 1, gmshEntities%Surface_numPhysicalTags(i)
               if (physicalTag .eq. gmshEntities%Surface_PhysicalTags(i,j)) &
                  Tag = i
            end do
         end do
      elseif (nsd .eq. 3) then
         do i = 1, gmshEntities%numVolumes
            do j = 1, gmshEntities%Volume_numPhysicalTags(i)
               if (physicalTag .eq. gmshEntities%Volume_PhysicalTags(i,j)) &
                  Tag = i
            end do
         end do
      end if

      ! connectivity
      j = 1
      do i = 1, gmshElements%numElementBlocks
         if (gmshElements%ElementDim(i) .eq. nsd) then
            if (gmshElements%EntityTag(i) .eq. Tag) then
               msh%nEl = gmshElements%numElementsinBlock(i)
               msh%eNoN = gmshElements%eNoN(i)
               allocate(msh%IEN(msh%eNoN,msh%nEl))
               msh%IEN = gmshElements%conn(1:msh%eNoN,j:(j+msh%nEl))
            end if
         end if
         j = j + gmshElements%numElementsInBlock(i)
      end do

      ! coordinates
      msh%nNo = gmshNodes%numNodes
      allocate(msh%x(nsd,msh%nNo))
      msh%x = gmshNodes%coord

      return
      end subroutine gmsh_bodymesh
!***********************************************************************
!     Assemble boundary meshes
      subroutine gmsh_boundarymesh(msh)
      use varmod
      use gmshMod

      implicit none

      Type(meshType), intent(inout) :: msh

      integer :: i, j, k, l, ii, ie, Tag, iFa
      integer,allocatable :: physicalTag(:), gIEN(:,:), g2l(:), l2g(:)
      integer,allocatable :: sharedelem(:,:), numshared(:)
      logical,allocatable :: flag(:)

      allocate(g2l(msh%nNo))

      ! element that share the same node
      ! this is for finding global element id
      allocate(sharedelem(msh%nNo,50), numshared(msh%nNo))
      sharedelem = 0
      numshared  = 0
      do i = 1, msh%nEl
         do j = 1, msh%eNoN
            k = msh%IEN(j,i)
            numshared(k) = numshared(k)+1
            if (numshared(k) .GT. 50) then
               write(stdout,ftab4) "ERROR: numshared is larger than 50."
               stop
            end if
            sharedelem(k,numshared(k)) = i
         end do
      end do
      
      msh%nFa = 0
      do i = 1, gmshPhysicalNames%num
         if (gmshPhysicalNames%dimension(i) .eq. nsd-1) then
            msh%nFa = msh%nFa+1
         end if
      end do
      allocate(msh%fa(msh%nFa), physicalTag(msh%nFa))
      j = 0
      do i = 1, gmshPhysicalNames%num
         if (gmshPhysicalNames%dimension(i) .eq. nsd-1) then
            j = j + 1
            physicalTag(j)  = gmshPhysicalNames%physicalTag(i)
            msh%fa(j)%fname = gmshPhysicalNames%name(i) 
         end if
      end do

      ! loop over faces to construct mesh files
      do iFa = 1, msh%nFa
         ! Find the corresponding EntityTags for boundary meshes
         ! ndim=2, CurveTag
         ! ndim=3, SurfaceTag
         if (nsd .eq. 2) then
            do j = 1, gmshEntities%numCurves
               do k = 1, gmshEntities%Curve_numPhysicalTags(j)
                  if (physicalTag(iFa) .eq. gmshEntities%Curve_physicalTags(j,k)) &
                     Tag = j
               end do
            end do
         elseif (nsd .eq. 3) then
            do j = 1, gmshEntities%numSurfaces
               do k = 1, gmshEntities%Surface_numPhysicalTags(j)
                  if (physicalTag(iFa) .eq. gmshEntities%Surface_physicalTags(j,k)) &
                     Tag = j
               end do
            end do
         end if

         ! coord and conn for boundary mesh
         j = 1
         do i = 1, gmshElements%numElementBlocks
            if (gmshElements%EntityTag(i) .eq. Tag) then
               msh%fa(iFa)%eNoN = gmshElements%eNoN(i)
               msh%fa(iFa)%nEl  = gmshElements%numElementsInBlock(i)
               if (allocated(gIEN)) deallocate(gIEN)
               allocate(gIEN(msh%fa(iFa)%eNoN,msh%fa(iFa)%nEl))
               gIEN = gmshElements%conn(1:msh%fa(iFa)%eNoN,j:(j+msh%fa(iFa)%nEl))
               exit
            end if
            j = j + gmshElements%numElementsInBlock(i)
         end do

         ! convert gIEN to local IEN and determine gN
         g2l = 0
         msh%fa(iFa)%nNo = 0
         do i = 1, msh%fa(iFa)%nEl
            do j = 1, msh%fa(iFa)%eNoN
               k = gIEN(j,i)
               if (g2l(k) .eq. 0) then
                  msh%fa(iFa)%nNo = msh%fa(iFa)%nNo + 1
                  g2l(k) = msh%fa(iFa)%nNo
               end if
            end do
         end do
         allocate(msh%fa(iFa)%gN(msh%fa(iFa)%nNo),msh%fa(iFa)%x(3,msh%fa(iFa)%nNo))
         do i = 1, msh%nNo
            if (g2l(i) .ne. 0) then
               j = g2l(i)
               msh%fa(iFa)%gN(j) = i
               msh%fa(iFa)%x(:,j) = msh%x(:,i)
            end if
         end do
         allocate(msh%fa(iFa)%IEN(msh%fa(iFa)%eNoN,msh%fa(iFa)%nEl))
         do i = 1, msh%fa(iFa)%nEl
            do j = 1, msh%fa(iFa)%eNoN
               msh%fa(iFa)%IEN(j,i) = g2l(gIEN(j,i))
            end do
         end do

         ! find global element id
         allocate(msh%fa(iFa)%gE(msh%fa(iFa)%nEl))
         if (allocated(flag)) deallocate(flag)
         allocate(flag(msh%fa(iFa)%eNoN))
         do i = 1, msh%fa(iFa)%nEl
            ii = gIEN(1,i)
            do k = 1, numshared(ii)
               flag = .FALSE.
               ie = sharedelem(ii,k)
               do j = 1, msh%fa(iFa)%eNoN               
                  do l = 1, msh%eNoN
                     if (gIEN(j,i) .eq. msh%IEN(l,ie)) then
                        flag(j) = .TRUE.
                        exit
                     end if
                  end do
               end do
               if (ALL(flag)) then
                  msh%fa(iFa)%gE(i) = ie
                  exit
               end if 
            end do
         end do
      end do

      return
      end subroutine gmsh_boundarymesh
!***********************************************************************