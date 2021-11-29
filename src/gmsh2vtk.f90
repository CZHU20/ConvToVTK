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

      integer :: fid

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
      ! Please note not all nodes are grid points, e.g. center of the circle.
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
      do fid = 1, msh%nNo
         write(11111,*) fid, msh%x(1:nsd,fid)
      end do

      do fid = 1, msh%nEl
         write(22222,"(5I8)") fid, msh%IEN(:,fid)
      end do

      ! stop

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

      return
      end subroutine gmsh_readentities
!***********************************************************************
!     Assemble body mesh
      subroutine gmsh_bodymesh(msh)
      use varmod
      use gmshMod

      implicit none

      Type(meshType), intent(inout) :: msh

      integer :: i, j, id, physicalTag, numEntityTags
      integer :: Tag, iTag, nEl
      integer,allocatable :: entityTag(:), gIEN(:,:)
      character(len=80) :: name

      j = 0
      physicalTag = -1
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
      if (physicalTag .eq. -1) then
         write(stdout,ftab4) "ERROR: cannot find corresponding entitiy."
         stop
      end if

      ! Using the physical tag to determine the 
      ! surface tag (entity tag) for 2D problem or 
      ! volume tag (entity tag) for 3D problem.
      ! Please note that there might be multiple entity tags
      ! for one physical tag.
      if (nsd .eq. 2) then
         allocate(entityTag(gmshEntities%numSurfaces))
         entityTag = 0
         numEntityTags = 0
         do i = 1, gmshEntities%numSurfaces
            do j = 1, gmshEntities%Surface_numPhysicalTags(i)
               if (physicalTag .eq. gmshEntities%Surface_PhysicalTags(i,j)) then
                  numEntityTags = numEntityTags + 1
                  entityTag(numEntityTags) = i
               end if
            end do
         end do
      elseif (nsd .eq. 3) then
         allocate(entityTag(gmshEntities%numVolumes))
         entityTag = 0
         numEntityTags = 0
         do i = 1, gmshEntities%numVolumes
            do j = 1, gmshEntities%Volume_numPhysicalTags(i)
               if (physicalTag .eq. gmshEntities%Volume_PhysicalTags(i,j)) then
                  numEntityTags = numEntityTags + 1
                  entityTag(numEntityTags) = i
               end if
            end do
         end do
      end if
      if (numEntityTags .eq. 0) then
         write(stdout,ftab4) "ERROR: cannot find corresponding entitiy."
         stop
      end if

      ! raw global connectivity
      nEl = 0
      do iTag = 1, numEntityTags
         Tag = entityTag(iTag)
         j = 0
         do i = 1, gmshElements%numElementBlocks
            if (gmshElements%ElementDim(i) .eq. nsd) then
               if (gmshElements%EntityTag(i) .eq. Tag) then
                  msh%eNoN = gmshElements%eNoN(i)
                  if (.not. allocated(gIEN)) & 
                     allocate(gIEN(msh%eNoN,gmshElements%numElements))
                  gIEN(1:msh%eNoN,nEl+1:nEl + gmshElements%numElementsinBlock(i)) = &
                  gmshElements%conn(1:msh%eNoN,j+1:j+gmshElements%numElementsinBlock(i))
                  nEl = nEl + gmshElements%numElementsinBlock(i)
                  exit
               end if
            end if
            j = j + gmshElements%numElementsInBlock(i)
         end do
      end do

      ! Because not all points are mesh points, e.g. center of the circule,
      ! we need to build global-to-local mapping to eliminate those points.
      allocate(g2l(gmshNodes%numNodes))
      g2l = 0
      id = 0
      do i = 1, nEl
         do j = 1, msh%eNoN
            Tag = gIEN(j,i)
            if (g2l(Tag) .eq. 0) then
               id = id + 1
               g2l(Tag) = id
            end if
         end do
      end do

      ! connectivity
      msh%nEl = nEl
      allocate(msh%IEN(msh%eNoN,msh%nEl))
      do i = 1, nEl
         do j = 1, msh%eNoN
            msh%IEN(j,i) = g2l(gIEN(j,i))
         end do
      end do

      ! coordinates
      msh%nNo = id
      allocate(msh%x(nsd,msh%nNo))
      do i = 1, gmshNodes%numNodes
         if (g2l(i) .ne. 0) then
            msh%x(:,g2l(i)) = gmshNodes%coord(1:nsd,i)
         end if
      end do

      ! deallocate
      deallocate(entityTag, gIEN)

      return
      end subroutine gmsh_bodymesh
!***********************************************************************
!     Assemble boundary meshes
      subroutine gmsh_boundarymesh(msh)
      use varmod
      use gmshMod

      implicit none

      Type(meshType), intent(inout) :: msh

      integer :: i, j, k, l, ii, ie, Tag, iFa, iTag
      integer :: nEl, eNoN, numEntityTags, numinBlocks
      integer,allocatable :: physicalTag(:), gIEN(:,:), lg2l(:)
      integer,allocatable :: sharedelem(:,:), numshared(:), entityTag(:)
      logical,allocatable :: flag(:)

      ! For boundary points, it requires two global-to-local mappings.
      ! g2l: since we need to remove some unnecessary points in gmshNodes, 
      !      this mapping is built to contruct the body mesh
      ! lg2l: when writing boundaries to vtp files, we need a local connectivity
      !      system, this mapping maps **body mesh**, i.e. after g2l, to local.
      allocate(lg2l(msh%nNo))

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

      ! Loop over faces to construct mesh files.
      ! Please note that there might be multiple entity tags
      ! for one physical tag
      do iFa = 1, msh%nFa

         ! Find the corresponding EntityTags for boundary meshes
         ! ndim=2, CurveTag
         ! ndim=3, SurfaceTag
         if (nsd .eq. 2) then
            if (.not. allocated(entityTag)) &
               allocate(entityTag(gmshEntities%numCurves))
            entityTag = 0
            numEntityTags = 0
            do j = 1, gmshEntities%numCurves
               do k = 1, gmshEntities%Curve_numPhysicalTags(j)
                  if (physicalTag(iFa) .eq. gmshEntities%Curve_physicalTags(j,k)) then
                     numEntityTags = numEntityTags + 1
                     entityTag(numEntityTags) = j
                  end if
               end do
            end do
         elseif (nsd .eq. 3) then
            if (.not. allocated(entityTag)) &
               allocate(entityTag(gmshEntities%numSurfaces))
            entityTag = 0
            numEntityTags = 0
            do j = 1, gmshEntities%numSurfaces
               do k = 1, gmshEntities%Surface_numPhysicalTags(j)
                  if (physicalTag(iFa) .eq. gmshEntities%Surface_physicalTags(j,k)) then
                     numEntityTags = numEntityTags + 1
                     entityTag(numEntityTags) = j
                  end if
               end do
            end do
         end if
         if (numEntityTags .eq. 0) then
            write(stdout,ftab4) "ERROR: cannot find corresponding entitiy."
            stop
         end if         

         ! raw connectivity for boundary mesh
         nEl = 0
         do iTag = 1, numEntityTags
            Tag = entityTag(iTag)
            j = 0
            do i = 1, gmshElements%numElementBlocks
               if (gmshElements%ElementDim(i) .eq. nsd-1) then
               if (gmshElements%EntityTag(i) .eq. Tag) then
                  eNoN = gmshElements%eNoN(i)
                  ! msh%fa(iFa)%eNoN = eNoN
                  if (.not. allocated(gIEN)) & 
                     allocate(gIEN(eNoN,gmshElements%numElements))
                  numinBlocks = gmshElements%numElementsinBlock(i)
                  gIEN(1:eNoN,nEl+1:nEl+numinBlocks) = gmshElements%conn(1:eNoN,j+1:j+numinBlocks)
                  nEl = nEl + numinBlocks
                  exit
               end if
               end if
               j = j + gmshElements%numElementsInBlock(i)
            end do
         end do
         
         ! first global to local mapping
         do i = 1, nEl
            do j = 1, eNoN
               gIEN(j,i) = g2l(gIEN(j,i))
            end do
         end do

         ! convert gIEN to local IEN and determine gN
         ! second global to local mapping
         lg2l = 0
         msh%fa(iFa)%nNo  = 0
         msh%fa(iFa)%eNoN = eNoN
         msh%fa(iFa)%nEl  = nEl
         do i = 1, msh%fa(iFa)%nEl
            do j = 1, msh%fa(iFa)%eNoN
               k = gIEN(j,i)
               if (lg2l(k) .eq. 0) then
                  msh%fa(iFa)%nNo = msh%fa(iFa)%nNo + 1
                  lg2l(k) = msh%fa(iFa)%nNo
               end if
            end do
         end do
         allocate(msh%fa(iFa)%gN(msh%fa(iFa)%nNo),msh%fa(iFa)%x(nsd,msh%fa(iFa)%nNo))
         do i = 1, msh%nNo
            if (lg2l(i) .ne. 0) then
               j = lg2l(i)
               msh%fa(iFa)%gN(j) = i
               msh%fa(iFa)%x(:,j) = msh%x(:,i)
            end if
         end do
         allocate(msh%fa(iFa)%IEN(msh%fa(iFa)%eNoN,msh%fa(iFa)%nEl))
         do i = 1, msh%fa(iFa)%nEl
            do j = 1, msh%fa(iFa)%eNoN
               msh%fa(iFa)%IEN(j,i) = lg2l(gIEN(j,i))
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

         ! ! debug
         ! if (msh%fa(iFa)%fname .eq. "bottom") then
         !    print *, msh%fa(iFa)%gE
         !    print *, msh%fa(iFa)%gN
         !    stop
         ! end if

      end do !iFa

      ! deallocate
      deallocate(physicalTag, gIEN, g2l, lg2l)
      deallocate(sharedelem, numshared, entityTag, flag)

      return
      end subroutine gmsh_boundarymesh
!***********************************************************************