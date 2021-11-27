!--------------------------------------------------------------------
!
! All the subroutines that handle p1 to p2 conversion. It takes two mesh
! formats: 1) P1 vtk mesh and 2) P1 text mesh (connectivity + connectivit)
! and convert them to P2 meshes.
! 
!--------------------------------------------------------------------

      module P1P2mod
      use varmod
      
      integer, allocatable :: xadj(:), adjncy(:), mdptID(:)

      contains 

         !==========================================

         ! Generate P2 mesh from P1 mesh -- body
         subroutine P1toP2(lM)
         implicit none 

         type(meshType), intent(inout) :: lM 

         integer :: a, e, i, j, is, ie, ii, jj
         integer :: npts, eNoN
         integer, allocatable :: conn(:,:), indx(:)
         real(kind=8), allocatable :: coord(:,:)

         ! Build adjacency structure use CSR
         call getAdjcncy(lM)
         write(stdout,ftab2) "nEg: "//TRIM(STR(xadj(lM%nNo+1) - 1))

         ! Add mid-point to each edge
         allocate(mdptID(xadj(lM%nNo+1) - 1))
         allocate(coord(nsd,(xadj(lM%nNo+1) - 1)/2+lM%nNo))
         mdptID = 0
         coord = 0D0
         coord(:,1:lM%nNo) = lM%x
         npts = lM%nNo
         do i = 1, lM%nNo
            is = xadj(i)
            ie = xadj(i+1) - 1
            do j = is, ie
               if (i .lt. adjncy(j)) then 
                  npts = npts + 1
                  coord(:,npts) = 5D-1*(coord(:,i)+coord(:,adjncy(j)))
                  mdptID(j) = npts          
               end if
            end do 
         end do 
         ! edge i-->j and edge j-->i should have the same mdptID
         do i = 1, lM%nNo
            is = xadj(i)
            ie = xadj(i+1) - 1
            do j = is, ie 
               if (i .gt. adjncy(j)) then 
                  ii = adjncy(j)
                  jj  = findindx(ii,i)
                  mdptID(j) = mdptID(jj)
               end if 
            end do
         end do  

         ! Build new connectivity matrix
         if (nsd .eq. 3) then 
            eNoN = 10 
            lM%vtkType = 24
         else if (nsd .eq. 2) then 
            eNoN = 6
            lM%vtkType = 22
         end if 
         allocate(conn(eNoN,lM%nEl),indx(lM%eNoN))
         do ie = 1, lM%nEl
            indx = lM%IEN(:,ie)
            conn(1:lM%eNoN,ie) = indx
            if (nsd .eq. 3) then 
               conn(5, ie) = mdptID(findindx(indx(1),indx(2)))
               conn(6, ie) = mdptID(findindx(indx(2),indx(3)))
               conn(7, ie) = mdptID(findindx(indx(1),indx(3)))
               conn(8, ie) = mdptID(findindx(indx(1),indx(4)))
               conn(9, ie) = mdptID(findindx(indx(2),indx(4)))
               conn(10,ie) = mdptID(findindx(indx(3),indx(4)))
            else if (nsd .eq. 2) then 
               conn(4, ie) = mdptID(findindx(indx(1),indx(2)))
               conn(5, ie) = mdptID(findindx(indx(2),indx(3)))
               conn(6, ie) = mdptID(findindx(indx(1),indx(3)))  
            end if          
         end do 

         ! Update body mesh information 
         deallocate(lM%x, lM%IEN)
         lM%nNo     = npts
         lM%eNoN    = eNoN
         allocate(lM%x(nsd,lM%nNo), lM%IEN(lM%eNoN,lM%nEl))
         lM%x       = coord 
         lM%IEN     = conn 

         deallocate(coord, conn, indx)
         end subroutine P1toP2

         !==========================================

         ! Generate P2 mesh from P1 mesh -- boundary
         subroutine P1toP2B(lM, lFa, insd)
         implicit none 

         type(meshType), intent(in)    :: lM 
         type(faceType), intent(inout) :: lFa 
         integer, intent(in) :: insd

         integer :: i, npts, eNoN, ie, ge
         integer, allocatable, dimension(:) :: indx, gN
         integer, allocatable, dimension(:,:)  :: conn
         integer :: min_val, max_val
         
         ! Build new connectivity matrix 
         if (insd .eq. 2) then 
            eNoN = 6 
            lFa%vtkType = 22
         else if (insd .eq. 1) then 
            eNoN = 3
            lFa%vtkType = 21
         end if 
         allocate(conn(eNoN,lFa%nEl), indx(lFa%eNoN))

         if (insd .eq. 2) then 
            do ie = 1, lFa%nEl
               indx = lFa%IEN(:,ie)
               conn(1:lFa%eNoN, ie) = indx
               conn(4, ie) = mdptID(findindx(indx(1),indx(2)))
               conn(5, ie) = mdptID(findindx(indx(2),indx(3)))
               conn(6, ie) = mdptID(findindx(indx(1),indx(3)))  
            end do               
         else if (insd .eq. 1) then 
            do ie = 1, lFa%nEl
               indx = lFa%IEN(:,ie)
               conn(1:lFa%eNoN, ie) = indx
               conn(3, ie) = mdptID(findindx(indx(1),indx(2)))
            end do 
         end if 

         ! Number of points and gN
         deallocate(indx)
         allocate(indx(eNoN*lFa%nEl), gN(eNoN*lFa%nEl))
         indx = RESHAPE(conn,(/eNoN*lFa%nEl/))
         min_val = MINVAL(indx) - 1
         max_val = MAXVAL(indx)
         npts = 0
         do while (min_val < max_val)
            npts = npts + 1
            min_val = MINVAL(indx, mask = indx>min_val)
            gN(npts) = min_val
         end do

         ! Update surface mesh information 
         deallocate(lFa%IEN)
         if (allocated(lFa%gN)) deallocate(lFa%gN)
         lFa%nNo  = npts
         lFa%eNoN = eNoN
         allocate(lFa%IEN(lFa%eNoN,lFa%nEl),lFa%gN(lFa%nNo))
         lFa%gN  = gN(1:npts)
         ! Map face IEN to local numbering to write VTP files
         deallocate(indx)
         allocate(indx(lM%nNo))
         indx = 0 
         do i = 1, npts
            indx(gN(i)) = i 
         end do 
         do ie = 1, lFa%nEl
            lFa%IEN(:,ie) = indx(conn(:,ie))
         end do 
         if (allocated(lFa%x)) deallocate(lFa%x)
         allocate(lFa%x(nsd,lFa%nNo))
         do i = 1, npts 
            lFa%x(:,i) = lM%x(:,lFa%gN(i))
         end do

         deallocate(conn, indx, gN)

         end subroutine P1toP2B

         !==========================================

         function findindx(ii,jj) result (indx)
         implicit none 
         integer :: ii, jj, indx
         integer :: is, ie, j

         is = xadj(ii)
         ie = xadj(ii+1) - 1
         indx = 0
         do j = is,ie
            if (adjncy(j) .eq. jj) then
               indx = j
               exit
            end if
         end do
         !print *, indx, FINDLOC(adjncy(is:ie),jj,DIM=1) + is - 1, is - 1
         end function findindx
         
         !==========================================
      
         subroutine checkIEN(lM,insd)
         implicit none 

         type(meshType), intent(inout) :: lM 
         integer, intent(in) :: insd

         integer :: ie
         integer, allocatable, dimension(:) :: indx

         real(kind=8) ::  vol
         real(kind=8), allocatable, dimension(:,:) :: xx
         real(kind=8), allocatable, dimension(:,:) :: rr

         allocate(xx(insd,lM%eNoN), indx(lM%eNoN))
         allocate(rr(insd,insd))

         if (insd .eq. 3) then
            do ie = 1, lM%nEl
               indx = lM%IEN(:,ie)
               xx   = lM%x(:,indx)
               rr(:,1)   = xx(:,2) - xx(:,1)
               rr(:,2)   = xx(:,3) - xx(:,1)
               rr(:,3)   = xx(:,4) - xx(:,1)
               vol  = NORM(rr(:,3),CROSS(rr(:,1:2)))
               if (vol .gt. 0) then 
                  lM%IEN(lM%eNoN-1,ie) = indx(lM%eNoN)
                  lM%IEN(lM%eNoN  ,ie) = indx(lM%eNoN-1)
               end if 
            end do
         else if (insd .eq. 2) then 
            do ie = 1, lM%nEl
               indx = lM%IEN(:,ie)
               xx   = lM%x(:,indx)
               rr(:,1)   = xx(:,2) - xx(:,1)
               rr(:,2)   = xx(:,3) - xx(:,1)
               vol  = rr(1,1)*rr(2,2) - rr(1,2)*rr(2,1)
               if (vol .lt. 0) then 
                  lM%IEN(lM%eNoN-1,ie) = indx(lM%eNoN)
                  lM%IEN(lM%eNoN  ,ie) = indx(lM%eNoN-1)
               end if 
            end do            
         end if

         deallocate(xx, indx, rr)
         end subroutine checkIEN

         !==========================================

         subroutine checkIENB(lM,lFa,insd)
         implicit none 

         type(meshType), intent(in)    :: lM 
         type(faceType), intent(inout) :: lFa 
         integer, intent(in) :: insd

         integer :: ie, ge
         integer, allocatable, dimension(:) :: indx

         real(kind=8) ::  vol
         real(kind=8), allocatable :: xx(:,:), rr(:,:), center(:)

         allocate(xx(nsd,lFa%eNoN), indx(lFa%eNoN))
         allocate(rr(nsd,nsd), center(nsd))

         ! Reorient the surface element so that the norm calculated  
         ! from cross(pt2-pt1,pt3-pt1) is pointing outward
         if (insd .eq. 2) then 
            do ie = 1, lFa%nEl
               indx = lFa%IEN(:,ie)
               ge   = lFa%gE(ie)
               xx   = lM%x(:,indx)
               center = SUM(lM%x(:,lM%IEN(:,ge)),DIM = 2)/lM%eNoN

               rr(:,1)= xx(:,2) - xx(:,1)
               rr(:,2)= xx(:,3) - xx(:,1) 
               rr(:,3)= center  - xx(:,1)
               vol  = NORM(rr(:,3),CROSS(rr(:,1:2)))
               if (vol .gt. 0) then 
                  lFa%IEN(lFa%eNoN-1,ie) = indx(lFa%eNoN)
                  lFa%IEN(lFa%eNoN  ,ie) = indx(lFa%eNoN-1)
               end if 
            end do
         else if (insd .eq. 1) then 
            do ie = 1, lFa%nEl
               indx = lFa%IEN(:,ie)
               ge   = lFa%gE(ie)
               xx   = lM%x(:,indx)
               center = SUM(lM%x(:,lM%IEN(:,ge)),DIM = 2)/lM%eNoN

               rr(:,1)   = xx(:,2) - xx(:,1)
               rr(:,2)   = center  - xx(:,1)
               vol  = rr(1,1)*rr(2,2) - rr(1,2)*rr(2,1)
               if (vol .lt. 0) then 
                  lFa%IEN(lFa%eNoN-1,ie) = indx(lFa%eNoN)
                  lFa%IEN(lFa%eNoN  ,ie) = indx(lFa%eNoN-1)
               end if 
            end do
         end if

         deallocate(xx, indx, rr, center)
         end subroutine checkIENB

         !==========================================
         
         subroutine getAdjcncy(lM)
         implicit none 
         
         type(meshType), intent(in) :: lM 
         
         integer :: maxedge = 40 
         integer :: a, b, e, Ac, Bc, i, j
         integer, allocatable :: ncount(:), edges(:,:)
         logical :: flag
         
         allocate(edges(maxedge,lM%nNo), ncount(lM%nNo))
         edges = 0 
         ncount = 0 

         ! Determine edge information
         do e = 1, lM%nEl 
            do a = 1, lM%eNoN
               Ac = lM%IEN(a,e)
               do b = 1, lM%eNoN 
                  Bc = lM%IEN(b,e)
                  if (Ac .eq. Bc) cycle
                  flag = .true.
                  j = 0
                  LOOP: do i = 1, ncount(Ac)
                     if (Bc .eq. edges(i,Ac)) then 
                        flag = .false.
                        exit LOOP
                     else if (Bc .gt. edges(i,Ac)) then 
                        j = i
                     end if
                  end do LOOP
                  if (flag) then 
                     if (ncount(Ac) .eq. maxedge) then 
                        write(stdout, ftab4) "ERROR: maxedges is too small."
                        stop 
                     end if 
                     edges(j+2:ncount(Ac)+1,Ac) = edges(j+1:ncount(Ac),Ac)
                     edges(j+1,Ac) = Bc
                     ncount(Ac) = ncount(Ac) + 1
                  end if
               end do 
            end do 
         end do 

         ! Build adjacency array using CSR
         allocate(xadj(lM%nNo+1))
         xadj = 0
         xadj(1) = 1
         do Ac = 2, lM%nNo+1
            xadj(Ac) = xadj(Ac-1) + ncount(Ac-1)
         end do
         allocate(adjncy(xadj(lM%nNo+1) - 1))
         adjncy = 0 
         do Ac = 1, lM%nNo
            a = xadj(Ac)
            b = xadj(Ac+1) - 1
            adjncy(a:b) = edges(1:ncount(Ac),Ac)
         end do

         deallocate(edges, ncount)
         end subroutine getAdjcncy

         !==========================================

         !     Returning number of data/words in a string
         pure function CheckNoNumbers(sTmp)
         use varmod
         implicit none
         character(len=strL), intent(in) :: sTmp

         logical isBnk
         integer CheckNoNumbers, i

         isBnk = .TRUE.
         CheckNoNumbers = 0
         do i=strL, 1, -1
            if (isBnk) then
               if (sTmp(i:i) .NE. " ") then
                  CheckNoNumbers = CheckNoNumbers + 1
                  isBnk = .FALSE.
               end if
            else
               if (sTmp(i:i) .EQ. " ") then
                  isBnk = .TRUE.
               end if
            end if
         end do

         return
         end function CheckNoNumbers
         
      end module P1P2mod
   
!***********************************************************************

      subroutine conv_vtk(msh)
      use P1P2mod 
      implicit none
      
      type(meshType), intent(inout) :: msh
      character(len=strL) :: fname

      integer :: ifa

      if ( (nsd .eq. 2 .and. msh%vtkType .ne. 5 ) .or. &
      &    (nsd .eq. 3 .and. msh%vtkType .ne. 10) ) then 
         write(stdout, ftab4) "ERROR: P2 conversion only works for triangle (2D) or tetrahedron (3D)."
         stop 
      end if 
      
      write(stdout, ftab1)
      write(stdout, ftab1) "Converting to P2 mesh.."

      ! Check orientation of the elements
      call checkIEN(msh,nsd)
      do ifa = 1, msh%nFa
         call checkIENB(msh,msh%fa(ifa),nsd-1)
      end do 

      !==============================================================
      !      Convert body mesh
      !==============================================================
      call P1toP2(msh)

      !==============================================================
      !      Convert surface mesh
      !==============================================================
      do ifa = 1, msh%nFa
         call P1toP2B(msh, msh%fa(ifa), nsd-1)
      end do 

      deallocate(xadj, adjncy, mdptID)
         
      end subroutine conv_vtk

!***********************************************************************

      subroutine readccne(msh)
      use P1P2mod
      implicit none 

      type(meshType), intent(inout) :: msh 

      character(len=strL) :: temp
      integer :: fid, i, e, iFa

      fid = 100 
      open(fid, file = trim(msh%fname)//'.connectivity')
      read(fid,'(A)') temp 
      msh%eNoN = CheckNoNumbers(temp) - 1
      call selectel(msh)
      

      !==============================================================
      !      Read body mesh
      !==============================================================      
      ! reading connectivity file
      msh%nEl = 1 
      do 
         read(fid,*,end=113)
         msh%nEl = msh%nEl + 1
      end do 
113   rewind(fid)

      allocate (msh%IEN(msh%eNoN,msh%nEl))
      do e=1, msh%nEl
         read (fid,*) i, msh%IEN(:,e)
      end do
      close (fid)

      ! reading coordinates file
      open(fid, file = trim(msh%fname)//'.coordinates')
      read (fid,"(A)") temp
      if (CheckNoNumbers(temp) .NE. nsd+1) THEN
         write(stdout, ftab4) "Error: nsd is not consistent with coordinate file"
      end if
      msh%nNo = 1
      do
         read (fid,*,end=111)
         msh%nNo = msh%nNo + 1
      end do
111   rewind(fid)

      allocate (msh%x(nsd,msh%nNo))
      do i=1, msh%nNo
         read (fid,*) e, msh%x(:,i)
      end do
      close (fid)


      !==============================================================
      !      Read surface mesh
      !==============================================================
      do iFa=1, msh%nFa
         ! reading face connectivity file
         open(fid, file = trim(msh%fa(iFa)%fname))
         read(fid,"(A)",end=110) temp
         msh%fa(iFa)%eNoN = CheckNoNumbers(temp) - 2
         msh%fa(iFa)%nEl = 1
         do
            read (fid,*,end=110)
            msh%fa(iFa)%nEl = msh%fa(iFa)%nEl + 1
         end do
110      allocate (msh%fa(iFa)%gE(msh%fa(iFa)%nEl))
         allocate (msh%fa(iFa)%IEN(msh%fa(iFa)%eNoN,msh%fa(iFa)%nEl))
         rewind(fid)
         do e=1, msh%fa(iFa)%nEl
            read (fid,*) msh%fa(iFa)%gE(e),i,msh%fa(iFa)%IEN(:,e)
         end do
         close (fid)
      end do

      write(stdout,ftab2) "Nsd: "//TRIM(STR(nsd))
      write(stdout,ftab2) "nNo: "//TRIM(STR(msh%nNo))
      write(stdout,ftab2) "nEl: "//TRIM(STR(msh%nEl))
      write(stdout,ftab2) "nFa: "//TRIM(STR(msh%nFa))

      end subroutine readccne

!***********************************************************************


   
