!--------------------------------------------------------------------
!
! All the subroutines that handle Gambit Neutral Mesh Files 
!
!--------------------------------------------------------------------

      subroutine conv_gmsh2vtk(msh)
      use varmod
      implicit none

      type(meshType), intent(out)  :: msh

      integer :: a, b, e, i, iFa, Ac, Ec, eOrd, fid, nToks, eNoNb
      character(len=strL) :: rLine, tokenList(maxToks)

      integer, allocatable :: ptr(:)

      fid = 100
      write(stdout,ftab1) "Loading file "//TRIM(msh%fname)
      open(fid, file=trim(msh%fname))
      call findKwrd(fid, "NUMNP")
      read(fid,*) msh%nNo, msh%nEl, i, msh%nFa, a, b

      nsd = max(a,b)
      write(stdout,ftab2) "Nsd: "//TRIM(STR(nsd))
      write(stdout,ftab2) "nNo: "//TRIM(STR(msh%nNo))
      write(stdout,ftab2) "nEl: "//TRIM(STR(msh%nEl))
      write(stdout,ftab2) "nFa: "//TRIM(STR(msh%nFa))

!     Read nodal coordinates
      allocate(msh%x(nsd,msh%nNo))
      call findKwrd(fid, "NODAL")
      do a=1, msh%nNo
         read(fid,*) i, msh%x(:,a)
      end do

!     Read element connectivity/IEN
      call findKwrd(fid, "ELEMENTS/CELLS")
      read(fid,'(A)') rLine
      call parseString(rLine, tokenList, nToks)
      if (nToks .eq. 0) then
         write(stdout,ftab4) &
            "ERROR: could not parse element connectivity"
         STOP
      end if
      read(tokenList(3),*) msh%eNoN

      allocate(msh%IEN(msh%eNoN, msh%nEl))
      if (nToks-3 .eq. msh%eNoN) then
         e = 1
         do a=1, msh%eNoN
            read(tokenList(3+a),*) msh%IEN(a,e)
         end do
         do e=2, msh%nEl
            read(fid,*) i, i, i, msh%IEN(:,e)
         end do
      else
         rewind(fid)
         call findKwrd(fid, "ELEMENTS/CELLS")
         do e=1, msh%nEl
            read(fid,*) i, i, i, msh%IEN(1:ntoks-3,e)
            read(fid,*) msh%IEN(ntoks-3+1:msh%eNoN,e)
         end do
      end if

      if (nsd.eq.3 .and. msh%eNoN.eq.8) then
         do e=1, msh%nEl
            a = 3
            b = 4
            i = msh%IEN(a,e)
            msh%IEN(a,e) = msh%IEN(b,e)
            msh%IEN(b,e) = i
            a = 7
            b = 8
            i = msh%IEN(a,e)
            msh%IEN(a,e) = msh%IEN(b,e)
            msh%IEN(b,e) = i
         end do
      end if

      allocate(msh%fa(msh%nFa))
      call selectel(msh)

!     Read face data
      do iFa=1, msh%nFa
         eNoNb = msh%fa(iFa)%eNoN
         call findKwrd(fid, "BOUNDARY")
         read(fid,'(A)') rLine
         call parseString(rLine, tokenList, nToks)
         read(tokenList(1),*) msh%fa(iFa)%fname
         write(stdout,ftab3) "Face <"//TRIM(msh%fa(iFa)%fname)//">"
         if (TRIM(tokenList(2)) .ne. "1") then
            write(stdout,ftab4) &
               "ERROR: element information not found on boundary"
            STOP
         end if
         read(tokenList(3),*) msh%fa(iFa)%nEl
         allocate(msh%fa(iFa)%gE(msh%fa(iFa)%nEl))
         allocate(msh%fa(iFa)%IEN(eNoNb, msh%fa(iFa)%nEl))
         do e=1, msh%fa(iFa)%nEl
            read(fid,*) msh%fa(iFa)%gE(e), i, eOrd
            Ec = msh%fa(iFa)%gE(e)
            do a=1, eNoNb
               i = msh%eOrd(a,eOrd)
               msh%fa(iFa)%IEN(a,e) = msh%IEN(i,Ec)
            end do
         end do
         read(fid,*)
      end do
      close(fid)

!     Convert to biquadratic elements if needed
      if (nsd.eq.2 .and. msh%eNoN.eq.4) then
         write(stdout,ftab1,advance='no') &
            "Detected a quad mesh. Convert to biquadratic? (y/n)  "
         read(*,'(A)') rLine
         rLine = adjustl(trim(rLine))
         if (rLine.eq.'y' .or. rLine.eq.'yes') then
            msh%nonL = .true.
            call conv_2biquad(msh)
         end if
      end if

!     Setup face data stucture
      allocate(ptr(msh%nNo))
      do iFa=1, msh%nFa
!     Get face gN and x data
         ptr = 0
         msh%fa(iFa)%nNo = 0
         do e=1, msh%fa(iFa)%nEl
            do a=1, msh%fa(iFa)%eNoN
               Ac = msh%fa(iFa)%IEN(a,e)
               if (ptr(Ac) .eq. 0) then
                  msh%fa(iFa)%nNo = msh%fa(iFa)%nNo + 1
                  ptr(Ac) = 1
               end if
            end do
         end do

         allocate(msh%fa(iFa)%gN(msh%fa(iFa)%nNo))
         allocate(msh%fa(iFa)%x(nsd,msh%fa(iFa)%nNo))
         a = 0
         do Ac=1, msh%nNo
            if (ptr(Ac) .gt. 0) then
               a = a + 1
               msh%fa(iFa)%gN(a) = Ac
               msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
            end if
         end do

!     Map face IEN to local numbering to write VTP files
         ptr = 0
         do a=1, msh%fa(iFa)%nNo
            Ac = msh%fa(iFa)%gN(a)
            ptr(Ac) = a
         end do

         do e=1, msh%fa(iFa)%nEl
            do a=1, msh%fa(iFa)%eNoN
               Ac = msh%fa(iFa)%IEN(a,e)
               msh%fa(iFa)%IEN(a,e) = ptr(Ac)
            end do
         end do
      end do

      end subroutine conv_gmsh2vtk

!***********************************************************************

      subroutine findKwrd(fid, kwrd)
      use varmod
      implicit none
      integer, intent(in) :: fid
      character(len=*), intent(in) :: kwrd

      integer :: itok, ntoks
      character(len=strL) :: rLine, tokenList(maxToks)

      MY_LOOP : do
         read(fid,'(A)',end=001) rLine
         call parseString(rLine, tokenList, nToks)
         if (nToks .gt. 0) then
            do itok=1, nToks
               if (TRIM(tokenList(itok)) .eq. TRIM(kwrd)) exit MY_LOOP
            end do
         end if
      end do MY_LOOP

      return

   001  write(stdout,ftab4) "ERROR: end of file reached.."
      write(stdout,ftab4) 'Failed processing for keyword "'// &
            trim(kwrd)//'"'
      STOP
      end subroutine findKwrd

!***********************************************************************

      subroutine conv_2biquad(lM)
      use varmod
      implicit none

      type(meshType), intent(inout) :: lM

      integer a, e, i, j, k, Ac, Ec, iFa, eNoNb

      integer, allocatable :: adj(:,:), tmpI(:,:)
      real(kind=8), allocatable :: Xtmp(:,:)

      write(stdout,ftab1)
      write(stdout,ftab1) "Converting to biquadratic mesh.."

      allocate(adj(lM%eNoF,lM%nEl))
      adj = 0
      call getAdjcncy(lM, adj)

!     Reset eNoN and vtkType for 9-node biquadratic elements
      lM%eNoN = 9
      lM%vtkType = 28

!     Adjust array sizes based on new connectivity
      allocate(tmpI(4,lM%nEl))
      do e=1, lM%nEl
         do a=1, 4
            tmpI(a,e) = lM%IEN(a,e)
         end do
      end do
      deallocate(lM%IEN)
      allocate(lM%IEN(lM%eNoN,lM%nEl))
      lM%IEN = 0
      do e=1, lM%nEl
         do a=1, 4
            lM%IEN(a,e) = tmpI(a,e)
         end do
      end do
      deallocate(tmpI)

!     Upper estimation for nNo
      allocate(Xtmp(nsd,lM%eNoN*lM%nEl))
      Xtmp = 0D0
      do a=1, lM%nNo
         do i=1, nsd
            Xtmp(i,a) = lM%x(i,a)
         end do
      end do
      deallocate(lM%x)

!     We will start adding intermediate nodes from faces. Here I assume
!     there is no overlapping between faces
      do iFa=1, lM%nFa
         lM%fa(iFa)%eNoN = 3
         lM%fa(iFa)%vtkType = 21
         eNoNb = lM%fa(iFa)%eNoN

!     Reset face IEN and gE structure
         allocate(tmpI(2,lM%fa(iFa)%nEl))
         do e=1, lM%fa(iFa)%nEl
            tmpI(1,e) = lM%fa(iFa)%IEN(1,e)
            tmpI(2,e) = lM%fa(iFa)%IEN(2,e)
         end do
         deallocate(lM%fa(iFa)%IEN)
         allocate(lM%fa(iFa)%IEN(eNoNb,lM%fa(iFa)%nEl))
         lM%fa(iFa)%IEN = 0
         do e=1, lM%fa(iFa)%nEl
            lM%fa(iFa)%IEN(1,e) = tmpI(1,e)
            lM%fa(iFa)%IEN(2,e) = tmpI(2,e)
         end do
         deallocate(tmpI)

         do e=1, lM%fa(iFa)%nEl
            i  = lM%fa(iFa)%IEN(1,e)
            j  = lM%fa(iFa)%IEN(2,e)
            Ec = lM%fa(iFa)%gE(e)
            call addNodesBetween(lM, i, j, Ec, .true., Xtmp)
            lM%fa(iFa)%IEN(3,e) = lM%nNo
            do k=1, lM%eNoF
               Ac = adj(k,Ec)
               if (Ac .eq. 0) cycle
               call addNodesBetween(lM, i, j, Ac, .false., Xtmp)
            end do
         end do
      end do

!     Now adding nodes for interior elements
      do e=1, lM%nEl
         do a=1, 4
            i = lM%IEN(a,e)
            j = lM%IEN(a+1,e)
            if (a .eq. 4) j = lM%IEN(1,e)

!           ignore boundary edges
            if (adj(a,e) .eq. 0) cycle

            call addNodesBetween(lM, i, j, e, .true., Xtmp)
            do k=1, lM%eNoF
               Ec = adj(k,e)
               if (Ec .eq. 0) cycle
               call addNodesBetween(lM, i, j, Ec, .false., Xtmp)
            end do
         end do
      end do

!     Add the node in the middle of the element
      do e=1, lM%nEl
         lM%nNo = lM%nNo + 1
         lM%IEN(9,e) = lM%nNo
         i = lM%IEN(1,e)
         j = lM%IEN(3,e)
         Xtmp(:,lM%nNo) = (Xtmp(:,i) + Xtmp(:,j))/2.0D0
      end do

      allocate(lM%x(nsd,lM%nNo))
      do a=1, lM%nNo
         lM%x(:,a) = Xtmp(:,a)
      end do
      deallocate(Xtmp, adj)

      write(stdout,ftab2) "nNo: "//TRIM(STR(lM%nNo))
      write(stdout,ftab2) "nEl: "//TRIM(STR(lM%nEl))

      contains

         !==========================================

         subroutine getAdjcncy(lM, adj)
         use varmod
         implicit none

         type(meshType), intent(in) :: lM
         integer, intent(inout) :: adj(lM%eNoF,lM%nEl)

         logical :: flag
         integer :: a, b, e, f, i, j, Ac, Bc, Ec, maxAdj

         integer, allocatable :: incNd(:), tmpI(:,:)

         allocate(incNd(lM%nNo))
         incNd = 0
         do e=1, lM%nEl
            do a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               incNd(Ac) = incNd(Ac) + 1
            end do
         end do

         maxAdj = MAXVAL(incNd)
         ALLOCATE(tmpI(maxAdj,lM%nNo))
         incNd = 0
         tmpI  = 0
         do e=1, lM%nEl
            do a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               incNd(Ac) = incNd(Ac) + 1
               tmpI(incNd(Ac), Ac) = e
            end do
         end do

         do e=1, lM%nEl
            do f=1, lM%eNoF
               a  = lM%eOrd(1,f)
               b  = lM%eOrd(2,f)
               Ac = lM%IEN(a,e)
               Bc = lM%IEN(b,e)
               flag = .false.
               LOOP: do i=1, incNd(Ac)
                  Ec = tmpI(i,Ac)
                  if (Ec .eq. e) cycle
                  do j=1, incNd(Bc)
                     if (tmpI(j,Bc) .eq. e) cycle
                     if (Ec .eq. tmpI(j,Bc)) then
                        flag = .true.
                        exit LOOP
                     end if
                  end do
               end do LOOP
               if (flag) adj(f,e) = Ec
            end do
         end do
         deallocate(tmpI, incNd)
         end subroutine getAdjcncy

         !==========================================

         subroutine addNodesBetween(lM, Ac1, Ac2, e, flag, Xtmp)
         use varmod
         implicit none

         type(meshType), intent(inout) :: lM
         integer, intent(in) :: Ac1, Ac2, e
         logical, intent(in) :: flag
         real(kind=8), intent(inout) :: Xtmp(nsd,lM%eNoN*lM%nEl)

         integer a, b1, b2

         b1 = 0
         b2 = 0
         do a=1, 4
            if (lM%IEN(a,e) .eq. Ac1) then
               b1 = a
            else if(lM%IEN(a,e) .eq. Ac2) then
               b2 = a
            end if
         end do
         if (b1 .gt. b2) then
            a  = b1
            b1 = b2
            b2 = a
         end if
         if (b1 .eq. 0) return

         if (b1.eq.1 .and. b2.eq.2) then
            a = 5
         else if (b1.eq.2 .and. b2.eq.3) then
            a = 6
         else if (b1.eq.3 .and. b2.eq.4) then
            a = 7
         else if (b1.eq.1 .and. b2.eq.4) then
            a = 8
         end if

         if (lM%IEN(a,e) .eq. 0) then
            if (flag) then
               lM%nNo = lM%nNo + 1
               Xtmp(:,lM%nNo) = (Xtmp(:,Ac1) + Xtmp(:,Ac2))/2.0D0
            end if
            lM%IEN(a,e) = lM%nNo
         end if
         end subroutine addNodesBetween

         !==========================================

      end subroutine conv_2biquad

!***********************************************************************