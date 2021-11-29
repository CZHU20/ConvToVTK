!--------------------------------------------------------------------
!
! This program is a simple utility to 
! 1. convert Gambit Neutral Mesh Files tp VTK based files
!  (VTU: Unstructured Grid / VTP: Polydata)
! 2. convert P1 mesh to P2 mesh
!
!--------------------------------------------------------------------

      module varmod
      use stdParams
      use genUtils

      type faceType
         integer :: nNo, nEl, eNoN, eNoE, vtkType
         integer, allocatable :: IEN(:,:), gN(:), gE(:)
         real(kind=8), allocatable :: x(:,:)
         character(len=strL) :: fname
      end type faceType

      type meshType
         logical :: nonL = .FALSE.
         integer :: nNo, nEl, eNoN, eNoF, vtkType, nFa
         integer, allocatable :: IEN(:,:), eOrd(:,:)
         real(kind=8), allocatable :: x(:,:)
         type(faceType), allocatable :: fa(:)
         character(len=strL) :: fname
      end type meshType

      integer :: nsd

      end module varmod

!***********************************************************************

      program convert_mesh
      use varmod
      implicit none 

      type(meshType) :: msh
      character(len=strL) :: fName, temp 

      integer :: i, j, k, fid, nlines, io

      i = IARGC()
      if (i .eq. 0) then
         write(stdout,ftab4) "ERROR: Input file name not specified"
         STOP
      else if (i .gt. 1) then
         write(stdout,ftab4) "ERROR: Too many arguments"
         STOP
      end if
      call getarg(1,fName)

      if ( endswith(trim(fName), '.neu') ) then 
         write(stdout, ftab1) "Converting Gambit neu file to vtk.."
         msh%fname = fName
         call conv_gneu2vtk(msh)
      elseif ( endswith(trim(fName), '.msh') ) then
         write(stdout,ftab1) "Converting GMSH msh file to vtk.."
         msh%fname = fName
         call conv_gmsh2vtk(msh)
      else
         write(stdout,ftab1,advance='no') &
            "Dimension of the body mesh (2/3):  "
         read(*,*) nsd
         fid = 100
         open(fid,file=trim(fName))
         read(fid,'(A)') temp 
         
         if ( endswith(trim(temp), '.vtu') ) then 
            ! mesh files are in vtk format
            rewind(fid)
            nlines = 0
            write(stdout, ftab1) "Reading list of vtk files.."
            do
               read(fid,'(A)',iostat=io) temp
               if (io/=0) exit
               if (len(trim(temp)) .gt. 0) nlines = nlines + 1
            end do 

            rewind(fid)
            read(fid,'(A)') msh%fname
            msh%nFa = nlines - 1
            allocate(msh%fa(msh%nFa))
            do i = 1, msh%nFa
               read(fid,'(A)') msh%fa(i)%fname
            end do 
            close(fid)

            ! read P1 mesh
            call readVTU(msh,msh%fname)
            call selectel(msh)
            do i = 1, msh%nFa
               call readVTP(msh%fa(i),msh%fa(i)%fname)
            end do 
            write(stdout,ftab2) "Nsd: "//TRIM(STR(nsd))
            write(stdout,ftab2) "nNo: "//TRIM(STR(msh%nNo))
            write(stdout,ftab2) "nEl: "//TRIM(STR(msh%nEl))
            write(stdout,ftab2) "nFa: "//TRIM(STR(msh%nFa))

            msh%nonL = .true.
           ! if read in vtk files, assuming they will be converted to P2 mesh
           ! vtk compatible IEN
           ! https://raw.githubusercontent.com/lorensen/VTKExamples/master/src/Testing/Baseline/Cxx/GeometricObjects/TestIsoparametricCellsDemo.png
            call conv_vtk(msh)

            ! prepare fname for output
            j = SCAN(msh%fname,'/',.True.)
            msh%fname = trim(msh%fname(j+1:))
            do i = 1, msh%nFa
               j = SCAN(msh%fa(i)%fname,'/',.True.)
               k = SCAN(msh%fa(i)%fname,'.',.True.)
               msh%fa(i)%fname = trim(msh%fa(i)%fname(j+1:k-1))
            end do 

         else if ( endswith(trim(temp), '.coordinates') .or. &
                 & endswith(trim(temp), '.connectivity')) then 
            ! mesh files are in plain text format
            write(stdout, ftab1) "Reading list of coordinate/connectivity based mesh files.."
            read(fid,'(A)') temp 
            i = SCAN(temp,'.',.True.)
            msh%fname = temp(:i-1)

            nlines = 0
            do
               read(fid,'(A)',iostat=io) temp
               if (io/=0) exit
               if (len(trim(temp)) .gt. 0) nlines = nlines + 1
            end do 
            rewind(fid)
            msh%nFa = nlines
            allocate(msh%fa(msh%nFa))
            read(fid,'(A)') temp
            read(fid,'(A)') temp
            do i = 1, msh%nFa
               read(fid,'(A)') msh%fa(i)%fname
            end do 
            close(fid)

            call readccne(msh)

            msh%nonL = .true.
            ! if read in vtk files, assuming they will be converted to P2 mesh
            ! vtk compatible IEN
            ! https://raw.githubusercontent.com/lorensen/VTKExamples/master/src/Testing/Baseline/Cxx/GeometricObjects/TestIsoparametricCellsDemo.png
             call conv_vtk(msh)

            ! prepare fname for output
             j = SCAN(msh%fname,'/',.True.)
             msh%fname = trim(msh%fname(j+1:))//'.vtu'
             do i = 1, msh%nFa
                j = SCAN(msh%fa(i)%fname,'/',.True.)
                k = SCAN(msh%fa(i)%fname,'.',.True.)
                msh%fa(i)%fname = trim(msh%fa(i)%fname(j+1:k-1))
             end do 

         end if
      end if 

!     Now, write VTU/VTP files
      write(stdout, ftab1) "Writing to vtk file system.."
      call VTK(msh)

      end program convert_mesh

!***********************************************************************

      subroutine selectel(lM)
      use varmod
      implicit none
      type(meshType), intent(inout) :: lM

      integer :: iFa, eNoNb

!     select mesh element type
      if (nsd .eq. 2) then
         select case (lM%eNoN)
         case (3)   ! tri !
            lM%eNoF = 3
            lM%vtkType = 5
         case (4)   ! quad !
            lM%eNoF = 4
            lM%vtkType = 9
         case default
            write(stdout,ftab4) &
               "ERROR: mesh element type not defined"
            STOP
         end select
      else
         select case (lM%eNoN)
         case (3)   ! tri !
            lM%eNoF = 3
            lM%vtkType = 5
         case (4)   ! tet !
            lM%eNoF = 4
            lM%vtkType = 10
         case (8)   ! hex !
            lM%eNoF = 6
            lM%vtkType = 12
         case default
            write(stdout,ftab4) &
               "ERROR: mesh element type not defined"
            STOP
         end select
      end if

!     select face element type
      do iFa=1, lM%nFa
         if (nsd .eq. 2) then
            select case (lM%eNoN)
            case (3, 4)  ! tri, quad !
               lM%fa(iFa)%eNoN = 2
               lM%fa(iFa)%eNoE = 1
               lM%fa(iFa)%vtkType = 3
            case default
               write(stdout,ftab4) &
                  "ERROR: face element type not defined"
               STOP
            end select
         else
            select case (lM%eNoN)
            case (3) ! tri !
               lM%fa(iFa)%eNoN = 2
               lM%fa(iFa)%eNoE = 1
               lM%fa(iFa)%vtkType = 3
            case (4) ! tet !
               lM%fa(iFa)%eNoN = 3
               lM%fa(iFa)%eNoE = 3
               lM%fa(iFa)%vtkType = 5
            case (8) ! hex !
               lM%fa(iFa)%eNoN = 4
               lM%fa(iFa)%eNoE = 4
               lM%fa(iFa)%vtkType = 9
            case default
               write(stdout,ftab4) &
                  "ERROR: face element type not defined"
               STOP
            end select
         end if
      end do ! iFa

      eNoNb = lM%fa(1)%eNoN
      allocate(lM%eOrd(eNoNb, lM%eNoF))
      if (nsd .eq. 2) then
         select case (lM%eNoN)
         case(3)
            lM%eOrd(1,1) = 1
            lM%eOrd(2,1) = 2
            lM%eOrd(1,2) = 2
            lM%eOrd(2,2) = 3
            lM%eOrd(1,3) = 3
            lM%eOrd(2,3) = 1
         case (4)  ! quad !
            lM%eOrd(1,1) = 1
            lM%eOrd(2,1) = 2
            lM%eOrd(1,2) = 2
            lM%eOrd(2,2) = 3
            lM%eOrd(1,3) = 3
            lM%eOrd(2,3) = 4
            lM%eOrd(1,4) = 4
            lM%eOrd(2,4) = 1
         end select
      else
         select case (lM%eNoN)
         case(3)
            lM%eOrd(1,1) = 1
            lM%eOrd(2,1) = 2
            lM%eOrd(1,2) = 2
            lM%eOrd(2,2) = 3
            lM%eOrd(1,3) = 3
            lM%eOrd(2,3) = 1
         case (4)  ! tet  !
            lM%eOrd(1,1) = 2
            lM%eOrd(2,1) = 1
            lM%eOrd(3,1) = 3
            lM%eOrd(1,2) = 1
            lM%eOrd(2,2) = 2
            lM%eOrd(3,2) = 4
            lM%eOrd(1,3) = 2
            lM%eOrd(2,3) = 3
            lM%eOrd(3,3) = 4
            lM%eOrd(1,4) = 3
            lM%eOrd(2,4) = 1
            lM%eOrd(3,4) = 4
         case (8)  ! brick/hex !
            lM%eOrd(1,1) = 1
            lM%eOrd(2,1) = 2
            lM%eOrd(3,1) = 6
            lM%eOrd(4,1) = 5
            lM%eOrd(1,2) = 2
            lM%eOrd(2,2) = 3
            lM%eOrd(3,2) = 7
            lM%eOrd(4,2) = 6
            lM%eOrd(1,3) = 3
            lM%eOrd(2,3) = 4
            lM%eOrd(3,3) = 8
            lM%eOrd(4,3) = 7
            lM%eOrd(1,4) = 1
            lM%eOrd(2,4) = 5
            lM%eOrd(3,4) = 8
            lM%eOrd(4,4) = 4
            lM%eOrd(1,5) = 1
            lM%eOrd(2,5) = 4
            lM%eOrd(3,5) = 3
            lM%eOrd(4,5) = 2
            lM%eOrd(1,6) = 5
            lM%eOrd(2,6) = 6
            lM%eOrd(3,6) = 7
            lM%eOrd(4,6) = 8
         end select
      end if

      return
      end subroutine selectel

!***********************************************************************

      subroutine VTK(lM)
      use varmod
      use vtkXMLMod
      implicit none

      type(meshType), intent(inout) :: lM
      integer :: iFa
      character(len=strL) :: fName
      logical :: flag

!     Write mesh vtu file
      write(fName,'(A)') "mesh-complete.mesh.vtu"
      lM%IEN(:,:) = lM%IEN(:,:) - 1
      call writeVTU(lM, fName)
      lM%IEN(:,:) = lM%IEN(:,:) + 1

!     Write face vtp files
      inquire(file="mesh-surfaces",exist=flag)
      if (.not.flag) call system("mkdir  mesh-surfaces")
      do iFa=1, lM%nFa
         write(fName,'(A)') "mesh-surfaces/"//TRIM(lM%fa(iFa)%fname)// &
            ".vtp"
         lM%fa(iFa)%IEN = lM%fa(iFa)%IEN - 1
         call writeVTP(lM%fa(iFa), fName)
         lM%fa(iFa)%IEN = lM%fa(iFa)%IEN + 1
      end do

      return
      end subroutine VTK

!***********************************************************************