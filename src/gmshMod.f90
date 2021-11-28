!--------------------------------------------------------------------
!
!     GMSH data structures interfaced.
!
!--------------------------------------------------------------------

      module gmshMod

      type gmshPhysicalNameType
         integer :: num
         integer,allocatable :: dimension(:)
         integer,allocatable :: physicalTag(:)
         character(len=80),allocatable :: name(:)
      end type gmshPhysicalNameType

      ! type gmshEntityType
         
      ! end type gmshEntityType
      
      type gmshNodeType
         integer :: numNodeBlocks 
         integer :: numNodes      
         integer :: minNodeTag    
         integer :: maxNodeTag             
         integer,allocatable :: numNodesInBlock(:)
         real(kind=8),allocatable :: coord(:,:)
      end type gmshNodeType

      type gmshElementType
         integer :: numElementBlocks
         integer :: numElements
         integer :: minElementTag
         integer :: maxElementTag

         integer,allocatable :: ElementDim(:)
         integer,allocatable :: EntityTag(:)
         integer,allocatable :: eNoN(:)
         integer,allocatable :: numElementsInBlock(:)
         integer,allocatable :: conn(:,:)
      end type gmshElementType

      ! Variables
      type(gmshPhysicalNameType) :: gmshPhysicalNames
      type(gmshNodeType)         :: gmshNodes
      type(gmshElementType)      :: gmshElements

      end module gmshMod