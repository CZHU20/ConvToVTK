!--------------------------------------------------------------------
!
!     GMSH data structures interfaced.
!
!--------------------------------------------------------------------

      module gmshMod
      
      type gmshNodeType
         integer :: numNodeBlocks 
         integer :: numNodes      
         integer :: minNodeTag    
         integer :: maxNodeTag             
         integer,allocatable :: numNodesInBlock(:)
      end type gmshNodeType

      end module gmshMod