      module global_node_displacement_orthotropic
      implicit none
!     
      integer,save:: kk=0
      real(8),save:: n_iter(19823)
        
      contains
!=============subtoutines==========================     
      subroutine n_init
      n_iter=0.d0
      end subroutine n_init
!=========================================================================                        
      end module global_node_displacement_orthotropic
      