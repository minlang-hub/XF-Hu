! Source code provided by:
! Dr. Xiaofei Hu (Email:hxf@dlut.edu.cn); Dr. Weian Yao; Mr. Peng Zhang (Ph.D candidate)
! Dalian University of Technology
!
! The theory and example parameters are referred to:
! Zhang P, Hu XF, Wang XY, Yao WA. An iteration scheme for phase field model for
! cohesive fracture and its implementation in Abaqus. Engineering Fracture Mechanics, 2018, 204: 268-287.
!
      include 'element_orthotropic.for'
      include 'global_node_displacement_orthotropic.for'
!           
      subroutine UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
! LOP:    =,0,1,2,3
!			0:at the beginning of program
!			1:at the beginning of iteration step
!			2:at the end of iteration step
!			3: indicates that the subroutine is being called at the end of the analysis
      use global_node_displacement_orthotropic
      INCLUDE 'ABA_PARAM.INC'  
      
      DIMENSION TIME(2)

      if (LOP==2) then
      elseif (LOP==1 .or. LOP==0) then
          call n_init
      end if
!
      end subroutine UEXTERNALDB

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     & KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     & LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)

      use global_node_displacement_orthotropic
      INCLUDE 'ABA_PARAM.INC'
!     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     & DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)  
!
!          
      real(8):: coord(2,4)
      integer:: k1,k2
      real(8):: pf_s, pf
      
      real(8):: pre_d
        
      n_iter(jelem)=n_iter(jelem)+1.d0
      pre_d=(1.d0+(-1.d0)**n_iter(jelem))/2.d0
 
      do k1=1,2
          do k2=1,4
              coord(k1,k2)=coords(k1,k2)
          enddo
      enddo
      
      call element_orthotropic(RHS,AMATRX,PROPS,coord,U,SVARS, pre_d)
      

      END subroutine uel
