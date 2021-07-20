      subroutine element_orthotropic(rhs,amatrix,props,coords,u,svars,
     &  pre_d)
      implicit none
!     props=[ E  mu  t  Gc ft lb]

      real(8):: rhs(12,1), amatrix(12,12), props(6), coords(2,4), u(12)
      real(8):: svars(4), pre_d
       
!     local varibles
      real(8):: guassp(2), gwight(2), d(3,3), jacb(2,2), inv_jacb(2,2)
      real(8):: det_jacb, b(3,8), uu(8,1), dd(4,1), nd(1,4), bd(2,4)
      real(8):: ppf(1,1), pf, omega, domega, ddomega
      real(8):: ru(8,1), rd(4,1), rr(12,1), kuu(8,8), kdd(4,4) 
      real(8):: kk(12,12), energy_plus
      real(8):: Gf, lb, QQ(12,12), dalpha, ddalpha, c0
      real(8):: rd1(4,1), kdd1(4,4),indicate
      integer:: numgus, i, j, indexq(12), k1, k2
!      
      real(8):: dalpha0, pf0, ddalpha0, c00
      real(8):: omega0, domega0, ddomega0
      real(8):: ph_energy, d_ph_energy
      
   
!     initialize varibles
      amatrix=0.d0; rhs=0.d0; d=0.d0; jacb=0.d0; inv_jacb=0.d0; b=0.d0;
      uu=0.d0;  dd=0.d0;  nd=0.d0;  bd=0.d0; ppf=0.d0
      pf=0.d0;  omega=0.d0;  domega=0.d0;  ddomega=0.d0;
      ru=0.d0; rd=0.d0; rr=0.d0; kuu=0.d0; kdd=0.d0;  kk=0.d0;
      energy_plus=0.d0; QQ=0.d0; dalpha=0.d0; ddalpha=0.d0; c0=0.d0
      rd1=0.d0; kdd1=0.d0;
      
      dalpha0=0.d0; pf0=0.d0; ddalpha0=0.d0; c00=0.d0;
      omega0=0.d0; domega0=0.d0; ddomega0=0.d0
      indicate=0.d0
      
      
      do i=1,4
          uu(2*i-1,1)=u(3*i-2)
          uu(2*i,1)=u(3*i-1)
          dd(i,1)=u(3*i)
      enddo
      
      Gf=props(4);  lb=props(6)
      
      numgus=2
      guassp(1)=-1.d0/dsqrt(3.d0); guassp(2)=1.d0/dsqrt(3.d0)
      gwight(1)=1.d0;  gwight(2)=1.d0
      
      do i=1,numgus
          do j=1,numgus      
              call jacob_matrix(jacb,inv_jacb,det_jacb,guassp(i),
     &                            guassp(j),coords)
              call b_matrix(nd,bd,b,inv_jacb,guassp(i),guassp(j))
              
              ppf=matmul(nd,dd)
              pf=ppf(1,1)
              
              
              call d_matrix(d,props)
          
     
              call gener_engeryp(energy_plus,b,uu,props)         
              if (svars(2*(i-1)+j)>energy_plus) then
                  energy_plus=svars(2*(i-1)+j)
              else
                  svars(2*(i-1)+j)=energy_plus
              endif

              
              call alpha_dd(dalpha,pf,ddalpha,c0)
              call degradfunc(omega,domega,ddomega,pf,props,c0)
                      
              pf0=0.d0
              call alpha_dd(dalpha0,pf0,ddalpha0,c00)
          call degradfunc(omega0,domega0,ddomega0,pf0,props,c00)
              if (energy_plus<-2.d0*Gf/c0/lb/domega0) then
                  energy_plus=-2.d0*Gf/c0/lb/domega0
              endif
     
             ph_energy=domega*energy_plus+Gf/c0*1.d0/lb*dalpha
             d_ph_energy=ddomega*energy_plus+Gf/(c0*lb)*ddalpha
             
              rd1=gwight(i)*gwight(j)*
     &            (ph_energy*transpose(nd)+Gf/c0*(2.d0*lb*
     &        matmul(matmul(transpose(bd),bd),dd)))*det_jacb*props(3)
              rd=rd+rd1
              
              kdd1=gwight(i)*gwight(j)*
     &      ((d_ph_energy)*
     &      matmul(transpose(nd),nd)+2.d0*lb*Gf/c0*
     &       matmul(transpose(bd),bd))*det_jacb*props(3)
              kdd=kdd+kdd1
     
     
             kuu=kuu+gwight(i)*gwight(j)*
     &       (omega*matmul(matmul(transpose(b),d),b))*det_jacb*props(3)
       
          enddo
      enddo
      ru=matmul(kuu,uu)
      rd=rd*pre_d
      
      do i=1,8
          do j=1,8
              kk(i,j)=kuu(i,j)
          enddo
      enddo
      
      do i=9,12
          do j=9,12
              kk(i,j)=kdd(i-8,j-8)
          enddo
      enddo
      
      do i=1,8
          rr(i,1)=ru(i,1)
      enddo
      
      do i=1,4
          rr(i+8,1)=rd(i,1)
      enddo
      
      data indexq/1,2,9,3,4,10,5,6,11,7,8,12/   
      do i=1,12
          QQ(indexq(i),i)=1.d0
      enddo
      
      rhs=-matmul(transpose(QQ),rr)
      amatrix=matmul(matmul(transpose(QQ),kk),QQ)
  
      end subroutine element_orthotropic
!================subroutines========================================================
!================traditional material propoty matrix================================
      subroutine d_matrix(d,props)
      implicit none
      real(8):: d(3,3), props(6)
      real(8):: e, mu, a1, a2
      
      d=0.d0
      e=props(1);  mu=props(2)
      a1=e/(1.d0-mu*mu)
      a2=(1.d0-mu)/2.d0
      
      d(1,1)=1.d0; d(1,2)=mu
      d(2,1)=d(1,2);   d(2,2)=1.d0
                                       d(3,3)=a2
      d=a1*d
      end subroutine d_matrix      
!==================shape function and its derivative with xi and eta======================    
      subroutine shapefuc(n,dn_xieta,xi,eta)
      
      implicit none      
      real(8):: n(4), dn_xieta(2,4),xi, eta

      n(1)=1.d0/4.d0*(1.d0-xi)*(1.d0-eta)
      n(2)=1.d0/4.d0*(1.d0+xi)*(1.d0-eta)
      n(3)=1.d0/4.d0*(1.d0+xi)*(1.d0+eta)
      n(4)=1.d0/4.d0*(1.d0-xi)*(1.d0+eta)
      
      dn_xieta(1,1)=-1.d0/4.d0*(1.d0-eta)
      dn_xieta(1,2)= 1.d0/4.d0*(1.d0-eta)
      dn_xieta(1,3)= 1.d0/4.d0*(1.d0+eta)
      dn_xieta(1,4)=-1.d0/4.d0*(1.d0+eta)
      
      dn_xieta(2,1)=-1.d0/4.d0*(1.d0-xi)
      dn_xieta(2,2)=-1.d0/4.d0*(1.d0+xi)
      dn_xieta(2,3)= 1.d0/4.d0*(1.d0+xi)
      dn_xieta(2,4)= 1.d0/4.d0*(1.d0-xi)
      end subroutine shapefuc
!======================jacob matrix=========================================================     
      subroutine jacob_matrix(jacb,inv_jacb,det_jacb,xi,eta,coords)
      implicit none
      real(8):: jacb(2,2), inv_jacb(2,2), det_jacb, xi, eta, coords(2,4)
      
!     local varibles
      real(8):: n(4), dn_xieta(2,4)
      integer:: i
      
      jacb=0.d0
      
      call shapefuc(n,dn_xieta,xi,eta)
      
      do i=1,4
          jacb(1,1)=jacb(1,1)+dn_xieta(1,i)*coords(1,i)
          jacb(1,2)=jacb(1,2)+dn_xieta(1,i)*coords(2,i)
          jacb(2,1)=jacb(2,1)+dn_xieta(2,i)*coords(1,i)
          jacb(2,2)=jacb(2,2)+dn_xieta(2,i)*coords(2,i)
      enddo
      
      det_jacb=jacb(1,1)*jacb(2,2)-jacb(1,2)*jacb(2,1)
      
      inv_jacb(1,1)=jacb(2,2);  inv_jacb(1,2)=-jacb(1,2)
      inv_jacb(2,1)=-jacb(2,1); inv_jacb(2,2)=jacb(1,1)
      inv_jacb=1.d0/det_jacb*inv_jacb
      end subroutine jacob_matrix
!===============traditional b matrix==============================================      
      subroutine b_matrix(nd,bd,b,inv_jacb,xi,eta)
      implicit none
      real(8):: nd(1,4), bd(2,4), b(3,8), inv_jacb(2,2), xi, eta
      
!     local varibles
      real(8):: n(4), dn_xieta(2,4), dn_x(4), dn_y(4)
      integer:: i, j
      
!     initialize varibles
      b=0.d0
      
      call shapefuc(n,dn_xieta,xi,eta)
      
      do i=1,4
          nd(1,i)=n(i)
      enddo
      
      do i=1,4
          dn_x(i)=inv_jacb(1,1)*dn_xieta(1,i)
     &             +inv_jacb(1,2)*dn_xieta(2,i)
          dn_y(i)=inv_jacb(2,1)*dn_xieta(1,i)
     &             +inv_jacb(2,2)*dn_xieta(2,i)
      enddo
      
      do j=1,4
          b(1,2*(j-1)+1)=dn_x(j)
                                   b(2,2*(j-1)+2)=dn_y(j)
          b(3,2*(j-1)+1)=dn_y(j);  b(3,2*(j-1)+2)=dn_x(j)
      enddo
      
      do j=1,4
          bd(1,j)=dn_x(j)
          bd(2,j)=dn_y(j)
      enddo
      
      end subroutine b_matrix    
!===============degradation function and its derivations========================
      subroutine degradfunc(omega,domega,ddomega,pf,props,c0)
      implicit none
      real(8):: omega, domega, ddomega, pf, props(6), c0
!     local varibles
      real(8):: g, dg, ddg, fai, dfai, ddfai, E, Gf, ft, lb
      real(8):: p, a1, a2, a3

      E=props(1); Gf=props(4); ft=props(5); lb=props(6)
      a1=4.d0/(c0*lb)*E*Gf/(ft*ft)
      a2=1.3868d0;  p=2.d0
      a3=0.6567d0;
      
      g=(1.d0-pf)**p; dg=-p*(1.d0-pf)**(p-1.d0); 
      ddg=p*(p-1.d0)*(1.d0-pf)**(p-2.d0)
      
      fai=g+a1*pf+a1*a2*pf**2.d0+a1*a2*a3*pf**3.d0
      dfai=dg+a1+2.d0*a1*a2*pf+3.d0*a1*a2*a3*pf**2.d0
      ddfai=ddg+2.d0*a1*a2+6.d0*a1*a2*a3*pf
      
      omega=g/fai
      
      domega=(dg*fai-g*dfai)/(fai**2.d0)
      ddomega=((ddg*fai-g*ddfai)*fai-2.d0*(dg*fai-g*dfai)*dfai)
     &        /(fai**3.d0)

      end subroutine degradfunc
!=====================r alpha===================================================
      subroutine alpha_dd(dalpha,pf,ddalpha,c0)
      implicit none
      real(8):: dalpha, pf, ddalpha, c0
      
      c0=3.1415926535897932384626433832795028841971693993751d0
      dalpha=2.d0-2.d0*pf
      ddalpha=-2.d0
      end subroutine alpha_dd     
!==================energy_plus==================================================
      subroutine gener_engeryp(energy_plus,b,uu,props)
      implicit none
      real(8):: energy_plus, b(3,8), uu(8,1), props(6)
      real(8):: epsi(3,1), strain(2,2)
      real(8):: E, mu, lamd, G
      real(8):: valu(2), en1, pn1, pn2
      
      E=props(1);   mu=props(2)
      lamd=E*mu/(1.d0-mu**2.d0); G=E/(2.d0*(1.d0+mu))
      
      epsi=matmul(b,uu)
      strain(1,1)=epsi(1,1);   strain(1,2)=0.5d0*epsi(3,1)
      strain(2,1)=strain(1,2); strain(2,2)=epsi(2,1);
      
      call eig(strain,valu)
      
      en1=valu(1)+valu(2)
      
      pn1=0.5d0*(dabs(valu(1))+valu(1))
      pn2=0.5d0*(dabs(valu(2))+valu(2))     
      energy_plus=lamd/2.d0*(0.5d0*(dabs(en1)+en1))**2.d0
     &                +G*(pn1*pn1+pn2*pn2)
      end subroutine gener_engeryp                 
!==================phase decomposite============================================
      subroutine eig(a,valu)
      implicit none
      real(8):: a(2,2),ca(2,2),valu(2),eigv(2,2)
      real(8):: r(2,2),eigg(2,2)
      real(8):: tol,t,coss,sins,eps
      integer:: i,j
      
      eps=1.d-15
      
      eigv=0.d0
      eigv(1,1)=1.d0; eigv(2,2)=1.d0
      do i=1,2
          do j=1,2
              ca(i,j)=a(i,j)
          enddo
      enddo
      
      do while (dabs(ca(1,2))>eps .or. dabs(ca(2,1))>eps)
          tol=(ca(1,1)-ca(2,2))/(2.d0*ca(1,2))
          t=tol/dabs(tol)/(dabs(tol)+dsqrt(1.d0+tol**2.d0))
          coss=(1.d0+t**2.d0)**(-0.5d0)
          sins=t*coss
          r(1,1)=coss; r(1,2)=-sins
          r(2,1)=sins; r(2,2)=coss
          ca=matmul(matmul(transpose(r),ca),r)
          eigv=matmul(eigv,r)
      enddo
      
      if (ca(1,1)>ca(2,2)) then
          valu(1)=ca(1,1)
          valu(2)=ca(2,2)
      else
          valu(1)=ca(2,2)
          valu(2)=ca(1,1)
          eigg(1,1)=eigv(1,2);  eigg(1,2)=eigv(1,1)
          eigg(2,1)=eigv(2,2);  eigg(2,2)=eigv(2,1)
          do i=1,2
              do j=1,2
                  eigv(i,j)=eigg(i,j)
              enddo
          enddo
      endif

      end subroutine eig
!==================end of the subroutines=======================================