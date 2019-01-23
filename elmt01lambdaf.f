      subroutine elmt01(d,ul,xl,ix,tl,s,r,ndf,ndm,nst,isw)
      
      IMPLICIT NONE
      
      ! CZE BY SHR
      
      !---  use common blocks
      include 'bdata.h'
      include 'cdata.h'
      include 'counts.h'
      include 'eldata.h'
      include 'elplot.h'
      include 'eltran.h'
      include 'hdata.h'
      include 'iofile.h'
      include 'prstrs.h'
      include 'tdata.h'
      include 'pointer.h'
      include 'comblk.h'
      include 'errchk.h'
      include 'rdata.h'
      include 'pdata6.h' 
      include 'strnum.h'
      include 'augdat.h'
      include 'qudshp.h'
      
      integer :: i,j,k
      integer :: ndf,ndm,nst,isw
      integer :: ix(4)
      integer :: href,hnum  
      integer :: nqp, iqp
      double precision :: s(nst,nst)
      double precision :: r(nst)
      double precision :: d(*)
      double precision :: ul(ndf,nen,6)
      double precision :: xl(ndm,nen) 
      double precision :: tl(1)
      double precision :: hist(2*4)
      double precision :: Xnode(ndm,nen)  ! same as 'xl'
      double precision :: Unode(ndf,nen)  ! same as 'ul'  
      double precision :: xcnode(ndm,nen)   
      logical :: flg_gnl   ! in case we implement nonlin CZ for large defo
      logical :: flg_tan
      double precision, allocatable :: S1(:),S2(:),T1(:),T2(:),Dg(:)

      
C	FEAP Abbruchbedingung für MATLAB Optimierung
      if (niter .gt. 45) then
        write(*,*) 'niter > 45, Rechnung abgebrochen'
        stop
      endif  
C	Ende      
      
      flg_gnl = .FALSE.       ! geom. nonlinearity
      flg_tan = .FALSE.       ! tangant stiffness
      
      hnum = 5 ! NUM OF HISTORY VAR PER GUASS POINT! gap(1), gap(2), tra(1), tra(2), d
      
      select case (isw)

      case(1)

         call dinput (d,14)
         if(errck)then
            write(*,*)'error1'
         else
            do i=1,14
            write(*,*)'INPUTs for CZE: i=',i,'d1=',d(i)
            end do
         endif

         write(*,*)'***HERE ARE YOUR INPUTS SUMMARY FOR THE CZE***'
         write(*,*)'beta',d(8)
         write(*,*)'K0',d(7)/d(6)
         write(*,*)'lam0',d(6)
C	Achtung hier lambda_f anstatt Gc
         write(*,*)'lamf',d(9)
         write(*,*)'Gc',0.5*d(7)*d(9)
C	Änderung Ende
         
         nqp = NINT(d(12))
         nh1 = nqp*hnum + 1 + 8  ! +1 we add the area AND current coordinate of the element
         write(*,*)'nqp = ',nqp
         write(*,*)'nst = ',nst
         
!          pause        

      case(3)      ! compute tangent and residual
         
         flg_tan = .TRUE.
      
         call SE_ResTan(
     &        ndm,       nen,
     &        ndf,       nst,
     &        xl,        ul,
     &        dt,        flg_tan,
     &        r,         s,
     &        d,         n,
     &        hr,        nh3,
     &        nh1,       nh2)

      case(4)
         
      case(6)
         
         flg_tan = .FALSE.
      
         call SE_ResTan(
     &        ndm,       nen,
     &        ndf,       nst,
     &        xl,        ul,
     &        dt,        flg_tan,
     &        r,         s,
     &        d,         n,
     &        hr,        nh3,
     &        nh1,       nh2)
     
      case(8)
        
      case(12)

         print*,'elemt',n
         write(*,*)'dt=',dt
         write(*,*)'ttim=',ttim

         call getdef(ndm,nen,ndf,xl,ul)   ! writes down the coordinates

         nqp = NINT(d(12))
         allocate( S1(nqp) )
         allocate( S2(nqp) )
         allocate( T1(nqp) )
         allocate( T2(nqp) )
         allocate( Dg(nqp) )
         
         do iqp = 1,nqp   
            S1(iqp) = hr(nh2+hnum*(iqp-1)+0)
            S2(iqp) = hr(nh2+hnum*(iqp-1)+1)
            T1(iqp) = hr(nh2+hnum*(iqp-1)+2)
            T2(iqp) = hr(nh2+hnum*(iqp-1)+3)   
            Dg(iqp) = hr(nh2+hnum*(iqp-1)+4)
         end do

         open (100,file = 'S1')
         call printvalue(100,'e10.3',S1,1,nqp)
         open (200,file = 'S2')
         call printvalue(200,'e10.3',S2,1,nqp)            
         open (300,file = 'T1')
         call printvalue(300,'e10.3',T1,1,nqp) 
         open (400,file = 'T2')
         call printvalue(400,'e10.3',T2,1,nqp)
         open (500,file = 'Dg')
         call printvalue(500,'e10.3',Dg,1,nqp)
         
         open (601,file = '2DArea')
         call printvalue(601,'e10.3',hr(nh2+nqp*hnum),1,1)
         
         open (701,file = 'lam')
         call printvalue(701,'e10.3',
     &   hr(nh2+nqp*hnum+1:nh2+nqp*hnum+2),1,2)
         
         deallocate(S1)
         deallocate(S2)
         deallocate(T1)
         deallocate(T2)
         deallocate(Dg)
         
!          write(*,*)'NORMAL damage! at elemt = ',eln
        
      case default

      end select
      
      return
      
      end subroutine elmt01
      
      
      
!##################################################################
!
!
!     
!##################################################################

      subroutine SE_ResTan(
     &   ndm,       nen,
     &   ndf,       nst,
     &   Xnode,     Unode,
     &   dtt,        flg_tan,
     &   r,         s,
     &   d,         eln,
     &   hr,        nh3,
     &   nh1,       nh2)
     
      implicit none
      
      include 'tdata.h'

      integer, intent(in) :: ndm
      integer, intent(in) :: nen
      integer, intent(in) :: ndf
      integer, intent(in) :: nst 
      double precision :: Xnode(ndm,nen)  ! same as 'xl'
      double precision :: Unode(ndf,nen)  ! same as 'ul'
      integer :: nqp
      double precision, allocatable :: qdp(:,:)
      logical, intent(in) :: flg_tan
      double precision :: s(nst,nst)
      double precision :: r(nst)
      double precision :: d(*)
            
      double precision :: XIvec(1)
      double precision :: xcnode(ndm,nen)
      double precision :: Nvec(nen/2)
      double precision :: NDvec(nen/2,ndm-1)
      double precision :: Jac(ndm,ndm)
      double precision :: detJ
      double precision :: F(ndm,ndm)
      double precision :: H(ndm,ndm)
      double precision :: JacInv(ndm,ndm)
      double precision :: NDXvec(nen,ndm)
      double precision :: xxx(ndm,nen)
      
      integer :: eln 
      
      double precision :: hr(*)
      integer :: nh1,nh2,nh3 
      
      ! CZE PRAMETERS
      double precision :: UnodeV(nst)
      double precision :: res(nst)
      double precision :: Lmat(2*ndm,nst)
      double precision :: Nmat(2,4)
      double precision :: B(2,nst)
      double precision :: XMid(ndm,2)
      double precision :: MidLine
      double precision :: Rotc(nst,nst)
      double precision :: CT(2,2)
      double precision :: STRE(2)
      double precision :: gap(2), gap1(2), gap2(2), go2(2),stro(2)
      double precision :: dtt, rx
      integer :: iqp,i,j, hnum, href
      double precision :: lam1,lam2,lamf,beta
      double precision :: lamo1,lamo2
      double precision :: ndata(6)
      
      nqp = NINT(d(12))
      allocate( qdp(nqp,2) )
      
      hnum = 5
      r    = 0.0D0
      s    = 0.0D0
      beta = d(8)

      xcnode = Xnode + Unode

      call VecX(ndm,nst,Unode,UnodeV)

      call SetLmat(ndm,nst,Lmat)
      
      call SetXMid(xcnode,XMid)
      
      call SetRotMat(XMid,MidLine,Rotc)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! this part is for updateing the rotation matrix
      xxx(1,1) = hr(nh1+nqp*hnum+1) 
      xxx(2,1) = hr(nh1+nqp*hnum+2) 
      xxx(1,2) = hr(nh1+nqp*hnum+3) 
      xxx(2,2) = hr(nh1+nqp*hnum+4)
      xxx(1,3) = hr(nh1+nqp*hnum+5)
      xxx(2,3) = hr(nh1+nqp*hnum+6)
      xxx(1,4) = hr(nh1+nqp*hnum+7)
      xxx(2,4) = hr(nh1+nqp*hnum+8)
!       
      if (ttim.EQ.dt) then
         call SetXMid(Xnode,XMid)
      else
         call SetXMid(xxx,XMid)
      endif
      
      call SetRotMat(XMid,MidLine,Rotc)
!       
      hr(nh2+nqp*hnum+1:nh2+nqp*hnum+8) 
     &   = [xcnode(1,1),xcnode(2,1),xcnode(1,2),xcnode(2,2),
     &      xcnode(1,3),xcnode(2,3),xcnode(1,4),xcnode(2,4)]
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! this part is for updateing the rotation matrix

      hr(nh2+nqp*hnum) = MidLine

      call SetQuadPt(nqp,qdp)            
      
      do iqp = 1,nqp

         XIvec = qdp(iqp,1)               
         
         call SetNvec(ndm,nen,XIvec,Nvec,Nmat)
         
         rx  = DOT_PRODUCT( XMid(1,:), Nvec )
         
         B   = 0.0D0
         gap = 0.0D0
         
         B = MATMUL(Nmat,Lmat)
         
         gap = MATMUL(B,MATMUL(Rotc,UnodeV))

         call mate(d,dt,iqp,eln,gap,STRE,CT,hr,nh1,nh2,ndata)
 
         if (flg_tan) then

            s = s +
     &          MATMUL
     &          (
     &          MATMUL( TRANSPOSE( B ), CT ),
     &          B
     &          ) * MidLine/2.0D0 * qdp(iqp,2)!*rx*2*3.1415               ! in case of axisymm BC
                       
         end if

         r = r -
     &       MATMUL( TRANSPOSE( B ), STRE ) * MidLine/2.0D0 * qdp(iqp,2)
!     &       *rx*2*3.1415                                                ! in case of axisymm BC
         

      end do
      
      s = MATMUL(TRANSPOSE(Rotc),MATMUL(s,Rotc))
      
      r = MATMUL(TRANSPOSE(Rotc),r)
      
      deallocate(qdp)
      
      return
      
      end subroutine SE_ResTan
!     
!##################################################################      
!     
!##################################################################
!       
      subroutine SetQuadPt(nqp,qdp)
      
      implicit none
      
      integer :: nqp
      double precision :: qdp(nqp,2)
      double precision :: psq13, msq13
      double precision :: psq35, msq35 
      double precision, parameter :: one = 1.0D0
      integer :: i,j   ! FOR FURTHER GENERALIZATION

      psq13 = DSQRT(1.0D0/3.0D0)
      msq13 = - psq13
      psq35 = DSQRT(3.0D0/5.0D0)
      msq35 = - psq35
      qdp   = 0.0D0     
      
      if (nqp.EQ.1) then       ! one point integration
      
         qdp(1,:) = [0.0D0,2.0D0]
         
      elseif (nqp.EQ.2) then   ! two point gauß integration / ! Newton-Cotes integration. 
         
         qdp(1,:) = [msq13,one]
         qdp(2,:) = [psq13,one]
         
!          qdp(1,:) = [-one,one]
!          qdp(2,:) = [one,one]
      
      elseif (nqp.EQ.3) then   ! three point gauß integration
         
!          qdp(1,:) = [msq35,5.0D0/9.0D0]
!          qdp(2,:) = [0.0D0,8.0D0/9.0D0]
!          qdp(3,:) = [psq35,5.0D0/9.0D0]
         
         qdp(1,:) = [-one ,1.0D0/3.0D0]
         qdp(2,:) = [0.0D0,4.0D0/3.0D0]
         qdp(3,:) = [one  ,1.0D0/3.0D0]
      
      elseif (nqp.EQ.4) then        
         
         qdp(1,:) = [-DSQRT(3.0D0/7.0D0+DSQRT(6.0D0/5.0D0)*2.0D0/7.0D0),
     &               (18.0D0-DSQRT(30.0D0))/36.0D0]
         qdp(2,:) = [-DSQRT(3.0D0/7.0D0-DSQRT(6.0D0/5.0D0)*2.0D0/7.0D0),
     &               (18.0D0+DSQRT(30.0D0))/36.0D0]
         qdp(3,:) = [+DSQRT(3.0D0/7.0D0-DSQRT(6.0D0/5.0D0)*2.0D0/7.0D0),
     &               (18.0D0+DSQRT(30.0D0))/36.0D0]
         qdp(4,:) = [+DSQRT(3.0D0/7.0D0+DSQRT(6.0D0/5.0D0)*2.0D0/7.0D0),
     &               (18.0D0-DSQRT(30.0D0))/36.0D0]
     
      elseif (nqp.EQ.5) then        
         
         qdp(1,:) = [-(1.0D0/3.0D0)*DSQRT(5+2.0D0*DSQRT(10.0D0/7.0D0)),
     &               (322.0D0-13.0D0*DSQRT(70.0D0))/900.0D0]
         qdp(2,:) = [-(1.0D0/3.0D0)*DSQRT(5-2.0D0*DSQRT(10.0D0/7.0D0)),
     &               (322.0D0+13.0D0*DSQRT(70.0D0))/900.0D0]
         qdp(3,:) = [0.0D0,128.0D0/225.0D0]
         qdp(4,:) = [+(1.0D0/3.0D0)*DSQRT(5-2.0D0*DSQRT(10.0D0/7.0D0)),
     &               (322.0D0+13.0D0*DSQRT(70.0D0))/900.0D0]
         qdp(5,:) = [+(1.0D0/3.0D0)*DSQRT(5+2.0D0*DSQRT(10.0D0/7.0D0)),
     &               (322.0D0-13.0D0*DSQRT(70.0D0))/900.0D0]
         
      endif

      return

      end subroutine SetQuadPt
!     
!##################################################################      
!     
!##################################################################
!     
      subroutine SetNvec(ndm,nen,XIvec,Nvec,Nmat)
      
      implicit none

      integer, intent(in) :: ndm
      integer, intent(in) :: nen
      double precision, intent(in) :: XIvec(ndm-1)
      double precision, intent(out) :: Nmat(2,4)
      double precision :: Nvec(nen/2)
      double precision :: r1(2)
      double precision :: g1(2), g2(2)
      double precision :: h1(2)
      double precision :: xi
      double precision, parameter :: p = +1.0D0/2.0D0
      double precision, parameter :: m = -1.0D0/2.0D0
      double precision, parameter :: z = 0.0D0
      integer :: i,j
      data r1 /p, p/
      data g1 /m, p/

      xi = XIvec(1)
      
      Nvec = 0.0D0

      Nvec(:) = r1
     &     + g1 * xi 

      Nmat(1,:) = [Nvec(1),z,Nvec(2),z]
      Nmat(2,:) = [z,Nvec(1),z,Nvec(2)]

      return

      end subroutine SetNvec
!     
!##################################################################
!     
      subroutine SetNDvec(ndm,nen,XIvec,NDvec)
      
      implicit none

      integer, intent(in) :: ndm
      integer, intent(in) :: nen
      double precision, intent(in) :: XIvec(ndm-1)
      double precision, intent(out) :: NDvec(nen/2,ndm-1)     
      double precision :: g1(2), g2(2)
      double precision :: h1(2)
      double precision :: xi
      double precision, parameter :: p = +1.0D0/2.0D0
      double precision, parameter :: m = -1.0D0/2.0D0
      data g1 /m, p/

      xi = XIvec(1)

      NDvec = 0.0D0

      NDvec(:,1) = g1 

      return

      end subroutine SetNDvec
! 
!##################################################################
!
!##################################################################
!     
      subroutine mate(d,dt,iqp,eln,gap,STRE,CT,hr,nh1,nh2,ndata)

      implicit none

      double precision, intent(in) :: d(*)
      double precision, intent(in) :: gap(2)
      double precision, intent(out) :: CT(2,2)
      double precision, intent(out) :: STRE(2)
      double precision :: W(4)   
      double precision :: DELn
      double precision :: DELt
      double precision :: SIGmn
      double precision :: SIGmt
      double precision :: hr(*)
      integer :: sepcase,href,hnum,eln,nh1,nh2,nh3
      double precision :: gapd
      double precision :: gapo(2)
      double precision :: stro(2)
      double precision :: gapU(2)
      double precision :: strU(2)
      double precision :: ULF(2)
      double precision :: contTOL 
      double precision :: PenValn 
      double precision :: PenValt
      double precision :: PenValS 
      double precision :: nu
      double precision :: dt
      double precision :: dmg
      double precision :: kmax
      double precision :: k0
      double precision :: ndata(6)
      integer :: i,j, iqp, nqp, nn
      
      sepcase = d(1)
      DELn    = d(2)
      DELt    = d(3)
      SIGmn   = d(4)
      SIGmt   = d(5)
      W(1)    = d(6)  
      W(2)    = d(7)  
      W(3)    = d(8)  
      W(4)    = d(9) 	!jetzt lambda_f statt Gc
      PenValn = d(10)
      nu      = d(11)
      nqp     = NINT(d(12))
      
      k0   = W(2)/W(1) 
      CT   = 0.0D0
      STRE = 0.0D0
      hnum = 5

      href = nh1+hnum*(iqp-1)-1
      gapo(1:2) = hr(href+1:href+2)
      stro(1:2) = hr(href+3:href+4)
      dmg = hr(href+5)
      
!       gapo(1:2) = ndata(3*(iqp-1)+1:3*(iqp-1)+2)
!       dmg = ndata(3*(iqp-1)+3)
!       write(*,*)'***gapo(1) = ',gapo(1)
!       write(*,*)'***gapo(2) = ',gapo(2)
!       write(*,*)'***dmg = ',dmg
!       pause

      contTOL = 0.0D0
      
      write(*,*)'gap2 = ',gap(2)
      
      if (gap(2).LT.contTOL) then

         call TangSeprLaw(eln,d,dt,gap,gapo,dmg,CT,STRE)
	     call NormSeprLaw(eln,d,dt,gap,gapo,dmg,CT,STRE)  ! to make it more smooth and more stable 
         
          write(*,*)'contact-------------- '
!          pause
!        nn = 2    ! set it as d input
!         
!         if (MOD(nn,2).EQ.1) then
!            
!            CT(2,2) = CT(2,2) + PenValn*nn*(gap(2)**(nn-1))  ! for odd nn
!                
!         else
!            
!            CT(2,2) = CT(2,2) - PenValn*nn*(gap(2)**(nn-1))  ! for even nn: 2
!            
!         endif
!         
!         STRE(2) = STRE(2) + gap(2)*CT(2,2)/nn

          CT(2,2) = PenValn
          CT(1,2) = 0
          CT(2,1) = 0
          STRE(2) = PenValn*gap(2)
          
                
      else      
      
         call NormSeprLaw(eln,d,dt,gap,gapo,dmg,CT,STRE)
         call TangSeprLaw(eln,d,dt,gap,gapo,dmg,CT,STRE)    
                  
      endif

      href = nh2+hnum*(iqp-1)-1     
      hr(href+1:href+hnum) 
     &   = [gap(1),gap(2),STRE(1),STRE(2),dmg]
      
      return

      end subroutine mate

!##################################################################
!           
!##################################################################
!     
      subroutine NormSeprLaw(eln,d,dt,gap,gapo,dmg,CT,STRE)

      implicit none

      double precision, intent(in) :: d(*)
      double precision, intent(in) :: gap(2),gapo(2)
      double precision, intent(out) :: CT(2,2)
      double precision, intent(out) :: STRE(2)
      double precision :: W(4)   
      double precision :: DELn
      double precision :: DELt
      double precision :: SIGmn
      double precision :: SIGmt
      double precision :: n,t
      double precision :: contTOL 
      double precision :: PenValn 
      double precision :: PenValt
      double precision :: nu
      double precision :: dt
      double precision :: lam,lam0,lamf,k0,beta,dmg,Ndmg,dddl,Gc
      double precision :: u(3),v(3)
      double precision :: Am(2,2),Bm(2,2),Cm(2,2)
      double precision :: FL,dF,ep,si,r0,RS,FS
      logical :: dl
      double precision :: alp,lfo,Gr,gam
      double precision :: PI = 3.1415927
      integer :: sepcase,eln
     
      sepcase = d(1)
      si      = d(2)
      ep      = d(3)
      Gr      = d(4)
      gam     = d(5)
      W(1)    = d(6)  
      W(2)    = d(7)  
      W(3)    = d(8)  
      W(4)    = d(9)	!jetzt lambda_f statt Gc
      PenValn = d(10)
      nu      = d(11)
      contTOL = d(12)
      
      lam  = DSQRT(((gap(2)+DABS(gap(2)))/2.0D0)**2 + (beta*gap(1))**2)

      select case (sepcase)

      
      case(1)        !  *********   L-J potential  *********
         
         
         STRE(2) = 0.0D0
         CT(2,2) = 0.0D0
         CT(2,1) = 0.0D0
!          pause
         
         
      case(2)        !  *********   Exponential TRAC MODEL -- Reza-paper! based on Needleman  *********
      
      case(3)        !  *********   NEW TRAC MODEL baaed on wrrigers paggi 2011  *********  
      
      case(4)        !  *********   NEW BI-LINEAR MODEL based on Samimi 2009  *********

C	hier Änderung W(4) jetzt lamf anstatt Gc
         beta = W(3)
         k0   = W(2)/W(1)   
         lam0 = W(1)
!         lamf = 2*W(4)/W(2)
	 lamf = W(4)
	 Gc   = 0.5*w(2)*w(4)
C	Änderung Ende

	 
         if (lam.LT.lam0) then
            Ndmg = 0.0D0
         elseif (lam.LT.lamf) then
         
            Ndmg = (lamf*(lam-lam0))/(lam*(lamf-lam0))
         
         if (Ndmg.GT.dmg) then 
            CT(2,2) = 
     &   -k0*((lamf*lam0)/((lamf-lam0)*lam**3))*gap(2)**2
            CT(2,1) = 
     &   -k0*((lamf*lam0)/((lamf-lam0)*lam**3))*(beta**2)*gap(2)*gap(1)
!             write(*,*)'NORMAL damage! at elemt = ',eln           
         endif
            
         else
            Ndmg = 1.0D0
         endif        
         
!          write(*,*)'Ndmg = ',Ndmg
!          write(*,*)'dmg = ',dmg
!          pause

              
         if (Ndmg.GT.dmg) then
            dmg = Ndmg
            nu  = nu*(1-dmg)**(1.0D0)
            STRE(2) = k0*(1.0D0-dmg)*gap(2)  
     &              + nu*(gap(2)-gapo(2))/dt
            CT(2,2) = k0*(1.0D0-dmg)
     &              + nu/dt
     &              + CT(2,2)
!             CT(2,1) = 0.d0
         else
!             dmg = Ndmg
            nu = nu*(1-dmg)**(1.0D0)
            STRE(2) = k0*(1.0D0-dmg)*gap(2)
     &              + nu*(gap(2)-gapo(2))/dt       
            CT(2,2) = k0*(1.0D0-dmg)
     &              + nu/dt
     &              + CT(2,2)
!             CT(2,1) = 0.d0
         endif
  
      case(5)       
      
      
      case(6)        !  *********   NEW FORMULATION SD TO LD  *********        
        
      case default        

         write(*,*)'sepcase = default'
         pause

      end select
           
      return

      end subroutine NormSeprLaw
!     
!##################################################################
!           
!##################################################################
!
      subroutine TangSeprLaw(eln,d,dt,gap,gapo,dmg,CT,STRE)

      implicit none

      double precision, intent(in) :: d(*)
      double precision, intent(in) :: gap(2),gapo(2)
      double precision, intent(out) :: CT(2,2)
      double precision, intent(out) :: STRE(2)     
      double precision :: W(4)   ! WEIGHT FACTORs
      double precision :: DELn
      double precision :: DELt
      double precision :: SIGmn
      double precision :: SIGmt      
      double precision :: n,t
      double precision :: contTOL 
      double precision :: PenValn 
      double precision :: PenValt
      double precision :: nu
      double precision :: dt
      double precision :: lam, lam0, lamf, k0, beta, dmg, Ndmg,dddl,Gc     
      double precision :: u(3),v(3)
      double precision :: Am(2,2),Bm(2,2),Cm(2,2)
      double precision :: A,b
      double precision :: PI = 3.1415927
      logical :: dl
      double precision :: alp,lfo,Gr,gam
      integer :: sepcase,eln
            
      sepcase = d(1)
      A       = d(2)
      b       = d(3)
      Gr      = d(4)
      gam     = d(5)
      W(1)    = d(6)  
      W(2)    = d(7)  
      W(3)    = d(8)  
      W(4)    = d(9)	!jetzt lambda_f statt Gc 
      PenValn = d(10)
      nu      = d(11)
      contTOL = d(12)
C	hier Änderung W(4) jetzt lamf anstatt Gc
         beta = W(3)
         k0   = W(2)/W(1)   
         lam0 = W(1)
!         lamf = 2*W(4)/W(2)
	 lamf = W(4)
	 Gc   = 0.5*w(2)*w(4)
C	Änderung Ende
      
      lam  = DSQRT(((gap(2)+DABS(gap(2)))/2.0D0)**2 + (beta*gap(1))**2)
      
      select case (sepcase)

      case(1)        !  *********   L-J potential TRAC MODEL  *********
      
      
         STRE(1) = A*DSIN(2*PI*gap(1)/b)
         CT(1,1) = A*2*PI/b*DCOS(2*PI*gap(1)/b)
         CT(1,2) = 0.0D0
         
!          write(*,*)'STRE(1)',STRE(1)
!          write(*,*)'gap(1)',gap(1)
        
         
         
      case(2)        !  *********   Exponential TRAC MODEL -- Reza-paper! based on Needleman  *********    
      
      case(3)        !  *********   NEW TRAC MODEL baaed on wrrigers paggi 2011  ********* 
                  
      case(4)        !  *********   NEW BI-LINEAR MODEL based on Samimi 2009  *********
         
C	hier Änderung W(4) jetzt lamf anstatt Gc
         beta = W(3)
         k0   = W(2)/W(1)   
         lam0 = W(1)
!         lamf = 2*W(4)/W(2)
	 lamf = W(4)
	 Gc   = 0.5*w(2)*w(4)
C	Änderung Ende
         
         if (lam.LT.lam0) then
            Ndmg = 0.0D0
         elseif (lam.LT.lamf) then
         
            Ndmg = (lamf*(lam-lam0))/(lam*(lamf-lam0))
         
            if (Ndmg.GT.dmg) then 
            CT(1,1) =  
     &   -k0*((lamf*lam0)/((lamf-lam0)*lam**3))*(beta**4)*(gap(1)**2)
            CT(1,2) = 
     &   -k0*((lamf*lam0)/((lamf-lam0)*lam**3))*(beta**2)*gap(2)*gap(1)
!             write(*,*)'TANGENTIAL damage! at elemt = ',eln
            endif
     
         else
            Ndmg = 1.0D0
         endif
         
         if (Ndmg.GT.dmg) then
            dmg = Ndmg
            nu  = nu*(1-dmg)**(1.0D0)
            STRE(1) = k0*(1.0D0-dmg)*(beta**2)*gap(1)  
     &              + nu*(gap(1)-gapo(1))/dt
            CT(1,1) = k0*(1.0D0-dmg)*beta**2
     &              + nu/dt
     &              + CT(1,1)
!             CT(1,2) = 0.d0
         else
!             dmg = Ndmg
            nu = nu*(1-dmg)**(1.0D0)
            STRE(1) = k0*(1.0D0-dmg)*(beta**2)*gap(1)
     &              + nu*(gap(1)-gapo(1))/dt       
            CT(1,1) = k0*(1.0D0-dmg)*beta**2
     &              + nu/dt
     &              + CT(1,1)
!             CT(1,2) = 0.d0
         endif
      
      case(5)       !  *********   NEW TRAC MODEL baaed on wrrigers paggi 2011  *********
       
      case(6)       !  *********   NEW FORMULATION SD TO LD  *********      
        
      case default        !     nothing to do

         write(*,*)'sepcase = default'
         
         pause

      end select
           
      return

      end subroutine TangSeprLaw
!        
!##################################################################
 
!##################################################################

      subroutine SetLmat(ndm,nst,Lmat)
      
      implicit none
      
      integer, intent(in) :: nst
      integer, intent(in) :: ndm
      double precision, parameter :: p = +1.0D0
      double precision, parameter :: m = -1.0D0
      double precision, parameter :: z = 0.0D0
      double precision, intent(out) :: Lmat(2*ndm,nst)
      
      integer i
        
      Lmat = 0.0D0
      
      Lmat(1,:) = [m,z,z,z,z,z,p,z]
      Lmat(2,:) = [z,m,z,z,z,z,z,p]
      Lmat(3,:) = [z,z,m,z,p,z,z,z]
      Lmat(4,:) = [z,z,z,m,z,p,z,z]
      
      return
      
      end
      
!##################################################################
!
!##################################################################

      subroutine VecX(ndm,nst,xcnode,xcnodeV)
      
      implicit none
      
      integer, intent(in) :: nst
      integer, intent(in) :: ndm
      double precision :: xcnode(ndm,4)
      double precision :: xcnodeV(nst)
      
      integer i,j
      
      do j = 1,nst/2
         do i = 1,ndm
            xcnodeV(i+ndm*(j-1)) = xcnode(i,j)
         enddo
      enddo
      
      return
      
      end
      
!##################################################################
!
!##################################################################
      
      subroutine SetXMid(xcnode,XMid)
      
      implicit none
      
      double precision :: xcnode(2,4)
      double precision :: XMid(2,2)
      
      integer :: i,j
      
      XMid(:,1) = (1.0D0/2.0D0)*(xcnode(:,1)+xcnode(:,4))
      XMid(:,2) = (1.0D0/2.0D0)*(xcnode(:,2)+xcnode(:,3))
      
      return
      
      end

!##################################################################
!
!##################################################################

      subroutine SetRotMat(XMid,MidLine,Rotc)
      
      implicit none
      
      double precision :: XMid(2,2)
      double precision :: MidLine
      double precision :: rot(2,2)
      double precision :: Teta
      double precision :: Rotc(8,8)
      
      integer :: i,j
      
      MidLine = NORM2(XMid(:,1)-XMid(:,2))
      
      Teta = DATAN((XMid(2,2)-XMid(2,1))/(XMid(1,2)-XMid(1,1)))
      
      if ((XMid(1,2)-XMid(1,1)).LT.0) then
      
         Teta = Teta + 4.D0*DATAN(1.D0)
         
      endif
      
      rot(1,:) = [+DCOS(Teta),+DSIN(Teta)] !  [1.0D0,0.0D0] !  [+DCOS(Teta),+DSIN(Teta)]
      rot(2,:) = [-DSIN(Teta),+DCOS(Teta)] !  [0.0D0,1.0D0] !  [-DSIN(Teta),+DCOS(Teta)]
      
!       print *, "rot-----------"
!       print *,Teta
!       DO j = 1,2
!       print *,(rot(i,j),i=1,2)
!       END DO
!       print *,""
!       print *,"/////////XD////////"
!       print *,XMid(1,2)-XMid(1,1)
!       pause
      
      Rotc = 0.0D0
      Do i = 1,4
         Rotc(2*i-1:2*i,2*i-1:2*i) = rot
      enddo
      
      return
      
      end
!
!##################################################################

!##################################################################
!                printm(0,'SHP_DX','e10.2',NDXvec,8,3)
      subroutine printm(un,name,fmt,ma,m,n)

      implicit none

      integer, intent(in) :: un
      integer, intent(in) :: m
      integer, intent(in) :: n
      character*(*), intent(in) :: name
      character*(*), intent(in) :: fmt
      double precision, intent(in) :: ma(m*n)

      integer :: lb
      integer :: ub
      integer :: i
      integer :: j
      character(len=2) nc
      character(len=100) ffmt

      write(nc,'(i2)') n

      ffmt = '(' // TRIM(ADJUSTL(nc)) // '(' // TRIM(fmt) // ',X))'
      
      write(un,*)

      write(un,*) name   

      do i=1,m

         write(un,TRIM(ffmt)) ( ma( (j-1)*m + i ), j=1, n )

      end do

      return

      end
!
!##################################################################
!     
!##################################################################
!              
      subroutine printvalue(un,fmt,ma,m,n)

      implicit none

      integer, intent(in) :: un
      integer, intent(in) :: m
      integer, intent(in) :: n
      character*(*), intent(in) :: fmt
      double precision, intent(in) :: ma(m*n)

      integer :: lb
      integer :: ub
      integer :: i
      integer :: j
      character(len=2) nc
      character(len=100) ffmt

      write(nc,'(i2)') n

      ffmt = '(' // TRIM(ADJUSTL(nc)) // '(' // TRIM(fmt) // ',X))'
      
!       write(un,*)  

      do i=1,m

         write(un,TRIM(ffmt)) ( ma( (j-1)*m + i ), j=1, n )

      end do

      return

      end
!
!##################################################################
      
      subroutine getdef(ndm,nen,ndf,Xnode,Unode)
      
      implicit none
      
      integer, intent(in) :: ndm
      integer, intent(in) :: nen
      integer, intent(in) :: ndf
      double precision :: Xnode(ndm,nen)  ! same as 'xl'
      double precision :: Unode(ndf,nen)  ! same as 'ul'
      double precision :: xcnode(ndm,nen)
      
      xcnode = Xnode + Unode
      open (700,file = 'xc')
      call printvalue(700,'e10.3',xcnode,ndm,nen)
      
      return
      
      end subroutine getdef 
