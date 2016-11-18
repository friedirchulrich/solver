      subroutine char_bc(rho,ru,rv,rw,e,tc,pinput,ma,cvinput,ginput,rinp
     &                   ut,dx,dy,mmax,nx,nxe,kz)
c-----------------------------------------------------------------------
c     characteristic boundary condition at the free stream (Paul Harris) 
c     09/04/2002
c     Copyright Christian Stemmer
C...|....1....|....2....|....3....|....4....|....5....|....6....|....7..
      implicit none
	
      real rho(mmax,nx,kz),ru(mmax,nx,kz),rv(mmax,nx,kz),rw(mmax,nx,kz)
      real e(mmax,nx,kz),cv(mmax,nx,kz),rmix(mmax,nx,kz)
      real t(mmax-2:mmax,nx,kz),u(mmax-2:mmax,nx,kz)
      real v(mmax-2:mmax,nx,kz),w(mmax-2:mmax,nx,kz)
      real gamma(mmax,nx,kz)
      real ginput,rinput,pinput(mmax,nx,1),cvinput
      real th(nx),rh(nx),uh(nx),vh(nx),wh(nx),p(mmax,nx),tc(mmax,nx,kz)

      real ma,dx,dy,dyx
      real fma,fmu,fx,fx1,fx2,f1,f2,r1,r2,u1,u2,v1,v2,w1,w2,t1,t2

      integer my,my1,my2,nx,nxe,mmax,kz
      integer m,n,n1,n2,i,k,j

c
      do i=1,mmax
      do j=1,nx
      do n=1,kz
      gamma(i,j,k)=ginput
      rmix(i,j,k)=rinput
      cv(i,j,k)=cvinput
      end do
      end do
      end do
c     print *,' inside char_bc.f'
c     print *,' mmax,nx,nxe,kz ',mmax,nx,nxe,kz
c     tan theta = dy/dx
      dyx=dy/dx
      my  = mmax
      my1 = mmax-1
      my2 = mmax-2
c
      n=1
      do 10 k=1,kz
      do 10 m=my2,my
        u(m,n,k)=ru(m,n,k)/rho(m,n,k)
        v(m,n,k)=rv(m,n,k)/rho(m,n,k)
        w(m,n,k)=rw(m,n,k)/rho(m,n,k)
        t(m,n,k)=(e(m,n,k)-0.5*(ru(m,n,k)*u(m,n,k)+rv(m,n,k)*v(m,n,k)+
     &             rw(m,n,k)*w(m,n,k)))/(rho(m,n,k)*cv(m,n,k))
  10  continue
C     write (0,*) '  e(my,1)= ',e(my,1),'  e(my-1,1)= ',e(my-1,1),
C    &        '  e(my-2,1)= ',e(my-2,1)
C     write (0,*) '  t(my,1)= ',t(my,1),'  t(my-1,1)= ',t(my-1,1),
C    &        '  t(my-2,1)= ',t(my-2,1)


      do 20 k=1,kz
      do 20 n=2,nxe
         
c --- lokale Machzahl (fma) und charakteristischer Winkel fmu)
c --- free-stream Ma number (ma) - my - upper boundary

         u(my2,n,k)=ru(my2,n,k)/rho(my2,n,k)
         u(my1,n,k)=ru(my1,n,k)/rho(my1,n,k)
         u(my ,n,k)=ru(my ,n,k)/rho(my ,n,k)
         v(my2,n,k)=rv(my2,n,k)/rho(my2,n,k)
         v(my1,n,k)=rv(my1,n,k)/rho(my1,n,k)
         v(my ,n,k)=rv(my ,n,k)/rho(my ,n,k)
         w(my2,n,k)=rw(my2,n,k)/rho(my2,n,k)
         w(my1,n,k)=rw(my1,n,k)/rho(my1,n,k)
         w(my ,n,k)=rw(my ,n,k)/rho(my ,n,k)
c     write (0,*) 'w = ',w(my,n,k),rw(my,n,k),rho(my,n,k),w(my1,n,k),
c    &                   w(my2,n,k),mmax
c     werner nimmt hier den Betrag - hilft aber wohl kaum
         t(my2,n,k)=( e(my2,n,k)
     &            -0.5*(ru(my2,n,k)*u(my2,n,k)+rv(my2,n,k)*v(my2,n,k)+
     &             rw(my2,n,k)*w(my2,n,k)))/rho(my2,n,k)/cv(my2,n,k)
         t(my1,n,k)=( e(my1,n,k)
     &            -0.5*(ru(my1,n,k)*u(my1,n,k)+rv(my1,n,k)*v(my1,n,k)+
     &             rw(my1,n,k)*w(my1,n,k)))/rho(my1,n,k)/cv(my1,n,k)
         t(my ,n,k)=( e(my,n,k)
     &            -0.5*(ru(my ,n,k)*u(my ,n,k)+rv(my ,n,k)*v(my ,n,k)+
     &             rw(my ,n,k)*w(my ,n,k)))/rho(my ,n,k)/cv(my ,n,k)
C        if (n.eq.2) then
C        print *,' nx=',n,' u=',u(mmax,n),' v=',v(mmax,n),
C    &           ' t=',t(mmax,n),' e=',e(mmax,n)
C        print *,' nx=',n,' u=',u(mmax-1,n),' v=',v(mmax-1,n),
C    &           ' t=',t(mmax-1,n),' e=',e(mmax-1,n)
C        print *,' nx=',n,' u=',u(mmax-2,n),' v=',v(mmax-2,n),
C    &           ' t=',t(mmax-2,n),' e=',e(mmax-2,n)
C        endif

C        write (0,1000) ' nx=',n,'  u=',u(mmax,n,1),'  v=',v(mmax,n,1),
C    &           '  t=',t(mmax,n,1),'  e=',e(mmax,n,1)
C        write (0,*) 'u,v,w,gamma,rmix,t  '
C        write (0,*) my,n,k
C        write (0,*) u(my,n,k)
C        write (0,*) v(my,n,k)
C        write (0,*) w(my,n,k)
C        write (0,*) gamma(my,n,k)
C        write (0,*) rmix(my,n,k)
C        write (0,*) t(my,n,k)
c        fma=ma * sqrt((u(my,n)**2+v(my,n)**2) / t(my,n))

         fma=sqrt((u(my,n,k)**2+v(my,n,k)**2+w(my,n,k)**2)/
     &             (gamma(my,n,k)*rmix(my,n,k)*t(my,n,k)))
     
         if (fma.le.1.0) then
            print *,' u=',u(mmax,n,1),' v=',v(mmax,n,1),
     &              ' nx=',n,'  z='
         end if
         fmu=asin(1/fma) + atan( v(mmax,n,k)/u(mmax,n,k) )
c        fmu=(u(mmax,n)*sqrt(fmu**2-1)-v(mmax,n))/
c    &       (u(mmax,n)+v(mmax,n)*sqrt(fmu**2-1))
C         print *,'nx,fma,fmu,fx  ',n,fma,fmu,fx
     
c
c --- Lage der charakteristischen Hilfspunkte

         fx = dyx / tan(fmu)
C         print *,'nx,fma,fmu,fx  ',n,fma,fmu,fx
         if (real(n)-fx*2. .lt. 1.) fx=0.
         if (fx .lt. 0.) fx=0.
         fx1=real(n) - fx
         fx2=real(n) - fx*2.
         n1=int(fx1)
         n2=int(fx2)
         f1=fx1 - real(n1)
         f2=fx2 - real(n2)

c        if (n.eq.9) then
c        write (0,*) '  n= ',n,' fx= ',fx,' fmu= ',fmu,
c    &           'n1= ',n1,'  n2= ',n2,'f1= ',f1,'  f2= ',f2
c        endif

c --- lineare Interpolation
         r1= rho(my1,n1,k) + f1*( rho(my1,n1+1,k) - rho(my1,n1,k) )
         u1= u  (my1,n1,k) + f1*( u  (my1,n1+1,k) - u  (my1,n1,k) )
         v1= v  (my1,n1,k) + f1*( v  (my1,n1+1,k) - v  (my1,n1,k) )
         w1= w  (my1,n1,k) + f1*( w  (my1,n1+1,k) - w  (my1,n1,k) )
         t1= t  (my1,n1,k) + f1*( t  (my1,n1+1,k) - t  (my1,n1,k) )

         r2= rho(my2,n2,k) + f2*( rho(my2,n2+1,k) - rho(my2,n2,k) )
         u2= u  (my2,n2,k) + f2*( u  (my2,n2+1,k) - u  (my2,n2,k) )
         v2= v  (my2,n2,k) + f2*( v  (my2,n2+1,k) - v  (my2,n2,k) )
         w2= w  (my2,n2,k) + f2*( w  (my2,n2+1,k) - w  (my2,n2,k) )
         t2= t  (my2,n2,k) + f2*( t  (my2,n2+1,k) - t  (my2,n2,k) )
c         
c --- charakteristische Bedingung
         rh  (n)=(4.*r1 - r2) / 3.
         uh  (n)=(4.*u1 - u2) / 3.
         vh  (n)=(4.*v1 - v2) / 3.
         wh  (n)=(4.*w1 - w2) / 3.
         th  (n)=(4.*t1 - t2) / 3.

C        if (n.ge.1.and.n.le.4) then
C          print *,'t1,t2 =',t1,t2,t(my1,n1),t(my1,n1+1),t(my1,n1),
C    &              t(my2,n2),t(my2,n2+1),t(my2,n2)
C        endif

         rho(my,n,k)=         rh(n)
         ru (my,n,k)=         rh(n)*uh(n)
         rv (my,n,k)=         rh(n)*vh(n)
         rw (my,n,k)=         rh(n)*wh(n)
         tc (my,n,k)=         th(n)
         e  (my,n,k)=cv(my,n,k)*rh(n)*th(n)+((ru(my,n,k)*uh(n)+
     &               rv(my,n,k)*vh(n)+rw(my,n,k)*wh(n))*0.5)
c        if (n.eq.9) then
c         write (0,*) 'rho,ru,rv,rw,wh,w1,w2=  ',rh(n),ru(my,n),
c    &              rv(my,n),rw(my,n,k),wh(n),w1,w2,w(my1,n1,k),
c    &              w(my1,n1+1,k),w(my1,n1,k)
c         write (0,*)'w2=  ',w2,w(my2,n2,k),w(my2,n2+1,k),rw(my2,n2,k),
c    &              rho(my2,n2,k)
c        endif
 20   continue

C     print *,'  e(my,1)= ',e(my,1),'  e(my-1,1)= ',e(my-1,1),
C    &        '  e(my-2,1)= ',e(my-2,1)
C     print *,'  t(my,1)= ',t(my,1),'  t(my-1,1)= ',t(my-1,1),
C    &        '  t(my-2,1)= ',t(my-2,1)

C     print *,'==========nach char-bc==================='
C     print *,'rho(mmax)= '
C     print "(5es20.12)",(rho(my,n),n=1,5)
C     print *,'ru(mmax)= '
C     print "(5es20.12)",(ru(my ,n),n=1,5)
C     print *,'rv(mmax)= '
C     print "(5es20.12)",(rv(my ,n),n=1,5)
C     print *,'v(mmax)= '
C     print "(5es20.12)",(rv(my ,n)/rho(my,n),n=1,5)
C     print *,'e(mmax)= '
C     print "(5es20.12)",(e (my ,n),n=1,5)
C     print *,'t(mmax)= '
C     print "(5es20.12)",(t (my ,n),n=1,5)
C     print *,'p(mmax)= '
C     print "(5es20.12)",(p (my ,n),n=1,5)


 1000 FORMAT(a4,2x,i4,3(a4,f13.5),(a4,f17.5))
      
C     stop 'char_bc'

      return
      end

