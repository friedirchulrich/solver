      program similarity
* similarity solution of compressible laminar boundary layer
*
* everything is dimensional
*
************************************************************************
* Tue Feb 27 17:14:12 MST 2001

c f77 -o similarity -r8 similarity.f
c
c f77 -o similarity -r8 -g -DEBUG:div_check=3:subscript_check=ON
c -DEBUG:trap_uninitialized=ON:verbose_runtime=ON similarity.f
c
c ftnchek similarity.f

      parameter (ny=10000) ! integration steps
      parameter (no=  10) ! output every no profile point
      dimension f(0:ny),fp(0:ny),fpp(0:ny),g(0:ny),gp(0:ny)

c length of boundary layer
c Re= 17,25E-6  Re=pe/te/rgas
      re_x     =  300790.

      ue       = 0.11495E+03 ! [m/s]5
      te       = 2.73E+02 ! [K]
      pe       = 101325.0 ! [N/m^2]
      re       = 1.29 ! [Kg/m^3]
c p=rho*R*T
c.....for air:
C     rgas     = 0.28702224E+03 ! [J/(kgK)]
      rgas     = pe/te/re
c.....for pure nitrogen:
c     rgas     = 0.28812511E+03 ! [J/(kgK)]
      gam      = 1.4
      prandtl  = 0.72

* maximum number of iterations, abortion criterion
      irelax   = 100
      tol      = 1.03e-06
* maximum eta
      etamax   = 20.

      rhoe     = pe/(rgas*te)
      vise     = computemu(1.,te,1.)
      print *,'vise= ',vise
      cp       = gam*rgas/(gam-1.)
      rhovise  = rhoe*vise
      rrhovise = 1./rhovise
      ue2rhe   = ue**2/(cp*te)
      reue     = rhoe*ue
      he       = rgas/(gam-1.)*te+0.5*ue**2+pe/rhoe
      xi       = re_x*vise**2
      deta     = etamax/float(ny)
      fac1     = sqrt(2.*xi)/ue*deta
      fac2     = reue/vise

* adiabatic wall flat plate BCs
      f (0)    = 0.
      fp(0)    = 0.
      gp(0)    = 0.
* initial guess (incompressible, laminar BL, adiabatic wall temperature condition)
      g(0)     = 1.02
      fpp(0)   =  0.47
C     fpp(0)   = 0.5
C     g(0)     = 1+sqrt(prandtl)*ue**2/(gam*rgas*te)*0.5*(gam-1.)
* deltas
      dfpp     = fpp(0)*1.e-5
      dg       =   g(0)*1.e-5
      do i=1,irelax
         call integrat(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &        rrhovise,te,pe,deta,ue2rhe,ny)
         fpe    = fp(ny)
         ge     =  g(ny)
         err=max(abs(fpe-1.),abs(ge-1.))
         print '(a,i5,a,e12.5)','it ',i,' err= ',err
         if (err.lt.tol) goto 1
C        if (i.eq.2) goto 1
         fpp(0) = fpp(0)+dfpp ! delta in fpp direction
         call integrat(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &        rrhovise,te,pe,deta,ue2rhe,ny)
         fpefp  = fp(ny)
         gefp   =  g(ny)
         fpp(0) = fpp(0)-dfpp
         g(0)   = g(0)+dg     ! delta in g direction
         call integrat(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &        rrhovise,te,pe,deta,ue2rhe,ny)
         fpegp  = fp(ny)
         gegp   =  g(ny)
         g(0)   = g(0)-dg
* new guess
         dfpdg   = (fpegp-fpe)/dg
         dfpdfpp = (fpefp-fpe)/dfpp
         dgdg    = ( gegp- ge)/dg
         dgdfpp  = ( gefp- ge)/dfpp
         rjac    = 1./(dfpdg*dgdfpp-dfpdfpp*dgdg)
         ddg     = rjac*( dgdfpp*(1.-fpe)-dfpdfpp*(1.-ge))
         ddfpp   = rjac*(-dgdg  *(1.-fpe)+dfpdg  *(1.-ge))
         fpp(0)  = fpp(0)+ddfpp
         g(0)    = g(0)  +ddg
         print *,'g(0),fpp(0)',g(0),fpp(0)
      enddo
      print '(a)','no convergence in simlarity'
      stop
 1    continue
      print '(a)','analytical BL solution to similarity.dat'
      open(99,file='Ma5_chris.dat',status='unknown',
     & form='formatted')
      print '(a)','output format: y,eta,p,rho,T,u,dudy'
      print '(a)','               1   2 3   4 5 6'
      write(99,'(A)') 'TITLE="BOUNDARY-LAYER"'
      write(99,'(A)') 'VARIABLES = x,y,z , Eta , rho , T, u, v'
      write(99,'(A)') 'ZONE T="T",I= 1,J=12,K=1'
      eta=0.
      y=0.
      adelta1 = 0.
      atheta  = 0.
      adelta3 = 0.
      addelta = 0.
      atw=g(0)*te
      arw=pe/(rgas*atw)
      hw=rgas/(gam-1.)*atw+pe/arw
      uulast=-1.e30
      do n=0,ny
         uu=fp(n)*ue
         if (abs((uu-uulast)/uu).le.tol) goto 2
         uulast=uu
         tt=g(n)*te
         rho=pe/(rgas*tt)
         if (mod(n,no).eq.0) then
          write(99,'(8(e12.5))')1.0,y,1.0,eta,rho,tt,uu
         endif
         eta=eta+deta
         dy=fac1/rho
         y=y+dy
         hh=rgas/(gam-1.)*tt+0.5*uu**2+pe/rho
         rruu=rho*uu
         adelta1 = adelta1 +           (1.-rruu/reue    )*dy
         atheta  = atheta  + rruu/reue*(1.-  uu/  ue    )*dy
         adelta3 = adelta3 + rruu/reue*(1.- (uu/  ue)**2)*dy
         addelta = addelta + rruu/reue*((hh-he)/(hw-he) )*dy
      enddo
 2    continue
      close(99)
      readelta1 = fac2*adelta1
      reatheta  = fac2*atheta
      readelta3 = fac2*adelta3
      readdelta = fac2*addelta
      avisw     = vise*computemu(atw,te,vise)
      acf=sqrt(2.)*arw*avisw/(rhovise*sqrt(re_x))*fpp(0)
      print '(a,e12.5)','T_w [K]   = ',atw
      print '(a,e12.5)','c_F       = ',acf
      print '(a,e12.5)','Re_delta1 = ',readelta1
      print '(a,e12.5)','Re_theta  = ',reatheta
      print '(a,e12.5)','Re_delta3 = ',readelta3
      print '(a,e12.5)','Re_Delta  = ',readdelta
      print '(a,e12.5)','H         = ',readelta1/reatheta
      print '(a,e12.5)','H_32      = ',readelta3/reatheta
      print '(a,e12.5)','deta      = ',deta

      end



      subroutine integrat(f,fp,fpp,g,gp,vise,rgas,prandtl,
     &     rrhovise,te,pe,deta,ue2rhe,ny)
      dimension f(0:ny),fp(0:ny),fpp(0:ny),g(0:ny),gp(0:ny)
      real p1,p2,p3

      p1=deta/24.
      p2=13.*p1
      p3=deta/3.

      t    = te*g(0)
      rho  = pe/(rgas*t)
      vis  = vise*computemu(t/te,te,vise)
      c    = rho*vis*rrhovise
      rc   = 1./c
      cold = c
      dcde = (c - cold)/deta
      rhs    = rc*(fpp(0)*dcde + f(0)*fpp(0))
      fpp(1) = fpp(0) - deta*rhs
      rhs    = rc*( gp(0)*dcde + prandtl*(f(0)*gp(0) + 
     &     c*ue2rhe*fpp(0)**2) )
      gp(1)  = gp(0)  - deta*rhs
      fp(1)=fp(0)+p1*(9.*fpp(0)+19.*fpp(1)-5.*fpp(2)+fpp(3))
      f (1)=f (0)+p1*(9.*fp (0)+19.*fp (1)-5.*fp (2)+fp (3))
      g (1)=g (0)+p1*(9.*gp (0)+19.*gp (1)-5.*gp (2)+gp (3))
      do n=2,ny-1
         t    = te*g(n-1)
         rho  = pe/(rgas*t)
         vis  = vise*computemu(t/te,te,vise)
         c    = rho*vis*rrhovise
         rc   = 1./c
         dcde = (c - cold)/deta
         cold = c
         rhs    = rc*(fpp(n-1)*dcde + f(n-1)*fpp(n-1))
         fpp(n) = fpp(n-1) - deta*rhs
         rhs    = rc*( gp(n-1)*dcde + prandtl*(f(n-1)*gp(n-1) + 
     &        c*ue2rhe*fpp(n-1)**2) )
         gp(n)  = gp(n-1)  - deta*rhs
         fp(n)  = fp(n-1) + deta*fpp(n)
         f (n)  = f (n-1) + deta*fp (n)
         g (n)  = g (n-1) + deta*gp (n)
      enddo
      t    = te*g(ny-1)
      rho  = pe/(rgas*t)
      vis  = vise*computemu(t/te,te,vise)
      c    = rho*vis*rrhovise
      rc   = 1./c
      dcde = (c - cold)/deta
      cold = c
      rhs    = rc*(fpp(ny-1)*dcde + f(ny-1)*fpp(ny-1))
      fpp(ny) = fpp(ny-1) - deta*rhs
      rhs    = rc*( gp(ny-1)*dcde + prandtl*(f(ny-1)*gp(ny-1) + 
     &     c*ue2rhe*fpp(ny-1)**2) )
      gp(ny)  = gp(ny-1)  - deta*rhs
      fp(ny)=fp(ny-2)+p3*(fpp(ny-2)+4.*fpp(ny-1)+fpp(ny))
      f (ny)=f (ny-2)+p3*(fp (ny-2)+4.*fp (ny-1)+fp (ny))
      g (ny)=g (ny-2)+p3*(gp (ny-2)+4.*gp (ny-1)+gp (ny))

      end



      function computemu(temp,temp_0,rmu_0)

      parameter (t1=40.,t2=110.4,c1=6.80689413e-8,c2=14.458e-7)

      t_dim=temp*temp_0
    
      if (t_dim.gt.t2) then
         computemu=c2*sqrt(t_dim)*t_dim/(t_dim+t2)
      elseif (t_dim.gt.t1) then
         computemu=c1*t_dim
      else
         computemu=c1*t1
      endif

      computemu=computemu/rmu_0

      end
