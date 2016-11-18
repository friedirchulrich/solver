      subroutine nin(mmax,nmax,f,y,y0,h,idir)
c---------------------------------------------------------------------
c                  numerische integration
c---------------------------------------------------------------------
c     mmax ... anzahl punkte
c     f   ... zu integrierende funktion
c     y   ... integralfunktion
c     y0  ... startwert
c     h   ... schrittweite
c     idir... integrationsrichtung   1: von unten nach oben
c                                   -1: von oben  nach unten
c---------------------------------------------------------------------
      implicit none

      integer mmax,nmax,idir,m,j,z
      real p1,p2,p3,h,y0
      real f(mmax,nmax)
      real y(mmax,nmax)

      p1=h/24.
      p2=13.*p1
      p3=h/3.
      if (idir.eq.1) then
        do z=1,nmax
          y(1,z)=y0
          y(2,z)=y(1,z)+p1*(9.*f(1,z)+19.*f(2,z)-5.*f(3,z)+f(4,z))
          j=mmax-1
          do 100 m=3,mmax-1
            y(m,z)=y(m-1,z)-p1*(f(m-2,z)+f(m+1,z))+p2*(f(m-1,z)+f(m,z))
 100      continue
          y(mmax,z)=y(j-1,z)+p3*(f(j-1,z)+4.*f(j,z)+f(mmax,z))
        enddo
      else
        do z=1,nmax
         y(mmax  ,z) = y0
         y(mmax-1,z) = y(1,z)+p1*(9.*f(1,z)+19.*f(2,z)-5.*f(3,z)+f(4,z))
         j=2
         do 200 m=mmax-2,2,-1
           y(m,z)=y(m+1,z)-p1*(f(m+2,z)+f(m-1,z))+p2*(f(m+1,z)+f(m,z))
 200     continue
         y(1,z)=y(j+1,z)+p3*(f(j+1,z)+4.*f(j,z)+f(1,z))
        enddo
      endif
      return
      end
