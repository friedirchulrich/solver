      subroutine lagrange_intpol(y,dy,dymax,u,uint,t,tint,r,rint,
     &                           ru,ruint,mmax)
C...|....1....|....2....|....3....|....4....|....5....|....6....|....7..
C     interpolate between non-aequidistant points through lagrange 
C     polynomials 
C.......................................................................
      implicit none

      integer dymax,mmax,m_l,j
      integer index
      real y(dymax),u(dymax),t(dymax),ru(dymax),r(dymax)
      real uint(mmax),tint(mmax),rint(mmax),ruint(mmax)
      real y1,y2,y3,y4,y5,hh,dy,fak1,fak2,fak3,fak4,fak5

      uint (1)=u (1)
      tint (1)=t (1)
      rint (1)=r (1)
      ruint(1)=ru(1)

      hh=0.
      do m_l=2,mmax
        hh=hh+dy
        index=2
        do j=2,dymax
          if (y(j).le.hh) then
           index=j
          endif
        enddo

        if (index.eq.dymax) then
         uint(m_l)=u(dymax)
         tint(m_l)=t(dymax)
         rint(m_l)=r(dymax)
         ruint(m_l)=ru(dymax)
        else
c         if (abs(hh-y(index+1)).lt.abs(hh-y(index)))
c     &    index=index+1
         if (index.eq.(dymax-1)) index=dymax-2
         if (index.eq.2) index=3
C        print "(a,f9.5,4es12.5)",'hh,y(i)',hh,y(index-1),y(index),
C    &            y(index+1),y(index+2)
         y1=y(index-1)
         y2=y(index  )
         y3=y(index+1)
         y4=y(index+2)
         y5=y(index-2)
         fak1=(hh-y2)*(hh-y3)*(hh-y4)*(hh-y5)/((y1-y2)*(y1-y3)*(y1-y4)
     &         *(y1-y5))
         fak2=(hh-y1)*(hh-y3)*(hh-y4)*(hh-y5)/((y2-y1)*(y2-y3)*(y2-y4)
     &         *(y2-y5))
         fak3=(hh-y1)*(hh-y2)*(hh-y4)*(hh-y5)/((y3-y1)*(y3-y2)*(y3-y4)
     &         *(y3-y5))
         fak4=(hh-y1)*(hh-y2)*(hh-y3)*(hh-y5)/((y4-y1)*(y4-y2)*(y4-y3)
     &         *(y4-y5))
         fak5=(hh-y1)*(hh-y2)*(hh-y3)*(hh-y4)/((y5-y1)*(y5-y2)*(y5-y3)
     &         *(y5-y4))
         uint(m_l)=u(index-1)*fak1+u(index)*fak2+u(index+1)*fak3+
     &           u(index+2)*fak4+u(index-2)*fak5
         tint(m_l)=t(index-1)*fak1+t(index)*fak2+t(index+1)*fak3+
     &           t(index+2)*fak4+t(index-2)*fak5
         rint(m_l)=r(index-1)*fak1+r(index)*fak2+r(index+1)*fak3+
     &           r(index+2)*fak4+r(index-2)*fak5
         ruint(m_l)=ru(index-1)*fak1+ru(index)*fak2+ru(index+1)*fak3+
     &           ru(index+2)*fak4+ru(index-2)*fak5
        endif
      enddo

c     stop 'test'

      return
      end
