      program interpolate
      implicit real*8 (a-h,o-z)
      character file*80,name*18
      dimension tg(2)
      common/gridsyn/ xlam(4000),y1(4,4000,80,19),teff(80),grav(19),
     .                y11(4000),dy11(2,4000),nf,nt,ng,ic
      common/coeffs/ coeff(4,2,4000,80,19)

      open(9,file='input_interpolate',status='old')
c
c     Read the ng input files
c
      read(9,*) ng
      if(ng.lt.4) stop 'ng.lt.4 cannot spline'
      if(ng.gt.19) stop 'ng.gt.19'
      ic=1
      do ig=1,ng
         read(9,100) file          
         open(2,file=file,status='old')
         write(*,*) file
         read(2,*) nf
         if(nf.gt.4000) stop 'nf.gt.4000'
         read(2,'(10f8.2)') (xlam(ij),ij = 1,nf)
         do it = 1,1000
            if(it.gt.80) stop 'nt.gt.80'
            read(2,101,end=3) teff(it),grav(ig)
            read(2,'(6e12.5)') (y1(ic,ij,it,ig), ij = 1,nf)
         enddo
 3       nt=it-1
         if(nt.lt.4) stop 'nt.lt.4 cannot spline'
         grav(ig)=dlog10(grav(ig))
      enddo

c
c     Calculate spline coefficients
c
      write(*,*) ' Calling splineg'
      call splineg
      write(*,*) ' Calling splinet'
      call splinet

      open(50,file='list',status='old')
c      do i=1,29
c         read(50,*)
c      enddo
      
      do i=1,2000
         read(50,102,end=99) name,tg(1),tg(2)
         call getf(tg(1),tg(2))
c         if(name(9:9).eq.'') then
            write(file,61) name(1:6),name(7:12)
 61         format(a6,a6,'_SED')
c         else
c            write(file,62) name(1:4),name(10:14)
c 62         format('J',a4,a5,'_synth')
c         endif
         open(60,file=file)
         write(60,'(f6.0,2x,f4.2)') tg(1),tg(2)
         write(60,'(i5)') nf
         do j=1,nf
            write(60,'(f8.2,2x,1pe12.5)') xlam(j),y11(j)
         enddo
         close(60)
      enddo

 100  format(a)
 101  format(24x,f11.1,11x,e11.3)
 102  format(a18,16x,f6.0,6x,f4.2)
      
 99   end

      subroutine getf(t,g)
      implicit real*8 (a-h,o-z)
      dimension f1(100),df1(100),tmp(100)
      common/gridsyn/ xlam(4000),y1(4,4000,80,19),teff(80),grav(19),
     .                y11(4000),dy11(2,4000),nf,nt,ng,ic

c     Teff indexes
      klot=1
      khit=nt
 1    if(khit-klot.gt.1) then
         k=(khit+klot)/2
         if(teff(k).gt.t)then
            khit=k
         else
            klot=k
         endif
         go to 1
      endif

c     log g indexes
      klog=1
      khig=ng
 2    if(khig-klog.gt.1) then
         k=(khig+klog)/2
         if(grav(k).gt.g)then
            khig=k
         else
            klog=k
         endif
         go to 2
      endif

      do ij=1,nf
         do it=1,nt
            call splintg(klog,khig,ij,it,g,f1(it),df1(it))
         enddo
         call spline(teff,f1,nt,tmp)
         call splint(klot,khit,teff,f1,tmp,nt,t,y11(ij))
         call spline(teff,df1,nt,tmp)
         call splint(klot,khit,teff,df1,tmp,nt,t,dy11(2,ij))

         do ig=1,ng
            call splintt(klot,khit,ij,ig,t,f1(ig),df1(ig))
         enddo
         call spline(grav,df1,ng,tmp)
         call splint(klog,khig,grav,df1,tmp,ng,g,dy11(1,ij))
      enddo
      return
      end


      subroutine splineg
      implicit real*8 (a-h,o-z)
      dimension u(19)
      common/gridsyn/ xlam(4000),y1(4,4000,80,19),teff(80),grav(19),
     .                y11(4000),dy11(2,4000),nf,nt,ng,ic
      common/coeffs/ coeff(4,2,4000,80,19)
c
c     Evaluate spline coefficients for gravity interpolation
c
      do 1 ij=1,nf
      do 1 it=1,nt
         coeff(ic,2,ij,it,1)=0.
         u(1)=0.
         do 2 i=2,ng-1
           sig=(grav(i)-grav(i-1))/(grav(i+1)-grav(i-1))
           p=sig*coeff(ic,2,ij,it,i-1)+2.
           coeff(ic,2,ij,it,i)=(sig-1.)/p
           u(i)=(6.*((y1(ic,ij,it,i+1)-y1(ic,ij,it,i))/
     .          (grav(i+1)-grav(i))
     .          -(y1(ic,ij,it,i)-y1(ic,ij,it,i-1))/(grav(i)-grav(i-1)))/
     .          (grav(i+1)-grav(i-1))-sig*u(i-1))/p
    2    continue
         qn=0.
         un=0.
         coeff(ic,2,ij,it,ng)=(un-qn*u(ng-1))/
     .        (qn*coeff(ic,2,ij,it,ng-1)+1.)
         do 3 k=ng-1,1,-1
    3    coeff(ic,2,ij,it,k)=coeff(ic,2,ij,it,k)*coeff(ic,2,ij,it,k+1)
     .           +u(k)
    1 continue
      return
      end

      subroutine splinet
      implicit real*8 (a-h,o-z)
      dimension u(80)
      common/gridsyn/ xlam(4000),y1(4,4000,80,19),teff(80),grav(19),
     .                y11(4000),dy11(2,4000),nf,nt,ng,ic
      common/coeffs/ coeff(4,2,4000,80,19)
c
c     Evaluate spline coefficients for gravity interpolation
c
      do 1 ij=1,nf
      do 1 ig=1,ng
         coeff(ic,1,ij,1,ig)=0.
         u(1)=0.
         do 2 i=2,nt-1
           sig=(teff(i)-teff(i-1))/(teff(i+1)-teff(i-1))
           p=sig*coeff(ic,1,ij,i-1,ig)+2.
           coeff(ic,1,ij,i,ig)=(sig-1.)/p
           u(i)=(6.*((y1(ic,ij,i+1,ig)-y1(ic,ij,i,ig))/
     .          (teff(i+1)-teff(i))
     .          -(y1(ic,ij,i,ig)-y1(ic,ij,i-1,ig))/(teff(i)-teff(i-1)))/
     .          (teff(i+1)-teff(i-1))-sig*u(i-1))/p
    2    continue
         qn=0.
         un=0.
         coeff(ic,1,ij,nt,ig)=(un-qn*u(nt-1))/
     .        (qn*coeff(ic,1,ij,nt-1,ig)+1.)
         do 3 k=nt-1,1,-1
    3    coeff(ic,1,ij,k,ig)=coeff(ic,1,ij,k,ig)*coeff(ic,1,ij,k+1,ig)
     .           +u(k)
    1 continue
      return
      end

      subroutine splint(klo,khi,xa,ya,y2a,n,x,y)
      implicit real*8 (a-h,o-z)
      dimension xa(n),ya(n),y2a(n)
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end 

      subroutine splintg(klo,khi,ij,it,x,y,dy)
      implicit real*8 (a-h,o-z)
      common/gridsyn/ xlam(4000),y1(4,4000,80,19),teff(80),grav(19),
     .                y11(4000),dy11(2,4000),nf,nt,ng,ic
      common/coeffs/ coeff(4,2,4000,80,19)
      h=grav(khi)-grav(klo)
      if (h.eq.0.) pause 'bad grav input.'
      a=(grav(khi)-x)/h
      b=(x-grav(klo))/h
      y=a*y1(ic,ij,it,klo)+b*y1(ic,ij,it,khi)+
     !     ((a**3-a)*coeff(ic,2,ij,it,klo)+
     !     (b**3-b)*coeff(ic,2,ij,it,khi))*(h**2)/6.
      dy=(y1(ic,ij,it,khi)-y1(ic,ij,it,klo))/h-(3.*a**2-1.)/6.*h*
     !     coeff(ic,2,ij,it,klo)+(3.*b**2-1.)/6.*h*coeff(ic,2,ij,it,khi)
      return
      end

      subroutine spline(x,y,n,y2)
      implicit real*8 (a-h,o-z)
      dimension x(n),y(n),y2(n),u(100)
      y2(1)=0.
      u(1)=0.
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      qn=0.
      un=0.
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
12    y2(k)=y2(k)*y2(k+1)+u(k)
      return
      end

      subroutine splintt(klo,khi,ij,ig,x,y,dy)
      implicit real*8 (a-h,o-z)
      common/gridsyn/ xlam(4000),y1(4,4000,80,19),teff(80),grav(19),
     .                y11(4000),dy11(2,4000),nf,nt,ng,ic
      common/coeffs/ coeff(4,2,4000,80,19)
      h=teff(khi)-teff(klo)
      if (h.eq.0.) pause 'bad teff input.'
      a=(teff(khi)-x)/h
      b=(x-teff(klo))/h
      y=a*y1(ic,ij,klo,ig)+b*y1(ic,ij,khi,ig)+
     !     ((a**3-a)*coeff(ic,1,ij,klo,ig)+
     !     (b**3-b)*coeff(ic,1,ij,khi,ig))*(h**2)/6.
      dy=(y1(ic,ij,khi,ig)-y1(ic,ij,klo,ig))/h-(3.*a**2-1.)/6.*h*
     !     coeff(ic,1,ij,klo,ig)+(3.*b**2-1.)/6.*h*coeff(ic,1,ij,khi,ig)
      return
      end
