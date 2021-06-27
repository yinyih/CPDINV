!***************************** CWT.F90 *****************************
!
!This Program reads grid data from MinCurV.F90,then do continuous wavelet transform(CWT) using a new, hybrid wavelet-the 'FAN wavelet',
!which is described in Kirby and Swain, 2004.
!wavelets at different azimuths and scales are convolved with the entire magnetic anomaly grid, 
!avoiding the segmentation of the signal into ﬁnite‐size windows. 
!The wavelet transform offers a better compromise between spatial and wavenumber resolution than moving windows 
!since the wavelet transform uses a combination of multiple wavelet scales to construct a power spectrum at each grid point
!
!
! Written by:
!  Yihong Yin, Chun-Feng LI,
!  with additional code by Kirby
!
! 
! References:
!  Kirby, J.F. (2005). Which wavelet best reproduces the Fourier power
!    spectrum?, Computers and Geosciences, 31(7): 846-864.
!  Gaudreau, É., Audet, P., & Schneider, D. A. (2019). Mapping Curie Depth Across Western Canada From a Wavelet Analysis of Magnetic Anomaly Data.
!    Journal of Geophysical Research: Solid Earth, 124(5), 4365-4385. doi:10.1029/2018jb016726
! 
! 
! usage: CWT input file name, center wavenumber, scale increment, moving distance, output file name
!
! For further information about program, use command "man CWT " in Terminal.
!

Program cwt2d
  
  implicit none 
 
  character*80 infile,outfile,steplenstr,wvn0str,voicestr
  integer::steplen
  real(kind=8)::wvn0,voice
   integer::i,j,m,n,l

 !Store the parameters of input file 
  character*4 flag
  integer::nunit_in,columns,rows
  real::xmin,xmax,ymin,ymax,fmin,fmax
  real(kind=8),allocatable::data(:,:)
   
  !Store the parameters of CWT
   real(kind=8),allocatable,dimension(:,:):: rwt,iwt
   character*3 mother
   real(kind=8)::m1,phic,range,lam0,ds,dx,dy
   integer::nx,ny,ns
   real(kind=8)::lammin,lammax
   real(kind=8),allocatable,dimension(:,:,:):: w3d
   real(kind=8)::s0,freq,kpeak,power,df,scale,sumpower
   
  
   integer:: nxyA0

 
  real(kind=8),allocatable,dimension(:,:)::coordinate
  integer::c,windowsize
  real(kind=8)::wsize,o
 
  integer::nrows,ncolumns,ngrid
  character(10)::ch1='',ch2=''
 
 integer::narg
 narg=IARGC()
   if (narg==5) then 
    call getarg(1,infile)
    call getarg(2,wvn0str)
    call getarg(3,voicestr)
    call getarg(4,steplenstr)
    call getarg(5,outfile) 
  else
    write(*,*) 'pls input correct parameter!'
  endif
  
  read(wvn0str,*) wvn0
  read(voicestr,*) voice
  read(steplenstr,*) steplen
  
  !infile='mag_test.grd'
  !outfile='spectra_2.668.dat'
  !ds=0.25
  !wvn0=2.668
  !steplen=15

  windowsize=1
  ds=1/voice
  
  call search_unit(10,nunit_in)
  call open_old_file(nunit_in,infile,'formatted')
 
  !read input files' headers :
  !DSAA
  !numbers of west direction, numbers of south direction
  !xmin,xmax
  !ymin,ymax
  !zmin,zmax
  !all  in km 
  read(nunit_in,*) flag
  read(nunit_in,*) columns,rows
  read(nunit_in,*) xmin,xmax
  read(nunit_in,*) ymin,ymax
  read(nunit_in,*) fmin,fmax
  allocate(data(rows,columns))
  
  !read space domain input data :
  do i=1,rows
     read(nunit_in,*) (data(i,j),j=1,columns)
  end do
  
  nx=columns;ny=rows
  call primefac(nx)
  call primefac(ny)
  
  allocate(rwt(ny,nx),iwt(ny,nx))
 

  dx=(xmax-xmin)/(columns-1)
  dy=(ymax-ymin)/(rows-1)

  close(nunit_in)
  write(*,*)'CWT is running, pls waiting ...'

  

 lammin=2.d0*min(dx,dy)
 lammax=max((columns-1)*dx,(rows-1)*dy)
 ns=1+log(lammax/lammin)/(ds*log(2.d0))
 nrows=0
 ngrid=0
 allocate(w3d(ns,rows,columns))

 call cwt(w3d,s0,ns,columns,rows,nx,ny,data,rwt,iwt,dx,dy,wvn0,voice)

  open(15,file=outfile,form='formatted')
  write(15,*)ch1,ch2
  close(15)
  
    
     open(15,file=outfile,form='formatted',position='append')
     
     do i=1,rows,steplen
      if(i+windowsize>rows) exit
      nrows=nrows+1
       do j=1,columns,steplen  
         if(j+windowsize>columns) exit
         ngrid=ngrid+1
         write(15,50)
         50   format(6x,' f',4x,'log mean |F|') 
         do l=1,ns 
            scale=s0*2.d0**(-(l-1)*ds)
            freq=wvn0/(scale*2*3.1415926)
            power=sqrt(w3d(l,i,j))
             write(15,13) freq,log(power)
             13    format(1x,2e13.5) 
          end do
       end do 
     end do
     close(15)
    
 ncolumns=ngrid/nrows
 wsize=windowsize*(dx+dy)/2
 o=float(windowsize)/steplen
 allocate(coordinate(ngrid,2))
 c=1
 do i=1,nrows
   do j=1,ncolumns
     if((j.lt.ncolumns).and.(i.lt.rows)) then 
       coordinate(c,1)=xmin+wsize/2+(j-1)*wsize/o
       coordinate(c,2)=ymin+wsize/2+(i-1)*wsize/o
     end if
     if((j.eq.ncolumns).and.(i.lt.rows)) then 
       coordinate(c,1)=xmin+(columns-1)*dx-wsize/2
       coordinate(c,2)=ymin+wsize/2+(i-1)*wsize/o
     end if
     if((j.lt.ncolumns).and.(i.eq.rows)) then 
       coordinate(c,1)=xmin+wsize/2+(j-1)*wsize/o
       coordinate(c,2)=ymin+(rows-1)*dy-wsize/2
     end if
     if((j.eq.ncolumns).and.(i.eq.rows)) then 
       coordinate(c,1)=xmin+(columns-1)*dx-wsize/2
       coordinate(c,2)=ymin+(rows-1)*dy-wsize/2
     end if
     c=c+1
   end do 
 end do

 open(15,file=outfile,access='stream')
  write(ch1,'(i10)') ncolumns
  write(ch2,'(i10)') nrows
  write(15) ch2,ch1
  close(15)

  open(16,file='coordinate.dat',status='unknown')
    do i=1,nrows*ncolumns
      write(16,*)(coordinate(i,j),j=1,2)
    end do
  close(16)

  write(*,*)'CWT Successfully'

 deallocate(w3d,data,coordinate,iwt,rwt)
  
end Program cwt2d




subroutine cwt(w3d,s0,ns,nxA,nyA,nxA0,nyA0,fg,rfg,ifg,dxA,dyA,wvn0,voice)
 
  implicit none
 
  real(kind=8),parameter:: flag = -9999.d0
 
 real(kind=8)::fg(nyA,nxA)
  real(kind=8)::rfg(nyA0,nxA0),ifg(nyA0,nxA0)
  real(kind=8)::w3d(ns,nyA,nxA)
 ! arrays :
 integer,allocatable,dimension(:):: ncoi,kx,mmmm,lll
 integer,allocatable,dimension(:,:,:):: coi3dA,coi3dB,coi3dP

 real(kind=8),allocatable,dimension(:):: scale,theta1,coi2

 real(kind=8),allocatable,dimension(:,:):: rwg,iwg,xwvn,ywvn,wvnsq,kpsi,fgc,fhc,sumkpsi,wps
 real(kind=8),allocatable,dimension(:,:):: agg_az,ahh_az,agg,ahh

 complex(kind=8),allocatable,dimension(:,:):: agh_az,agh,coy,adm, &
 &erradm,errscy,aerradm,aerrscy,xpsi,xgwt,xhwt

 complex(kind=8),allocatable,dimension(:,:,:):: wgk,whk,q3d,c3d,eq3d,ec3d

 complex(kind=8),allocatable,dimension(:,:,:):: wg4d


 integer:: nxA,nyA,nxA0,nyA0,nxyA0,i,j,nxhalf2,nyhalf2,nfan,ntheta,m,ns,k,n, &
  &npts0,npts
 integer:: ncx,ncy,nxB,nyB,ncxy,i1,i2,j1,j2,ix1,ix2,jx1,jx2,ii,jj,iiy,jjx, &
 &iii,jjj
 integer:: ip,iypt,jxpt,nbx,nby,kk
 integer:: iucom,iug,iuh,iupts,ouwav,ouwg,ouwh,ouQ,ouQe,ouC,ouCe,oucoi
 integer(kind=8):: nsize

 real(kind=8):: twopi,deg2rad,extent,dfand,voice,ds,west,south,dxA,dyA, &
 &tapr,dkx,dky,dfan,wvn0,dtheta0,dtheta,lammin,lammax,s0,lambda,efw, &
 &u0,v0,theta,thetac,theta2,kmax,kpeak,kdelta,efwfac,qdtheta,thetalast
 real(kind=8):: xd,yd,gauss,rxpsi,ixpsi
 real(kind=8):: coifac,eastA,northA,coi,coiL,coiR,coiB,coiT,x,y,coitmp,scl
 real(kind=8):: dxB,dyB,xpt,ypt,aggp,ahhp,coyR,coyI,errscypR,errscypI,admR, &
 &admI,erradmpR,erradmpI
 real(kind=8):: gbytes,gigs
 complex(kind=8):: aghp,cohp,admp,errscyp,erradmp,xpsiij

 character(100):: filemag,filewave,filewg
 character(1):: havewc,freeboug,eqtopo,msup,aniso,zeros,detr,mir,tap,grid
 character(1):: coyQ,grdpts,coiopt,coitype,conv,cQerr,lgpmode
  
 twopi = 8 * atan(1.d0)
 deg2rad = twopi / 360
 
  
  ds=1/voice
  dfand=0!
  extent=180!
  if (extent >= 179.5d0) then
     aniso = 'i' 
  else
     aniso = 'a'
  end if

  nxyA0 = nxA0 * nyA0
  
  rfg=0.d0; ifg=0.d0;

  rfg(1:nyA,1:nxA)=fg

  call fft(rfg,ifg,nxyA0,nyA0,nyA0,1)
  call fft(rfg,ifg,nxyA0,nxA0,nxyA0,1)

  qdtheta = exp(-1.d0)
  dtheta0 = 2 * acos( (log(qdtheta + exp(-wvn0**2)) / wvn0**2) + 1 )
  if (aniso == 'i') then
    nfan = 1
    dfand = 0.d0               ! never used in isotropic inversion
  else
    nfan = 180 / dfand
  end if
  extent = extent * deg2rad
  dfan = dfand * deg2rad
 

  dtheta= 2.d0*sqrt(-2.d0*log(0.75))/wvn0
  ntheta = int(extent/dtheta)

  thetalast = (ntheta-1) * dtheta

  if(extent-thetalast<dtheta/2) ntheta=ntheta-1
  if (aniso == 'a') then
    if (mod(ntheta,2) == 0) ntheta = ntheta + 1
  end if
  if (ntheta < 2) ntheta = 2

  allocate(theta1(nfan)) ; theta1=0.d0

  if (aniso == 'i') then
    theta1(1) = 0.d0
  else
    do m=1,nfan
      thetac = (m-1)*dfan
      theta1(m) = thetac - (ntheta-1)*dtheta/2
    end do
   end if
  lammin = min(dxA , dyA)
  lammax = min( (nxA-1)*dxA , (nyA-1)*dyA )
  kmax = twopi / lammax

  kpeak = wvn0
  kdelta = (wvn0 + sqrt(wvn0**2 + 4)) / 2
  efwfac = kpeak

  s0 = efwfac / kmax
  
  allocate(scale(ns))


  do k=1,ns
    scale(k) = s0 * 2.d0**(-(k-1)*ds)
    efw = efwfac / scale(k)
    lambda = twopi / efw
 
  end do

 
  allocate(xwvn(nyA0,nxA0),ywvn(nyA0,nxA0),wvnsq(nyA0,nxA0))
  xwvn=0.d0; ywvn=0.d0; wvnsq=0.d0
  
  nxhalf2 = nxA0/2 + 2
  nyhalf2 = nyA0/2 + 2
  dkx = twopi / (nxA0*dxA)
  dky = twopi / (nyA0*dyA)
  
  do i=1,nyA0
   do j=1,nxA0
    xwvn(i,j) = (j-(j/nxhalf2)*nxA0-1)*dkx
    ywvn(i,j) = (i-(i/nyhalf2)*nyA0-1)*dky
   end do
  end do
  wvnsq = xwvn**2 + ywvn**2

   allocate(kpsi(nyA0,nxA0),sumkpsi(nyA0,nxA0))
   allocate(rwg(nyA0,nxA0),iwg(nyA0,nxA0))

  allocate(wps(nyA0,nxA0))

  fansA: do m=1,nfan
 
   
    w3d=0.d0;
  
   !=======================================================================
   scalesA: do k=1,ns

      wps=0.d0

     write(*,'("wavenumber domain multiplication at scale ",i4)') k
 
    !-----------------------------------------------------------------------
     azimk: do n=1,ntheta
 
      kpsi=0.d0; rwg=0.d0; iwg=0.d0 ;sumkpsi=0.d0 
 
      theta = theta1(m) + (n-1)*dtheta
 
      u0 = wvn0 * cos(theta)
      v0 = wvn0 * sin(theta)
 
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! complete Morlet wavelet :
      kpsi = exp( -((scale(k)*xwvn-u0)**2 + (scale(k)*ywvn-v0)**2) / 2)
      kpsi = kpsi - exp( -( scale(k)**2*wvnsq + wvn0**2 ) / 2)
      where (kpsi < 0.d0) kpsi = 0.d0

      ! energy (broadband) normalisation :
      kpsi = kpsi * sqrt(nxyA0 / sum(kpsi**2))
      
      rwg = rfg * kpsi
      iwg = ifg * kpsi
    
      !inverse Fourier transform
      call fft(rwg,iwg,nxyA0,nyA0,nyA0,-1)
      call fft(rwg,iwg,nxyA0,nxA0,nxyA0,-1)
      rwg=rwg/nxyA0 ; iwg=iwg/nxyA0
     
      wps =wps+(rwg*rwg+iwg*iwg)
      
     end do azimk
     
     w3d(k,:,:) = wps(1:nyA,1:nxA)/ntheta
    
   end do scalesA
 
  end do fansA
  
  deallocate(xwvn,ywvn,wvnsq)
  deallocate(theta1,scale)
  deallocate(kpsi,rwg,iwg,sumkpsi,wps)

end subroutine cwt
!**************************************************
!the function of FFT
!*************************************************
subroutine fft(a,b,ntot,n,nspan,isn)
      parameter (nmaxf=67,nmaxp=301)
      implicit double precision (a-h,o-z)
      double precision c72,s72
      double precision a(2*ntot),b(2*ntot-1)
      dimension nfac(20),np(nmaxp)
      dimension at(nmaxf),ck(nmaxf),bt(nmaxf),sk(nmaxf)
      equivalence (i,ii)

 ! array storage in nfac for a maximum of 20 prime factors of n.
 ! if n has more than one square-free factor, the product of the
 ! square-free factors must be less than or equal to 301.

 ! array storage for maximum prime factor of 67.

      maxf=nmaxf
      maxp=nmaxp
      ierr=0
      if (n.lt.2) return
      inc=isn
      c72=0.30901699437495d0
      s72=0.95105651629515d0
      s120=0.86602540378444d0
      rad=6.2831853071796d0
      if (isn.ge.0) go to 10
      s72=-s72
      s120=-s120
      rad=-rad
      inc=-inc
 10    nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*float(jc)*0.5d0
      i=0
      jf=0
 !determine the factors of n
      m=0
      k=n
      go to 20
 15    m=m+1
      nfac(m)=4
      k=k/16
 20    if (k-(k/16)*16.eq.0) go to 15
      j=3
      jj=9
      go to 30
 25    m=m+1
      nfac(m)=j
      k=k/jj
 30    if (mod(k,jj).eq.0) go to 25
      j=j+2
      jj=j**2
      if (jj.le.k) go to 30
      if (k.gt.4) go to 40
      kt=m
      nfac(m+1)=k
      if (k.ne.1) m=m+1
      go to 80
 40    if (k-(k/4)*4.ne.0) go to 50
      m=m+1
      nfac(m)=2
      k=k/4
 50    kt=m
      j=2
 60    if (mod(k,j).ne.0) go to 70
      m=m+1
      nfac(m)=j
      k=k/j
 70    j=((j+1)/2)*2+1
      if (j.le.k) go to 60
 80    if (kt.eq.0) go to 100
      j=kt
 90    m=m+1
      nfac(m)=nfac(j)
      j=j-1
      if (j.ne.0) go to 90
 !compute fourier transform
 100   sd=radf/float(kspan)
      cd=dsin(sd)
      cd=2.0d0*cd*cd
      sd=dsin(sd+sd)
      kk=1
      i=i+1
      if (nfac(i).ne.2) go to 400
 !transform for factor of 2 ( including rotation factor)
      kspan=kspan/2
      k1=kspan+2
 210   k2=kk+kspan
      ak=a(k2)
      bk=b(k2)
      a(k2)=a(kk)-ak
      b(k2)=b(kk)-bk
      a(kk)=a(kk)+ak
      b(kk)=b(kk)+bk
      kk=k2+kspan
      if (kk.le.nn) go to 210
      kk=kk-nn
      if (kk.le.jc) go to 210
      if (kk.gt.kspan) go to 800
 220   c1=1.0d0-cd
      s1=sd
 230   k2=kk+kspan
      ak=a(kk)-a(k2)
      bk=b(kk)-b(k2)
      a(kk)=a(kk)+a(k2)
      b(kk)=b(kk)+b(k2)
      a(k2)=c1*ak-s1*bk
      b(k2)=s1*ak+c1*bk
      kk=k2+kspan
      if (kk.lt.nt) go to 230
      k2=kk-nt
      c1=-c1
      kk=k1-k2
      if (kk.gt.k2) go to 230
      ak=cd*c1+sd*s1
      s1=(sd*c1-cd*s1)+s1
      c1=c1-ak
      kk=kk+jc
      if (kk.lt.k2) go to 230
      k1=k1+inc+inc
      kk=(k1-kspan)/2+jc
      if (kk.le.jc+jc) go to 220
      go to 100
 !transform for factors of 3 (optional code)
 320   k1=kk+kspan
      k2=k1+kspan
      ak=a(kk)
      bk=b(kk)
      aj=a(k1)+a(k2)
      bj=b(k1)+b(k2)
      a(kk)=ak+aj
      b(kk)=bk+bj
      ak=-0.5d0*aj+ak
      bk=-0.5d0*bj+bk
      aj=(a(k1)-a(k2))*s120
      bj=(b(k1)-b(k2))*s120
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      kk=k2+kspan
      if (kk.lt.nn) go to 320
      kk=kk-nn
      if (kk.le.kspan) go to 320
      go to 700
 !transform for factor of 4
 400   if (nfac(i).ne.4) go to 600
      kspnn=kspan
      kspan=kspan/4
 410   c1=1.0d0
      s1=0.0d0
 420   k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      akp=a(kk)+a(k2)
      akm=a(kk)-a(k2)
      ajp=a(k1)+a(k3)
      ajm=a(k1)-a(k3)
      a(kk)=akp+ajp
      ajp=akp-ajp
      bkp=b(kk)+b(k2)
      bkm=b(kk)-b(k2)
      bjp=b(k1)+b(k3)
      bjm=b(k1)-b(k3)
      b(kk)=bkp+bjp
      bjp=bkp-bjp
      if (isn.lt.0) go to 450
      akp=akm-bjm
      akm=akm+bjm
      bkp=bkm+ajm
      bkm=bkm-ajm
      if (s1.eq.0.0d0) go to 460
 430   a(k1)=akp*c1-bkp*s1
      b(k1)=akp*s1+bkp*c1
      a(k2)=ajp*c2-bjp*s2
      b(k2)=ajp*s2+bjp*c2
      a(k3)=akm*c3-bkm*s3
      b(k3)=akm*s3+bkm*c3
      kk=k3+kspan
      if (kk.le.nt) go to 420
 440   c2=cd*c1+sd*s1
      s1=(sd*c1-cd*s1)+s1
      c1=c1-c2
      c2=c1*c1-s1*s1
      s2=2.0d0*c1*s1
      c3=c2*c1-s2*s1
      s3=c2*s1+s2*c1
      kk=kk-nt+jc
      if (kk.le.kspan) go to 420
      kk=kk-kspan+inc
      if (kk.le.jc) go to 410
      if (kspan.eq.jc) go to 800
      go to 100
 450   akp=akm+bjm
      akm=akm-bjm
      bkp=bkm-ajm
      bkm=bkm+ajm
      if (s1.ne.0) go to 430
 460   a(k1)=akp
      b(k1)=bkp
      a(k2)=ajp
      b(k2)=bjp
      a(k3)=akm
      b(k3)=bkm
      kk=k3+kspan
      if (kk.le.nt) go to 420
      go to 440
 !transform for factor of 5 (optional code)
 510   c2=c72*c72-s72*s72
      s2=2.0d0*c72*s72
 520   k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      k4=k3+kspan
      akp=a(k1)+a(k4)
      akm=a(k1)-a(k4)
      bkp=b(k1)+b(k4)
      bkm=b(k1)-b(k4)
      ajp=a(k2)+a(k3)
      ajm=a(k2)-a(k3)
      bjp=b(k2)+b(k3)
      bjm=b(k2)-b(k3)
      aa=a(kk)
      bb=b(kk)
      a(kk)=aa+akp+ajp
      b(kk)=bb+bkp+bjp
      ak=akp*c72+ajp*c2+aa
      bk=bkp*c72+bjp*c2+bb
      aj=akm*s72+ajm*s2
      bj=bkm*s72+bjm*s2
      a(k1)=ak-bj
      a(k4)=ak+bj
      b(k1)=bk+aj
      b(k4)=bk-aj
      ak=akp*c2+ajp*c72+aa
      bk=bkp*c2+bjp*c72+bb
      aj=akm*s2-ajm*s72
      bj=bkm*s2-bjm*s72
      a(k2)=ak-bj
      a(k3)=ak+bj
      b(k2)=bk+aj
      b(k3)=bk-aj
      kk=k4+kspan
      if (kk.lt.nn) go to 520
      kk=kk-nn
      if (kk.le.kspan) go to 520
      go to 700
 !transform  for  odd  factors
 600   k=nfac(i)
      kspnn=kspan
      kspan=kspan/k
      if (k.eq.3) go to 320
      if (k.eq.5) go to 510
      if (k.eq.jf) go to 640
      jf=k
      s1=rad/float(k)
      c1=dcos(s1)
      s1=dsin(s1)
      if (jf.gt.maxf) go to 996
      ck(jf)=1.0d0
      sk(jf)=0.0d0
      j=1
 630   ck(j)=ck(k)*c1+sk(k)*s1
      sk(j)=ck(k)*s1-sk(k)*c1
      k=k-1
      ck(k)=ck(j)
      sk(k)=-sk(j)
      j=j+1
      if (j.lt.k) go to 630
 640   k1=kk
      k2=kk+kspnn
      aa=a(kk)
      bb=b(kk)
      ak=aa
      bk=bb
      j=1
      k1=k1+kspan
 650   k2=k2-kspan
      j=j+1
      at(j)=a(k1)+a(k2)
      ak=at(j)+ak
      bt(j)=b(k1)+b(k2)
      bk=bt(j)+bk
      j=j+1
      at(j)=a(k1)-a(k2)
      bt(j)=b(k1)-b(k2)
      k1=k1+kspan
      if (k1.lt.k2) go to 650
      a(kk)=ak
      b(kk)=bk
      k1=kk
      k2=kk+kspnn
      j=1
 660   k1=k1+kspan
      k2=k2-kspan
      jj=j
      ak=aa
      bk=bb
      aj=0.0d0
      bj=0.0d0
      k=1
 670   k=k+1
      ak=at(k)*ck(jj)+ak
      bk=bt(k)*ck(jj)+bk
      k=k+1
      aj=at(k)*sk(jj)+aj
      bj=bt(k)*sk(jj)+bj
      jj=jj+j
      if(jj.gt.jf) jj=jj-jf
      if(k.lt.jf) go to 670
      k=jf-j
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      j=j+1
      if (j.lt.k) go to 660
      kk=kk+kspnn
      if (kk.le.nn) go to 640
      kk=kk-nn
      if (kk.le.kspan) go to 640
 ! multiply by rotation factor (except for factors of 2 and 4)
 700   if (i.eq.m) go to 800
      kk=jc+1
 710   c2=1.0-cd
      s1=sd
 720   c1=c2
      s2=s1
      kk=kk+kspan
 730   ak=a(kk)
      a(kk)=c2*ak-s2*b(kk)
      b(kk)=s2*ak+c2*b(kk)
      kk=kk+kspnn
      if (kk.le.nt) go to 730


      ak=s1*s2
      s2=s1*c2+c1*s2
      c2=c1*c2-ak
      kk=kk-nt+kspan
      if (kk.le.kspnn) go to 730
      c2=c1-(cd*c1+sd*s1)
      s1=s1+(sd*c1-cd*s1)
      kk=kk-kspnn+jc
      if (kk.le.kspan) go to 720
      kk=kk-kspan+jc+inc
      if (kk.le.jc+jc) go to 710
      go to 100
 !permute results to normal order--- done in two stages
 !permutation for square factors of n
 800   np(1)=ks
      if (kt.eq.0) go to 890
      k=kt+kt+1
      if (m.lt.k) k=k-1
      j=1
      np(k+1)= jc
 810   np(j+1)=np(j)/nfac(j)
      np(k)=np(k+1)*nfac(j)
      j=j+1
      k=k-1
      if (j.lt.k) go to 810
      k3=np(k+1)
      kspan=np(2)
      kk=jc+1
      k2=kspan+1
      j=1
      if(n.ne.ntot) goto 850
 !permutation for single-variate transform (optional code)
 820   ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=kspan+k2
      if (k2.lt.ks) go to 820
 830   k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if (k2.gt.np(j)) go to 830
      j=1
 840   if (kk.lt.k2) go to 820
      kk=kk+inc
      k2=kspan+k2
      if (k2.lt.ks) go to 840
      if (kk.lt.ks) go to 830
      jc=k3
      go to 890
 ! permutation for multivariate transform
 850   k=kk+jc
 860   ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=k2+inc
      if (kk.lt.k) go to 860
      kk=kk+ks-jc
      k2=k2+ks-jc
      if (kk.lt.nt) go to 850
      k2=k2-nt+kspan
      kk=kk-nt+jc
      if (k2.lt.ks) go to 850
 870   k2=k2-np(j)
      j=j+1
      k2=np(j+1) +k2
      if (k2.gt.np(j)) go to 870
      j=1
 880   if (kk.lt.k2) go to 850
      kk=kk+jc
      k2=kspan+k2
      if (k2.lt.ks) go to 880
      if (kk.lt.ks) go to 870
      jc=k3
 890   if (2*kt+1.ge.m) return
      kspnn=np(kt+1)
 !permutation for square free factors of n
      j=m-kt
      nfac(j+1)=1
 900   nfac(j)=nfac(j)*nfac(j+1)
      j=j-1
      if (j.ne.kt) go to 900
      kt=kt+1
      nn=nfac(kt)-1
      if (nn.gt.maxp) go to 998
      jj=0
      j=0
      go to 906
 902   jj=jj-k2
      k2=kk
      k=k+1
      kk=nfac(k)
 904   jj=kk+jj
      if (jj.ge.k2) go to 902
      np(j)=jj
 906   k2=nfac(kt)
      k=kt+1
      kk=nfac(k)
      j=j+1
      if (j.le.nn) go to 904
 !determine the permutation cycles of length greter than 1
      j=0
      go to 914
 910   k=kk
      kk=np(k)
      np(k)=-kk
      if (kk.ne.j) go to 910
      k3=kk
 914   j=j+1
      kk=np(j)
      if (kk.lt.0) go to 914
      if (kk.ne.j) go to 910
      np(j)=-j
      if (j.ne.nn) go to 914
      maxf=inc*maxf
 !reorder a and b, following the permutation cycles
      go to 950
  924   j=j-1
      if (np(j).lt.0) go to 924
      jj=jc
 926   kspan=jj
      if (jj.gt.maxf) kspan=maxf
      jj=jj-kspan
      k=np(j)
      kk=jc*k+ii+jj
      k1=kk+kspan
      k2=0
 928   k2=k2+1
      at(k2)=a(k1)
      bt(k2)=b(k1)
      k1=k1-inc
      if (k1.ne.kk) go to 928
 932   k1=kk+kspan
      k2=k1-jc*(k+np(k))
      k=-np(k)
 936   a(k1)=a(k2)
      b(k1)=b(k2)
      k1=k1-inc
      k2=k2-inc
      if (k1.ne.kk) go to 936
      kk=k2
      if (k.ne.j) go to 932
      k1=kk+kspan
      k2=0
 940   k2=k2+1
      a(k1)=at(k2)
      b(k1)=bt(k2)
      k1=k1-inc
      if (k1.ne.kk) go to 940
      if (jj.ne.0) go to 926
      if (j.ne.1) go to 924
 950   j=k3+1
      nt=nt-kspnn
      ii=nt-inc+1
      if (nt.ge.0) go to 924
      return

 ! error finish, insufficient array storage
  996 continue
      write (6,997)
      ierr=1
  997 format(' the array dimension parameter maxf is insufficient')

  998 continue
      write(6,999)
      ierr=1
  999 format(' the array dimension parameter maxp is insufficient')

      return
      end

subroutine primefac(n)
	implicit none
	integer,parameter:: nfac=19,maxf=67,maxp=301
	integer:: pr(nfac),fact(maxp),kpr,new,kfct,k,prod
	integer,intent(inout):: n

	pr=(/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67/)
	fact=0

	kpr=1
	kfct=1
	new=n

 ! loop to determine all the square-free prime factors of n :

	all: do

 !loop to determine all the prime factors of n :

	 factors: do

	  do k=1,nfac
	   if (new==pr(k)) exit factors
	  enddo

 !if n can be divided by a prime no. without remainder :

	  if (mod(new,pr(kpr))==0) then
	    new=new/pr(kpr)
	    fact(kfct)=pr(kpr)
	    kfct=kfct+1
	    do k=1,nfac
	     if (new==pr(k)) fact(kfct)=new
	    enddo
	    kpr=1                              ! start from pr(1) again

 ! if n cannot be divided by a prime no. :

	  else
	    kpr=kpr+1                          ! try the next prime no.

 ! if no more primes, n -> n+1 :

	    if (kpr>nfac) then
	      n=n+1
	      kpr=1
	      new=n
	      fact=0
	      kfct=1
	      cycle factors
	    endif

	  endif

	 enddo factors

 !now determine the square-free factors of n :

	 prod=1
	 k=1
	 squarefree: do
	  if (fact(k)/=fact(k+1)) then
	    prod=prod*fact(k)
	    k=k+1
	  else
	    k=k+2
	  endif
	  if (k>kfct) exit squarefree
	 enddo squarefree

 ! n -> n+1 if product of square-free factors > maxp

	 if (prod>maxp) then
	   n=n+1
	   kpr=1
	   new=n
	   fact=0
	   kfct=1
	   cycle
	 else
	   exit all
	 endif

	enddo all

	return
end subroutine primefac

!**********************************************
!Open an  existing file. 
!*************************************************
subroutine open_old_file(nunit,filename,form)
 character*(*) filename,form
 integer nunit
 logical file_exist

 inquire(file=filename,exist=file_exist)
 if(file_exist) then
  open(nunit,file=filename,status='old',form=form)
 else
  write(*,*) 'file:',filename,'doesnot exist'
  write(*,*) 'please input filename='
  read(*,*) filename
  open(nunit,file=filename,status='old',form=form)
 end if
end subroutine open_old_file
!**********************************************
!From file channel number Nunit_ Start starts 
!to search for an unopened file channel number Nunit_ in
!*********************************************************
subroutine search_unit(nunit_start,nunit_in)
    logical unit_open
    integer nunit_in,nunit_start

    nunit_in=nunit_start
    unit_open=.TRUE.
    do while(unit_open)
      nunit_in=nunit_in+1
      inquire(unit=nunit_in,opened=unit_open)
    end do
end subroutine search_unit
