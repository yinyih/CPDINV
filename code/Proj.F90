
!***************************** Proj.F90 *****************************
!
!This program reads (longitude, latitude) positions from standard input and computes (x,y) coordinates
!using the ellipsoidal transverse mercator map projection. Optionally, it can read (x,y) positions
!and compute (longitude, latitude) values in the inverse transformation.
!
!
! Credits:  Yihong Yin and Chun-Feng Li @ Zhejiang University
!
! 
! References:
!
!  Yin Y.H, Li C.-F., Lu Y., 2021, Estimating Curie-point depths using both wavelet-based and Fourier spectral centroid methods in the western Pacific marginal seas. 
!    Geophysical Journal International,  DOI: 10.1093/gji/ggab257
! 
! usage: Proj input file name, -1/1, central longitude, base latitude, output file name
!
! For further information about program, type "man Project " in Terminal.
!

Program Proj

 use,intrinsic ::ISO_C_BINDING
 implicit none
  !Command line read parameters: 
 !
 
 character*80 infile,file_type,outfile,cmstr,blstr,pro_type
 integer::status1,narg,type
 integer::access,status!test if the file exists
 real(kind=4)::cm,bl
  

  
  real,allocatable::gridData(:,:)


 integer::grd=2,xyz=1,nvals=0,i,j
 real,allocatable::val(:,:)
 real,allocatable::lon(:),lat(:),xkm(:),ykm(:)
 
  character*4 flag
  integer::ny,nx
  real::xmin,xmax,ymin,ymax,fmin,fmax,dx,dy
  real,allocatable::data(:,:)
  
 narg=IARGC()
 if (narg==5) then 
  call getarg(1,infile)
  call getarg(2,pro_type)
  call getarg(3,cmstr)
  call getarg(4,blstr)
  call getarg(5,outfile)
 else
  write(*,*) 'pls input correct parameter'
 endif

 read(pro_type,*) type
 read(cmstr,*)cm
 read(blstr,*)bl

 !infile='z_pf.dat'
 !outfile='pf.dat'
 !cm=116.156
 !bl=13.259
 !type=-1

 i=index(infile,'.')
 file_type=infile(i+1:)
 write(*,*)'Project is running, pls wait...'
 !****************Open file and Read data Start***************************

  !*************Read grid data******************************
 if(file_type.eq.'grd') then
   open(grd,file=infile,status ='OLD',access="stream",form="unformatted")
   read(grd,*) flag
  read(grd,*) ny,nx
  read(grd,*) xmin,xmax
  read(grd,*) ymin,ymax
  read(grd,*) fmin,fmax
  
   allocate (gridData(nx,ny))
   allocate (val(ny*nx,3))
   write(*,*) 'The head file information is as follows:'
   write(*,*) 'flag=',flag
   write(*,*) 'rows/columns=',ny,nx
   write(*,*) 'X_Min=',xmin,'X_Max=',xmax
   write(*,*) 'Y_Min=',ymin,'Y_Max=',ymax
   write(*,*) 'Z_Min=',fmin,'Z_Max=',fmax

   do i=1,ny
    read(grd,*) (gridData(i,j),j=1,nx)
 end do
 close(grd)

   dx=(xmax-xmin)/(nx-1)
   dy=(ymax-ymin)/(ny-1)
   do j=1,ny
     do i=1,nx
      val(i+(j-1)*nx,1)=xmin+(i-1)*dx
      val(i+(j-1)*nx,2)=ymax-(j-1)*dy
      val(i+(j-1)*nx,3)=gridData(i,ny-j+1)
      end do
    end do
   nvals=ny*nx
   deallocate (gridData)
   write(*,*) 'Read Grd Successfully'
   close(grd)
  
 !*************Read xyz data******************************
 else
   open(xyz,file = infile,status='OLD',position='rewind')
    do while(.true.)
     read(xyz,*,iostat=status1)
     if(status1 /= 0) exit
     nvals=nvals+1
    enddo
    allocate (val(nvals,3))
    rewind(xyz)
    do i=1,nvals,1
      read(xyz,*,iostat=status1) (val(i,j),j=1,3)
      if(status1/=0) exit
    end do
    write(*,*) 'Read XYZ Successfully'
   close(xyz)
 endif
   
 !***************************Project Start***********************
 ! Start in 2019,11,10 by Yin Yihong
 allocate (lon(nvals),lat(nvals),xkm(nvals),ykm(nvals))

 
 
 if(type.eq.1) then
   lon=val(:,1)
   lat=val(:,2)
   do i=1,nvals,1
     call etmftw(lon(i),lat(i),xkm(i),ykm(i),cm,bl)
   end do
   
   val(:,1)=xkm
   val(:,2)=ykm
    write(*,*) 'Forward Projection done successfully'
 end if
  if(type.eq.-1) then
   xkm=val(:,1)
   ykm=val(:,2)
   do i=1,nvals,1
     call etminv(xkm(i),ykm(i),lon(i),lat(i),cm,bl)
   end do
    
   val(:,1)=lon
   val(:,2)=lat
    write(*,*) 'Inverse Projection done successfully'
 end if

 deallocate (lon,lat,xkm,ykm)

 !!!! end in 2019,11,13 by Yin Yihong
 !***************************Project End***************************\

 open(4,file=outfile,status='unknown')
 do i=1,nvals
   write(4,*) (val(i,j), j=1,3)
 end do
 close(4)
 deallocate(val)
  
  
end Program Proj

!******************************************************
! Orthographic projection: Project (longitude,latitude) to the (x,y)
!****************************Transverse Mercator Projection******************
subroutine etmftw(xdeg,ydeg,xkm,ykm,cm,bl)
    
  real k0,n,mfcn,m,m0,lam,lam0
  parameter (a=6378.2064)
  parameter (e=0.0822719,esq=0.676866E-2)
  parameter (k0=0.9996)
  parameter (deg2rad=0.17453293E-1,rad2deg=57.2957795)
  common/etmconsts/ iwarn,epsq,m0,c0,c2,c4,c6,lam0,phi0,cmdeg  

  c0=1.-e**2/4.-3.*e**4/64.-5.*e**6/256.
  c2=3.*e**2/8.+3.*e**4/32.+45.*e**6/1024.
  c4=15.*e**4/256.+45.*e**6/1024.
  c6=35.*e**6/3072.

  lam0=cm*deg2rad
  phi0=bl*deg2rad
  m0=a*(c0*phi0-c2*sin(2.*phi0)+c4*sin(4.*phi0)-c6*sin(6.*phi0)) 
  epsq=esq/(1.0-esq)
      
  lam=xdeg*deg2rad
  phi=ydeg*deg2rad

  n=a/sqrt(1.-esq*sin(phi)**2)
  tnphi=tan(phi)
  t=tnphi**2
 
  csphi=cos(phi)
  c=epsq*csphi**2
  biga=csphi*(lam-lam0)
  m=a*(c0*phi-c2*sin(2.*phi)+c4*sin(4.*phi)-c6*sin(6.*phi))

  xkm=k0*n*(biga+(1.-t+c)*biga**3/6.+(5.-18.*t+t**2+72.*c-58.*epsq)*biga**5/120.)
  ykm=k0*(m-m0+n*tnphi*(.5*biga**2 +(5.-t+9.*c+4.*c**2)*biga**4/24. + (61.-58.*t+t**2+600.*c-330.*epsq)*biga**6/720.))

end subroutine etmftw

!******************************************************
!Back projectionProject (x,y) to  (longitude,latitude)
!****************************Transverse Mercator Projection******************
subroutine etminv(xkm,ykm,xdeg,ydeg,cm,bl)

    
  real k0,n,mu,m,m0,lam,lam0
  parameter (a=6378.2064)
  parameter (e=0.0822719,esq=0.676866E-2)
  parameter (k0=0.9996)
  parameter (deg2rad=0.17453293E-1,rad2deg=57.2957795)
  common/etmconsts/ iwarn,epsq,m0,c0,c2,c4,c6,lam0,phi0,cmdeg  

  c0=1.-e**2/4.-3.*e**4/64.-5.*e**6/256.
  c2=3.*e**2/8.+3.*e**4/32.+45.*e**6/1024.
  c4=15.*e**4/256.+45.*e**6/1024.
  c6=35.*e**6/3072.

  lam0=cm*deg2rad
  phi0=bl*deg2rad
  m0=a*(c0*phi0-c2*sin(2.*phi0)+c4*sin(4.*phi0)-c6*sin(6.*phi0))

  m=m0 + ykm/k0
  e1=(1.-sqrt(1.-esq))/(1.+sqrt(1.-esq))
  mu=m/(a*(1.-esq/4.-3.*esq**2/64.-5.*esq**3/256.))
  phi1=mu+(3.*e1/2.-27.*e1**3/32.)*sin(2.*mu) + (21.*e1**2/16.-55.*e1**4/32.)*sin(4.*mu) +(151.*e1**3/96.)*sin(6.*mu)
  c1=epsq*cos(phi1)**2
  t1=tan(phi1)**2
  n1=a/sqrt(1.-esq*sin(phi1)**2)
  r1=a*(1-esq)/((1-esq*sin(phi1)**2)*sqrt(1-esq*sin(phi1)**2))
  d=xkm/(n1*k0)
   
  phi=phi1-(n1*tan(phi1)/r1)*(d**2/2.-(5.+3.*t1+10.*c1-4.*c1**2-9.*epsq)*d**4/24.&
      &+ (61.+90*t1+298*c1+45*t1**2-252.*epsq-3.*c1**2)*d**6/720.)
  lam=lam0+(d-(1.+2.*t1+c1)*d**3/6.+(5.-2.*c1+28.*t1-3.*c1**2+8.*epsq+24.*t1**2)*d**5/120.) /cos(phi1)

  xdeg=lam*rad2deg
  ydeg=phi*rad2deg


end subroutine etminv
