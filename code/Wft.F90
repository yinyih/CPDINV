
  
!***************************** FFT.F90 *****************************
!
!reads grid data from MinCurV,then do fast fourier transform and calculate area-weighted radial spectrumin dw band in different windowsize.
!It will output a spectrun file and a new coordinate file based on different windowsize and steplength.
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
!  Phillips, J.D., 1997. Potential-field geophysical software for the PC, version 2.2. U.S.
!    Geological Survey, Open-File, Report, 97-725.	
! 
! 
! usage: FFT input file, window size, step length, output file
!
! For further information about program, use command "man FFT " in Terminal.
!


program Wft
 implicit none
 
 
  character*80 infile,outfile,size,len,dwh,w
  integer::windowsize,steplen
  real::dw,wmax
 
 !Store the parameters of input file 
  character*4 flag
  integer::nunit_in,columns,rows
  real::xmin,xmax,ymin,ymax,fmin,fmax,dx,dy
  real,allocatable::data(:,:)
  
  
  !Store the parameters of FFT
  real,allocatable,dimension(:,:):: rfg,ifg,afg
  integer:: nxyA0
  integer::i,j,m,n
  
  
    !Store the parameters of power spectrum
  real,dimension(1000)::sum
  integer,dimension(1000)::number
  integer::narg

  
  integer::nrows,ncolumns,ngrid
  

  character(10)::ch1='',ch2=''
  
  
  real,allocatable,dimension(:,:)::coordinate
  integer::c
  real::wsize,o
 


  narg=IARGC()
   if (narg==6) then 
    call getarg(1,infile)
    call getarg(2,size)
    call getarg(3,len)
    call getarg(4,dwh)
    call getarg(5,w)
    call getarg(6,outfile) 
  else
    write(*,*) 'pls input correct parameter!'
  endif

  read(size,*) windowsize
  read(len,*) steplen
  read(dwh,*) dw
  read(w,*) wmax

  !infile='mag.grd'
  !windowsize=64
 ! steplen=32
  !dw=0.006
  !wmax=0.15
  !outfile='spectra.dat'
  write(*,*)'FFT is running, pls waiting ...'

  open(15,file=outfile,form='formatted')
  write(15,*)ch1,ch2
  close(15)

 call search_unit(10,nunit_in)
 call open_old_file(nunit_in,infile,'formatted')

 read(nunit_in,*) flag
 read(nunit_in,*) columns,rows
 allocate(data(rows,columns))
 read(nunit_in,*) xmin,xmax
 read(nunit_in,*) ymin,ymax
 read(nunit_in,*) fmin,fmax
 do i=1,rows
    read(nunit_in,*) (data(i,j),j=1,columns)
 end do
 close(nunit_in)
  
 allocate(rfg(windowsize,windowsize),ifg(windowsize,windowsize),afg(windowsize,windowsize))
 
 dx=((xmax-xmin)/(columns-1))
 dy=((ymax-ymin)/(rows-1))
 
 nxyA0=windowsize*windowsize
 nrows=0
 ncolumns=0
 ngrid=0
 

  do m=1,rows,steplen
    if(m+windowsize>rows) exit
    nrows=nrows+1
   do n=1,columns,steplen
     rfg=0.d0
     ifg=0.d0
     if(n+windowsize>columns) exit
     ngrid=ngrid+1
    !Start FFT at each window
	do i=1,windowsize
	 do j=1,windowsize
         rfg(i,j)=data(m+i-1,n+j-1)
      end do
	end do
     
     call fft(rfg,ifg,nxyA0,windowsize,windowsize,1)
     call fft(rfg,ifg,nxyA0,windowsize,nxyA0,1)
     !Finish FFT at each window
     call f_sras(outfile,rfg,-ifg,sum,number,windowsize,windowsize,windowsize*2,dy,dx,dw,wmax)    
   
    end do
  end do
 
 ncolumns=ngrid/nrows
 wsize=windowsize*(dx+dy)/2
 o=float(windowsize/steplen)
 allocate(coordinate(ngrid,2))
 
 ! Output coordinate file
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
 
  ! Output coordinate file
  open(16,file='coordinate.dat',status='unknown')
    do i=1,nrows*ncolumns
      write(16,*)(coordinate(i,j),j=1,2)
    end do
  close(16)
  
  ! Output spectrum file
  open(15,file=outfile,access='stream')
  write(ch1,'(i10)') nrows
  write(ch2,'(i10)') ncolumns
 
  write(15) ch2,ch1
  
  close(15)
  
  
  write(*,*)'FFT Successfully'

 deallocate(rfg,ifg,afg,data)
end	program Wft	



!**************************************************
!the function of FFT
!*************************************************
subroutine fft(a,b,ntot,n,nspan,isn)

 implicit none
 
 integer,parameter:: nmaxf=67, nmaxp=301
 
 integer:: nfac(20), np(nmaxp)
 
 real(kind=4)::a(2*ntot), b(2*ntot-1)
 real(kind=4),dimension(nmaxf):: at,bt,ck,sk
 
 integer:: ntot,n,nspan
 integer:: isn
 integer:: maxf,maxp,ierr,inc,nt,ks,kspan,nn,jc,i,jf,m,k,j,jj,kt,kk,k1,k2, &
  &kspnn,k3,k4,ii
 
 real(kind=4):: c72,s72,s120,rad,radf,sd,cd,ak,bk,c1,s1,aj,bj,akp,akm, &
  &ajm,ajp,bkp,bkm,bjp,bjm,c2,s2,c3,s3,aa,bb
 
 equivalence (i,ii)
 
 ! array storage in nfac for a maximum of 20 prime factors of n.
 ! if n has more than one square-free factor, the product of the
 ! square-free factors must be less than or equal to 301.
 
 ! array storage for maximum prime factor of 67.
 
 maxf = nmaxf
 maxp = nmaxp
 ierr = 0
 if (n<2) return
  inc = isn
  c72 = 0.30901699437495d0
  s72 = 0.95105651629515d0
  s120 = 0.86602540378444d0
  rad = 6.2831853071796d0
 if (isn>=0) go to 10
  s72 = -s72
  s120 = -s120
  rad = -rad
  inc = -inc
  10 nt = inc * ntot
  ks = inc * nspan
  kspan = ks
  nn = nt - inc
  jc = ks / n
  radf = rad * jc / 2
  i = 0
  jf = 0
 !  determine the factors of n
  m = 0
  k = n
 go to 20
  15 m = m + 1
  nfac(m) = 4
  k = k / 16
 20 if (k-(k/16)*16==0) go to 15
 j = 3
 jj = 9
 go to 30
 25 m = m + 1
 nfac(m) = j
 k = k / jj
 30 if (mod(k,jj)==0) go to 25
 j = j + 2
 jj = j**2
 if (jj<=k) go to 30
 if (k>4) go to 40
 kt = m
 nfac(m+1) = k
 if (k/=1) m = m + 1
 go to 80
 40 if (k-(k/4)*4/=0) go to 50
 m = m + 1
 nfac(m) = 2
 k = k / 4
 50 kt = m
 j = 2
 60 if (mod(k,j)/=0) go to 70
 m = m + 1
 nfac(m) = j
 k = k / j
 70 j = ((j+1)/2)*2 + 1
 if (j<=k) go to 60
 80 if (kt==0) go to 100
 j = kt
 90 m = m + 1
 nfac(m) = nfac(j)
 j = j - 1
 if (j/=0) go to 90
 !  compute fourier transform
 100 sd = radf / kspan
 cd = sin(sd)
 cd = 2 * cd * cd
 sd = sin(sd+sd)
 kk = 1
 i = i + 1
 if (nfac(i)/=2) go to 400
 !  transform for factor of 2 ( including rotation factor)
 kspan = kspan / 2
 k1 = kspan + 2
 210 k2 = kk + kspan
 ak = a(k2)
 bk = b(k2)
 a(k2) = a(kk) - ak
 b(k2) = b(kk) - bk
 a(kk) = a(kk) + ak
 b(kk) = b(kk) + bk
 kk = k2 + kspan
 if (kk<=nn) go to 210
 kk = kk - nn
 if (kk<=jc) go to 210
 if (kk>kspan) go to 800
 220 c1 = 1 - cd
 s1 = sd
 230 k2 = kk + kspan
 ak = a(kk) - a(k2)
 bk = b(kk) - b(k2)
 a(kk) = a(kk) + a(k2)
 b(kk) = b(kk) + b(k2)
 a(k2) = c1*ak - s1*bk
 b(k2) = s1*ak + c1*bk
 kk = k2 + kspan
 if (kk<nt) go to 230
 k2 = kk - nt
 c1 = -c1
 kk = k1 - k2
 if (kk>k2) go to 230
 ak = cd*c1 + sd*s1
 s1 = (sd*c1-cd*s1) + s1
 c1 = c1 - ak
 kk = kk + jc
 if (kk<k2) go to 230
 k1 = k1 + inc + inc
 kk = (k1-kspan)/2 + jc
 if (kk<=jc+jc) go to 220
 go to 100
 !  transform for factors of 3 (optional code)
 320 k1 = kk + kspan
 k2 = k1 + kspan
 ak = a(kk)
 bk = b(kk)
 aj = a(k1) + a(k2)
 bj = b(k1) + b(k2)
 a(kk) = ak + aj
 b(kk) = bk + bj
 ak = -aj/2 + ak
 bk = -bj/2 + bk
 aj = (a(k1)-a(k2)) * s120
 bj = (b(k1)-b(k2)) * s120
 a(k1) = ak - bj
 b(k1) = bk + aj
 a(k2) = ak + bj
 b(k2) = bk - aj
 kk = k2 + kspan
 if (kk<nn) go to 320
 kk = kk - nn
 if (kk<=kspan) go to 320
 go to 700
 !  transform for factor of 4
 400 if (nfac(i)/=4) go to 600
 kspnn = kspan
 kspan = kspan/4
 410 c1 = 1.d0
 s1 = 0.d0
 420 k1 = kk + kspan
 k2 = k1 + kspan
 k3 = k2 + kspan
 akp = a(kk) + a(k2)
 akm = a(kk) - a(k2)
 ajp = a(k1) + a(k3)
 ajm = a(k1) - a(k3)
 a(kk) = akp + ajp
 ajp = akp - ajp
 bkp = b(kk) + b(k2)
 bkm = b(kk) - b(k2)
 bjp = b(k1) + b(k3)
 bjm = b(k1) - b(k3)
 b(kk) = bkp + bjp
 bjp = bkp - bjp
 if (isn<0) go to 450
 akp = akm - bjm
 akm = akm + bjm
 bkp = bkm + ajm
 bkm = bkm - ajm
 if (s1==0.d0) go to 460
 430 a(k1) = akp*c1 - bkp*s1
 b(k1) = akp*s1 + bkp*c1
 a(k2) = ajp*c2 - bjp*s2
 b(k2) = ajp*s2 + bjp*c2
 a(k3) = akm*c3 - bkm*s3
 b(k3) = akm*s3 + bkm*c3
 kk = k3 + kspan
 if (kk<=nt) go to 420
 440 c2 = cd*c1 + sd*s1
 s1 = (sd*c1-cd*s1) + s1
 c1 = c1 - c2
 c2 = c1*c1 - s1*s1
 s2 = 2 * c1 * s1
 c3 = c2*c1 - s2*s1
 s3 = c2*s1 + s2*c1
 kk = kk - nt + jc
 if (kk<=kspan) go to 420
 kk = kk - kspan + inc
 if (kk<=jc) go to 410
 if (kspan==jc) go to 800
 go to 100
 450 akp = akm + bjm
 akm = akm - bjm
 bkp = bkm - ajm
 bkm = bkm + ajm
 if (s1/=0) go to 430
 460 a(k1) = akp
 b(k1) = bkp
 a(k2) = ajp
 b(k2) = bjp
 a(k3) = akm
 b(k3) = bkm
 kk = k3 + kspan
 if (kk<=nt) go to 420
 go to 440
 !  transform for factor of 5 (optional code)
 510 c2 = c72*c72 - s72*s72
 s2 = 2 * c72 * s72
 520 k1 = kk + kspan
 k2 = k1 + kspan
 k3 = k2 + kspan
 k4 = k3 + kspan
 akp = a(k1) + a(k4)
 akm = a(k1) - a(k4)
 bkp = b(k1) + b(k4)
 bkm = b(k1) - b(k4)
 ajp = a(k2) + a(k3)
 ajm = a(k2) - a(k3)
 bjp = b(k2) + b(k3)
 bjm = b(k2) - b(k3)
 aa = a(kk)
 bb = b(kk)
 a(kk) = aa + akp + ajp
 b(kk) = bb + bkp + bjp
 ak = akp*c72 + ajp*c2 + aa
 bk = bkp*c72 + bjp*c2 + bb
 aj = akm*s72 + ajm*s2
 bj = bkm*s72 + bjm*s2
 a(k1) = ak - bj
 a(k4) = ak + bj
 b(k1) = bk + aj
 b(k4) = bk - aj
 ak = akp*c2 + ajp*c72 + aa
 bk = bkp*c2 + bjp*c72 + bb
 aj = akm*s2 - ajm*s72
 bj = bkm*s2 - bjm*s72
 a(k2) = ak - bj
 a(k3) = ak + bj
 b(k2) = bk + aj
 b(k3) = bk - aj
 kk = k4 + kspan
 if (kk<nn) go to 520
 kk = kk - nn
 if (kk<=kspan) go to 520
 go to 700
 !  transform  for  odd  factors
 600 k = nfac(i)
 kspnn = kspan
 kspan = kspan / k
 if (k==3) go to 320
 if (k==5) go to 510
 if (k==jf) go to 640
 jf = k
 s1 = rad / k
 c1 = cos(s1)
 s1 = sin(s1)
 if (jf>maxf) go to 996
 ck(jf) = 1.d0
 sk(jf) = 0.d0
 j = 1
 630 ck(j) = ck(k)*c1 + sk(k)*s1
 sk(j) = ck(k)*s1 - sk(k)*c1
 k = k - 1
 ck(k) = ck(j)
 sk(k) = -sk(j)
 j = j + 1
 if (j<k) go to 630
 640 k1 = kk
 k2 = kk + kspnn
 aa = a(kk)
 bb = b(kk)
 ak = aa
 bk = bb
 j = 1
 k1 = k1 + kspan
 650 k2 = k2 - kspan
 j = j + 1
 at(j) = a(k1) + a(k2)
 ak = at(j) + ak
 bt(j) = b(k1) + b(k2)
 bk = bt(j) + bk
 j = j + 1
 at(j) = a(k1) - a(k2)
 bt(j) = b(k1) - b(k2)
 k1 = k1 + kspan
 if (k1<k2) go to 650
 a(kk) = ak
 b(kk) = bk
 k1 = kk
 k2 = kk + kspnn
 j = 1
 660 k1 = k1 + kspan
 k2 = k2 - kspan
 jj = j
 ak = aa
 bk = bb
 aj = 0.d0
 bj = 0.d0
 k = 1
 670 k = k + 1
 ak = at(k)*ck(jj) + ak
 bk = bt(k)*ck(jj) + bk
 k = k + 1
 aj = at(k)*sk(jj) + aj
 bj = bt(k)*sk(jj) + bj
 jj = jj + j
 if (jj>jf) jj = jj - jf
 if (k<jf) go to 670
 k = jf - j
 a(k1) = ak - bj
 b(k1) = bk + aj
 a(k2) = ak + bj
 b(k2) = bk - aj
 j = j + 1
 if (j<k) go to 660
 kk = kk + kspnn
 if (kk<=nn) go to 640
 kk = kk - nn
 if (kk<=kspan) go to 640
 !  multiply by rotation factor (except for factors of 2 and 4)
 700 if (i==m) go to 800
 kk = jc + 1
 710 c2 = 1 - cd
 s1 = sd
 720 c1 = c2
 s2 = s1
 kk = kk + kspan
 730 ak = a(kk)
 a(kk) = c2*ak - s2*b(kk)
 b(kk) = s2*ak + c2*b(kk)
 kk = kk + kspnn
 if (kk<=nt) go to 730
 ak = s1 * s2
 s2 = s1*c2 + c1*s2
 c2 = c1*c2 - ak
 kk = kk - nt + kspan
 if (kk<=kspnn) go to 730
 c2 = c1 - (cd*c1+sd*s1)
 s1 = s1 + (sd*c1-cd*s1)
 kk = kk - kspnn + jc
 if (kk<=kspan) go to 720
 kk = kk - kspan + jc + inc
 if (kk<=jc+jc) go to 710
 go to 100
 !  permute results to normal order--- done in two stages
 !  permutation for square factors of n
 800 np(1) = ks
 if (kt==0) go to 890
 k = kt + kt + 1
 if (m<k) k = k - 1
 j = 1
 np(k+1) = jc
 810 np(j+1) = np(j) / nfac(j)
 np(k) = np(k+1) * nfac(j)
 j = j + 1
 k = k - 1
 if (j<k) go to 810
 k3 = np(k+1)
 kspan = np(2)
 kk = jc + 1
 k2 = kspan + 1
 j = 1
 if(n/=ntot) goto 850
 !  permutation for single-variate transform (optional code)
 820 ak = a(kk)
 a(kk) = a(k2)
 a(k2) = ak
 bk = b(kk)
 b(kk) = b(k2)
 b(k2) = bk
 kk = kk + inc
 k2 = kspan + k2
 if (k2<ks) go to 820
 830 k2 = k2 - np(j)
 j = j + 1
 k2 = np(j+1) + k2
 if (k2>np(j)) go to 830
 j = 1
 840 if (kk<k2) go to 820
 kk = kk + inc
 k2 = kspan + k2
 if (k2<ks) go to 840
 if (kk<ks) go to 830
 jc = k3
 go to 890
 !  permutation for multivariate transform
 850 k = kk + jc
 860 ak = a(kk)
 a(kk) = a(k2)
 a(k2) = ak
 bk = b(kk)
 b(kk) = b(k2)
 b(k2) = bk
 kk = kk + inc
 k2 = k2 + inc
 if (kk<k) go to 860
 kk = kk + ks - jc
 k2 = k2 + ks - jc
 if (kk<nt) go to 850
 k2 = k2- nt + kspan
 kk = kk- nt + jc
 if (k2<ks) go to 850
 870 k2 = k2 - np(j)
 j = j + 1
 k2 = np(j+1) + k2
 if (k2>np(j)) go to 870
 j = 1
 880 if (kk<k2) go to 850
 kk = kk + jc
 k2 = kspan + k2
 if (k2<ks) go to 880
 if (kk<ks) go to 870
 jc = k3
 890 if (2*kt+1>=m) return
 kspnn = np(kt+1)
 !  permutation for square free factors of n
 j = m - kt
 nfac(j+1) = 1
 900 nfac(j) = nfac(j) * nfac(j+1)
 j = j - 1
 if (j/=kt) go to 900
 kt = kt + 1
 nn = nfac(kt) - 1
 if (nn>maxp) go to 998
 jj = 0
 j = 0
 go to 906
 902 jj = jj - k2
 k2 = kk
 k = k + 1
 kk = nfac(k)
 904 jj = kk + jj
 if (jj>=k2) go to 902
 np(j) = jj
 906 k2 = nfac(kt)
 k = kt + 1
 kk = nfac(k)
 j = j + 1
 if (j<=nn) go to 904
 !  determine the permutation cycles of length greter than 1
 j = 0
 go to 914
 910 k = kk
 kk = np(k)
 np(k) = -kk
 if (kk/=j) go to 910
 k3 = kk
 914 j = j + 1
 kk = np(j)
 if (kk<0) go to 914
 if (kk/=j) go to 910
 np(j) = -j
 if (j/=nn) go to 914
 maxf = inc * maxf
 !  reorder a and b, following the permutation cycles
 go to 950
 924 j = j - 1
 if (np(j)<0) go to 924
 jj = jc
 926 kspan = jj
 if (jj>maxf) kspan = maxf
 jj = jj - kspan
 k = np(j)
 kk = jc*k + ii + jj
 k1 = kk + kspan
 k2 = 0
 928 k2 = k2 + 1
 at(k2) = a(k1)
 bt(k2) = b(k1)
 k1 = k1 - inc
 if (k1/=kk) go to 928
 932 k1 = kk + kspan
 k2 = k1 - jc*(k+np(k))
 k = -np(k)
 936 a(k1) = a(k2)
 b(k1) = b(k2)
 k1 = k1 - inc
 k2 = k2 - inc
 if (k1/=kk) go to 936
 kk = k2
 if (k/=j) go to 932
 k1 = kk + kspan
 k2 = 0
 940 k2 = k2 + 1
 a(k1) = at(k2)
 b(k1) = bt(k2)
 k1 = k1 - inc
 if (k1/=kk) go to 940
 if (jj/=0) go to 926
 if (j/=1) go to 924
 950 j = k3 + 1
 nt = nt - kspnn
 ii = nt - inc + 1
 if (nt>=0) go to 924
 return
 
 ! error finish, insufficient array storage
 996 continue
 write(*,'("fft.f95: the array dimension parameter maxf is insufficient")')
 ierr=1
 
 998 continue
 write(*,'("fft.f95: the array dimension parameter maxp is insufficient")')
 ierr=1
 
 return
end subroutine fft

!*********************************************
!Calculation of weighted radial amplitude spectrum
!*********************************************
subroutine f_sras(outfile,rfg,ifg,a,n,n2,n1,n22,dy,dx,dw,wmax)

      
     character*80 outfile
     real,dimension(n1,n2)::rfg,ifg
     dimension a(1000),n(1000)

      ncp1=n2+1
      cn=1./(n2*dy)
      rn=1./(n1*dx)

      scale=float(n1)*float(n2)
      nrnq=float(n1)/2.+1.0000001
      ncnq=float(n2)/2.+1.0000001
      pi=3.1415927e0


      dwh=0.5*dw
      rdx=1./dx
      rdy=1./dy
      rmin=rdx
      if(rdy.lt.rdx)rmin=rdy
      if(2.*wmax.ge.rmin)go to 900
      nrowmax=1+wmax/rn
      ncolmax=1+wmax/cn
    

      open(15,file=outfile,status='unknown',form='formatted',position='append')   
      write(15,50)
 50   format(6x,' w',4x,'log mean |F|')
      do 11 kw=1,1000
      a(kw)=0.0
 11    n(kw)=0
  
 ! Read fftfil.coef file and modify.
 !  First row, DC term treated separately.
     
      a(1)=a(1)+sqrt(rfg(1,1)*rfg(1,1)+ifg(1,1)*ifg(1,1))
      n(1)=n(1)+1

      do 1 i=2,n2
      !ic=2*i
      ii=i-1
      if(i.gt.ncnq) ii=-(ncp1-i)
      v2=float(ii)*cn
      ! (If necessary, here consider effect of sign of v2.)
      w=abs(v2)
      Gr=rfg(1,i)
      Gi=ifg(1,i)
      if(w.gt.wmax)go to 1
      kw=1+w/dw
      a(kw)=a(kw)+sqrt(Gr*Gr+Gi*Gi)
      n(kw)=n(kw)+1
  1     continue

     ! Read rest of fftfil.cof file, row-wise, and modify.
      do 3 jr=2,nrowmax
      jj=jr-1
      v1=float(jj)*rn
      v1sq=v1*v1
      !read(14,rec=jr)g
      do 2 i=1,n2
      !ic=2*i
      ii=i-1
      if(i.gt.ncnq) ii=-(ncp1-i)
      v2=float(ii)*cn
      v2sq=v2*v2
      w=sqrt(v1sq+v2sq)
     !(If necessary, consider effect of sgn of v1, v2.)
      Gr=rfg(jr,i)
      Gi=ifg(jr,i)
      if(w.gt.wmax)go to 2
      kw=1+w/dw
      a(kw)=a(kw)+sqrt(Gr*Gr+Gi*Gi)
     
      n(kw)=n(kw)+1
  2     continue
  3     continue

      c=1./(n1*n2)
      kwmax=1+wmax/dw
      do 12 kw=1,kwmax
      rkw=kw
      w=rkw*dw-dwh
      aras=a(kw)*c
      
      rasm=-1
      rasmlog=-1
      if(n(kw).ne.0)rasm=aras/n(kw)
      if(rasm.gt.0) rasmlog=log(rasm)
      write(15,13)w,rasmlog
 13    format(1x,2e13.5,i5,2e13.5)
 12    continue
      go to 999
 900   write(6,901)
 901   format(' Cant handle it.')
 999   close(15)
      
      return
end subroutine
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
!!From file channel number Nunit_ Start starts 
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
