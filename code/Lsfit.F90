!***************************** LS_Fitting.F90 *****************************
!
!This program reads the spectra file from FFT,then calculate the Z0 and Zt in low wavenumber domain and high wavenumber domain 
!which you select respectively by least square fitting method.
!
!
! Written by:
!  Yihong Yin, Chun-Feng LI,
!
! 
! usage: LS_Fitting input file name, Zt/Z0,Fractal exponent,line of each spectrum,kmin,kmax output file name
!
! For further information about program, use command "man LS_Fitting " in Terminal.
!

program Lsfit


 implicit none
  

  
  real ,parameter::PI=3.141592654
  character*80 infile,outfile,Z,linestr,Bstr,index1str,index2str,outfile2
  integer::line
  real::B!B is Fractal exponent
  integer::index1,index2!index of kmin and kmax
 
  integer::nunit_in,columns,rows
  real,allocatable,dimension(:,:)::data
  
 
  real,allocatable,dimension(:)::k,rasmlog 
  real,allocatable,dimension(:)::k0,rasmlog0,k1,rasmlog1 
  real,allocatable,dimension(:)::zt,uncertainty
  integer::i,j,m,n,length
  integer*8 i1
  character*80::outfile1

  integer::narg
  narg=IARGC()

 if(narg==8) then 

    call getarg(1,infile)
    call getarg(2,Z)
    call getarg(3,Bstr)
    call getarg(4,linestr)
    call getarg(5,index1str)
    call getarg(6,index2str)
    call getarg(7,outfile)
    call getarg(8,outfile2)
  else
   write(*,*) 'pls input correct parameter!'
  endif
  
  read(Bstr,*) B
  read(linestr,*) line
  read(index1str,*) index1
  read(index2str,*) index2

  !infile='spectra_cwt_10.673.dat'
  !Z='Zt'
  !line=43
  !index1=32
  !index2=37
  !outfile='Zt_10.673.txt'
  !outfile2='Zt_10.673_un.txt'
  !B=3

    !length=repeat(outfile,'.')
     !if(length.eq.q) then
     !  i1=index(outfile,'.')
      ! outfile1=outfile(1:i1-1)
    !else
      ! 
     !end if

  write(*,*) 'The least square fitting is running,pls waiting...'
  
  call search_unit(10,nunit_in)
  call open_old_file(nunit_in,infile,'formatted')
  
  
  read(nunit_in,*) columns,rows
 
  
  allocate(zt(rows*columns))
  allocate(uncertainty(rows*columns))

  do m=1,rows*columns
  
     allocate(data(line-1,2))
     allocate(k(line-1),rasmlog(line-1))
     
     do i=1,line
       if (i.eq.1) then
         read(nunit_in,*)
       else
         read(nunit_in,*) (data(i-1,j),j=1,2)
       end if
      end do 
     k=data(:,1)! read wavenumber
     rasmlog=data(:,2)!read Amplitude spectrum
     
    
      allocate(k0(index1:index2),rasmlog0(index1:index2),k1(index1+1:index2+1),rasmlog1(index1+1:index2+1))
      k0=k(index1:index2)
      rasmlog0=rasmlog(index1:index2)

     
      k1=k(index1+1:index2+1)
      rasmlog1=rasmlog(index1+1:index2+1)

      if(Z.eq.'Zt') then
        call least_square(k0*2*PI,rasmlog0+(B-1)/2*log(k0*2*PI),index1,index2,zt(m),uncertainty(m))
      else if(Z.eq.'Z0') then
        if((rasmlog(index1+1)+(B-1)/2*log(k(index1+1)*2*PI)-log(k(index1+1)*2*PI))&
            .gt.(rasmlog(index1)+(B-1)/2*log(k(index1)*2*PI)-log(k(index1)*2*PI))) then
		!the first point value more than the following points
         call least_square(k1*2*PI,rasmlog1+(B-1)/2*log(k1*2*PI)-log(k1*2*PI),index1+1,index2+1,zt(m),uncertainty(m))
        else
		!the first point value more than the following points 
         call least_square(k0*2*PI,rasmlog0+(B-1)/2*log(k0*2*PI)-log(k0*2*PI),index1,index2,zt(m),uncertainty(m))
        end if
      else
        write(*,*) 'please enter the second correct parameters'
      end if
      deallocate(k0,rasmlog0,k1,rasmlog1,k,rasmlog,data)

    
  end do 

  open(13,file=outfile,status='unknown')
  write(13,*) columns,rows

  open(14,file=outfile2,status='unknown')
  write(14,*) columns,rows
 

  do i=1,rows*columns
    write(13,*) (abs(zt(i)))
    write(14,*) (abs(uncertainty(i)))
  end do
  close(13)
  close(14)
  
  deallocate(zt,uncertainty)

  write(*,*) 'The least square fitting has done successfully'
 
  
end program Lsfit


!*****************************************************
! The least square fitting
!*****************************************************
subroutine  least_square(x,y,index1,index2,z,s)
 

 real::A,B,C,D,E,F,y_new,y_tmp,x_tmp,s,Lxx
 real::x(index2-index1+1),y(index2-index1+1)
 integer::i,j,n
 real::z
 real::temp,tempSumXX,tempSumYY,Xmean,Ymean,b_tmp
 

 tempSumYY=0
 n=index2-index1+1

 A=0;B=0;C=0;D=0;E=0;F=0
 y_tmp=0;x_tmp=0;

 do i=1,n
   A=A+x(i)*x(i)
   B=B+x(i)
   C=C+x(i)*y(i)
   D=D+y(i)
 end do

 temp=n*A-B*B
 if(temp.ne.0) then
   z=(n*C-B*D)/temp !calculate slope
 else
   z=0
 end if
  if(z.gt.0) z=0

  Xmean=B/n
  Ymean=D/n
  b_tmp= Ymean-z*Xmean

 !calculate uncertainty
  
 do i=1,n
   
    y_new= z*x(i)+b_tmp
    y_tmp=y_tmp+(y_new-y(i))*(y_new-y(i))
    x_tmp=x_tmp+(Xmean-x(i))*(Xmean-x(i))
 end do

 if(x_tmp.ne.0) then
  s=sqrt((y_tmp/x_tmp)/(n-2))
 else
  s=0
  end if 

 
 
end subroutine least_square

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


