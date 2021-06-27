!*****************************Depth.F90 *****************************
!
!This Program reads Z0 file and Zt file and their uncertainty name from LS_Fitting and calculate Zb according format Zb=2*Z0-Zt.
!then read coordinate file from FFT to generate the xyz data including (x,y,Zb)
!
!
!  Written by:
!  Yihong Yin, Chun-Feng LI,
!
!
! usage: DepthZ0 input file, Zt input file, Z0 uncertainty input file, Zt uncertainty input file, Zb output file, Zb uncertainty output file
!
! For further information about program, use command "man Depth " in Terminal.
!


Program Depth
 implicit none

  
  character*80 infile1,infile2,infile3,infile4,outfile1,outfile2

  integer::columns,rows
  real,allocatable,dimension(:)::data1,data2,data3,data4
  real,allocatable,dimension(:,:)::coordinate
  real,allocatable,dimension(:)::data,uncertainty
  integer::i,j,nunit_in,i1,i2,i3

  integer::narg

  

  narg=IARGC()
  if(narg==6) then 
   call getarg(1,infile1)
   call getarg(2,infile2)
   call getarg(3,infile3)
   call getarg(4,infile4)
   call getarg(5,outfile1)
   call getarg(6,outfile2)
  else
   write(*,*) 'pls input correct parameter!'
  endif

   !i1=index(infile1,'.')
   !infile11=infile1(1:i1-1)
   !i2=index(infile1,'.')
   !infile22=infile2(1:i2-1)
   !i3=index(outfile,'.')
   !outfile11=outfile(1:i3-1)
  !infile1='Z0_my.txt'
  !infile2='Zt_my.txt'
  !infile3='coordinate.dat'
  !outfile='o.txt'

  write(*,*)'Depth is running, pls waiting ...'
  call search_unit(10,nunit_in)
  call open_old_file(nunit_in,infile1,'formatted')
  call open_old_file(nunit_in+1,infile2,'formatted')
  call open_old_file(nunit_in+2,infile3,'formatted')
  call open_old_file(nunit_in+3,infile4,'formatted')
  call open_old_file(nunit_in+4,'coordinate.dat','formatted')
  


  read(nunit_in,*) columns,rows
  
  read(nunit_in+1,*) 
 

  allocate(data1(rows*columns),data2(rows*columns),data3(rows*columns),data4(rows*columns),uncertainty(rows*columns))
  allocate(data(rows*columns),coordinate(rows*columns,2))

  do i=1,rows*columns
    read(nunit_in,*) data1(i)
    read(nunit_in+1,*)data2(i)
    read(nunit_in+2,*) data3(i)
    read(nunit_in+3,*) data4(i)
    read(nunit_in+4,*)(coordinate(i,j),j=1,2)
  end do
  
  do i=1,rows*columns
    data(i)=2*data1(i)-data2(i)
     if(data(i).lt.0) data(i)=0
    uncertainty(i)=sqrt(4*data3(i)*data3(i)+ data4(i)*data4(i));
  end do
  
  deallocate(data1,data2,data3,data4)

  open(13,file=outfile1,status='unknown')
  open(14,file=outfile2,status='unknown')
  
  !write(13,"(a4)") 'DSAA'
  !write(13,*) columns,rows
  !write(13,*) xmin,xmax
  !write(13,*) ymin,ymax
  !write(13,*) minval(data),maxval(data)

  do i=1,rows*columns
    write(13,*) (coordinate(i,j),j=1,2),data(i)
    write(14,*) (coordinate(i,j),j=1,2),uncertainty(i)
  end do

  close(13)
  close(14)
  deallocate(data,coordinate,uncertainty)

  write(*,*)'Depth has done successfully...'
end Program Depth

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
   
   
   
