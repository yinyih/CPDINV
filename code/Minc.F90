!***************************** MinCurV.F90 *****************************
!
!This program generates a 2-dimensional grid, equally incremente in x and y, from randomly placed data points. 
!The minimum curvature algorithm  produces a smooth grid by iteratively solving a set of difference equations which minimize the total
!2nd horizontal derivative and attempt to honor input data.,
!the output file header includes 'DSAA';columns rows ;xmin xmax;ymin ymax;zmin zmax
!
! Written by:
!  Yihong Yin, Chun-Feng LI,
! 
!
! ReferencesS
!  Wang W Y. (2010). Minimum Curvature method and Fortran programming in data processing of potential field (in Chinese). Geological Publishing House.
! 
!
! usage: MinCurV input file name, xmin/xmax/ymin/ymax, dx,dy,output file name
!
! For further information about program, use command "man MinCurV" in Terminal.
!

Program Minc
 
 implicit none

 
 character*80 input_file,range,dx_str,dy_str,output_file !Necessary parameters
 character*80  duplicate_str,edge_str,restrict_str,iteration_max_str,eps_max_str!Optional parameters
 character*80 xmin_str,xmax_str,ymin_str,ymax_str
 
 integer::column_x,column_y,column_u,column,mpoint_random,mpoint,line
 integer::iteration_max
 integer::duplicate,method_edge,method_restrict,method_initial
 integer::i1,i2,i3
 real::xmin,xmax,ymin,ymax,dx,dy
 real::eigval,eps_max
 integer*2::number_arg,number
 logical flage !
 real,allocatable::FIELD(:,:),FIELD_x(:),FIELD_y(:),FIELD_u(:)
 integer::ierr=0

 integer::narg
 narg=IARGC()
  if (narg==5) then 
   call getarg(1,input_file) 
   call getarg(2,range)
   call getarg(3,dx_str)
   call getarg(4,dy_str)
   call getarg(5,output_file)
   
   write(*,*)'The default parameters are as followed'
   write(*,*)'the iteration_max is 100000'
   write(*,*) 'the eps_max is 0.005'
   write(*,*)'the duplicate is 0 '
   write(*,*)'the method_edge is 1'
   write(*,*) 'the method_restrict is 2 '
   
    iteration_max=100000
    eps_max=0.005
    duplicate=0
    method_edge=1
    method_restrict=2
  else if(narg==10) then
   call getarg(1,input_file) 
   call getarg(2,range)
   call getarg(3,dx_str)
   call getarg(4,dy_str)
   call getarg(5,iteration_max_str)
   call getarg(6,eps_max_str)
   call getarg(7,duplicate_str)
   call getarg(8,edge_str)
   call getarg(9,restrict_str)
   call getarg(10,output_file)
  read(iteration_max_str,*) iteration_max
  read(eps_max_str,*) eps_max
  read(duplicate_str,*) duplicate
  read(edge_str,*) method_edge
  read(restrict_str,*) method_restrict
 
  else
    write(*,*)'pls input correct parameters'
  endif
  

 i1=index(range,'/')
 xmin_str=range(1:i1-1)
 i2=index(range(i1+1:),'/')+i1
 xmax_str=range(i1+1:i2-1)
 i3=index(range(i2+1:),'/')+i2
 ymin_str=range(i2+1:i3-1)
 ymax_str=range(i3+1:)
  
  read(xmin_str,*) xmin
  read(xmax_str,*) xmax
  read(ymin_str,*) ymin
  read(ymax_str,*) ymax
  read(dx_str,*) dx
  read(dy_str,*) dy
  
 write(*,*)'Minimum Curvature is running, pls waiting ...'
 
 
 eigval=1.70141e+38
 column_x=1
 column_y=2
 column_u=3
 method_initial=4
 mpoint=(xmax-xmin)/dx+1
 line=(ymax-ymin)/dy+1

 column=max(column_x,column_y,column_u)
 call get_xyz_mpoint(input_file,mpoint_random,column)
  IF (mpoint_random.gt.0) then
  
   allocate (FIELD_x(1:mpoint_random),FIELD_y(1:mpoint_random),FIELD_u(1:mpoint_random),STAT=ierr)
   allocate (FIELD(1:mpoint,1:line),STAT=ierr)
   IF(ierr.eq.0) then
   call input_xyz_column(input_file,mpoint_random,FIELD_x,FIELD_y,FIELD_u,column_x,column_y,column_u)
 
   call MinCurV_2D_random_grid_sub(eigval,eps_max,iteration_max,duplicate,method_edge,method_restrict,method_initial,&
                                mpoint_random,FIELD_x,FIELD_y,FIELD_u,xmin,xmax,mpoint,ymin,ymax,line,FIELD,flage)
   call OUTPUT_GRD(FIELD,output_file,mpoint,line,xmin,xmax,ymin,ymax,1,mpoint,1,line,eigval)
   end IF
   deallocate(FIELD_x,FIELD_y,FIELD_u,FIELD)

   write(*,*)'Minimum Curvature has done successfully'
   
  end IF

end Program Minc


!*******************************************************
!Read parameters for 2D discrete data gridding
!********************************************************

subroutine READCMD_2D_grid(CMDFILE,input_file,output_file,xmin,xmax,mpoint,ymin,ymax,line,eigval,eps_max,iteration_max,&
                            duplicate, method_edge,method_restrict,column_x,column_y,column_u)
    
    character*(*) cmdfile,input_file,output_file
    integer mpoint,line,duplicate,column_x,column_y,column_u
    real xmin,xmax,ymin,ymax,eigval,eps_max
    
    integer nunit_in
   
    call search_unit(0,nunit_in)   
    call open_old_file(nunit_in,cmdfile,'formatted')
    read(nunit_in,*)
    call read_char(nunit_in,input_file)  
    call read_char(nunit_in,output_file) 
    call read_real_1(nunit_in,xmin)
    call read_real_1(nunit_in,xmax)
    call read_integer(nunit_in,mpoint)
    call read_real_1(nunit_in,ymin)
    call read_real_1(nunit_in,ymax)
    call read_integer(nunit_in,line)
    read(nunit_in,*)
    call read_integer(nunit_in,duplicate)
    call read_integer(nunit_in,method_edge)
    call read_integer(nunit_in,method_restrict)
    
    call read_real_1(nunit_in,eigval)
    call read_integer(nunit_in,iteration_max)
    call read_real_1(nunit_in,eps_max)
    call read_integer(nunit_in,column_x)
    call read_integer(nunit_in,column_y)
    call read_integer(nunit_in,column_u)
  
    close(nunit_in)


end subroutine READCMD_2D_grid

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

!********************************************
!In the open file channel number nunit_ In, read a string 'char'after = i
!********************************************
subroutine read_char(nunit_in,char) 
 
 character*(*)char
 character*200 ctemp

 indexposi=0
 do while(indexposi==0)
   read(nunit_in,'(a)') ctemp
  indexposi=index(ctemp,'=')
 end do !while
 char=ctemp(indexposi+1:)
 char=trim(char)
 
end subroutine read_char

!********************************************
!In the open file channel number nunit_ In, read a real number 'r1' after = 
!********************************************
subroutine read_real_1(nunit_in,r1) 

 character*200 ctemp,rr1
 indexposi=0

 do while(indexposi==0)
   read(nunit_in,'(a)') ctemp
  indexposi=index(ctemp,'=')
 end do !while
 rr1=ctemp(indexposi+1:)
 rr1=trim(rr1)
 read(rr1,*) r1
 
end subroutine read_real_1

!**********************************************
!In the open file channel number nunit_ In, read a integer number 'm'after = i
!**********************************************
subroutine read_integer(nunit_in,m)
 
 
 character*200 ctemp,rr1
 indexposi=0

 do while(indexposi==0)
   read(nunit_in,'(a)') ctemp
  indexposi=index(ctemp,'=')
 end do !while
 rr1=ctemp(indexposi+1:)
 rr1=trim(rr1)
 read(rr1,*) m  

end subroutine read_integer

!********************************************
!Get the number of discrete pointsmpoint_random
!********************************************
subroutine get_xyz_mpoint(input_file,mpoint_random,column)
 
 implicit none
 character*(*) input_file
 integer::nunit_in,column,status1,mpoint_random,i
 real temp
 mpoint_random=0

 call search_unit(10,nunit_in)
 call open_old_file(nunit_in,input_file,'formatted')

 do while(.true.)
  read(nunit_in,*,iostat = status1)(temp,i=1,column)
  if(status1 /= 0) exit
  mpoint_random=mpoint_random+1
 end do 
       
 close(nunit_in)  
end subroutine get_xyz_mpoint

!************************************************
!Read three columns of XYZ data
!***********************************************
subroutine input_xyz_column(input_file,mpoint,x_point,y_point,z_point,column_x,column_y,column_z)
 
 
 integer::column_x,column_y,column_z,mpoint
 character*(*) input_file
 real x_point(mpoint),y_point(mpoint),z_point(mpoint)
 integer nunit_in
 real,allocatable::temp(:)
 integer::column,i,k
 
 column=max(column_x,column_y,column_z)
 allocate(temp(column))

 call search_unit(10,nunit_in)
 call open_old_file(nunit_in,input_file,'formatted')
 do i=1,mpoint
  read(nunit_in,*) (temp(k),k=1,column)
  x_point(i)=temp(column_x)
  y_point(i)=temp(column_y)
  z_point(i)=temp(column_z)
 
 
 end do

  close(nunit_in)
  deallocate(temp)

end subroutine input_xyz_column

!************************************************
!Minimum curvature gridding of two dimensional discrete data
!The sub modules used are as follows：
!duplicate_2d_eigval_sub:eliminate 2-D repeated point 
!method_inital_2d_random_sub:determinate initial value and constraint point selection of 2-D grid nodes
!method-edge_2d_eigval_sub:determinate endpoint of 2-D grid data
!method_inital_2d_sub:Determinate initial value of blank point in 2-D grid data
!MinCurV_2D_random_sub:Minimum curvature iteration for constrained data on 2-D meshes
!***********************************************
subroutine MinCurV_2D_random_grid_sub(eigval,eps_max,iteration_max,duplicate,method_edge,method_restrict,method_initial,&
                                mpoint_random,FIELD_x,FIELD_y,FIELD_u,xmin,xmax,mpoint,ymin,ymax,line,FIELD,flage)
  
  real eigval,eps_max,xmin,xmax,ymin,ymax
  integer iteration_max,mpoint_random,mpoint,line
  integer method_edge,method_restrict,method_initial,duplicate
  real FIELD(1:mpoint,1:line)
  real FIELD_x(1:mpoint_random),FIELD_y(1:mpoint_random),FIELD_u(1:mpoint_random)
  logical flage

  real,allocatable::DELT_x(:,:),DELT_y(:,:),DELT_u(:,:)
  integer,allocatable::M_eigval(:),N_eigval(:)
  
 
  !Deal with repeated points, eliminate eigenvalues and data outside the grid range
  call duplicate_2d_eigval_sub(eigval,duplicate,mpoint_random,FIELD_x,FIELD_y,FIELD_u,mpoint,xmin,xmax,line,ymin,ymax)
  
  
  
  !Determine the initial values of grid nodes and constraint points
  allocate(DELT_x(1:mpoint,1:line),DELT_y(1:mpoint,1:line),DELT_u(1:mpoint,1:line),STAT=IERR)
  call method_inital_2d_random_sub(eigval,mpoint_random,FIELD_x,FIELD_y,FIELD_u,xmin,xmax,mpoint,&
                                    ymin,ymax,line,FIELD,DELT_x,DELT_y,DELT_u)

  !Find the number of blank points number_eigval
  number_eigval=0
  do j=1,line
   do i=1,mpoint
    if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
     number_eigval=number_eigval+1
    end if
   end do
  end do
  
  allocate(M_eigval(number_eigval),N_eigval(number_eigval),STAT=ierr)
  
  !Give the value of the end of the edge,
  !and find the blank point to save it in the array
  call METHOD_edge_2D_eigval_sub(method_edge,method_restrict,eigval,mpoint,line,FIELD,number_eigval,M_eigval,N_eigval)

  !Convert the maximum relative error (eps_max) to maximum  absolute error (eps_abs)
  fmin=huge(fmin)
  fmax=-huge(fmax)
  DO j=1,line
   DO i=1,mpoint
    if(abs(FIELD(i,j)).lt.abs(eigval*0.9)) then
     fmin=min(fmin,FIELD(i,j))
     fmax=max(fmax,FIELD(i,j))
    end if
   end DO
  end DO
  eps_abs=abs(eps_max*abs(fmax-fmin))

  !Assign initial value to blank point
  call method_inital_2d_sub(method_initial,eigval,mpoint,line,FIELD)

  !Do the constrained minimum curvature iteration of 2-D grid data
  call MinCurV_2D_random_sub(eps_abs,iteration_max,eigval,DELT_x,DELT_y,DELT_u,&
                              xmin,xmax,mpoint,ymin,ymax,line,FIELD,flage,&
                               number_eigval,M_eigval,N_eigval)
  deallocate(DELT_x,DELT_y,DELT_u,M_eigval,N_eigval)


end subroutine MinCurV_2D_random_grid_sub

!*****************************************************
!eliminate 2-D repeated point and eigenvalue
!*****************************************************
subroutine duplicate_2d_eigval_sub(eigval,duplicate,mpoint_random,FIELD_x,FIELD_y,FIELD_u,mpoint,xmin,xmax,line,ymin,ymax)
 
 real eigval,xmin,xmax,ymin,ymax
 integer duplicate,mpoint_random,mpoint,line
 real FIELD_x(1:mpoint_random),FIELD_y(1:mpoint_random),FIELD_u(1:mpoint_random)

 !eliminate points with function values as eigenvalues
 number=0
 do i=1,mpoint_random
  if(abs(FIELD_x(i)).lt.abs(eigval*0.9)) then 
   if(abs(FIELD_y(i)).lt.abs(eigval*0.9)) then 
     if(abs(FIELD_u(i)).lt.abs(eigval*0.9)) then 
     number=number+1
     FIELD_x(number)=FIELD_x(i)
     FIELD_y(number)=FIELD_y(i)
     FIELD_u(number)=FIELD_u(i)
     end if
   end if
   end if
 end do
 mpoint_random=number

 !Remove the data out of the gridding range
  dx=(xmax-xmin)/(mpoint-1)
  dy=(ymax-ymin)/(line-1)
  dxy=max(dx,dy)
  xxmin=xmin-dxy
  xxmax=xmax+dxy
  yymin=ymin-dxy
  yymax=ymax+dxy
  number=0

  do i=1,mpoint_random
   if(FIELD_x(i).lt.xxmax.and.(FIELD_x(i).gt.xxmin)) then
    if(FIELD_y(i).lt.yymax.and.(FIELD_y(i).gt.yymin)) then
     number=number+1
     FIELD_x(number)=FIELD_x(i)
     FIELD_y(number)=FIELD_y(i)
     FIELD_u(number)=FIELD_u(i)
    end if
   end if
  end do
  mpoint_random=number

  !Sort the data according to the coordinates from large to small
  call shell_3real(FIELD_x,FIELD_y,FIELD_u,mpoint_random,1,mpoint_random)
  
  !Eliminate duplicate data
  eps_x=epsilon(eps_x)*dx
  eps_y=epsilon(eps_y)*dy
  k=1
  do while((k.lt.mpoint_random).and.(FIELD_u(k).lt.eigval*0.9))
   L=K+1
    number=1
    duplicate_x=FIELD_x(k)
    duplicate_y=FIELD_y(K)
    duplicate_u=FIELD_u(K)
    do while((l.le.mpoint_random).and.(abs(FIELD_x(k)-FIELD_x(L)).le.eps_x))
     if(abs(FIELD_y(k)-FIELD_y(L)).le.eps_y) then
      if(duplicate.eq.-1) then
       number=1
       duplicate_x=duplicate_x+FIELD_x(L)
       duplicate_y=duplicate_x+FIELD_y(L)
       duplicate_u=min(duplicate_u,FIELD_u(L))
      else if(duplicate.eq.1) then
       number=1
       duplicate_x=duplicate_x+FIELD_x(L)
       duplicate_y=duplicate_x+FIELD_y(L)
       duplicate_u=min(duplicate_u,FIELD_u(L))
      else if(duplicate.eq.0) then
       number=number+1
       duplicate_x=duplicate_x+FIELD_x(L)
       duplicate_y=duplicate_x+FIELD_y(L)
       duplicate_u=min(duplicate_u,FIELD_u(L))
      else
      end if

      FIELD_x(L)=eigval
      FIELD_y(L)=eigval
      FIELD_u(L)=eigval
     end if
     l=l+1

    end do
    FIELD_x(k)=duplicate_x/number
    FIELD_u(k)=duplicate_u/number
    k=k+1
  end do

 ! Eliminate points as eigenvalues
  number=0
  do i=1,mpoint_random
   if(abs(FIELD_u(i)).lt.abs(eigval*0.9)) then
    number=number+1
    FIELD_x(number)=FIELD_x(i)
    FIELD_y(number)=FIELD_y(i)
    FIELD_u(number)=FIELD_u(i)
   end if
  end do
  mpoint_random=number
end subroutine duplicate_2d_eigval_sub

!*****************************************
!Sorting real arrays
!************************************************
subroutine shell_3real(A,A1,A2,n,mm,nn)
 dimension A(1:n),A1(1:n),A2(1:n)

 b=nn-mm+1
 b=alog(b)/alog(2.0)
 m=b
 l=1
 do i=1,m
  l=l*2
 end do
 do i=1,m
  k=l-1
  l=l/2
  do j=k+mm,nn
   t=A(j)
   t1=A1(j)
   t2=A2(j)
   is=j-k
   do while((is.gt.mm-1).and.(A(is).gt.t))
    A(is+k)=A(is)
    A1(is+k)=A1(is)
    A2(is+k)=A2(is)
    is=is-k
   end do 
  A(is+k)=t
  A1(is+k)=t1
  A2(is+k)=t2
  end do
 end do

end subroutine shell_3real

!*********************************************************
!Determine the initial value and constraint point of 2-D grid data node
!*************************************************************
subroutine method_inital_2d_random_sub(eigval,mpoint_random,FIELD_x,FIELD_y,FIELD_u,xmin,xmax,mpoint,&
                                       ymin,ymax,line,FIELD,DELT_x,DELT_y,DELT_u)

 
 real eigval,xmin,xmax,ymin,ymax
 integer mpoint_random,mpoint,line
 real FIELD(1:mpoint,1:line)
 real FIELD_x(1:mpoint_random),FIELD_y(1:mpoint_random),FIELD_u(1:mpoint_random)
 real DELT_x(1:mpoint,1:line),DELT_y(1:mpoint,1:line),DELT_u(1:mpoint,1:line)

 dx=(xmax-xmin)/(mpoint-1)
 dy=(ymax-ymin)/(line-1)
 eps_x=epsilon(eps_x)*dx
 eps_y=epsilon(eps_y)*dy

 !Assign the coordinate values of grid constraint points as eigenvalues
 DELT_x=eigval
 DELT_y=eigval
 DELT_u=eigval
 DO k=1,mpoint_random
  xk=(FIELD_x(k)-xmin)/dx!Deal with points in X direction
  i1=int(xk)+1
  i2=i1+1
  xi1=FIELD_x(k)-(xmin+(i1-1)*dx)
  xi2=FIELD_x(k)-(xmin+(i2-1)*dx)
  if(abs(xi1).lt.eps_x) then ! Handle special cases (at nodes)
   i2=i1
   xi2=xi1
  end if
  yk=(FIELD_y(k)-ymin)/dy!Deal with points in Y direction
  j1=int(yk)+1
  j2=j1+1
  yj1=FIELD_y(k)-(ymin+(j1-1)*dy)
  yj2=FIELD_y(k)-(ymin+(j2-1)*dy)
  if(abs(yj1).lt.eps_y) then !Handle special cases (at nodes)
   j2=j1
   yj2=yj1
  end if

  IF((i1.ge.1).and.(i1.le.mpoint).and.(j1.ge.1).and.(j1.le.line)) then
   if(abs(xi1*xi1+yj1*yj1).lt.abs(DELT_x(i1,j1)**2+DELT_y(i1,j1)**2)) then
    DELT_x(i1,j1)=xi1
    DELT_y(i1,j1)=yj1
    DELT_u(i1,j1)=FIELD_u(k)
   end if
  end IF
  IF((i1.ge.1).and.(i1.le.mpoint).and.(j2.ge.1).and.(j2.le.line)) then
   if(abs(xi1*xi1+yj2*yj2).lt.abs(DELT_x(i1,j2)**2+DELT_y(i1,j2)**2)) then
    DELT_x(i1,j2)=xi1
    DELT_y(i1,j2)=yj2
    DELT_u(i1,j2)=FIELD_u(k)
   end if
  end IF
  IF((i2.ge.1).and.(i2.le.mpoint).and.(j1.ge.1).and.(j1.le.line)) then
   if(abs(xi2*xi2+yj1*yj1).lt.abs(DELT_x(i2,j1)**2+DELT_y(i2,j1)**2)) then
    DELT_x(i2,j1)=xi2
    DELT_y(i2,j1)=yj1
    DELT_u(i2,j1)=FIELD_u(k)
   end if
  end IF
  IF((i2.ge.1).and.(i2.le.mpoint).and.(j2.ge.1).and.(j2.le.line)) then
   if(abs(xi2*xi2+yj2*yj2).lt.abs(DELT_x(i2,j2)**2+DELT_y(i2,j2)**2)) then
    DELT_x(i2,j2)=xi2
    DELT_y(i2,j2)=yj2
    DELT_u(i2,j2)=FIELD_u(k)
   end if
  end IF
 end DO

 FIELD=DELT_u
   
end subroutine method_inital_2d_random_sub

!*********************************************************
!Determine the 2-D grid data endpoint values and blank points
!*************************************************************
subroutine METHOD_edge_2D_eigval_sub(method_edge,method_restrict,eigval,mpoint,line,FIELD,number_eigval,M_eigval,N_eigval)
 real eigval
 integer method_edge,method_restrict,mpoint,line,number_eigval
 real FIELD(1:mpoint,1:line)
 integer M_eigval(1:number_eigval),N_eigval(1:number_eigval)
 logical flage_number
 integer radius

 IF(method_restrict.eq.2) then !The end point is not used as a constraint
  !Find the blank points and save them in the array
  number_eigval=0
  do j=1,line
   do i=1,mpoint
    if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
     number_eigval=number_eigval+1
     M_eigval(number_eigval)=i
     N_eigval(number_eigval)=j
    end if
   end do
  end do
 end IF
  !When the blank point is the endpoint value
  IF((number_eigval.gt.0).and.(number_eigval.lt.mpoint*line)) then
   IF(method_edge.eq.1) then!The endpoint value is 0
    field_mean=0.0
   else IF(method_edge.eq.3) then!The endpoint value is the average of all the numbers
    number=0
    field_mean=0.0
    do j=1,line
     do i=1,mpoint
      if(abs(FIELD(i,j)).lt.abs(eigval*0.9)) then
       number=number+1
       field_mean=field_mean+FIELD(i,j)
      end if
     end do
    end do
    field_mean=field_mean/number
    else IF(method_edge.eq.4) then!The endpoint value is the value of the nearest actual data
     i=1 !The first point
     DO while((i.eq.1).or.(i.eq.mpoint))
      Do j=1,line
       if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
        number=0
        field_mean=0.0
        radius=1!Initial value of circulation radius
        DO while(number.eq.0)
         do k=max(1+1,i-radius),min(mpoint-1,i+radius)
          do m=max(1+1,j-radius),min(line-1,j+radius)
           if(abs(FIELD(k,m)).lt.abs(eigval*0.9)) then
            field_mean=field_mean+FIELD(k,m)
            number=number+1
           end if
          end do
         end do
         radius=radius+1
        end DO
        FIELD(i,j)=field_mean/number
       end if
      end DO
      i=(i-1)+mpoint!The last point
     end DO
    
    j=1
    DO while((j.eq.1).or.(j.eq.line)) 
     DO i=1,mpoint
      if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
       number=0
       field_mean=0.0
       radius=1
       DO while(number.eq.0)
        do k=max(1+1,i-radius),min(mpoint-1,i+radius)
          do m=max(1+1,j-radius),min(line-1,j+radius)
           if(abs(FIELD(k,m)).lt.abs(eigval*0.9)) then
            field_mean=field_mean+FIELD(k,m)
            number=number+1
           end if
          end do
         end do
         radius=radius+1
       end DO
       FIELD(i,j)=field_mean/number
      end if
     end DO
    j=(j-1)+line
    end DO
    else!The endpoint value is taken as the average value of the boundary data
    m2=mpoint
    m3=1
    n2=line
    n3=1
    DO j=1,line
     DO i=1,mpoint
      if(abs(FIELD(i,j)).lt.abs(eigval*0.9)) then
       m2=min(m2,i)
       m3=max(m3,i)
       n2=min(n2,j)
       n3=max(n3,j)
      end if
     end DO
    end DO
    
    number=0
    field_mean=0.0
    DO j=n2,n3
     if(abs(FIELD(m2,j)).lt.abs(eigval*0.9)) then
      number=number+1
      field_mean=field_mean+FIELD(m2,j)
     end if
     if(abs(FIELD(m3,j)).lt.abs(eigval*0.9)) then
      number=number+1
      field_mean=field_mean+FIELD(m3,j)
     end if
    end DO

    DO i=m2+1,m3-1
     if(abs(FIELD(i,n2)).lt.abs(eigval*0.9)) then
      number=number+1
      field_mean=field_mean+FIELD(i,n2)
     end if
     if(abs(FIELD(i,n3)).lt.abs(eigval*0.9)) then
      number=number+1
      field_mean=field_mean+FIELD(i,n3)
     end if
    end DO
    field_mean=field_mean/number
   end IF

   IF (method_edge.eq.4) then!Endpoint assignment
    DO j=1,line
     if(abs(FIELD(1,j)).gt.abs(eigval*0.9)) FIELD(1,j)=field_mean
     if(abs(FIELD(mpoint,j)).gt.abs(eigval*0.9)) FIELD(mpoint,j)=field_mean
    end DO

    DO i=1,mpoint
     if(abs(FIELD(i,1)).gt.abs(eigval*0.9)) FIELD(i,1)=field_mean
     if(abs(FIELD(i,line)).gt.abs(eigval*0.9)) FIELD(i,line)=field_mean
    end DO
   end IF
  end IF
  
  IF(method_restrict.NE.2)then!End point as constraint point
  !Find the blank point and save it in the array
  number_eigval=0
  DO j=1,line
   DO i=1,mpoint
    if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
     number_eigval=number_eigval+1
     M_eigval(number_eigval)=i
     N_eigval(number_eigval)=j
    end if
   end DO
  end DO
  end IF
end subroutine METHOD_edge_2D_eigval_sub

!******************************************************************
!Initial value of blank point in 2-D grid data
!Submodule：
!method_inital_2d_1Dbi:Determinate initial value of 1-D bidirectional interpolation
!method_inital_2d_distance：Determinate initial value of 2-D bidirectional interpolation
!*********************************************************************
subroutine method_inital_2d_sub(method_initial,eigval,mpoint,line,FIELD)
 real eigval
 integer method_initial,mpoint,line
 real FIELD(1:mpoint,1:line)

 !Find the number of blank points:number_eigval
 number_eigval=0
 do j=1,line
  do i=1,mpoint
   if(abs(FIELD(i,j)).gt.abs(eigval*0.9))then
    number_eigval=number_eigval+1
   end if
  end do
 end do

 !Determine the initial value
 IF((number_eigval.gt.0).and.(number_eigval.lt.mpoint*line)) then
  IF(method_initial.eq.1) then !The initial value is 0
   do j=1,line
    do i=1,mpoint
     if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
      FIELD(i,j)=0
     end if
    end do
   end do
  ELSE IF((method_initial.eq.2).or.(method_initial.eq.3)) then !Bidirectional interpolation
    call method_inital_2d_1Dbi(method_initial,eigval,mpoint,line,FIELD)
  ELSE
    call method_inital_2d_distance(eigval,mpoint,line,FIELD)

  end IF
 end IF

end subroutine method_inital_2d_sub

!******************************************************************
!Use 1-D bidirectional interpolation to complete 2-D interpolation
!submodule：
!method_inital_1d_sub:determinate initial value of 1-D bidirectional interpolation
!*********************************************************************
subroutine method_inital_2d_1Dbi(method_initial,eigval,mpoint,line,FIELD)
 real eigval
 integer method_initial,mpoint,line
 real FIELD(1:mpoint,1:line)
 real,allocatable::FIELD_y(:,:),FIELD_x(:,:)
 real,allocatable::FIELD_tmp(:)

 allocate(FIELD_x(1:mpoint,1:line),FIELD_y(1:mpoint,1:line),STAT=ierr)
 allocate(FIELD_tmp(1:max(line,mpoint)),STAT=ierr)

 !Find the number of blank points:number_eigval
 number_eigval=0
 do j=1,line
  do i=1,mpoint
   if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
    number_eigval=number_eigval+1
   end if 
  end do
 end do
 
 m3=1
 n2=line
 n3=1
 do j=1+1,line-1
  do i=1+1,mpoint-1
   if(abs(FIELD(i,j)).lt.abs(eigval*0.9)) then
    m2=min(m2,i)
    m3=max(m3,i)
    n2=min(n2,j)
    n3=max(n3,j)
   end if
  end do
 end do
 
 if(m2.eq.2) m2=1
 if(m3.eq.mpoint-1) m3=mpoint
 if(n2.eq.2) n2=1
 if(n3.eq.line-1) n3=line

 FIELD_x=FIELD
 FIELD_y=FIELD

 !Determine the initial value
  nstep=max(1,line/10)
  mstep=max(1,mpoint/10)
  DO while((number_eigval.gt.0).and.(number_eigval.lt.mpoint*line))
   ! the initial value is determined along the X (m) direction firstly, 
   ! and then the initial value is determined along the Y (n) direction
    DO jj=n2,n3+(nstep-1),nstep
     j=min(jj,n3)
     do i=1,mpoint
      FIELD_tmp(i)=FIELD_x(i,j)
     end do
     call method_inital_1d_sub(method_initial,eigval,mpoint,FIELD_tmp)
     do i=1,mpoint
      FIELD_x(i,j)=FIELD_tmp(i)
     end do
    end DO
    DO i=1,mpoint,mstep
     do j=1,line
      FIELD_tmp(j)=FIELD_x(i,j)
     end do
     call method_inital_1d_sub(method_initial,eigval,line,FIELD_tmp)
     do j=1,line
      FIELD_x(i,j)=FIELD_tmp(j)
     end do
    end DO

    !the initial value is determined along the Y (n) direction firstly, 
   ! and then the initial value is determined along the X (m) direction
    DO ii=m2,m3+(mstep-1),mstep
     i=min(ii,m3)
     do j=1,line
      FIELD_tmp(j)=FIELD_y(i,j)
     end do
     call method_inital_1d_sub(method_initial,eigval,line,FIELD_tmp)
     do j=1,line
      FIELD_y(i,j)=FIELD_tmp(j)
     end do
    end DO
    DO j=1,line,nstep
     do i=1,mpoint
      FIELD_tmp(i)=FIELD_y(i,j)
     end do
     call method_inital_1d_sub(method_initial,eigval,mpoint,FIELD_tmp)
     do i=1,mpoint
      FIELD_y(i,j)=FIELD_tmp(i)
     end do
    end DO
	
    !merge the interpolation results of two directions
    DO j=1,line
     do i=1,mpoint
      if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
       if((FIELD_y(i,j).lt.eigval*0.9).and.(FIELD_x(i,j).lt.eigval*0.9)) then
        FIELD(i,j)=(FIELD_y(i,j)+FIELD_x(i,j))/2.0
       end if
      end if
     end do
    end DO
	
   !Find the number of eigenvalues
   number_eigval=0
   DO j=1,line
    DO i=1,mpoint
     if(abs(FIELD(i,j)).gt.abs(eigval*0.9)) then
      number_eigval=number_eigval+1
     end if
    end DO
   end DO
   mstep=1
   nstep=1

  end DO 

  deallocate(FIELD_x,FIELD_y,FIELD_tmp,STAT=ierr)

end subroutine method_inital_2d_1Dbi

!******************************************************************
!Use 1-D grid data to determine the initial value of blank value
!*********************************************************************
subroutine method_inital_1d_sub(method_initial,eigval,mpoint,FIELD)
 real eigval
 integer method_initial,mpoint
 real FIELD(1:mpoint)
 
 pi=3.141592654
 i1=1
 i2=1

 DO while(i2.le.mpoint)

  !determinate the blank data'S starting point (i1) and ending point(i2)
   i=i2
   do while((i.le.mpoint).and.(abs(FIELD(i)).lt.abs(eigval*0.9)))
    i=i+1
   end do
   !i1=max(1,i-1)
   i1=i-1
   i=i1+1

   do while((i.le.mpoint).and.(abs(FIELD(i)).gt.abs(eigval*0.9)))
    i=i+1
   end do
   !i2=min(mpoint,i)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
   i2=i

   
   IF(((i2-i1).gt.0).and.((i2-i1).lt.(mpoint-1))) then
    if(method_initial.eq.1) then! The blank point is 0
     do i=i1+1,i2-1
      FIELD(i)=0.0
     end do
    else if(method_initial.eq.3) then!The blank point is determinated by distance weighting method
     do i=i1+1,i2-1
      if((i1.ge.1).and.(i2.le.mpoint)) then
        distance_i1=1.0/(i1-i)**2
        distance_i2=1.0/(i2-i)**2
        distance=distance_i1+distance_i2
        sum_dis=FIELD(i1)*distance_i1+FIELD(i2)*distance_i2
        FIELD(i)=sum_dis/distance
      else if((i1.ge.1).and.(i2.gt.mpoint).and.(FIELD(i1-1).lt.(eigval*0.9))) then
      !else if((i1.ge.1).and.(i2.gt.mpoint)) then
        distance_i1=1.0/(i1-i)**2
        distance_i2=1.0/(i1-1-i)**2
        distance=distance_i1+distance_i2
        sum_dis=FIELD(i1)*distance_i1+FIELD(i1-1)*distance_i2
        FIELD(i)=sum_dis/distance
      else if((i1.lt.1).and.(i2.le.mpoint).and.(FIELD(i2+1).lt.(eigval*0.9))) then
      !else if((i1.lt.1).and.(i2.le.mpoint)) then
        distance_i1=1.0/(i2+1-i)**2
        distance_i2=1.0/(i2-i)**2
        distance=distance_i1+distance_i2
        sum_dis=FIELD(i2+1)*distance_i1+FIELD(i2)*distance_i2
        FIELD(i)=sum_dis/distance
      else
        FIELD(i)=0.0
      end if
     end do
    else!The blank point is determinated by cosine decreasing method
     x0=pi/(i2-i1)/2.0

     do i=i1+1,i2-1
       x=x0*(i2-i)
       if((i1.ge.1).and.(i2.le.mpoint)) then
         FIELD(i)=(FIELD(i2)-FIELD(i1))*cos(x)+FIELD(i1)
       else if((i1.ge.1).and.(i2.gt.mpoint).and.(FIELD(i1-1).lt.(eigval*0.9))) then
         FIELD(i)=(FIELD(i1-1)-FIELD(i1))*cos(x)+FIELD(i1)
       else if((i1.lt.1).and.(i2.le.mpoint).and.(FIELD(i2+1).lt.(eigval*0.9))) then
         FIELD(i)=(FIELD(i2)-FIELD(i2+1))*cos(x)+FIELD(i2+1)
       else
         FIELD(i)=0.0
       end if
     end do 
    end if
   ELSE
     do i=i1+1,i2-1
      FIELD(i)=0.0
     end do
   end IF
   i2=i2+1
 end DO



end subroutine method_inital_1d_sub

!******************************************************************
!Distance weighted interpolation of 2-D grid data
!*********************************************************************
subroutine method_inital_2d_distance(eigval,mpoint,line,FIELD)
 real eigval
 integer mpoint,line
 real FIELD(1:mpoint,1:line)
 real,allocatable::FIELD_dis(:,:)
 real,allocatable::FIELD_real(:)
 integer,allocatable::M_real(:),N_real(:)
 
 !Find the number of eigenvalues:number_eigval
 number_eigval=0
 DO j=1,line
  DO i=1,mpoint
   if(abs(FIELD(i,j)).gt.abs(eigval*0.9))then
    number_eigval=number_eigval+1
   end if
  end DO
 end DO

 IF((number_eigval.gt.0).and.(number_eigval.lt.mpoint*line)) then
  allocate(FIELD_dis(1:mpoint,l:line),STAT=ierr)
  FIELD_dis=FIELD
  do j=1+1,line-1 
   do i=1+1,mpoint-1
    if(abs(FIELD(i,j)).lt.abs(eigval*0.9)) then
     number_real=0
     if(FIELD(i,j-1).lt.eigval*0.9) number_real=number_real+1
     if(FIELD(i,j+1).lt.eigval*0.9) number_real=number_real+1
     if(FIELD(i-1,j).lt.eigval*0.9) number_real=number_real+1
     if(FIELD(i+1,j).lt.eigval*0.9) number_real=number_real+1
     if(number_real.eq.4) then 
       FIELD_dis(i,j)=eigval
     end if
    end if
   end do
  end do

  !Find the number of non blank points
  number_real=0
  do j=1,line
   do i=1,mpoint
    if(abs(FIELD_dis(i,j)).lt.abs(eigval*0.9)) then
     number_real=number_real+1
    end if
   end do 
  end do 

  !Put the coordinates and values of non blank points into the array:field_real，M_real,N_real
  allocate(FIELD_real(1:number_real),M_real(1:number_real),N_real(1:number_real),STAT=ierr)
  number=0
  do j=1,line
   do i=1,mpoint
    if(FIELD_dis(i,j).lt.eigval*0.9) then
     number=number+1
     FIELD_real(number)=FIELD_dis(i,j)
     M_real(number)=i
     N_real(number)=j
    end if
   end do
  end do

 !Some points are weighted and others are determined by cosine interpolation
  
  DO j=1,line,max(1,(line/10))
   DO i=1,mpoint,max(1,(mpoint/10))
    if(FIELD(i,j).gt.eigval*0.9) then
     distance=0.0
     sum_dis=0.0
     do k=1,number_real
      distanceij=1.0/((M_real(k)-i)**2+(N_real(k)-j)**2)
      distance=distance+distanceij
      sum_dis=sum_dis+FIELD_real(k)*distanceij
     end do
     FIELD(i,j)=sum_dis/distance
    end if
   end DO
  end DO

  deallocate(FIELD_dis,STAT=ierr)
  deallocate(FIELD_real,M_real,N_real,STAT=ierr)
  call method_inital_2d_1Dbi(2,eigval,mpoint,line,FIELD)
 end IF
end subroutine method_inital_2d_distance

!******************************************************************
! Iterate constrained minimum curvature iteration for 2-D grid data
!Submodule：mincurv_2D_boundary:Minimum curvature boundary processing for 2-D grid data
!*********************************************************************
subroutine MinCurV_2D_random_sub(eps_abs,iteration_max,eigval,DELT_x,DELT_y,DELT_u,&
                              xmin,xmax,mpoint,ymin,ymax,line,FIELD,flage,&
                               number_eigval,M_eigval,N_eigval)
  real eps_abs,eigval,xmin,xmax,ymin,ymax
  integer iteration_max,mpoint,line,number_eigval
  logical flage
  real FIELD(1:mpoint,1:line)
  real DELT_x(1:mpoint,1:line),DELT_y(1:mpoint,1:line),DELT_u(1:mpoint,1:line)
  integer  M_eigval(1:number_eigval),N_eigval(1:number_eigval)

  integer,allocatable::M_flage(:,:),S(:,:),T(:,:)
  real,allocatable::U(:,:)
  real,allocatable::A0(:,:),A1(:,:),A2(:,:),A3(:,:),A4(:,:),A5(:,:)

  allocate(U(1-2:mpoint+2,1-2:line+2),STAT=ierr)
  allocate(M_flage(1:mpoint,1:line),STAT=ierr)
  allocate(A0(1:mpoint,1:line),A1(1:mpoint,1:line),A2(1:mpoint,1:line),&
              A3(1:mpoint,1:line),A4(1:mpoint,1:line),A5(1:mpoint,1:line),stat=ierr)
  allocate(S(1:mpoint,1:line),T(1:mpoint,1:line),STAT=ierr)

 !Assign a flag value to the node, and the default value is 0, that is, no iteration is performed
  M_flage=0
  do k=1,number_eigval
    M_flage(M_eigval(k),N_eigval(k))=-1 !Mark of unconstrained point
  end do
  do j=1,line
    do i=1,mpoint
        if((abs(DELT_x(i,j)).lt.eigval*0.9).and.(abs(DELT_y(i,j)).lt.eigval*0.9))then
            M_flage(i,j)=1!Mark of constrained point
        end if

    end do
  end do 

  dx=(xmax-xmin)/(mpoint-1)
  dy=(ymax-ymin)/(line-1)
  alph=dx/dy
  beta=dy/dx
  alph0=-1.0/(2*(3+4*alph*alph+3*(alph**4)))
  alph1=(alph**4)
  alph2=2*(alph*alph)
  alph3=-4*(1+alph*alph)
  alph4=-4*(1+alph*alph)*alph*alph
  alph5=-2*(1+alph*alph)
  alph6=-2*(1+alph*alph)*alph*alph

  do j=1,line
   do i=1,mpoint
    IF((abs(DELT_x(i,j)).lt.eigval*0.9).and.(abs(DELT_y(i,j)).lt.eigval*0.9)) then
     if(DELT_x(i,j).gt.0.0)then
      S(i,j)=1
      else
      S(i,j)=-1
     end if
     if(DELT_y(i,j).gt.0.0)then
      T(i,j)=1
      else
      T(i,j)=-1
     end if
     DELT_x(i,j)=abs(DELT_x(i,j)/dx)
     DELT_y(i,j)=abs(DELT_y(i,j)/dy)
     beta1=2*(1+alph*alph)**2
     beta2=(1+alph**4)
     denom=beta2*(DELT_x(i,j)+DELT_y(i,j))*(1+DELT_x(i,j)+DELT_y(i,j))&
           +beta1*(DELT_x(i,j)*DELT_y(i,j)+DELT_x(i,j)+DELT_y(i,j)+1)
     denom=denom*2.0
     A0(i,j)=-(DELT_x(i,j)+DELT_y(i,j))*(1+DELT_x(i,j)+DELT_y(i,j))/denom
     A1(i,j)=beta1*DELT_x(i,j)*(1+DELT_x(i,j)+2*DELT_y(i,j))/denom
     A2(i,j)=-2*beta1*DELT_x(i,j)*(DELT_x(i,j)+DELT_y(i,j))/denom
     A3(i,j)=-2*beta1*DELT_y(i,j)*(DELT_x(i,j)+DELT_y(i,j))/denom
     A4(i,j)=-beta1*DELT_x(i,j)*(1+DELT_x(i,j))/denom
     A5(i,j)=2*beta1/denom
    end IF
   end do
  end do

  !Assign the input data to a newly defined temporary array:U(1:MPOINT)
  do j=1,line
   do i=1,mpoint
    U(i,j)=FIELD(i,j)
   end do
  end do
  
  !Give the relevant initial values
  eps_error=2*abs(eps_abs)!Initial value of residual maximum absolute value 
  iteration=0!Initial value of iteration number

  DO while((eps_error.gt.eps_abs).and.(iteration.le.iteration_max))

   
   !Dealing with boundaries
   call MinCurV_2D_boundary(U,mpoint,line,alph)
   !Start iteration
   eps_error=0.0
   DO j=1,line
    DO i=1,mpoint
     IF(M_flage(i,j).eq.-1)then !The boundary points are unconstrained points
     temp=(U(i+2,j)+U(i-2,j))+alph1*(U(i,j+2)+U(i,j-2))&
                            +alph2*(U(i+1,j+1)+U(i-1,j+1)+U(i+1,j-1)+U(i-1,j-1))&
                            +alph3*(U(i+1,j)+U(i-1,j))&
                            +alph4*(U(i,j+1)+U(i,j-1))
      temp=temp*alph0
      eps_error=max(abs(temp-U(i,j)),eps_error)
      U(i,j)=temp
     ELSE IF(M_flage(i,j).eq.1)then !The boundary points are constrained points
      temp=(U(i+2,j)+U(i-2,j))+alph1*(U(i,j+2)+U(i,j-2))&
                            +alph2*(U(i+1,j+1)+U(i-1,j+1)+U(i+1,j-1)+U(i-1,j-1))&
                            +alph5*(U(i+1,j)+U(i-1,j))&
                            +alph6*(U(i,j+1)+U(i,j-1))
      temp=temp*A0(i,j)
      temp=temp+A0(i,j)*alph5*(-U(i-S(i,j),j+T(i,j))+2*U(i-S(i,j),j)+alph2*U(i,j-T(i,j))+U(i+S(i,j),j-T(i,j)))
      temp=temp+A1(i,j)*U(i-S(i,j),j+T(i,j))&
               +A2(i,j)*U(i-S(i,j),j)&
               +A3(i,j)*U(i,j-T(i,j))&
               +A4(i,j)*U(i+S(i,j),j-T(i,j))&
               +A5(i,j)*DELT_u(i,j)
      eps_error=max(abs(temp-U(i,j)),eps_error)
      U(i,j)=temp
     ELSE
     end IF
    end DO
   end DO
   iteration=iteration+1
   
  end DO

  if(iteration.le.iteration_max) then
   flage=.true.
  else
   flage=.false.
  end if
  
  do j=1,line
   do i=1,mpoint
    FIELD(i,j)=U(i,j)
   end do
  end do

  deallocate(M_flage,S,T,U,A0,A1,A2,A3,A4,A5)
end subroutine MinCurV_2D_random_sub

!******************************************************************
!2-D unconstrained minimum curvature iterative boundary treatment
!*********************************************************************
subroutine MinCurV_2D_boundary(U,mpoint,line,alph)
 real U(1-2:mpoint+2,1-2:line+2)

 !The first boundary condition is processed at the i-1:mpoint point
  do j=1,line
    U(1-1,j)=2*U(1,j)-U(1+1,j)!The left point
    U(mpoint+1,j)=2*U(mpoint,j)-U(mpoint-1,j)!The rightt point
  end do

 !The first boundary condition is processed on the j-1:line line
 do i=1,mpoint
  U(i,1-1)=2*U(i,1)-U(i,1+1)
  U(i,line+1)=2*U(i,line)-U(i,line-1)
 end do

 !Processing corner point
 U(1-1,1-1)=U(1+1,1-1)+U(1-1,1+1)-U(1+1,1+1) 
 U(mpoint+1,1-1)=U(mpoint+1,1+1)+U(mpoint-1,1-1)-U(mpoint-1,1+1) 
 U(1-1,line+1)=U(1+1,line+1)+U(1-1,line-1)-U(1+1,line-1) 
 U(mpoint+1,line+1)=U(mpoint+1,line-1)+U(mpoint-1,line+1)-U(mpoint-1,line-1) 

 !The value of the second boundary condition at the second point in the i direction
 i=1!left
 alph1=alph*alph
 alph2=-2*(1+alph1)
 do j=1,line
    Pij=alph1*(U(i+1,j+1)-U(i-1,j+1)+U(i+1,j-1)-U(i-1,j-1))&
        +alph2*(U(i+1,j)-U(i-1,j))
    U(i-2,j)=U(i+2,j)+Pij
 end do
 i=mpoint!right
 do j=1,line
    Pij=alph1*(U(i+1,j+1)-U(i-1,j+1)+U(i+1,j-1)-U(i-1,j-1))&
        +alph2*(U(i+1,j)-U(i-1,j))
    U(i+2,j)=U(i-2,j)-Pij
 end do

 !The value of the second line in the j direction of the second boundary condition
 j=1!bottom
 beta=1.0/alph
 beta1=beta*beta
 beta2=-2*(1+beta1)
 do i=1,mpoint 
    Qij=beta1*(U(i+1,j+1)-U(i+1,j-1)+U(i-1,j+1)-U(i-1,j-1))&
        +beta2*(U(i,j+1)-U(i,j-1))
    U(i,j-2)=U(i,j+2)+Qij
 end do
 j=line!top
 do i=1,mpoint 
    Qij=beta1*(U(i+1,j+1)-U(i+1,j-1)+U(i-1,j+1)-U(i-1,j-1))&
        +beta2*(U(i,j+1)-U(i,j-1))
    U(i,j+2)=U(i,j-2)-Qij
 end do

end subroutine MinCurV_2D_boundary
!******************************************************************
!Grid data output
!*********************************************************************
subroutine OUTPUT_GRD(A,filename,m,n,xmin,xmax,ymin,ymax,m1,m2,n1,n2,eigval)
 real A(m,n)
 character*(*)filename
 integer nunit_out

 zmin=A(1,1)
 zmax=A(1,1)

 do j=n1,n2
    do i=m1,m2
        if(A(i,j).lt.eigval/2.0) then
            zmin=min(zmin,A(i,j))
            zmax=max(zmax,A(i,j))
        end if
    end do
 end do
 call search_unit(10,nunit_out)
 open (nunit_out,file=filename,status='unknown',form='formatted')
 close (nunit_out,status='delete')
 open (nunit_out,file=filename,status='unknown',form='formatted')
 write(nunit_out,'(a)') 'DSAA'
 write(nunit_out,*) m2-m1+1,n2-n1+1
 write(nunit_out,*) xmin,xmax
 write(nunit_out,*) ymin,ymax
 write(nunit_out,*) zmin,zmax
 do j=n1,n2
    write(nunit_out,*)(A(i,j),i=m1,m2)
 end do
 close (nunit_out)
end subroutine OUTPUT_GRD


!end in 12.19 2019 by Yin Yihong
