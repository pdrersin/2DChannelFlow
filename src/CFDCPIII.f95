!Final course project MAE 6263 - Dr. Jayaraman
!2-D Pressure Driven Channel Flow
!Explicit Solver
!Romit Maulik - PhD Student - Computational Fluid Dynamics Laboratory
!Oklahoma State University - Stillwater

program NS2D
implicit none


integer :: nx,ny,fostep,ns,i,j,outpar,k,psit
real*8 :: xl,yl,regl,ft,gm,pedx,nu,ucl,dy,dx,dt,delp,rho
real*8 :: t1,t2,cflu,cflv
real*8,allocatable :: u(:,:),v(:,:),p(:,:),hx(:,:),hy(:,:)

common/Reynolds/ regl
common/PoissonIter/psit

open(15,file='CFDCP3Input.txt')
read(15,*)pedx		!Cell Peclet number
read(15,*)regl		!Global Reynolds number
read(15,*)nu		!Kinematic Viscosity of fluid
read(15,*)rho		!Density of fluid
read(15,*)gm		!Viscous courant
read(15,*)xl		!Total Length in x direction
read(15,*)yl		!Total Length in y direction
read(15,*)ft		!Final Time
read(15,*)fostep	!fostep;File output every this percent of total timesteps(Choose multiples of ten)
read(15,*)psit		!Number of iterations in Poisson Solver for Pressure
close(15)

!Calculating U at centerline
ucl = nu*regl/yl

!Calculating grid discretizations
dx = pedx*nu/ucl
dy = dx

!Calculating number of timesteps
ns = nint(ft/dt)

!Calculating number of grid discretizations
nx = nint(xl/dx)
ny = nint(yl/dy)

!Calculating number of timesteps
ns = nint(ft/dt)

!Calculating number of grid discretizations
nx = nint(xl/dx)
ny = nint(yl/dy)


!Initial Condition & Boundary Condition Setup
allocate(u(0:nx,0:ny))
allocate(v(0:nx,0:ny))
allocate(p(0:nx,0:ny))
allocate(hx(0:nx,0:ny))
allocate(hy(0:nx,0:ny))


do i = 0,nx
  do j = 0,ny
    u(i,j) = 0
    v(i,j) = 0
    p(i,j) = 0
  end do
end do

!Calculating pressure drop using Darcy-Weisbach relation
!Note that P2 (right buondary) is kept at reference zero
do j = 0,ny
  p(nx,j) = 0
end do

delp = xl*64/(regl)*(ucl**2)/(2d0*yl)*rho

do j = 0,ny
  p(0,j) = delp
end do

!IC File Ouput
open(20,file="InitialField.plt")
write(20,*)'Title="IC Data set"'
write(20,*)'variables ="x","y","u","v"'
close(20)

open(20,file="InitialField.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.3)',advance="no")dfloat(i)/dfloat(ny),dfloat(j)/dfloat(ny),u(i,j),v(i,j)
      write(20,*) ''
    end do
  end do
close(20)


!Output file setup
open(20,file="ContourPlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="x","y","u","v"'
close(20)

!Output file setup at y=0.5
open(20,file="LinePlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="y","u"'
close(20)


call cpu_time(t1)
!Time integration - 
do k = 1,ns

!Stability Check
cflu = maxval(u)
cflv = maxval(v) 
if ((cflu*dt/(dx)+cflv*dt/(dy))>1d0) then
  print*,'Unstable - Reduce Timestep'
  print*,dt
!  call exit(10)
end if

call updateH(u,v,hx,hy,nx,ny,dx,dy)

call updateP(hx,hy,nx,ny,dx,dy)

call updatefield(u,v,p,hx,hy,nx,ny,dx,dy,dt)

outpar = ns*fostep/100

  
!Output to transient .plt file for Tecplot  
if (mod(k,outpar)==0) then
open(20,file="ContourPlots.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',k
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.3,1600F14.3,1600F14.3,1600F14.3)',advance="no")dfloat(i)/dfloat(ny),dfloat(j)/dfloat(ny),u(i,j),v(i,j)
      write(20,*) ''
    end do
  end do
close(20)

open(20,file="LinePlots.plt",position="append")
write(20,"(a,i8,a)")'Zone I = ',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',k
    do j = 0,ny
      write (20, '(1600F14.3,1600F14.3)',advance="no")dfloat(j)/dfloat(ny),u(nx/2,j)
      write(20,*) ''
    end do
close(20)

end if


 
end do

call cpu_time(t2)

open(4,file='cpu.txt')
write(4,*)"cpu time (sec)=",(t2-t1)
close(4)


end


!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!Subroutine for Convection Diffusion Term Updation
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine updateH(u,v,hx,hy,nx,ny,dx,dy)
implicit none

integer :: nx,ny,i,j
real*8::dx,dy,regl
real*8,dimension(0:nx,0:ny):: u,v,hx,hy
real*8,dimension(-1:nx+1,0:ny)::utemp,vtemp
real*8,dimension(0:nx,0:ny)::udx,udy,uddx,uddy,vdx,vdy,vddx,vddy

common/Reynolds/ regl

!Dummy variables and periodic boundary conditions
do i = 0,nx
  do j = 0,ny
    utemp(i,j) = u(i,j)
    vtemp(i,j) = v(i,j)
  end do
end do

do j = 0,ny
  utemp(-1,j) = u(nx-1,j)
  vtemp(-1,j) = v(nx-1,j)
end do

!Calculating first order derivatives  
do i = 0,nx
  	do j = 1,ny-1
  		udx(i,j) = (utemp(i,j)-utemp(i-1,j))/dx
        vdx(i,j) = (vtemp(i,j)-vtemp(i-1,j))/dx
	end do
end do

do i = 0,nx
  	do j = 1,ny-1
  		udy(i,j) = (utemp(i,j)-utemp(i,j-1))/dy
        vdy(i,j) = (vtemp(i,j)-vtemp(i,j-1))/dy
	end do
end do

!Calculating second order derivatives
do i = 0,nx
  	do j = 1,ny-1
  		uddx(i,j) = (utemp(i+1,j)+utemp(i-1,j)-2.*utemp(i,j))/(2d0*dx)
        vddx(i,j) = (vtemp(i+1,j)+vtemp(i-1,j)-2.*vtemp(i,j))/(2d0*dx)
	end do
end do

do i = 0,nx
  	do j = 1,ny-1
  		uddy(i,j) = (utemp(i,j+1)+utemp(i,j-1)-2d0*utemp(i,j))/(2d0*dy)
        vddy(i,j) = (vtemp(i,j+1)+vtemp(i,j-1)-2d0*vtemp(i,j))/(2d0*dy)
	end do
end do

!Formulation of Hx,Hy arrays

do i = 0,nx
	do j = 1,ny-1
    	hx(i,j) = -utemp(i,j)*(udx(i,j))-vtemp(i,j)*(udy(i,j))+1/regl*(uddx(i,j)+uddy(i,j))
    	hy(i,j) = -utemp(i,j)*(vdx(i,j))-vtemp(i,j)*(vdy(i,j))+1/regl*(vddx(i,j)+vddy(i,j))
	end do
end do


return
end



!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!Subroutine for Pressure Updation
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine updateP(hx,hy,nx,ny,dx,dy)
implicit none

common/PoissonIter/psit

integer :: nx,ny,i,j,k,psit,q
real*8::dy,dx,ta,tb
real*8,dimension(0:nx,0:ny):: p,hx,hy,hxdx,hydy
real*8,dimension(-1:nx+1,0:ny)::hxtemp,hytemp,ptemp



!Calculating divergence of H terms

do i = 0,nx
  do j = 0,ny
    hxtemp(i,j) = hx(i,j)
    hytemp(i,j) = hy(i,j)
    ptemp(i,j)	= p(i,j)
  end do
end do
    
do i = 0,nx
  do j = 1,ny-1
    hxdx(i,j) = (hxtemp(i+1,j)-hxtemp(i-1,j))/(2d0*dx)
    hydy(i,j) = (hytemp(i,j+1)-hxtemp(i,j-1))/(2d0*dy)
  end do
end do    


do k = 1,psit

q = 0
    
do i = 1,nx-1
  do j = 1,ny-1
    if (j>1.and.j<ny-1) then
      ta = p(i+1,j)+ptemp(i-1,j)+p(i,j+1)+ptemp(i,j-1)
	  tb = hxdx(i,j) + hydy(i,j)
      ptemp(i,j) = ta/4-tb/4*(dx**2)
	else if (j==1) then
      ta = p(i+1,j)+ptemp(i-1,j)+p(i,j+1)
      tb = hxdx(i,j) + hydy(i,j)
      ptemp(i,j) = ta/3-tb/3*(dx**2)
    else if (j==ny-1) then
      ta = p(i+1,j)+ptemp(i-1,j)+p(i,j-1)
      tb = hxdx(i,j) + hydy(i,j)
      ptemp(i,j) = ta/3-tb/3*(dx**2)
	end if
  end do

	  ptemp(i,ny) = ptemp(i,ny-1)
      ptemp(i,0) = ptemp(i,1)
		  
end do

call l1normcheck(ptemp,p,nx,ny,q)

if (q.ne.0) then
    do j = 0,ny
  	do i = 0,nx
    	p(i,j) = ptemp(i,j)
    end do
	end do
    print*,k    
  exit
end if  

end do



return
end

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!Subroutine for L1 Norm Check
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine l1normcheck(utemp,u,nx,ny,q)
implicit none

common/contol/ tol

integer::q,nx,ny,i,j
real*8,dimension(-1:nx+1,0:ny)::utemp,u
real*8::sumu,tol

sumu = 0.

do j = 0,ny
	do i = 0,nx
  		sumu = sumu + (utemp(i,j)-u(i,j))
	end do
end do


if (sumu<tol) then
  q = 1
end if

return
end

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!Subroutine for velocity update
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine updatefield(u,v,p,hx,hy,nx,ny,dx,dy,dt)
implicit none

real*8::dx,dy,dt
integer::nx,ny,i,j
real*8,dimension(0:nx,0:ny)::u,v,p,hx,hy
real*8,dimension(0:nx,0:ny)::gradpx,gradpy
real*8,dimension(-1:nx+1,0:ny)::ptemp

!Defining dummys for pressure

do i = 0,nx
  do j = 0,ny
    ptemp(i,j) = p(i,j)
  end do
end do

do j = 0,ny
  ptemp(-1,j) = ptemp(nx-1,j)
  ptemp(nx+1,j) = ptemp(1,j)
end do


!Calculation of pressure gradient
do i = 1,nx-1
  do j = 1,ny-1
    gradpx(i,j) = (ptemp(i,j)-ptemp(i-1,j))/dx
    gradpy(i,j) = (ptemp(i,j)-ptemp(i,j-1))/dy    
  end do

  gradpx(i,ny)=0
  gradpy(i,0)=0
end do

!Updating the field
do i = 0,nx
  do j = 1,ny-1
	u(i,j) = u(i,j) + dt*(hx(i,j)-gradpx(i,j))
  	v(i,j) = v(i,j) + dt*(hy(i,j)-gradpy(i,j))
  end do
end do



return
end