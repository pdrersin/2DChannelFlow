!Final course project MAE 6263 - Dr. Jayaraman
!2-D Pressure Driven Channel Flow
!Explicit Solver
!Romit Maulik - PhD Student - Computational Fluid Dynamics Laboratory
!Oklahoma State University - Stillwater

program NS2D
implicit none


integer :: nx,ny,fostep,ns,i,j,outpar,k,psit
real*8 :: xl,yl,regl,ft,gm,nu,ucl,dy,dx,dt,delp,rho
real*8 :: t1,t2,cflu,cflv,tol
real*8,allocatable :: u(:,:),v(:,:),p(:,:),hx(:,:),hy(:,:)

common/Reynolds/ regl
common/PoissonIter/psit
common/contol/ tol
common/tint/k,outpar
common/density/rho

open(15,file='CFDCP3Input.txt')
read(15,*)nx		!Divisions in X direction
read(15,*)ny		!Divisions in y direction
read(15,*)ns		!Number of timesteps
read(15,*)regl		!Global Reynolds number
read(15,*)nu		!Kinematic Viscosity of fluid
read(15,*)rho		!Density of fluid
read(15,*)gm		!Viscous courant
read(15,*)xl		!Total Length in x direction
read(15,*)yl		!Total Length in y direction
read(15,*)ft		!Final Time
read(15,*)fostep	!fostep;File output every this percent of total timesteps(Choose multiples of ten)
read(15,*)tol		!Tolerance for Pressure Poisson Solver
close(15)

!Calculating U at centerline
ucl = 2.0d0

!Calculating grid discretizations
dx = xl/dfloat(nx)
dy = yl/dfloat(ny)

!Calculating timestep
dt = ft/dfloat(ns)

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
  p(0,j) = 2.0
end do

outpar = ns*fostep/100

!IC File Ouput
open(20,file="InitialField.plt")
write(20,*)'Title="IC Data set"'
write(20,*)'variables ="x","y","u","v","p"'
close(20)

open(20,file="InitialField.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.3)',advance="no")dfloat(i)/dfloat(ny),dfloat(j)/dfloat(ny),u(i,j),v(i,j),p(i,j)
      write(20,*) ''
    end do
  end do
close(20)

!Output file setup
open(20,file="ContourPlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="x","y","u","v","p"'
close(20)

open(20,file="GradPContourPlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="x","y","gradpx","gradpy"'
close(20)


!Output file setup at y=0.5
open(20,file="LinePlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="x","p"'
close(20)

call cpu_time(t1)
!Time integration - 
do k = 1,ns

!Stability Check
cflu = maxval(u)
cflv = maxval(v) 
if ((cflu*dt/(dx)+cflv*dt/(dy))>1d0) then
  print*,'Unstable - Reduce Timestep'
!  print*,dt
!  call exit(10)
end if

call updateH(u,v,hx,hy,nx,ny,dx,dy)

call updateP(p,hx,hy,nx,ny,dx,dy)

call updatefield(u,v,p,hx,hy,nx,ny,dx,dy,dt)



  
!Output to transient .plt file for Tecplot  
if (mod(k,outpar)==0) then
open(20,file="ContourPlots.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',k
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.3,1600F14.3,1600F14.3,1600F14.3)',advance="no")dfloat(i)/dfloat(ny),&
      &dfloat(j)/dfloat(ny),u(i,j),v(i,j),p(i,j)
      write(20,*) ''
    end do
  end do
close(20)

open(20,file="LinePlots.plt",position="append")
write(20,"(a,i8,a)")'Zone I = ',nx+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',k
    do i = 0,nx
      write (20, '(1600F14.3,1600F14.3)',advance="no")dfloat(i)/dfloat(ny),p(i,ny/2)
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

do j = 1,ny-1
  utemp(-1,j) = u(nx-1,j)
  utemp(nx+1,j) = u(1,j)
  vtemp(nx+1,j) = v(1,j)
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

	hx(i,ny)= 0.
    hx(i,0)= 0.

end do


return
end



!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!Subroutine for Pressure Updation
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine updateP(p,hx,hy,nx,ny,dx,dy)
implicit none

common/PoissonIter/psit

integer :: nx,ny,i,j,psit,k,q
real*8::dy,dx,ta,tb
real*8,dimension(0:nx,0:ny):: p,hx,hy,hxdx,hydy,ptemp1,hytemp,ptemp2
real*8,dimension(-1:nx+1,0:ny)::hxtemp,pd

psit = 500

!Calculating divergence of H terms


do i = 0,nx
  do j = 0,ny
    hxtemp(i,j) = hx(i,j)
    hytemp(i,j) = hy(i,j)
    ptemp1(i,j)	= p(i,j)
    ptemp2(i,j)	= p(i,j)    
    pd(i,j)		= p(i,j)
    hxdx(i,j)=0.0
    hydy(i,j)=0.0    
  end do
end do

!Adding Dummy terms in x direction for hxdx

do j = 0,ny
  hxtemp(-1,j) = hx(nx-1,j)
  hxtemp(nx+1,j) = hx(1,j)
  pd(-1,j) = p(nx-1,j)
  pd(nx+1,j)=p(1,j)
end do


    
do i = 0,nx
  do j = 1,ny-1
    hxdx(i,j) = (hxtemp(i+1,j)-hxtemp(i-1,j))/(2d0*dx)
    hydy(i,j) = (hytemp(i,j+1)-hytemp(i,j-1))/(2d0*dy)
  end do
end do    

q = 0
  
do k = 1,psit

do i = 1,nx-1
  do j = 1,ny-1
    if (j>1.and.j<ny-1) then
      ta = pd(i+1,j)+pd(i-1,j)+pd(i,j+1)+p(i,j-1)
	  tb = hxdx(i,j) + hydy(i,j)
      ptemp1(i,j) = ta/4-tb/4*(dx**2)
	else if (j==1) then
      ta = pd(i+1,j)+pd(i-1,j)+pd(i,j+1)
      tb = hxdx(i,j) + hydy(i,j)
      ptemp1(i,j) = ta/3-tb/3*(dx**2)
    else if (j==ny-1) then
      ta = pd(i+1,j)+pd(i-1,j)+pd(i,j-1)
      tb = hxdx(i,j) + hydy(i,j)
      ptemp1(i,j) = ta/3-tb/3*(dx**2)
	end if
  end do

	ptemp1(i,ny) = ptemp1(i,ny-1)
	ptemp1(i,0) = ptemp1(i,1)  
end do		  

call l1normcheck(ptemp1,ptemp2,nx,ny,q)

if (q.ne.1) then

do j = 0,ny
  do i = 0,nx
    	p(i,j) = ptemp1(i,j)
  end do
end do

print*,k

exit

end if

do j = 0,ny
  do i = 0,nx
    	ptemp2(i,j) = ptemp1(i,j)
  end do
end do





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
real*8,dimension(0:nx,0:ny)::utemp,u
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

common/tint/k,outpar
common/density/rho

real*8::dx,dy,dt,rho
integer::nx,ny,i,j,k,outpar
real*8,dimension(0:nx,0:ny)::u,v,p,hx,hy
real*8,dimension(0:nx,0:ny)::gradpx,gradpy
real*8,dimension(-1:nx+1,0:ny)::ptemp

!Defining dummys for pressure

do i = 0,nx
  do j = 0,ny
    ptemp(i,j) = p(i,j)
  end do
end do

do i=0,nx
  do j=0,ny
    gradpx(i,j)=0.0
    gradpy(i,j)=0.0
  end do
end do  

!Periodic BCs for gradpx

do j= 0,ny
  ptemp(-1,j) = ptemp(nx-1,j)
  ptemp(nx+1,j) = ptemp(1,j)
end do

!Calculation of pressure gradient
do i = 1,nx-1
  do j = 0,ny
    gradpx(i,j) = (ptemp(i+1,j)-ptemp(i-1,j))/(2.*dx)
  end do
end do

do i = 0,nx
  	do j = 1,ny-1
    	gradpy(i,j) = (ptemp(i,j+1)-ptemp(i,j-1))/(2d0*dy)
	end do
end do

  	do j = 0,ny
    	gradpx(0,j) = gradpx(1,j)
    	gradpx(nx,j) = gradpx(nx-1,j)        
	end do

  
if (mod(k,outpar)==0) then
open(20,file="GradPContourPlots.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',k
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.3,1600F14.3,1600F14.3,1600F14.3)',advance="no")dfloat(i)/dfloat(ny),dfloat(j)/dfloat(ny)&
      &,gradpx(i,j),gradpy(i,j)
      write(20,*) ''
    end do
  end do
close(20)
end if

!Updating the field
do i = 0,nx
  do j = 1,ny-1
	u(i,j) = u(i,j) + dt*(hx(i,j)-1/rho*gradpx(i,j))
  	v(i,j) = v(i,j) + dt*(hy(i,j)-1/rho*gradpy(i,j))
  end do
end do



return
end