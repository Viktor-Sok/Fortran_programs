!Sokolow V.S.
!33413/1
module grids
implicit none
real,parameter:: Rmin=1, Rmax=11, time=10  !пространственно-временные данные задачи<<<<<<<<<<<
integer,parameter:: n=160                 !число разбиений по r <<<<<<<<<
end module grids

program lab1
use grids
use numerical_libraries
implicit none
integer, parameter ::mxparm=50            
integer::k,i,s,ido,istep,nout

real::u,h,yprime,r1,t,fi,param(mxparm),tol,tend,ut,d,md,utoch
dimension u(n+1),yprime(n+1),r1(n+1),d(n+1),utoch(n+1)
external fcn
s=n+1                !число дифф. уравнений
h=(Rmax-Rmin)/n      !шаг

r1(1)=Rmin
do i=1,n             !вычисление основной сетки
r1(i+1)=h+r1(i)
end do

open (1,file='result.txt')
write(1,3) (r1(i),i=1,n+1)
3 format ('istep',3X,'time',7X,'r=',<n+1>(1X,F6.2)) 

do i=1,n+1           !задание начальных условий аппроксимирующей системы
u(i)=fi(r1(i))
end do
tol=1.0e-4           !допустимая локальная погрешность
param=0.0            !действуем по умолчанию
param(10)=1.0e-5     !допустимая глобальная погрешность
param(4)=1000000     !допустимое чилсо шагов  
istep=0              
ido=1

t=0.0
do while (istep<time)
istep=istep+1
tend=istep
call ivprk (ido,s,fcn,t,tend,tol,param,u)
write (1,4) istep,t,u
4 format (I2,F12.3,7X,<n+1>(1X,F6.2)) 
end do
  
do i=1,n+1
utoch(i)=ut(r1(i),t)                 !модуль невязки точного и приближенного решения
d(i)=abs(utoch(i)-u(i))
end do
   
md=d(1)                    !Чёбышевская норма (максимум модуля невязки на отрезке интегрирования)  
do i=2,n+1
if (md<d(i)) md=d(i)
end do 
write(1,5) utoch
5 format (1X,'Точное решение=',6X,<n+1>(1X,F6.2)) 
write (1,*) 'Чёбышевская норма=',md
write (1,*) 'Колличесвто шагов= ', param(34)
write (1,*) 'Колличество вычесленных производных= ', param(35)
write (1,*) 'Величина последнего шага интегрирования ', param(31)

end program lab1



subroutine fcn(s,t,u,yprime)   
use grids          
implicit none
integer:: s,i
real :: h,r1,r2,f,k,u,q,t,yprime,pv1,pv2
dimension r1(s),r2(s-1),u(s),yprime(s)
!заполнение массивов основной и вспомогательной сетки
h=(Rmax-Rmin)/n     !шаг

r1(1)=Rmin
do i=1,n            !вычисление основной сетки
r1(i+1)=h+r1(i)
end do

r2(1)=Rmin+h/2
do i=1,n-1          !вычисление вспомогательной сетки
r2(i+1)=h+r2(i)
end do

yprime(1)=pv1(t)  !вычисление правой части системы
yprime(s)=pv2(t)
do i=2,s-1
yprime(i)=r2(i)/r1(i)*k(r2(i),t)*(u(i+1)-u(i))/h**2-r2(i-1)/r1(i)*k(r2(i-1),t)*(u(i)-u(i-1))/h**2-q(r1(i),t)*u(i)+f(r1(i),t)
end do
end subroutine fcn



real function k(r,t)                 ! Функция k
real:: r,t
k=t+r
end  function k


real function  q(r,t)                 ! Функция q
real::r,t
q=t+r
end function  q


real function  f(r,t)                 ! Функция f
real::r,t
f=t*(2-9*r+r**3+t**2+r*t)-12*r**2+r**4
end function f


real function pv1(t)                 ! Функция pv1 - производная от функции на границе Rmin
real::t
pv1=2*t
end function pv1


real function pv2(t)                 ! Функция pv2 - производная от функции на границе Rmax
real::t
pv2=2*t
end function pv2


real function fi(r)                 ! Функция fi (начальное условие)
real::r
fi=r**3
end function fi
 

real function ut(r,t)               ! точная функция-решение
real::r,t
ut=r**3+t**2
end function ut








