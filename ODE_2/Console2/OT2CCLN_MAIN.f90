 real*16 function f1(x1,x2,x3,t) 
 real*16 x1,x2,x3,t
 real*16 xr 
xr=-(t+1)*x1+x2
!xr=2*x1+4*x2+exp(4*t)
f1=xr
end function

 real*16 function f2(x1,x2,x3,t) 
 real*16 x1,x2,x3,t
 real*16 xr 
xr=-1*x1*cos(t)+x3
!xr=x1+2*x2;
f2=xr
end function

 real*16 function f3(x1,x2,x3,t) 
 real*16 x1,x2,x3,t
 real*16 xr 
xr=t*t +x1*(x1*sin(x1+t)+(t*t-1))
f3=xr
end function



subroutine solve_euler(x_1_i, x_2_i, x_3_i,t,h)
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o
real*16, INTENT(IN):: t,h
 
   x_1_o=h*f1(x_1_i,x_2_i,x_3_i,t);
   x_2_o=h*f2(x_1_i,x_2_i,x_3_i,t);
   x_3_o=h*f3(x_1_i,x_2_i,x_3_i,t);
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o

end subroutine


subroutine solve_rk4(x_1_i, x_2_i, x_3_i,t1,h)
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o
real*16 c, b, a
real*16, INTENT(IN):: t1,h
real*16 hhalf,t

hhalf=h/2.0Q0;
t=t1-h;
 a=f1(x_1_i,x_2_i,x_3_i,t);
 b=f1(x_1_i+hhalf*a,x_2_i+hhalf*a,x_3_i+hhalf*a,t+hhalf);
 c=f1(x_1_i+hhalf*b,x_2_i+hhalf*b,x_3_i+hhalf*b,t+hhalf);
 d=f1(x_1_i+c,x_2_i+c,x_3_i+c,t+h);
 x_1_o=1.0Q0/6.0Q0*h*(a+b+c+d)

 a=f2(x_1_i,x_2_i,x_3_i,t);
 b=f2(x_1_i+hhalf*a,x_2_i+hhalf*a,x_3_i+hhalf*a,t+hhalf);
 c=f2(x_1_i+hhalf*b,x_2_i+hhalf*b,x_3_i+hhalf*b,t+hhalf);
 d=f2(x_1_i+c,x_2_i+c,x_3_i+c,t+h);
 x_2_o=1.0Q0/6.0Q0*h*(a+b+c+d)
 
  a=f3(x_1_i,x_2_i,x_3_i,t);
 b=f3(x_1_i+hhalf*a,x_2_i+hhalf*a,x_3_i+hhalf*a,t+hhalf);
 c=f3(x_1_i+hhalf*b,x_2_i+hhalf*b,x_3_i+hhalf*b,t+hhalf);
 d=f3(x_1_i+c,x_2_i+c,x_3_i+c,t+h);
 x_3_o=1.0Q0/6.0Q0*h*(a+b+c+d)
 
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o

end subroutine


subroutine writer(i,x_1_i, x_2_i, x_3_i,t)
10 format ((f16.6),',',(f16.6),',',(f16.6),',',(f16.6))
 real*16 x_1_i, x_2_i, x_3_i, t
 integer*8 i;
if (mod(i,5000) .eq. 0) then
write(100,10)t,x_1_i,x_2_i,x_3_i
end if 
end subroutine

PROGRAM ODE


 real*16 h,t
 real*16 etol
 real*16 x_1_i, x_2_i, x_3_i
 real*16 x_1_o, x_2_o, x_3_o
 real*16 x_1, x_2, x_3

integer*8 i



t=0
x_0=0.0Q0
x_1=1.1Q0  !1.1
x_2=2.2Q0 !2.2
x_3=3.3Q0

x_0_i=x_0
x_1_i=x_1
x_2_i=x_2
x_3_i=x_3

h=.000001Q0

open (100, FILE='datos.csv')
do i=1,10000000
  t=h*i
!  call solve_euler(x_1_i, x_2_i, x_3_i,t,h)
  call solve_rk4(x_1_i, x_2_i, x_3_i,t,h)
 
  call writer(i,x_1_i, x_2_i, x_3_i,t)
end do
close (100,STATUS='KEEP')
WRITE (*,*)'DONE'
!pause
end program ODE



