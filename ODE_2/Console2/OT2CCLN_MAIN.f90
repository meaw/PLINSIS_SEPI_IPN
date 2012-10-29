 real*16 function f1(x1,x2,x3,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t
 real*16 xr 
xr=-(t+1)*x1+x2
!xr=2*x1+4*x2+exp(4*t)
f1=xr
end function

 real*16 function f2(x1,x2,x3,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t
 real*16 xr 
xr=-1*x1*cos(t)+x3
!xr=x1+2*x2;
f2=xr
end function

 real*16 function f3(x1,x2,x3,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t
 real*16 xr 
xr=t*t +x1*(x1*sin(x1+t)+(t*t-1))
f3=xr
end function



subroutine solve_euler_trap(x_1_i, x_2_i, x_3_i,t,h)
IMPLICIT NONE
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o
real*16 x_1_t, x_2_t, x_3_t
real*16, INTENT(IN):: t,h
real*16 hhalf
hhalf=h/2.0Q0
 
   x_1_t=x_1_i+hhalf*f1(x_1_i,x_2_i,x_3_i,t);
   x_2_t=x_2_i+hhalf*f2(x_1_i,x_2_i,x_3_i,t);
   x_3_t=x_3_i+hhalf*f3(x_1_i,x_2_i,x_3_i,t);

 
   x_1_o=h*f1(x_1_t,x_2_t,x_3_t,t+hhalf);
   x_2_o=h*f2(x_1_t,x_2_t,x_3_t,t+hhalf);
   x_3_o=h*f3(x_1_t,x_2_t,x_3_t,t+hhalf);
   
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o

end subroutine

subroutine solve_euler(x_1_i, x_2_i, x_3_i,t,h)
IMPLICIT NONE
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

subroutine solve_euler_mod(x_1_i, x_2_i, x_3_i,t,h)   
IMPLICIT NONE
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o
real*16 x_1_t1, x_2_t1, x_3_t1
real*16, INTENT(IN):: t,h
 x_1_t1=f1(x_1_i,x_2_i,x_3_i,t);
 x_2_t1=f2(x_1_i,x_2_i,x_3_i,t);
  x_3_t1=f3(x_1_i,x_2_i,x_3_i,t);

 x_1_o=1.0Q0/2.0Q0*h*(x_1_t1+f1(x_1_i+x_1_t1*h ,x_2_i+x_2_t1*h,x_3_i+x_3_t1*h,t+h))
 

 x_2_o=1.0Q0/2.0Q0*h*(x_2_t1+f1(x_1_i+x_1_t1*h ,x_2_i+x_2_t1*h*h,x_3_i+x_3_t1*h,t+h)*h)
 

 x_3_o=1.0Q0/2.0Q0*h*(x_3_t1+f1(x_1_i+x_1_t1*h ,x_2_i+x_2_t1*h,x_3_i+x_3_t1*h,t+h)*h)
 
 

   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o

end subroutine




subroutine solve_trap(x_1_i, x_2_i, x_3_i,t,h)
IMPLICIT NONE
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp,y_k,xval
integer*8 max_iter,itern
max_iter=2000;

        tol=h/10000Q0;
        y_k=x_1_i;
        y_kp=x_1_i;
        do itern=1 , max_iter
            xval=y_k+h/2.0Q0*(f1(y_kp,x_2_i,x_3_i,t + h)+f1(x_1_i,x_2_i,x_3_i,t))
             if (abs(xval - y_kp)<tol) then
                    x_1_o=xval
                   goto 2
             end if
            y_kp=xval
        end do
        x_1_o=xval;
    2   y_k=x_2_i;
        y_kp=x_2_i;
        do itern=1 , max_iter
            xval=y_k+h/2.0Q0*(f2(x_1_i,y_kp,x_3_i,t + h)+f2(x_1_i,x_2_i,x_3_i,t))
             if (abs(xval - y_kp)<tol) then
                    x_2_o=xval
                   goto 3
             end if
            y_kp=xval
        end do
        x_2_o=xval
    3   y_k=x_3_i;
        y_kp=x_3_i;
        do itern=1 , max_iter
            xval=y_k+(h/2.0Q0)*(f3(x_1_i,x_2_i,y_kp,t + h)+f3(x_1_i,x_2_i,x_3_i,t))
             if (abs(xval - y_kp)<tol) then
                    x_3_o=xval
                   goto 4
             end if
            y_kp=xval
        end do
        x_3_o=xval
    4   x_1_i=x_1_o;
     x_2_i=x_2_o;
      x_3_i=x_3_o;
end subroutine

subroutine solve_euler_implicit(x_1_i, x_2_i, x_3_i,t,h)
IMPLICIT NONE
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp_1,y_kp_3,y_kp_2,y_k,xval,y_kp_1_o,y_kp_3_o,y_kp_2_o
integer*8 max_iter,itern,iternj,j
max_iter=2000;

        tol=h/10000Q0;
        
          y_kp_1=x_1_i;
        y_kp_2=x_2_i;
        y_kp_3=x_3_i;
        
        y_kp_1_o=y_kp_1+1
         y_kp_2_o=y_kp_1+1
          y_kp_3_o=y_kp_1+1
        do iternj=1 , max_iter
        
        
        
         if( (abs(y_kp_1_o - y_kp_1)<tol) .and. (abs(y_kp_2_o - y_kp_2)<tol) .and. (abs(y_kp_3_o - y_kp_3)<tol)) then
         
         j=j
         exit
         
         end if
        
         y_kp_1_o=y_kp_1
         y_kp_2_o=y_kp_2
         y_kp_3_o=y_kp_3
         
        
        y_k=x_1_i;
        do itern=1 , max_iter
            xval=y_k+h*(f1(y_kp_1,y_kp_2,y_kp_3,t + h))
             if (abs(xval - y_kp_1)<tol) then
             y_kp_1=xval
                   goto 2
             end if
            y_kp_1=xval
        end do
    2   y_k=x_2_i;
        do itern=1 , max_iter
            xval=y_k+h*(f2(y_kp_1,y_kp_2,y_kp_3,t + h))
             if (abs(xval - y_kp_2)<tol) then
              y_kp_2=xval
                   goto 3
             end if
            y_kp_2=xval
        end do
    3   y_k=x_3_i;
        do itern=1 , max_iter
            xval=y_k+h*(f3(y_kp_1,y_kp_2,y_kp_3,t + h))
             if (abs(xval - y_kp_3)<tol) then
              y_kp_3=xval
                   goto 4
             end if
            y_kp_3=xval
        end do
    4  j=j
    end do
     x_1_i=y_kp_1;
     x_2_i=y_kp_2;
      x_3_i=y_kp_3;
end subroutine












subroutine solve_am4(x_1_i, x_2_i, x_3_i,t,h,size)
IMPLICIT NONE
 real*16 f1,f2,f3
 integer size
real*16, INTENT(INOUT)::  x_1_i(size), x_2_i(size), x_3_i(size)
real*16 x_1_o, x_2_o, x_3_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp_1,y_kp_3,y_kp_2,y_k,xval,y_kp_1_o,y_kp_3_o,y_kp_2_o
integer*8 max_iter,itern,iternj,j
max_iter=2000;


        tol=h/10000Q0;
      
        y_kp_1=x_1_i(4);
        y_kp_2=x_2_i(4);
        y_kp_3=x_3_i(4);
        
        y_kp_1_o=y_kp_1+1
         y_kp_2_o=y_kp_1+1
          y_kp_3_o=y_kp_1+1
        do iternj=1 , max_iter
        
        
        
         if( (abs(y_kp_1_o - y_kp_1)<tol) .and. (abs(y_kp_2_o - y_kp_2)<tol) .and. (abs(y_kp_3_o - y_kp_3)<tol)) then
         
         j=j
         exit
         
         end if
        
         y_kp_1_o=y_kp_1
         y_kp_2_o=y_kp_2
         y_kp_3_o=y_kp_3
         
           y_k=x_1_i(4);
        do itern=1 , max_iter
            xval=y_k+h/24.0Q0*(9.0Q0*f1(y_kp_1,y_kp_2,y_kp_3,t + h)+19.0Q0*f1(x_1_i(4),x_2_i(4),x_3_i(4),t )-5.0Q0*f1(x_1_i(3),x_2_i(3),x_3_i(3),t - h)+1.0Q0*f1(x_1_i(2),x_2_i(2),x_3_i(2),t -h-h))
             if (abs(xval - y_kp_1)<tol) then
             y_kp_1=xval
                    x_1_o=xval
                   goto 2
             end if
            y_kp_1=xval
        end do
        x_1_o=xval;
    2   y_k=x_2_i(4);
        do itern=1 , max_iter
           xval=y_k+h/24.0Q0*(9.0Q0*f2(y_kp_1,y_kp_2,y_kp_3,t + h)+19.0Q0*f2(x_1_i(4),x_2_i(4),x_3_i(4),t )-5.0Q0*f2(x_1_i(3),x_2_i(3),x_3_i(3),t - h)+1.0Q0*f2(x_1_i(2),x_2_i(2),x_3_i(2),t -h-h))
            
             if (abs(xval - y_kp_2)<tol) then
             y_kp_2=xval
                    x_2_o=xval
                   goto 3
             end if
             y_kp_2=xval
        end do
        x_2_o=xval
    3   y_k=x_3_i(4);
        do itern=1 , max_iter
        
        
            xval=y_k+h/24.0Q0*(9.0Q0*f3(y_kp_1,y_kp_2,y_kp_3,t + h)+19.0Q0*f3(x_1_i(4),x_2_i(4),x_3_i(4),t )-5.0Q0*f3(x_1_i(3),x_2_i(3),x_3_i(3),t - h)+1.0Q0*f3(x_1_i(2),x_2_i(2),x_3_i(2),t -h-h))
       
             if (abs(xval - y_kp_3)<tol) then
             y_kp_3=xval
                    x_3_o=xval
                   goto 4
             end if
            y_kp_3=xval
        end do
        x_3_o=xval
    4 j=j
    
      end do
     call rotate(x_1_i, x_2_i, x_3_i,size)
    
     x_1_i(4)=x_1_o;
     x_2_i(4)=x_2_o;
      x_3_i(4)=x_3_o;
end subroutine









subroutine solve_rk4(x_1_i, x_2_i, x_3_i,t1,h)
IMPLICIT NONE
 real*16 f1,f2,f3
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i
real*16 x_1_o, x_2_o, x_3_o
real*16 c, b, a,d
real*16, INTENT(IN):: t1,h
real*16 hhalf,t

hhalf=h/2.0Q0;
t=t1;
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








subroutine rotate(x_1_i, x_2_i, x_3_i,size)
IMPLICIT NONE
integer size
real*16, INTENT(INOUT)::  x_1_i(size), x_2_i(size), x_3_i(size)
integer i;

do i=1,size-1
x_1_i(i)=x_1_i(i+1)
x_2_i(i)=x_2_i(i+1)
x_3_i(i)=x_3_i(i+1)
end do


end subroutine





subroutine writer(i,x_1_i, x_2_i, x_3_i,t)
IMPLICIT NONE
10 format ((f16.6),',',(f16.6),',',(f16.6),',',(f16.6))
 real*16 x_1_i, x_2_i, x_3_i, t
 integer*8 i;
if (mod(i,500) .eq. 0) then
write(100,10)t,x_1_i,x_2_i,x_3_i
end if 
end subroutine

PROGRAM ODE


 real*16 h,t
 real*16 etol
 real*16 x_1_i(4), x_2_i(4), x_3_i(4)
 real*16 x_1, x_2, x_3

integer*8 i



t=0
x_1=1.1Q0  !1.1
x_2=2.2Q0 !2.2
x_3=3.3Q0


x_1_i(4)=x_1
x_2_i(4)=x_2
x_3_i(4)=x_3

h=.00001Q0

open (100, FILE='datos.csv')
do i=1,1000000
  t=h*(i-1);
!  call solve_euler(x_1_i, x_2_i, x_3_i,t,h)
!  call solve_rk4(x_1_i, x_2_i, x_3_i,t,h)
 !call solve_euler_mod(x_1_i, x_2_i, x_3_i,t,h)    !!verificarlo
 call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4),t,h)
! call solve_euler_trap(x_1_i, x_2_i, x_3_i,t,h)
 !call solve_trap(x_1_i(4), x_2_i(4), x_3_i(4),t,h) 
 
 call writer(i,x_1_i(4), x_2_i(4), x_3_i(4),t)
end do

!  t=h*(0);
!  call rotate(x_1_i, x_2_i, x_3_i,4)
!call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4),t,h)
!
!  t=h*(1);
!  call rotate(x_1_i, x_2_i, x_3_i,4)
!call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4),t,h)
!
!  t=h*(2);
!  call rotate(x_1_i, x_2_i, x_3_i,4)
!call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4),t,h)
!
!
!do i=4,1000000
!  t=h*(i-1);
!call solve_am4(x_1_i, x_2_i, x_3_i,t,h,4) 
! 
!call writer(i,x_1_i(4), x_2_i(4), x_3_i(4),t)
!end do






close (100,STATUS='KEEP')
WRITE (*,*)'DONE'
!pause
end program ODE



