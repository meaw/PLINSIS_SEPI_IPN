
 real*16 function f1(x1,x2,x3,x4,x5,x6,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t,x4,x5,x6
 real*16 xr 
 
 integer*4 solving_set
common /solver_data/solving_set
xr=0;
 if (solving_set .eq. 1) then 
xr=x1*x1
end if

 if (solving_set .eq. 2) then 
xr=-(t+1.0Q0)*x1+x2
end if

 if (solving_set .eq. 3) then 
xr=-x1-2.0Q0*x2+t*t
end if

 if (solving_set .eq. 4) then 
xr=100.0Q0*(sin(t)-x1)
end if

 if (solving_set .eq. 5) then 
xr=x2
end if

!xr=2*x1+4*x2+exp(4*t)
f1=xr
end function

 real*16 function f2(x1,x2,x3,x4,x5,x6,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t,x4,x5,x6
 real*16 xr 
 
 integer*4 solving_set
common /solver_data/solving_set
xr=0;
 if (solving_set .eq. 2) then 
xr=-x1*Qcos(t)+x3
end if
 if (solving_set .eq. 3) then 
xr=1.0Q0/2.0Q0*x1+2.0Q0*x2+t
end if

 if (solving_set .eq. 5) then 
xr=-1.001Q0*x2-1.0Q0*x1
end if
!xr=x1+2*x2;
f2=xr
end function

 real*16 function f3(x1,x2,x3,x4,x5,x6,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t,x4,x5,x6
 real*16 xr 
 integer*4 solving_set
common /solver_data/solving_set
 xr=0;
 
  if (solving_set .eq. 2) then 
xr=t*t +x1*(x1*Qsin(x1+t)+(t*t-1.0Q0))
end if
f3=xr
end function

real*16 function f4(x1,x2,x3,x4,x5,x6,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t,x4,x5,x6
 real*16 xr 
 integer*4 solving_set
common /solver_data/solving_set
 xr=0;
 

f4=xr
end function

real*16 function f5(x1,x2,x3,x4,x5,x6,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t,x4,x5,x6
 real*16 xr 
 integer*4 solving_set
common /solver_data/solving_set
 xr=0;
 
f5=xr
end function

real*16 function f6(x1,x2,x3,x4,x5,x6,t) 
 IMPLICIT NONE
 real*16 x1,x2,x3,t,x4,x5,x6
 real*16 xr 
 integer*4 solving_set
common /solver_data/solving_set
 xr=0;
 

f6=xr
end function


subroutine solve_euler_trap(x_1_i, x_2_i, x_3_i,x_4_i,x_5_i,x_6_i,t,h)
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i, x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o
real*16 x_1_t, x_2_t, x_3_t, x_4_t, x_5_t, x_6_t
real*16, INTENT(IN):: t,h
real*16 hhalf
hhalf=h/2.0Q0
 
   x_1_t=x_1_i+hhalf*f1(x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i,t);
   x_2_t=x_2_i+hhalf*f2(x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i,t);
   x_3_t=x_3_i+hhalf*f3(x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i,t);
   x_4_t=x_4_i+hhalf*f4(x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i,t);
   x_5_t=x_5_i+hhalf*f5(x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i,t);
   x_6_t=x_6_i+hhalf*f6(x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i,t);
 
   x_1_o=h*f1(x_1_t,x_2_t,x_3_t,x_4_t,x_5_t,x_6_t,t+hhalf);
   x_2_o=h*f2(x_1_t,x_2_t,x_3_t,x_4_t,x_5_t,x_6_t,t+hhalf);
   x_3_o=h*f3(x_1_t,x_2_t,x_3_t,x_4_t,x_5_t,x_6_t,t+hhalf);
   x_4_o=h*f4(x_1_t,x_2_t,x_3_t,x_4_t,x_5_t,x_6_t,t+hhalf);
   x_5_o=h*f5(x_1_t,x_2_t,x_3_t,x_4_t,x_5_t,x_6_t,t+hhalf);
   x_6_o=h*f6(x_1_t,x_2_t,x_3_t,x_4_t,x_5_t,x_6_t,t+hhalf);
   
   
   
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o
   x_4_i=x_4_i+x_4_o
   x_5_i=x_5_i+x_5_o
   x_6_i=x_6_i+x_6_o
end subroutine

subroutine solve_euler(x_1_i, x_2_i, x_3_i, x_4_i, x_5_i, x_6_i,t,h)
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i, x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o, x_4_o, x_5_o, x_6_o
real*16, INTENT(IN):: t,h
 
   x_1_o=h*f1(x_1_i,x_2_i,x_3_i, x_4_i, x_5_i, x_6_i,t);
   x_2_o=h*f2(x_1_i,x_2_i,x_3_i, x_4_i, x_5_i, x_6_i,t);
   x_3_o=h*f3(x_1_i,x_2_i,x_3_i, x_4_i, x_5_i, x_6_i,t);
   x_4_o=h*f4(x_1_i,x_2_i,x_3_i, x_4_i, x_5_i, x_6_i,t);
   x_5_o=h*f5(x_1_i,x_2_i,x_3_i, x_4_i, x_5_i, x_6_i,t);
   x_6_o=h*f6(x_1_i,x_2_i,x_3_i, x_4_i, x_5_i, x_6_i,t);
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o
   x_4_i=x_4_i+x_4_o
   x_5_i=x_5_i+x_5_o
   x_6_i=x_6_i+x_6_o

end subroutine

subroutine solve_euler_modified(x_1_i, x_2_i,x_4_i, x_5_i, x_6_i, x_3_i,t,h)   
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o
real*16 x_1_t1, x_2_t1, x_3_t1, x_4_t1, x_5_t1, x_6_t1
real*16 x_1_t2, x_2_t2, x_3_t2,x_4_t2, x_5_t2, x_6_t2
real*16, INTENT(IN):: t,h
 x_1_t1=h*f1(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 x_2_t1=h*f2(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
  x_3_t1=h*f3(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
   x_4_t1=h*f4(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 x_5_t1=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
  x_6_t1=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
  

 x_1_t2=h*f1(x_1_i+x_1_t1,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+h);
 x_2_t2=h*f2(x_1_i,x_2_i+x_2_t1,x_3_i,x_4_i, x_5_i, x_6_i,t+h);
 x_3_t2=h*f3(x_1_i,x_2_i,x_3_i+x_3_t1,x_4_i, x_5_i, x_6_i,t+h);
 x_4_t2=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+x_4_t1, x_5_i, x_6_i,t+h);
 x_5_t2=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+x_5_t1, x_6_i,t+h);
 x_6_t2=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+x_6_t1,t+h);
  
  
 x_1_o=1.0Q0/2.0Q0*(x_1_t1+x_1_t2) 
 x_2_o=1.0Q0/2.0Q0*(x_2_t1+x_2_t2)
  x_3_o=1.0Q0/2.0Q0*(x_3_t1+x_3_t2)
  x_4_o=1.0Q0/2.0Q0*(x_4_t1+x_4_t2) 
 x_5_o=1.0Q0/2.0Q0*(x_5_t1+x_5_t2)
  x_6_o=1.0Q0/2.0Q0*(x_6_t1+x_6_t2)

   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o
   x_4_i=x_4_i+x_4_o
   x_5_i=x_5_i+x_5_o
   x_6_i=x_6_i+x_6_o

end subroutine


subroutine solve_trap(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t,h)


IMPLICIT NONE


 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o, x_4_o, x_5_o, x_6_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp,y_k,xval,y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,y_kp_1_o,y_kp_2_o,y_kp_3_o,y_kp_4_o,y_kp_5_o,y_kp_6_o,toll
integer*8 max_iter,itern,iternj,j


 integer*4 solving_set
common /solver_data/solving_set


max_iter=2000;
        tol=1e-9;
        toll=1e-7;
        y_kp_1=x_1_i
        y_kp_2=x_2_i
        y_kp_3=x_3_i
        y_kp_4=x_4_i
        y_kp_5=x_5_i
        y_kp_6=x_6_i
        
        y_kp_1_o=y_kp_1+1
         y_kp_2_o=y_kp_2+1
          y_kp_3_o=y_kp_3+1
          
            y_kp_4_o=y_kp_4+1
         y_kp_5_o=y_kp_5+1
          y_kp_6_o=y_kp_6+1
          
          
        do iternj=1 , 10
        
        
        
         if( (abs(y_kp_1_o - y_kp_1)<toll) .and. (abs(y_kp_2_o - y_kp_2)<toll) .and. (abs(y_kp_3_o - y_kp_3)<toll) .and. (abs(y_kp_4_o - y_kp_4)<toll) .and. (abs(y_kp_5_o - y_kp_5)<toll) .and. (abs(y_kp_6_o - y_kp_6)<toll)) then
         
         j=j
         exit
         
         end if
        
         y_kp_1_o=y_kp_1
         y_kp_2_o=y_kp_2
         y_kp_3_o=y_kp_3
         y_kp_4_o=y_kp_4
         y_kp_5_o=y_kp_5
         y_kp_6_o=y_kp_6


        y_k=x_1_i;
       
        do itern=1 , max_iter
            xval=y_k+h/2.0Q0*(f1(y_kp_1,y_kp_2,y_kp_3,y_kp_4, y_kp_5, y_kp_6,t + h)+f1(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t))
            if (solving_set==4) then 
          y_kp_1=(50*h*sin(t+h)+50*h*sin(t)-y_k*(50*h-1))/(50*h+1)
            xval=y_kp_1;
            end if
            
             if (abs(xval - y_kp_1)<tol) then
                  y_kp_1=xval
                   
                   goto 2
             end if
            y_kp_1=xval
        end do

        
    2   y_k=x_2_i;
       
        do itern=1 , max_iter
            xval=y_k+h/2.0Q0*(f2(y_kp_1,y_kp_2,y_kp_3,y_kp_4, y_kp_5, y_kp_6,t + h)+f2(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t))
             if (abs(xval - y_kp_2)<tol) then
              y_kp_2=xval
                   goto 3
             end if
             y_kp_2=xval
        end do

    3   y_k=x_3_i;
        do itern=1 , max_iter
            xval=y_k+(h/2.0Q0)*(f3(y_kp_1,y_kp_2,y_kp_3,y_kp_4, y_kp_5, y_kp_6,t + h)+f3(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t))
             if (abs(xval - y_kp_3)<tol) then
                    y_kp_3=xval
                   goto 4
             end if
            y_kp_3=xval
        end do
    4    y_k=x_4_i;
        do itern=1 , max_iter
            xval=y_k+(h/2.0Q0)*(f4(y_kp_1,y_kp_2,y_kp_3,y_kp_4, y_kp_5, y_kp_6,t + h)+f4(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t))
             if (abs(xval - y_kp_4)<tol) then
                    y_kp_4=xval
                   goto 5
             end if
            y_kp_4=xval
        end do
 5      y_k=x_5_i;
        do itern=1 , max_iter
            xval=y_k+(h/2.0Q0)*(f5(y_kp_1,y_kp_2,y_kp_3,y_kp_4, y_kp_5, y_kp_6,t + h)+f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t))
             if (abs(xval - y_kp_5)<tol) then
                    y_kp_5=xval
                   goto 6
             end if
            y_kp_5=xval
        end do
 6      y_k=x_6_i;
        do itern=1 , max_iter
            xval=y_k+(h/2.0Q0)*(f6(y_kp_1,y_kp_2,y_kp_3,y_kp_4, y_kp_5, y_kp_6,t + h)+f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t))
             if (abs(xval - y_kp_6)<tol) then
                    y_kp_6=xval
                   goto 7
             end if
            y_kp_6=xval
        end do
   7   j=j
   end do
     x_1_i=y_kp_1;
     x_2_i=y_kp_2;
      x_3_i=y_kp_3;
      x_4_i=y_kp_4;
     x_5_i=y_kp_5;
      x_6_i=y_kp_6;
end subroutine

subroutine solve_euler_implicit(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t,h)
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp_1,y_kp_3,y_kp_2,y_kp_4,y_kp_5,y_kp_6,y_k,xval,y_kp_1_o,y_kp_3_o,y_kp_2_o,y_kp_4_o,y_kp_5_o,y_kp_6_o,toll,tmp
integer*8 max_iter,itern,iternj,j

 integer*4 solving_set
common /solver_data/solving_set
max_iter=2000;
     tol=1e-8;
             toll=1e-5;
          y_kp_1=x_1_i;
        y_kp_2=x_2_i;
        y_kp_3=x_3_i;
        y_kp_4=x_4_i;
        y_kp_5=x_5_i;
        y_kp_6=x_6_i;
        
        y_kp_1_o=y_kp_1+1
         y_kp_2_o=y_kp_2+1
          y_kp_3_o=y_kp_3+1
          
            y_kp_4_o=y_kp_4+1
         y_kp_5_o=y_kp_5+1
          y_kp_6_o=y_kp_6+1
          
          
        do iternj=1 , 10
        
        
        
         if( (abs(y_kp_1_o - y_kp_1)<toll) .and. (abs(y_kp_2_o - y_kp_2)<toll) .and. (abs(y_kp_3_o - y_kp_3)<toll) .and. (abs(y_kp_4_o - y_kp_4)<toll) .and. (abs(y_kp_5_o - y_kp_5)<toll) .and. (abs(y_kp_6_o - y_kp_6)<toll)) then
         
         j=j
         exit
         
         end if
        
         y_kp_1_o=y_kp_1
         y_kp_2_o=y_kp_2
         y_kp_3_o=y_kp_3
         y_kp_4_o=y_kp_4
         y_kp_5_o=y_kp_5
         y_kp_6_o=y_kp_6
         
        
        y_k=x_1_i;
        do itern=1 , max_iter
           xval=y_k+h*(f1(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h))
         if (solving_set==4) then 
        y_kp_1=(100*h*sin(t+h)+y_k)/(100*h+1)
         xval=y_kp_1;
        end if
        
        
       
             if (abs(xval - y_kp_1)<tol) then
             y_kp_1=xval
                   goto 2
             end if
            y_kp_1=xval
        end do
    2   y_k=x_2_i;
        do itern=1 , max_iter
            xval=y_k+h*(f2(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h))
             if (abs(xval - y_kp_2)<tol) then
              y_kp_2=xval
                   goto 3
             end if
            y_kp_2=xval
        end do
    3   y_k=x_3_i;
        do itern=1 , max_iter
            xval=y_k+h*(f3(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h))
             if (abs(xval - y_kp_3)<tol) then
              y_kp_3=xval
                   goto 4
             end if
            y_kp_3=xval
        end do
    4  y_k=x_4_i;
        do itern=1 , max_iter
            xval=y_k+h*(f4(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h))
             if (abs(xval - y_kp_4)<tol) then
              y_kp_4=xval
                   goto 5
             end if
            y_kp_4=xval
        end do
    5  y_k=x_5_i;
        do itern=1 , max_iter
            xval=y_k+h*(f5(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h))
             if (abs(xval - y_kp_5)<tol) then
              y_kp_5=xval
                   goto 6
             end if
            y_kp_5=xval
        end do  
    6y_k=x_6_i;
        do itern=1 , max_iter
            xval=y_k+h*(f6(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h))
             if (abs(xval - y_kp_6)<tol) then
              y_kp_6=xval
                   goto 7
             end if
            y_kp_6=xval
        end do
    7 j=j
    end do
     x_1_i=y_kp_1;
     x_2_i=y_kp_2;
      x_3_i=y_kp_3;
      x_4_i=y_kp_4;
     x_5_i=y_kp_5;
      x_6_i=y_kp_6;
      
end subroutine


subroutine solve_am4(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t,h,size)
IMPLICIT NONE
 real*16 f1,f2,f3,f6,f4,f5
 integer size
real*16, INTENT(INOUT)::  x_1_i(size), x_2_i(size), x_3_i(size),x_4_i(size), x_5_i(size), x_6_i(size)
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp_1,y_kp_3,y_kp_2,y_kp_4,y_kp_5,y_kp_6,y_k,xval,y_kp_1_o,y_kp_3_o,y_kp_2_o,y_kp_4_o,y_kp_5_o,y_kp_6_o
integer*8 max_iter,itern,iternj,j

 integer*4 solving_set
common /solver_data/solving_set

max_iter=2000;


            tol=1e-9;
      
        y_kp_1=x_1_i(4);
        y_kp_2=x_2_i(4);
        y_kp_3=x_3_i(4);
        y_kp_4=x_4_i(4);
        y_kp_5=x_5_i(4);
        y_kp_6=x_6_i(4);
            
        y_kp_1_o=y_kp_1+1
         y_kp_2_o=y_kp_2+1
          y_kp_3_o=y_kp_3+1
          y_kp_4_o=y_kp_4+1
         y_kp_5_o=y_kp_5+1
          y_kp_6_o=y_kp_6+1
        do iternj=1 , 10
        
        
        
          if( (abs(y_kp_1_o - y_kp_1)<tol) .and. (abs(y_kp_2_o - y_kp_2)<tol) .and. (abs(y_kp_3_o - y_kp_3)<tol) .and. (abs(y_kp_4_o - y_kp_4)<tol) .and. (abs(y_kp_5_o - y_kp_5)<tol) .and. (abs(y_kp_6_o - y_kp_6)<tol)) then
        
         j=j
         exit
         
         end if
        
         y_kp_1_o=y_kp_1
         y_kp_2_o=y_kp_2
         y_kp_3_o=y_kp_3
         
         y_kp_4_o=y_kp_4
         y_kp_5_o=y_kp_5
         y_kp_6_o=y_kp_6
         
         
           y_k=x_1_i(4);
        do itern=1 , max_iter
            xval=y_k+h/24.0Q0*(9.0Q0*f1(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+19.0Q0*f1(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-5.0Q0*f1(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)+1.0Q0*f1(x_1_i(2),x_2_i(2),x_3_i(2),x_4_i(2),x_5_i(2),x_6_i(2),t -h-h))
         
         if (solving_set==4) then
                  y_kp_1=-(125*h*sin(t-h)-25*h*sin(t-2*h)-225*h*sin(t+h)-475*h*sin(t)+y_k*(475*h-6)-25*(5*x_1_i(3)-x_1_i(2))*h)/(3*(75*h+2));
                 
                    xval=y_kp_1
                    end if
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
           xval=y_k+h/24.0Q0*(9.0Q0*f2(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+19.0Q0*f2(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-5.0Q0*f2(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)+1.0Q0*f2(x_1_i(2),x_2_i(2),x_3_i(2),x_4_i(2),x_5_i(2),x_6_i(2),t -h-h))
            
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
            xval=y_k+h/24.0Q0*(9.0Q0*f3(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+19.0Q0*f3(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-5.0Q0*f3(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)+1.0Q0*f3(x_1_i(2),x_2_i(2),x_3_i(2),x_4_i(2),x_5_i(2),x_6_i(2),t -h-h))
       
             if (abs(xval - y_kp_3)<tol) then
             y_kp_3=xval
                    x_3_o=xval
                   goto 4
             end if
            y_kp_3=xval
        end do
        x_3_o=xval
    4  y_k=x_4_i(4);
        do itern=1 , max_iter
        
        
            xval=y_k+h/24.0Q0*(9.0Q0*f4(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+19.0Q0*f4(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-5.0Q0*f4(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)+1.0Q0*f4(x_1_i(2),x_2_i(2),x_3_i(2),x_4_i(2),x_5_i(2),x_6_i(2),t -h-h))
       
             if (abs(xval - y_kp_4)<tol) then
             y_kp_4=xval
                    x_4_o=xval
                   goto 5
             end if
            y_kp_4=xval
        end do
        x_4_o=xval
    5 y_k=x_5_i(4);
        do itern=1 , max_iter
            xval=y_k+h/24.0Q0*(9.0Q0*f5(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+19.0Q0*f5(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-5.0Q0*f5(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)+1.0Q0*f5(x_1_i(2),x_2_i(2),x_3_i(2),x_4_i(2),x_5_i(2),x_6_i(2),t -h-h))
             if (abs(xval - y_kp_5)<tol) then
             y_kp_5=xval
                    x_5_o=xval
                   goto 6
             end if
            y_kp_5=xval
        end do
        x_5_o=xval
      6  y_k=x_6_i(4);
        do itern=1 , max_iter
            xval=y_k+h/24.0Q0*(9.0Q0*f6(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+19.0Q0*f6(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-5.0Q0*f6(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)+1.0Q0*f6(x_1_i(2),x_2_i(2),x_3_i(2),x_4_i(2),x_5_i(2),x_6_i(2),t -h-h))
       
             if (abs(xval - y_kp_6)<tol) then
             y_kp_6=xval
                    x_6_o=xval
                   goto 7
             end if
            y_kp_6=xval
        end do
        x_6_o=xval 
        7 j=j
    
      end do
     call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,size)
    
     x_1_i(4)=x_1_o;
     x_2_i(4)=x_2_o;
      x_3_i(4)=x_3_o;
           x_4_i(4)=x_4_o;
     x_5_i(4)=x_5_o;
      x_6_i(4)=x_6_o;
end subroutine


subroutine solve_gear2(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t,h,size)
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
 integer size
real*16, INTENT(INOUT)::  x_1_i(size), x_2_i(size), x_3_i(size),x_4_i(size), x_5_i(size), x_6_i(size)
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o

real*16, INTENT(IN):: t,h
real*16 tol,y_kp_1,y_kp_3,y_kp_2,y_kp_4,y_kp_5,y_kp_6,y_k,xval,y_kp_1_o,y_kp_3_o,y_kp_2_o,y_kp_4_o,y_kp_5_o,y_kp_6_o
integer*8 max_iter,itern,iternj,j
max_iter=2000;

      
        tol=h/100000Q0;
      
         y_kp_1=x_1_i(4);
        y_kp_2=x_2_i(4);
        y_kp_3=x_3_i(4);
        y_kp_4=x_4_i(4);
        y_kp_5=x_5_i(4);
        y_kp_6=x_6_i(4);
            
        y_kp_1_o=y_kp_1+1
         y_kp_2_o=y_kp_2+1
          y_kp_3_o=y_kp_3+1
          y_kp_4_o=y_kp_4+1
         y_kp_5_o=y_kp_5+1
          y_kp_6_o=y_kp_6+1
        do iternj=1 , max_iter
        
        
        
          if( (abs(y_kp_1_o - y_kp_1)<tol) .and. (abs(y_kp_2_o - y_kp_2)<tol) .and. (abs(y_kp_3_o - y_kp_3)<tol) .and. (abs(y_kp_4_o - y_kp_4)<tol) .and. (abs(y_kp_5_o - y_kp_5)<tol) .and. (abs(y_kp_6_o - y_kp_6)<tol)) then
        
         j=j
         exit
         
         end if
        
        
        

        
         y_kp_1_o=y_kp_1
         y_kp_2_o=y_kp_2
         y_kp_3_o=y_kp_3
         y_kp_4_o=y_kp_4
         y_kp_5_o=y_kp_5
         y_kp_6_o=y_kp_6
         
           y_k=x_1_i(4);
        do itern=1 , max_iter
            xval=h*(2.0Q0/3.0Q0)*f1(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+(4.0Q0/3.0Q0)*f1(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-(1.0Q0/3.0Q0)*f1(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)
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
            xval=h*(2.0Q0/3.0Q0)*f2(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+(4.0Q0/3.0Q0)*f2(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-(1.0Q0/3.0Q0)*f2(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)
            
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
        
        
            xval=(2.0Q0/3.0Q0)*h*f3(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+(4.0Q0/3.0Q0)*f3(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-(1.0Q0/3.0Q0)*f3(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)
       
             if (abs(xval - y_kp_3)<tol) then
             y_kp_3=xval
                    x_3_o=xval
                   goto 4
             end if
            y_kp_3=xval
        end do
        x_3_o=xval
    4 y_k=x_4_i(4);
        do itern=1 , max_iter
        
        
            xval=(2.0Q0/3.0Q0)*h*f4(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+(4.0Q0/3.0Q0)*f4(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-(1.0Q0/3.0Q0)*f4(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)
       
             if (abs(xval - y_kp_4)<tol) then
             y_kp_4=xval
                    x_4_o=xval
                   goto 5
             end if
            y_kp_4=xval
        end do
        x_4_o=xval
     5   y_k=x_5_i(4);
        do itern=1 , max_iter
        
        
            xval=(2.0Q0/3.0Q0)*h*f5(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+(4.0Q0/3.0Q0)*f5(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-(1.0Q0/3.0Q0)*f5(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)
       
             if (abs(xval - y_kp_5)<tol) then
             y_kp_5=xval
                    x_5_o=xval
                   goto 6
             end if
            y_kp_5=xval
        end do
        x_5_o=xval
      6 y_k=x_6_i(4);
        do itern=1 , max_iter
        
        
            xval=(2.0Q0/3.0Q0)*h*f6(y_kp_1,y_kp_2,y_kp_3,y_kp_4,y_kp_5,y_kp_6,t + h)+(4.0Q0/3.0Q0)*f6(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),t )-(1.0Q0/3.0Q0)*f6(x_1_i(3),x_2_i(3),x_3_i(3),x_4_i(3),x_5_i(3),x_6_i(3),t - h)
       
             if (abs(xval - y_kp_6)<tol) then
             y_kp_6=xval
                    x_6_o=xval
                   goto 7
             end if
            y_kp_6=xval
        end do
        x_6_o=xval
     7   j=j
    
      end do
     call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,size)
    
     x_1_i(4)=x_1_o;
     x_2_i(4)=x_2_o;
      x_3_i(4)=x_3_o;
           x_4_i(4)=x_4_o;
     x_5_i(4)=x_5_o;
      x_6_i(4)=x_6_o;
end subroutine



subroutine solve_rkm(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t1,h)
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o
real*16 c, b, a,d,e
real*16, INTENT(IN):: t1,h
real*16 half,t

half=1/2.0Q0;
t=t1;
 a=h*f1(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f1(x_1_i+1.0Q0/3.0Q0*a,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 c=h*f1(x_1_i+1.0Q0/6.0Q0*(a+b),x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 d=h*f1(x_1_i+1.0Q0/8.0Q0*(a+3.0Q0*c),x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+1.0Q0/2.0Q0*h);
 e=h*f1(x_1_i+1.0Q0/2.0Q0*(a-3.0Q0*c+4*d),x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+h);
 
 x_1_o=1.0Q0/6.0Q0*(a+4.0Q0*d+e)

 a=h*f2(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f2(x_1_i,x_2_i+1.0Q0/3.0Q0*a,x_3_i,x_4_i, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 c=h*f2(x_1_i,x_2_i+1.0Q0/6.0Q0*(a+b),x_3_i,x_4_i, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 d=h*f2(x_1_i,x_2_i+1.0Q0/8.0Q0*(a+3.0Q0*c),x_3_i,x_4_i, x_5_i, x_6_i,t+1.0Q0/2.0Q0*h);
 e=h*f2(x_1_i,x_2_i+1.0Q0/2.0Q0*(a-3.0Q0*c+4*d),x_3_i,x_4_i, x_5_i, x_6_i,t+h);
  x_2_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
 
  a=h*f3(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f3(x_1_i,x_2_i,x_3_i+1.0Q0/3.0Q0*a,x_4_i, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 c=h*f3(x_1_i,x_2_i,x_3_i+1.0Q0/6.0Q0*(a+b),x_4_i, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 d=h*f3(x_1_i,x_2_i,x_3_i+1.0Q0/8.0Q0*(a+3.0Q0*c),x_4_i, x_5_i, x_6_i,t+1.0Q0/2.0Q0*h);
 e=h*f3(x_1_i,x_2_i,x_3_i+1.0Q0/2.0Q0*(a-3.0Q0*c+4*d),x_4_i, x_5_i, x_6_i,t+h);
  x_3_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
  
  
  a=h*f4(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+1.0Q0/3.0Q0*a, x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 c=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+1.0Q0/6.0Q0*(a+b), x_5_i, x_6_i,t+1.0Q0/3.0Q0*h);
 d=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+1.0Q0/8.0Q0*(a+3.0Q0*c), x_5_i, x_6_i,t+1.0Q0/2.0Q0*h);
 e=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+1.0Q0/2.0Q0*(a-3.0Q0*c+4*d), x_5_i, x_6_i,t+h);
 
 x_4_o=1.0Q0/6.0Q0*(a+4.0Q0*d+e)
 
 a=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+1.0Q0/3.0Q0*a, x_6_i,t+1.0Q0/3.0Q0*h);
 c=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+1.0Q0/6.0Q0*(a+b), x_6_i,t+1.0Q0/3.0Q0*h);
 d=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+1.0Q0/8.0Q0*(a+3.0Q0*c), x_6_i,t+1.0Q0/2.0Q0*h);
 e=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+1.0Q0/2.0Q0*(a-3.0Q0*c+4*d), x_6_i,t+h);
 
 x_5_o=1.0Q0/6.0Q0*(a+4.0Q0*d+e)
 
 a=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+1.0Q0/3.0Q0*a,t+1.0Q0/3.0Q0*h);
 c=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+1.0Q0/6.0Q0*(a+b),t+1.0Q0/3.0Q0*h);
 d=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+1.0Q0/8.0Q0*(a+3.0Q0*c),t+1.0Q0/2.0Q0*h);
 e=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+1.0Q0/2.0Q0*(a-3.0Q0*c+4*d),t+h);
 
 x_6_o=1.0Q0/6.0Q0*(a+4.0Q0*d+e)
 
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o
   
   x_4_i=x_4_i+x_4_o
   x_5_i=x_5_i+x_5_o
   x_6_i=x_6_i+x_6_o

end subroutine



subroutine solve_gear1(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t1,h)
IMPLICIT NONE
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i
real*16, INTENT(IN):: t1,h
call solve_euler_implicit(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t1,h)
end subroutine



subroutine solve_rk4(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t1,h)
IMPLICIT NONE
 real*16 f1,f2,f3,f4,f5,f6
real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i
real*16 x_1_o, x_2_o, x_3_o,x_4_o, x_5_o, x_6_o
real*16 c, b, a,d
real*16, INTENT(IN):: t1,h
real*16 half,t

half=1/2.0Q0;
t=t1;
 a=h*f1(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f1(x_1_i+half*a,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+half*h);
 c=h*f1(x_1_i+half*b,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+half*h);
 d=h*f1(x_1_i+c,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t+h);
 x_1_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)

 a=h*f2(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f2(x_1_i,x_2_i+half*a,x_3_i,x_4_i, x_5_i, x_6_i,t+half*h);
 c=h*f2(x_1_i,x_2_i+half*b,x_3_i,x_4_i, x_5_i, x_6_i,t+half*h);
 d=h*f2(x_1_i,x_2_i+c,x_3_i,x_4_i, x_5_i, x_6_i,t+h);
 x_2_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
 
 a=h*f3(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f3(x_1_i,x_2_i,x_3_i+half*a,x_4_i, x_5_i, x_6_i,t+half*h);
 c=h*f3(x_1_i,x_2_i,x_3_i+half*b,x_4_i, x_5_i, x_6_i,t+half*h);
 d=h*f3(x_1_i,x_2_i,x_3_i+c,x_4_i, x_5_i, x_6_i,t+h);
 x_3_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
 
 a=h*f4(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+half*a, x_5_i, x_6_i,t+half*h);
 c=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+half*b, x_5_i, x_6_i,t+half*h);
 d=h*f4(x_1_i,x_2_i,x_3_i,x_4_i+c, x_5_i, x_6_i,t+h);
 x_4_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
 
  a=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+half*a, x_6_i,t+half*h);
 c=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+half*b, x_6_i,t+half*h);
 d=h*f5(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i+c, x_6_i,t+h);
 x_5_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
 
  a=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i,t);
 b=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+half*a,t+half*h);
 c=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+half*b,t+half*h);
 d=h*f6(x_1_i,x_2_i,x_3_i,x_4_i, x_5_i, x_6_i+c,t+h);
 x_6_o=1.0Q0/6.0Q0*(a+2.0Q0*b+2.0Q0*c+d)
 
   x_1_i=x_1_i+x_1_o
   x_2_i=x_2_i+x_2_o
   x_3_i=x_3_i+x_3_o
   x_4_i=x_4_i+x_4_o
   x_5_i=x_5_i+x_5_o
   x_6_i=x_6_i+x_6_o
end subroutine








subroutine rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,size)
IMPLICIT NONE
integer size
real*16, INTENT(INOUT)::  x_1_i(size), x_2_i(size), x_3_i(size),x_4_i(size), x_5_i(size), x_6_i(size)
integer i;

do i=1,size-1
x_1_i(i)=x_1_i(i+1)
x_2_i(i)=x_2_i(i+1)
x_3_i(i)=x_3_i(i+1)
x_4_i(i)=x_4_i(i+1)
x_5_i(i)=x_5_i(i+1)
x_6_i(i)=x_6_i(i+1)
end do


end subroutine


subroutine set_ini(x_1_i, x_2_i, x_3_i, x_4_i, x_5_i, x_6_i,x_1, x_2, x_3,x_4, x_5, x_6)
IMPLICIT NONE

real*16, INTENT(INOUT)::  x_1_i, x_2_i, x_3_i, x_4_i, x_5_i, x_6_i
real*16, INTENT(INOUT)::  x_1, x_2, x_3,x_4, x_5, x_6
x_1_i=x_1
x_2_i=x_2
x_3_i=x_3
x_4_i=x_4
x_5_i=x_5
x_6_i=x_6

end subroutine


subroutine writer(tk,i,x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t)
IMPLICIT NONE
10 format ((f16.6),',',(f16.6),',',(f16.6),',',(f16.6),',',(f16.6),',',(f16.6),',',(f16.6))

 real*16 x_1_i, x_2_i, x_3_i, t,x_4_i, x_5_i, x_6_i
 integer*8 i,titer,ttime,tk;
 if (tk<1000) then
 write(100,10)t,x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i
 else
 
 
if (mod(i,tk/1000) .eq. 0) then
write(100,10)t,x_1_i,x_2_i,x_3_i,x_4_i,x_5_i,x_6_i
end if 
end if
end subroutine

PROGRAM ODE
IMPLICIT NONE

 real*16 h,t
 real*16 etol
 real*16 x_1_i(4), x_2_i(4), x_3_i(4),x_4_i(4), x_5_i(4), x_6_i(4)
 real*16 x_1, x_2, x_3,x_4, x_5, x_6

integer*8 i,titer,ttime
integer*4 solving_set
character  root *20

common /solver_data/solving_set


t=0

x_1=0.0Q0  !1.1
x_2=0.0Q0 !2.2
x_3=0.0Q0
x_4=0.0Q0  !1.1
x_5=0.0Q0 !2.2
x_6=0.0Q0
ttime=0.0Q0
!!Cambiar Estos datos de acuerdo al problema a resolver
solving_set=5
h=0.0175Q0
root='P5\0,0175\'



if (solving_set==2) then
x_1=1.1Q0  !1.1
x_2=2.2Q0 !2.2
x_3=3.3Q0
ttime=10.0Q0
end if

if (solving_set==1) then
x_1=0.5Q0  !1.1
ttime=1.0Q0
end if


if (solving_set==3) then
x_1=-2.0Q0  !1.1
x_2=1.0Q0  !1.1
ttime=2.0Q0
end if
if (solving_set==4) then
x_1=0.0Q0  !1.1}
ttime=4.0Q0
end if
if (solving_set==5) then
x_1=1.0Q0  !1.1
x_2=-1.0Q0  !1.1
ttime=3.0Q0
end ifr


titer=ttime/h;

Write(*,*) 'Solving Euler'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_euler.csv')
do i=1,titer
  t=h*(i-1);

  call solve_euler(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
 call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')

Write(*,*) 'Solving RK4'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_rk4.csv')
do i=1,titer
  t=h*(i-1);
call solve_rk4(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
 call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')


Write(*,*) 'Solving Euler Mod'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_euler_mod.csv')
do i=1,titer
  t=h*(i-1);
call solve_euler_modified(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)   
 call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')

Write(*,*) 'Solving Euler Implicit'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_euler_implicit.csv')
do i=1,titer
  t=h*(i-1);
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
 call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')

Write(*,*) 'Solving Euler Pre- Trap Corr'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_euler_trap.csv')
do i=1,titer
  t=h*(i-1);
 call solve_euler_trap(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
 call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')

Write(*,*) 'Solving Trap'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_trap.csv')
do i=1,titer
  t=h*(i-1);
call solve_trap(x_1_i(4), x_2_i(4) , x_3_i(4) ,x_4_i(4), x_5_i(4), x_6_i(4),t,h) 
 call writer(titer,i,x_1_i(4), x_2_i(4),   x_3_i(4),x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')

Write(*,*) 'Solving RKM'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_rkm.csv')
do i=1,titer
  t=h*(i-1);
call solve_rkm(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
 call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')



Write(*,*) 'Solving AM4'
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_am4.csv')
  t=h*(0);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
t=h*(1);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
t=h*(2);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
do i=4,titer
  t=h*(i-1);
call solve_am4(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,t,h,4)  
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')




Write(*,*) 'Solving Gear1 =Euler Impl '
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_Gear1.csv')
  t=h*(0);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
t=h*(1);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
t=h*(2);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
do i=4,titer
  t=h*(i-1);
 call solve_gear1(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h) 
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')

Write(*,*) 'Solving Gear2  '
call set_ini(x_1_i(4),x_2_i(4),x_3_i(4),x_4_i(4),x_5_i(4),x_6_i(4),x_1,x_2,x_3,x_4,x_5,x_6)
open (100, FILE=trim(root)//'datos_Gear2.csv')
  t=h*(0);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
t=h*(1);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
t=h*(2);
call rotate(x_1_i, x_2_i, x_3_i,x_4_i, x_5_i, x_6_i,4)
call solve_euler_implicit(x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t,h)
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
do i=4,titer
  t=h*(i-1);
 !call solve_gear2(x_1_i, x_2_i, x_3_i,t,h,4) 
call writer(titer,i,x_1_i(4), x_2_i(4), x_3_i(4)  ,x_4_i(4), x_5_i(4), x_6_i(4),t)
end do
close (100,STATUS='KEEP')















WRITE (*,*)'DONE'
!pause
end program ODE



