integer function get_solver_type(sol_name)
character sol_name *(*)

if (trim(sol_name) .eq. "EULER-IMPLICIT") then
get_solver_type=1
end if

if (trim(sol_name) .eq. "TRAPEZOIDAL") then

get_solver_type=2
end if

if (trim(sol_name) .eq. "SIMPSON") then
get_solver_type=3
end if

if (trim(sol_name) .eq. "ADAMS-BASHFORTH") then
get_solver_type=4
end if


Write(*,*)"Solver:",trim(sol_name),get_solver_type
end function


REAL*8 function fode_eval(x,y)
real*8 y,x
real*8 val
val=y*y
fode_eval=val
end function



real*8 function funceval( x)
    !funciones
    real*8 fode_eval
    !Argumenos de entrada
    real*8 x
    !Argumenos locales
    real*8 y_k,y_kp,t,x1

    !Argumenos Globales
    integer degree,solv_type
    real*8 h,tol,initial_val,values(10)
    COMMON /config/ h,tol,initial_val,values,degree,solv_type
        !Algunas condiciones
        y_k=initial_val
        y_kp=x
        t=0
        x1=.1
        !Solvers
        if (solv_type .eq. 1) then
            funceval=y_k+h*(fode_eval(x1,y_kp))-y_kp
        end if
        if (solv_type .eq. 2) then
            funceval=y_k+h/2*(fode_eval(x1,y_k)+fode_eval(x1,y_kp))-y_kp
        end if
        
        if (solv_type .eq. 3) then
            funceval=values(2)+h/2*(3*fode_eval(x1,values(2))+fode_eval(x1,values(1)) )-y_kp
        end if
end function

real*8 function funcderiv( x)
    real*8 x
    real*8 y0,y1,y2
    real*8 step
    real*8 funceval

    !Argumenos Globales
    integer degree,solv_type
    real*8 h,tol,initial_val,values(10)
    COMMON /config/ h,tol,initial_val,values,degree,solv_type
        step=.000001
        y0=funceval( x-step)
        y1=funceval( x)
        y2=funceval( x+step)
        y0=(y1-y0)/(+step)
        y2=(y2-y1)/(+step)
        y1=(y0+y2)/2
        funcderiv=y1
end function

real*8 function nsolver_svar(xini_,max_iter)
    !funciones
    real*8 funceval
    real*8 funcderiv
    !variables
    real*8 xini_,xini,xval 
    integer max_iter
    integer itern
    !Argumenos Globales
    integer degree,solv_type
    real*8 h,tol,initial_val,values(10)
    COMMON /config/ h,tol,initial_val,values,degree,solv_type
    
        xini=xini_
        do itern=1 , max_iter
            xval=xini-(funceval(xini)/ funcderiv(xini))
            if (abs(xval-xini)<tol) then
                nsolver_svar=xval
                return 
            end if
            xini=xval
        end do
        write(*,*)'Numero de iteraciones Excedido' , max_iter
        nsolver_svar=xval
end function


real*8 function solver_svar(xini,max_iter)
    !funciones
    real*8 fode_eval
    !variables
    real*8 xini,xval 
    integer max_iter
    integer itern

    real*8 y_k,x1,y_kp
    !Argumenos Globales
    integer degree,solv_type
    real*8 h,tol,initial_val,values(10)
    COMMON /config/ h,tol,initial_val,values,degree,solv_type
    
        y_k=initial_val
        y_kp=y_k
        x1=.1
        do itern=1 , max_iter
            xval=y_k+h*(fode_eval(x1,y_kp))
             if (abs(xval-y_kp)<tol) then
                    solver_svar=xval
                    return
            end if
            y_kp=xval
        end do
        write(*,*)'Numero de iteraciones Excedido' , max_iter
        solver_svar=xval
end function



subroutine nsolver_seed(xini_)
    !funciones
    real*8 nsolver_svar
    
    !variables
    real*8 xini_ 
    integer i
    !Argumenos Globales
    integer degree,solv_type
    real*8 h,tol,initial_val,values(10)
    COMMON /config/ h,tol,initial_val,values,degree,solv_type
    values(1)=xini_
    do i=2,degree
       values(i)=nsolver_svar(xini_,10)
    end do
end subroutine



program ODE

implicit none
integer i;
real*8 initial_guess
real feval
real sol

!Funciones
real*8  nsolver_svar
real*8  solver_svar
integer get_solver_type


    !Argumenos Globales
    integer degree,solv_type
    real*8 h,tol,initial_val,values(10)
    COMMON /config/ h,tol,initial_val,values,degree,solv_type
    
h=.00500000000000000000000
tol=.000000001
initial_val=1.0


solv_type=get_solver_type("EULER-IMPLICIT")
    initial_val=1.0
    do i=1,20
        initial_guess=initial_val;
        sol=nsolver_svar(initial_guess,10)
        write(*,*)h*i, sol
        initial_val=sol
    end do
    initial_val=1.0
    pause
    do i=1,20
        initial_guess=initial_val;
        sol=solver_svar(initial_guess,100)
        write(*,*)h*i, sol
        initial_val=sol
    end do

 solv_type=get_solver_type("TRAPEZOIDAL")
initial_val=1.0
    do i=1,20
        initial_guess=initial_val;
        sol=nsolver_svar(initial_guess,10)
        write(*,*)h*i, sol
        initial_val=sol
    end do
    pause
    
     solv_type=get_solver_type("ADAMS-BASHFORTH")
    
     initial_val=1.0
     degree=2
     call nsolver_seed(initial_val)
     initial_val=1.0
    
    do i=1,20
        initial_guess=initial_val;
        sol=nsolver_svar(initial_guess,10)
        write(*,*)h*i, sol
        initial_val=sol
    end do
    pause
    
    
!crea las cantidades previas
!nsolver_seed(4)



pause
end program ODE



