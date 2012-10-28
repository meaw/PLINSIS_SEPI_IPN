
  
  subroutine matmul4(A,B, C )
implicit none

   real*4, dimension(2,2), INTENT(IN) :: A,B
   real*4, dimension(2,2)  :: C
   integer i,j,k
   c=0;
   do i = 1, 2
    do j = 1, 2
     do k= 1, 2
      C(i, j) =  C(i, j)+A(i,k)*B(k,j)
    end do
  end do
  end do
   end subroutine
   
     subroutine matmul16(A,B, C )
implicit none

   real*16, dimension(2,2), INTENT(IN) :: A,B
   real*16, dimension(2,2)  :: C
   integer i,j,k
   c=0Q0;
   do i = 1, 2
    do j = 1, 2
     do k= 1, 2
      C(i, j) =  C(i, j)+A(i,k)*B(k,j)
    end do
  end do
  end do
   end subroutine
     subroutine matmul8(A,B, C )
implicit none

   real*8, dimension(2,2)  :: A
   real*8, dimension(2,2)  :: B
   real*8, dimension(2,2)  :: C
   integer i,j,k
   c=0;
   do i = 1, 2
    do j = 1, 2
     do k= 1, 2
      C(i, j) =  C(i, j)+A(i,k)*B(k,j)
    end do
  end do
  end do
   end subroutine
  
  
  ! ==================================================================

real*16 function e_tol16(A)
implicit none

   real*16, dimension(2,2), INTENT(IN) :: A
   integer i,j
   real*16 sum_e
   
   sum_e=0Q0;
   do i=1,2
   do j=1,2
       sum_e=sum_e+abs(A(i,j))
   end do
   end do
   e_tol16=sum_e
   end function
  

! ==================================================================

real*8 function e_tol8(A)
implicit none

   real*8, dimension(2,2), INTENT(IN) :: A
   integer i,j
   real*8 sum_e
   
   sum_e=0;
   do i=1,2
   do j=1,2
       sum_e=sum_e+abs(A(i,j))
   end do
   end do
   e_tol8=sum_e
   end function
   
! ==================================================================

real*4 function e_tol4(A)
implicit none

   real*4, dimension(2,2), INTENT(IN) :: A
   integer i,j
   real*4 sum_e
   
   sum_e=0;
   do i=1,2
   do j=1,2
       sum_e=sum_e+abs(A(i,j))
   end do
   end do
   e_tol4=sum_e
   end function
subroutine  diag8(A_mres,n)
    real*8, dimension( 2, 2 ) :: A_mres
    integer n,i
    A_mres=0
    do i=1,n
    A_mres(i,i)=1
    end do
end subroutine
subroutine  diag4(A_lres,n)
    real*4, dimension( 2, 2 ) :: A_lres
    integer n,i
    A_lres=0
    do i=1,n
    A_lres(i,i)=1
    end do
end subroutine

   
   subroutine  diag16(A_hres,n)
    real*16, dimension( 2, 2 ) :: A_hres
    integer n,i
    A_hres=0
    do i=1,n
    A_hres(i,i)=1Q0
    end do
end subroutine
   
      subroutine fp4( pnum,ret)
        real*4  pnum
        character (LEN=130) pstr
        character (LEN=40) ret
        INTEGER BEGIN,I,DLIMIT,j
        DLIMIT=7+1 !P.DEC
        1000 format (f120.80)
        write(pstr,1000)pnum
        pstr = ADJUSTL (pstr)    !RESTRINGE  A LA CANTIDAD DE DIGITOS REALES
        BEGIN=0
        IF (PSTR(1:1) .NE. '-') THEN
        pstr='+'//pstr
        END IF
        DO I=2,130
          IF (PSTR(i:I) .EQ. '.' ) THEN
          DLIMIT=DLIMIT-1;
          END IF
          
            IF (.not.((PSTR(i:I) .EQ. '0') .or.  (PSTR(i:I) .EQ. '.')) ) THEN
                DO J=I+dlimit,130
                if (PSTR(J:J) .ne.  '.') then
                    PSTR(J:J)='0' 
                    end if
                END DO
                exit
            END IF
        END DO
        ret=PSTR(1:40)
    end subroutine
   
   subroutine fp8( pnum,ret)
        real*8  pnum
        character (LEN=130) pstr
        character (LEN=40) ret
        INTEGER BEGIN,I,DLIMIT,J
        DLIMIT=16+1 !P.DEC
        1000 format (f120.80)
        write(pstr,1000)pnum
        pstr = ADJUSTL (pstr)    !RESTRINGE  A LA CANTIDAD DE DIGITOS REALES
        BEGIN=0
        IF (PSTR(1:1) .NE. '-') THEN
        pstr='+'//pstr
        END IF
        DO I=2,130
          IF (PSTR(i:I) .EQ. '.' ) THEN
          DLIMIT=DLIMIT-1;
          END IF
          
            IF (.not.((PSTR(i:I) .EQ. '0') .or.  (PSTR(i:I) .EQ. '.')) ) THEN
                DO J=I+dlimit,130
                    PSTR(J:J)='0'
                END DO
                exit
            END IF
        END DO
        ret=PSTR(1:40)
    end subroutine
   
    subroutine  solver_EXACT8()
        implicit none
        100 format (A12,i4,f55.35)
        101 format (f39.34,f39.34)
        1010 format (A39,x,A39)
        real*8, dimension( 2, 2 ) :: A_mres,p,p1,A_mres_T
        CHARACTER (len=40)fp_print(10)
        A_mres=0
        A_mres(1,1)=exp(-1.0)
        A_mres(2,2)=exp(-17.0)

        P(1,1)=1
        P(1,2)=3
        P(2,1)=2
        P(2,2)=4

        P1(1,1)=-2
        P1(1,2)=1.5
        P1(2,1)=1
        P1(2,2)=-.5

call        Matmul8(A_mres,p1,A_mres_T)
call        Matmul8(p,A_mres_T,A_mres)
            write (*,100) 'Exact (8)'
                call FP8(A_mres(1,1),fp_print(1))
                call FP8(A_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP8(A_mres(2,1),fp_print(1))
                call FP8(A_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
    end subroutine
   
    subroutine fp16( pnum,ret)
        real*16  pnum
        character (LEN=130) pstr
        character (LEN=40) ret
        INTEGER BEGIN,I,DLIMIT,J
        DLIMIT=34+1 !P.DEC
        1000 format (f120.80)
        write(pstr,1000)pnum
        pstr = ADJUSTL (pstr)    !RESTRINGE  A LA CANTIDAD DE DIGITOS REALES
        BEGIN=0
        IF (PSTR(1:1) .NE. '-') THEN
        pstr='+'//pstr
        END IF
        DO I=2,130
            IF (PSTR(i:I) .EQ. '.' ) THEN
                DLIMIT=DLIMIT-1;
            END IF
          
            IF (.not.((PSTR(i:I) .EQ. '0') .or.  (PSTR(i:I) .EQ. '.')) ) THEN
                DO J=I+DLIMIT,130
                    PSTR(J:J)='0'
                END DO
                exit
            END IF
        END DO
        ret=PSTR(1:40)
    end subroutine

    subroutine  solver_EXACT16()
        implicit none
        100 format (A12,i4,f55.35)
        101 format (f39.34,f39.34)
        1010 format (A39,x,A39)
        102 format (f38.23,f38.23)
        real*16, dimension( 2, 2 ) :: A_mres,p,p1,A_mres_T
        CHARACTER (len=40)fp_print(10)
        A_mres=0Q0 
        A_mres(1,1)=Qexp(-1.0Q0 )
        A_mres(2,2)=Qexp(-17.0Q0 )

        P(1,1)=1Q0 
        P(1,2)=3Q0 
        P(2,1)=2Q0 
        P(2,2)=4Q0 

        P1(1,1)=-2Q0 
        P1(1,2)=1.5Q0 
        P1(2,1)=1Q0 
        P1(2,2)=-.5Q0 

        call Matmul16(A_mres,p1,A_mres_T)
        call Matmul16(p,A_mres_T,A_mres)
        write (*,100) 'Exact (16)'
            call FP16(A_mres(1,1),fp_print(1))
            call FP16(A_mres(1,2),fp_print(2))
        WRITE (*,1010)fp_print(1),fp_print(2)
            call FP16(A_mres(2,1),fp_print(1))
            call FP16(A_mres(2,2),fp_print(2))
        WRITE (*,1010)fp_print(1),fp_print(2)
    end subroutine
       
   
   
   
   
   
   
   
   
   
   
   
   subroutine  solver_trad16(A_mres)
    implicit none
    100 format (A12,i4,x,x,A39)
     1010 format (A39,x,A39)
    real*16, dimension( 2, 2 ) :: A_mres,AT_mres,AS_mres,AT2_mres,AT_mres_T
    real*16  fact16,e_tol16
    CHARACTER (len=40)fp_print(10)
    integer i
    call diag16(AT_mres,2)
    AS_mres=AT_mres;
    fact16=1
    do i=1,200
        fact16=fact16*i
        call matmul16(AT_mres,A_mres,AT_mres_T)
        AT_mres=AT_mres_T
        AT2_mres=AT_mres/fact16
        AS_mres=AS_mres+AT2_mres
        if (mod(i,10) .eq. 0) then
            call FP16(e_tol16(AT2_mres),fp_print(1))
            write (*,100) 'Iteracion',i,fp_print(1)
                call FP16(AS_mres(1,1),fp_print(1))
                call FP16(AS_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP16(AS_mres(2,1),fp_print(1))
                call FP16(AS_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
        end if
        if (e_tol16(AT2_mres)<1e-34) then
        call FP16(e_tol16(AT2_mres),fp_print(1))
            write (*,100) 'BIteracion',i,fp_print(1)
                call FP16(AS_mres(1,1),fp_print(1))
                call FP16(AS_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP16(AS_mres(2,1),fp_print(1))
                call FP16(AS_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
    end do
end subroutine
  
subroutine  solver_fpaware16(A_mres)
     implicit none
   100 format (A12,i4,x,x,A39)
    1010 format (A39,x,A39)
    !real*4, dimension( 2, 2 ) :: A_lres,AT_lres,AS_lres,AT2_lres
    real*16, dimension( 2, 2 ) :: A_mres, A_E,A_F,A_T
    real*16  fact16,e_tol16
    integer i
    CHARACTER (len=40)fp_print(10)
  A_E = 0;
 call diag16(A_F,2) 
fact16 = 1;
do i=1,200




   A_E = A_E + A_F;

   
     if ((mod(i-1,10) .eq. 0) .and. i>1) then
        call FP16(e_tol16(A_F),fp_print(1))
            write (*,100) 'Iteracion',i-1,fp_print(1)
                call FP16(A_E(1,1),fp_print(1))
                call FP16(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP16(A_E(2,1),fp_print(1))
                call FP16(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
     end if
        if (e_tol16(A_F)<1e-34) then
           call FP16(e_tol16(A_F),fp_print(1))
            write (*,100) '>BIteracion',i-1,fp_print(1)
                call FP16(A_E(1,1),fp_print(1))
                call FP16(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP16(A_E(2,1),fp_print(1))
                call FP16(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
        
call        matmul16(A_F,A_mres,A_T)
           A_F = A_T/fact16;
   fact16 = fact16+1;
    end do
end subroutine
   
   
   
   
   
   subroutine  solver_fpaware_hump16(A_mres)
     implicit none

    1010 format (I4,x,A39)
    !real*4, dimension( 2, 2 ) :: A_lres,AT_lres,AS_lres,AT2_lres
    real*16, dimension( 2, 2 ) :: A_mres, A_E,A_F,A_T
    real*16  fact16,e_tol16
    integer i
    CHARACTER (len=40)fp_print(10)
  A_E = 0;
 call diag16(A_F,2) 
fact16 = 1;
do i=1,200
        A_E = A_E + A_F;
        call FP16(e_tol16(A_F),fp_print(1))
        write (*,1010) i,fp_print(1)

        call        matmul16(A_F,A_mres,A_T)
        A_F = A_T/fact16;
        fact16 = fact16+1;
  end do
end subroutine
   
   
subroutine  solver_trad8(A_mres)
    implicit none
   100 format (A12,i4,x,x,A39)
     1010 format (A39,x,A39)
    !real*4, dimension( 2, 2 ) :: A_lres,AT_lres,AS_lres,AT2_lres
    real*8, dimension( 2, 2 ) :: A_mres,AT_mres,AS_mres,AT2_mres,AT_mres_T
    real*8  fact8,e_tol8
    CHARACTER (len=40)fp_print(10)
    integer i
    call diag8(AT_mres,2)
    AS_mres=AT_mres;
    fact8=1
    do i=1,100
        fact8=fact8*i
        call matmul8(AT_mres,A_mres,AT_mres_T)
        AT_mres=AT_mres_T
        AT2_mres=AT_mres/fact8
        AS_mres=AS_mres+AT2_mres
        if (mod(i,10) .eq. 0) then
           call FP8(e_tol8(AT2_mres),fp_print(1))
            write (*,100) 'Iteracion',i,fp_print(1)
            
                call FP8(AS_mres(1,1),fp_print(1))
                call FP8(AS_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP8(AS_mres(2,1),fp_print(1))
                call FP8(AS_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
        end if
        if (e_tol8(AT2_mres)<1e-15) then
           call FP8(e_tol8(AT2_mres),fp_print(1))
            write (*,100) '>BIteracion',i,fp_print(1)
                call FP8(AS_mres(1,1),fp_print(1))
                call FP8(AS_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP8(AS_mres(2,1),fp_print(1))
                call FP8(AS_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
    end do
end subroutine
  
subroutine  solver_fpaware8(A_mres)
     implicit none
     100 format (A12,i4,x,x,A39)
    1010 format (A39,x,A39)
    !real*4, dimension( 2, 2 ) :: A_lres,AT_lres,AS_lres,AT2_lres
    real*8, dimension( 2, 2 ) :: A_mres, A_E,A_F,A_T
    real*8  fact8,e_tol8
    integer i
 
     
    CHARACTER (len=40)fp_print(10)
  A_E = 0;
 call diag8(A_F,2) 
fact8 = 1;
do i=1,100
   A_E = A_E + A_F;
   
   
     if ((mod(i-1,10) .eq. 0) .and. i>1) then
     

 
          call FP8(e_tol8(A_F),fp_print(1))
            write (*,100) 'Iteracion',i-1,fp_print(1)
        
                call FP8(A_E(1,1),fp_print(1))
                call FP8(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP8(A_E(2,1),fp_print(1))
                call FP8(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
     end if
        if (e_tol8(A_F)<1e-15) then
           call FP8(e_tol8(A_F),fp_print(1))
            write (*,100) '>BIteracion',i-1,fp_print(1)
                call FP8(A_E(1,1),fp_print(1))
                call FP8(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP8(A_E(2,1),fp_print(1))
                call FP8(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
        
       call matmul8(A_F,A_mres,A_T)
        A_F = A_T/fact8;
   fact8 = fact8+1;
   
    end do
end subroutine
     
     
     integer function  check_matrix(A_hres)
    real*4, dimension( 2, 2 ) :: A_hres
    real*4 a;
    integer n,i,j
    a=huge(a)*huge(a)
 n=0
    do i=1,2
    do j=1,2
   if (( A_hres(i,j) .eq. a) .or. ( A_hres(i,j) .eq. -1*a)) then
   n=1
   end if
    end do
    end do
    check_matrix= n
end function
   
   
   
   subroutine  solver_trad4(A_mres)
    implicit none
   100 format (A12,i4,x,x,A39)
     1010 format (A39,x,A39)
    !real*4, dimension( 2, 2 ) :: A_lres,AT_lres,AS_lres,AT2_lres
    real*4, dimension( 2, 2 ) :: A_mres,AT_mres,AS_mres,AT2_mres,AT_mres_T
    real*4  fact4,e_tol4
    CHARACTER (len=40)fp_print(10)
    integer i,check_matrix
    call diag4(AT_mres,2)
    AS_mres=AT_mres;
    fact4=1
    do i=1,100
        fact4=fact4*i
        call matmul4(AT_mres,A_mres,AT_mres_T)
        if (check_matrix(AT_mres_T) .ne. 0) then
        Write(*,*) 'Error: Un elemento de la matriz es mas grande que el tamaño del tipo flotante'
        exit
        end if
        
        
        AT_mres=AT_mres_T
        AT2_mres=AT_mres/fact4
        AS_mres=AS_mres+AT2_mres
        if (mod(i,10) .eq. 0) then
           call FP4(e_tol4(AT2_mres),fp_print(1))
            write (*,100) 'Iteracion',i,fp_print(1)
            
                call FP4(AS_mres(1,1),fp_print(1))
                call FP4(AS_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP4(AS_mres(2,1),fp_print(1))
                call FP4(AS_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
        end if
        if (e_tol4(AT2_mres)<1e-7) then
           call FP4(e_tol4(AT2_mres),fp_print(1))
            write (*,100) '>BIteracion',i,fp_print(1)
                call FP4(AS_mres(1,1),fp_print(1))
                call FP4(AS_mres(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP4(AS_mres(2,1),fp_print(1))
                call FP4(AS_mres(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
    end do
end subroutine
   
  subroutine  solver_fpaware4(A_Lres)
     implicit none
      100 format (A12,i4,x,x,A39)
    1010 format (A39,x,A39)
    real*4, dimension( 2, 2 ) :: A_Lres, A_E,A_F,A_T
    real*4  fact4,e_tol4
    integer i
    CHARACTER (len=40)fp_print(10)
  A_E = 0;
 call diag4(A_F,2) 
fact4 = 1;
do i=1,100
   A_E = A_E + A_F;
  
     if ((mod(i-1,10) .eq. 0) .and. (i>1)) then
          call FP4(e_tol4(A_F),fp_print(1))
            write (*,100) 'Iteracion',i-1,fp_print(1)
        
                call FP4(A_E(1,1),fp_print(1))
                call FP4(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP4(A_E(2,1),fp_print(1))
                call FP4(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
     end if
        if (e_tol4(A_F)<.0000001) then
                call FP4(e_tol4(A_F),fp_print(1))
            write (*,100) '>BIteracion',i-1,fp_print(1)
            
                call FP4(A_E(1,1),fp_print(1))
                call FP4(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP4(A_E(2,1),fp_print(1))
                call FP4(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
call     matmul4(A_F/fact4,A_lres,A_T)  
         A_F =  A_T;
   fact4 = fact4+1;
   
    end do
end subroutine


subroutine  solver_fpaware4b(A_Lres)
     implicit none
      100 format (A12,i4,x,x,A39)
    1010 format (A39,x,A39)
    real*4, dimension( 2, 2 ) :: A_Lres, A_E,A_F,A_T
    real*4  fact4,e_tol4
    integer i
    CHARACTER (len=40)fp_print(10)
  A_E = 0;
 call diag4(A_F,2) 
fact4 = 1;
do i=1,100
   A_E = A_E + A_F;
  
     if ((mod(i-1,10) .eq. 0) .and. (i>1)) then
          call FP4(e_tol4(A_F),fp_print(1))
            write (*,100) 'Iteracion',i-1,fp_print(1)
        
                call FP4(A_E(1,1),fp_print(1))
                call FP4(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP4(A_E(2,1),fp_print(1))
                call FP4(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
     end if
        if (e_tol4(A_F)<.0000001) then
                call FP4(e_tol4(A_F),fp_print(1))
            write (*,100) '>BIteracion',i-1,fp_print(1)
            
                call FP4(A_E(1,1),fp_print(1))
                call FP4(A_E(1,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
                call FP4(A_E(2,1),fp_print(1))
                call FP4(A_E(2,2),fp_print(2))
            WRITE (*,1010)fp_print(1),fp_print(2)
            exit 
        end if
call     matmul4(A_F,A_lres,A_T)  
         A_F =  A_T/fact4;
   fact4 = fact4+1;
   
    end do
end subroutine
   
   
  program Main
   implicit none


real*4, dimension( 2, 2 ) :: A_lres
real*8, dimension( 2, 2 ) :: A_mres
real*16, dimension( 2, 2 ) :: A_hres


A_hres(1,1)=-49Q0
A_hres(1,2)=24Q0
A_hres(2,1)=-64Q0
A_hres(2,2)=31Q0
A_mres=A_hres
A_lres=A_mres
call  solver_EXACT16()
call  solver_EXACT8()

write(*,*) 'COMMON 16'
call solver_trad16(A_hres)
write(*,*) 'SOLVER FP AWARE 16'
call solver_fpaware16(A_hres)


write(*,*) 'COMMON 8'
call solver_trad8(A_mres)
write(*,*) 'SOLVER FP AWARE 8'
call solver_fpaware8(A_mres)
write(*,*) 'SOLVER FP AWARE 4'
call solver_fpaware4(A_lres)
write(*,*) 'SOLVER FP AWARE 4b'
call solver_fpaware4b(A_lres)
write(*,*) 'COMMON 4'
call solver_trad4(A_lres)


pause
call solver_fpaware_hump16(A_hres)
 pause
  end program
