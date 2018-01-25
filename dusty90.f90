program dusty 

! The converted dusty deck to Fortran 90
! Prof. A.J. Pounds
! Departments of Chemistry and Computer Science
! Mercer University
! Spring 2008
!
! Modified to use the timing libraries in Spring 2016

integer, parameter :: MAXDIM = CPPFLAGS 

integer, dimension(MAXDIM) :: IA
integer :: N
real (kind=8), dimension(MAXDIM) :: AV, BV, CV
real (kind=8), dimension(MAXDIM,MAXDIM) :: OP, ID, AM, BM, CM, DM 
real (kind=8) :: check, BOT, TOP, HOLDA, HOLDB, TRACE3
real (kind=4) :: start, finish
real (kind=8) :: seed

real (kind=8) :: wall, cpu
double precision trig
external trig


interface
    real (kind=8) function conrand(seed)
    real (kind=8)  :: seed
    end function conrand
end interface

interface
    real(kind=8) function cputime()
    end function cputime 
end interface

interface
    real(kind=8) function walltime()
    end function walltime 
end interface

N = MAXDIM

!call cpu_time(start) 

wall = walltime()
cpu  = cputime()

seed = 1.0D0

!     Fill arrays

! Loop 10 Series -- Filling Arrays

do i = 1, N
   AV(i) = bessel_jn(0,dble(conrand(seed) * &
               (-1)**(mod(int(10*conrand(seed)),N))))
enddo

do i = 1, N
   BV(i) = bessel_jn(1,dble(conrand(seed) * & 
               (-1)**(mod(int(10*conrand(seed)),N))))
enddo

check = 0.0
do i = 1, N
  ival = N
  check = check + AV(i) * BV(i)
  call idcheck(ival,check,AV,BV,ID)
enddo  

! Compute |AV><BV|

do  i = 1, N
  do  j = 1, N
    call idcheck(N,check,AV,BV,ID)
    if ( check .gt. 0.5 ) then
       OP(i,j) = AV(i) * BV(j) / BV(i)
    else
       OP(i,j) = AV(j) * BV(i) / BV(j)
    endif
  enddo 
  IA(i) = i
enddo


do i = 1, N 
  do  j = 0, i, 8  
       IA(I) = mod(mod(i+j,N),N)+1 
  enddo 
enddo

do i = 1, N
enddo

! Loop 20 

! Added to check values of arrays
do i = 1, N
enddo

do  i = 1, N
   call idcheck(N,check,AV,BV,ID)
   CV(IA(I)) = (AV(IA(I)) + BV(IA(I))) / check
enddo



! Loop 30 

do i = 2, N
   call idcheck(N,check,AV,BV,ID)
   AV(i) = AV(i-1) * BV(i) + CV(i)
enddo



! Loop 40 

do i = 1, N
   call idcheck(N,check,AV,BV,ID)
   do j = 1, N
      if ( check .gt. 0.5 ) then
         BOT = OP(i,j) 
         TOP = AV(j) * BV(j)
         HOLDA = AV(j)
         AV(j) = BV(j) + CV(j) / (TOP-BOT) * ID(i,i)
         BV(j) = HOLDA + CV(j) / (TOP-BOT) * ID(j,j)
        ! stop
         AM(i,j) = AV(j) * trig(IA(i),IA(j)) 
         BM(i,j) = BV(j) * trig(IA(j),IA(i)) 
      else
         BOT = OP(i,j) 
         TOP = AV(j) * BV(j)
         HOLDA = AV(j)
         AV(j) = BV(j) - CV(j) / (TOP-BOT) * ID(j,j)
         BV(j) = HOLDA - CV(j) / (TOP-BOT) * ID(i,i)
         AM(i,j) = AV(j) / trig(IA(i),IA(j))
         BM(i,j) = BV(j) / trig(IA(j),IA(i)) 
      endif
   enddo 
enddo

! Loop 50

do  i = 1, N
   do  j = 1, N
      CM(i,j) = 0.0
      do k = 1, N
         if ( i .lt. j ) then
            CM(i,j) = CM(i,j) - AM(i,k) * BM(k,j) / check 
         else
            CM(i,j) = CM(i,j) + AM(i,k) * BM(k,j) / check 
         endif
      enddo 
  enddo
enddo


! Loop 60

do i = 1, N
   do j = 1, N
      sum = 0.0
      do k = 1, N
         sum = sum + CM(i,k) * AM (j,k)
      enddo 
      DM(i,j) = sum
   enddo 
enddo

do i = 1, N
  do  j = 1, N
     CM(i,j) = DM(i,j)
  enddo 
enddo

! Loop 70

do i = 1, N
   do j = 1, N
      sum = 0.0
      do k = 1, N
         sum = sum - CM(i,k) * BM (j,k)
      enddo 
      DM(i,j) = sum
   enddo 
enddo


HOLDA = abs(AM(1,1))
HOLDB = abs(BM(1,1)) 
do i = 1, N
  do j = 1, N
    HOLDA = max(HOLDA,abs(AM(i,j))) 
    HOLDB = max(HOLDB,abs(BM(i,j))) 
  enddo 
enddo
         
TRACE3 = 0.0
       
! Loop 80

do i = 1, N
  TRACE3 = TRACE3 + (AM(IA(i),IA(i)) + BM(IA(i),IA(i)) & 
                  - DM(IA(i),IA(i))) / (HOLDA * HOLDB)
enddo

!call cpu_time(finish) 

cpu = cputime() - cpu
wall = walltime() - wall


#ifndef DO_TIMING 
  print *, 'Final trace = ', trace3, ' and IDCHECK ', check
  print *, '-- RUNTIME -> ', cpu, ' seconds'
#else
  print "(I5,5X,F9.4)", MAXDIM, cpu
#endif 

end
        
         
real (kind=8) function trig (i,j)
real (kind=8) ::  x, y, z
pi = acos(-1.0)
x = dble(i) - dble(j)
y = dble(i) + dble(j) 
z = exp ( sin(sqrt(x**2+y**2)*pi  ) )  
trig = x + y + log10(abs(1+z+(x*y*z)))/ (abs(x)+abs(y))
return
end 

subroutine idcheck(N,check,AV,BV,ID)

real (kind=8) :: AV(*), BV(*), ID(N,*)
real (kind=8) :: l2
real (kind=8) :: check, check2
real (kind=8) :: a, b, c, d 

real (kind=8) :: tempa, tempb, pi;

do i = 1, N  
  do j = 1, N
    if ( i .eq. j ) then 
       if (( AV(i) .lt. 0 ) .and. ( BV(j) .lt. 0 )) then
         ID(i,j) = 1.0
       elseif (( AV(i) .lt. 0 ) .and. ( BV(j) .gt. 0 )) then
         ID(i,j) = -1.0
       elseif (( AV(i) .gt. 0 ) .and. ( BV(j) .lt. 0 )) then
         ID(i,j) = -1.0
       else
         ID(i,j) = 1.0
       endif
    elseif ( i .ne. j ) then 
       ID(i,j) =  cos(check+2.0*i*acos(-1.0)/N)+  &
                 2.0*sin(check+ 2.0*j*acos(-1.0)/N)
    endif
  enddo 
enddo

l2 = 0.0
do i = 1, N
  l2 = l2 + AV(i)**2
enddo


l2 = sqrt(l2)
do i = 1, N
  AV(i) = AV(i) / l2
enddo


l2 = 0.0
do i = 1, N
  l2 = l2 + BV(i)**2
enddo


l2 = sqrt(l2)
do i = 1, N
  BV(i) = BV(i) / l2
enddo

      a = 0.0D0
      b = 0.0D0
      c = 0.0D0
      d = 0.0D0
      do 70 i = 1, N
        do 80 j = 1, N
           do 90 k = 1, N
               goto ( 200, 300, 400, 500 ) int(mod(i+j+k,4)+1) 
200            a  = a +  AV(i) * BV(j) * ID(j,k) 
               check = check + a
               goto 100
300            b  = b +  AV(j) * BV(i) * ID(k,j) 
               check = check - b 
               goto 100
400            c  = c -  AV(i) * BV(j) * ID(k,j) 
               check = sqrt(b**2 + c**2)
               goto 100
500            d  = d -  AV(j) * BV(i) * ID(j,k) 
               check2 = a + b + c + d
100            continue
90         continue
80      continue
70    continue


      check = min(abs(check2),abs(check))/max(abs(check2),abs(check))           


      return
      end


real (kind=8) function conrand(seed)
!
! Function to generate a sequence of random numbers.
! Adapted from the  "Minimal Standard Method, real version 1 in Pascal"
! Park, S, Miller, K. "Random Number Generators: Good Ones are
! Hard to Find".  Communications of the ACM. vol 31, number 10,
! October 1988. pp. 1192-1201.
!
! Fortran 2003 Version tested on 64 Bit Linux, gfortran compiler
! Andrew J. Pounds, Ph.D.
! Departments of Chemistry and Computer Science
! Mercer University
! Fall 2011
!
real (kind=8) :: seed
real (kind=8) :: a, m
real (kind=8) :: temp
a = 16807.0D0
m = 2147483647.0D0
temp = a*seed 
seed = temp - m * int(temp/m)
conrand = seed / m
return
end function conrand
