program toto
  use moindre_carree
  implicit none
  real*8,dimension(3,3) :: A
  real*8,dimension(3) :: x , y , z
  integer :: i
  real*8 :: tmp
  A = 0
  A(1,1) = 1
  A(2,2) = 1
  A(3,3) = 1
  print*,'On va entrer le vecteur x'
 1 do i = 1 , 3
     print*,'element',i,'de x'
     read*,x(i)
  end do
  print*,'On va entrer le vecteur y'
  do i = 1 , 3
     print*,'element',i,'de y'
     read*,y(i)
  end do
end program toto
