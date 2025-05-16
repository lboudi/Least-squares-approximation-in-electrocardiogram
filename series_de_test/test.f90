program test
  implicit none
  integer :: i , N
  real*8 :: x
  open(unit = 15 , file = 'donnees_test.txt', status = 'new')
  print*,'Entrez entier N nombre de points'
  read*,N
  do i = -N , N , 1
     x = i**3
     write(15,*) i ,x
  end do
  close(15)
end program test

     
