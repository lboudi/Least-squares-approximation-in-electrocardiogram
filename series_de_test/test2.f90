program resoudre
  use moindre_carree
  implicit none
  real*8,dimension(10,10) :: A
  real*8,dimension(10) :: b , x , tmp 
  integer :: i , j
!!$  open(unit = 17 , file = 'tete.txt',status = 'old')
!!$  open(unit = 18 , file = 'titi.txt' , status = 'old')
!!$  do i = 1 , 10
!!$     read(18,*) b(i)
!!$  end do
!!$  do i = 1 , 10
!!$     do j = 1 , 10
!!$        read(17,*) A(i,j)
!!$     end do
!!$  end do
  do i = 1 , 10
     do j = 1 , 10
        A(i,j) = 0.5
        A(j,i) = A(i,j)
     end do
  end do
  do i = 1 , 10
     A(i,i) = 1.0
  end do
  print*,'On va rentrer b'
  do i = 1 , 10
     print*,'element',i,'de b'
     read*,b(i)
  end do

  print*,'Voici A'
  do i = 1 , 10
     print*,A(i,:)
  end do
  print*,'Voici b'
  print*,b

  print*,'Maintenant on va appliquer le gradient'
  
  call gradient(A,x,b)
  print*,'Voici x'
  print*,x
  print*,'Voici le produit Ax'
  tmp = 0.0
  do i = 1 , 10
     do j = 1  , 10
        tmp(i) = tmp(i) + A(i,j)*x(j)
     end do
  end do
  print*,tmp
  print*,'Voici b'
  print*,b
     
end program resoudre
