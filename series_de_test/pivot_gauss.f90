program ex4
  implicit none

  !Declarations
  real*8,dimension(:,:),allocatable::A,A2 !la matrice A
  real*8,dimension(:),allocatable::x,b,b2 !les vecteurs x et b
  integer::n
  integer::i,j

  !Instructions  
  n=5

  !On alloue la mémoire nécessaire à représenter la matrice et les vecteurs
  allocate(A(n,n))
  allocate(A2(n,n))
  allocate(x(n))
  allocate(b(n))
  allocate(b2(n))

  A=0
  A(1,1)=1
  A(1,3)=2
  A(1,5)=10
  A(2,2)=1.D-16
  A(3,3)=1
  A(2,4)=4
  A(3,1)=2
  A(3,2)=5
  A(3,4)=4
  A(4,2)=4
  A(4,3)=4
  A(5,4)=10
  A(5,5)=2

  A2=A
  
  !On initialise le vecteur b aux valeurs fournies dans l'énoncé
  do i=1,n
     b(i)=i
  end do

  print*,'XXXXXX Valeurs initiales XXXXXX'
  print*,'----A----'
  do i=1,n
     print*,A(i,:)
  end do
  print*,'----b----'
  print*,b
  
  print*,'XXXXXX Appel de la fonction XXXXXX'
  call pivot_gauss(A,b)
  
  print*,'XXXXXX Matrice triangulaire XXXXXX'
  print*,'----A----'
  do i=1,n
     print*,A(i,:)
  end do
  print*,'----b----'
  print*,b
  
  call remontee(A,b,x)
  
  print*,'XXXXXX Solution XXXXXX'
  print*,'----x----'
  print*,x

  b2=0
  do i=1,n
     do j=1,n
        b2(i)=b2(i)+A2(i,j)*x(j)
     end do
  end do

  print*,'----b2=Ax----'
  print*,b2
     
  !On libère la mémoire occupée par la matrice et les vecteurs
  deallocate(A)
  deallocate(x)
  deallocate(b)
  deallocate(b2)
       
contains
  !------------- Pivot de Gauss -----------------
  subroutine pivot_gauss(A,b)
    implicit none
    real*8,dimension(:,:),intent(inout)::A
    real*8,dimension(:),intent(inout)::b
    real*8,dimension(:),allocatable::tmp_vec
    real*8::tmp_scal
    integer::n,k,i
    real*8::a_ik_max
    integer::idc_a_ik_max
    real*8::coef

    n=size(b)

    allocate(tmp_vec(n))
    
    do k=1,n-1
       ! Recherche du plus grand pivot en valeur absolue
       a_ik_max=A(k,k)
       idc_a_ik_max=k
       do i=k+1,n
          if (abs(A(i,k))>abs(a_ik_max)) then
             a_ik_max=A(i,k)
             idc_a_ik_max=i
          end if
       end do
       !On inverse les lignes i et k
       !En pratique, il n'est pas nécessaire de faire la permutaion
       !mais on peut stocker les changements d'indices (non fait ici)
       if (idc_a_ik_max /= k) then
          tmp_vec=A(idc_a_ik_max,:)
          A(idc_a_ik_max,:)=A(k,:)
          A(k,:)=tmp_vec
          tmp_scal=b(idc_a_ik_max)
          b(idc_a_ik_max)=b(k)
          b(k)=tmp_scal
          print*,'Echange des lignes',k,'et',idc_a_ik_max
          print*,'A='
          do i=1,n
             print*,A(i,:)
          end do
          print*,'b='
          print*,b
       end if
       print*,'Ligne du pivot:',k,'Pivot:',A(k,k)
       !On elimine A(i,k) dans les equations k+1 a n
       do i=k+1,n
          coef=A(i,k)/A(k,k)
          do j=k,n
             A(i,j)=A(i,j)-A(k,j)*coef
          end do
          b(i)=b(i)-b(k)*coef
       end do
       print*,'A='
       do i=1,n
          print*,A(i,:)
       end do
       print*,'b='
       print*,b
    end do

    deallocate(tmp_vec)

  end subroutine pivot_gauss

  !------------- Remontee -----------------
  subroutine remontee(A,b,x)
    implicit none
    real*8,dimension(:,:),intent(in)::A
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out)::x
    real*8::somme
    integer::n,i

    n=size(b)
    
    !Solution pour l'element n: trivial
    x(n)=b(n)/A(n,n)
    !On part du bas vers le haut
    do i=n-1,1,-1
       somme=0
       do j=i+1,n
          somme=somme+A(i,j)*x(j)
       end do
       x(i)=(b(i)-somme)/A(i,i)
    end do

  end subroutine remontee

end program ex4
