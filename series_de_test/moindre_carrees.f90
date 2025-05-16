module moindre_carree
  implicit none



contains

  subroutine lire_y(n_cas,N,y,te)
    implicit none
    real*8,dimension(:),intent(inout) :: y , te
    integer,intent(in) :: n_cas ,  N
    integer :: t 
    real*8 :: a , b , c , d , e , f , g , h , i , j , k
    select case(n_cas)
    case(1)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , y(t)  
          t = t + 1
       end do
    case(2)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a , y(t) 
          t = t + 1
       end do
    case(3)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a , b , y(t) 
          t = t + 1
       end do
    case(4)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a , b , c , y(t) 
          t = t + 1
       end do
    case(5)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,  b , c , d , y(t) 
          t = t + 1
       end do
    case(6)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,  b , c , d , e , y(t) 
          t = t + 1
       end do
    case(7)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a , b , c , d , e , f , y(t) 
          t = t + 1
       end do
    case(8)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,  b , c , d , e , f , g , y(t) 
          t = t + 1
       end do
    case(9)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,   b , c , d , e , f , g , h , y(t) 
          t = t + 1
       end do
    case(10)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,  b , c , d , e , f , g , h , i , y(t) 
          t = t + 1
       end do
    case(11)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,  b , c , d , e , f , g , h , i , j , y(t)
          t = t + 1
       end do
    case(12)
       t = 1
       do while ( t <= N )
          read(15,*) te(t) , a ,  b , c , d , e , f , g , h , i , j , k , y(t)
          t = t + 1
       end do

    end select

  end subroutine lire_y







  subroutine resolution(A,L,D,b,z,x,y1,n_case)
    implicit none
    real*8,dimension(:,:),intent(inout) :: A , L
    real*8,dimension(:),intent(inout) :: b , D , z , y1,x
    integer :: n_case
    x = 0;
    y1 = 0;
    z = 0;
    D = 0;
    L = 0;
    select case(n_case)
    case(0)
       call factorisation_cholesky(A,L,D)
       call descente(L,b,z)
       call diagonal(D , z , y1)
       call remontee( L , y1 , x)
    case(1)
       call gradient(A,x,b)
    end select
  end subroutine resolution

  subroutine lit_parametre(fichier1,fichier2,deg,sigma,M,T1,T2)
    implicit none
    character(len = 50),intent(inout) ::fichier1 , fichier2
    integer,intent(inout) :: deg , M
    real*8,intent(inout) :: T1 , T2 , sigma
    print*,'Entrez le fichier  d entree'
    read*,fichier1
    print*,'Entrez le fichier de sortie'
    read*,fichier2
    print*,'Entrez le degre du polynome'
    read*,deg
    print*,'Entrez le sigma la variance choisir un grand nombre'                            
    read*,sigma
    print*,'Entrez T1 temps de départ  entre [1:60000]' 
    read*,T1
    print*,'T2 temps d arrivee entre [1:60000]'
    read*,T2
    print*,'Entrez l entier M pour avoir le pas de temps'
    read*,M
  end subroutine lit_parametre


  subroutine allouer(A,Pi,L,y1,D,b,z,x,N,deg)
    implicit none
    real*8,dimension(:,:),intent(inout),allocatable :: A,Pi,L
    real*8,dimension(:),intent(inout),allocatable :: y1,D,b,z,x
    integer :: N , deg
    allocate(Pi(N,deg))
    allocate(A(deg,deg))
    allocate(L(deg,deg))
    allocate(y1(deg))
    allocate(D(deg))
    allocate(b(deg))
    allocate(z(deg))
    allocate(x(deg))
  end subroutine allouer


  subroutine desallouer(A,Pi,L,y1,D,b,z,x)
    implicit none
    real*8,dimension(:,:),intent(inout),allocatable :: A,Pi,L
    real*8,dimension(:),intent(inout),allocatable :: y1,D,b,z,x
    deallocate(Pi)
    deallocate(A)
    deallocate(L)
    deallocate(y1)
    deallocate(D)
    deallocate(b)
    deallocate(z)
    deallocate(x)
  end subroutine desallouer


  subroutine systeme_lineaire(A,Pi,xi,y,te,b,sigma)
    implicit none
    real*8,dimension(:,:),intent(inout) :: A , Pi 
    real*8,dimension(:),intent(inout) :: b , y , te
    real*8 :: xi , sigma
    integer ::i , j , N , d , k
    N = size(Pi(:,1))
    d = size(Pi(1,:))
    Pi = 0
    do i = 1 , N
       do j = 1 , d
          Pi(i,j) = (( te(i) - xi )**(j-1))
       end do
    end do
    A = 0
    do i = 1 , d
       do k = 1 , N
          do j = 1 , d
             A(i,j) = A(i,j) + Pi(k,i)*exp(((-1)/2)*((te(k) - xi)**2)/(sigma**2))*Pi(k,j)
          end do
       end do
    end do
    b = 0
    do i = 1 , d
       do k = 1 , N
          b(i) = b(i) + Pi(k,i)*y(k)
       end do
    end do
  end subroutine systeme_lineaire



  !------------- Factorisation de Cholesky -----------------
  subroutine factorisation_cholesky(A,L,D)
    implicit none
    real*8,dimension(:,:),intent(in)::A
    real*8,dimension(:,:),intent(out)::L
    real*8,dimension(:),intent(out)::D
    integer::n,i,j,k

    n=size(D)

    L=0
    do i=1,n
       do j=1,i-1
          !Calcul de a_ij
          L(i,j)=A(i,j)
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*D(k)*L(j,k)
          end do
          L(i,j)=L(i,j)/D(j)
       end do
       L(i,i)=1
       !Calcul de d_ii
       D(i)=A(i,i)
       do k=1,i-1
          D(i)=D(i)-L(i,k)*L(i,k)*D(k)
       end do
    end do
  end subroutine factorisation_cholesky
  !------------- Descente -----------------
  subroutine descente(A,b,x)
    implicit none
    real*8,dimension(:,:),intent(in)::A
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out)::x
    real*8::somme
    integer::n,i,j

    n=size(b)

    !Solution pour l'element 1: trivial
    x(1)=b(1)/A(1,1)
    !On part du bas vers le haut
    do i=2,n
       somme=0
       do j=1,i-1
          somme=somme+A(i,j)*x(j)
       end do
       x(i)=(b(i)-somme)/A(i,i)
    end do
  end subroutine descente
  !------------- Diagonal -----------------
  subroutine diagonal(D,b,x)
    implicit none
    real*8,dimension(:),intent(in)::D,b
    real*8,dimension(:),intent(out)::x
    real*8::somme
    integer::n,i

    n=size(b)

    do i=1,n
       x(i)=b(i)/D(i)
    end do

  end subroutine diagonal

  !------------- Remontee_t -----------------
  subroutine remontee(A,b,x)
    implicit none
    real*8,dimension(:,:),intent(in)::A
    real*8,dimension(:),intent(in)::b
    real*8,dimension(:),intent(out)::x
    real*8::somme
    integer::n,i,j

    n=size(b)

    !Solution pour l'element n: trivial
    x(n)=b(n)/A(n,n)
    !On part du bas vers le haut
    do i=n-1,1,-1
       somme=0
       do j=i+1,n
          somme=somme+A(j,i)*x(j)
       end do
       x(i)=(b(i)-somme)/A(i,i)
    end do

  end subroutine remontee


  subroutine gradient(A,x,b)
    implicit none
    real*8,dimension(:,:),intent(in) :: A
    real*8,dimension(:),intent(in) :: b
    real*8,dimension(:),intent(inout) :: x
    real*8,dimension(size(x)) :: x1 , tmp , r
    real*8 :: somme, sum   , critere_arret , eps , alpha 
    integer :: i , j , n , k = 1
    eps = 1d-6

    n = size(x)
    ! INITIALISATION DE r = b - Ax0

    r = 0
    do i = 1 , n
       do j = 1 , n
          r(i) = r(i) + A(i,j)*x(j)
       end do
    end do
    r = b - r
    ! INITIALISATION DE ALPHA alpha = (rk,rk)/(Ark,rk)
    somme = 0.0

    do i = 1 , n
       somme = somme + r(i)**2
    end do
    tmp = 0.0

    do i = 1 , n
       do j = 1 , n
          tmp(i) = tmp(i) + A(i,j)*r(j)
       end do
    end do
    sum = 0.0

    do i = 1 , n
       sum = sum + tmp(i)*r(i)
    end do
    alpha = somme/sum
    
    ! xk+1 = xk + alphak*rk

    x = x + alpha*r

    ! rk+1 = rk - alphak*Ark

    tmp = 0.0

    do i = 1 , n
       do j = 1 , n
          tmp(i) = tmp(i) + A(i,j)*r(j)
       end do
    end do

    r = r - alpha*tmp

    ! Calcul du premier critere_arret ||Ax - b||


    tmp = 0.0

    do i = 1 , n
       do j = 1 , n
          tmp(i) = tmp(i) + A(i,j)*x(j)
       end do
    end do


    tmp = tmp - b


    somme = 0.0

    do i = 1 , n
       somme = somme + tmp(i)**2
    end do

    critere_arret  = sqrt(somme)





    do while ( critere_arret > eps )

       ! CALCUL DE ALPHA alpha = (rk,rk)/(Ark,rk)


       somme = 0.0

       do i = 1 , n
          somme = somme + r(i)**2
       end do
       tmp = 0.0

       do i = 1 , n
          do j = 1 , n
             tmp(i) = tmp(i) + A(i,j)*r(j)
          end do
       end do

       sum = 0.0

       do i = 1 , n
          sum = sum + tmp(i)*r(i)
       end do
       alpha = somme/sum
       
       ! xk+1 = xk + alphak*rk

       x = x + alpha*r

       ! rk+1 = rk - alphak*Ark

       tmp = 0.0
       

       do i = 1 , n
          do j = 1 , n
             tmp(i) = tmp(i) + A(i,j)*r(j)
          end do
       end do

       r = r - alpha*tmp

       ! CALCUL DU critere_arret ||Ax - b||

       tmp = 0.0

       do i = 1 , n
          do j = 1 , n
             tmp(i) = tmp(i) + A(i,j)*x(j)
          end do
       end do


       tmp = tmp - b

       somme = 0.0

       do i = 1 , n
          somme = somme + tmp(i)**2
       end do

       critere_arret  = sqrt(somme)

       k = k + 1
    end do
  end subroutine gradient


end module moindre_carree



