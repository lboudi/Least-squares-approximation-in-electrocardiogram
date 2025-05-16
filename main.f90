program moindre_carre
  use moindre_carree
  implicit none
  real*8,dimension(:,:),allocatable :: A , Pi , L 
  real*8,dimension(:),allocatable :: te , D , b ,  z , y , y1 , x
  real*8 :: T1 , T2 , h , xi , sigma , somme , time1 , time2
  integer :: t , N , M , k , deg , i , n_case , n_cas
  character(len = 50) ::fichier1 , fichier2 , temp

  call  lit_parametre(fichier1,fichier2,deg,sigma,M,T1,T2)

  open(unit = 15 , file = fichier1 , status = 'old')

  open(unit = 16 , file = fichier2 , status = 'new')


  t = 1

  do  while( t <= 4 )
     read(15,*) temp
     t = t + 1
  end do

  N = T2 - T1 + 1
  allocate(te(N))
  allocate(y(N))
  t = 1
  do  while (t < T1)
     read(15,*) h
     t = t + 1
  end do
  print*,'Entrez un entier pour lire les donnees en y de 1:12'
  read*,n_cas
  call   lire_y(n_cas,N,y,te)
  h = ( T2 - T1)/M
  i = 0
  xi = T1 + i*h
  print*,'Entrez entier n_case 0 pour cholesky et 1 pour gradient'
  read*,n_case
  call cpu_time(time1)
  do while( xi <= T2 )
     print*,'On est a l iteration',i
     i = i + 1    

     call allouer(A,Pi,L,y1,D,b,z,x,N,deg)
     call systeme_lineaire(A,Pi,xi,y,te,b,sigma)


     call  resolution(A,L,D,b,z,x,y1,n_case)

     write(16,*) xi ,x(1),x(2),x(3)
     call desallouer(A,Pi,L,y1,D,b,z,x)
     xi = T1 + i*h
  end do
  call cpu_time(time2)
  time1 = time2 - time1;
  print*,'Voici le temps pris pour le calcul'
  print*,time1

  deallocate(y)
  deallocate(te)
  close(15)
  close(16)
end program moindre_carre

