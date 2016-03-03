

!******* ALL **********************************************************
  REAL target_lat,target_lon
  INTEGER points_num, levels_num, moments_num,  i
  
!******* ALLOCATABLE **************************************************

!  REAL(8) Levels[ALLOCATABLE](:), PData(:,:,:)
  real(8),allocatable,dimension(:) :: Levels
  real(8),allocatable,dimension(:,:,:) :: PData
  real(8),allocatable,dimension(:,:) :: Pcoor


  
!**********************************************************************
!***** INPUT PARAMETERS ***********************************************
!**********************************************************************

! файл температуры воздуха - всего 17 точек по высоте в NCEP
    open(1, file='air.dat')
    
! считываем гео координаты точки на которую нужно интерполировать
    read(1,*) target_lat, target_lat
! считываем кол-во точек покрывающего участка, кол-во точек по высоте
! и количество моментов времени
    read(1,*) points_num, levels_num, moments_num

! выдел€ем массив давлений [мб]
    ALLOCATE(Levels(1:levels_num),stat=ier)

! выдел€ем массив данных покрывающего участка
    ALLOCATE(PData(1:points_num, 1:moments_num, 1:levels_num))

! выдел€ем массив координат точек покрывающего участка    
      ALLOCATE(Pcoor(1:points_num, 1:2))

! считываем из файла массив давлений [мб]      
    read(1,*) Levels(1:levels_num)

! считываем из координаты каждой геоточки и затем все данные по времени и по высоте     
    read(1,*) (Pcoor(i, 1:2), ((PData(i, j, k), k=1,levels_num,1), j=1, moments_num,1), i=1, points_num, 1)
    close(1)

! выводим  считанные значени€ давлени€ и темп по первой точкe air в 
! файл дл€ примера    
    open(5,file='meteoint.out.dat')
    
    write(5,'(f8.3), 1x') Levels(1:levels_num)
  write(5,*) (Pcoor(i, 1:2),ACHAR(13),ACHAR(10), ((PData(i, j, k), k=1,levels_num,1),ACHAR(13),ACHAR(10),ACHAR(13),ACHAR(10), j=1, moments_num,1),ACHAR(13),ACHAR(10), i=1, 1, 1)
  close(5)

! файл скорости ветра u - всего 17 точек по высоте в NCEP
      open(2, file='uwnd.dat')  
      
! повтор€ем теже команды считывани€ как дл€ air     
        
      close(2)
      
! файл скорости ветра v  - всего 17 точек по высоте в NCEP
      open(3, file='vwnd.dat')  

! повтор€ем теже команды считывани€ как дл€ air
  
      close(3)

! файл влажности - всего восемь точек по высоте в NCEP      
!      open(1, file='rhum.dat') 
!   close(1)
    

  
    DEALLOCATE(Levels, PData, Pcoor)
!**********************************************************************
!***** PROCESSING ****************************************************
!**********************************************************************

      END


