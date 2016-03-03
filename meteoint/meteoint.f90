

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

! ���� ����������� ������� - ����� 17 ����� �� ������ � NCEP
    open(1, file='air.dat')
    
! ��������� ��� ���������� ����� �� ������� ����� ���������������
    read(1,*) target_lat, target_lat
! ��������� ���-�� ����� ������������ �������, ���-�� ����� �� ������
! � ���������� �������� �������
    read(1,*) points_num, levels_num, moments_num

! �������� ������ �������� [��]
    ALLOCATE(Levels(1:levels_num),stat=ier)

! �������� ������ ������ ������������ �������
    ALLOCATE(PData(1:points_num, 1:moments_num, 1:levels_num))

! �������� ������ ��������� ����� ������������ �������    
      ALLOCATE(Pcoor(1:points_num, 1:2))

! ��������� �� ����� ������ �������� [��]      
    read(1,*) Levels(1:levels_num)

! ��������� �� ���������� ������ �������� � ����� ��� ������ �� ������� � �� ������     
    read(1,*) (Pcoor(i, 1:2), ((PData(i, j, k), k=1,levels_num,1), j=1, moments_num,1), i=1, points_num, 1)
    close(1)

! �������  ��������� �������� �������� � ���� �� ������ ����e air � 
! ���� ��� �������    
    open(5,file='meteoint.out.dat')
    
    write(5,'(f8.3), 1x') Levels(1:levels_num)
  write(5,*) (Pcoor(i, 1:2),ACHAR(13),ACHAR(10), ((PData(i, j, k), k=1,levels_num,1),ACHAR(13),ACHAR(10),ACHAR(13),ACHAR(10), j=1, moments_num,1),ACHAR(13),ACHAR(10), i=1, 1, 1)
  close(5)

! ���� �������� ����� u - ����� 17 ����� �� ������ � NCEP
      open(2, file='uwnd.dat')  
      
! ��������� ���� ������� ���������� ��� ��� air     
        
      close(2)
      
! ���� �������� ����� v  - ����� 17 ����� �� ������ � NCEP
      open(3, file='vwnd.dat')  

! ��������� ���� ������� ���������� ��� ��� air
  
      close(3)

! ���� ��������� - ����� ������ ����� �� ������ � NCEP      
!      open(1, file='rhum.dat') 
!   close(1)
    

  
    DEALLOCATE(Levels, PData, Pcoor)
!**********************************************************************
!***** PROCESSING ****************************************************
!**********************************************************************

      END


