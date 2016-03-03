*	PARAMETER( nf=110001 )

******** ALL **********************************************************
	REAL target_lat,target_lon
	INTEGER points_num, levels_num, i
	
******** ALLOCATABLE **************************************************
	REAL(8) Levels[ALLOCATABLE](:), PData[ALLOCATABLE] (:,:)
	REAL Pcoor[ALLOCATABLE](:,:)

*	CHARACTER*256 pyt_refl,file_loc
	
***********************************************************************
****** INPUT PARAMETERS ***********************************************
***********************************************************************

* ���� ����������� ������� - ����� 17 ����� �� ������ � NCEP
	  open(1, file='air.dat')
	  
* ��������� ��� ���������� ����� �� ������� ����� ���������������
	  read(1,*) target_lat, target_lat
* ��������� ���-�� ����� ������������ ������� � ���-�� ����� �� ������
	  read(1,*) points_num, levels_num

* �������� ������ �������� [��]
	  ALLOCATE(Levels(1:levels_num),stat=ier)

* �������� ������ ������ ������������ �������
	  ALLOCATE(PData(1:points_num, 1:levels_num))

* �������� ������ ��������� ����� ������������ �������	  
      ALLOCATE(Pcoor(1:points_num, 1:2))

* ��������� �� ����� ������ �������� [��]      
	  read(1,*) Levels(1:levels_num)

* ��������� �� ���������� ������ �������� � ����� ��� ������ �� ������ 	  
	  read(1,*) ((Pcoor(i, 1:2),PData(i, 1:levels_num)), i=1, points_num)
	  close(1)

* ������� ��� ��������� �������� � ���� ��� �������	  
	  open(5,file='air.out.dat')
	  
	  write(5,'(f8.3), 1x') Levels(1:levels_num)
	  write(5,*) ((Pcoor(i, 1:2),PData(i, 1:levels_num)), i=1, points_num)
	  
	  close(5)

* ���� �������� ����� u	- ����� 17 ����� �� ������ � NCEP
      open(2, file='uwnd.dat')	
      
* ��������� ���� ������� ���������� ��� ��� air     
      	
      close(2)
      
* ���� �������� �����	v  - ����� 17 ����� �� ������ � NCEP
      open(3, file='vwnd.dat')	

* ��������� ���� ������� ���������� ��� ��� air
	
      close(3)

* ���� �������� ��������� - ����� ������ ����� �� ������ � NCEP	    
      open(1, file='rhum.dat')	
      
      DEALLOCATE(Levels, PData, Pcoor)
      read(1,*) target_lat, target_lat
	  read(1,*) points_num, levels_num
	  
	  ALLOCATE(Levels(1:levels_num),stat=ier)
	  ALLOCATE(PData(1:points_num, 1:levels_num))
      ALLOCATE(Pcoor(1:points_num, 1:2))
      
	  read(1,*) Levels(1:levels_num)
	  read(1,*) ((Pcoor(i, 1:2),PData(i, 1:levels_num)), i=1, points_num)
	  close(1)
	  
	  open(5,file='meteoint.out.dat')
	  
	  write(5,'(f8.3, 2x)') Levels(1:levels_num)
	  write(5,*) ((Pcoor(i, 1:2),PData(i, 1:levels_num)), i=1, points_num)
	  
	  close(5)
      	
      close(1)
    

	
		DEALLOCATE(Levels, PData, Pcoor)
***********************************************************************
****** PROCESSING ****************************************************
***********************************************************************

      END


