c   Program to spatial interpolation, using Gandin's optimum algorithm  
c  Subroutines:
c <regress> -  Procedure of the finding to exponential regression 
c              for dependency correlation from distances
c <Gandin> - Gandin optimal spatial interpolation procedure
c <inter_japan_time_step> - a forecast of some meteoparameter 
c      INPUT PARAMETERS:
c    Lat_0,Lon_0 - latitude and longitude of a target point (from initial data file)
c    N_point - a number base points (from initial data file)
c    N_time  - a number temporary counting out (from initial data file)
c    N_h     - a number pressure levels (from initial data file)
c    Delta   - a target time (from initial data file) 
c    kkk     - amount affecting station 
c      OUTPUT DATA:
c    VALUES OF METEOPARAMETER IN THE TARGET POINT AND TARGET TIME (in file RESULT.DAT)
c      ADDITION INFORMATION:
c      standard deviation  relative nearest point(in file SQ.DAT)
c      result od interpolation for average values at basic points (in file AVER.DAT)
c      time series at target point (in file HT.DAT)  
c      d(TEMP)/dt  at basic points 

      program main_japan
	use MSFLIB

      Real f0,f_in,f_fin,L0,L_in,L_fin
	Real srv(8),Lat_0,Lon_0
      character*12 f_ini(4)

	data f_ini/'air.dat','uwnd.dat','vwnd.dat','rhum.dat'/ 
	data file_xy/'xy.dat'/
      Character*20 Regress_Err,Inter_Error,Time_Error,Time_Err,file_xy

          REAL(8) P[ALLOCATABLE](:), T[ALLOCATABLE] (:,:,:)
	    REAL(8) U[ALLOCATABLE](:), V[ALLOCATABLE] (:,:,:)
	    REAL Lat[ALLOCATABLE](:),Lon[ALLOCATABLE](:)
          REAL tt[ALLOCATABLE](:),tr[ALLOCATABLE](:)
          REAL xc[ALLOCATABLE](:),yc[ALLOCATABLE](:)
	    REAL r[ALLOCATABLE](:,:),T_out[ALLOCATABLE](:)
          REAL data_out[ALLOCATABLE](:,:),t_work[ALLOCATABLE](:,:)
          REAL d_st[ALLOCATABLE](:)
100   format(1x,16F8.1)
101   format(1x,I3,2x,I3,2x,'Pressure= ',F7.1) 
102   format(20x,'Pressure= ',F7.1,4x,'time (in hour)',F7.1)   
104   format(1x,f7.1,2x,28F8.2)
105   format(2x,17F7.1)
106   format(2x,'SQD: ',8F9.3)
c ***************************************************************
C          print*,'    Input file_read number'
C	    read*,N_read_file
C     IF N_read_file=1 THEN TEMPERATURE 
C     IF N_read_file=2 THEN WIND - U
C     IF N_read_file=1 THEN WIND - V
C     IF N_read_file=1 THEN HUMIDITY

       N_read_file=1  
       KKK=3 
C	 N_read_file=2 
C	 N_read_file=3 
C	 N_read_file=4 

c       read coord of point interpolating to and number nearest station
c
       CALL READ_TIME(DELTA_0)
C      open(1,file='japan.ini')
C	 read(1,*) delta_0    !time in hours
C       read(1,*) kkk        !Number points taking into account for interpolation
C	close(1)

C     FILES TO WRITE ADDITION INFORMATION

	open(18,file='ht.dat') 
	open(12,file='sq.dat')
	open(13,file='aver.dat') 
	open(21,file='Result.dat') 
c      считывание исходного массива        

	  open(1,file=f_ini(N_read_file))
	  read(1,*) Lat_0,Lon_0 ! Latitude and Longitude for interpolation
	  read(1,*) N_point,N_h,N_time  !number of points (Lat,Lon), Press and time

c      выделяем массив давлений [мб]
	  ALLOCATE(P(1:N_h),stat=ier)
c      выделяем массив данных покрывающего участка
	  ALLOCATE(T(1:N_point, 1:N_time, 1:N_h))
c      выделяем массив координат точек покрывающего участка	  

        ALLOCATE(Lat(1:N_point))
        ALLOCATE(Lon(1:N_point))
	  ALLOCATE(xc(1:N_point))
	  ALLOCATE(yc(1:N_point))
	  ALLOCATE(d_st(1:N_point))
	  ALLOCATE(tt(1:N_point))
	  ALLOCATE(tr(1:N_point))
	  ALLOCATE(t_out(1:N_time))
	  ALLOCATE(data_out(1:N_h,1:N_time))
	  ALLOCATE(r(1:N_point,1:N_point))
	  ALLOCATE(t_work(1:N_time,1:N_point))

	  read(1,*) (P(i),i=1,N_h) ! Pressure

      do kf=1,N_point
         	read(1,*) lat(kf),Lon(kf)
	 do k=1,N_time
		read(1,*) (T(kf,k,j),j=1,N_h)	
	 enddo
      enddo 
	 close(1)

c      interpolating point is surrounded on the 0,0
       nx = sqrt(N_point*1.0)
       ny = nx
	 X_0=0
	 Y_0=0
	
        if(N_read_file.EQ.4) N_h=8
c      Строим координатную сетку для базисных станций

       call xy_coord(N_point,Lat_0,Lon_0,Lat,Lon,xc,yc,r,file_xy)

c      Находим точку сети, ближайшую к интерполируемой точке  
       call number_nearest_point (n_point,Lat_0,Lon_0,Lat,Lon,xc,yc,
     &	                        d_st,jst)
110    format(10x,'number nearest station: ',I3,2x,'distance=',F9.3)
111    format(2x,'lat=',F7.2,4x,'lon=',F7.2,2x,'Dtime=',F9.4,1x,'hour')
112    format(2x,'x=',F9.3,2x,'y=',F9.3,2x)
       
	 write(12,111) Lat(jst),Lon(jst),DELTA_0
	 write(12,112) xc(jst),yc(jst)
       write(12,110) jst,d_st(jst)
	 
	 write(*,111) Lat(jst),Lon(jst),DELTA_0
	 write(*,112) xc(jst),yc(jst)
	write(*,110) jst,d_st(jst)
c      Интерполируем данные базисных станций на заданное время 
c       simple interpolation

	open(15,file='time_d.dat')
        N_ti=N_time-1
	do k=1,N_h 
  	 do i=1,N_point 
        tt(i) = T(i,N_ti,k)
	  tr(i) = T(i,N_time,k)
       enddo 
       
      call inter_japan_time_step(n_point,tr,tt,delta_0)  
c         Заменяем последнее данное для сеточного времени на интерполированное
      
	 do i=1,N_Point 
        T(i,N_ti,k) = T(i,N_ti,k)+tt(i) 
	 enddo 
      enddo !k  
      close(15)

c ********* цикл по высотам *********
      do kh=1,N_h
c ************************************
	write(*,*) '       RESULT for pressure ',P(kh)

c       Формируем массив метеопараметров для высоты P(kf)
c	   open(2,file=file_t)  
c	      write(2,101) N_time,N_point,P(kh)
c	       do i=1,N_time
c              write(2,100) (T(j,i,kh),j=1,N_point)
c	       enddo
c         close(2)
         do i = 1,n_time
	    do j = 1,n_point
           t_work(i, j) = T(j,i,kh)
	    enddo  
    	   enddo  
c      Определяем параметры экспоненты, интерполирующей 
c       зависимость корреляций от расстояния: y=a*exp(-br)        
c
       call regress(N_time,N_point,t_work,r,a_exp,b_exp,Regress_Err)
c
c      Интерполируем заданную точку 
c
        
       Call inter_Gandin(x_0,y_0,N_time,N_point,t_work,jst,xc,yc,
     &	 a_exp,b_exp,t_out,srv,kkk,inter_Error)

	   do i=1,N_time
	 	 data_out(kh,i)=t_out(i)
         enddo
      enddo !kh

      write(18,105) (P(i),i=1,N_h)
	do i=1,N_time
       write(18,105) (data_out(k,i),k=1,N_h)
      enddo  
C       write result to file
       write(21,'(17F9.2)') (data_out(k,N_ti),k=1,N_h)
c      close(17)
	close(18)
	close(12)
	close(13)
	close(20)
	 DEALLOCATE(P,T,xc,yc,lat,lon,tt,tr,t_work,t_out,d_st)

	stop
	end


	subroutine xy_coord(N_point,Lat_0,Lon_0,Lat,Lon,x,y,r,file_xy)
	  parameter (r0 = 6400,d = 2.5)

c     Строит координатную сетку для заданной сетки широт и догот 
	real  f,L,L1,Lat_0,Lon_0,pi
	real Lat(N_point),Lon(N_point),x(N_point),y(N_point)
	real r(N_point,N_point)
	character*20 file_xy

	pi = 4.* ATAN(1.0)  
	L1 = pi * r0 / 180.

      open(10,file=file_xy)
	 
	do i=1,N_point 
	   f= lat(i)
	   L=lon(i)
	   y(i) = l1 * (f - Lat_0)
	   x(i) = (l - Lon_0) * L1 * COS(f * pi / 180.)
         write(10,1) x(i), y(i), f, L
	  enddo
1     format(4F10.3)
	CLOSE(10)
c     Расстояния между точками сетки        
       do i=1,n_point 
        do j=i,n_point
	    r(i,j)=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2)
        enddo
	 enddo  
	 do i = 1, n_point
	  do j = i, n_point 
	    r(j, i)= r(i, j)
        enddo
       enddo 

	return
	end subroutine xy_coord

	

	 

	




