
C     Time-interpolating procedure
c       A linear time interpolation are using to describe the change of a meteoparameter.		
	subroutine inter_japan_time_step(ns,tr,t,delta)
                      
	Real t(ns),tr(ns),d
	
	 do i=1,ns
         d = tr(i)-t(i)
	    t(i)=d*delta/6.
       enddo 
	 write(15,'(9F9.2)') (t(i),i=1,ns)
	write(15,*) ' '
	
       return
	end subroutine inter_japan_time_step
c
       subroutine number_nearest_point(n,Lat_0,Lon_0,Lat,Lon,x,y,d,jst)
     
       Real Lat_0,Lon_0,Lat(n),Lon(n),x(n),y(n),d(n)
	  do i=1,n
	   d(i)=sqrt(x(i)**2+y(i)**2)
        enddo
	  amin=d(1)
	   do i=2,n
	    if(d(i).LT.amin) then
		 amin = d(i)
           num=i
          endif
	   enddo
       	 jst=num
       return
	 end subroutine number_nearest_point

       subroutine READ_TIME(DELTA)

	 CHARACTER *80 STR, SHOUR*2, SMIN*2, SSEC*2 
	 CHARACTER *2 ST,ST1,ST2  
       INTEGER  HOUR,MIN,SEC
	
	 OPEN(1,FILE='TIME.DAT')
	  READ(1,'(A)') STR
         sHOUR=str(31:32)
         sMIN = str(34:35)
	   sSEC = str(37:38) 
	  write(ST, '(A2)') sHOUR
        write(ST1,'(A2)') sMIN
	  write(ST2,'(A2)') sSEC
	  READ(ST, '(i2)') HOUR
	  READ(ST1,'(i2)') MIN
	  READ(ST2,'(i2)') SEC
       CLOSE(1)
	 IF(HOUR.GE.1.AND.HOUR.LE.7) HOUR=HOUR-1
	 IF(HOUR.GE.7.AND.HOUR.LE.13) HOUR=HOUR-7
	 IF(HOUR.GE.13.AND.HOUR.LE.19) HOUR=HOUR-13
	 IF(HOUR.GE.19.AND.HOUR.LE.24) HOUR=HOUR-19
	 IF(HOUR.GE.0.AND.HOUR.LT.1) HOUR=HOUR+6
	 
	 DELTA = HOUR + MIN/60. + SEC / 3600.
	  
	  RETURN

	end subroutine READ_TIME
