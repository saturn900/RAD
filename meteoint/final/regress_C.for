c
c       Procedure of the finding to exponential regression 
c          for dependency correlation from distances
c
	subroutine regress(nr,nz,t,r0,a_exp,b_exp,Regress_Err)
	use MSFLIB

c    INPUT PARAMETERS:
c	t - array of meteoparameters
c     r0 - distances between station 
c     nr - amount members in time series
c     nz - amount points in a net
c    OUTPUT DATA:
c     a_exp, b_exp - factors of the equation to exponential regression  a*exp(-b*r)
c     Cdet =R**2 (R - factor to correlations)
c     Er - a root-mean-square mistake at interpolations correlation
c     Regress_Err - quality to interpolations 

       Parameter ( aa = 999)
       Integer z,i,j,kk
	 Real  a_exp,b_exp,xx,yy
       Real  L0(12),t(nr, nz), r0(nz, nz)
       Character*20 d$,Regress_Err,mes(5)

	INTEGER NN[ALLOCATABLE] (:),MM[ALLOCATABLE] (:)
	REAL    p0[ALLOCATABLE](:,:)
	REAL x[ALLOCATABLE](:),y[ALLOCATABLE](:),E[ALLOCATABLE](:)
	
	DATA mes/'- exellent','- good','- satisfact','- bad','- very bad'/

      d$ = 'inf_regress.dat'
      Regress_Err ='OK'

	  ALLOCATE(NN(1:(nz ** 2- nz)/2),stat=ier)
	  ALLOCATE(MM(1:(nz ** 2- nz)/2),stat=ier)
	  ALLOCATE(p0(1:nz,1:nz))
	  ALLOCATE(x(1:(nz ** 2- nz)/2),stat=ier)
	  ALLOCATE(y(1:(nz ** 2- nz)/2),stat=ier)
	  ALLOCATE(E(1:(nz ** 2- nz)/2),stat=ier)
        n=nr
	  z=nz 

	do j = 1 , z 
	 do k = j , z
	  a = 0
	  b = 0
	  a1 = 0
	  b1 = 0
	  d = 0
	   do i = 1, n
		xx = t(i, j)
		yy = t(i, k)
		a = a + xx
		a1 = a1 + xx ** 2
		b = b + yy
		b1 = b1 + yy ** 2
		d = d + xx * yy
	   enddo
		a = a / n
		b = b / n
		a1 = a1 / n - a ** 2
		b1 = b1 / n - b ** 2
		d = d / n
		if(a1.NE.0.AND.b1.NE.0) then
		 p0(j, k) = (d - a * b) / SQRT(a1 * b1)
	    else
           p0(j, k) =0.01
		endif
	  enddo
	enddo	
c         определ€ем величину максимального и минимального элементов массива данных
c                              p0-минимум;p9-максимум
		p00 = p0(1, 1)
		P9 = p0(1, 1)
		do i = 1 ,z
           do j = i ,z
         IF (p0(i, j).EQ.aa ) go to 190
         IF (p0(i, j).LE.0  ) go to 190
		  IF (p0(i, j).LE.p00 ) p00 = p0(i, j)
            IF (p0(i, j).GE.P9 ) P9 = p0(i, j)
 190      continue
           enddo
		enddo 
200     continue
		r00 = .5
		r9 = .5
		do i = 1, z
		 do j = i, z
			IF (r0(i, j).LE.0 ) go to  199
			IF (r0(i, j).LT.r00) r00 = r0(i, j) 
			IF (r0(i, j).GT.r9) r9 = r0(i, j)
199       continue
            enddo
	    enddo
c _________________________________
       Rm = r9 * 0.9
	 CC=0
c _________________________________

       i1 = 0
       N1 = 0
        do kk = 1,z*z/2-z
	   x(kk)=0
	   y(kk)=0 
        enddo
        do i2 = 1,z
	   do j2 = i2,z
           IF (i2.EQ.j2) go to 310
	     IF (i2.EQ.CC) go to 310
	     IF (j2.eq.CC) go to 310
		 IF (p0(i2, j2).LE.0) go to 310
c           IF p0(i, j) = AA THEN 310
           IF (RM.LT.r0(i2, j2))  go to 310 
c		    go to 310
c		  ELSE 
		  i1 = i1 + 1
		  N1 = N1 + 1
c	     Endif
            x(i1) = r0(i2, j2)
		  y(i1) = p0(i2, j2)
		  NN(i1) = i2
		  MM(i1) = j2
310     continue
	   enddo
	  enddo 
          j = N1

	      y0 = p00
		  Y9 = P9
		  X0 = 0
		  X9 = RM
	      y0 = 0

c        Regression by exponent function 

        x1 = 0
	  X2 = 0
	  Y1 = 0
	  Y2 = 0
	  p = 0
       do i = 1, j
	    xx = x(i)
		yy = y(i)
		if(xx.LE.0.OR.yy.LE.0.or.xx.LE.1.E-10.OR.yy.LE.1.E-10) then
	write(*,*) 'Program Interpol(regression): Error! Negative values 
     &for distances'
		 write(*,*) i, x(i), y(i)
	     endif

		yy = ALOG(yy)
		x1 = x1 + xx
		Y1 = Y1 + yy
		X2 = X2 + xx * xx
		Y2 = Y2 + yy * yy
		p = p + xx * yy
	 enddo
	 V = (x1 * Y1 - j * p) / (x1 * x1 - j * X2)
	 W = (Y1 - V * x1) / j
	 R = (p - x1 * Y1 / j) / SQRT((X2-x1*x1/j)*(Y2-Y1*Y1/j))
	 a = W
	 b = V
	 a = EXP(W)

471	Y3 = 0
	Y4 = 0
	P1 = 0
	Y5 = 0
	Y1 = 0
	Y2 = 0
472	do i = 1,N1

	 xx = x(i)
	 yy = y(i)
	 Y1 = Y1 + yy
	 Y2 = Y2 + yy * yy
	 YR = a * EXP(b * xx)
	 Y3 = Y3 + YR
	 Y4 = Y4 + YR * YR
	 P1 = P1 + yy * YR 
	 Y5 = Y5 + (yy - YR) * (yy - YR) 
	 E(i) = ABS(yy - YR)
	enddo 
	 RR = (P1 - Y1 * Y3 / N1) / SQRT((Y2-Y1*Y1/N1)*(Y4-Y3*Y3/N1))
	 CDet=RR*RR
c	  Multi. correlation :  RR
c	  Coeff. determination: CDet
c	  Average square error:  ER

	 ER = SQRT(Y5 / N1)

	If(CDET.LE.1.AND.CDET.GT.0.85)    Regress_Err=mes(1)
	If(CDET.LE.0.85.AND.CDET.GT.0.70) Regress_Err=mes(2)
	If(CDET.LE.0.70.AND.CDET.GT.0.55) Regress_Err=mes(3)
	If(CDET.LE.0.55.AND.CDET.GT.0.40) Regress_Err=mes(4)
	If(CDET.LE.0.40.AND.CDET.GT.0.)   Regress_Err=mes(5)

	OPEN (4,file=d$)
45    format(45F8.1)
46    format(2I3,F8.3)
39    format('      »спользуетс€ экспоненциальна€ регресси€')
38    format('      –езультаты расчета программы GRAPH_EXP')
40    format('радиус отбора= ',F9.3)
41    format('x0=', F9.3,2x,'x9= ',F9.3,2x,'y0 = ',F9.3,2x,'y9= ',F9.3)
42    format('a=', F12.5,2x,'b= ', F12.5,2x, 'число точек = ', I4)
43    format('коэфф. коррел€ции R=',F12.5,' радиус коррел€ции= ',F9.3)
44    format('коэфф. детерминации CDET= ',F12.4,2x,A20)
47    format('средн€€ квадрат.ошибка = ',F12.4)
      write(4,*) '                     '   
	write(4,39) 
	write(4,38) 
	write(4,40) RM
	write(4,41) X0, X9, y0,  Y9
	write(4,42) a, b, N1
	write(4,43) R, 1/b
	write(4,44) cdet,Regress_Err
	write(4,47) ER

       a_exp=a
 	 b_exp=b
       do i = 1,N1
	 IF( E(i).GT.2 * ER) THEN 
	   Regress_Err='Bad'
	   ELSE 
	 go to 1580
1580   continue
       endif
	 enddo 
	 DEALLOCATE(x,y,E,p0,MM,NN)  
	return
	End subroutine regress




