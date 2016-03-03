c
c     Gandin optimal spatial interpolation procedure
c   
      Subroutine inter_Gandin(x_0,y_0,m0,ns,t,jst,cx,cy,
     * a_exp,b_exp,t_out,srv,kkk,inter_Error)
	use MSFLIB
c
c    INPUT PARAMETERS:
c      t - array of meteoparameters
c      m0 - amount members in time series
c      ns - amount points in a net
c      jst - number station nearest to target point (with coordinate x,y;if ns=81 then jst=41)
c      x_0,y_0 - coordinates target point(x_0=0, y_0=0)
c      cx,cy - coord. all stations 
c      a_exp,b_exp - factors of the equation to exponential regression 
c      kkk - amount affecting station 
c    OUTPUT DATA:
c      t_out - output array containing interpolated parameters for each time
c      srv - standard deviation
c      inter_Error - quality to interpolations
c
  

	Real t(m0,ns+1),cx(ns),cy(ns),srv(kkk),t_out(m0),dds(15)
	Real a_exp,b_exp

	Character*20 read_t,read_xy,inter_Error

      Integer nns[ALLOCATABLE](:),nsk[ALLOCATABLE](:),n[ALLOCATABLE](:) 
	REAL tr[ALLOCATABLE](:,:),tt[ALLOCATABLE](:,:)
	REAL tw[ALLOCATABLE](:,:),rd[ALLOCATABLE](:,:)
	REAL t_res[ALLOCATABLE](:,:),r[ALLOCATABLE](:,:)
	REAL ds[ALLOCATABLE](:,:)
	REAL c[ALLOCATABLE](:,:),c1[ALLOCATABLE](:,:)
	REAL a[ALLOCATABLE](:,:),ao[ALLOCATABLE](:,:)
      REAL p[ALLOCATABLE](:,:)
	REAL rd0[ALLOCATABLE](:),rd1[ALLOCATABLE](:),xv[ALLOCATABLE](:)
      REAL s_int[ALLOCATABLE](:),sred[ALLOCATABLE](:)
	REAL xk[ALLOCATABLE](:),yk[ALLOCATABLE](:),h[ALLOCATABLE](:)
	REAL b[ALLOCATABLE](:),a1[ALLOCATABLE](:),b1[ALLOCATABLE](:)
	REAL rs[ALLOCATABLE](:),sr[ALLOCATABLE](:),s[ALLOCATABLE](:)
	REAL g[ALLOCATABLE](:),dis[ALLOCATABLE](:),x[ALLOCATABLE](:)
          
      inter_Error='OK'

      ALLOCATE(tr(1:m0,1:ns+1))  
	ALLOCATE(tt(1:m0,1:ns+1))
	ALLOCATE(tw(1:m0,1:ns+1))
      ALLOCATE(rd(1:ns,1:ns)) 
      ALLOCATE(rd0(1:ns))
	ALLOCATE(rd1(1:ns))
	ALLOCATE(xv(1:ns))
	ALLOCATE(s_int(1:ns))
	ALLOCATE(sred(1:ns))
      ALLOCATE(t_res(1:m0,1:kkk))  
      ALLOCATE(r(1:ns,1:ns))
      ALLOCATE(xk(1:ns))
	ALLOCATE(yk(1:ns))  
	ALLOCATE(nns(1:ns))
	ALLOCATE(nsk(1:ns))
      ALLOCATE(h(1:ns))
	ALLOCATE(a1(1:ns))
	ALLOCATE(b1(1:ns))
	ALLOCATE(rs(1:ns))
	ALLOCATE(s(1:ns))
	ALLOCATE(sr(1:ns))
	ALLOCATE(g(1:ns))
	ALLOCATE(dis(1:ns))
	ALLOCATE(x(1:ns))
	ALLOCATE(b(1:ns))
	ALLOCATE(ds(1:ns,1:15))
c	ALLOCATE(n(1:12,1:ns))
	ALLOCATE(c(1:ns,1:ns))
	ALLOCATE(c1(1:ns,1:ns))
	ALLOCATE(a(1:ns,1:ns)) 
	ALLOCATE(ao(1:ns,1:ns))
	ALLOCATE(P(1:ns,1:ns))
	 
c	open(2,file=read_t)  
c	   read(2,*) m0,ns 
c       	do i = 1,m0
c		  read(2,*) (t(i, j),j = 1,ns)
c	      write(*,*) (t(i, j),j = 1,ns)
c	    enddo
c	close(2)

c       tr- initial array 
	 	do i = 1,m0
           do j = 1,ns
	       tr(i,j)=t(i,j)
		 enddo
		enddo    

	do  i = 1, ns 
	 nsk(i) = i
	enddo 

c          pасчет средних значений по мес€цам

	do j = 1, ns
	 s(j)=0
	 do i = 1, m0
	  s(j) = s(j) + t(i, j)
       enddo
	enddo 

	do j = 1, ns
	 s(j) = s(j) / m0
	enddo

c             вычисление отклонений

	do j = 1, ns
	 do i = 1, m0
	 t(i, j) = t(i, j) - s(j)
	 enddo
	enddo

c           считываем координаты станций
    
c	OPEN (2,file=read_xy)
c	do i = 1, ns
c	 read(2,*)  cx(i),cy(i),lat(i),long(i)
C		 write(*,*)  cx(i),cy(i),lat(i),long(i)
c	enddo
c	CLOSE (2)

c       –ассто€ни€ между станци€ми  
	do k = 1,ns
	 do j = 1, k
		rd(k, j) = SQRT((cx(k) - cx(j))**2 + (cy(k) - cy(j))**2)
		rd(j, k) = r(k, j)
	 enddo !j
	enddo !k

c	  calculates distances between all stations and needed
	do k = 1,ns
	 rd0(k)=SQRT((cx(k) - x_0)**2 + (cy(k) - y_0)**2)
      enddo  

C расчет ковариационной матрицы дл€ определени€ вли€ющих станций

	do k = 1, ns
	 do k1 = k, ns
	  ss = 0
	   aa = 0 
	    bb = 0
		 aa1 = 0
		  bb1 = 0
	     do i = 1, m0
			aa = aa + t(i, k)
		    bb = bb + t(i, k1)
		    aa1 = aa1 + t(i, k)** 2
		    bb1 = bb1 + t(i, k1)** 2
			ss = ss + t(i, k) * t(i, k1)
		 enddo !i
		aa = aa / m0
		aa1 = aa1 / m0
		bb = bb / m0
		bb1 = bb1 / m0
		 ss = ss / m0
	 if(aa1.NE.0.AND.bb1.NE.0) then
	  r(k, k1) = (ss - aa * aa1) / SQRT((aa1 - aa**2)*(bb1 - bb**2))
	 else
        r(k, k1) =0.01
	 endif
	r(k1, k) = r(k, k1) 
       enddo !k1
	enddo  !k

c            Ќачало цикла по  числу соседей n0

	do n0 = 3, kkk

            ns_1=1
C                       определение вли€ющих станций
          j5=JST        ! ЅЋ»∆ј…Ўјя   »Ќ“≈–ѕќЋ»–”≈ћќ… “ќ„ ≈ —“јЌ÷»я 
	I1 = 0
55    do J1 = 1,ns
        IF (J1.EQ.j5) go to  60
          I1 = I1 + 1
		rs(I1) = r(j5, J1)
          nns(I1) = nsk(J1)
60    enddo
       nn = I1

	
c      ”пор€дочивание массива коррел€ций

	call YPOR_2(1,Nn,ns,rs,nns)

c          формирование нового массива

	do i = 1, m0
	 do jj = 1, n0
		J0 = nns(jj)
		tt(i, jj) = t(i, J0)
	    rd1(jj)=rd0(j0)
	 enddo ! jj
	  tt(i, n0 + 1) = t(i, j5)
	enddo ! i

	do k = 1,n0 + 1
	 do k1 = k, n0 + 1
		ss = 0 
		do i = 1, m0
			ss = ss + tt(i, k) * tt(i, k1)
	    enddo !i
	   a(k, k1) = ss 
	   a(k1, k) = a(k, k1) 
       enddo !k1
      enddo !k
	 
	do i = 1, n0 + 1
	 do jj = 1,n0 + 1
		P(i, jj) = a(i, jj)
       enddo !jj
      enddo !i
c                   ÷икл по J2
	do j2 = 1, n0
	 do i = 1, n0
	  b(i) = 0
	 enddo ! i
	    b(j2) = a_exp*exp(b_exp*rd1(j2)) ! интерпол€ционные параметры из graph
c		b(j2) = 1.
	 do J3 = 1,n0
	  do J4 = 1, n0
		a(J3, J4) = P(J3, J4)
	  enddo !j4
       enddo !j3

1440  N1 = n0 - 1
	do k = 1, N1
1450   IF (ABS(a(k, k)).GT.0) go to 1530
1460    k1 = k + 1
        do m = k1,n0
1470     IF(ABS(a(m, k)).GT.0) go to 1490
1480      GOTO 1510
1490     do l = 1,n0
           U = a(k, l)
		 a(k, l) = a(m, l)
1500       a(m, l) = U
         enddo ! l
1510    continue
	 enddo ! m
1520    U = b(k)
        b(k) = b(m)
	  b(m) = U
1530    g(k) = b(k) / a(k, k)
         k1 = k + 1
1540    do i = k1,n0
         b(i) = b(i) - a(i, k) * g(k)
1550     do J1 = k,n0
           j = n0 - J1 + k
		 c(k, j) = a(k, j) / a(k, k)
1560       a(i, j) = a(i, j) - a(i, k) * c(k, j)
1570     enddo !J1
        enddo ! i
	enddo ! k
1580     m = n0
         X(m) = b(m) / a(m, m)
1590     m = m - 1
         ss = 0
	    do l = m,N1
1600       ss = ss + c(m, l + 1) * X(l + 1)
          enddo ! l
1610     X(m) = g(m) - ss
         IF( m.GT.1) go to 1590

	 do i = 1,n0
	  ao(j2, i) = X(i)
	 enddo ! i
	enddo !j2

c                онец ÷икла по J2


c         PRINT *,'form array free members'

	 do i = 1,n0
	  b(i) = a(i, n0 + 1)
c	PRINT *,B(I)
	 enddo ! i
c          vesovye coefficients

	do i = 1,n0
	 ss = 0
	  do j = 1,n0
		ss = ss + ao(i, j) * b(j)
	  enddo !j
	   X(i) = ss
	 enddo ! i

c                 регрессионный прогноз
c        нормировка и отбор весовых коэффициентов 
c _______________________________________________________	    
       jx=0  
	do i = 1,n0
       Xv(i)=0
       if(X(i).GT.0) then
	  jx=jx+1
	   Xv(jx)=X(i)
       endif
	enddo 
	ss=0 
      do jt=1,n0
	 ss=ss +Xv(jt)
	enddo
	do jt=1,n0
	 xv(jt)=xv(jt)/ss
	 X(jt)=xv(jt)
	enddo

c	    _______________________________________________________	    
	
	do i = 1,m0
	 ss = 0
	 do jj = 1,n0
	  ss = ss + X(jj) * tt(i, jj)
	 enddo ! jj

c     рассчитываем значени€ температуры (метеопараметра)
          t_out(i)=ss
	    t_res(i,n0)=ss 
c	go to 777
	enddo ! i


c      интерполированное среднее значение
c       средние значение станций, участвующих в интерпол€ции         
	   s_int(n0)=0
           do jj = 1,n0
	      s_int(n0) = s_int(n0) + X(jj) * s(nns(jj))
           enddo
c	  print *, 'n0= ',n0, ' average = ',s_int(n0)

            
c          	s_int(n0) - интерполированное среднее значение искомой станции 
	do i = 1,m0
       t_out(i)=t_out(i)+S_int(n0)
	 t_res(i,n0)=t_res(i,n0) + S_int(n0)
      enddo   
	 
      dss = 0
	ssr = 0
	sred=0   
         do i=1,m0
	    del = tr(i, JST)- t_res(i,n0) 
		dss = dss + del**2
         enddo
	    dss = dss / m0
	    srv(n0)=dss
		srv(n0)=sqrt(srv(n0))
	enddo ! n0

45    format(45F8.3)
16    format(16F8.1)		    
100   format(F10.3,4x,10F9.3)
101   format(45F9.3)
102   format(8x,6F9.2,5x,'real ',2x,F9.3,4x,F7.2)
103   format('number point= ',I2,2x,' x= ',F9.2,2x,'y= ',F9.2)
104   format( ' sq.err= ',6F9.3)
105   format( ' true average= ',F9.3,2x,' average= ',6F9.3)

c       Writing Results

	do i = 1,m0
	 delta=t_res(i,kkk)+s_int(kkk)-tr(i,jst)
	enddo 	
	 write(12,104) (srv(n0),n0=3,kkk)
	 write(13,105) s(jst),(s_int(n0),n0=3,kkk)

777    continue
        do i = 1,m0
          t_out(i)=t_res(i,kkk)
	  enddo   
c		PRINT*, 'End of work'
		return
       
	END Subroutine inter_Gandin

	 SUBROUTINE YPOR_2(N1,N2,NPW2,XK,Yk)
      DIMENSION XK(NPW2),YK(NPW2)
      Integer Yk,YKB
      l=0
      MMM=N2-1
    1 IND=0
      DO 2 I=N1,MMM
       XKB=XK(I)
	 YKB=YK(I)
        IF( XKB.GE.XK(I+1) ) GO TO 2
        UBUB=XK(I+1)
	  IU= YK(I+1)
         XK(I)=UBUB
         XK(I+1)=XKB
         YK(I)=IU
         YK(I+1)=YKB 
       IND=1
    2 CONTINUE
      IF(IND.EQ.1)  Then
        l = l+1
c        print *,'l= ',l
       GO TO 1
                    else
      endif
      RETURN
      end  SUBROUTINE YPOR_2
