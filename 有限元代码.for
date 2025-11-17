!     定义随时间变化的外界温度 （第二类边界条件）
      SUBROUTINE FILM(H,SINK,TEMP,KSTEP,KINC,TIME,NOEL,NPT,
     1 COORDS,JLTYP,FIELD,NFIELD,SNAME,NODE,AREA)

C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION H(2),TIME(2),COORDS(3), FIELD(NFIELD)
      CHARACTER*80 SNAME

	v=2.6                         !v为日平均风速( m/s),根据实际条件确定
	hc=3.7*v+9.4                  !hc为对流系数 (W/m^2 C) W=J/s
	H(1)=3600*hc                  !对流系数（J/（h m^2 C））
	H(2)=0
	Tamax=35.6                    !!日最高气温（C）,根据实际条件确定
	Tamin=22.8                    !!日最低气温（C）,根据实际条件确定
	Ta=(Tamax+Tamin)/2            !日平均气温（C）
	Tm=(Tamax-Tamin)/2            !日气温变化幅度（C）
	w=0.2618                      !频率（pi/12）
	t0=9                          !气温变化时间差(h)
	SN1=SIN(w*(TIME(1)-t0))
	SN2=SIN(2*w*(TIME(1)-t0))
	SINK=Ta+Tm*(0.96*SN1+0.14*SN2) !日气温变化（两个正弦函数的组合表达）
C
      RETURN
      END
      
!     定义随时间变化的热流（第一类边界条件）
      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,
     1 JLTYP,TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION FLUX(2), TIME(2), COORDS(3)
      CHARACTER*80 SNAME

C      user coding to define FLUX(1) and FLUX(2)
	
      FLUX(2)=0
      Qs=26.3E6           !!日太阳辐射总量(J/m^2),根据实际条件确定
      c=10.7              !!日照时间(h)
      m=12.0/c
      q0=0.131*m*Qs       !中午最大辐射(J/s m^2)
      w=0.2618            !频率（pi/12）
      pi=3.14159265
      as=0.90             !!路面的太阳辐射吸收率
	t=TIME(1)-12
      q=q0/(m*pi)
      ak0=q0/pi
      sa0=pi/(2*m)
      do k=1,30         
	  if(k==m)then
	      ak=ak0*( sin((m+k)*sa0)/(m+k)+sa0 )
	  else
            ak=ak0*( sin((m+k)*sa0)/(m+k)+sin((m-k)*sa0)/(m-k) )
        end if
        q=q+ak*cos(k*pi*t/12)   !太阳辐射（Fourier级数表达式）
      end do     
      FLUX(1)=as*q      !进入路面的热流量
      
      RETURN
      END
      
