      Program SCX_Loss
C
C  Simulate the concentration of aluminum, oxygen and aluminum oxide in bond and substrate.
C
      implicit none
C
C  Declare variables
C
      integer NDIM                     ! Number of points throughout thickness
      parameter (NDIM=401)             ! Sets dimension of array
      real*8 Ox(0:NDIM)                ! Oxygen concentration
      real*8 Al(0:NDIM)                ! Aluminum concentration
      real*8 AlOx(0:NDIM)              ! Aluminum oxide concentration
      real*8 Oxu(0:NDIM)               ! Updated oxygen concentration
      real*8 Alu(0:NDIM)               ! Updated aluminum concentration
      real*8 AlOxu(0:NDIM)             ! Updated aluminum oxide concentration
      real*8 dx                        ! X direction spacing
      real*8 dt                        ! Time step
      real*8 DAl(0:NDIM)               ! Aluminum diffusion coefficient
      real*8 DAl0(0:NDIM)              ! Aluminum diffusion - 0 mole fraction Al
      real*8 DAl1(0:NDIM)              ! Aluminum diffusion - 1 mole fraction Al
      real*8 DOx(0:NDIM)               ! Oxygen coefficient
      real*8 testAl(0:NDIM)            ! Time steps for Aluminum
      real*8 testOx(0:NDIM)            ! Time steps for Oxygen
      real*8 Al100                     ! Molar concentration of 100% aluminum
      real*8 DAl0_0                    ! Aluminum diffusion - 0 mole fraction Al
      real*8 DAl1_0                    ! Aluminum diffusion - 1 mole fraction Al
      real*8 DOx_0                     ! Oxygen coefficient
      real*8 DAl0_0InAL2O3             ! Aluminum diffusion - 0 mole fraction Al
      real*8 DAl1_0InAL2O3             ! Aluminum diffusion - 1 mole fraction Al
      real*8 DAl0InAL2O3(0:NDIM)        ! Aluminum diffusion - 0 mole fraction Al
      real*8 DAl1InAL2O3(0:NDIM)       ! Aluminum diffusion - 1 mole fraction Al
	  real*8 DAlInAl2O3(0:NDIM)        ! Aluminum diffusion in Al2O3
      real*8 DOx_0InAL2O3              ! Oxygen coefficient
      real*8 DOxInAL2O3(0:NDIM)        ! Oxygen coefficient
      real*8 R0                        ! Reaction coefficient
	  real*8 RP                        ! Reference reaction coefficient
	  real*8 nAl                       ! ALuminum order (power)
	  real*8 nO                         ! Oxygen order (power)
      real*8 R(0:NDIM)                 ! Reaction rate
      real*8 QAl                       ! Aluminum diffusion energy
      real*8 QOx                       ! Oxygen diffusion energy
      real*8 DG0                       ! Reaction Gibbs Free energy
      real*16 TAl                      ! Aluminum reference temperature
      real*16 TOx                      ! Oxygen reference temperature
      real*16 TG0                      ! Gibbs Free Energy reference temperature
      real*8 RT                        ! Gas Constant (8.3144621 J/mol-K)
      real*8 h                         ! Bond layer thickness
      real*8 AlB                       ! Inital aluminum concentration in bond
      real*8 AlSCX                     ! Inital aluminum concentration in SCX
      real*8 AlOcr                     ! Critical concentration of Al2O3
      real*8 Ox0                       ! Oxygen concentration at surface
      real*8 x(0:NDIM)                 ! Distance from surface
      real*16 Temp(0:NDIM)             ! Temperature
      real*8 t                         ! Current time
      real*8 tf                        ! Final time
      real*8 interpolate               ! Interpolation value
      real*8 xa(7),ya(7)               ! Tabular data
      real*8 signum                    ! Signum function
      logical writeit                  ! Flag to write output
	  integer NoWrite                  ! Number of output times
	  real*8 time_interval             ! Time intervals between output
	  real*8 last_write                ! Last write time  
	  real*8 next_write                ! Next write time
      integer i                        ! Counter
C
C  Define fixed parameters
C
      RT=8.3144621D3                   ! J/mol-K
      Al100=0.1D0                      ! moles/cm^3
      xa=(/0.00, 0.07, 0.19, 0.29, 0.36, 0.40, 0.48/)				! Aluminum Atomic Fraction
      ya=(/2.E-10, 2.E-10, 3.E-11, 3.E-11, 3.E-11, 1.E-11, 3.E-12/)	! Aluminum Diffusion Coefficient in Ni @ 1323K
C
C  Get parameters
C
      open (unit=505,file='SCX_Parameters.inp', status='old')
      read (505,*) DAl0_0,DAl1_0,DOx_0
      read (505,*) DAl0_0InAl2O3,DAl1_0InAl2O3,DOx_0InAl2O3
	  read (505,*) RP,nAl,nO
      read (505,*) AlB,AlSCX,Ox0
      read (505,*) QAl,QOx,DG0
      read (505,*) TAl,TOx,TG0
      read (505,*) h,AlOcr
      close (Unit=505)
C
C  Convert from reference reaction to actual reaction coefficient
C
      R0=RP*AlB**(2.-nAl)*Ox0**(3.-nO)
C
C  Get and echo the spacing, final time and number of tables
C
      write (6,'(A,$)') 'Please enter x interval: '
      read (5,*) dx
	  write (6,*) dx
      write (6,'(A,$)') 'Please enter final time: '
      read (5,*) tf
	  write (6,*) tf
      write (6,'(A,$)') 'Please enter number of output times: '
      read (5,*) NoWrite
	  write (6,*) NoWrite
	  write (6,'(/)')
	  time_interval=tf/float(NoWrite) 
C
C  Echo parameters
C
      write (6,'(A)') 'Diffusion in SCX'
      write (6,'(3(A,1PE11.3)/)')
     $   'DAl0_0 =',DAl0_0,'  DAl1_0 = ',DAl1_0,'  DOx_0 = ',DOx_0
      write (6,'(A)') 'Diffusion in AlOx'
      write (6,'(3(A,1PE11.3)/)')
     $   'DAl0_0 =',DAl0_0InAL2O3,'  DAl1_0 = ',DAl1_0InAL2O3
     $ ,                          '  DOx_0 = ',DOx_0InAL2O3
      write (6,'(3(A,1PE11.3)/)')
     $   'RP     =',RP    ,'  nAl    = ',nAl   ,'  nO    = ',nO
      write (6,'(3(A,1PE11.3)/)')
     $   'AlB    =',AlB   ,'  AlSCX  = ',AlSCX ,'  OxO   = ',Ox0
      write (6,'(3(A,1PE11.3)/)')
     $   'QAl    =',QAl   ,'  QOx    = ',QOx   ,'  DG0   = ',DG0
      write (6,'(3(A,1PE11.3)/)')
     $   'TAl    =',TAl   ,'  TOx    = ',TOx   ,'  TG0   = ',TG0
      write (6,'(3(A,1PE11.3)/)')
     $   'h      =',h     ,'  AlOcr  = ',AlOcr
	  write (6,'(/)')
C
C  Find time step
C
      dt=dx**2/4/max(DAl0_0,DAl1_0,DOx_0)
C
C  Initialize variables
C
      t=0.
      x(0)=0.
      Ox(0)=Ox0
      Al(0)=AlB*Al100
      AlOx(0)=0
      do i=1,NDIM
         x(i)=x(i-1)+dx
         Ox(i)=0.
         AlOx(i)=0.
         if (x(i).le.h) then
            Al(i)=AlB*Al100
         else
            Al(i)=AlSCX*Al100
        end if
      end do
C
C  Get the temperature at each point
C
      call GetTemp(x,Temp,t,NDIM)
C
C  Find the base aluminum and oxygen diffusion coefficients
C

      call GetDiffusion
     $     (DAl0,DAl0InAl2O3,Temp,DAl0_0,DAl0_0InAL2O3,QAl,TAl,RT,NDIM)
      call GetDiffusion
     $     (DAl1,DAl1InAl2O3,Temp,DAl1_0,DAl1_0InAL2O3,QAl,TAl,RT,NDIM)
	  DAlInAl2O3=DAl0InAl2O3+Al/Al100*(DAl1InAl2O3-DAL0InAl2O3)
      call GetDiffusion
     $     (DOx,DOxInAl2O3,Temp,DOx_0,DOx_0InAL2O3,QOx,TOx,RT,NDIM)
C
C  Find aluminum diffusion coefficients in Nickel
C
      do i=0,NDIM
	    DAl(i)=interpolate(xa,ya,Al(i)/Al100,7)*DAl0(i)/ya(1)
	  end do
C
C  Add effects of Al2O3
C
      DAl=DAl/(1.+DAl/DAlinAl2O3*(1.-Al/(Al+AlOx)))
      DOx=DOx/(1.+DOx/DOxinAl2O3*(1.-Ox/(Ox+AlOx)))
C
C  Find the aluminum and oxygen reaction coefficients
C
      call GetReaction(R,Temp,R0,DG0,TG0,RT,NDIM)
C
C  Update the time step
C
      testAl=1./(2.*DAl/dx**2+abs(nAl*R*Al**(NAl-1)*Ox**NO))
      testOx=1./(2.*DOx/dx**2+abs(nO*R*Al**NAl*Ox**(nO-1)))
      dt=0.5*min(minval(testAl),minval(testOx))
C
C  Write the initial output
C 
      write (6,'( /// 2(A,1PE12.3))') '   time =',t,
     $                                '   Surface temperature =',Temp(0)
      write (6,'(3X,A,11X,A,4X,A,6X,A)') 
     $      'x','Aluminum','Oxygen','AlOx'
      do i=0,NDIM
            write (6,'(1P8E12.3,I5)') x(i),Al(i),Ox(i),AlOx(i)
      end do
C
C  Loop over time
C
      do while (t.le.tf.and.t.ge.0.)
C
C  Update the time 
C
         t=t+dt
C
C  Is this a write interval?
C
         last_write=int(t/time_interval)*time_interval
		 next_write=last_write+time_interval
         writeit=((t.gt.next_write-dt)
     $          .and.(t.le.next_write))
     $          .or.(t.eq.tf)
C
C  Write heading
C
         if (writeit) then
            write (6,'(/// 3X,A,1PE12.3)') 'time =',t
            write (6,'(3X,A,11X,A,4X,A,6X,A)') 
     $            'x','Aluminum','Oxygen','AlOx'
         end if
C
C  Get the temperature at each point
C
         call GetTemp(x,Temp,t,NDIM)
C
C  Find the base aluminum and oxygen diffusion coefficients
C
         call GetDiffusion
     $     (DAl0,DAl0InAl2O3,Temp,DAl0_0,DAl0_0InAL2O3,QAl,TAl,RT,NDIM)
         call GetDiffusion
     $     (DAl1,DAl1InAl2O3,Temp,DAl1_0,DAl1_0InAL2O3,QAl,TAl,RT,NDIM)
	     DAlInAl2O3=DAl0InAl2O3+Al/Al100*(DAl1InAl2O3-DAL0InAl2O3)
         call GetDiffusion
     $     (DOx,DOxInAl2O3,Temp,DOx_0,DOx_0InAL2O3,QOx,TOx,RT,NDIM)
C
C  Find aluminum diffusion coefficients
C
      do i=0,NDIM
	    DAl(i)=interpolate(xa,ya,Al(i)/Al100,7)*DAl0(i)/ya(1)
	  end do
C
C  Add effects of Al2O3
C
         DAl=DAl/(1.+DAl/DAlinAl2O3*(1.-Al/(Al+AlOx)))
         Where (Ox+AlOx.gt.0.) 
     $          DOx=DOxinAl2O3*(1.-Ox/(Ox+AlOx))+DOx*Ox/(Ox+AlOx)
C	            where (DOx.gt.5.*DOxInAl2O3) DOx=-DOx_0
C
C  Find the aluminum and oxygen diffusion coefficients
C
         call GetReaction(R,Temp,R0,DG0,TG0,RT,NDIM)
C
C  Update the time step
C
         testAl=1./(2.*DAl/dx**2+nAl*R*abs(Al**(NAl-1)*Ox**nO))
         testOx=1./(2.*DOx/dx**2+nO*R*abs(Al**NAl*Ox**(nO-1)))
         dt=0.5*min(minval(testAl),minval(testOx))
C
C  Upddate aluminum concentration
C
         if (AlOx(0).ge.AlOcr) then
		    Alu(0)=0.
         else
            Alu(0)=Al(0)+dt/2./dx**2*((DAl(1)+DAl(0))*(Al(1)-Al(0))
     $                            -(DAl(0)+DAl(1))*(Al(0)-Al(1)))
     $                            -nAl*R(0)*Al(0)**nAl*Ox(0)**nO*dt
         end if
         do i=1,NDIM-1
		    if (AlOx(i).ge.AlOcr) then
			   Alu(i)=0.
            elseif (AlOx(i-1).lt.AlOcr) then
               Alu(i)=Al(i)+dt/2./dx**2
     $            *((DAl(i+1)+DAl(i))*(Al(i+1)-Al(i))
     $            -(DAl(i)+DAl(i-1))*(Al(i)-Al(i-1)))*signum(Al(i-1))
     $            -nAl*R(i)*Al(i)**nAl*Ox(i)**nO*dt
	        else
               Alu(i)=Al(i)+dt/2./dx**2
     $            *((DAl(i+1)+DAl(i))*(Al(i+1)-Al(i))
     $            -(DAl(i)+DAl(i+1))*(Al(i)-Al(i+1)))
     $            -nAl*R(i)*Al(i)**nAl*Ox(i)**nO*dt
            end if
         end do
            if (AlOx(NDIM).ge.AlOcr) then
               Alu(NDIM)=0.
            else
               Alu(NDIM)=Al(NDIM)+dt/2./dx**2*
     $            ((DAl(NDIM-1)+DAl(NDIM))*(Al(NDIM-1)-Al(NDIM))
     $            -(DAl(NDIM)+DAl(NDIM-1))*(Al(NDIM)-Al(NDIM-1)))
     $            *signum(Al(NDIM-1))
     $            -nAl*R(NDIM)*Al(NDIM)**nAl*Ox(NDIM)**nO*dt
	    end if
C
C  Update aluminum concentration
C		 
         Oxu(0)=Ox0
         do i=1,NDIM-1
            if (AlOx(i).ge.AlOcr) then
               Oxu(i)=Ox0
			else if (DOx(i).lt.0.) then
			   Oxu(i)=Ox0
            else
               Oxu(i)=Ox(i)+dt/dx**2*(DOx(i)*(Ox(i+1)-Ox(i))
     $               -DOx(i)*(Ox(i)-Ox(i-1)))
     $               -nO*R(i)*Al(i)**nAl*Ox(i)**nO*dt
            end if
         end do
            if (AlOx(NDIM).ge.AlOcr) then
               Oxu(NDIM)=Ox0
			else if (DOx(i).lt.0.) then
			   Oxu(i)=Ox0
            else
               Oxu(NDIM)=Ox(NDIM)+dt/dx**2
     $                  *(DOx(NDIM)*(Ox(NDIM-1)-Ox(NDIM))
     $                  -DOx(NDIM)*(Ox(NDIM)-Ox(NDIM-1)))
     $                  *signum(Ox(NDIM-1))
     $                  -nO*R(NDIM)*Al(NDIM)**nAl*Ox(NDIM)**nO*dt
            end if
C
C  Update aluminum oxide concentration
C	     
            where (AlOx.ge.AlOcr)
               AlOxu=AlOcr
            elsewhere
               AlOxu=AlOx+R*Al**nAl*Ox**nO*dt
            end where
C
C  Move updated values to current values
C
         Al=Alu
         Ox=OXu
         AlOx=AlOxu
C
C  Write the output
C
         if (writeit) then
            do i=0,NDIM
               write (6,'(1P8E12.3,I5)') x(i),Al(i),Ox(i),AlOx(i)
            end do
         end if
      end do
C
C Write final interval
C 
      write (6,'( /// 2(A,1PE12.3))') '   time =',t,
     $                                '   Surface temperature =',Temp(0)
         write (6,'(3X,A,11X,A,4X,A,6X,A)') 
     $         'x','Aluminum','Oxygen','AlOx'
         do i=0,NDIM
            write (6,'(1P8E12.3,I5)') x(i),Al(i),Ox(i),AlOx(i)
         end do
C
C  Normal stop
C
      stop 
      end


      subroutine GetTemp(x,Temp,t,NDIM)
C
C  Get values for the temperature
C
      implicit none
C
C  Declare variables
C
      real*8 x(0:NDIM)         ! Depth
      real*16 Temp(0:NDIM)     ! Temperature - K
      real*16 F                ! Uniform Temperature - degF
      real*8 t                 ! Time
      integer NDIM             ! No. of points
C
C  Read the temperature
C
      open (unit=222,file='SCX_TemperatureHistory.inp',status='old')
	  read (222,*) F           ! Read temperature in deg F
	  close(unit=222)
C
C  Convert to K
C
      Temp=5./9.*(F+459.6)     ! K
C
C  Normal return
C
      return
      end


      Subroutine GetDiffusion(Diff1,Diff2,Temp,D1,D2,Q0,T0,RT,NDIM)
C
C  Get the diffusion coefficient
C
      implicit none
C
C Declare variables
C
      real*8 Diff1(0:NDIM)     ! Diffusion constant
      real*8 Diff2(0:NDIM)     ! Diffusion constant
      real*16 Temp(0:NDIM)     ! Temperature
      real*8 D1                ! Reference Diffusion
      real*8 D2                ! Reference Diffusion
      real*8 Q0                ! Energy
      real*16 T0               ! Reference temperature
      real*8 RT                ! Gas Constant
      integer NDIM             ! Array lenght
C
C  Find diffusion constant
C
      Diff1=-Q0/RT*(T0-Temp)/(Temp*T0)
      Diff1=D1*exp(Diff1)
      Diff2=D2/D1*Diff1
C
C  Normal return
C
      return
      end


      Subroutine GetReaction(Reac,Temp,R0,DG,TG,RT,NDIM)
C
C  Get the diffusion coefficient
C
      implicit none
C
C Declare variables
C
      real*8 Reac(0:NDIM)      ! Reaction rate constant
      real*16 Temp(0:NDIM)     ! Temperature
      real*8 R0                ! Reference rate
      real*8 DG                ! Gibbs Free energy
      real*16 TG               ! Reference temperature
      real*8 RT                ! Gas Constant
      integer NDIM             ! Array lenght
C
C  Find diffusion constnat
C
      Reac=-DG/RT*(TG-Temp)/(Temp*TG)
      Reac=R0*exp(Reac)
C
C  Normal return
C
      return
      end


      real*8 function signum(x)
C                                           /  1 if x>0
C  Calculate the signum function signum(x)=<   0 if x=0
C                                           \ -1 if x<0
      real*8 x        ! Argument
      if (x.ne.0.) then
         signum=x/abs(x)
      else
         signum=0.
      end if
C
C  Normal return
C
      return
      end



      real*8 Function interpolate(xa, ya, x, n)
C
C  Interpolate using x to find the corresponding y
C  from a table of x values in xa and y values in ya
C
      implicit none
C
C  Declare variables
C

      integer i                          ! Counter
      integer n                          ! Number of values
      real*8 xa(*), ya(*)                ! Table values
	  real*8 x                           ! Independent variable
C
C  Check to see if x is between table ends or beyond and set appropriately
C
      if (x.le.xa(1)) then               ! x is beyond the lowest table entry
         interpolate = ya(1)
      else if (x.ge.xa(n)) then          ! x is beyond the highest table entry
         interpolate = ya(n)
      else                               ! x is between the end points
         i = 1
         do while (xa(i+1).lt.x) ! Look for table entry just below x
            i = i + 1
         end do
C
C  Interpolate
C
         interpolate =
     $      ya(i+1) * (x - xa(i)) / (xa(i+1) - xa(i))+
     $      ya(i) * (xa(i+1) - x) / (xa(i+1) - xa(i))
      end if
C
C  Normal return
C      return
       end


