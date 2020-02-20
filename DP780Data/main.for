      include 'modules.for'
      
      program main

      use mvoigt, only : initialize_voigt,invten,profac
      use meshelby, only : eshelbyb
      use sc_estimate, only : av_modulus,sc_new
      use miscellaneous_sub, only : read_main,reorient_grain
      use diffract, only : dif_planes
      use output
      use state_var
      use back_stress
      
      use const
      use flags
      use files
      use mphase_props 
      use grain_props
      use sample_props_v
      use mphase_state
      use grain_state
      use sample_state 
      use mphase_rate
      use grain_rate
      use sample_rate_v
      use bc
      use bc_v
      use solve_t
      use hard_law1, only : hl_read1=>hl_read,hl_ini_tau1=>hl_ini_tau
     #    ,hl_update_statv1=>hl_update_statv,hl_twin_mfp1=>hl_twin_mfp
     #    ,hl_pts_damp1=>hl_pts_damp
      use hard_law2, only : hl_read2=>hl_read,hl_ini_tau2=>hl_ini_tau
     #    ,hl_update_statv2=>hl_update_statv
      use back_stress_v, only : iphBsCtrl
      use phase_transf
      use twinning
      
      parameter (NSTATV=1000000)
      
      real :: aux3333(3,3,3,3),aux3(3),aux6(6),ass2_old(6,6)
     #   ,error_mod_org(2),a_guess(2),e2(6,6),spin_cor_s(3,3)
     #   ,etss_ini(6),stss_ini(6),STATEV(NSTATV)
      character(len=150) :: label
      
      type (sample) :: s  
      
    3 FORMAT(1h ,78('*'))
    5 FORMAT(1h ,'START THE PROCESS: ',a)
    6 FORMAT(1h ,'STEP:',i4,'/',i4,'    PROC:',i4,'/',i4)
     
      !set hard coded parameters
      a_guess=2.0/3.0
      iDeTwOpt=1
      iDtwMfp=0

      !read main input file
      call read_main(nph,axis,eulerph,nproc,nproc_cycle,ncycle,filecrys
     #    ,filesamp,fileproc,fileprev,filediff,filetemp,error_mod
     #    ,itmax_mod,itmax_grain,label,xmmin_in)
      error_mod_org=error_mod
      
      !open output files and write header
      if (iOutput.eq.1) call open_file

      !initializes voigt variables
      call initialize_voigt

      !initializes points and weights for Gaussian integration of Eshelby
      call eshelbyb (aux3,aux3333,0.,aux3333
     #    ,aux3333,0)           

      !add phase and grain properties from data in SX file and texture file
      do iph=1,nph
          call data_crystal(filecrys(iph),icrysym(iph),iph)
          if (kCL.eq.1) call hl_read1(iph)
          if (kCL.eq.2) call hl_read2(iph)
          if (iBackStress.eq.1) call bs_read(iph)
          call add_grains(iph,filesamp(iph))
          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2=SUM(nphngr(1:iph))
          wgt_ph(iph)=SUM(wgt(ng1:ng2)) 
      enddo
      !renormalize weights
      ngrain=SUM(nphngr(1:nph))
      totwgt=0.0
      do ng=1,ngrain
          totwgt=totwgt+wgt(ng)
      enddo
      do ng=1,ngrain
          wgt(ng)=wgt(ng)/totwgt
      enddo    
      do iph=1,nph
          wgt_ph(iph)=wgt_ph(iph)/totwgt      
      enddo
      !if only one crystal is present
      if(ngrain.eq.1)iSingleCry=1

      !initialize weight of grains for twinning
      if (itwinning.eq.1) then
          call initialize_tw(ngrain)
      endif
      
      !phtr initiate arrays
      if (iPhTr.eq.1) call phtr_read(filetemp(1))
      
      !read diffraction file and calculate initial lattice strain
      if(i_diff_dir.eq.1) then
          do iph=1,nph
              call dif_planes(filediff(iph),0.0,istep,0,iph)
          enddo
      endif

      !initialize shape
      do iph=1,nph
          call initialize_phase_shape(iph)
          if (ishape.ge.2) then
              ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
              ng2=SUM(nphngr(1:iph))
              call initialize_grain_shape(iph)
              call update_grain_shape(ng1,ng2)
              if (iPhTr.eq.1) call initialize_fijgr_pt(iph)
          endif
          !initiate array axisgr
          call update_phase_shape(iph)
      enddo

      !find Hill estimate for sample elastic stiffness and complience
      call av_modulus(1,1,ngrain,wgt,scs2,ccs2,ass2,css2,sss2)

      !set tengent stiffness to elastic stiffness for grains
      acs2(:,:,1:ngrain)=ccs2(:,:,1:ngrain)

      !find sc estimate for elastic stiffness 
      interaction=1 !(for linear elasticity interaction=1 is corect)
      call sc_new(0,liter,iverify,e2,interaction,0,1,ngrain
     #     ,etrss,alfass,ass2,1,0,axisgr,axisph
     #     ,css2,einvsa,einvsagr,escr4,escr4gr,aef
     #     ,aefgr,acs2,alfacs,aloc2,itmax_mod,ngrain,wgt,a_guess
     #     ,error_mod,meffc,etrcs_eig,etrss_eig,xmmin_in)  
     
      !find complience
      call invten (ass2,sss2)

      !previous procedure
      if (i_prev_proc.ge.1) then
          !initialize constants for hardening law for each phase
          do iph=1,nph
              ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
              if (kCL.eq.1) call hl_ini_tau1(ng1,tau,1)
          enddo
          if (kCL.eq.2) call hl_ini_tau2(1,tau,1)
          call read_statv(s,fileprev,NSTATV,STATEV)
          call get(s)
          !initialize twinning variables
          if(itwinning.eq.1) then 
              !find mfp
              call twin_barriers
              !set ptsdamp_old
              PTSdamp_old(:,1:ngrain)=PTSdamp(:,1:ngrain)
          endif 
          if (i_prev_proc.eq.2) then
              !error message
              if(ngParent.ne.ngrain)then
                  write(*,*)'ERROR: cannot evaluate CRSS in twins 
     #created in previous procedure'
                  pause
                  stop
              endif
          endif
          do iph=1,nph
              ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
              ng2=SUM(nphngr(1:iph))
              !update crystal properties 
              if (ng1.le.ng2) then
                  call cr_to_sa(ng1,ng2,0)
                  call cr_to_sa(ng1,ng2,1)    
              endif 
          enddo         
      endif
      
      !write header to output file
      if(iOutput.eq.1) call write_header(label)
      
      !set up initial values and constants
      step=1.0
      time_acum=0.0     
      
      !loop over processes
      iproc=0 
      icycle=1
      do while(iproc.lt.nproc)         
          iproc=iproc+1
          call cpu_time(start_time)      
          
          !write to terminal
          write(*,3)
          write(*,3)
          write(*,5) fileproc(iproc)
          
          !add bc from process file
          if (ivarBC.eq.0) then
              call add_bc(fileproc(iproc),i_temp_cij,i_ref_et,i_ref_st
     #        ,i_bc_mode,etss,stss,etssref)
          elseif (ivarBC.eq.1) then
              call add_bc_var(fileproc(1),iproc,i_temp_cij,i_ref_et
     #        ,i_ref_st,i_bc_mode,etss,stss,etssref)
          endif

          !initializes critical stress and accumulated shear in each grain.
          if (i_prev_proc.ne.1.and.iproc.eq.1.and.icycle.eq.1) then
              do ng=1,ngrain
                  if(i_prev_proc.ne.3)then
                      if (kCL.eq.1) call hl_ini_tau1(ng,tau,i_prev_proc)
                      if (kCL.eq.2) call hl_ini_tau2(ng,tau,i_prev_proc)
                  endif
                  if(i_prev_proc.eq.3)then
                      if (kCL.eq.1) call hl_ini_tau1(ng,tau,0)
                      if (kCL.eq.2) call hl_ini_tau2(ng,tau,0)
                  endif
              enddo
          endif
          
          !define initial temperature
          temp=temp_s
          
          !set up initial state (used to define bc)
          etss_ini=etss
          stss_ini=stss
          
          !write initial state
          if(iOutput.eq.1.and.ivarBC.eq.0) 
     #        call write_file(istep,iproc,step,temp)
     
          !loop over steps of process
          istep=0
          do while (istep.lt.nsteps)  
              !correct bc in case of rotations and volume change
              if(itwinning.eq.1.or.irot.ge.0.or.iPhTr.eq.1) then
                  call correct_bcr(i_bc_mode,istep,etss_ini,etss
     #                           ,stss_ini,stss)
              endif
              
              istep=istep+1  
              iTotStep=iTotStep+1
              
              if (i_temp_cij.eq.1) then
                  call temp_dep(temp)
                  call cr_to_sa(1,ngrain,1)
                  call cr_to_sa(1,ngrain,0)
              endif
              
              !write to terminal
              write(*,3)
              write(*,6) istep,nsteps,iproc,nproc
              
              !solve for props at t
              idiv=1  
              do while (idiv.eq.1)
                  idiv=0    
                  ass2_old=ass2  
                  
                  !define range of grains
                  ng1=1
                  ng2=ngrain 
                  
                  error_max = 0.0
                  error_max_out = 0.0
                  
                  !solves for rates at time t     
                  call main_calc(iproc,istep,ng1,ng2,idiv ,iact,nact
     #                ,jran,itmax_mod,ngrain,ng_update,nout_old,ietbc
     #                ,istbc,etrss,strss,alfass,ass2,tau,tau_bcst,gamtot
     #                ,gamd,ccs2,mcs,nmcs,stcs,alfacs,deltemp,e2,aef
     #                ,aefgr,acs2,aloc2,taud,taud_bcst,etrcs,strcs
     #                ,axisph,axisgr,css2,wgt,a_guess,error_mod,einvsa
     #                ,einvsagr,escr4,escr4gr,meffc,tau_update,etrbc
     #                ,strbc,etrcs_eig,etrss_eig,xmmin_in) 
                  
                  !increase tolerance if convergence is not reached
                  if (idiv.eq.1) call inc_toler(1,i_converge
     #                ,error_mod,error_mod_org,ass2,ass2_old,a_guess)
                  if (idiv.eq.0) error_mod(1) = error_mod_org(1)
                  if (idiv.eq.0) a_guess=2.0/3.0
              enddo          
              
              !update total variables
              !update gamtot and tau
              call update_grain_gam(step,gamd)
              call update_grain_tau(step,taud)
              
              !updates hard vars
              do ng=1,ngrain
                  if (kCL.eq.1)call hl_update_statv1(ng,step)
                  if (kCL.eq.2)call hl_update_statv2(ng)
                  if(iBackStress.eq.1.and
     #              .iphBsCtrl(ngrnph(ng)).ne.-1)then
                      call bs_update_statv(ng)
                  endif
              enddo
              
              if (iBackStress.eq.1.and.
     #              any(mask=iphBsCtrl(1:nph).eq.-1))then
                  call ph_bcst
              endif
              !update phase transf
              if(iPhTr.eq.1)then              
                  !find escgr for each slip system and grain
                  call SFW_calc(1,nphngr(1))
                  !finde lode parameter and stress triaxiality for grains
                  call lode_triax(1,nphngr(1))
                  !calculate dwgt array, create new grains and update ngrain
                  call wgtd_calc(wgtd)
                  !update weights
                  do ng=ng1,ngrain
                      wgt(ng)=wgt(ng)+wgtd(ng)
                  enddo
                  do iph=1,nph
                      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
                      ng2=SUM(nphngr(1:iph))
                      wgt_ph(iph)=SUM(wgt(ng1:ng2))
                  enddo
                  !update phase transformation strain (etrcs_eig)
                  call eig_calc
              endif
              
              !update twinning
              if (itwinning.eq.1) then
                  !find mfp
                  call twin_barriers
                  !find twin fraction
                  call twin_fraction(step)
                  !find wgtd
                  call wgtd_calc_tw(wgtd)
                  !update weights
                  do ng=ng1,ngrain
                      wgt(ng)=wgt(ng)+wgtd(ng)
                  enddo
                  !correct parent grain strss
                  call corr_p_st(stcs)
                  !change CRSS in twins due to mfp
                  if (ngParent.ne.ngrain) then
                      if (kCL.eq.1) call hl_twin_mfp1(tau)
                  endif
                  !rise CRSS of twin sys if volume fraction is large
                  if (kCL.eq.1.or.kCL.eq.1) call hl_pts_damp1
     #                (iproc,istep,tau)   
              endif
              
              !update shape of inclusion (array AXISPH or AXISGR)
              if(ishape.eq.1) then
                  call update_phase_etrss(etrcs,etrss_ph) !MZ_2ph
                  do iph=1,nph    
                      call update_phase_fij(iph,step,etrss_ph) 
                      call update_phase_shape(iph) 
                  enddo
              elseif (ishape.ge.2) then
                  ! update deformation gradient (fijgr)
                  if (iPhTr.eq.1) call update_fijgr_pt !update phase transf def grad
                  call update_grain_fij(step,etrcs)
                  ! update ellipsoid shape (axisgr)
                  if (itwinning.eq.1) then
                      call update_grain_shape(1,ngParent)
                      do ng=ngParent+1,ngrain
                          ngp=iParentGrain(ng) ! identify parent grain
                          IPTS=iParentSystem(ng) !indentify parent system
                          call twin_shape(ng,ngp,IPTS,axisgr)
                      enddo
                  else
                      call update_grain_shape(1,ngrain)
                  endif
              endif      
              
              if (istep.lt.nsteps) then
                  !old update
                  IF (irot.EQ.1) THEN
                      !updates the texture (array r)
                      call update_grain_spin(step,omegabcr,etrss,ESCR4
     #                    ,EINVSA)
                      call update_grain_orientatin(step,ngrain,omegag)
                      !recalculates m,q,burg,alpha,ccs2,scs2
                      do iph=1,nph 
                          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
                          ng2=SUM(nphngr(1:iph))
                          if (ng1.le.ng2) then
                              call cr_to_sa(ng1,ng2,0)
                              call cr_to_sa(ng1,ng2,1)    
                          endif
                      enddo  
                      !updates stresses in sample and grains, to stay on yield surface  
                      if (istep.ne.0) then
                          call update_grain_corot_st(ngrain,step
     #                       ,spin_cor_s,omegag)
                          call update_sample_corot_st(omegabcr
     #                       ,spin_cor_s)
                      endif
                      if (iSingleCry.eq.1) then
                          stss=stcs(:,1)
                      endif             
                  ENDIF
                  !new update
c                  IF (irot.EQ.1) THEN
c                      !update omegag
c                      call update_grain_spin(step,omegabcr,etrss,ESCR4
c     #                    ,EINVSA)
c                      !update drotcs
c                      do ng=1,ngrain
c                          call REORIENT_GRAIN(drotcs(:,:,ng)
c     #                        ,omegag(:,:,ng))
c                      enddo
c                      !update variables 
c                      call update_grain_config(1,ngrain)
c                      !sample
c                      if (iSingleCry.eq.1) omegabcr=omegag(:,:,1)
c                      call REORIENT_GRAIN (drotss,omegabcr)
c                      call update_sample_config
c                  ENDIF
              end if
              
              !update the overall stresses and strains (tot and el)
              call update_sample_state(etrss,strss)
              
              !update grain stresses and strains
              call update_grain_state(step,deltemp,etrcs,strcs
     #            ,etrcs_eig)
	 
              !update phase state
              do iph=1,nph
                  call update_phase_state(iph)
              enddo
              
              !update temperature
              temp=temp+deltemp
              
              !calculate averages and deviations
              call g_average(etrss,strss,etrcs,strcs,etrav)
              
              !correct total stress and strain to be volume average
              if (itwinning.eq.1.or.iPhTr.eq.1) then !.or.irot.eq.1
                  do i=1,6
                      stss(i)=stav(i)
                      etss(i)=etav(i) 
                  enddo
                  !elastic strain
                  do i=1,6
                      etelss(i)=0.0
                      do j=1,6
                          etelss(i)=etelss(i)+sss2(i,j)*stss(j)
     #                        *profac(j)
                      enddo
                  enddo
              endif
              
              !write
              if (iOutput.eq.1) then 
                  call write_file(istep,iproc,step,temp)
                  call write_temp
              endif
          enddo
          !cycle_proc   when one cycle is completed, the counter gets 0, so the process files will be read for the next(icycle+1) cycle
          if(iproc.eq.nproc_cycle.and.icycle.lt.ncycle)then
              iproc=0
              icycle=icycle+1
          endif
      enddo
      
      !close output file
      if (iOutput.eq.1) call close_file
      
      !store state variables
      call set(s)      
      
      !write state variables
      if (i_prev_proc.eq.0) then
          call write_statv(s,fileprev,NSTATV,STATEV)
      endif
      
      !measure time
      call cpu_time (end_time)
      dum_time=end_time-start_time
      time_acum=time_acum+dum_time    
      write(*,'('' TOTAL TIME='',f8.2,'' secs'')') time_acum
      
      end program main