c      
c***********************************************************************
c *** const ************************************************************
c***********************************************************************
      module const
      
      PARAMETER (NGR= 10000)  !Number of grains.
      PARAMETER (NMOD=10)     !Number of deformation modes.
      PARAMETER (NSLS=100)     !Number of deformation systems.
      PARAMETER (NDIFFX=60)   !Number of difraction directions.
      PARAMETER (NPROCX=10000)    !Number of thermomechanical processes.
      PARAMETER (KPT=30)      !Number of diffraction output detail points.
      PARAMETER (NPHM=5)	!Number of phases MZ_2ph
      
      REAL(8), PARAMETER :: NaN = 
     #TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8)      
      
      end module const
c      
c***********************************************************************
c *** flags ************************************************************
c***********************************************************************
      module flags
      
      use const
      
      integer :: ishape,irot,ipileup,kSM,iPoleFigFlag,i_diff_dir,iDiag
     #    ,kCL,iSingleCry,iTwinLaw,i_prev_proc,iDetwOpt,iDtwMfp,ilatBS
     #    ,iBackStress,iPhTr,itwinning,iOutput,itexskip,nCoatedPh
     #    ,nCoatingPh,ivarBC,inonSch,iTwPh(NPHM),iUW,iUmat,iDiagBs
      
      end module flags
c      
c***********************************************************************
c *** files ************************************************************
c***********************************************************************
      module files
      
      use const
      
      character(len=150) :: filecrys(NPHM),filesamp(NPHM)
     #             ,fileproc(NPROCX),fileprev,filediff(NPHM),filetemp(5)
      
      end module files
c      
c***********************************************************************
c *** mphase_props_v ***************************************************
c***********************************************************************
      module mphase_props_v
      
      use const
      !crydatci
      integer :: nsm(NMOD,NPHM),itw(NMOD),nmodes(NPHM)
     #    ,nsys(NPHM),nslmod(NPHM),nslsys(NPHM),ntwmod(NPHM)
     #    ,ntwsys(NPHM),nphngr(NPHM),icrysym(NPHM),ISECTW(NMOD)
     #    ,ngrnph(NGR)
      !crydatcr
      real :: ccc2(6,6,NPHM),scc2(6,6,NPHM),alfacc(6,NPHM)
     #       ,qc2(3,3,nsls,NPHM),mc2(3,3,NSLS,NPHM),ccc2p0(6,6,NPHM)
     #       ,ccc2dp(6,6,NPHM),wgt_ph(NPHM),bcc(3,NSLS,NPHM)
     #       ,ncc(3,NSLS,NPHM),stw(NMOD,NPHM),TWTHRES(4,NMOD,NPHM)
     #       ,nmc2(3,3,NSLS,NPHM),nscc(4,NMOD,NPHM)
      !flow law
      integer :: mode_slip(NMOD,NSLS,NPHM),iSysMode(NSLS,NPHM)     
      !twini
      integer iTwinSys(NSLS) 
      
      end module mphase_props_v
c      
c***********************************************************************
c *** sample_props_v ***************************************************
c***********************************************************************   
      module sample_props_v
      
      use const
      !polycry
      real :: css2(6,6),sss2(6,6)
      !samdatai
      integer :: ngrain,nph,ngParent
      
      end module sample_props_v
c      
c***********************************************************************
c *** grain_props_v ******************************************************
c***********************************************************************
      module grain_props_v
      use const
      !crydatsr
      real :: ccs2(6,6,NGR),scs2(6,6,NGR),mcs(6,NSLS,NGR)
     #        ,alfacs(6,NGR),qcs(6,nsls,ngr)
     #        ,phi(NGR),the(NGR),ome(NGR),wgt(NGR),r(3,3,NGR)
     #        ,bcs(3,NSLS,NGR),ncs(3,NSLS,NGR),nmcs(6,NSLS,NGR)  
      
      end module grain_props_v      
c      
c***********************************************************************
c *** mphase_state_v ***************************************************
c***********************************************************************
      module mphase_state_v
      
      use const
      !total
      real :: stss_ph(6,NPHM),etss_ph(6,NPHM),etelss_ph(6,NPHM) 
     #        ,etss_all(6)
      !shape
      real :: axis(3,NPHM),axisph(0:3,3,NPHM),eulerph(3,NPHM)
     #        ,fijph(3,3,NPHM)      
      
      end module mphase_state_v
c
c***********************************************************************
c *** grain_state_v ****************************************************
c***********************************************************************
      module grain_state_v

      use const
      
      !total
      real :: stcs(6,NGR),etcs(6,NGR),etelcs(6,NGR),stcsref(6,NGR)
     #    ,etelhycs(6,NGR),tau(NSLS,NGR),gamtot(NGR),etthcs(6,ngr)
     #    ,etcs_eig(6,NGR)
      !shape
      real :: fijgr(3,3,NGR),axisgr(0:3,3,NGR),fijgr_pt(3,3,NGR)
      !logical
      integer :: nact(NGR),iact(NSLS,NGR),ng_update(NGR),nactSlip(NGR)
     #    ,iactSlip(NSLS,NGR)
      !twin/phase transf
c      integer :: iParentGrain(NGR),iParentSystem(NGR),iParentMode(NGR)
c     #    ,iChildGrain(NSLS,NGR)
      !other
      real :: tau_update(NSLS,NGR),actwgt(0:20)   
      
      end module grain_state_v  
c      
c***********************************************************************
c *** mphase_rate ******************************************************
c***********************************************************************
      module mphase_rate_v
      
      use const
      !rate
      real :: etrss_ph(6,NPHM),strss_ph(6,NPHM),etelrss_ph(6,NPHM)   
      !eshelby
      real :: ESCR4(3,3,3,3,NPHM),EINVSA(3,3,3,3,NPHM)
     #        ,aef(6,6,NPHM),auxsample(6,NPHM)       
      end module mphase_rate_v
c
c***********************************************************************
c *** grain_rate_v *****************************************************
c***********************************************************************
      module grain_rate_v
      
      use const
      
      !rate
      real :: strcs(6,NGR),etrcs(6,NGR),omegag(3,3,ngr),acs2(6,6,NGR)
     #        ,gamd(NSLS,NGR),taud(NSLS,NGR),taud_bcst(NSLS,NGR,2)
     #        ,wgtd(NGR),etrcs_eig(6,NGR),drotcs(3,3,NGR)
      !eshelby
      real :: EINVSAGR(3,3,3,3,NGR),ESCR4GR(3,3,3,3,NGR),aloc2(6,6,ngr)
     #        ,meffc(NGR),as(3,3,3,3,ngr),aefgr(6,6,NGR)      
      
      end module grain_rate_v
c      
c***********************************************************************
c *** sample_state_v ***************************************************
c***********************************************************************
      module sample_state_v
      !total
      real :: stss(6),etss(6),etelss(6),etss_eig(6)
      !other
      real :: etssref(6),stav(6),etav(6),actav       
      !seed for random number generator
      integer :: jran
      
      end module sample_state_v
c      
c***********************************************************************
c *** sample_rate_v ****************************************************
c***********************************************************************
      module sample_rate_v
      
      !rate
      real :: strss(6),etrss(6),ass2(6,6),alfass(6),strav(6),etrav(6)
     #    ,etrss_eig(6),drotss(3,3)
     
      end module sample_rate_v
c      
c***********************************************************************
c *** twinning_v *******************************************************
c***********************************************************************
      module twinning_v
      
      use const
      !twini
      integer :: ktwtag(NSLS,NGR),link(NGR)
     #               ,MaxTwins,IncludeTS
     #               ,iTwinLevel(NGR)
     #               ,iMaxTwinLevel,iParentGrain(NGR)
     #               ,iParentSystem(NGR),iParentMode(NGR)
     #               ,iChildGrain(NSLS,NGR)
      !twinr
      real :: vfrac_mod_acum(NMOD),vfrac_mod(NMOD)
     #               ,offset(6),TwinFrac(NMOD,NPHM),TwinCRSS(NMOD,NPHM)
     #               ,TVF,CTVF,wgtx(NGR)     
      !twin CG
      integer :: IPTSGR(NSLS,NGR),IPTSGRC(NGR)
      real :: sysmfp(NSLS,NSLS),tau_tw(NSLS),TwFrSy(NSLS,NGR)
     #    ,WGTCG_all(NGR),TWFRGR(NMOD,NGR),PTS_THRES(NMOD,NPHM)
     #    ,TWFR_THRES(NMOD,NPHM),PTSDAMP(NSLS,NGR),PTSDAMP_OLD(NSLS,NGR)
     #    ,PTVFM(NMOD),CTVFM(NMOD),newGR(NGR),newGRSh(NGR),axis3(NGR)
      
      !state variables definition
      type state_tw
          real :: TwFrSy(NSLS),wgtcg_all,iParentGrain,iParentMode
     #        ,iParentSystem,iPTSgr(NSLS),iPTSgrc,iChildGrain(NSLS)
     #        ,PTSdamp(NSLS),newGR,newGRSh,axis3 !,tau0_tw(NSLS)
      end type state_tw   
      
      end module twinning_v
c      
c***********************************************************************
c *** hard_law1_v *******************************************************
c***********************************************************************
      module hard_law1_v
      
      use const

      integer :: iDirSys(NSLS,NPHM),iOpposite(NSLS),iact_sys(NSLS,NGR)
      
      !flags 
      integer :: iFCC,iGrShCRSS,irevlaw(NPHM)
      !hard law state variables
      real :: rho_rev(NSLS,NGR),rho_forw(NSLS,NGR)
     #  ,rho_tot(NSLS,NGR),rho_tot_max(NSLS,NGR),rho_deb(NGR)
     #  ,rho_act(NSLS,NGR),tau_act(NSLS,NGR),gam_acc(NSLS,NGR)
     #  ,rssmin(NSLS,NGR)
      !hard law constants (reading)
      real :: chi_inter(NPHM),q_rate(NPHM)
     #  ,grsze,BURG(NMOD,NPHM),ACTENER(NMOD,NPHM),aK1(NMOD,NPHM)
     #  ,DRAG(NMOD,NPHM),edot_zero(NMOD,NPHM),rho_ini_for(NMOD,NPHM)
     #  ,rho_ini_deb(NMOD,NPHM),tau0_mode(NMOD,NPHM)
     #  ,tau0_mode_a(NMOD,NPHM),tau0_mode_b(NMOD,NPHM)
     #  ,tau0_mode_c(NMOD,NPHM),TLATENT(NMOD,NMOD,NPHM)
     #  ,TLATENT1(NMOD,NMOD,NPHM),HPfac(NMOD,NPHM),HPK0(NMOD,NPHM)
     #  ,HPKCG(NMOD,NMOD,NPHM),tauHP(NMOD,NPHM)
     #  ,a_deb(NMOD,NPHM),a_deb_a(NMOD,NPHM),a_deb_b(NMOD,NPHM)
     #  ,a_deb_c(NMOD,NPHM),p_rev(NPHM),aM_par(NPHM)
     #  ,d_mod(NMOD,NMOD+1,NPHM),g_mod(NMOD,NMOD+1,NPHM),d0,d1,d2,d3,d4
     #  ,d5,g0,g1,g2,g3,g4,g5,rev_coeff(2,NPHM)
     #  ,tau_crit(NMOD,NPHM),tau_crit_a(NMOD,NPHM)
     #  ,tau_crit_b(NMOD,NPHM),tau_crit_c(NMOD,NPHM)
     #  ,tau_prop(NMOD,NPHM),tau_prop_a(NMOD,NPHM)
     #  ,tau_prop_b(NMOD,NPHM),tau_prop_c(NMOD,NPHM)
     #  ,tau_crit_de(NMOD,NPHM)
     #  ,tau_crit_de_a(NMOD,NPHM),tau_crit_de_b(NMOD,NPHM)
     #  ,tau_crit_de_c(NMOD,NPHM),tau_prop_de(NMOD,NPHM)
     #  ,tau_prop_de_a(NMOD,NPHM),tau_prop_de_b(NMOD,NPHM)
     #  ,tau_prop_de_c(NMOD,NPHM),alpha_reg,rev_par(NPHM)
      !hard law constants (derived)
      real :: alatent(NSLS,NSLS,NPHM),alatentDD(NSLS,NSLS,NPHM)
     #  ,shearmod(NMOD,NPHM),xmfp(NSLS,NGR),Boltz
      !temp arrays
      integer :: itemp_act_sys(NSLS,NGR)
      real :: temp_rho_tot_0(NSLS,NGR),rho_tot_0(NSLS,NGR)
     #  ,Drho_deb_Dgamma(NSLS),aK2(NMOD,NPHM)
     #  ,Drho_forw_Dgamma(NSLS),Drho_rev_Dgamma(NSLS)
     #  ,Drho_tot_Dgamma(NSLS),Drho_tot_Dgamma_hd(NSLS)
      
      !state variables definition
      type state_hl
          real :: tau(NSLS),rho_rev(NSLS),rho_forw(NSLS),rho_tot(NSLS)
     #     ,rho_tot_max(NSLS),rho_act(NSLS),tau_act(NSLS),iact_sys(NSLS)
     #     ,itemp_act_sys(NSLS),rho_deb,rssmin(NSLS)
     #     ,temp_rho_tot_0(NSLS)
      end type state_hl      
      
      end module hard_law1_v
c      
c***********************************************************************
c *** back_stress_v ****************************************************
c***********************************************************************
      module back_stress_v
      
      use const
      
      integer :: ibs, iphBsCtrl(NPHM)
      
      real :: tau_bcst(NSLS,NGR,2),aL_bs(NSLS,NSLS,NPHM)
     #    ,tau_sat(NMOD,NPHM),ni(NMOD,NPHM),gam_b(NMOD,NPHM)
     #    ,fact_1(NMOD,NPHM),vol_frac,taud_bcst_sys(NSLS,NGR)
     #    ,tau_bcst_tot(NGR),alpha_reg_bs
      
      type state_bs
          real :: gam_acc(NSLS),tau_bcst(NSLS,2)
      end type state_bs
      
      end module back_stress_v
c      
c***********************************************************************
c *** phase_transf_v ***************************************************
c***********************************************************************
      module phase_transf_v
      
      use const
      
      real :: etcs_pt(6,NGR),alpha0,alphaK,beta0,betaK,anexp,wgtcr
     #    ,triax(NGR),escgr(NGR),sfe,gamtot_phtr(NGR),hn_a(3),e1_a(3)
     #    ,e2_a(3),rsm_e(3,3,NGR),O_cubic_sym(3,3,24)
      
      type state_pt
          real :: etrcs_eig(6),etcs_pt(6),iChildGrain,iParentGrain
     #        ,gamtot
      end type state_pt
      
      end module phase_transf_v
c      
c***********************************************************************
c *** bc_v *************************************************************
c***********************************************************************
      module bc_v
      
      !prodatai
      integer :: nsteps,nproc,itmax_mod,itmax_grain,istbc(6),ietbc(6)
     #    ,iTotStep,nproc_cycle,ncycle,nsteps_min
     
      !prodatar
      real :: temp_s,temp_f,deltemp,stbc(6),strbc(6),etbc(6)
     #               ,etrbc(6),error_mod(2),ctrl_incr,edot_macro
     #               ,error_max,thres_et

      !JAW
      real :: etbc_sym(3,3),omegabc(3,3),omegabcr(3,3)
     #               ,fulletbc(3,3)
      integer :: ifulletbc(3,3),icvx,i_control_var
      
      !varBC
      integer :: ietbc_var(6),istbc_var(6)
      real :: etbc_var(6,100000),stbc_var(6,100000)
      
      end module bc_v
c      
c***********************************************************************
c *** diffract_v *******************************************************
c***********************************************************************
      module diffract_v
      
      use const
      
      !diffract
      real :: wgtset(NDIFFX),wgtsetini(NDIFFX),wgtgrset(NDIFFX,ngr)
     #    ,RAND_WGT(NDIFFX),vs(NDIFFX,3,NPHM),vc(ndiffx,3,24,NPHM)
     #    ,DEC(3,3,NDIFFX),als(NDIFFX)
      integer :: ngrset(NDIFFX),igrset(NDIFFX,NGR),ndiff(NPHM)
     #    ,nfamily(NDIFFX,NPHM)
      real :: toler
      
      end module diffract_v
c      
c***********************************************************************
c *** output_v *********************************************************
c***********************************************************************
      module output_v
      
      use const
      
      real :: etsseq,etsspleq,etsseigeq,stsseq,etrsseq
     #           ,etrsspleq,strsseq,volume,pressure,wtotal,wplastic
     #           ,stet,stetav
      
      real :: TEMP(10000)
      
      end module output_v
c      
c***********************************************************************
c *** epsc_var *********************************************************
c***********************************************************************
      module epsc_var

      use const
      use flags  
      use twinning_v
      use mphase_props_v
      use mphase_state_v
      use mphase_rate_v
      use sample_props_v
      use sample_state_v
      use sample_rate_v
      use grain_props_v
      use grain_state_v
      use grain_rate_v
      use diffract_v
      use bc_v
      use hard_law1_v
c      use mvoigt
      
      !ResPole
      real :: chiPoleFig(NDIFFX),etaPoleFig(NDIFFX)
      integer :: nPoleFig
      
      integer :: jact(0:NSLS)
     #               ,iwrite9
      
     
      !split 
      real :: ass2_old(6,6),error_mod_org(2),a_guess(2)
      integer :: iseed
      
      end module epsc_var
c      
c***********************************************************************
c *** miscellaneous_sub ************************************************
c***********************************************************************
      module miscellaneous_sub
      !jacobi :
      !eigsrt :
      !det :Calculates determinant of 3x3 matrix      
      !lubksbc :Solves the set of N linear equations
      !ludcmpc :Performs the LU decomposition of a matrix 
      !tmismatch :Calculates scalar difference between tensors
      !tnorm :Calculates norm of a tensor
      !ran2 :Generates random number given seed
      !euler :Calc. transform matrix for each grain 
      !reorient_grain :Calc. rotation matrix from incrmental spin
      !rodrigues :Calc. rotation matrix from incrmental spin
      contains

      subroutine read_main(nph,axis,eulerph,nproc,nproc_cycle,ncycle
     #    ,filecrys,filesamp,fileproc,fileprev,filediff,filetemp
     #    ,error_mod,itmax_mod,itmax_grain,label,xmmin_in)
      use const
      use flags
      
    1 FORMAT(a)      
    2 FORMAT(1h ,a)
    3 FORMAT(1h ,78('*'))
  100 FORMAT(1h ,14('*'),' SELF-CONSISTENT THERMO-ELASTOPLASTIC CODE'
     #,' "EPSC" ',14('*'))
      !out
      integer, intent(out) :: nph,nproc
      real, intent(out) :: axis(3,NPHM),eulerph(3,NPHM),error_mod(2)
     #    ,xmmin_in
      character(len=150), intent(out) :: filecrys(NPHM),filesamp(NPHM)
     #             ,fileproc(NPROCX),fileprev,filediff(NPHM),filetemp(5)
     #             ,label
     
      CHARACTER*78 prosa

      !Open and read main control file: "EPSC4.IN"   -> unit=1
      OPEN(unit=1,file='epsc4.in',status='old')       
      
      read(1,1) label          ! simulation label
      write(*,3)
      write(*,100)
      write(*,2) label
      write(*,3)
      READ(1,*) nph           
      READ(1,1) prosa
      READ(1,*) ishape          
      do iph=1,nph !MZ_2ph          
          READ(1,*) (axis(i,iph),i=1,3)  !MZ_pseudo read axis and eulerph for each phase
          READ(1,*) (eulerph(i,iph),i=1,3) !initial ellipsoid orientation angles
      enddo
      !Name of TEXTURE file
      read(1,1) prosa
      do iph=1,nph !MZ_2ph       
          read(1,1) filesamp(iph)
      enddo !MZ_2ph 
      READ(1,*) irot            
      !Name of MATERIAL file
      read(1,1) prosa
      do iph=1,nph !MZ_2ph      
          read(1,1) filecrys(iph) 
      enddo !MZ_2ph    
      !Hardening law, backstress MZ_bs      
      read(1,1) prosa !MZ_bs
      read(1,*) kCL,iTwinLaw,iBackStress,iPhTr,iOutput,nCoatedPh
     #    ,nCoatingPh,ivarBC,inonSch,xmmin_in
      !Precision Settings
      read(1,1) prosa     
      read(1,*) itmax_mod
      READ(1,*) error_mod(1)
      read(1,*) itmax_grain
      !Flag for previous procedure and the file name 
      read(1,1) prosa
      read(1,*) i_prev_proc
      read(1,1) fileprev
      READ(1,*) itexskip
      !Flag for difraction calculation and the file with diff directions
      read(1,*) i_diff_dir
      do iph=1,nph   
          read(1,1) filediff(iph)
      enddo
      !Flag for strain pole figure calculation
      read(1,*) i_strpf
      !Number of thermomechanical process in this simulation 
      read(1,1) prosa
      read(1,*) nproc,nproc_cycle,ncycle
      if (nproc.gt.NPROCX) then
        write(*,'(1h ,''ERROR: Number of processes greater''
     #  ,'' than code dimension !!!'',/,1h ,''DIMENSION in code = ''
     #  ,i3)') NPROCX
        write(*,*)
        write(*,'(1h ,''STOP IN MAIN '')')
        stop
      endif
      !Names of PROCESS files 
      read(1,1) prosa
      do n=1,nproc
        read(1,1) fileproc(n)
      enddo
      read(1,1) prosa    
      read(1,1) filetemp(1)
      
      CLOSE(unit=1)         
      
      end subroutine read_main
      
      SUBROUTINE jacobi(a,n,np,d,v,nrot,ier)

      !Copied from VPSC6 by JN, 8-7-07
c **********************************************************************

      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do ip=1,n
        do iq=1,n
          v(ip,iq)=0.
        enddo
        v(ip,ip)=1.
      enddo
      do ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
      enddo
      nrot=0
      do i=1,50
        sm=0.
        do ip=1,n-1
          do iq=ip+1,n
            sm=sm+abs(a(ip,iq))

          enddo
        enddo
c
        if(sm.eq.0.)then
        ier=0
        return
        endif
c
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do ip=1,n-1
          do iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     #        g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
              enddo
              do j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
              enddo
              do j=1,n
                g=v(j,ip)
                h=v(j,iq)

                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
              enddo
              nrot=nrot+1
            endif
          enddo
        enddo
        do ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
        enddo
      enddo
c      pause 'too many iterations in jacobi'
c
      ier=1
c
      return
      END SUBROUTINE

      SUBROUTINE eigsrt(d,v,n,np)

      !Copied from VPSC6 by JN, 8-7-07
c **********************************************************************

      INTEGER n,np
      REAL d(np),v(np,np)
      INTEGER i,j,k
      REAL p
      do i=1,n-1
        k=i
        p=d(i)
        do j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
        enddo
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
          enddo
        endif
      enddo
      return
      END SUBROUTINE

      function det(a)

      !Copied from VPSC6 by JN, 8-7-07
c **********************************************************************

      dimension a(3,3)

      det=a(1,1)*a(2,2)*a(3,3)
     #   +a(1,2)*a(2,3)*a(3,1)
     #   +a(1,3)*a(2,1)*a(3,2)
     #   -a(1,3)*a(2,2)*a(3,1)
     #   -a(2,3)*a(3,2)*a(1,1)
     #   -a(1,2)*a(2,1)*a(3,3)
      return
      end function

      SUBROUTINE LUDCMPC(A,N,NP,INDX,D)

      PARAMETER (TINY = 1.0e-20)
      PARAMETER (NMAX = 120)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      d=1.0
      do i=1,n
         aamax=0.0
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo
         if (aamax.eq.0.0) then
            PRINT *, 'PROGRAM STOP --> ZERO MATRIX IN SUBROUT LUDCMPC'
            STOP
         endif
         vv(i)=1.0/aamax
      enddo
      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
         enddo
         aamax=0.0
         do i=j,n
            sum=a(i,j)
            do k=1,j-1
               sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
            dum = vv(i)*abs(sum)
            if (dum.ge.aamax) then
               aamax=dum
               imax=i
            endif
         enddo
         if (j.ne.imax) then
            do k=1,n
               dum=a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
            enddo
            d = -d
            vv(imax) = vv(j)
         endif
         indx(j) = imax
         if (a(j,j).eq.0.0) a(j,j) = TINY
         if (j.ne.n) then
            dum = 1.0/(a(j,j))
            do i=j+1,n
               a(i,j) = a(i,j)*dum
            enddo
         endif
      enddo
      return
      end subroutine

      SUBROUTINE LUBKSBC(A,N,NP,INDX,D)

      DIMENSION A(NP,NP),INDX(N),D(N)

      ii=0
      do i=1,n
         ll=indx(i)
         sum=d(ll)
         d(ll) = d(i)
         if (ii.ne.0) THEN
            do j=ii,i-1
               sum = sum - a(i,j) * d(j)
            enddo
         else if (sum.ne.0.0) then
            ii=i
         endif
         d(i) = sum
      enddo
      do i=n,1,-1
         sum = d(i)
         if (i.lt.n) then
            do j=i+1,n
               sum = sum - a(i,j) * d(j)
            enddo
         endif
         d(i) = sum/a(i,i)
      enddo
      return
      end SUBROUTINE

      FUNCTION tmismatch (v1,v2,nrows,ncols)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TMISMATCH   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES RELATIVE DIFFERENCE BETWEEN TWO NROWSxNCOLS MATRICES
C     THE DIFFERENCE IS RELATIVE TO THE NORM OF THE ARITHMETIC AVERAGE
C     OF BOTH DATA.
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      DIMENSION v1(36),v2(36)
      DIMENSION v_dif(36),v_ave(36)

      do i=1,nrows*ncols
        v_dif(i)=v1(i)-v2(i)
        v_ave(i)=0.5*(v1(i)+v2(i))
      enddo
      tmismatch=tnorm(v_dif,nrows,ncols)/tnorm(v_ave,nrows,ncols)

      RETURN
      END FUNCTION

      FUNCTION tnorm(v,nrows,ncols)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     FUNCTION TNORM   ---->   VERSION OF 27/DEC/98
C
C     CALCULATES THE NORM OF A NROWSxNCOLS-MATRIX (NROWS,NCOLS =< 6)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION v(36)

      tnorm=0.0
      do i=1,nrows*ncols
        tnorm=tnorm+v(i)*v(i)
      enddo
      tnorm=sqrt(tnorm)

      RETURN
      END FUNCTION

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END FUNCTION

      subroutine euler(iopt,ph,th,tm,a)
c
c     CALCULATE THE EULER ANGLES ASSOCIATED WITH THE TRANSFORMATION
c     MATRIX A(I,J) IF IOPT=1 AND VICEVERSA IF IOPT=2
c     A(i,j) TRANSFORMS FROM SYSTEM sa TO SYSTEM ca.
c     ph,th,om ARE THE EULER ANGLES (in degrees) OF ca REFERRED TO sa.
c *****************************************************************************

      dimension a(3,3)
      pi=4.*atan(1.0)

      if(iopt.eq.1) then
        th=acos(a(3,3))
        if(abs(a(3,3)).ge.0.9999) then
          tm=0.
          ph=atan2(a(1,2),a(1,1))
        else
          sth=sin(th)
          tm=atan2(a(1,3)/sth,a(2,3)/sth)
          ph=atan2(a(3,1)/sth,-a(3,2)/sth)
        endif
        th=th*180./pi
        ph=ph*180./pi
        tm=tm*180./pi
      else if(iopt.eq.2) then
        sph=sin(ph*pi/180.)
        cph=cos(ph*pi/180.)
        sth=sin(th*pi/180.)
        cth=cos(th*pi/180.)
        stm=sin(tm*pi/180.)
        ctm=cos(tm*pi/180.)
        a(1,1)=ctm*cph-sph*stm*cth
        a(2,1)=-stm*cph-sph*ctm*cth
        a(3,1)=sph*sth
        a(1,2)=ctm*sph+cph*stm*cth
        a(2,2)=-sph*stm+cph*ctm*cth
        a(3,2)=-sth*cph
        a(1,3)=sth*stm
        a(2,3)=ctm*sth
        a(3,3)=cth
      endif

      return
      end subroutine

      SUBROUTINE REORIENT_GRAIN (AROT,C)
      
C *****************************************************************************
C     SUBROUTINE REORIENT_GRAIN  (ex-ORIENT)   --->   VERSION 19/01/01
c                                     (copied from VPSC)
C     BUILDS INCREMENTAL ROTATION MATRIX 'AROT' BASED ON RODRIGUES FORMULA.
C     'C' IS THE INCREMENTAL LATTICE SPIN.
C     'AROT' TRANSFORMS FROM INITIAL TO FINAL ORIENTATION.
C *****************************************************************************

      dimension c(3,3),th(3,3),th2(3,3),v(3),vbar(3),arot(3,3)

      v(1)=c(3,2)
      v(2)=c(1,3)
      v(3)=c(2,1)
      snorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      snorm1=tan(snorm/2.)
      IF(SNORM.LT.1.e-06) SNORM=1
      do 20 i=1,3
      vbar(i)=snorm1*v(i)/snorm
20    continue
      snorm=vbar(1)*vbar(1)+vbar(2)*vbar(2)+vbar(3)*vbar(3)
      th(3,2)= vbar(1)
      th(1,3)= vbar(2)
      th(2,1)= vbar(3)
      th(2,3)=-vbar(1)
      th(3,1)=-vbar(2)
      th(1,2)=-vbar(3)
      do 40 i=1,3
40    th(i,i)=0.
      do 30 i=1,3
      do 30 j=1,3
      th2(i,j)=0.
      do 50 k=1,3
50    th2(i,j)=th2(i,j)+th(i,k)*th(k,j)
30    continue
      do 60 i=1,3
      do 60 j=1,3
60    arot(i,j)=(i/j)*(j/i)+2.*(th(i,j)+th2(i,j))/(1.+snorm)

      return
      end subroutine

      SUBROUTINE RODRIGUES (C,AROT)
C *****************************************************************************
C     SUBROUTINE RODRIGUES (ex-REORIENT_GRAIN)   --->   VERSION 03/NOV/09

C     BUILDS INCREMENTAL ROTATION MATRIX 'AROT' BASED ON RODRIGUES FORMULA.
C             --> AROT=I+sin(phi)*C+(1-cos(phi))*C^2
C     'C'   :INCREMENTAL LATTICE SPIN (ANTI-SYMMETRIC) (normalized)
C     'AROT':INCREMENTAL TRANSFORMATION FROM INITIAL TO FINAL ORIENTATION.
C *****************************************************************************
      dimension c(3,3),c2(3,3),v(3),xid33(3,3),arot(3,3)
      data xid33/1.,0.,0.,0.,1.,0.,0.,0.,1./

c ** v(i) is the Rodrigues spin axis. The norm is the rotation angle phi.
      v(1)=c(3,2)
      v(2)=c(1,3)
      v(3)=c(2,1)
      vnorm=sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      if(vnorm.lt.1.e-06) then
cx        write(*,*) ' zero rotation inside Rodrigues'
cx        write(*,*) ' press RETURN to continue'
cx        pause
        do i=1,3
        do j=1,3
          arot(i,j)=xid33(i,j)
        enddo
        enddo
        return
      endif

      coef1=sin(vnorm)/vnorm
      coef2=(1.-cos(vnorm))/vnorm**2

      do i=1,3
      do j=1,3
        c2(i,j)=0.
        do k=1,3
          c2(i,j)=c2(i,j)+c(i,k)*c(k,j)
        enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
        arot(i,j)=xid33(i,j)+coef1*c(i,j)+coef2*c2(i,j)
      enddo
      enddo

      return
      end subroutine
      
      FUNCTION pythag(a,b)
      REAL a,b,pythag
c     Computes (a2 + b2)1=2 without destructive underflow or overflow.      
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
        pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END FUNCTION
      
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      REAL a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500) !Maximum anticipated value of n.
c     USES pythag
c     Given a matrix a(1:m,1:n), with physical dimensions mp by np, this routine computes its
c     singular value decomposition, A = U W V^T. The matrix U replaces a on output. The
c     diagonal matrix of singular values W is output as a vector w(1:n). The matrix V (not the
c     transpose V T ) is output as v(1:n,1:n).
      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX) !,pythag
      g=0.0 !Householder reduction to bidiagonal form.
      scale=0.0
      anorm=0.0
      do i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if(i.le.m)then
          do k=i,m
            scale=scale+abs(a(k,i))
          enddo 
          if(scale.ne.0.0)then
            do  k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
            enddo 
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do  j=l,n
              s=0.0
              do  k=i,m
                s=s+a(k,i)*a(k,j)
              enddo 
              f=s/h
              do k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
              enddo 
      
      
            enddo 
            do k=i,m
              a(k,i)=scale*a(k,i)
            enddo 
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if((i.le.m).and.(i.ne.n))then
          do k=l,n
            scale=scale+abs(a(i,k))
          enddo 
          if(scale.ne.0.0)then
            do k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
            enddo 
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do  k=l,n
              rv1(k)=a(i,k)/h
            enddo 
            do  j=l,m
              s=0.0
              do  k=l,n
               s=s+a(j,k)*a(i,k)
              enddo 
              do  k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
              enddo 
            enddo 
            do  k=l,n
              a(i,k)=scale*a(i,k)
            enddo 
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
      enddo 
      do  i=n,1,-1 !Accumulation of right-hand transformations.
        if(i.lt.n)then
          if(g.ne.0.0)then
            do  j=l,n !Double division to avoid possible underflow.
              v(j,i)=(a(i,j)/a(i,l))/g
            enddo 
            do  j=l,n
              s=0.0
              do  k=l,n
                s=s+a(i,k)*v(k,j)
              enddo 
              do  k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
              enddo 
            enddo 
          endif
          do  j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
          enddo 
        endif
        v(i,i)=1.0
      
      
        g=rv1(i)
        l=i
      enddo 
      do  i=min(m,n),1,-1 !Accumulation of left-hand transformations.
        l=i+1
        g=w(i)
        do  j=l,n
          a(i,j)=0.0
        enddo 
        if(g.ne.0.0)then
          g=1.0/g
          do  j=l,n
            s=0.0
            do  k=l,m
              s=s+a(k,i)*a(k,j)
            enddo 
            f=(s/a(i,i))*g
            do  k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
            enddo 
          enddo 
          do  j=i,m
            a(j,i)=a(j,i)*g
          enddo 
        else
          do  j= i,m
            a(j,i)=0.0
          enddo 
        endif
        a(i,i)=a(i,i)+1.0
      enddo 
      do  k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
        do its=1,30 
          do  l=k,1,-1 !Test for splitting.
            nm=l-1 !Note that rv1(1) is always zero.
            if((abs(rv1(l))+anorm).eq.anorm) goto 2
            if((abs(w(nm))+anorm).eq.anorm) goto 1
          enddo 
    1     c=0.0 !Cancellation of rv1(l), if l > 1.
          s=1.0
          do  i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do  j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
            enddo 
          enddo 
    2     z=w(k)
          if(l.eq.k)then !Convergence.
            if(z.lt.0.0)then !Singular value is made nonnegative.
              w(k)=-z
              do  j=1,n
                v(j,k)=-v(j,k)
              enddo       
          
          
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l) !Shift from bottom 2-by-2 minor.
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0 !Next QR transformation:
          s=1.0
          do  j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do  jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
            enddo 
            z=pythag(f,h)
            w(j)=z !Rotation can be arbitrary if z = 0.
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do  jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
            enddo 
          enddo 
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
        enddo 
    3   continue
      enddo 
      return
      END SUBROUTINE
            
      SUBROUTINE cond_numb(np,n,b,cond_number) !AAanis
      
      dimension a(np,np),b(np,np),ai(np,np),w(np),v(np,np)
      
      a=b
      ! normalize matrix
c      ANORM=TNORM(A,21,21)
c      DO I=1,6
c      DO J=1,6
c        A(I,J)=A(I,J)/ANORM
c      ENDDO
c      ENDDO      
      
      ! SVD
      call svdcmp(a,n,n,np,np,w,v)

      ! condition number
      cond_number=maxval(w(1:n))/minval(w(1:n))
c      if(cond_number.gt.1e6) pause

      ! inverse using SDV
      ai=0.0
      do i=1,6
          if (w(i).gt.1e-6) then
              ai(i,i)=1.0/w(i)
          else
              ai(i,i)=0.0
          endif
      enddo
      ai=matmul(ai,transpose(a))
      ai=matmul(v,ai)
c      
c      DO I=1,6
c        DO J=1,6
c          AI(I,J)=AI(I,J)*INVFAC(I,J)/ANORM
c        ENDDO
c      ENDDO      
c      
      ! replace the input matrix with its inverse
c      b = ai
      
      END SUBROUTINE
      
      SUBROUTINE NON_SCHMID(bcs,xncs,nscc,p)
      !The non-shmid tensor calculation 
      ! inptus are burgers vector and normal of each slip plane. the output is the symmetric part of the non-schmid tensor for each slip plane
      
      DIMENSION bcs(3), xncs(3), tcs(3), tbd(3,3), tnd(3,3),
     # xnnd(3,3), ttd(3,3),bbd(3,3), p(3,3), nscc(4)
     
      Real nscc
      
      
      call CROSS_PRODUCT(xncs,bcs,tcs)
      
      ! First term of the non-shmid tensor (the symmetric part)

      do i=1,3
          do j=1,3
              tbd(i,j)= (0.5)*(tcs(i)*bcs(j)+bcs(i)*tcs(j))
          enddo
      enddo

       ! Second term of the non-shmid tensor (the symmetric part)
      do i=1,3
          do j=1,3
              tnd(i,j)= (0.5)*(tcs(i)*xncs(j)+xncs(i)*tcs(j))
          enddo
      enddo

       ! Third term of the non-shmid tensor

      do i=1,3
          do j=1,3
              xnnd(i,j)= xncs(i)*xncs(j)
          enddo
      enddo

        
       ! Fourth term of the non-shmid tensor

      do i=1,3
          do j=1,3
              ttd(i,j)= tcs(i)*tcs(j)
          enddo
      enddo

       ! Fifth term of the non-shmid tensor

      do i=1,3
          do j=1,3
              bbd(i,j)= bcs(i)*bcs(j)
          enddo
      enddo
      
      p = (nscc(1))*tbd+(nscc(2))*tnd+(nscc(3))*xnnd+(nscc(4))
     #    *ttd-(nscc(3)+nscc(4))*bbd 
      
      Return 
      End subroutine 
        
      SUBROUTINE CROSS_PRODUCT (A, B, C)     !cross product (right-handed)

      IMPLICIT NONE                                  ! no default typing

      DOUBLE PRECISION, DIMENSION(3), INTENT (IN)    :: A   ! multiplicand 3-vector
      DOUBLE PRECISION, DIMENSION(3), INTENT (IN)    :: B   ! multiplier 3-vector
      DOUBLE PRECISION, DIMENSION(3), INTENT (OUT)   :: C   ! result: 3-vector cross product


      C(1) = A(2)*B(3) - A(3)*B(2)    ! compute cross product components
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)

      RETURN

      END SUBROUTINE CROSS_PRODUCT 
      
      end module miscellaneous_sub      
c      
c***********************************************************************
c *** voigt ************************************************************
c***********************************************************************
      module mvoigt

      !arrays
      integer :: ijv(6,2),i6(6)
      real :: invfac(6,6),profac(6),id2(6,6)
      real :: dummy1(6),dummy2(3,3),dummy3(6,6),dummy4(3,3,3,3)
      
      private
      public :: INVTEN,VOIGT,initialize_voigt,ijv,invfac,profac,id2,i6
     #          ,dummy1,dummy2,dummy3,dummy4,tens_mult,tens_mult_1
      
      contains 
      !invten       :Finds the inverse of 6x6 matrix         
      !voigt        :Assigns the vectors and matrix in Voigt 
      
      SUBROUTINE INVTEN(A,AI)
c **********************************************************************
c *** CALCULATES THE INVERSE 'AI' OF THE TENSOR 'A' USING THE 6X6 VOIGT
c     REPRESENTATION.
c *** NORMALIZE/RENORMALIZE TO WORK WITH TENSORS OF ORDER UNIT.
c **********************************************************************
c *** USES:    LUDCMPC    LUBKSBC                                    ***
c **********************************************************************
c     VERSION: 01/nov/99
c     Option of singular value decomposition was eliminated on july/99
c **********************************************************************

      use miscellaneous_sub, only : tnorm, ludcmpc, lubksbc

      DIMENSION A(6,6),AI(6,6),AX(6,6)
      DIMENSION ID(6,6),INDX(6)
c      real INVFAC,ID2,ID
      real ID
c      COMMON/VOIG/IJV(6,2),INVFAC(6,6),PROFAC(6),ID2(6,6)
      ANORM=TNORM(A,6,6)
      DO I=1,6
      DO J=1,6
        AX(I,J)=A(I,J)/ANORM
      ENDDO
      ENDDO

      CALL LUDCMPC(AX,6,6,INDX,DET)
C      DO I=1,6
C        DET=DET*AX(I,I)
C      ENDDO
      DO I=1,6
        DO J=1,6
          ID(I,J)=(I/J)*(J/I)
        ENDDO
      ENDDO
      DO I=1,6
        CALL LUBKSBC(AX,6,6,INDX,ID(1,I))
      ENDDO
      DO I=1,6
        DO J=1,6
          AI(I,J)=ID(I,J)*INVFAC(I,J)/ANORM
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE VOIGT(T1,T2,C2,C4,IOPT)
C **********************************************************************
C *** TRANSFORMS 6X1 MATRIX T1 INTO 3x3 TENSOR T2 IF IOPT=1          ***
C *** AND VICEVERSA IF IOPT=2.                                       ***
C *** TRANSFORMS 6X6 MATRIX C2 INTO 3x3x3x3 TENSOR C4 IF IOPT=3      ***
C *** AND VICEVERSA IF IOPT=4.                                       ***
c **********************************************************************
c *** VERSION 30/JUL/99                                              ***
c **********************************************************************
      DIMENSION T1(6),T2(3,3),C2(6,6),C4(3,3,3,3)
c      real INVFAC,ID2
c      COMMON/VOIG/IJV(6,2),INVFAC(6,6),PROFAC(6),ID2(6,6)
C
      IF(IOPT.EQ.1) THEN
        DO I=1,6
          I1=IJV(I,1)
          I2=IJV(I,2)
          T2(I1,I2)=T1(I)
          T2(I2,I1)=T1(I)
        ENDDO
      ENDIF
C
      IF(IOPT.EQ.2) THEN
        DO I=1,6
          I1=IJV(I,1)
          I2=IJV(I,2)
          T1(I)=T2(I1,I2)
        ENDDO
      ENDIF
C
      IF(IOPT.EQ.3) THEN
        DO I=1,6
          I1=IJV(I,1)
          I2=IJV(I,2)
          DO J=1,6
            J1=IJV(J,1)
            J2=IJV(J,2)
            C4(I1,I2,J1,J2)=C2(I,J)
            C4(I2,I1,J1,J2)=C2(I,J)
            C4(I1,I2,J2,J1)=C2(I,J)
            C4(I2,I1,J2,J1)=C2(I,J)
          ENDDO
        ENDDO
      ENDIF
C
      IF(IOPT.EQ.4) THEN
        DO I=1,6
          I1=IJV(I,1)
          I2=IJV(I,2)
          DO J=1,6
            J1=IJV(J,1)
            J2=IJV(J,2)
            C2(I,J)=C4(I1,I2,J1,J2)
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE

      subroutine initialize_voigt
      integer :: ijvx(6,2)
      DATA ((ijvx(i,j),j=1,2),i=1,6)/1,1,2,2,3,3,2,3,1,3,1,2/
      
      !Assigns values to Voigt notation matrix (ijv), identity matrix ***
      !(id2), product factor (profac) and inverse factor (invfac)     ***
      do i=1,6
        i6(i)=0
        if(i.le.3) i6(i)=1
        ijv(i,1)=ijvx(i,1)
        ijv(i,2)=ijvx(i,2)
        profac(i)=1.0+(i/4)
        do j=1,6
          id2(i,j)=(i/j)*(j/i)*(1.0-0.50*(i/4))
          invfac(i,j)=1.0/((1+i/4)*(1+j/4))
        enddo
      enddo      
      
      return
      end subroutine
      
      SUBROUTINE TENS_MULT(res,t1,t2)
      
      DIMENSION t1(6,6),t2(6,6),res(6,6)
   
      res=0.0      

      do i=1,6
        do j=1,6
          do k=1,6
            res(i,j) = res(i,j) + t1(i,k)*t2(k,j)*profac(k)
          enddo
        enddo
      enddo
      
      RETURN
      END SUBROUTINE
   
      SUBROUTINE TENS_MULT_1(res,t1,t2)
      
      DIMENSION t1(6,6),t2(6),res(6)
      
      res=0.0
      
      do i=1,6
        do j=1,6
          res(i) = res(i) + t1(i,j)*t2(j)*profac(j)
        enddo
      enddo
       
      RETURN
      END SUBROUTINE 
      
      end module mvoigt
c      
c***********************************************************************
c *** meshelby *********************************************************
c***********************************************************************
      module meshelby

      private 
      public eshelbyb
      
      contains
      !av_modulus   :Calculates Voigt & Reuss averages       
      !eshelby      :Calculates Eshelby tensor               
      !esh_inv      :Tensor inversion in Voigt notation      
      !esh_mult     :Tensor multiplication in Voigt notation 

      SUBROUTINE ESHELBYB (axis,c4,keff,esim,escr,ioption)
C ***********************************************************************
C     SUBROUTINE ESHELBY      --->      VERSION 15/NOV/07
C
C     IOPTION=0: Initialize arrays assoc. with Gauss integration points.
C     IOPTION=1: Calculate elastic Eshelby tensor for elastic inclusion.
C     IOPTION=2: Calculate incompressible Eshelby tensors ESIM (strain-
C                rate) & ESCR (spin-rate) associated with the visco-
C                plastic inclusion.
C     IOPTION=3: Calculate incompressible and hydrostatic Eshelby tensors
C                PESH (deviatoric pressure), PDIL (spherical pressure) &
C                DESH (dilatation) for visco-plastic inclusion.
C     IOPTION=4: Calculates 1st term in d(S)/d(M)=d(T)/d(M)*L+T*d(L)/d(M)
C     IOPTION=5: Calculates 2nd term in d(S)/d(M)=d(T)/d(M)*L+T*d(L)/d(M)
C
C     Options 2-3-4-5 admit a non-zero visco-plastic bulk modulus KEFF.
C
C     Algorithms are based in Lebensohn et al, MSMSE 6 (1998) p.447.
C
C     Uses explicit matrix inversion and explicit Voigt notation (when
C     possible) to optimize computing time.
C
C     Modified oct/2005 to adapt number of integration points to the shape
C     of the ellipsoid in order to keep Eshelby tensor within a certain
C     tolerance (based on analysis done by Gwenaelle Proust).
C     Aspect ratio criterion was adopted for the case AXIS(2)>AXIS(1)>AXIS(3)
C
C     Modified oct/2007 to use a Gauss-Lobatto integration with fix (ten)
C     integration points and weights, chosen to optimize the integration
C     within (eleven) different domains of the ellipsoid aspect ratios.
C     (based on analysis made by Laurent Capolungo).
C        if IGAUSSLEG=0 uses Gauss Lobatto  (cases 1 to 11)
C        if IGAUSSLEG=1 uses Gauss Legendre (case 12) (kept as an option)
C ***********************************************************************
      use mvoigt

      DIMENSION c4(3,3,3,3),esim(3,3,3,3),escr(3,3,3,3)
      DIMENSION p(3,3),pesh(3,3),desh(3,3)
      DIMENSION c2(6,6),gamma2(6,6),gamma4(3,3,3,3)
      DIMENSION axis(3),x1(10),a1(10),a1inv(10)
      DIMENSION aa1x(6),aa2x(3,3),aaww1x(6),aaww2x(3,3)
CFEB
CFEE
      PARAMETER (ngaumx=16,ngaumx2=256)
      DIMENSION xph(ngaumx),xth(ngaumx),wph(ngaumx),wth(ngaumx)
      COMMON/ESHELBY3/ngaussph(12),ngaussth(12)
      COMMON/ESHELBY4/alpha(12,3,ngaumx2),aa1(12,6,ngaumx2),
     # aww(12,3,ngaumx2),aaww1(12,6,ngaumx2),ww(12,ngaumx2)
C *** Integration points and weights, and limits of aspect ratios
      COMMON/ESHELBY5/punti(10,11),pesi(10,11),dte(0:10),
     #                puntigl(16),pesigl(16)

      REAL      keff
      INTEGER   case

      pi=4.*atan(1.0)

      IGAUSSLEG=0     ! hardwires the Gauss-Lobatto option
ccc   IGAUSSLEG=1     ! hardwires the Gauss-Legendre option

C ***********************************************************************c
C     INITIALIZATION RUN
C     Calculates Gauss-Legendre integration points and weights in the
c     interval [0,pi]. Initializes Gauss-Lobatto points and weights.
C     Initializes arrays associated with each point to avoid repeating
C     its calculation at every call. All values of the points and
C     weights were calculated to optimize the error.
C ***********************************************************************

      if(ioption.eq.0) then

        do i=1,11
          ngaussph(i)=10      ! Gauss-Lobatto
          ngaussth(i)=10
        end do
        ngaussph(12)=16       ! Gauss-Legendre
        ngaussth(12)=16

c****************************
c  CASE 1
c****************************
        punti(1,1)=4.71236594E-02
        punti(2,1)=0.241774723
        punti(3,1)=0.565131843
        punti(4,1)=0.968887568
        punti(5,1)=1.37937832
        punti(6,1)=1.76221442
        punti(7,1)=2.17270517
        punti(8,1)=2.57646084
        punti(9,1)=2.89981818
        punti(10,1)=3.09446883

        pesi(1,1)=0.120191820
        pesi(2,1)=0.264987558
        pesi(3,1)=0.373805553
        pesi(4,1)=0.420841277
        pesi(5,1)=0.390970200
        pesi(6,1)=0.390970260
        pesi(7,1)=0.420841366
        pesi(8,1)=0.373805553
        pesi(9,1)=0.264987499
        pesi(10,1)=0.120192111

c****************************
c  CASE 2
c****************************
        punti(1,2)=1.57080423E-02
        punti(2,2)=0.144995824
        punti(3,2)=0.425559640
        punti(4,2)=0.829968274
        punti(5,2)=1.31460333
        punti(6,2)=1.82698941
        punti(7,2)=2.31162453
        punti(8,2)=2.71603298
        punti(9,2)=2.99659705
        punti(10,2)=3.12588477

        pesi(1,2)=5.41692823E-02
        pesi(2,2)=0.207461149
        pesi(3,2)=0.348739326
        pesi(4,2)=0.452716887
        pesi(5,2)=0.507709801
        pesi(6,2)=0.507709682
        pesi(7,2)=0.452716798
        pesi(8,2)=0.348738998
        pesi(9,2)=0.207461327
        pesi(10,2)=5.41692935E-02

c****************************
c  CASE 3
c****************************
        punti(1,3)=3.76990959E-02
        punti(2,3)=0.198626831
        punti(3,3)=0.483041346
        punti(4,3)=0.871647120
        punti(5,3)=1.32964790
        punti(6,3)=1.81194484
        punti(7,3)=2.26994562
        punti(8,3)=2.65855122
        punti(9,3)=2.94296598
        punti(10,3)=3.10389376

        pesi(1,3)=9.68142375E-02
        pesi(2,3)=0.224478707
        pesi(3,3)=0.341134071
        pesi(4,3)=0.430180043
        pesi(5,3)=0.478189558
        pesi(6,3)=0.478189170
        pesi(7,3)=0.430180043
        pesi(8,3)=0.341134191
        pesi(9,3)=0.224478647
        pesi(10,3)=9.68143344E-02

c****************************
c  CASE 4
c****************************
        punti(1,4)=3.45576368E-02
        punti(2,4)=0.187556863
        punti(3,4)=0.468425453
        punti(4,4)=0.859980166
        punti(5,4)=1.32527423
        punti(6,4)=1.81631863
        punti(7,4)=2.28161263
        punti(8,4)=2.67316723
        punti(9,4)=2.95403576
        punti(10,4)=3.10703516

        pesi(1,4)=8.95763785E-02
        pesi(2,4)=0.217725381
        pesi(3,4)=0.341026783
        pesi(4,4)=0.435772508
        pesi(5,4)=0.486694932
        pesi(6,4)=0.486695170
        pesi(7,4)=0.435772508
        pesi(8,4)=0.341026902
        pesi(9,4)=0.217725128
        pesi(10,4)=8.95764604E-02

c****************************
c  CASE 5
c****************************
        punti(1,5)= 3.14158052E-02
        punti(2,5)=0.177928671
        punti(3,5)= 0.457155794
        punti(4,5)= 0.851592362
        punti(5,5)= 1.32222414
        punti(6,5)= 1.81936860
        punti(7,5)=2.29000044
        punti(8,5)=2.68443704
        punti(9,5)=2.96366405
        punti(10,5)=3.11017680

        pesi(1,5)=8.26927349E-02
        pesi(2,5)=0.213228315
        pesi(3,5)=0.342008322
        pesi(4,5)=0.440196186
        pesi(5,5)=0.492670894
        pesi(6,5)=0.492670983
        pesi(7,5)=0.440195888
        pesi(8,5)=0.342008322
        pesi(9,5)=0.213227972
        pesi(10,5)=8.26930404E-02

c****************************
c  CASE 6
c****************************
        punti(1,6)= 2.98452154E-02
        punti(2,6)=0.173592165
        punti(3,6)=0.452448040
        punti(4,6)=0.848216832
        punti(5,6)=1.32101476
        punti(6,6)=1.82057810
        punti(7,6)= 2.29337597
        punti(8,6)=2.68914461
        punti(9,6)=2.96800065
        punti(10,6)=3.11174774

        pesi(1,6)=7.93928578E-02
        pesi(2,6)=0.211627841
        pesi(3,6)=0.342669785
        pesi(4,6)=0.442057431
        pesi(5,6)=0.495048553
        pesi(6,6)=0.495048642
        pesi(7,6)=0.442057490
        pesi(8,6)=0.342670023
        pesi(9,6)=0.211627468
        pesi(10,6)=7.93929026E-02

c****************************
c  CASE 7
c****************************
        punti(1,7)=2.67036632E-02
        punti(2,7)=0.165752888
        punti(3,7)=0.444431901
        punti(4,7)=0.842614472
        punti(5,7)=1.31902647
        punti(6,7)= 1.82256627
        punti(7,7)=2.29897833
        punti(8,7)=2.69716072
        punti(9,7)=2.97583985
        punti(10,7)=3.11488938

        pesi(1,7)=7.30879456E-02
        pesi(2,7)=0.209402516
        pesi(3,7)=0.344104946
        pesi(4,7)=0.445234656
        pesi(5,7)=0.498966068
        pesi(6,7)= 0.498966306
        pesi(7,7)=0.445234746
        pesi(8,7)= 0.344104946
        pesi(9,7)=0.209402665
        pesi(10,7)=7.30878562E-02

c****************************
c  CASE 8
c****************************
        punti(1,8)=2.67036632E-02
        punti(2,8)=0.165752888
        punti(3,8)=0.444431901
        punti(4,8)=0.842614472
        punti(5,8)=1.31902647
        punti(6,8)=1.82256627
        punti(7,8)=2.29897833
        punti(8,8)=2.69716072
        punti(9,8)=2.97583985
        punti(10,8)=3.11488938

        pesi(1,8)=7.30879456E-02
        pesi(2,8)=0.209402516
        pesi(3,8)=0.344104946
        pesi(4,8)=0.445234656
        pesi(5,8)=0.498966068
        pesi(6,8)=0.498966306
        pesi(7,8)=0.445234746
        pesi(8,8)=0.344104946
        pesi(9,8)=0.209402665
        pesi(10,8)= 7.30878562E-02

c****************************
c  CASE 9
c****************************
        punti(1,9)=2.43473575E-02
        punti(2,9)=0.160516247
        punti(3,9)=0.439386278
        punti(4,9)=0.839168847
        punti(5,9)=1.31781363
        punti(6,9)=1.82377899
        punti(7,9)=2.30242372
        punti(8,9)=2.70220637
        punti(9,9)=2.98107672
        punti(10,9)=3.11724544

        pesi(1,9)=6.86219111E-02
        pesi(2,9)=0.208388865
        pesi(3,9)=0.345189095
        pesi(4,9)=0.447236270
        pesi(5,9)=0.501360059
        pesi(6,9)=0.501359940
        pesi(7,9)=0.447236151
        pesi(8,9)=0.345189214
        pesi(9,9)=0.208388969
        pesi(10,9)=6.86219335E-02

c****************************
c  CASE 10
c****************************
        punti(1,10)=2.19910536E-02
        punti(2,10)=0.155757755
        punti(3,10)=0.434985727
        punti(4,10)=0.836206555
        punti(5,10)=1.31677616
        punti(6,10)= 1.82481658
        punti(7,10)=2.30538607
        punti(8,10)=2.70660710
        punti(9,10)=2.98583508
        punti(10,10)=3.11960149

        pesi(1,10)=6.43825606E-02
        pesi(2,10)=0.207786217
        pesi(3,10)=0.346235514
        pesi(4,10)=0.448981822
        pesi(5,10)=0.503410578
        pesi(6,10)= 0.503410578
        pesi(7,10)=0.448981792
        pesi(8,10)=0.346235693
        pesi(9,10)=0.207785636
        pesi(10,10)= 6.43827692E-02

c****************************
c  CASE 11
c****************************
        punti(1,11)=2.04204638E-02
        punti(2,11)=0.152822554
        punti(3,11)=0.432348520
        punti(4,11)=0.834448099
        punti(5,11)=1.31616223
        punti(6,11)=1.82543063
        punti(7,11)=2.30714464
        punti(8,11)=2.70924401
        punti(9,11)=2.98877001
        punti(10,11)=3.12117243

        pesi(1,11)=6.16818815E-02
        pesi(2,11)=0.207559645
        pesi(3,11)=0.346902698
        pesi(4,11)=0.450027168
        pesi(5,11)=0.504624724
        pesi(6,11)= 0.504624426
        pesi(7,11)=0.450027317
        pesi(8,11)=0.346902847
        pesi(9,11)=0.207559645
        pesi(10,11)=6.16819337E-02

c**************************************
c  CASE 12: GAULEG generates the points
c**************************************
        call gauleg(0.0,pi,puntigl,pesigl,ngaussph(12))

C *****************************************************************
C *** Calculates and saves arrays that depend on integration points

       do case=1,12

         if (case.eq.12) then
           do i=1,16
             xph(i)=puntigl(i)
             xth(i)=puntigl(i)
             wph(i)= pesigl(i)
             wth(i)= pesigl(i)
           end do
         else
           do i=1,10
             xph(i)=punti(i,case)
             xth(i)=punti(i,case)
             wph(i)= pesi(i,case)
             wth(i)= pesi(i,case)
           end do
         end if
c--------------------------------------------------------------
c *** integration [0,pi][0,pi] adds a factor 2 in Eqs. B11 & B14.

         do ith=1,ngaussth(case)
           sinth=sin(xth(ith))
           costh=cos(xth(ith))
           simbtet=wth(ith)*sinth/(2.0*pi)

           do iph=1,ngaussph(case)
             ny=iph+(ith-1)*ngaussph(case)
             ww(case,ny)=simbtet*wph(iph)
             alpha(case,1,ny)=sinth*cos(xph(iph))
             alpha(case,2,ny)=sinth*sin(xph(iph))
             alpha(case,3,ny)=costh

             do i=1,3
             do j=1,3
               aa2x(i,j)  =alpha(case,i,ny)*alpha(case,j,ny)
               aaww2x(i,j)=aa2x(i,j)*ww(case,ny)
             enddo
             enddo
             call voigt(aa1x  ,aa2x  ,c2,c4,2)
             call voigt(aaww1x,aaww2x,c2,c4,2)
             do i=1,6
               aa1(case,i,ny)  =aa1x(i)
               aaww1(case,i,ny)=aaww1x(i)
             enddo

c *** Array AWW is used only if ICAUCHY=1.
             do i=1,3
               aww(case,i,ny)=alpha(case,i,ny)*ww(case,ny)
             enddo
           enddo
         enddo

        enddo      ! end of do case=1,12

      endif      ! ENDIF FOR IOPTION=0

C************************************************************************
C     End of initialization
C************************************************************************

C ***********************************************************************
C     CALCULATION OF ESHELBY TENSORS FOR STIFFNESS 'C4' AND ELLIPSOID
C     AXES 'AXIS'
C     ASSUMES: AXIS2 > AXIS1 > AXIS3  --> RATIO1 > RATIO2 > 1
C ***********************************************************************

      if(ioption.ge.1) then

        abc=axis(1)*axis(2)*axis(3)
        ratio1=axis(2)/axis(3)
        ratio2=axis(1)/axis(3)

        if (igaussleg.eq.1) then
          case=12
        else
          dte(0) = 0.0
          dte(1) =-0.7*ratio1+7
          dte(2) =-ratio1+17
          dte(3) =-ratio1+23
          dte(4) =-ratio1+26
          dte(5) =-ratio1+29.3
          dte(6) =-ratio1+32
          dte(7) =-ratio1+34.85
          dte(8) =-ratio1+37
          dte(9) =-ratio1+41.9
          dte(10)=-ratio1+44.5
          case=11
          do i=1,10
            if(ratio2.ge.dte(i-1) .and. ratio2.lt.dte(i)) case=i
          enddo
        endif

c--------------------------------------------------------

        npoints=ngaussph(case)*ngaussth(case)

        pdil=0.
        do j=1,3
        do i=1,3
          p(i,j)=0.
        enddo
        enddo
        do j=1,6
        do i=1,6
          gamma2(i,j)=0.
        enddo
        enddo

        call voigt(aa1x,aa2x,c2,c4,4)

        do ny=1,npoints

c   Compute Voigt components A1(1)-A(6) of tensor A(3,3) defined by Eq.B3:
c   --->  A(i,j)=L(i,j,k,l)*a(j)*a(l)

          do i=1,6
            aa1x(i)=aa1(case,i,ny)
          enddo
          call esh_mult_voigt(c2,aa1x,a1)
cw
      IF(IOPTION.EQ.1) THEN

c   If solving an elastic inclusion invert the system
c   --> A(3,3) x X(3,3) = C(3,3)
c   Inverts A(3,3) using explicit Voigt notation.
c   Uses explicit form of C(3,3) to calculate solution in Voigt notation.

            call esh_inv3_voigt(a1,a1inv)
            do i=1,6
              x1(i)=a1inv(i)
            enddo

      ENDIF

          ro3=((alpha(case,1,ny)*axis(1))**2+
     #         (alpha(case,2,ny)*axis(2))**2+
     #         (alpha(case,3,ny)*axis(3))**2)**1.5
          abcoro3=abc/ro3

c   Compute the Eshelby integral Eq.B11 defining:
c         Gamma(m,j,n,i)=T(m,n,i,j)=a(m)*a(j)*G(n,i)
c   with the property:
c         Gamma(m,j,n,i)=Gamma(j,m,n,i)=Gamma(m,j,i,n)

          do i=1,6
          do j=1,6
            gamma2(i,j)=gamma2(i,j)+aaww1(case,i,ny)*x1(j)*abcoro3
          enddo
          enddo

        end do   ! end of loop over double integration

c ********************************************************************
c   Go back to the 3*3*3*3 notation
        call voigt(aa1x,aa2x,gamma2,gamma4,3)

c   Compute symmetric (distortion) Eshelby tensor from Eq.B9.
c       esim(n,m,k,l)=0.5*(gamma(m,j,n,i)+gamma(n,j,m,i))*c4(i,j,k,l)
c   Compute anti-symmetric (rotation) Eshelby tensor from Eq.B9.
c       escr(n,m,k,l)=0.5*(gamma(m,j,n,i)-gamma(n,j,m,i))*c4(i,j,k,l)

        do l=1,3
        do k=1,3
        do m=1,3
        do n=1,3
c
          dumsim=0.
          dumscr=0.
c
        do j=1,3
        do i=1,3
c
          dumsim=dumsim+(gamma4(m,j,n,i)+gamma4(n,j,m,i))*c4(i,j,k,l)
          dumscr=dumscr+(gamma4(m,j,n,i)-gamma4(n,j,m,i))*c4(i,j,k,l)
c
        enddo
        enddo
c
            esim(n,m,k,l)=0.5*dumsim
            escr(n,m,k,l)=0.5*dumscr
c
        enddo
        enddo
        enddo
        enddo


      endif      !  endif for IOPTION.GE.1

      RETURN
      END SUBROUTINE

      SUBROUTINE ESH_INV3_VOIGT (A,AINV)
C
C *************************************************************************
C     SUBROUTINE ESH_INV3_VOIGT   --->   version 23/jul/01
C
C     Inverts the 3x3 symmetric matrix 'A' using explicit Voigt notation:
C     11->1, 22->2, 33->3, 23=32->4, 31=13->5, 12=21->6
C *************************************************************************
      DIMENSION A(10),AINV(10)

      DET = A(1)*A(2)*A(3) + 2*A(4)*A(5)*A(6) - A(1)*A(4)*A(4)
     #     - A(2)*A(5)*A(5) - A(3)*A(6)*A(6)

      AINV(1) = ( A(2)*A(3) - A(4)*A(4))/DET
      AINV(2) = ( A(1)*A(3) - A(5)*A(5))/DET
      AINV(3) = ( A(1)*A(2) - A(6)*A(6))/DET
      AINV(4) = (-A(1)*A(4) + A(5)*A(6))/DET
      AINV(5) = ( A(4)*A(6) - A(2)*A(5))/DET
      AINV(6) = (-A(3)*A(6) + A(4)*A(5))/DET

      RETURN
      END SUBROUTINE

      SUBROUTINE ESH_INV4_VOIGT (A,AINV)
C
C **********************************************************************
C     SUBROUTINE ESH_INV4_VOIGT   --->   VERSION 20/JUL/01

C     Inverts the 4*4 symmetric matrix 'A' using explicit Voigt notation:
C     11-->1, 22-->2, 33-->3, 23=32-->4, 31=13-->5, 12=21-->6
C     14-->7, 24-->8, 34-->9, 44-->10.
C **********************************************************************
      DIMENSION A(10),AINV(10)


      ainv(1) = a(2)*a(3)*a(10)+2*a(4)*a(8)*a(9) -
     #          a(2)*a(9)*a(9)-a(3)*a(8)*a(8)-a(4)*a(4)*a(10)

      ainv(2) = a(1)*a(3)*a(10)+2*a(5)*a(7)*a(9) -
     #          a(1)*a(9)*a(9)-a(3)*a(7)*a(7)-a(5)*a(5)*a(10)

      ainv(3) = a(1)*a(2)*a(10)+2*a(6)*a(7)*a(8) -
     #          a(1)*a(8)*a(8)-a(2)*a(7)*a(7)-a(6)*a(6)*a(10)

      ainv(4) = a(1)*a(4)*a(10)+a(5)*a(7)*a(8)+a(6)*a(7)*a(9) -
     #          a(1)*a(8)*a(9)-a(4)*a(7)*a(7)-a(5)*a(6)*a(10)
      ainv(4) =-ainv(4)

      ainv(5) = a(4)*a(6)*a(10)+a(2)*a(7)*a(9)+a(5)*a(8)*a(8) -
     #          a(4)*a(7)*a(8)-a(6)*a(8)*a(9)-a(2)*a(5)*a(10)

      ainv(6) = a(3)*a(6)*a(10)+a(5)*a(8)*a(9)+a(4)*a(7)*a(9) -
     #          a(3)*a(7)*a(8)-a(6)*a(9)*a(9)-a(4)*a(5)*a(10)
      ainv(6) =-ainv(6)

      ainv(7) = a(4)*a(6)*a(9)+a(4)*a(5)*a(8)+a(2)*a(3)*a(7) -
     #          a(4)*a(4)*a(7)-a(2)*a(5)*a(9)-a(3)*a(6)*a(8)
      ainv(7) =-ainv(7)

      ainv(8) = a(1)*a(4)*a(9)+a(5)*a(5)*a(8)+a(3)*a(6)*a(7) -
     #          a(4)*a(5)*a(7)-a(5)*a(6)*a(9)-a(1)*a(3)*a(8)

      ainv(9) = a(1)*a(2)*a(9)+a(5)*a(6)*a(8)+a(4)*a(6)*a(7) -
     #          a(2)*a(5)*a(7)-a(6)*a(6)*a(9)-a(1)*a(4)*a(8)
      ainv(9) =-ainv(9)

      ainv(10)=a(1)*a(2)*a(3)+2*a(4)*a(5)*a(6) -
     #         a(1)*a(4)*a(4)-a(2)*a(5)*a(5)-a(3)*a(6)*a(6)

      det=   a(1)*ainv(1)+   a(2)*ainv(2)+   a(3)*ainv(3)+
     #    2.*a(4)*ainv(4)+2.*a(5)*ainv(5)+2.*a(6)*ainv(6)+
     #    2.*a(7)*ainv(7)+2.*a(8)*ainv(8)+2.*a(9)*ainv(9)+
     #       a(10)*ainv(10)
      det=   det/4.

      do i=1,10
        ainv(i)=ainv(i)/det
      enddo

      return
      end subroutine

      SUBROUTINE ESH_MULT_VOIGT(B,C,A)

C     Performs the multiplication:
C        A(i,k)=B(i,j,k,l)*C(j,l) using Voigt's notation
C        B is a 6*6 symmetric matrix
C        C is a 3*3 symmetric tensor
C        A will be a 3*3 symmetric tensor

      DIMENSION B(6,6),C(6),A(9)

      A(1)=B(1,1)*C(1)+B(6,6)*C(2)+B(5,5)*C(3)
     #    +2*(B(5,6)*C(4)+B(1,5)*C(5)+B(1,6)*C(6))

      A(2)=B(6,6)*C(1)+B(2,2)*C(2)+B(4,4)*C(3)
     #    +2*(B(2,4)*C(4)+B(4,6)*C(5)+B(2,6)*C(6))

      A(3)=B(5,5)*C(1)+B(4,4)*C(2)+B(3,3)*C(3)
     #    +2*(B(3,4)*C(4)+B(3,5)*C(5)+B(4,5)*C(6))

      A(4)=B(5,6)*C(1)+B(2,4)*C(2)+B(3,4)*C(3)
     #      +(B(2,3)+B(4,4))*C(4)
     #      +(B(3,6)+B(4,5))*C(5)
     #      +(B(4,6)+B(2,5))*C(6)

      A(5)=B(1,5)*C(1)+B(4,6)*C(2)+B(3,5)*C(3)
     #      +(B(3,6)+B(4,5))*C(4)
     #      +(B(1,3)+B(5,5))*C(5)
     #      +(B(1,4)+B(5,6))*C(6)

      A(6)=B(1,6)*C(1)+B(2,6)*C(2)+B(4,5)*C(3)
     #      +(B(4,6)+B(2,5))*C(4)
     #      +(B(1,4)+B(5,6))*C(5)
     #      +(B(1,2)+B(6,6))*C(6)

      RETURN
      END SUBROUTINE
      
      subroutine gauleg(x1,x2,x,w,n)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SUBROUTINE GAULEG(X1,X2,X,W,N)
C
C     Given lower and upper limits of integration (x1 & x2) returns arrays
C     x(n) & w(n) containing abcisas and weights of Gauss-Legendre quadrature
C     formula.
C     Scales interval to (-1,+1) to find the roots. Roots are symmetric with
C     respect to zero. Scales back to the original interval at the end.
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      
      dimension x(n),w(n)
      parameter(eps=3.0e-14)
      pi=4.0*atan(1.0)
      m=(n+1)/2
      xm=0.50*(x1+x2)
      xl=0.50*(x2-x1)
      xn=n
      do 12 i=1,m
        xi=i
        z=cos(pi*(xi-.250)/(xn+0.50))
c       BJCL RELEASE version will not run without this dummy print statement
c        print*,''
c        write(*,*)'z= ',z,xm
1       continue
        p1=1.0
        p2=0.0
        do 11 j=1,n
          xj=j
          p3=p2
          p2=p1
          p1=((2.0*j-1.0)*z*p2-(xj-1.0)*p3)/xj
11      continue
        pp=n*(z*p1-p2)/(z*z-1.0)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.eps) go to 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.0*xl/((1.0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      end subroutine

      end module meshelby
c      
c***********************************************************************
c *** mphase_props *****************************************************
c***********************************************************************
      module mphase_props
        
      use const      
      use mphase_props_v
      
      contains
      !crystal_sym  :Applies symmetry ops to system or diff vectors
      !data_crystal :Reads SX elastic & plastic information
      
      subroutine crystal_symmetry (ioption,ur1,icrysym,
     #                             isn,sn,sneq,isb,sb,nequiv)
c
c ***********************************************************************
c     subroutine crystal_symmetry   --->   version 09/JAN/2009
c
c *** If IOPTION=1:
c     Reads crystal symmetry 'icrysym' and unit cell parameters.
c     Generates vectors 'cvec(i,n)' of the unit cell.
c     Generates symmetry operators 'hh(i,j,nsymop)' for all crystal symmetries.
c *** If IOPTION=2:
c     Reads Miller indices of systems in 3 or 4-index notation 'isn(i)'
c     & 'isb(i)'. Calculates normal & burgers vectors 'sn(i)' & 'sb(i)'
c *** If IOPTION=3:
c     Generates 'nequiv' crystallographically equivalent orientations sneq(i,n)
c     of normal vector sn(i) by applying all the symmetry operations to it.
c     Discards repeated orientations and defines 'nequiv'.
c *** Simmetry parameter ICRYSYM:
c        1: CUBIC
c        2: HEXAGONAL
c        3: TRIGONAL
c        4: TETRAGONAL
c        5: ORTHORHOMBIC
c        6: MONOCLINIC
c        7: TRICLINIC
c
c ***********************************************************************
      dimension hh(3,3,24),hx(3,3,6),itag(24)
      dimension isn(4),sn(3),sneq(3,24),isb(4),sb(3)
      dimension cdim(3),cang(3),cvec(3,3)
      integer ur1
      character crysym*5
      save h,nsymop,cvec
      data pi /3.1415926535898/

c ****************************************************************************

      if(ioption.eq.1) then

        read(ur1,*)
        read(ur1,'(a)') crysym
        icrysym=0
        if(crysym.eq.'cubic' .or. crysym.eq.'CUBIC') icrysym=1
        if(crysym.eq.'hexag' .or. crysym.eq.'HEXAG') icrysym=2
        if(crysym.eq.'trigo' .or. crysym.eq.'TRIGO') icrysym=3
        if(crysym.eq.'tetra' .or. crysym.eq.'TETRA') icrysym=4
        if(crysym.eq.'ortho' .or. crysym.eq.'ORTHO') icrysym=5
        if(crysym.eq.'monoc' .or. crysym.eq.'MONOC') icrysym=6
        if(crysym.eq.'tricl' .or. crysym.eq.'TRICL') icrysym=7
        if(icrysym.eq.0) then
          write(*,*) ' *** CANNOT RECOGNIZE THE CRYSTAL SYMMETRY'
          stop
        endif

        READ(UR1,*) (CDIM(i),i=1,3),(CANG(i),i=1,3)
        DO I=1,3
          CANG(I)=CANG(I)*PI/180.
        ENDDO

c *** assumes 'a' coincident with 'x' and 'b' in the plane 'xy'
c       CVEC(1,1)=1.
c       CVEC(2,1)=0.
c       CVEC(3,1)=0.
c       CVEC(1,2)=COS(CANG(3))
c       CVEC(2,2)=SIN(CANG(3))
c       CVEC(3,2)=0.
c       CVEC(1,3)=COS(CANG(2))
c       CVEC(2,3)=(COS(CANG(1))-COS(CANG(2))*COS(CANG(3)))/SIN(CANG(3))
c       CVEC(3,3)=SQRT(1.-CVEC(1,3)**2-CVEC(2,3)**2)

c *** assumes 'c' coincident with 'z' and 'a' in the plane 'xz'
        CVEC(1,1)=SIN(CANG(2))
        CVEC(2,1)=0.
        CVEC(3,1)=COS(CANG(2))
        CVEC(1,2)=(COS(CANG(3))-COS(CANG(1))*COS(CANG(2)))/SIN(CANG(2))
        CVEC(3,2)=COS(CANG(1))
        CVEC(2,2)=SQRT(1.-CVEC(1,2)**2-CVEC(3,2)**2)
        CVEC(1,3)=0.
        CVEC(2,3)=0.
        CVEC(3,3)=1.

        DO J=1,3
        DO I=1,3
          CVEC(I,J)=CDIM(J)*CVEC(I,J)
        ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          DO M=1,6
            HX(I,J,M)=0.d0
          ENDDO
          DO N=1,24
            hh(I,J,N)=0.d0
          ENDDO
        ENDDO
        ENDDO

c *** identity operation ---> triclinic & all symmetries
      do i=1,3
        hh(i,i,1)=1.d0
      enddo
      nsymop=1

c *** 180 deg rotation around (001) ---> orthorhombic, monoclinic
      if(icrysym.eq.5 .or. icrysym.eq.6) then
        hh(1,1,2)= cos(pi)
        hh(2,2,2)= cos(pi)
        hh(3,3,2)= 1.d0
        hh(1,2,2)=-sin(pi)
        hh(2,1,2)= sin(pi)
        nsymop=2
      endif

c *** x-mirror & y-mirror ---> orthorhombic
      if(icrysym.eq.5) then
        hh(1,1,3)=-1.d0
        hh(2,2,3)= 1.d0
        hh(3,3,3)= 1.d0

        hh(1,1,4)= 1.d0
        hh(2,2,4)=-1.d0
        hh(3,3,4)= 1.d0
        nsymop=4
      endif

c *** cubic symmetry
      if(icrysym.eq.1) then

c *** rotations of (pi/3) & (2*pi/3) around <111>
        hx(1,3,1)= 1.d0
        hx(2,1,1)= 1.d0
        hx(3,2,1)= 1.d0

        hx(1,2,2)= 1.d0
        hx(2,3,2)= 1.d0
        hx(3,1,2)= 1.d0

        do m=1,2
          do n=1,nsymop
            mn=m*nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              hh(i,j,mn)=hh(i,j,mn)+hx(i,k,m)*hh(k,j,n)
            enddo
            enddo
            enddo
          enddo
        enddo
        nsymop=mn

c *** mirror across the plane (110)
        hx(1,2,3)= 1.d0
        hx(2,1,3)= 1.d0
        hx(3,3,3)= 1.d0

        do n=1,nsymop
          mn=nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              hh(i,j,mn)=hh(i,j,mn)+hx(i,k,3)*hh(k,j,n)
            enddo
            enddo
            enddo
        enddo
        nsymop=mn

c *** rotations of 90, 180, 270 around x3

        do m=1,3
          ang=pi/2.*float(m)
          hx(1,1,m)= cos(ang)
          hx(2,2,m)= cos(ang)
          hx(3,3,m)= 1.0
          hx(1,2,m)=-sin(ang)
          hx(2,1,m)= sin(ang)
          hx(1,3,m)= 0.0
          hx(3,1,m)= 0.0
          hx(2,3,m)= 0.0
          hx(3,2,m)= 0.0
        enddo

        do m=1,3
          do n=1,nsymop
            mn=m*nsymop+n
              do i=1,3
              do j=1,3
              do k=1,3
                hh(i,j,mn)=hh(i,j,mn)+hx(i,k,m)*hh(k,j,n)
              enddo
              enddo
              enddo
          enddo
        enddo
        nsymop=mn

      endif                    !end of condition for icrysym=1

c *** hexagonal, trigonal and tetragonal symmetry

      if(icrysym.ge.2 .and. icrysym.le.4) then
        if(icrysym.eq.2) nrot=6
        if(icrysym.eq.3) nrot=3
        if(icrysym.eq.4) nrot=4

c *** mirror plane at 30 deg or 60 deg or 45 deg with respect to x1
        ang=pi/float(nrot)
        hh(1,1,2)= cos(ang)**2-sin(ang)**2
        hh(2,2,2)=-hh(1,1,2)
        hh(3,3,2)= 1.d0
        hh(1,2,2)= 2.*cos(ang)*sin(ang)
        hh(2,1,2)= hh(1,2,2)
        nsymop=2

c *** rotations of 2*pi/6 around axis <001> for hexagonals.
c *** rotations of 2*pi/3 around axis <001> for trigonals.
c *** rotations of 2*pi/4 around axis <001> for tetragonals.
        do nr=1,nrot-1
          ang=nr*2.*pi/nrot
          hx(1,1,nr)= cos(ang)
          hx(2,2,nr)= cos(ang)
          hx(3,3,nr)= 1.d0
          hx(1,2,nr)=-sin(ang)
          hx(2,1,nr)= sin(ang)
        enddo

        do m=1,nrot-1
          do n=1,nsymop
            mn=m*nsymop+n
            do i=1,3
            do j=1,3
            do k=1,3
              hh(i,j,mn)=hh(i,j,mn)+hx(i,k,m)*hh(k,j,n)
            enddo
            enddo
            enddo
          enddo
        enddo
        nsymop=mn

      endif               !end of condition for icrysym= 2,3,4

c     write(10,*)
c     write(10,'(''  # of symmetry operations='',i4)') nsymop
c     write(10,'(''  symmetry matrices'')')
c     write(10,'(i3,9f7.3)') (n,((hh(i,j,n),j=1,3),i=1,3),n=1,nsymop)

      endif               !end of condition for ioption=1

c **************************************************************************
c   Converts Miller-Bravais indices of plane normal and slip direction
c   into normalized vectors sn(i) and sb(i), respectively.
c   Indices for cubic (1), tetragonal (4), orthorhombic (5), monoclinic (6)
c   & triclinic (7) systems are in 3-index notation.
c   For hexagonal (2) & trigonal (3) systems uses 4-index notation.
c **************************************************************************

      if (ioption.eq.2 .or. ioption.eq.3) then    ! iopt=3 necessary for EPSC

        if(icrysym.eq.2 .or. icrysym.eq.3) then
          isn(3)=isn(4)
          isb(1)=isb(1)-isb(3)
          isb(2)=isb(2)-isb(3)
          isb(3)=isb(4)
        endif

c *** assumes 'a' coincident with 'x' and 'b' in the plane 'xy'
c       sn(1)= isn(1)/cvec(1,1)
c       sn(2)=(isn(2)-cvec(1,2)*sn(1))/cvec(2,2)
c       sn(3)=(isn(3)-cvec(1,3)*sn(1)-cvec(2,3)*sn(2))/cvec(3,3)
c *** assumes 'c' coincident with 'z' and 'a' in the plane 'xz'
        sn(3)= isn(3)/cvec(3,3)
        sn(1)=(isn(1)-cvec(3,1)*sn(3))/cvec(1,1)
        sn(2)=(isn(2)-cvec(1,2)*sn(1)-cvec(3,2)*sn(3))/cvec(2,2)

        snnor=sqrt(sn(1)**2+sn(2)**2+sn(3)**2)
        do j=1,3
          sn(j)=sn(j)/snnor
          if(abs(sn(j)).lt.1.e-03) sn(j)=0.
        enddo

c *** Burgers vector calculation necessary for EPSC & VPSC
      if (ioption.eq.2) then
        do i=1,3
          sb(i)=isb(1)*cvec(i,1)+isb(2)*cvec(i,2)+isb(3)*cvec(i,3)
        enddo
        sbnor=sqrt(sb(1)**2+sb(2)**2+sb(3)**2)
        do j=1,3
          sb(j)=sb(j)/sbnor
          if(abs(sb(j)).lt.1.e-03) sb(j)=0.
        enddo
      endif      ! end of if(ioption.eq.2)

      endif      ! end of if(ioption.eq. 2 or 3)

c **************************************************************************
c *** generates all symmetry related vectors sneq(i,n) with z>0.
c *** eliminates redundant poles: coincidents and opposites
c **************************************************************************

      if(ioption.eq.3) then

        do n=1,nsymop
          itag(n)=0
          do i=1,3
          sneq(i,n)=0.d0
            do j=1,3
              sneq(i,n)=sneq(i,n)+hh(i,j,n)*sn(j)
            enddo
          enddo
        enddo

        if(icrysym.ne.7) then      ! nsymop=1 for trigonal
          do m=1,nsymop-1
            if(itag(m).eq.0) then
              do n=m+1,nsymop
                sndif=abs(sneq(1,m)-sneq(1,n))+abs(sneq(2,m)-sneq(2,n))
     #               +abs(sneq(3,m)-sneq(3,n))
                if(sndif .le. 0.0001) itag(n)=1
                sndif=abs(sneq(1,m)+sneq(1,n))+abs(sneq(2,m)+sneq(2,n))
     #               +abs(sneq(3,m)+sneq(3,n))
                if(sndif .le. 0.0001) itag(n)=1
              enddo
            endif
          enddo
        endif

        nequiv=0
        do n=1,nsymop
          if(itag(n).eq.0) then
            nequiv=nequiv+1
            isign=1
            if(sneq(3,n).lt.0.) isign=-1
            sneq(1,nequiv)=isign*sneq(1,n)
            sneq(2,nequiv)=isign*sneq(2,n)
            sneq(3,nequiv)=isign*sneq(3,n)
          endif
        enddo

      endif            !end of if(ioption=3)
c **************************************************************************

      return
      end subroutine

      SUBROUTINE data_crystal (filecrys,icrysym,iph)
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     SUBROUTINE data_crystal      -->      version JUN/16/2010
c **********************************************************************
c *** Reads CRYSTAL data (elastic, thermal, hardening) from "filecrys"
c *** Calculates Schmid tensors and initializes hardening parameters.
c **********************************************************************
c *** USES:    invten, crystal_symmetry                              ***
c **********************************************************************
      use mvoigt
      use flags, only : kSM,itwinning,inonSch,iTwPh
      use miscellaneous_sub, only : non_schmid
      
      CHARACTER prosa*78,filecrys*150,namesys*7,crysym*5
      DIMENSION isn(4),sn(3),sneq(3,24),isb(4),sb(3),ipole(4)
      DIMENSION mode(NMOD)
      DIMENSION fijx(3,3),fnew(3,3),aaa(3,3),ccc4(3,3,3,3)
      DIMENSION aux6(6),aux33(3,3),ap(3,3)
c ______________________________________________________________________

    1 FORMAT(a)
    2 FORMAT(1h ,78('*'))
    3 FORMAT(1h ,7('*'),' CRYSTAL data - File: ',a,8('*'))
    4 FORMAT(1h ,6E12.4)
    5 FORMAT(1h ,'Plastic modes used: ',i3)
    6 FORMAT(1h ,a,2i5,/1h ,4E12.4,/1h ,4E12.4,/1h ,8E12.4)
    8 FORMAT(1h ,4i3,3x,4i3)
   10 FORMAT(1h ,3i3,3x,3i3)
   12 FORMAT(1h ,a)
c ______________________________________________________________________

c *** Opens the file with the CRYSTAL data "filecrys"
      OPEN(unit=1,file=filecrys,status='old')

C -----------------------------------------------------------------

c *** Reads crystal symmetry and unit cell parameters.
c *** Generates all symmetry operations associated with CRYSYM.

      call crystal_symmetry (1,1,icrysym,isn,sn,sneq,isb,sb,npoles)
      nind=3
      if(icrysym.eq.2 .or. icrysym.eq.3) nind=4

C *** READS SINGLE CRYSTAL ELASTIC STIFFNESS
      read(1,1)  prosa
      read(1,*)  ((ccc2(i,j,iph),j=1,6),i=1,6) !MZ_2ph ccc2 -> extra dimension

      call invten(ccc2(:,:,iph),scc2(:,:,iph))   ! Calculates elastic compliance matrix !MZ_2ph ccc2,scc2 -> extra dimension

C Merkel: 05/2010 --> Large elastic strain and pressure dependence
      read(1,1)   prosa
      read(1,*)   kSM
      if (kSM.eq.1) then
C *** READS FIRST PRESSURE DERIVATIVE OF SINGLE CRYSTAL ELASTIC STIFFNESS
        read(1,1)  prosa
        read(1,*)  ((ccc2dp(i,j,iph),j=1,6),i=1,6) !MZ_2ph ccc2dp -> extra dim
C *** Save zero pressure values of elastic constants
        do i=1,6
         do j=1,6
           ccc2p0(i,j,iph)=ccc2(i,j,iph) !MZ_2ph ccc2,ccc2p0 -> extra dimension
         enddo
        enddo
      endif

C *** READS SINGLE CRYSTAL THERMAL EXPANSION COEFFICIENTS
      read(1,1)    prosa
      read(1,*)    (alfacc(i,iph),i=1,6)

C *** READS CRYSTALLOGRAPHIC MODE PARAMETERS
      read(1,1) prosa
      read(1,*) nmodx
      read(1,*) nmodes(iph) 
      mode(:)=0 !MZ_2ph bug fix (restart array for each SX file)
      read(1,*) (mode(i),i=1,nmodes(iph)) 

      if (nmodes(iph).gt.NMOD) then 
        write(*,'(1h ,''ERROR: Number of plastic modes greater than
     #  code dimension !!!'',/,1h ,''DIMENSION in code = '',i3)') NMOD
        write(*,*)
        write(*,'(1h ,''STOP IN ROUTINE *** data_crystal ***'')')
        stop
      endif

      im=1                 ! counter for number of modes
      isys=0               ! counter for number of systems
      nslsys(iph)=0             ! counter for number of slip systems
      ntwsys(iph)=0             ! counter for number of twin systems
      nmodes(iph)=nmodes(iph)        ! number of modes used 
      nslmod(iph)=0             ! counter for number of slip modes
      ntwmod(iph)=0             ! counter for number of twin modes

      do iloop=1,nmodx   ! loops over all modes in input file

        read(1,1) namesys
        read(1,*) modex,nsmx,nrsx,iopsysx,itwx
        if (inonSch.eq.1) then !inonSch
            READ(1,*)  c1, c2, c3, c4
        endif
c     ADDEDMZ TWSH - read threshold parameters of twin        
        if(itwx.eq.1) read(1,*) stwx,ISECTWX,THRES1X,THRES2X  !characteristic twin shear 
c     END
        nrsx=nrsx                                !just to fool the compiler

        if(modex.ne.iloop) then
          write(*,*) ' WARNING !!!'
          write(*,*) ' MODE NUMBERS MUST BE SEQUENTIAL IN CRYSTAL FILE'
          STOP
        endif

        if (iloop.ne.mode(im)) then
          do is=1,nsmx
            read(1,*)
          enddo
        endif

        if (iloop.eq.mode(im)) then

          if(iopsysx.eq.1 .and. itwx.eq.1) then
c            write(*,*) ' WARNING !!!'
c            write(*,*) ' IOPSYSX MUST BE =0 WHEN ITWX=1'
c            STOP
          endif

          nsm(im,iph)=(iopsysx+1)*nsmx
          itw(im)=itwx             ! itw=0 for slip , itw=1 for twin mode
          if(itw(im).eq.0) then
            stw(im,iph)=0
            nslmod(iph)=nslmod(iph)+1
            nslsys(iph)=nslsys(iph)+nsm(im,iph)
          endif
          if(itw(im).eq.1) then
            stw(im,iph)=stwx
            itwinning=1            ! Setting flag for active twin systems
            iTwPh(iph)=1
            ntwmod(iph)=ntwmod(iph)+1
            ntwsys(iph)=ntwsys(iph)+nsm(im,iph)       
            ISECTW(im)=ISECTWX ! not used, should go over system
            TWTHRES(1,im,iph)=THRES1X !TWTHRES(1:2,ntwmod_absolute number)
            TWTHRES(2,im,iph)=THRES2X    
          endif
          
          if (inonSch.eq.1) then !inonSch
              nscc(1,im,iph) = c1
              nscc(2,im,iph) = c2
              nscc(3,im,iph) = c3
              nscc(4,im,iph) = c4              
          endif

c *** Reads Miller indices of systems and transforms to cartesian vectors.
c *** Defines opposite system for slip modes. Odd=direct, even=opposite.
c *** Stores Burgers vector for twinning reorientation calculations (aug/98).
c *** Calculates Schmid tensor.

          do is=1,nsmx

            read(1,*) (isn(i),i=1,nind),(isb(i),i=1,nind)
            call crystal_symmetry
     #                     (2,1,icrysym,isn,sn,sneq,isb,sb,npoles)

            prod=sn(1)*sb(1)+sn(2)*sb(2)+sn(3)*sb(3)
            if(prod.ge.0.000001) then
              WRITE(*,'('' SYSTEM IS NOT ORTHOGONAL !!'')')
              WRITE(*,'('' ISN='',3I7)') (ISN(J),J=1,3)
              WRITE(*,'('' ISB='',3I7)') (ISB(J),J=1,3)
              WRITE(*,'(''   N='',3F7.3)') (SN(J),J=1,3)
              WRITE(*,'(''   B='',3F7.3)') (SB(J),J=1,3)
              STOP
            endif

            isys=isys+1
            do i=1,3
              bcc(i,isys,iph)= sb(i)
              ncc(i,isys,iph)= sn(i)
              do j=1,3
                mc2(i,j,isys,iph)= 0.50*(sb(i)*sn(j)+sb(j)*sn(i))
                qc2(i,j,isys,iph)= 0.50*(sb(i)*sn(j)-sb(j)*sn(i))
              enddo
            enddo
            
            if (inonSch.eq.1) then !inonSch
                !n-s effect, defining the non schmid tensor
                call NON_SCHMID(sb,sn,nscc(:,im,iph),ap)
                nmc2(:,:,isys,iph)=ap(:,:)            
            endif
            
            if(iopsysx.eq.1) then
              do i=1,3
                bcc(i,isys+1,iph)=-sb(i)
                ncc(i,isys+1,iph)= sn(i)
                do j=1,3
                  mc2(i,j,isys+1,iph)=-mc2(i,j,isys,iph)
                  qc2(i,j,isys+1,iph)=-qc2(i,j,isys,iph)
                enddo
              enddo          
              
              if (inonSch.eq.1) then !inonSch
                  !n-s effect, defining the non schmid tensor
                  call NON_SCHMID(-sb,sn,nscc(:,im,iph),ap)
                  nmc2(:,:,isys+1,iph)=ap(:,:)            
              endif 
              
              isys=isys+1
            endif

          enddo      ! end of loop over systems in the mode 'nsmx'

          if (isys.gt.NSLS) then
            write(*,'(1h ,''ERROR: Number of systems greater than code''
     #      ,'' dimension !!!'',/,1h ,''Code DIMENSION = '',i3)') NSLS
            write(*,*)
            write(*,'(1h ,''STOP IN ROUTINE *** data_crystal ***'')')
            stop
          endif
          if (im.gt.NMOD) then
            write(*,'(1h ,''ERROR: Number of modes greater than code''
     #      ,'' dimension !!!'',/,1h ,''Code DIMENSION = '',i3)') NMOD
            write(*,*)
            write(*,'(1h ,''STOP IN ROUTINE *** data_crystal ***'')')
            stop
          endif
          im=im+1

        endif    ! end of if(iloop.eq.mode(im))

      enddo      ! end of loop over total number of modes 'nmodx'

      nsys(iph)=isys

c     Set up an array that shows whether a system is a twin system or not BJCL
      nst=0
      do im=1,nmodes(iph) 
        do is=1,nsm(im,iph)
          nst=nst+1
          iSysMode(nst,iph)=im
          mode_slip(im,is,iph)=nst !MZ_2ph bug fix
          if(itw(im).eq.1) then
            iTwinSys(nst)=1
          else
            iTwinSys(nst)=0
          endif
        enddo
      enddo

cc ______________________________________________________________________
cc   Initialize values for fijph --> CT: need to move to appropriate subr !MZ_2ph commented and moved to InitializeShape(iph)
cc_______________________________________________________________________
c      da=eulerph(1,1)
c      db=eulerph(2,1)
c      dc=eulerph(3,1)
c
c      DO I=1,3
c        AXISPH(0,I,iph)=axis(I,iph) !MZ_pseudo set up array defining ellipsoid per phase
c        do j=1,ngr
c          axisgr(0,i,j)=axis(i,ngrnph(j))
c        end do
c      END DO
c
c      call euler(2,da,db,dc,aaa)
c
c      do i=1,3
c      do j=1,3
c        fijx(i,j)=(i/j)*(j/i)*AXISPH(0,I,iph)
c      enddo
c      enddo
c
c      do j=1,3
c      do i=1,3
c        fnew(i,j)=0.
c        do m=1,3
c          fnew(i,j)=fnew(i,j)+aaa(m,i)*fijx(m,j)
c        enddo
c      enddo
c      enddo
c
c      do i=1,3
c      do j=1,3
c        fijph(i,j)=fnew(i,j) !MZ_pseudo set up fijph per phase since shape of ellipsoids can be different
c      enddo
c      enddo
c
c      if (ishape.ge.2) then
c        do k=1,ngr
c          do i=1,3
c            do j=1,3
c              fijgr(i,j,k)=fnew(i,j)
c            end do
c          enddo
c        end do
c      end if
c
cc________________________________________________________________________

      return
      END SUBROUTINE
      
      SUBROUTINE temp_dep(temp)
c **********************************************************************
c *** Calculates stiffness moduli in units of GPa and thermal        ***
c *** expansion coefficients as a function of temperature            ***
c ***                                                                ***
c **********************************************************************
c *** USES:    invten                                                ***
c **********************************************************************
c *** VERSION 12/OCT/94                                              ***
c **********************************************************************
      use mvoigt, only : invten
      use sample_props_v, only : nph

      PARAMETER (C_REF=0.00001)      ! Sets units of stiffness

      tempf=temp
c ______________________________________________________________________
c
c *** These coefficients are for ZIRCONIUM                           ***
c *** Fits the curves obtained by:                                   ***
c *** FISHER,E.S. and RENKEN,C.J. Phys.Rev. 135 2A (1964) A482-A494. ***
c      ccc2(1,1,iph)=(159430-58.133*tempf+1.447d-2*tempf**2  
c     #           -4.099d-6*tempf**3)/C_REF
c      ccc2(2,2,iph)=ccc2(1,1,iph)
c      ccc2(1,2,iph)=(61357+49.009*tempf-4.1198d-2*tempf**2
c     #           +1.396d-5*tempf**3)/C_REF
c      ccc2(2,1,iph)=ccc2(1,2,iph)
c      ccc2(1,3,iph)=(64912+1.8131d-2*tempf+4.4831d-3*tempf**2
c     #           -3.704d-6*tempf**3)/C_REF
c      ccc2(2,3,iph)=ccc2(1,3,iph)
c      ccc2(3,1,iph)=ccc2(1,3,iph)
c      ccc2(3,2,iph)=ccc2(1,3,iph)
c     ccc2(3,3,iph)=(174080-30.996*tempf-1.3754d-3*tempf**2
c     #           -6.997d-9*tempf**3)/C_REF
c      ccc2(4,4,iph)=(37290-20.79*tempf+1.1433d-2*tempf**2
c     #           -5.0173d-6*tempf**3)/C_REF
c      ccc2(5,5,iph)=ccc2(4,4,iph)
c      ccc2(6,6,iph)=(ccc2(1,1,iph)-ccc2(1,2,iph))/2
      ccc2=0.0
      alfacc=0.0
      do iph=1,nph
c *** Elastic coefficients for URANIUM                             ***
          ccc2(1,1,iph)=(2.15753+9.73266e-05*tempf-3.52983e-07*tempf**2
     #           -2.56109e-10*tempf**3)/C_REF
          ccc2(2,2,iph)=(2.16106-0.000675579*tempf+4.44855e-07*tempf**2 
     #           -4.8623e-10*tempf**3)/C_REF
          ccc2(3,3,iph)=(2.9554-0.000890293*tempf-1.63953e-07*tempf**2
     #           -1.12126e-10*tempf**3)/C_REF
          ccc2(4,4,iph)=(1.39999-0.00039389*tempf-4.67451e-07*tempf**2
     #           +1.21138e-10*tempf**3)/C_REF
          ccc2(5,5,iph)=(0.897582-0.000356345*tempf-7.55658e-07*tempf**2
     #           +3.64044e-10*tempf**3)/C_REF
          ccc2(6,6,iph)=(0.895228-0.00060743*tempf+4.68913e-07*tempf**2
     #           -4.99329e-10*tempf**3)/C_REF
          ccc2(1,2,iph)=(0.37151+0.00047945*tempf-7.18435e-07*tempf**2
     #           +5.51519e-10*tempf**3)/C_REF
          ccc2(1,3,iph)=(0.227897-0.000125654*tempf+3.03661e-07*tempf**2
     #           +6.28745e-13*tempf**3)/C_REF
          ccc2(2,3,iph)=(0.744436+0.00280592*tempf-7.95823e-06*tempf**2
     #           +8.667e-09*tempf**3-3.33214e-12*tempf**4)/C_REF
          ccc2(2,1,iph)=ccc2(1,2,iph)
          ccc2(3,1,iph)=ccc2(1,3,iph)
          ccc2(3,2,iph)=ccc2(2,3,iph)
      
c *** Thermal expansion coefficients for Uranium  (Lloyd & Barrett)               ***
c      alfacc(1,iph)=(24.226e-6)-(9.832e-9)*(tempf)
c     #         +(46.021e-12)*(tempf**2)
c      alfacc(2,iph)=3.078e-6+(3.476e-9)*(tempf)
c     #        -(38.451e-12)*(tempf**2)
c      alfacc(3,iph)=8.721e-6+(37.043e-9)*(tempf)
c     #        +(9.081e-12)*(tempf**2)
      
c *** Thermal expansion coefficients for Uranium  (Toul)               ***
          alfacc(1,iph)=(28.6e-6)-(44.1e-9)*(tempf)+(83.7e-12)
     #        *(tempf**2)
          alfacc(2,iph)=(1.37e-6)+(5.37e-9)*(tempf)-(31.39e-12)
     #        *(tempf**2)
          alfacc(3,iph)=(3.46e-6)+(50.7e-9)*(tempf)-(5.71e-12)
     #        *(tempf**2)

          call invten(ccc2(:,:,iph),scc2(:,:,iph))
      
      enddo
      
      return
      END SUBROUTINE
      
      end module mphase_props
c      
c***********************************************************************
c *** grain_props_v ****************************************************
c***********************************************************************
      module grain_props
      
      use grain_props_v
      
      contains
      !cr_to_sa     :Transforms from Crystal to Sample system 
      !data_sample  :Reads sample data                        
     
      SUBROUTINE cr_to_sa(ng1,ng2,iopt)
c **********************************************************************
c *** Rotates MATERIAL properties from crystal to sample for grains  ***
c *** ng1 to ng2                                                     ***
c *** Option: 0 Rotates Schmid tensors & thermal coefficients        ***
c ***         1 Rotates the stiffness and compliance                 ***
c **********************************************************************
c *** USES:    voigt                                                 ***
c **********************************************************************
c *** VERSION: 12/aug/98                                             ***
c **********************************************************************
c

      use mvoigt
      use mphase_props
      use flags, only : inonSch
      
c ______________________________________________________________________
c
      DIMENSION ccc4(3,3,3,3),scc4(3,3,3,3),alfacc2(3,3)
c ______________________________________________________________________
c
c *** Rotates Schmid tensor, thermal tensor, Burgers vector from CR-->SA
      if (iopt.eq.0) then
c
c *** Rotation of Schmid tensor
        do 10 ng=ng1,ng2
          do 10 ns=1,nsys(ngrnph(ng))
            do 10 ij=1,6
              mcs(ij,ns,ng)=0.0
              if(inonSch.eq.1) nmcs(ij,ns,ng)=0.0 !inonSch
              qcs(ij,ns,ng)=0.0
              i=ijv(ij,1)
              j=ijv(ij,2)
              do 10 i1=1,3
                do 10 j1=1,3
                  mcs(ij,ns,ng)=mcs(ij,ns,ng)
     #                          +r(i1,i,ng)*r(j1,j,ng)
     #                          *mc2(i1,j1,ns,ngrnph(ng))
                  if(inonSch.eq.1) then !inonSch
                      nmcs(ij,ns,ng)=nmcs(ij,ns,ng)
     #                          +r(i1,i,ng)*r(j1,j,ng)
     #                          *nmc2(i1,j1,ns,ngrnph(ng))
                  endif
                  qcs(ij,ns,ng)=qcs(ij,ns,ng)
     #                          +r(i1,i,ng)*r(j1,j,ng)
     #                          *qc2(i1,j1,ns,ngrnph(ng)) 
   10   continue
c
c *** Burgers vector BCS(i,ns,ng) required for twinning reorientation
        do 15 ng=ng1,ng2
          do 15 ns=1,nsys(ngrnph(ng))
            do 15 i=1,3
              bcs(i,ns,ng)=0.0
                ncs(i,ns,ng)=0.0
              do 15 j=1,3
                bcs(i,ns,ng)=bcs(i,ns,ng)+r(j,i,ng)*bcc(j,ns,ngrnph(ng))
                ncs(i,ns,ng)=ncs(i,ns,ng)+r(j,i,ng)*ncc(j,ns,ngrnph(ng))
   15   continue
c
c *** Go from 6-vector to 3x3 tensor and rotate thermal expansion tensor
        call voigt(alfacc(:,ngrnph(ng1)),alfacc2,dummy3,dummy4,1)
        do 20 ng=ng1,ng2
          do 20 ij=1,6
            i=ijv(ij,1)
            j=ijv(ij,2)
            alfacs(ij,ng)=0.0
            do 20 i1=1,3
              do 20 j1=1,3
                alfacs(ij,ng)=alfacs(ij,ng)
     #                        +r(i1,i,ng)*r(j1,j,ng)*alfacc2(i1,j1)
   20   continue
c
      else
c
c *** Go from (6x6) matrix to (3x3x3x3) tensor and rotate elastic moduli
        call voigt(dummy1,dummy2,ccc2(:,:,ngrnph(ng1)),ccc4,3) !MZ_2ph ccc2 -> extra dimension
        call voigt(dummy1,dummy2,scc2(:,:,ngrnph(ng1)),scc4,3) !MZ_2ph scc2 -> extra dimension
        do 30 ng=ng1,ng2
          do 30 ij=1,6
            i=ijv(ij,1)
            j=ijv(ij,2)
            do 30 kl=1,6
              k=ijv(kl,1)
              l=ijv(kl,2)
              ccs2(ij,kl,ng)=0.0
              scs2(ij,kl,ng)=0.0
              do 30 i1=1,3
                do 30 j1=1,3
                  do 30 k1=1,3
                    do 30 l1=1,3
                       ccs2(ij,kl,ng)=ccs2(ij,kl,ng)+r(i1,i,ng)
     #                                *r(j1,j,ng)*r(k1,k,ng)
     #                                *r(l1,l,ng)*ccc4(i1,j1,k1,l1)
                       scs2(ij,kl,ng)=scs2(ij,kl,ng)+r(i1,i,ng)
     #                                *r(j1,j,ng)*r(k1,k,ng)
     #                                *r(l1,l,ng)*scc4(i1,j1,k1,l1)
   30   continue
c
      endif
c ______________________________________________________________________
c
      return
      END SUBROUTINE
     
      SUBROUTINE data_sample(filesamp,ngrain,iph) !MZ_2ph
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     SUBROUTINE data_sample      -->      version JUNE/06/2010
c **********************************************************************
c *** Reads from the file "filesamp" the SAMPLE data                 ***
c *** Renormalizes weights. Calculates rotation matrix for each grain **
c *** Reads state from previous procedure if i_prev_proc=1           ***
c **********************************************************************
c *** USES:    euler                                                 ***
c **********************************************************************
      
      use miscellaneous_sub, only : euler
      use mvoigt
      use mphase_props
      use sample_props_v, only : nph,ngParent
      use flags, only : ishape,iOutput
      
      integer, intent(in) :: iph
      integer, intent (out) :: ngrain

      CHARACTER*78 prosa,filesamp*150,fileprev*150,eul_conv*1

      DIMENSION aux33(3,3),aux(3,3),bur(3),hd(NSLS,NSLS)

c ______________________________________________________________________

    1 FORMAT(a)
    2 FORMAT(1h ,78('*'))
    3 FORMAT(1h ,8('*'),' SAMPLE data - File: ',a,9('*'))
    4 FORMAT(1h ,3f7.2)
    9 FORMAT(1h ,i5,4f12.5)
   10 FORMAT(a1,2i5)
   11 FORMAT(1h ,a1,2i5)
   12 FORMAT(1h ,a)
c ______________________________________________________________________
c *** Opens the file with the SAMPLE data "filesamp"

      OPEN(unit=1,file=filesamp,status='old')


      read(1,1) prosa
      read(1,1) prosa
      read(1,1) prosa
      read(1,*) eul_conv,ngrain
      if (ngrain.gt.NGR) then
        write(*,'(1h ,''ERROR: Number of grains greater than code''
     #,'' dimension !!!'',/,1h ,''DIMENSION in code = '',i4)') NGR
        write(*,*)
        write(*,'(1h ,''STOP IN ROUTINE *** data_sample ***'')')
        stop
      endif
        ngParent=ngrain
        nphngr(iph)=ngrain !MZ_2ph define array nphngr 
      if(iOutput.eq.1) write(12,11) eul_conv,ngrain
      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
      ng2=SUM(nphngr(1:iph))
      read(1,*) (phi(i),the(i),ome(i),wgt(i),
     #     i=ng1,ng2) !MZ_2ph

c     print *
c     print *, 'RANDOM SHIFT OF EULER ANGLES WITHIN 3 degs ACTIVATED'
c     print *, 'TO AVOID FULLY SYMMETRIC ORIENTATIONS'
c     print *
c     jran=-1
c      do i=1,ngrain
c        phi(i)= phi(i)+ran2(jran)*3.
c        the(i)= the(i)+ran2(jran)*3.
c        ome(i)= ome(i)+ran2(jran)*3.
c      enddo
c     write(12,9) (i,phi(i),the(i),ome(i),wgt(i),i=1,ngrain)
c
      CLOSE(UNIT=1)  
      
c      if (iph.eq.nph) then !MZ_2ph only when it reads all orientations it should calculate weight and rotation matrix etc.
      !MZ_2ph correct ngrain and ngParent to grain num rom all phases          
c          ngrain=SUM(nphngr)
c          ngParent=ngrain            
          
c ______________________________________________________________________
c *** Calculates rotation matrix for each grain:
c *** matrix r(i,j,ngr) transforms from sample to crystal
c         totwgt=0.0
c         do ng=1,ngrain
c           totwgt=totwgt+wgt(ng)
c         enddo
         ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
         ng2=SUM(nphngr(1:iph))
         do ng=ng1,ng2
           ngrnph(ng)=iph
           call euler(2,phi(ng),the(ng),ome(ng),aux33)
           wgt(ng)=wgt(ng) !/totwgt
           do i=1,3
             do j=1,3
               r(i,j,ng)=aux33(i,j)
             enddo
           enddo
         
         enddo

         CLOSE(unit=1)

c      endif !MZ_2ph
      
      return
      END SUBROUTINE
      
      SUBROUTINE update_grain_orientatin(STEP,ngrain,omegag)
      
      use miscellaneous_sub, only : reorient_grain
      
      !in
      integer, intent(in) :: ngrain
      real, intent(in) :: step,omegag(3,3,NGR)
      
      DIMENSION ROT(3,3),AROTG(3,3),rg(3,3)
      
      
      DO NG=1,NGRAIN
C *** CREATE A ROTATION MATRIX BASED UPON THE RESULT

        DO I=1,3
          DO J=1,3
            ROT(I,J)=(omegag(i,j,ng))*STEP
          ENDDO
        ENDDO

        CALL REORIENT_GRAIN (AROTG,ROT)
c       CALL RODRIGUES (ROT,AROTG)

        DO I=1,3
          DO J=1,3
            rg(j,i)=0.0
            DO K=1,3
              rg(j,i)=rg(j,i)+AROTG(I,K)*r(j,k,ng)
            END DO
          END DO
        END DO
        r(:,:,ng)=rg(:,:)
      END DO
        
      RETURN
      END SUBROUTINE  
      
      subroutine add_grains(iph,str)
      use mphase_props, only : nphngr
      character(len=150) :: str
      
      !read texture file
      call data_sample (str,ngrain,iph)     
      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
      ng2=SUM(nphngr(1:iph))
      !populate crystal properties 
      if (ng1.le.ng2) then
        call cr_to_sa(ng1,ng2,0)
        call cr_to_sa(ng1,ng2,1)    
      endif     
      
      end subroutine add_grains
      
      end module grain_props
c      
c***********************************************************************
c *** mphase ***********************************************************
c***********************************************************************
      module mphase_state
      
      use mphase_props_v
      use mphase_state_v
           
      contains
      !update_phase_shape    :updates average ellipsoid for Eshelby calc    
      !update_phase_fij      :updates deformation state for sub update_shape
      
      SUBROUTINE initialize_phase_shape(iph)
      
      use const
      use flags
      use miscellaneous_sub, only : euler

      DIMENSION fijx(3,3),fnew(3,3),aaa(3,3)
c ______________________________________________________________________
c   Initialize values for fijph --> CT: need to move to appropriate subr
c_______________________________________________________________________
      da=eulerph(1,iph)
      db=eulerph(2,iph)
      dc=eulerph(3,iph)    

      DO I=1,3
        AXISPH(0,I,iph)=axis(I,iph) !MZ_pseudo set up array defining ellipsoid per phase
      END DO

      call euler(2,da,db,dc,aaa)

      do i=1,3
      do j=1,3
        fijx(i,j)=(i/j)*(j/i)*AXISPH(0,I,iph)
      enddo
      enddo

      do j=1,3
      do i=1,3
        fnew(i,j)=0.
        do m=1,3
          fnew(i,j)=fnew(i,j)+aaa(m,i)*fijx(m,j)
        enddo
      enddo
      enddo

      do i=1,3
      do j=1,3
        fijph(i,j,iph)=fnew(i,j) !MZ_pseudo set up fijph per phase since shape of ellipsoids can be different
      enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE update_phase_shape(iph)

      !Copied from VPSC6 in abbreviated/modified form by JN, 8-7-07
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      use miscellaneous_sub, only : jacobi,eigsrt,det,euler


      DIMENSION W(3),BX(3,3),B(3,3),BT(3,3)


C *** IF iPH=0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN ELEMENT
C *** IF iPH>0 ELLIPSOID REPRESENTS AVERAGE DEFORMATION IN PHASE 'IPH'
C *** CALCULATES EIGENVALUES, EIGENVECTORS & EULER ANGLES OF ELEMENT GRAIN,
C     PHASE GRAIN, OR INDIVIDUAL GRAIN
C *** 'AXISPH' TRANSFORMS FROM ELLIPSOID TO SAMPLE AXES.
      DO I=1,3
      DO J=1,3
        BX(I,J)=0.
        DO K=1,3
          BX(I,J)=BX(I,J)+FIJPH(I,K,iph)*FIJPH(J,K,iph) 
        ENDDO
      ENDDO
      ENDDO 

      CALL JACOBI(BX,3,3,W,B,NROT,IER)
      CALL EIGSRT(W,B,3,3)
      IF (IER.EQ.1) THEN
        WRITE(*,*) 'ERROR IN UPDATE_SHAPE FOR PHASE ELLIPSOID'
        STOP
      ENDIF

C *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
C *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR.
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING 1 AND 2.

      SIGN=-1.
      IF(DET(B).LE.0.) SIGN=1.
      DO I=1,3
        EXCHANGE=B(I,1)
        B(I,1)=B(I,2)
        B(I,2)=EXCHANGE*SIGN
      ENDDO
      EXCHANGE=W(1)
      W(1)=W(2)
      W(2)=EXCHANGE

      DO I=1,3
        AXISPH(0,I,iph)=SQRT(W(I))
        DO J=1,3
          AXISPH(I,J,iph)=B(I,J)
          BT(I,J)=B(J,I)
        ENDDO
      ENDDO

      CALL EULER(1,ANG1,ANG2,ANG3,BT)
      EULERPH(1,iph)=ANG1
      EULERPH(2,iph)=ANG2
      EULERPH(3,iph)=ANG3    

      end subroutine

      SUBROUTINE update_phase_fij(iph,step,etrss_ph)
C *****************************************************************************
C     SUBROUTINE UPDATE_FIJ      --->      VERSION OF NOV/28/2005
C
C     USES THE VELOCITY GRADIENT (AVERAGE, PHASE or GRAIN) IN THE STEP
C     TO UPDATE INCREMENTALLY THE CORRESPONDING DEFORMATION TENSOR 'FIJ'
C
C     !Copied from VPSC6, highly abbreviated and modified by JN, 8-7-07
C *****************************************************************************
      use mvoigt
      use const
      
      !in
      real, intent(in) :: etrss_ph(6,NPHM)

      DIMENSION FNEW(3,3)
      DIMENSION XID3(3,3), etrss_33(3,3), etrcs6(6), etrcs_33(3,3)

C *** UPDATES THE DEFORM GRAD IN THE ELEMENT 'FIJPH(i,j)' USING THE
C     MACROSCOPIC DEFORMATION STEP 'DEL_ETSS_33(I,J)'


      DO I=1,3
        DO J=1,3
          XID3(I,J)=0
          IF (I.EQ.J) THEN
            XID3(I,J)=1
          END IF
        END DO
      END DO

      CALL VOIGT(etrss_ph(:,iph),etrss_33,dummy3,dummy4,1) !MZ_pseudo etrss_ph(:,iph) instead of etrss

      DO I=1,3
        DO J=1,3
          FNEW(I,J)=0.0
          DO K=1,3
            FNEW(I,J)=FNEW(I,J)+(etrss_33(I,K)*step+XID3(I,K))
     #      *FIJPH(K,J,iph)
          ENDDO
        ENDDO
      ENDDO
      DO I=1,3
        DO J=1,3
          FIJPH(I,J,iph)=FNEW(I,J) !MZ_pseudo fijph goes over phase and etrss_33 should go over phase (average value over phase grains of etrcs)
        ENDDO
      ENDDO

      END SUBROUTINE
  
      SUBROUTINE update_phase_etrss(etrcs,etrss_ph)
      
      use const
      use mphase_props, only : wgt_ph,ngrnph
      use sample_props_v, only : ngrain,nph
      use grain_props_v, only : wgt
      
      !in 
      real, intent(in) :: etrcs(6,NGR)
      !out
      real, intent(out) :: etrss_ph(6,NPHM)
      
      etrss_ph(:,:)=0.0

      do ng=1,ngrain
        do i=1,6
          etrss_ph(i,ngrnph(ng))=etrss_ph(i,ngrnph(ng))+
     #            etrcs(i,ng)*wgt(ng)
        enddo
      enddo
      do iph=1,nph
        do i=1,6
          etrss_ph(i,iph)=etrss_ph(i,iph)/wgt_ph(iph)
        enddo
        if(wgt_ph(iph).eq.0.0) then !in case phase has 0 weight don't update ellipsoid
          etrss_ph(1:6,iph)=0.0
        endif
      enddo
      !check if <etrcs>=etrss
c      do i=1,6
c          etrss_ph(i,3)=etrss_ph(i,1)*0.5+etrss_ph(i,2)*0.5
c      enddo
      
      RETURN
      END SUBROUTINE       
     
      subroutine update_phase_state(iph)
      
      use grain_state_v, only : etcs,stcs
      use grain_props_v, only : wgt
      
      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
      ng2=SUM(nphngr(1:iph))
      wgt_ph(iph)=SUM(wgt(ng1:ng2)) 
      
      etss_ph(:,iph)=0.0
      stss_ph(:,iph)=0.0
      do ng=ng1,ng2
          etss_ph(:,iph)=etss_ph(:,iph)+etcs(:,ng)*wgt(ng)
          stss_ph(:,iph)=stss_ph(:,iph)+stcs(:,ng)*wgt(ng)
      enddo
      etss_ph(:,iph)=etss_ph(:,iph)/wgt_ph(iph)
      stss_ph(:,iph)=stss_ph(:,iph)/wgt_ph(iph)
      
      end subroutine
     
      end module mphase_state
c
c***********************************************************************
c *** grain_state ******************************************************
c***********************************************************************
      module grain_state

      use grain_state_v
     
      contains
      !update_grain_shape    :updates average ellipsoid for Eshelby calc    
      !update_grain_fij      :updates deformation state for sub update_shape
      
      SUBROUTINE initialize_grain_shape(iph)
      
      use flags
      use miscellaneous_sub, only : euler
      use mphase_props, only : nphngr
      use mphase_state, only : fijph,axis
      
      !in
      integer, intent(in) :: iph

      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
      ng2=SUM(nphngr(1:iph))      
      
      DO I=1,3
        do j=ng1,ng2
          axisgr(0,i,j)=axis(i,iph)
        end do
      END DO
      
      do k=ng1,ng2
        do i=1,3
          do j=1,3
            fijgr(i,j,k)=fijph(i,j,iph)
          end do
        enddo
      end do

c________________________________________________________________________
      
      RETURN
      END SUBROUTINE

      SUBROUTINE update_grain_gam(step,gamd)
      
      use const
      use sample_props_v      
      use grain_props_v
      use mphase_props      
      
      !in
      real, intent(in) :: gamd(NSLS,NGR)
      do ng=1,ngrain
        if(nact(ng).ne.0) then            
          do ns1=1,nsys(ngrnph(ng))
            gamtot(ng) =gamtot(ng) +gamd(ns1,ng)*step
          enddo
        endif
      enddo      
        
      RETURN
      END SUBROUTINE
      
      SUBROUTINE update_grain_tau(step,taud)
      
      use const
      use sample_props_v      
      use grain_props_v
      use mphase_props      
      
      !in
      real, intent(in) :: taud(NSLS,NGR)
      
      !Updates CRSS's, Gamma's in systems and grains
      do ng=1,ngrain
          if(nact(ng).ne.0) then            
            do ns1=1,nsys(ngrnph(ng))
	    	 tau(ns1,ng)=tau(ns1,ng)+taud(ns1,ng)*step
            enddo
          endif
      enddo 

      RETURN
      END SUBROUTINE
      
      SUBROUTINE update_grain_state(step,deltemp,etrcs,strcs,etrcs_eig)
      
      use grain_props_v, only : alfacs
      use sample_props_v, only : ngrain        
      !in
      real, intent(in) :: strcs(6,NGR),etrcs(6,NGR),etrcs_eig(6,NGR)
      
      do i=1,6
        do ng=1,ngrain
          stcs(i,ng)=stcs(i,ng)+strcs(i,ng)*step
          etcs(i,ng)=etcs(i,ng)+etrcs(i,ng)*step
          etcs_eig(i,ng)=etcs_eig(i,ng)+etrcs_eig(i,ng)*step
          etthcs(i,ng)=etthcs(i,ng)+alfacs(i,ng)*deltemp*step
        enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE update_grain_shape(ng1,ng2)
      
      use flags
      use sample_props_v      
      use grain_props_v
      use mphase_props   
      use miscellaneous_sub, only : jacobi,eigsrt,det

      !Copied from VPSC6 in abbreviated/modified form by JN, 8-7-07
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      DIMENSION W(3),BX(3,3),B(3,3),BT(3,3),fijgr_tot(3,3)
      
      do igr=ng1,ng2
        if (iPhTr.eq.1) then
          fijgr_tot=matmul(FIJGR(:,:,igr),fijgr_pt(:,:,igr))
        else 
          fijgr_tot=FIJGR(:,:,igr)
        endif
        DO I=1,3
          DO J=1,3
            BX(I,J)=0.
            DO K=1,3
              BX(I,J)=BX(I,J)+fijgr_tot(I,K)*fijgr_tot(J,K)
            ENDDO
          ENDDO
        ENDDO

        CALL JACOBI(BX,3,3,W,B,NROT,IER)
        CALL EIGSRT(W,B,3,3)
        IF (IER.EQ.1) THEN
          WRITE(*,*) 'Error in update_shape for grain',igr
          STOP
        ENDIF

C *** EIGENVALUES (AND ASSOC EIGENVECTORS) ARE ORDERED FROM LARGER TO SMALLER.
C *** REDEFINE AXIS(2) TO BE THE LARGEST IN ORDER TO IMPROVE ACCURACY IN THE
C     CALCULATION OF THE ESHELBY TENSOR.
C *** IF DET(B)<0 MEANS THAT THE SYSTEM IS LEFT HANDED. IT IS MADE RIGHT
C     HANDED BY EXCHANGING 1 AND 2.

        SIGN=-1.
        IF(DET(B).LE.0.) SIGN=1.
        DO I=1,3
          EXCHANGE=B(I,1)
          B(I,1)=B(I,2)
          B(I,2)=EXCHANGE*SIGN
        ENDDO
        EXCHANGE=W(1)
        W(1)=W(2)
        W(2)=EXCHANGE


        DO I=1,3
          AXISGR(0,I,IGR)=SQRT(W(I))
          DO J=1,3
            AXISGR(I,J,IGR)=B(I,J)
          ENDDO
        ENDDO


      end do

c      write (*,*) (axisph(0,i),i=1,3)
c      do ig=1,5
c        write (*,*) ig
c        write (*,*) (axisgr(0,i,ig),i=1,3)
c      enddo
c      read (*,*)

      end subroutine

      SUBROUTINE UPDATE_GRAIN_FIJ(step,etrcs)
      
C *****************************************************************************
C     SUBROUTINE UPDATE_FIJ      --->      VERSION OF NOV/28/2005
C
C     USES THE VELOCITY GRADIENT (AVERAGE, PHASE or GRAIN) IN THE STEP
C     TO UPDATE INCREMENTALLY THE CORRESPONDING DEFORMATION TENSOR 'FIJ'
C
C     !Copied from VPSC6, highly abbreviated and modified by JN, 8-7-07
C *****************************************************************************

      use sample_props_v      
      use grain_props_v
      use mphase_props
      use mvoigt
      
      !in
      real, intent(in) :: etrcs(6,NGR)

      DIMENSION FNEW(3,3)
      DIMENSION XID3(3,3), etrss_33(3,3), etrcs6(6), etrcs_33(3,3)

C *** UPDATES THE DEFORM GRAD IN THE ELEMENT 'FIJPH(i,j)' USING THE
C     MACROSCOPIC DEFORMATION STEP 'DEL_ETSS_33(I,J)'
        do igr=1,ngrain
          
          DO I=1,3
            DO J=1,3
              XID3(I,J)=0
              IF (I.EQ.J) THEN
                XID3(I,J)=1
              END IF
            END DO
          END DO        
        
          do i=1,6
            etrcs6(i)=etrcs(i,igr)
          end do

          CALL VOIGT(etrcs6,etrcs_33,dummy3,dummy4,1)

          DO I=1,3
            DO J=1,3
              FNEW(I,J)=0.0
              DO K=1,3
                FNEW(I,J)=FNEW(I,J)+(etrcs_33(I,K)*step+XID3(I,K))
     #          *FIJGR(K,J,IGR)
              ENDDO
            ENDDO
          ENDDO
          DO I=1,3
            DO J=1,3
              FIJGR(I,J,IGR)=FNEW(I,J)
            ENDDO
          ENDDO
        enddo

      END SUBROUTINE

      SUBROUTINE update_grain_corot_st(ngrain,step,spin_cor_s,omegag) !rot_ferr
C **********************************************************************
C     SUBROUTINE COROTATION
C
C     co-rotational stress calculations
C **********************************************************************

      use mvoigt
      use grain_props_v, only : wgt
      
      !in 
      integer, intent(in) :: ngrain
      real, intent(in) :: step,omegag(3,3,NGR)
      
      !out
      real, intent(out) :: spin_cor_s(3,3) !MZ vol_avg_spin 

      DIMENSION stcs6(6),stcs33(3,3),tmp33(3,3),tmp6(6)
      DIMENSION stss6(6),stss33(3,3),tmps33(3,3),tmps6(6)
      DIMENSION aux66(6,6),aux3333(3,3,3,3)


          spin_cor_s(:,:) = 0.0    !MZ vol_avg_spin 		
          DO NG=1,ngrain !rot_ferr
            DO i=1,6
              stcs6(i)=stcs(i,ng)           !take the stress for the crystal
            ENDDO
            call voigt(stcs6,stcs33,aux66,aux3333,1)     !convert to 3x3 tensor
            DO I=1,3
              DO J=1,3
                tmp33(i,j)=0.0
                DO k=1,3
                  tmp33(i,j)=tmp33(i,j)+(omegag(i,k,ng)*stcs33(k,j)   !find rate terms specific to Jaumann derivative
     #                       -omegag(k,j,ng)*stcs33(i,k))*step
                END DO
              END DO
            END DO
            call voigt(tmp6,tmp33,aux66,aux3333,2)    !convert to 6 vector
            do i=1,6
              stcs(i,ng)=stcs(i,ng)+tmp6(i)           !update the stress in the grain to keep it on the yield surface
            enddo
            !MZ vol_avg_spin 
            spin_cor_s = spin_cor_s + tmp33*wgt(ng)			
          enddo

      END SUBROUTINE
       
      SUBROUTINE update_grain_config(ng1,ng2)
      ! updates crystal variables to current config at t+dt
      
      use mvoigt
      use grain_props_v, only : r,wgt
      use sample_props_v, only : nph
      use mphase_props_v, only : nphngr,wgt_ph
      use grain_props, only : cr_to_sa
      use grain_rate_v, only : drotcs,omegag
      
      DIMENSION aux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3)
      
      !update rotation matrix 
      do ng=ng1,ng2
        r(:,:,ng)=MATMUL(r(:,:,ng),transpose(drotcs(:,:,ng)))
      enddo
          
      !update all variables (or include phase division !???)
      do iph=1,nph
          ng1ph=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2ph=SUM(nphngr(1:iph))
          wgt_ph(iph)=SUM(wgt(ng1ph:ng2ph))          
          if (ng1ph.le.ng2ph) then
              call cr_to_sa(ng1ph,ng2ph,0)
              call cr_to_sa(ng1ph,ng2ph,1)
          endif
      enddo
      
      !update total stress (crystal)
      do ng=ng1,ng2
        !use drotcs
        call VOIGT(stcs(:,ng),aux33,aux66,aux3333,1)
        aux33=MATMUL(aux33,transpose(drotcs(:,:,ng)))
        aux33=MATMUL(drotcs(:,:,ng),aux33)
        call VOIGT(stcs(:,ng),aux33,aux66,aux3333,2)
        !use spin
c        call VOIGT(stcs(:,ng),aux33,aux66,aux3333,1)
c        aux33=MATMUL(aux33,omegag(:,:,ng))-MATMUL(omegag(:,:,ng),aux33)
c        call VOIGT(aux6,aux33,aux66,aux3333,2)
c        stcs(:,ng)=stcs(:,ng)+aux6
      enddo
      
      !update total strain (crystal)
      do ng=ng1,ng2
        call VOIGT(etcs(:,ng),aux33,aux66,aux3333,1)
        aux33=MATMUL(aux33,transpose(drotcs(:,:,ng)))
        aux33=MATMUL(drotcs(:,:,ng),aux33)
        call VOIGT(etcs(:,ng),aux33,aux66,aux3333,2)
      enddo     
      
      
      !update total eigen strain (crystal)
      do ng=ng1,ng2
        call VOIGT(etcs_eig(:,ng),aux33,aux66,aux3333,1)
        aux33=MATMUL(aux33,transpose(drotcs(:,:,ng)))
        aux33=MATMUL(drotcs(:,:,ng),aux33)
        call VOIGT(etcs_eig(:,ng),aux33,aux66,aux3333,2)
      enddo        

      RETURN
      END SUBROUTINE
      
      end module grain_state
c      
c***********************************************************************
c *** sample_state *****************************************************
c***********************************************************************
      module sample_state
     
      use sample_state_v
      
      contains
      
      SUBROUTINE update_sample_state(etrss,strss)
      
      use mvoigt, only:profac
      use sample_props_v
      
      real, intent(in) :: etrss(6),strss(6)
      
      !update total stress and strain
      stss=stss+strss
      etss=etss+etrss
      
      !update elastic strain
      do i=1,6
        etelss(i)=0.0
        do j=1,6
          etelss(i)=etelss(i)+sss2(i,j)*stss(j)*profac(j)
        enddo
      enddo   
      
      RETURN
      END SUBROUTINE     
      
      SUBROUTINE update_sample_corot_st(omegabcr,spin_cor_s)
      
      use mvoigt
      use flags, only : iSingleCry
      
      !in
      real, intent(in) :: omegabcr(3,3),spin_cor_s(3,3)      

      DIMENSION stss6(6),stss33(3,3),tmps33(3,3),tmps6(6)
      DIMENSION aux66(6,6),aux3333(3,3,3,3)
      
c *** Enforcing Rigid body rotations on sample. Added by JAW 10/09/09
      DO i=1,6
        stss6(i)=stss(i)
      ENDDO
      call voigt(stss6,stss33,aux66,aux3333,1)
      DO I=1,3
        DO J=1,3
          tmps33(i,j)=0.0
          DO k=1,3
            tmps33(i,j)=tmps33(i,j)+(omegabcr(i,k)*stss33(k,j)
     #                 -omegabcr(k,j)*stss33(i,k))*step
          END DO
        END DO
      END DO
c      call voigt(tmps6,tmps33,aux66,aux3333,2) !MZ vol_avg_spin - commented old spin correction
      call voigt(tmps6,spin_cor_s,aux66,aux3333,2) !MZ vol_avg_spin - use new spin correction
      do i=1,6
        stss(i)=stss(i)+tmps6(i)
      enddo
c *** End of Additions by JAW
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE g_average(etrss,strss,etrcs,strcs,etrav)
c **********************************************************************
c *** Calculates the averages and deviations for the stress rate,the ***
c *** strain rate, the stress, the strain and the pressure.          ***
c **********************************************************************
c *** VERSION: 24/APR/97
c **********************************************************************
c
      use mvoigt, only : profac
      use sample_props_v
      use grain_props_v
      use grain_state
      use flags, only : iOutput
      use grain_rate_v, only : acs2
      use sample_rate_v, only : ass2
      use miscellaneous_sub, only : tmismatch

      !in
      real, intent(in) :: etrss(6),strss(6),etrcs(6,NGR),strcs(6,NGR) 
      !out
      real, intent(out) :: etrav(6)

c
      DIMENSION strav(6),strdev(6),etrdev(6),stdev(6)
      DIMENSION etdev(6),stPGav(6),EtPGav(6),stCHav(6),etCHav(6)
c ______________________________________________________________________
c
    1 FORMAT(1h ,6e12.4)
    2 FORMAT(1h ,25e12.4)
    3 FORMAT(1h ,6e12.4)
c ______________________________________________________________________
c
      presss =(stss(1)+stss(2)+stss(3))/3.0
      presav =0.0
      presdev=0.0
        do ng=1,ngrain
          prescs =(stcs(1,ng)+stcs(2,ng)+stcs(3,ng))/3.0
          presav =presav +prescs   *wgt(ng)
          presdev=presdev+prescs**2*wgt(ng)
        enddo
      presdev=sqrt(abs(presdev-presav**2))
c
      do i=1,6
        strav(i)=0.0
        etrav(i)=0.0
        stav(i)=0.0
        etav(i)=0.0
        stPGav(i)=0.0
        etPGav(i)=0.0
        stCHav(i)=0.0
        etCHav(i)=0.0
        wgtPA=0.0
        wgtCH=0.0
        do ng=1,ngrain
          strav(i)=strav(i)+strcs(i,ng)*wgt(ng)
          etrav(i)=etrav(i)+etrcs(i,ng)*wgt(ng)
          stav(i)=stav(i)+stcs(i,ng)*wgt(ng)
          etav(i)=etav(i)+etcs(i,ng)*wgt(ng)
        enddo
        do ng=1,ngParent
          stPGav(i)=stPGav(i)+stcs(i,ng)*wgt(ng)
          etPGav(i)=etPGav(i)+etcs(i,ng)*wgt(ng)
        enddo
        do ng=ngParent+1,ngrain
          stCHav(i)=stCHav(i)+stcs(i,ng)*wgt(ng)
          etCHav(i)=etCHav(i)+etcs(i,ng)*wgt(ng)
        enddo
      enddo
      wgtPA=0.0
      do ng=1,ngParent
        wgtPA=wgtPA+wgt(ng)
      enddo
      if(wgtPA.gt.0) then
        do i=1,6
          stPGav(i)=stPGav(i)/wgtPA
          etPGav(i)=etPGav(i)/wgtPA
        enddo
      endif
      wgtCH=0.0
      do ng=ngParent+1,ngrain
        wgtCH=wgtCH+wgt(ng)
      enddo
      if(wgtCH.gt.0) then
        do i=1,6
          stCHav(i)=stCHav(i)/wgtCH
          etCHav(i)=etCHav(i)/wgtCH
        enddo
      endif
c
      strnorm=0.0
      etrnorm=0.0
      stnorm =0.0
      etnorm =0.0
      do i=1,6
        strnorm=strnorm+strav(i)**2 * profac(i)
        etrnorm=etrnorm+etrav(i)**2 * profac(i)
        stnorm =stnorm +stav(i)**2  * profac(i)
        etnorm =etnorm +etav(i)**2  * profac(i)
      enddo
      presdev=presdev/sqrt(stnorm)
      do i=1,6
        strdev(i)=0.0
        etrdev(i)=0.0
        stdev(i)=0.0
        etdev(i)=0.0
        do ng=1,ngrain
          strdev(i)=strdev(i)+(strcs(i,ng)-strav(i))**2*wgt(ng)
          etrdev(i)=etrdev(i)+(etrcs(i,ng)-etrav(i))**2*wgt(ng)
          stdev(i)=stdev(i)+(stcs(i,ng)-stav(i))**2*wgt(ng)
          etdev(i)=etdev(i)+(etcs(i,ng)-etav(i))**2*wgt(ng)
        enddo
      enddo
      do i=1,6
        strdev(i)=sqrt(strdev(i))/sqrt(strnorm)
        etrdev(i)=sqrt(etrdev(i))/sqrt(etrnorm)
        stdev(i) =sqrt(stdev(i))  /sqrt(stnorm)
        etdev(i) =sqrt(etdev(i))  /sqrt(etnorm)
      enddo
      aux=0.0
      do ng=1,ngrain
        aux=aux+tmismatch(acs2(:,:,ng),ass2,6,6)*wgt(ng)
      enddo
      
      if(iOutput.eq.1) then
c
        write(11,*)
        write(11,*)'Bound. Cond., Av. and Dev. STRESS RATE (normalized)'
        write(11,1) strss
        write(11,1) strav
        write(11,3) strdev
        write(11,*)'Bound. Cond., Av. and Dev. STRAIN RATE (normalized)'
        write(11,1) etrss
        write(11,1) etrav
        write(11,3) etrdev
        write(11,*)
        write(11,*) 'Bound. Cond., Av. and Dev. STRESS (normalized)'
        write(11,1) stss
        write(11,1) stav
        write(11,3) stdev
        write(11,1) stPGav
        write(11,1) stCHav
        write(11,*) 'Bound. Cond., Av. and Dev. STRAIN (normalized)'
        write(11,1) etss
        write(11,1) etav
        write(11,3) etdev
        write(11,1) etPGav
        write(11,1) etCHav
        write(11,*) 'Bound. pressure, Av. Press. and Dev. (normalized)'
        write(11,1) presss,presav,presdev
c       
        write(15,2) etrss,etrdev,etss,etdev,aux
        write(16,2) strss,strdev,stss,stdev,aux
      endif
c ______________________________________________________________________
c
      return
      END SUBROUTINE
      
      SUBROUTINE update_sample_config
      ! updates sample variables to current config at t+dt
      
      use mvoigt
      use sample_rate_v, only : drotss
      use bc_v, only : omegabcr
      
      DIMENSION aux6(6),aux33(3,3),aux66(6,6),aux3333(3,3,3,3)
      
      !update total stress (macroscopic)      
      !use drotss
      call VOIGT(stss,aux33,aux66,aux3333,1)
      aux33=MATMUL(aux33,transpose(drotss))
      aux33=MATMUL(drotss,aux33)
      call VOIGT(stss,aux33,aux66,aux3333,2)
      !use spin
c      call VOIGT(stss,aux33,aux66,aux3333,1)
c      aux33=MATMUL(aux33,omegabcr)-MATMUL(omegabcr,aux33)
c      call VOIGT(aux6,aux33,aux66,aux3333,2)
c      stss=stss+aux6
      
      !update total strain (macroscopic)
      call VOIGT(etss,aux33,aux66,aux3333,1)
      aux33=MATMUL(aux33,transpose(drotss))
      aux33=MATMUL(drotss,aux33)
      call VOIGT(etss,aux33,aux66,aux3333,2)        
      
      RETURN
      END SUBROUTINE
      
      end module sample_state
c      
c***********************************************************************
c *** mphase_rate ******************************************************
c***********************************************************************
      module mphase_rate
      
      use mphase_props   
      use mphase_rate_v
           
      !update_phase_shape    :updates average ellipsoid for Eshelby calc    
      !update_phase_fij      :updates deformation state for sub update_shape
     
      end module mphase_rate
c
c***********************************************************************
c *** grain_rate *******************************************************
c***********************************************************************
      module grain_rate

      use grain_rate_v
     
      contains

      SUBROUTINE update_grain_spin(step,omegabcr,etrss,escr4,einvsa) !added by JN 8-28-2007 !rot_ferr
c
      use miscellaneous_sub, only : reorient_grain
      use flags, only : ishape,iSingleCry
      use mphase_props
      use grain_props_v
      use sample_props_v, only : ngrain
      use mvoigt
      
      !in
      real, intent(in) :: step,omegabcr(3,3),etrss(6)
     #    ,ESCR4(3,3,3,3,NPHM),EINVSA(3,3,3,3,NPHM)

c
      DIMENSION dnsa(3,nsls),dbsa(3,nsls),rotslip(3,3),AROTG(3,3),
     #          ROTLOC(3,3),rg(3,3),rot(3,3),DEV(6),DEV33(3,3)
c     #          AS(3,3,3,3)
      REAL      LIJGR0(3,3)

C *** FOR EVERY GRAIN...
      DO NG=1,ngrain
C *** FINDS CRYSTALLOGRAPHIC ROTATION COMPONENT DUE TO LOCAL ROTATION
C *** ROTLOC=PI*S**(-1)*(DG-DAV)

        DO I=1,3
        DO J=1,3
        DO K=1,3
        DO L=1,3
          DUMMY=0.0
          DO K1=1,3
          DO L1=1,3
            if (ishape.eq.0.or.ishape.eq.1) then
              DUMMY=DUMMY+ESCR4(I,J,K1,L1,ngrnph(ng))
     #           *EINVSA(K1,L1,K,L,ngrnph(ng))!MZ_pseudo EINVSA and ESCR4 are per phase and ROTLOC will depend on phase as well but its temporary 
            ELSEIF (ishape.ge.2) THEN
              DUMMY=DUMMY+ESCR4GR(I,J,K1,L1,NG)*EINVSAGR(K1,L1,K,L,NG)
            end if
          ENDDO
          ENDDO
          AS(I,J,K,L,ng)=DUMMY
        ENDDO
        ENDDO
        ENDDO
        ENDDO

        DO I=1,6
          DEV(I)=ETRCS(I,NG)-ETRSS(I)
        END DO

        CALL VOIGT(DEV,DEV33,dummy3,dummy4,1)

        DO I=1,3
          DO J=1,3
          ROTLOC(I,J)=0.0
            DO K=1,3
              DO L=1,3
                ROTLOC(I,J)=ROTLOC(I,J)+AS(I,J,K,L,ng)*DEV33(K,L)
              END DO
            END DO
          END DO
        END DO

      ! CPFE      
        if (iSingleCry.eq.1) then
            ROTLOC=0.0
        endif  
            
C *** FINDS CRYSTALLOGRAPHIC ROTATION COMPONENT DUE TO PLASTIC SLIP ACTIVITY

        do i=1,3
          do j=1,3
            rg(j,i)=r(j,i,NG)
            LIJGR0(I,J)=0.0
          enddo
        enddo

        do is=1,nsys(ngrnph(ng))
          do i=1,3
            dnsa(i,is)=0.0
            dbsa(i,is)=0.0
            do j=1,3
              dnsa(i,is)=dnsa(i,is)+rg(j,i)*ncc(j,is,ngrnph(ng))
              dbsa(i,is)=dbsa(i,is)+rg(j,i)*bcc(j,is,ngrnph(ng))
            enddo
          enddo
          do i=1,3
            do j=1,3
           LIJGR0(i,j)=LIJGR0(i,j)+dbsa(i,is)*dnsa(j,is)*gamd(is,ng)
            enddo
          enddo
        enddo

        DO I=1,3
          DO J=1,3
            ROTSLIP(I,J)=(LIJGR0(I,J)-LIJGR0(J,I))/2.0
          ENDDO
        ENDDO

C *** SUMMING OF THE CRYSTALLOGRAPHIC GRAIN ROTATION COMPONENTS (RIGID minus PLASTIC)
c *** ROT and omega equations changed to include rigid rotation (omegabcr) by JAW 10/09/09
c *** 'omegag' is used inside COROTATION to calculate Jauman derivative
        DO I=1,3
          DO J=1,3
            ROT(I,J)=(omegabcr(i,j)+ROTLOC(I,J)-ROTSLIP(I,J))*STEP
            omegag(i,j,ng)=omegabcr(i,j)+rotloc(i,j)-rotslip(i,j)
          ENDDO
        ENDDO
c *** End of Changes by JAW
      ENDDO

      END SUBROUTINE
      
      end module grain_rate
c      
c***********************************************************************
c *** bc ***************************************************************
c***********************************************************************
      module bc

      use bc_v
      
      contains
      !data_process :Reads process data
      !load_conditions :interprets boundary conditions and finds   
      !                 symmetric and antisymmetric imposed strains      

      SUBROUTINE data_process (fileproc,i_temp_cij,i_ref_et,i_ref_st,
     #                         i_bc_mode,ivarBC)
     
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     SUBROUTINE DATA_PROCESS      -->      VERSION JUN/06/2010
c
c *** Reads from the file "fileproc" the PROCESS data                ***
c **********************************************************************   
      !out      
      integer, intent(out) :: i_temp_cij,i_ref_et,i_ref_st,i_bc_mode

      CHARACTER*78 prosa,fileproc*150
c ______________________________________________________________________
c
    1 FORMAT(a)
    2 FORMAT(1h ,78('*'))
    3 FORMAT(1h ,8('*'),' PROCESS data - File: ',a,8('*'))
    4 FORMAT(1h ,6(6x,i6))
    5 FORMAT(1h ,6d12.4)
    6 FORMAT(1h ,'ITMAX OVER MODULUS= ',i3,' ERROR= ',d12.4)
    7 FORMAT(1h ,'ITMAX OVER GRAINS= ',i3)
    8 FORMAT(1h ,d12.4)
   12 FORMAT(1h ,a)
c ______________________________________________________________________

c *** Opens the file with the PROCESS data "fileproc"                ***
      OPEN(unit=1,file=fileproc,status='old')

      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,*) nsteps
      READ(1,*) i_control_var
      READ(1,*) i_bc_mode
      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,1) PROSA

c *** Reading 9 flags for full strain tensor. Added by JAW 10/9/09
      READ(1,*) ifulletbc(1,1),ifulletbc(1,2),ifulletbc(1,3)
      READ(1,*) ifulletbc(2,1),ifulletbc(2,2),ifulletbc(2,3)
      READ(1,*) ifulletbc(3,1),ifulletbc(3,2),ifulletbc(3,3)
c *** Updating old 6 component flags
      ietbc(1)=ifulletbc(1,1)
      ietbc(2)=ifulletbc(2,2)
      ietbc(3)=ifulletbc(3,3)
      ietbc(4)=ifulletbc(2,3)*ifulletbc(3,2)
      ietbc(5)=ifulletbc(1,3)*ifulletbc(3,1)
      ietbc(6)=ifulletbc(1,2)*ifulletbc(2,1)

      READ(1,1) PROSA

c *** Reading full strain tensor. Added by JAW 10/9/09
      READ(1,*) fulletbc(1,1),fulletbc(1,2),fulletbc(1,3)
      READ(1,*) fulletbc(2,1),fulletbc(2,2),fulletbc(2,3)
      READ(1,*) fulletbc(3,1),fulletbc(3,2),fulletbc(3,3)

      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,*) istbc(1),istbc(6),istbc(5)
      READ(1,*) istbc(2),istbc(4)
      READ(1,*) istbc(3)
      READ(1,1) PROSA
      READ(1,*) stbc(1),stbc(6),stbc(5)
      READ(1,*) stbc(2),stbc(4)
      READ(1,*) stbc(3)
      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,*) temp_s
      READ(1,1) PROSA
      READ(1,*) deltemp
      READ(1,1) PROSA
      READ(1,*) i_temp_cij
      READ(1,1) PROSA
      READ(1,*) i_ref_et
      READ(1,1) PROSA
      READ(1,*) i_ref_st
      READ(1,1) PROSA
      READ(1,1) PROSA
      READ(1,*) edot_macro

      if (ivarBC.eq.1) then
          call varBC_read      
      endif
      
      CLOSE(unit=1)
      
      return
      END SUBROUTINE

      SUBROUTINE load_conditions


      DIMENSION tempa(6)

C     Added 8-9-07 by JN in order to allow for imposing stress control or strain
c     control using off-diagonal boundary condition elements.
c
c     i_control_var is fed in from main routine - values 0-6
c     icvx is returned - values 0-12

C     Checks imposed load conditions for the process to make sure they are
c     properly defined, and determines wether to impose stress or strain control
c     for the process based upon the ietbc(6) and istbc(6) arrays.  Allows
c     control based on diagonal (i=j) or off diagonal stress or strain components
c     of the symmetric stress or strain control tensors, etbc(6) and stbc(6).

c *** commented out by JAW 10/09/09
c      DO I=1,6
c        IF (ietbc(I).NE.0.AND.ietbc(I).NE.1) THEN
c          write (*,*)
c          WRITE (*,*) 'ERROR! IETBC COMPONENT ',I,' CAN ONLY BE 0 OR 1.'
c          STOP
c        END IF
c      END DO
      DO I=1,3
        DO J=1,3
          IF (ifulletbc(I,J).NE.0.AND.ifulletbc(I,J).NE.1) THEN
            write (*,*)
            WRITE (*,*) 'ERROR! IFBC COMPONENT (',I,',',J,
     #        ') CAN ONLY BE 0 OR 1.'
            STOP
          END IF
        END DO
      END DO
c *** End of edit by JAW
      DO I=1,6
        IF (istbc(I).NE.0.AND.istbc(I).NE.1) THEN
          write (*,*)
          WRITE (*,*) 'ERROR! ISTBC COMPONENT ',I,' CAN ONLY BE 0 OR 1.'
          STOP
        END IF
      END DO

C *** Checks for under-constraint or over-constraint of boundary conditions.

      DO I=1,6
        TEMPA(I)=ietbc(I)+ISTBC(I)
        IF (TEMPA(I).EQ.2) THEN
          write (*,*)
          WRITE (*,*) 'CANNOT CONSTRAIN STRAIN COMPONENT (',I,')'
          WRITE (*,*) 'AND STRESS COMPONENT (',I,') - MUST RELAX ONE.'
          STOP
        END IF
        IF (TEMPA(I).EQ.0) THEN
          write (*,*)
          WRITE (*,*) 'CANNOT RELAX STRAIN COMPONENT (',I,') AND'
          WRITE (*,*) 'STRESS COMPONENT (',I,') - MUST CONSTRAIN ONE.'
          STOP
        END IF
      END DO
c *** Full strain tensor checks Added by JAW 10/9/09 Commented out by BJCL 2010/02/04
c      do i=1,2
c      do j=i+1,3
c        if(ifulletbc(i,j)+ifulletbc(j,i).eq.0) then
c          write(*,*) 'CHECK OFF-DIAGONAL STRAIN BOUNDARY CONDITIONS'
c          write(*,*) 'CANNOT RELAX BOTH OFF-DIAGONAL COMPONENTS'
c          stop
c        ELSEIF (ifulletbc(i,j).ne.ifulletbc(j,i)) Then
c          write(*,*) 'CHECK OFF-DIAGONAL STRIAN BOUNDARY CONDITIONS'
c          write(*,*) 'YOUR LAB FRAME IS ROTATING'
c          stop
c        endif
c      enddo
c      enddo
c *** End of Additions by JAW
      !Determines stress or strain control component - looks for imposed BC in the specified control slot.

      IF (i_control_var.eq.0) THEN                 !Temperature Control if true
        icvx=0
      ELSE
        IF (ietbc(i_control_var).eq.1) THEN        !Incremental Strain Control if true
          icvx=i_control_var
        ELSE IF (istbc(i_control_var).eq.1) THEN   !Incremental Stress Control if true
          icvx=i_control_var+6
        ELSE
          WRITE (*,*)
          WRITE (*,*) 'ERROR-CANNOT FIND CONTROLLING BOUNDARY CONDITION'
          WRITE (*,*) 'FOR ICTRL =',I_CONTROL_VAR
          STOP
        END IF
      END IF

      END SUBROUTINE
      
      SUBROUTINE define_bc
      
      use mvoigt
      
      !Added by JAW 10/9/09
      do i=1,3
        do j=1,3
          etbc_sym(i,j)=0.5*(fulletbc(i,j)+fulletbc(j,i))
          omegabc(i,j)=0.5*(fulletbc(i,j)-fulletbc(j,i))
        end do
      end do
      call voigt(etbc,etbc_sym,dummy3,dummy4,2)
      !End of addition by JAW      
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE define_bcr(etss,stss,i_bc_mode,ivarBC)
      
      integer, intent(in) :: i_bc_mode
      real, intent(in) :: etss(6),stss(6)
      
      !Relative bc (mz)
      if(i_bc_mode.eq.0) then          
          do i=1,6
            strbc(i)=0.0
            if (istbc(i).eq.1) strbc(i)=stbc(i)/nsteps
            etrbc(i)=0.0
            if (ietbc(i).eq.1)  etrbc(i)=etbc(i)/nsteps
          enddo
      endif
      !Absolute bc (mz)
      if(i_bc_mode.eq.1) then
          do i=1,6
            strbc(i)=0.0
            if (istbc(i).eq.1) strbc(i)=(stbc(i)-stss(i))/nsteps
            etrbc(i)=0.0
            if (ietbc(i).eq.1)  etrbc(i)=(etbc(i)-etss(i))/nsteps
          enddo
      endif
      
      !omega bcr
      do i=1,3
        do j=1,3
          omegabcr(i,j)=omegabc(i,j)/nsteps
           if (omegabcr(i,j).gt.0) then
            write (*,*)
            WRITE (*,*) 'YOU HAVE IMPOSED A SPIN!!'
            write (*,*)
           endif
        end do
      end do      

      if (ivarBC.eq.0) then
        !Aadded to prevent process control by a zero component.
        if (icvx.ge.1.and.icvx.le.6) then
          if (etrbc(icvx).eq.0.0) then
            write (*,*) "WARNING - Strain increment for strain control"
            write (*,*) "component should not be zero."
            stop
          end if
        end if
        if (icvx.ge.7) then
          if (strbc(icvx-6).eq.0.0) then
            write (*,*) "WARNING - Stress increment for stress control"
            write (*,*) "component should not be zero."
            stop
          end if
        end if
      endif

      !Defines Incremental Temperature Step Consideration
      temp=temp_s
      temp_f=temp_s+deltemp*nsteps 
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE correct_bcr(i_bc_mode,istep,etss_ini
     #                       ,etss,stss_ini,stss)
      
      integer, intent(in) :: i_bc_mode
      real, intent(in) :: etss_ini(6),etss(6),stss_ini(6),stss(6)
      
      do i=1,6
        strbc(i)=0.0
        if(istbc(i).eq.1) then
          if(i_bc_mode.eq.1) then
            strbc(i)=(stbc(i)-stss_ini(i))*(istep+1)/nsteps
     #               +stss_ini(i)-stss(i)
          else
            strbc(i)=(stbc(i)-0.0)*(istep+1)/nsteps
     #               +stss_ini(i)-stss(i)
          endif
        endif
        etrbc(i)=0.0
        if(ietbc(i).eq.1) then
          if(i_bc_mode.eq.1) then
            etrbc(i)=(etbc(i)-etss_ini(i))*(istep+1)/nsteps
     #               +etss_ini(i)-etss(i)
          else
            etrbc(i)=(etbc(i)-0.0)*(istep+1)/nsteps
     #               +etss_ini(i)-etss(i)
          endif
        endif
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE et_st_ref(i_ref_et,i_ref_st,etssref)
      
      use sample_state, only : etss
      use grain_state, only : etthcs,stcsref,stcs
      use grain_props_v, only : wgt
      
      real, intent(out) ::etssref(6)
      
      integer, intent(in) :: i_ref_et,i_ref_st
      
      real  aux6(6)
      
      !If i_ref_et=1 resets total strain in sample and crystals to zero.
      !This is done for redefining the strain origin when plotting, but
      !does not alter the stress states or the elastic strains.
      !Total and grain stress-strains are updated incrementally at each
      !step using the strain rates calculated in the step.
      !'stcsref' introduced on 15/07/98 permits to refer 'elastic' strain
      !in grains to the values of the previous process (for comparing with
      !X-ray or neutron measurements). Done inside SUBROUTINE DIF_PLANES.
      !
      ! Merkel, 05/2010, removed initialization of etcsref: not used anywhere in the code
      if (i_ref_et.eq.1) then
        do i=1,6
          etssref(i)=etss(i)
        enddo
c        do ng=1,ngrain
c          do i=1,6
c            etcsref(i,ng)=etcs(i,ng)
c          enddo
c        enddo
      else if (i_ref_et.eq.-1) then
        do i=1,6
          etssref(i)=0.0
        enddo
c        do ng=1,ngrain
c          do i=1,6
c            etcsref(i,ng)=0.0
c          enddo
c        enddo
      endif
      if (i_ref_st.eq.1) then
        do ng=1,ngrain
          do i=1,6
            stcsref(i,ng)=stcs(i,ng)
          enddo
        enddo
      else if (i_ref_st.eq.-1) then
        do ng=1,ngrain
          do i=1,6
            stcsref(i,ng)=0.0
          enddo
        enddo
      endif
        do i=1,6
        aux6(i)=0.0
        do ng=1,ngrain
            aux6(i)=aux6(i)+etthcs(i,ng)*wgt(ng)
          enddo
        enddo
 
      RETURN
      END SUBROUTINE
      
      subroutine add_bc(fileproc,i_temp_cij,i_ref_et,i_ref_st,i_bc_mode
     #     ,etss,stss,etssref)
      
      !in
      character(len=150), intent(in) :: fileproc
      real, intent(in) :: stss(6),etss(6)
      !out
      integer, intent(out) :: i_temp_cij,i_ref_et,i_ref_st,i_bc_mode
      real, intent(out) :: etssref(6)
      
      !Calls routine for input of thermo-mechanical process conditions.
      call data_process(fileproc,i_temp_cij,
     #     i_ref_et,i_ref_st,i_bc_mode,0)

      call load_conditions
      
      !Calculates symmetric and antisymmetric portions of strain tensor
      call define_bc
      call define_bcr(etss,stss,i_bc_mode,0)

      !reference state for macro and grain state prior to beginning of process
      call et_st_ref(i_ref_et,i_ref_st,etssref)
      
      end subroutine add_bc
      
      subroutine varBC_read
      
      real :: varBC_all(100000,6)
      
      !read file
      read(1,*) (ietbc_var(i),i=1,6)
      read(1,*) (istbc_var(i),i=1,6)
      read(1,*) nproc
      do i=1,nproc
          read(1,*) varBC_all(i,1:sum(ietbc_var)+sum(istbc_var))
      enddo
      
      !Calculates symmetric and antisymmetric portions of strain tensor
      call define_bc
      
      !populate etbc_var and stbc_var
      k=1
      do i=1,6
          if (ietbc_var(i).eq.1) then
              etbc_var(i,:)=varBC_all(:,k)
              k=k+1
          elseif (istbc_var(i).eq.1) then
              stbc_var(i,:)=varBC_all(:,k)
              k=k+1
          endif
      enddo
      
c      do i=1,6
c          if (ietbc_var(i).eq.0.and.ietbc(i).eq.1) then
c              etbc_var(i,:)=etbc
c          elseif (istbc_var(i).eq.0.and.istbc(i).eq.1) then
c              stbc_var(i,:)=stbc
c          endif
c      enddo
      
      end subroutine
      
      subroutine add_bc_var(fileproc,iproc,i_temp_cij,i_ref_et,i_ref_st
     #     ,i_bc_mode,etss,stss,etssref)
      
      !in
      integer, intent(in) :: iproc
      character(len=150), intent(in) :: fileproc
      real, intent(in) :: stss(6),etss(6)
      !out
      integer, intent(out) :: i_temp_cij,i_ref_et,i_ref_st,i_bc_mode
      real, intent(out) :: etssref(6)
      
      real :: aux61(6),det_full(3,3),aux66(6,6),aux3333(3,3,3,3)
     #    ,aux31(3,3),drot_full(3,3)
      
      !Calls routine for input of thermo-mechanical process conditions.
      if (iproc.eq.1) then
          call data_process(fileproc,i_temp_cij,
     #        i_ref_et,i_ref_st,i_bc_mode,1)
      endif
      
      etbc=etbc_var(:,iproc)
      stbc=stbc_var(:,iproc)
      
      call load_conditions
      
      !Calculates symmetric and antisymmetric portions of strain tensor
      call define_bcr(etss,stss,i_bc_mode,1)

      !reference state for macro and grain state prior to beginning of process
      call et_st_ref(i_ref_et,i_ref_st,etssref)
      
      end subroutine add_bc_var
      
      end module bc
c      
c***********************************************************************
c *** diffract *********************************************************
c***********************************************************************
      module diffract
      
      use diffract_v
      
      contains
      !dif_planes   :Calculates residual strain along a direction
      
      SUBROUTINE dif_planes(filediff,temp,istep,iopt,iph)
c **********************************************************************
c *** iopt=0 - Read dif file and crystal values for hkl sets of interest
c
c *** iopt=1 - Determine grains in each hkl set
c            - calculates lattice strains and intensity for the set
c **********************************************************************
c *** VERSION 3/aug/08
c **********************************************************************
      use mvoigt, only : ijv,profac
      
      use const
      use flags    
      use twinning_v
      use mphase_props
      use mphase_state
      use sample_props_v
      use sample_state
      use grain_props_v
      use grain_state
      use bc

      CHARACTER prosa*78,filediff*150
      LOGICAL   heading

      DIMENSION pc(3,24),ps(3),etelcsx(6),n4(4),aux(3),rg(3,3)
      dimension isn(4),sn(3),sneq(3,24),isb(4),sb(3)                   !CNT
      DIMENSION para_w(NDIFFX),ipol(6,3)
c      DIMENSION vc(ndiffx,3,24,NPHM),nfamily(NDIFFX,NPHM)
c      DIMENSION vs(NDIFFX,3,NPHM)
      DIMENSION para_sq(NDIFFX)
c Merkel 05/2010, pressure in the grain, changes in hydrostatic strain tensors
      REAL pressure(6) , p, etelhycsx(6), etelcsy(6), p_incr(6)
c Merkel 05/2010, extra parameters for deviatoric strains
      DIMENSION para_w_dev(NDIFFX), para_sq_dev(NDIFFX)
c
      pi=4.0*atan(1.0)
c
c **********************************************************************
    1 FORMAT(a)
    2 FORMAT(1h ,a)
    3 FORMAT(1h ,i5,f6.2)
    4 FORMAT(1h ,4i3,2(3x,f8.2))
    5 FORMAT(1h ,3i3,3x,2(3x,f8.2))
    6 FORMAT(1h ,'  SET #:',i6,'  NGRSET:',i6,'  VOLFRAC:',f10.6)
    7 FORMAT(1h ,15i6)
    8 FORMAT(1h ,3I6,7E13.5,e15.5,3F12.2,F12.8)

      !define grain range belonging to phase iph
      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
      ng2=SUM(nphngr(1:iph))    
    
c **********************************************************************
c *** Reads indices of diffracting plane and angles of diffraction direction.
c *** Generates crystallographically equivalent plane normals 'pc'.
c *** Calculates diffraction direction 'ps'
      if (iopt.eq.0) then
        OPEN(unit=1,file=filediff,status='old')

        iunit=190+iph
        read(1,1) prosa
        if (iOutput.eq.1) write(iunit,2) prosa
        read(1,1) prosa
        if (iOutput.eq.1) write(iunit,2) prosa
        read(1,*) ndiff(iph),angdetector
        if (iOutput.eq.1) write(iunit,*) ndiff(iph),angdetector
        if (ndiff(iph).gt.NDIFFX) then
          write(*,'(1h ,''ERROR: Number of diffraction directions''
     #    '' greater than code dimension !!!'',/,1h ,''DIMENSION in''
     #    '' code = '',i3)') NDIFFX
          write(*,*)
          write(*,'(1h ,''STOP IN ROUTINE *** dif_planes ***'')')
          stop
        endif
        read(1,1) prosa
        if (iOutput.eq.1) write(iunit,2) prosa
        read(1,1) prosa
        if (iOutput.eq.1) write(iunit,2) prosa
        prosa=prosa      ! only to fool the compiler

        nind=3                                     !CNT
        if(icrysym(iph).eq.2 .or. icrysym(iph).eq.3) nind=4

        do n=1,ndiff(iph)

          read(1,*) (isn(i),i=1,nind),chi,eta                          !CNT
          if(nind.eq.3.and.iOutput.eq.1) write(iunit,'(3i4,2f10.1)') 
     #                    (isn(i),i=1,3),chi,eta
          if(nind.eq.4.and.iOutput.eq.1) write(iunit,'(4i4,2f10.1)') 
     #                    (isn(i),i=1,4),chi,eta
          eta=eta*pi/180.
          chi=chi*pi/180.
          ps(1)=cos(eta)*sin(chi)
          ps(2)=sin(eta)*sin(chi)
          ps(3)=         cos(chi)
          do i=1,nind
            n4(i)=isn(i)
          enddo

          call crystal_symmetry (3,1,icrysym(iph),isn,sn,pc,isb,sb
     #        ,nfamilyx)!CNT
          nfamily(n,iph)=nfamilyx
          toler=cos(angdetector*pi/180.0)
          rand_wgt(n)=nfamilyx*(1-toler)
c          rand_wgt(n)=rand_frac(n,nfamilyx*2,angdetector,angpole)
c          WRITE(21,*) N,'  ',nfamilyx*2,'   ',RAND_WGT(N)

c         write(19,'('' EQUIVALENT PLANES IN SET'',i5)') n
c         write(19,'(3(3x,3f7.3))') ((pc(i,nf),i=1,3),nf=1,nfamily)
c
          do i=1,3
            vs(n,i,iph)=ps(i)
            do j=1,nfamily(n,iph)
              vc(n,i,j,iph)=pc(i,j)
            end do
          enddo
        enddo
c		
      endif   !end of initialization run for iopt = 0
c
c *********************************************************************
c *** If twinning or rotations due to slip are active, find grains
c *** contributing to each diffraction set, ndiff
c *** OR, if this is the initialization step and no twinning or rotations
c *** will be occurring, then find the sets of grains this time only.

      isw=0
      if (iopt.eq.1.and.itwinning+irot.ge.1.or.nph.ne.1) isw=1  ! convoluted print control
      if (iopt.eq.0.and.itwinning+irot.eq.0) isw=1  ! improve it...

      if (isw.eq.1) then
        do n=1,ndiff(iph)
          wgtset(n)=0.0
            ngrset(n)=0

          do ng=ng1,ng2
            do ipl=1,nfamily(n,iph)
              do i=1,3
                ps(i)=0.0
                do j=1,3
                  ps(i)=ps(i)+r(j,i,ng)*vc(n,j,ipl,iph)
                enddo
              enddo

              prodesc=0.0
              do i=1,3
                prodesc=prodesc+ps(i)*vs(n,i,iph)
              enddo
               if (abs(prodesc).ge.toler) then
                ngrset(n)=ngrset(n)+1
                igrset(n,ngrset(n))=ng
                wgtset(n)=wgtset(n)+wgt(ng)
              endif
           enddo
          enddo
        enddo      ! end of do n=1,ndiff
c **********************************************************************
        if (itwinning+irot.ne.0) then
          do n=1,ndiff(iph)
            write(iunit,*)
            write(iunit,6) n,ngrset(n),wgtset(n)
            if (iwrite9.eq.1) then
              write(iunit,7) (igrset(n,n1),n1=1,ngrset(n))
            endif
          enddo
        endif

        heading=.true.
        close(unit=1)
      endif

      if (iopt.eq.1.and.itwinning+irot.eq.0.and.nph.eq.1) then
        do n=1,ndiff(iph)
          write(iunit,*)
          write(iunit,6) n,ngrset(n),wgtset(n)
          if (iwrite9.eq.1) then
            write(iunit,7) (igrset(n,n1),n1=1,ngrset(n))
          endif
        enddo
        heading=.true.
      endif

c **********************************************************************
c *** calculates average strain for a crystallographic family along
c *** the diffraction direction.

      if (iopt.eq.1) then

c Merkel 05/2010: need to keep trak of the total elastic strain in the grain.
c    Need to update the tensor in all grains, even if not seen in diffraction
c    Need to update the hydrostatic strain tensor in each grain

        if (kSM.eq.1) then

          do ng=ng1, ng2
            do i=1,6
              etelcsx(i)=0.0
              etelhycsx(i)=0.0
              p = (stcs(1,ng) + stcs(2,ng) + stcs(3,ng))/3.
              pref = (stcsref(1,ng) + stcsref(2,ng) + stcsref(3,ng))/3.
              p_incr = (/ p-pref, p-pref, p-pref, 0., 0., 0. /)
              do j=1,6
                etelhycsx(i)=etelhycsx(i)+scs2(i,j,ng)*
     #               (p_incr(j))*profac(j)
                etelcsx(i)=etelcsx(i)+scs2(i,j,ng)*
     #               (stcs(j,ng)-stcsref(j,ng))*profac(j)
              enddo
              etelcs(i,ng) = (1.+etelcs(i,ng))*etelcsx(i) + etelcs(i,ng)
              etelhycs(i,ng)=(1.+etelhycs(i,ng))*etelhycsx(i)
     #                          +etelhycs(i,ng)
            enddo
          enddo
        endif

        do n=1,ndiff(iph)
c Merkel, 05/2010: we now have two strains, the full elastic strain and
c   the elastic strain relative to the compressed state
          para_w(n)=0.0
          para_sq(n)=0.0
          para_w_dev(n)=0.0
          para_sq_dev(n)=0.0
          do ng=1,ngrset(n)
            igset=igrset(n,ng)
            eps=0.0
            eps_dev=0.0
            if (kSM.eq.1) then
              tmpEpsDev1 = 0.
              tmpEpsDev2 = 0.
              do ij=1,6
                i=ijv(ij,1)
                j=ijv(ij,2)
                tmpEpsDev1 = tmpEpsDev1 + vs(n,i,iph)*vs(n,j,iph)*
     #             (etelcs(ij,igset)-etelhycs(ij,igset))*profac(ij)
                tmpEpsDev2 = tmpEpsDev2 + vs(n,i,iph)*vs(n,j,iph)*
     #                       etelhycs(ij,igset)*profac(ij)

                eps=eps+vs(n,i,iph)*vs(n,j,iph)*etelcs(ij,igset)
     #                       *profac(ij) !calculates elastic strain of crystal along Q if a nearby hkl qualifies
              enddo
              eps_dev =  tmpEpsDev1 / (1.+tmpEpsDev2)
c Regular case, no large elastic strains
            else
              do i=1,6
                etelcsx(i)=0.0
                do j=1,6
                  etelcsx(i)=etelcsx(i)+scs2(i,j,igset)*
     #                 (stcs(j,igset)-stcsref(j,igset))*profac(j)
                enddo
c Comment the line below out to exclude thermal strains from diff output data
c                etelcsx(i)=etelcsx(i)+etthcs(i,igset)
              enddo
              do ij=1,6
                i=ijv(ij,1)
                j=ijv(ij,2)
c Calculates elastic strain of crystal along Q for given hkl
                eps=eps+vs(n,i,iph)*vs(n,j,iph)*etelcsx(ij)*profac(ij)
              enddo
            endif
            para_w(n)=para_w(n)+eps*wgt(igset)
            para_sq(n)=para_sq(n)+eps*eps*wgt(igset)
c Merkel, 05/2010, we now have two strains, the full elastic strain and the
c   elastic strain relative to the compressed state
            para_w_dev(n)=para_w_dev(n)+eps_dev*wgt(igset)
            para_sq_dev(n)=para_sq_dev(n)+eps_dev*eps_dev*wgt(igset)
          enddo      ! end of DO loop over grains in subset

          if(wgtset(n).ne.0.0) then
            para_w(n)=para_w(n)/wgtset(n)
            tmp = para_sq(n)/wgtset(n)-para_w(n)*para_w(n)
            if (abs(tmp).lt.1.E-10) tmp = 0
            para_sq(n)=sqrt(tmp)
c Merkel, 05/2010, we now have two strains, the full elastic strain and the elastic strain
c relative to the compressed state
            para_w_dev(n)=para_w_dev(n)/wgtset(n)
            tmp = para_sq_dev(n)/wgtset(n)-
     #                               para_w_dev(n)*para_w_dev(n)
            if (abs(tmp).lt.1.E-10) tmp = 0
            para_sq_dev(n)=sqrt(tmp)
          endif

        enddo      ! END OF LOOP OVER DIFFRACTING PLANES

c Merkel, 05/2010: Reset stress references for each grain
        if (kSM.eq.1) then

          do ng=ng1,ng2
            do i=1,6
               stcsref(i,ng) = stcs(i,ng)
            enddo
          enddo
        endif

        !calculate diffracting elastic constants
c        call diff_el_const


c **************************************************************************
c *** Write the results in "EPSC9.OUT" -> unit=19

        if (icvx.ge.7) then
          eref=etss(icvx-6)-etssref(icvx-6)*1e2
          sref=stss(icvx-6)
        else if (icvx.lt.7 .and. icvx.ge.1) then
          eref=etss(icvx)-etssref(icvx)*1e2
          sref=stss(icvx)
          else
            eref=temp
        endif

        if (heading) then
        write(iunit,*)
          write(iunit,'('' Temp      Eref      Sref   '',
     #    ''       dif1       dif2       dif3       dif4       dif5'',
     #    ''       dif6       dif7       dif8       dif9      dif10'')')

          heading=.false.
        endif

c **********************************************************************
c *** Standard output for DEC least squares routine.
c        if(itwinning.eq.1) then
          iunit=200+iph
          write(iunit,'(1400E13.5)')
     #      eref,
     #      ((etss(i)-etssref(i))*1e2,i=1,3),
     #      (stss(i)*1e0,i=1,3),
     #      (etav(i),i=1,3),
     #      (stav(i),i=1,3),
     #      (para_w(nd)*1e6,nd=1,ndiff(iph)),
c     #      (real(ngrset(nd)),nd=1,ndiff),
     #      (para_sq(nd)*1e6,nd=1,ndiff(iph)),
     #      (para_w_dev(nd)*1e6,nd=1,ndiff(iph)),
     #      (para_sq_dev(nd)*1e6,nd=1,ndiff(iph)),
     #      (wgtset(nd)/rand_wgt(nd),nd=1,ndiff(iph)),
     #      (dec(1,1,nd),nd=1,ndiff(iph)),
     #      (dec(2,2,nd),nd=1,ndiff(iph))
c        endif
        iunit=190+iph
        write(iunit,'(1400E13.5)') temp,
     #    ((etss(i)-etssref(i)),i=1,3),
     #    (stss(i),i=1,3),
     #    (para_w(nd)*1e6,nd=1,ndiff(iph))
c        WRITE(23,'(I6,I6,2F8.3,2(6E12.4,3X),200E12.4,3X)')
c     #    istep,ngrain,actav,temp,
c     #    ((etss(i)-etssref(i)),i=1,6),
c     #    (stss(i),i=1,6),(para_w(nd)*1e6,nd=1,ndiff),
c     #    ((wgtset(nd)/rand_wgt(nd)),nd=1,ndiff)
        
      endif            !END OF IOPTION=1


      RETURN
      END SUBROUTINE

      SUBROUTINE diff_el_const
      !Main ref:"Microstresses in textured polycrystals studied by multireflection
      !diffraction method and self-consistent model"
      !Summury: 
      !Due to inhomogenous plastic deformation residual stresses are present after unloading
      !of material. These can be macro (sample level) and micro (grain lelvel). Average residual
      !strains coming from these stresses of subset of specifically oriented grains 
      !(diffracting grains) can be measured (XRD).
      !<et_hkl> = F*st_s + q*<Q*Q*S_c*st_c>
      !F - diffracting elastic constants
      !S_c - crystal complience in sample frame
      !Q - transformation matrix from sample to laboratory (where residual strains are measured)
      !<Q*Q*S_c*st_c> - this should be calculated lattice strain by EPSC
      !This equation is rewritten in terms of interplanar spacing and fitted with experimentally 
      !measured values. 
      
      use diffract_v
      use mphase_props
      use grain_props, only : wgt
      use grain_state, only : nact,etcs
      use grain_rate, only : aloc2
      use mvoigt
      use sample_props_v, only : ngrain
      use sample_state_v, only : stss
      use sample_rate_v, only : ass2
      
      DIMENSION comp2_s(6,6),comp4_s(3,3,3,3)
      DIMENSION aloc4(3,3,3,3),stss_2(3,3)
      DIMENSION etcs_2(3,3)
           
      ! DEC needs to be based on elastic constants and loc tensor
      ! F_pq = sum_g( n_i*n_j*A_ijkl(g)*S_klpq*Vf(g) ) / sum_g( Vf(g) )      
      ! it works only for elasticity sum_ng(nact(ng)) = 0
      i_el = sum(nact(1:ngrain)) !0 for elasticity
      
c      if(iproc.eq.nproc-1.and.istep.eq.nsteps.and.ilog_res_st.eq.1)then
c          tau(:,:)=10000.0 !no plasticity in process 4          
c      endif
c      
c      if(iproc.eq.nproc-2.and.istep.eq.nsteps.and.ilog_res_st.eq.0)then
c          stcs = stcs*q_res_st
c          ! restart CRSS
c          if (irs.eq.1) then
c              ! read from file for t=293 K
c              open(1,file='MG4_1.SX')
c              do i=1,75
c                  read(1,*)
c              enddo
c              read(1,*)a,b,c,d,e,f1,g1,h1
c              call CRSS_disloc_dens1(0,ng)
c              
c              do ng=1,ngrain
c                 call CRSS_disloc_dens1(1,ng)
c              enddo
c              close(1)
c          endif
c      endif

      if (i_el.eq.0) then
      
      ! find complience in sample frame
      call invten(ass2,comp2_s)
      call voigt(dummy1,dummy2,comp2_s,comp4_s,3)
      
      ! find F (should there be multiplication with constant product factor)
      if (nph.gt.1) then
          write(*,*) 'ERROR: diffracting elastic constants work with one
     # phase'
          pause
      endif
      do n=1,ndiff(1) !go over scattering vectors
      DEC(:,:,n)=0.0    
      do n1=1,ngrset(n) !go over number of diffracting grains for each scattering vector
      ng=igrset(n,n1) !go over diffracting grains for each scattering vector

        ! find localization tensor in tensor form, sample frame
        call voigt(dummy1,dummy2,aloc2(:,:,ng),aloc4,3)
        
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do ip=1,3
                  do iq=1,3
                    DEC(ip,iq,n) = DEC(ip,iq,n) + vs(n,i,1)
     #               *vs(n,j,1)*aloc4(i,j,k,l)*comp4_s(k,l,ip,iq)
     #               *wgt(ng)/wgtset(n)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

      enddo !n1
      enddo !n
      
      
      ! check if DEC multiplied with stss gives lattice strain
      call voigt(stss,stss_2,dummy3,dummy4,1)
      do n=1,ndiff(1) !go over scattering vectors
         als(n)=0.0   
c         do n1=1,ngrset(n) !go over number of diffracting grains for each scattering vector
c           ng=igrset(n,n1) !go over diffracting grains for each scattering vector
         
         call  voigt(etcs(:,ng),etcs_2,dummy3,dummy4,1)
         
         do ip=1,3
           do iq=1,3
             als(n)=als(n) + DEC(ip,iq,n)*stss_2(ip,iq)!*1e6
c             als(n)=als(n) + etcs_2(ip,iq)*vs(n,ip)*vs(n,iq)
c     #        *wgt(ng)/wgtset(n)*1e6
           enddo
         enddo
           
c       enddo !n1
      enddo !n      
      
c      write(178,'(2i10 x 10000e20.5)') iproc,istep,DEC(1,1,1:12)
c     # ,DEC(2,2,1:12)
c      write(179,'(10000e20.5)') als(1:ndiff)
      else 
        DEC(:,:,:)=0.0
      endif
      
      RETURN
      END SUBROUTINE
      
      
      end module diffract
c      
c***********************************************************************
c *** hard_law1 ********************************************************
c***********************************************************************
      module hard_law1
      
      use hard_law1_v
      
      private
      
      public hl_read,hl_ini_tau,hl_update_statv,hl_hd,state_hl
     #    ,set_hl,get_hl,write_statv_hl,read_statv_hl,hl_twin_mfp
     #    ,hl_pts_damp,make_alat,make_alat_FCC
      contains
      
      !hl_read: reads parameters for hardening law
      !hl_ini_tau: calculates initial CRSS and assigns it to tau
      !hl_update_statv: updates state variables of hardening law
      !hl_hd: calculates hardening matrix
      !set_hl: sets statv object
      !get_hl: gets values from statv object
      !write_statv_hl: writes statv object to statv array
      !read_statv_hl: reads statv object from statv array
      
      subroutine set_hl(ng,hl)
      
      use grain_state, only : tau
      
      integer, intent(in) :: ng
      type(state_hl),intent(out) :: hl
      
      hl%tau=tau(:,ng)
      hl%rho_rev=rho_rev(:,ng)
      hl%rho_forw=rho_forw(:,ng)
      hl%rho_tot=rho_tot(:,ng)
      hl%rho_tot_max=rho_tot_max(:,ng)
      hl%rho_act=rho_act(:,ng)
      hl%tau_act=tau_act(:,ng)
      hl%iact_sys=iact_sys(:,ng)
      hl%itemp_act_sys=itemp_act_sys(:,ng)
      hl%rho_deb=rho_deb(ng)
      hl%temp_rho_tot_0=temp_rho_tot_0(:,ng)
      hl%rssmin=rssmin(:,ng)
      
      end subroutine set_hl
      
      subroutine get_hl(ng,hl)
      
      use grain_state, only : tau
      
      integer, intent(in) :: ng
      type(state_hl),intent(in) :: hl
      
      tau(:,ng)=hl%tau
      rho_rev(:,ng)=hl%rho_rev
      rho_forw(:,ng)=hl%rho_forw
      rho_tot(:,ng)=hl%rho_tot
      rho_tot_max(:,ng)=hl%rho_tot_max
      rho_act(:,ng)=hl%rho_act
      tau_act(:,ng)=hl%tau_act
      iact_sys(:,ng)=hl%iact_sys
      itemp_act_sys(:,ng)=hl%itemp_act_sys
      rho_deb(ng)=hl%rho_deb
      temp_rho_tot_0(:,ng)=hl%temp_rho_tot_0
      rssmin(:,ng)=hl%rssmin
      
      end subroutine get_hl
      
      subroutine write_statv_hl(hl,ng,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_hl),intent(in) :: hl
      integer,intent(in) :: ng,NSTATV
      integer,intent(inout) :: ns
      real,intent(out) :: STATEV(NSTATV)
      
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%tau(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%rho_rev(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%rho_forw(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%rho_tot(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%rho_tot_max
     #    (1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%rho_act(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))  
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%tau_act(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%iact_sys(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%itemp_act_sys
     #    (1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns)=hl%rho_deb
      ns=ns+1
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%temp_rho_tot_0
     #    (1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%rssmin(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
            
      end subroutine write_statv_hl
      
      subroutine read_statv_hl(hl,ng,iph,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_hl),intent(out) :: hl
      integer,intent(in) :: ng,NSTATV,iph
      real,intent(in) :: STATEV(NSTATV)
      integer,intent(inout) :: ns
      
      hl%tau(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%rho_rev(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%rho_forw(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%rho_tot(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%rho_tot_max(1:nsys(iph))=STATEV
     #    (ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%rho_act(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)  
      hl%tau_act(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%iact_sys(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%itemp_act_sys(1:nsys(iph))=STATEV
     #    (ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%rho_deb=STATEV(ns)
      ns=ns+1      
      hl%temp_rho_tot_0(1:nsys(iph))=STATEV(ns:ns+nsys(ngrnph(ng))-1)
      ns=ns+nsys(ngrnph(ng))
      hl%rssmin(1:nsys(ngrnph(ng)))=STATEV(ns:ns+nsys(ngrnph(ng))-1)
      ns=ns+nsys(ngrnph(ng))
      
      end subroutine read_statv_hl
      
      SUBROUTINE hl_read(iph)
      
      use const
      use mphase_props
      use flags, only : kCL,iDiag,itwinning,iTwinLaw
      use mvoigt, only : voigt
      use twinning_v, only : PTS_THRES,TwinFrac,TwinCRSS,TWFR_THRES
      use files, only : fileproc
      use bc, only : edot_macro
      
      DIMENSION ccc4(3,3,3,3),aux6(6),aux33(3,3)
      DIMENSION DRAG_a(NMOD,NPHM),DRAG_b(NMOD,NPHM),DRAG_c(NMOD,NPHM)
      DIMENSION cgr(NMOD)

      !Boltzmann's constant
      BOLTZ=1.380622e-23  

      ! Reads hardening law parameters
      READ(1,*)
      READ(1,*) iDiag !,alpha_reg
      READ(1,*) chi_inter(iph),q_rate(iph),grsze
      DO is=1,nslmod(iph)
          READ(1,*)
          READ(1,*) BURG(is,iph),ACTENER(is,iph)
          READ(1,*) aK1(is,iph),DRAG(is,iph)
c          READ(1,*) aK1(is,iph),DRAG_a(is,iph),DRAG_b(is,iph) !marcel
c                ,DRAG_c(is,iph)
          READ(1,*) edot_zero(is,iph)
          READ(1,*) rho_ini_for(is,iph),rho_ini_deb(is,iph)
          READ(1,*) tau0_mode_a(is,iph),tau0_mode_b(is,iph),
     #    tau0_mode_c(is,iph)
          if (itwinning.eq.1) then !MZ_2ph
              READ(1,*) (TLATENT(is,mo,iph),mo=1,ntwmod(iph))
              READ(1,*) (TLATENT1(is,mo,iph),mo=1,ntwmod(iph))
              !read detwinning slip twin interaction
c              READ(1,*) (TLATENT_de(is,mo),mo=1,ntwmod)
c              READ(1,*) (TLATENT1_de(is,mo),mo=1,ntwmod)  
          endif
          READ(1,*) HPK0(is,iph),(HPKCG(is,mo,iph),mo=1,ntwmod(iph))
          READ(1,*) a_deb_a(is,iph),a_deb_b(is,iph),a_deb_c(is,iph)
c      if (iGrShCRSS.eq.1) READ(1,*) cgr(is)     
      ENDDO
      READ(1,*) p_rev(iph),aM_par(iph),rev_par(iph),irevlaw(iph)
      !latent hardening
      READ(1,*) iFCC
      if (iFCC.eq.0) then
          do is1=1,nslmod(iph)
              READ(1,*) (d_mod(is1,is2,iph),is2=1,(nslmod(iph)+1))
          enddo
          do is1=1,nslmod(iph)
              READ(1,*) (g_mod(is1,is2,iph),is2=1,(nslmod(iph)+1))
          enddo
      elseif (iFCC.eq.1) then
          READ(1,*) d0,d1,d2,d3,d4,d5
          READ(1,*) g0,g1,g2,g3,g4,g5
      elseif (iFCC.eq.2) then
          do is1=1,nslmod(iph)
            READ(1,*) (d_mod(is1,is2,iph),is2=1,(nslmod(iph)+1))
          enddo
          do is1=1,nslmod(iph)
            READ(1,*) (g_mod(is1,is2,iph),is2=1,(nslmod(iph)+1))
          enddo
          READ(1,*) d0,d1,d2,d3,d4,d5
          READ(1,*) g0,g1,g2,g3,g4,g5
        endif


      DO it=1,ntwmod(iph)
          READ(1,*)
          READ(1,*) tau_crit_a(it,iph),tau_crit_b(it,iph)
     #      ,tau_crit_c(it,iph)
          READ(1,*) tau_prop_a(it,iph),tau_prop_b(it,iph)
     #      ,tau_prop_c(it,iph)
          READ(1,*) HPK0(it+nslmod(iph),iph),
     #      (HPKCG(it+nslmod(iph),mo,iph),mo=1,ntwmod(iph))
          READ(1,*) BURG(it+nslmod(iph),iph)
          if (iDetwOpt.eq.1) then
              READ(1,*) tau_crit_de_a(it,iph),tau_crit_de_b(it,iph)
     #           ,tau_crit_de_c(it,iph)
              READ(1,*) tau_prop_de_a(it,iph),tau_prop_de_b(it,iph)
     #           ,tau_prop_de_c(it,iph)
          endif
          itt=it+nslmod(iph)
          if (iTwinLaw.eq.1) then
              READ(1,*) TwinFrac(itt,iph),TwinCRSS(itt,iph)   !BJCL twin model parameters
          endif
          if (iTwinLaw.eq.2) READ(1,*) TwinFrac(itt,iph)
     #       ,TWFR_THRES(it,iph)!BJCL twin model parameters
          if (iTwinLaw.eq.2) READ(1,*) PTS_THRES(it,iph)
      ENDDO
      
      !get edot and temp_s from process file
c      call READ_PROC_PAR(fileproc(1),edot_macro,temp_s)
      
      !initialize arrays
      call voigt(aux6,aux33,ccc2(:,:,iph),ccc4,3)
      is=0
      do im=1,nmodes(iph)
          is=is+nsm(im,iph)
          shearmod(im,iph)=0.
          do i=1,3
              do j=1,3
                  do k=1,3
                      do l=1,3
                          shearmod(im,iph)=shearmod(im,iph)+
     #                    mc2(i,j,is,iph)*ccc4(i,j,k,l)
     #                    *mc2(k,l,is,iph)
                      enddo
                  enddo
              enddo
          enddo
      enddo
      
      !initialize iact arrays
      call initialize_iact(iph)
      
      !initialize latent hardening matrix
      if (iFCC.eq.0) then
          call make_alat(iph)
      elseif (iFCC.eq.1) then
          call make_alat_FCC(iph)
      elseif (iFCC.eq.2) then    
          call make_alat(iph)
          call make_alat_FCC(iph)
      endif
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE hl_ini_tau(ng,tau,iopt)
      
      use const
      use files, only : fileproc
      use mphase_props
      use bc, only : edot_macro,temp_s
      use grain_state, only : gamtot
      
      integer,intent (in) :: ng
      real, intent(out) :: tau(NSLS,NGR)

      iph1=ngrnph(ng)
      gamtot(ng)=0.0

      !set up initial state variables
      if (iopt.eq.0) then
          rho_deb(ng)=0.0
          is=0
          do imo=1,nslmod(iph1)
              do isy=1,nsm(imo,iph1)
                  is=is+1
                  !Defines the initial dislocation densities
                  rho_tot(is,ng)=rho_ini_for(imo,iph1)
                  rho_forw(is,ng)=rho_ini_for(imo,iph1)
                  !Finding the array of active dislocation desities.
                  rho_act(iDirSys(is,iph1),ng)=rho_ini_for(imo,iph1)
              enddo
          enddo
          do imo=1,nslmod(iph1)
              rho_deb(ng) = rho_deb(ng) + rho_ini_deb(imo,iph1)
          enddo
      endif
      
      !find tau based on state variables
      is=0
      do imo=1,nslmod(iph1)
          do isy=1,nsm(imo,iph1)
              is=is+1

              !Defines the initial critical stresss
              tau0_mode(imo,iph1)=tau0_mode_a(imo,iph1)+
     #        tau0_mode_b(imo,iph1)*exp(-temp_s
     #            /tau0_mode_c(imo,iph1))
c            tau0_mode(imo)=tau0_mode_a(imo,iph1)/
c     #        (1.0 + exp(-tau0_mode_b(imo,iph1)*
c     #       (temp_s - tau0_mode_c(imo,iph1)))) !marcel

                  !Define dependence of Drag with temperature
c            DRAG(imo,iph1) = DRAG_a(imo,iph1) + DRAG_b(imo,iph1)
c     #                       *exp(DRAG_c(imo,iph1)*temp_s)         !marcel
                  !Hall-Petch term
              HPfac(imo,iph1)=sqrt(burg(imo,iph1))*HPK0(imo,iph1)*
     #          shearmod(imo,iph1)
              if (iGrShCRSS.eq.0) then
                  tauHP(imo,iph1)=HPfac(imo,iph1)/sqrt(grsze)
              endif
              if (iGrShCRSS.eq.1) then
                  tau_GrSh = HPfac(imo,iph1)/sqrt(xmfp(is,ng)) !iGrShCRSS
              endif

              ! Defines the yield point for the slip mode considered
              sqrt_root = sqrt(sum_latent(1,is,ng))
              tau(is,ng)=tau0_mode(imo,iph1)+tauHP(imo,iph1)+chi_inter
     #          (iph1)*burg(imo,iph1)*shearmod(imo,iph1)*sqrt_root
     #          +0.086*shearmod(imo,iph1)*burg(imo,iph1)
     #          *sqrt(rho_deb(ng))*log(1.0/( burg(imo,iph1)
     #          *sqrt(rho_deb(ng)) ) )
              if (iGrShCRSS.eq.1) tau(is,ng) = tau(is,ng) + tau_GrSh
              a_deb(imo,iph1)=a_deb_a(imo,iph1)+a_deb_b(imo,iph1)*
     #          log(1.+temp_s/a_deb_c(imo,iph1))
              !BJCL: If the temperature is changed during the process, should this be re-calculated at each step?
              enddo
          enddo
      ! Consider the twin modes
      do imo=1,ntwmod(iph1)
          do isy=1,nsm(imo+nslmod(iph1),iph1)
              is=is+1
              tau(is,ng)=0.0
              js=0
              do jmo=1,nslmod(iph1)
                  do jsy=1,nsm(jmo,iph1)
                      js=js+1
                      !Equation 3.28 in Beyerlein and Tome 2008:
                      tau(is,ng)=tau(is,ng)+shearmod(imo+nslmod(iph1)
     #                  ,iph1)*burg(imo+nslmod(iph1),iph1)*burg(jmo
     #                  ,iph1)*(tlatent(jmo,imo,iph1)-tlatent1(jmo,imo
     #                  ,iph1)*log(edot_macro))*rho_tot(js,ng)/2. !/2. since there are two times more slip systems
                  enddo
              enddo

              !Hall-Petch term
              HPfac(imo+nslmod(iph1),iph1)=HPK0(imo+nslmod(iph1),iph1)
     #          *sqrt(burg(imo+nslmod(iph1),iph1))
     #          *shearmod(imo+nslmod(iph1),iph1)!martin_BUG
              tauHP(imo+nslmod(iph1),iph1)=HPfac(imo+nslmod(iph1),iph1)
     #          /sqrt(grsze)
              !Equation 3.26 in Beyerlein and Tome 2008 - and the text following 3.26:
              tau(is,ng)=tau(is,ng)+tau_crit_a(imo,iph1)
     #          +tau_crit_b(imo,iph1)*exp(-temp_s/tau_crit_c(imo,iph1))
     #          +tauHP(imo+nslmod(iph1),iph1)!martin_BUG
              !BJCL: If the temperature is changed during the process, should this be re-calculated at each step?
              !detwinning - calculation of tau_crit and tau_prop for detwinning
              tau_crit_de(imo,iph1)=tau_crit_de_a(imo,iph1)+
     #          tau_crit_de_b(imo,iph1)*exp(-temp_s
     #          /tau_crit_de_c(imo,iph1))
              tau_prop_de(imo,iph1)=tau_prop_de_a(imo,iph1)+
     #          tau_prop_de_b(imo,iph1)*exp(-temp_s
     #          /tau_prop_de_c(imo,iph1))
          enddo
      enddo
      !Finding the array of active CRSS.
      if (iopt.eq.0) then
          do is=1,nsys(iph1)
              tau_act(iDirSys(is,iph1),ng)=tau(is,ng)
          enddo
      endif
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE hl_update_statv(ng,step)
      
      use const
      use mphase_props
      use grain_state, only : iact,nact
      use grain_rate, only : gamd
      
      integer, intent(in) :: ng
      real, intent(in) :: step
      
      iph1=ngrnph(ng)
      
      if(sum(gamd(:,ng)).ne.0.0) then
      do ns=1,nslsys(ngrnph(ng))
          gam_acc(ns,ng)=gam_acc(ns,ng)+gamd(ns,ng)
      enddo      
      
c     Set up array of active systems.
      do is=1,nact(ng)
          n1=iact(is,ng)
          if (iTwinSys(n1).EQ.0) then
              if (mod(n1,2).EQ.1) then
                  iact_sys(iDirSys(n1,ngrnph(ng)),ng)=1
              else
                  iact_sys(iDirSys(n1,ngrnph(ng)),ng)=2
              endif
          else
              iact_sys(iDirSys(n1,ngrnph(ng)),ng)=1
          endif
      enddo  
      
      !update drho_dgamma derivatives
      call Drho_Dgamma(ng)
      
      !comensate for variation of active slip systems during iterations
      temp_rho_tot_0(1:nslsys(ngrnph(ng)),ng)=
     #  rho_tot_0(1:nslsys(ngrnph(ng)),ng) !store the converged values

      !find maximum value of tot disloc dens up to this point
      do is=1,nslsys(ngrnph(ng))
          if(rho_tot(is,ng).gt.rho_tot_max(is,ng)) then
              rho_tot_max(is,ng)=rho_tot(is,ng)
          endif
      enddo
      !update rev disloc
      !first update active systems (gamd>0)
      do is=1,nslsys(ngrnph(ng))
          rho_rev(is,ng)=rho_rev(is,ng)+Drho_rev_Dgamma(is)
     #    *gamd(is,ng)
      enddo
      !update opposite to active systems
      do is1=1,nact(ng)
          n1=iact(is1,ng)
          iop=iOpposite(n1)
          if (iTwinSys(n1).eq.0) then !TW
              if (rho_tot_0(iop,ng).gt.0.0) then
                  rho_rev(iop,ng)=rho_rev(iop,ng)+Drho_rev_Dgamma
     #            (iop)*gamd(n1,ng)
              endif
          endif
      enddo
      !update forward dislocations
      do is=1,nslsys(ngrnph(ng))
          iop=iOpposite(is)
          rho_forw(is,ng)=rho_forw(is,ng)+Drho_forw_Dgamma(is)
     #    *gamd(is,ng)
          rho_forw(iop,ng)=rho_forw(iop,ng)+Drho_forw_Dgamma(is)
     #    *gamd(is,ng)
      enddo
      !make sure that rev and forward disloc is greater than 0
      do is=1,nslsys(ngrnph(ng))
          if(rho_rev(is,ng).lt.0.0) then
              rho_rev(is,ng)=0.0
          endif
          if(rho_forw(is,ng).lt.0.0) then
              rho_forw(is,ng)=0.0
          endif
      enddo
      
      !update rho_for,rho_rev1,rho_rev2 by adding corresponding derivative times dgamma
      !rho_tot is the sum of rho_for,rho_rev1,rho_rev2
      
      do imo=1,nslmod(ngrnph(ng))
          do i=1,nsm(imo,ngrnph(ng))
              is=mode_slip(imo,i,ngrnph(ng))
              !MZ_revd update total disloc densitey as sum of rev and forward (2nd way)
              rho_tot(is,ng)=rho_tot(is,ng)+Drho_tot_Dgamma(is) !MZ_revd
     #          *gamd(is,ng)
              rho_tot(iOpposite(is),ng)=rho_tot(iOpposite(is),ng) !MZ_revd
     #          +Drho_tot_Dgamma(is)*gamd(is,ng)
              !Forces rho_for to never get smaller than its initial value --> necessary?
              if(rho_tot(is,ng).lt.rho_ini_for(imo,iph1))
     #          rho_tot(is,ng) =  rho_ini_for(imo,iph1)
              rho_deb(ng)=rho_deb(ng)+Drho_deb_Dgamma(is)
     #          *gamd(is,ng)             
              !Forces rho_deb to never get smaller than its initial value --> necessary?
              if(rho_deb(ng).lt.rho_ini_deb(imo,iph1))
     #          rho_deb(ng)  = rho_ini_deb(imo,iph1)
          enddo
      enddo
      
      !write set of active slip systems to memory
      do is=1,nslsys(ngrnph(ng))/2+ntwsys(ngrnph(ng))
          itemp_act_sys(is,ng)=iact_sys(is,ng)
      enddo  
      
      !update rho_rev during loading
      if (rev_par(ngrnph(ng)).gt.0.0) then
          call hl_dstatv_drss(ng,step)
      endif
	  
      endif
      RETURN
      END SUBROUTINE
      
      subroutine hl_dstatv_drss(ng,step)
      
      use const
      use mvoigt
      use flags, only : iLatBS,iBackStress
      use mphase_props
      use grain_state_v, only : iact,nact,stcs,tau
      use grain_rate_v, only : gamd,strcs
      use grain_props_v, only : nmcs,mcs
      use back_stress_v, only : tau_bcst
      
      integer, intent(in) :: ng
      real, intent(in) :: step
      real :: drho_rev(NSLS),drho_tot(NSLS),dtau(NSLS)
      
      !restart rssmin for active slip systems
      do ns1=1,nact(ng)
          n1=iact(ns1,ng)
          rssmin(n1,ng)=tau(n1,ng)
      enddo
      
      iph=ngrnph(ng)
      drho_tot(1:nslsys(ngrnph(ng)))=0.0
      do is=1,nslsys(ngrnph(ng))
          imo=iSysMode(is,ngrnph(ng))
          drho_rev(is)=0.0
          !find rss
          rss_old=rss
          rss=0.0
          drss=0.0
          do i=1,6 
              rss=rss+(mcs(i,is,ng)+nmcs(i,is,ng))
     #               *stcs(i,ng)*profac(i)  !inonSch
              drss=drss+(mcs(i,is,ng)+nmcs(i,is,ng))
     #               *strcs(i,ng)*profac(i)  
          enddo
          if(iBackStress.eq.1) rss=rss-tau_bcst(is,ng
     #       ,ilatBS)
          if (rss.lt.rssmin(is,ng)) rssmin(is,ng)=rss
          !rho_rev update
          if (rho_rev(is,ng).gt.0.0.and.rss.le.rssmin(is,ng).and.
     #        rss.lt.0.0.and..not.
     #        ANY(mask=iact(1:nact(ng),ng).eq.is).and.drss.lt.0.0)then
              drho_rev(is)=rev_par(iph)*(rho_rev(is,ng)
     #            +rho_rev(iOpposite(is),ng))/tau0_mode(imo,iph)*drss
              if (rho_rev(is,ng)+drho_rev(is).gt.0.0) then
                  rho_rev(is,ng)=rho_rev(is,ng)+drho_rev(is)
              endif
              !erase opposite sign disloc as well
              if (rho_rev(iOpposite(is),ng)+drho_rev(is).gt.0.0) then
                  rho_rev(iOpposite(is),ng)=rho_rev(iOpposite(is),ng)
     #                +drho_rev(is)
              endif
              drho_tot(is)=drho_tot(is)+drho_rev(is)
              drho_tot(iOpposite(is))=drho_tot(iOpposite(is))
     #            +drho_rev(is)
          endif
      enddo
      
      do is=1,nslsys(ngrnph(ng))
          rho_tot(is,ng)=rho_tot(is,ng)+drho_tot(is)
          rho_tot_max(is,ng)=rho_tot_max(is,ng)+drho_tot(is)
      enddo
      
      do ns1=1,nslsys(ngrnph(ng))
          dtau(ns1)=0.0
          do ns2=1,nslsys(ngrnph(ng))
             dtau(ns1)=dtau(ns1)+dtau_drho(ns1,ns2,ng)*drho_tot(ns2)/2.0
          enddo
      enddo
      do ns1=1,nsys(ngrnph(ng))
	    tau(ns1,ng)=tau(ns1,ng)+dtau(ns1)*step
      enddo
      
c      if (ng.eq.1) write(195,'(1000e25.5)') rho_tot(1:48,1),tau(1:48,1)
c     #    ,rssmin(1:48,1),gamd(1:48,1),dtau(1:48)
      
      return
      end subroutine
      
      SUBROUTINE Drho_Dgamma(ng)
      
      use const     
      use mphase_props
      use grain_state, only : nact,iact
      use bc,only : edot_macro,temp_s
      
      !in
      integer, intent(in) :: ng

      iph1=(ngrnph(ng))
      !compensate for variation of active slip systems during iterations
      rho_tot_0(1:nslsys(ngrnph(ng)),ng)=
     #  temp_rho_tot_0(1:nslsys(ngrnph(ng)),ng) !restore to previously defined value
      ! identification of system that is reversed and defining rho_tot_0 for it
      do is1=1,nact(ng)
          n1=iact(is1,ng)
          id=iDirSys(n1,iph1)
          if (iTwinSys(n1).eq.0) then
              if (itemp_act_sys(id,ng).NE.0) then !system was active at some point
                  iDifference=itemp_act_sys(id,ng)-iact_sys(id,ng)
              endif
              if (iDifference.ne.0) then
                  rho_tot_0(iOpposite(n1),ng)=rho_tot(n1,ng)
              endif
          endif
      enddo

      do is=1,nslsys(ngrnph(ng))
          Drho_deb_Dgamma(is)=0.
      enddo
      is=0
      do imo=1,nslmod(ngrnph(ng))
          !Equation 3.12 in Beyerlein and Tome 2008
          !The drag is given in MPa (N/mm^2) in the SX file and it needs to enter in the equation in Pa (N/m^2)
          aK2(imo,iph1)=( aK1(imo,iph1)*chi_inter(iph1)*burg(imo,iph1)
     #      /ActEner(imo,iph1) ) *( 1.-( Boltz*temp_s/( drag(imo,iph1)
     #      *1.e6*burg(imo,iph1)**3 ) )*log(edot_macro
     #      /edot_zero(imo,iph1)) )
c      if(ng.eq.1) write(*,*) 'k2/k1 = ',aK2/aK1(imo)

          do i=1,nsm(imo,ngrnph(ng))
              is=is+1
              !Equation 3.14 and 3.15 in Beyerlein and Tome 2008
              Drho_deb_Dgamma(is)=aK2(imo,iph1)*rho_tot(is,ng)*burg
     #        (imo,iph1)*sqrt(rho_deb(ng))*a_deb(imo,iph1)*q_rate(iph1)
              !The two lines below limits the change in rho to be
              !non-negative... Is that required?
              if(Drho_deb_Dgamma(is).LT.0.) Drho_deb_Dgamma(is)=0.
          enddo
      enddo

      !evolving reversable dislocations - defining derivative
      !active and opposite to active have deined evolution, rest is 0
      do is=1,nslsys(ngrnph(ng))
          Drho_rev_Dgamma(is)=0.0
          Drho_forw_Dgamma(is)=0.0
          Drho_tot_Dgamma(is)=0.0
      enddo
      !find active and opposite to currently active system and define derivative for it
      do is1=1,nact(ng)
          n1=iact(is1,ng)
          if (iTwinSys(n1).ne.1) then !TW
              iop=iOpposite(n1)
              imo=iSysMode(n1,ngrnph(ng))

              if (rho_tot_0(iop,ng).gt.0.0) then
                  Drho_rev_Dgamma(n1) = p_rev(iph1)*aK1(imo,iph1)*
     #            sqrt(sum_latent(2,n1,ng) !rev_iso (rho_forw(n1,ng)+rho_rev(n1,ng))
     #            )*f_rev(n1,ng)-aK2(imo,iph1)*rho_rev(n1,ng)
              else
                  Drho_rev_Dgamma(n1) = p_rev(iph1)*aK1(imo,iph1)*
     #            sqrt(sum_latent(2,n1,ng))-aK2(imo,iph1) !rho_forw(n1,ng)+rho_rev(n1,ng)
     #            *rho_rev(n1,ng) !rev_iso
              endif

              if (rho_tot_0(iop,ng).gt.0.0) then
                  Drho_rev_Dgamma(iop) = -aK1(imo,iph1)*sqrt(
     #            sum_latent(2,iop,ng) ) !rho_forw(n1,ng)+rho_rev(n1,ng)
     #            *(rho_rev(iop,ng)/rho_tot_0(iop,ng))**aM_par(iph1)
              endif
          endif
      enddo

      !evolve forward dislocations define derivatives
      do is1=1,nact(ng)
          is=iact(is1,ng)
          iop=iOpposite(is)
          imo=iSysMode(is,ngrnph(ng))
          if (iTwinSys(is).ne.1) then !TW
              Drho_forw_Dgamma(is) = (1.0-p_rev(iph1))*aK1(imo,iph1)
     #          *sqrt(sum_latent(2,is,ng)) !rho_forw(is,ng)+rho_rev(is,ng)
     #          -aK2(imo,iph1)*rho_forw(is,ng)
              Drho_forw_Dgamma(iop) = Drho_forw_Dgamma(is)
c      if(Drho_forw_Dgamma(is).LT.0.) Drho_forw_Dgamma(is)=0.
          endif
      enddo

      !define derivative of total disloc dens
      do is1=1,nact(ng)
          is=iact(is1,ng)
          if (iTwinSys(is).ne.1) then !TW
              Drho_tot_Dgamma(is)=Drho_forw_Dgamma(is)
     #          +Drho_rev_Dgamma(is)+Drho_rev_Dgamma(iOpposite(is))
              Drho_tot_Dgamma(iOpposite(is))=Drho_tot_Dgamma(is)
          endif
      enddo
      !define derivative used for hd matrix (>0) and disloc dens needs
      !to be higher than the previous one at point of reversal
      do is=1,nslsys(ngrnph(ng))
          if (Drho_tot_Dgamma(is).gt.0.0.and.rho_tot(is,ng).ge.
     #    rho_tot_max(is,ng)) then
              Drho_tot_Dgamma_hd(is)=Drho_tot_Dgamma(is)
          else
              Drho_tot_Dgamma_hd(is)=0.0
          endif
      enddo      
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE hl_hd(ng,hd)
      
      use const     
      use flags, only : iDtwMfp
      use mphase_props
      use grain_props_v, only : wgt
      use twinning_v
      use grain_state, only : nact,iact,gamtot,tau
      use grain_rate, only : gamd
      use bc,only : edot_macro
      use sample_props_v, only : ngParent
      
      !in
      integer, intent(in) :: ng
      
      !out
      real, intent(out) :: hd(NSLS,NSLS)
      
      iph1=(ngrnph(ng))
      
      !update drho_dgamma derivatives
      call Drho_Dgamma(ng)

      hd(:,:)=0.
      do n1=1,nsys(ngrnph(ng)) !goes over all sys, not only active
          ns1=iact(n1,ng)
          !Identify to which mode the slip sys is related
          imode1=iSysMode(n1,ngrnph(ng))
          imode2=0
          do ns2=1,nact(ng)
              n2=iact(ns2,ng)
              imode2=iSysMode(n2,ngrnph(ng))
              !slip-slip interaction
              if (iTwinsys(n1).NE.1.and.iTwinsys(n2).NE.1) then

                  !DEBRIS DISLOCATION (Derivative of 3.20)
                  hd(n1,n2)= hd(n1,n2)-0.086*burg(imode1,iph1)*shearmod
     #              (imode1,ngrnph(ng))*(  log( burg(imode1,iph1)*
     #              sqrt(rho_deb(ng)) )+1.  )*Drho_deb_Dgamma(n2)/
     #              (2.*sqrt(rho_deb(ng)))
                  
                  !FOREST DISLOCATION (Derivative of 3.19)
                  hd(n1,n2)= hd(n1,n2)+Drho_tot_Dgamma_hd(n2)
     #              *dtau_drho(n1,n2,ng)

              !twin-slip interaction
              else if (iTwinsys(n1).EQ.1.and.iTwinsys(n2).NE.1) then
              !find if tau_tw is large enough for systems that have a child
                  iflag=0
                  if (iChildGrain(n1,ng).ne.0.and.
     #            tau_tw(n1).gt.tau(n1,ng)) iflag=1
                  !EFFECT OF DISLOCATIONS ON SLIP SYSTEMS
                  hd(n1,n2)=shearmod(imode1,ngrnph(ng))*burg
     #            (imode1,iph1)*burg(imode2,iph1)*
     #            (tlatent(imode2,imode1-nslmod(ngrnph(ng)),iph1)-
     #            tlatent1(imode2,imode1-nslmod(ngrnph(ng)),iph1)
     #            *log(edot_macro))*Drho_tot_Dgamma_hd(n2)

                  !slip-twin interaction
              else if (iTwinsys(n1).NE.1.and.iTwinsys(n2).EQ.1) then
                  nt2=n2-nslsys(ngrnph(ng))
                  imodet2=imode2-nslmod(ngrnph(ng))

                  if (ng.le.ngParent) then !----------------CG for parents

                      if (IPTSGRC(ng).eq.n2) then
                      !calculate tauHP due to mean free path
                          tauHP_CG=HP_CG(ng,iph1,n1)
                          if (tauHP_CG.ge.tauHP(imode1,iph1)) then
                              d_c=grsze*TWTHRES(1,imode2,iph1)
                              HP1fac=HPKCG(imode1,imodet2,iph1)*sqrt(
     #                          burg(imode1,iph1))*shearmod(imode1,iph1)
                              hd(n1,n2)=hd(n1,n2)+
     #                          (  1.0/((TWFR_THRES(imodet2,iph1)
     #                          -PTS_THRES(imodet2,iph1)))*HP1fac/((1.0
     #                          -TWFRSY(nt2,ng))*d_c*sysmfp(n1,nt2))**
     #                          0.5 + 0.5*(TWFRSY(nt2,ng)-PTS_THRES
     #                          (imodet2,iph1))/(TWFR_THRES(imodet2,iph1
     #                          )-PTS_THRES(imodet2,iph1))*HP1fac/((1.0
     #                          -TWFRSY(nt2,ng))*d_c*sysmfp(n1,nt2))**
     #                          1.5*d_c*sysmfp(n1,nt2)   )*wgt(ng)
     #                          /wgtcg_all(ng)/stw(imode2,iph1)

                          endif
                      endif

                          !for twins during detwinning
                  elseif (iDtwMfp.eq.1) then
                      ngp=iParentGrain(ng)
                      ntp=iParentSystem(ng)-nslsys(ngrnph(ng))
                      ntp_abs=iParentSystem(ng)
                      imtp=iParentMode(ng)
                  
                      iDtwSys=iParentSystem(ng)
                      if (n2.eq.iDtwSys) then
                          if (IPTSGRC(ngp).gt.0) then
                              d_c=grsze*TWTHRES(1,imode2,iph1)
                              HP1fac=HPKCG(imode1,imodet2,iph1)
     #                          *sqrt(burg(imode1,iph1))
     #                          *shearmod(imode1,ngrnph(ng))
                              !derivaive
                              !COMMENT TO EXCLUDE DETWINNING EFFECT ON CRSS IN TWINS
                              hd(n1,n2)=hd(n1,n2)+
     #                          (  0.5*HP1fac/(TWFRSY(ntp,ngp)*d_c
     #                          *sysmfp(n1,ntp))**1.5*d_c*sysmfp(n1,ntp)
     #                            )*wgt(ng)/wgtcg_all(ng)/stw(imtp,iph1)
                              if (n1.eq.n2) hd(n1,n2)=0.0
                          endif
                      endif
                  endif

                  if (hd(n1,n2).lt.0.0) hd(n1,n2)=0.0
              !twin-twin interaction
              else if (iTwinsys(n1).EQ.1.and.iTwinsys(n2).EQ.1) then
                  nt1=n1-nslsys(ngrnph(ng))
                  imodet1=imode1-nslmod(ngrnph(ng))
                  nt2=n2-nslsys(ngrnph(ng))
                  imodet2=imode2-nslmod(ngrnph(ng))
                  !CG for parents
                  if (ng.le.ngParent) then
                      if (IPTSGRC(ng).eq.n2) then
                      !calculate tauHP due to mean free path
                          tauHP_CG=HP_CG(ng,iph1,n1)
                          if (tauHP_CG.ge.tauHP(imode1,iph1)) then
                              d_c=grsze*TWTHRES(1,imode2,iph1)
                              HP1fac=HPKCG(imode1,imodet2,iph1)
     #                          *sqrt(burg(imode1,iph1))
     #                          *shearmod(imode1,ngrnph(ng))
                              hd(n1,n2)=hd(n1,n2)+
     #                          (  1.0/((TWFR_THRES(imodet2,iph1)
     #                          -PTS_THRES(imodet2,iph1)))*HP1fac/((1.0
     #                          -TWFRSY(nt2,ng))*d_c*sysmfp(n1,nt2))
     #                          **0.5 + 0.5*(TWFRSY(nt2,ng)
     #                          -PTS_THRES(imodet2,iph1))/
     #                          (TWFR_THRES(imodet2,iph1)-
     #                          PTS_THRES(imodet2,iph1))*HP1fac
     #                          /((1.0-TWFRSY(nt2,ng))*d_c*
     #                          sysmfp(n1,nt2))**1.5*d_c*sysmfp(n1,nt2)
     #                          )*wgt(ng)/wgtcg_all(ng)/stw(imode2,iph1)
                          endif
                      endif
                  elseif (iDtwMfp.eq.1) then !for twins during detwinning
              
                      ngp=iParentGrain(ng)
                      ntp=iParentSystem(ng)-nslsys(ngrnph(ng))
                      ntp_abs=iParentSystem(ng)
                      imtp=iParentMode(ng)
              
                      iDtwSys=iParentSystem(ng)
                      if (n2.eq.iDtwSys) then
                          if (IPTSGRC(ngp).gt.0) then
                              d_c=grsze*TWTHRES(1,imode2,iph1)
                              HP1fac=HPKCG(imode1,imodet2,iph1)
     #                          *sqrt(burg(imode1,iph1))*shearmod
     #                          (imode1,ngrnph(ng))
                              !derivaive
                              !COMMENT TO EXCLUDE DETWINNING EFFECT ON CRSS IN TWINS
                              hd(n1,n2)=hd(n1,n2)+
     #                          (  0.5*HP1fac/(TWFRSY(ntp,ngp)*d_c*
     #                          sysmfp(n1,ntp))**1.5*d_c*sysmfp(n1,ntp)
     #                          )*wgt(ng)/wgtcg_all(ng)/stw(imtp,iph1)
                              if (n1.eq.n2) hd(n1,n2)=0.0
                          endif
                      endif
                  endif
                  if (hd(n1,n2).lt.0.0) hd(n1,n2)=0.0
              endif
          enddo      ! end of do ns2
      enddo        ! end of do ns1

      !Ensure all components of the hardening matrix are larger than 1.0
      do ns1=1,nact(ng)
        n1=iact(ns1,ng)
        do ns2=1,nact(ng)
          n2=iact(ns2,ng)
            if(hd(n1,n2).LT.1.0) then 
                hd(n1,n2)=1.0
                if (iTwinSys(n1).NE.1) then
                    hd(iOpposite(n1),n2)=1.0
                endif
            endif
c          if(ng.eq.1) write(*,'(E8.2,'' '')',advance='no') hd(n1,n2)              
        enddo
      enddo     
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE initialize_iact(iph)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
c      SUBROUTINE ACTIVE_SYS(ng,iopt)      
c      initiates arrays used for determination of currently active sys, disloc
c         dens, CRSS.
c         iDirSys(is,nph) - for given index of duble system (is) returns value of direct system
c         iOpposite(is) - for given value of double slip system (is) returns value of 
c                         opposite system (system with positive/negative shear)        
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
      use mphase_props

c     Set up direct system from double.          
      is=0
      do is1=1,nslsys(iph)/2
          do is2=1,2
              is=is+1
              iDirSys(is,iph)=is1
          enddo    
      enddo
      do it=1,ntwsys(iph)
          is=is+1
          iDirSys(is,iph)=nslsys(iph)/2+it
      enddo    
c     Set up opposite system.
      do is=1,nslsys(iph)
          if (mod(is,2).EQ.1) then
              iOpposite(is)=is+1
          else 
              iOpposite(is)=is-1
          endif    
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE make_alat(iph)
      
      use mphase_props, only : iSysMode,nslsys
      
      do is1=1,nslsys(iph)
          im1=iSysMode(is1,iph)
          do is2=1,nslsys(iph)
              im2=iSysMode(is2,iph)
              if (is1.eq.is2.or.is1.eq.iOpposite(is2)) then
                  alatentDD(is1,is2,iph)=g_mod(im1,1,iph)
                  alatent(is1,is2,iph)=d_mod(im1,1,iph)
              else
                  alatentDD(is1,is2,iph)=g_mod(im1,im2+1,iph)
                  alatent(is1,is2,iph)=d_mod(im1,im2+1,iph)
              endif
          enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE make_alat_FCC(iph)
      
      real :: aLat(12,12)
      
      !CRSS latent
        aLat(1,1)= d0      !A2
        aLat(1,2)= d1      !A3
        aLat(1,3)= d1      !A6
        aLat(1,4)= d3      !B2
        aLat(1,5)= d4      !B4
        aLat(1,6)= d4      !B5
        aLat(1,7)= d2      !C1
        aLat(1,8)= d4      !C3
        aLat(1,9)= d5      !C5
        aLat(1,10)=d2      !D1
        aLat(1,11)=d5      !D4
        aLat(1,12)=d4      !D6
        
        !A3!A3
        aLat(2,2)= d0      !A3
        aLat(2,3)= d1      !A6
        aLat(2,4)= d4      !B2
        aLat(2,5)= d2      !B4
        aLat(2,6)= d5      !B5
        aLat(2,7)= d4      !C1
        aLat(2,8)= d3      !C3
        aLat(2,9)= d4      !C5
        aLat(2,10)=d5      !D1
        aLat(2,11)=d2      !D4
        aLat(2,12)=d4      !D6        
        
        !A6!A6
        aLat(3,3)= d0      !A6
        aLat(3,4)= d5      !B2
        aLat(3,5)= d5      !B4
        aLat(3,6)= d2      !B5
        aLat(3,7)= d5      !C1
        aLat(3,8)= d4      !C3
        aLat(3,9)= d2      !C5
        aLat(3,10)=d4      !D1
        aLat(3,11)=d4      !D4
        aLat(3,12)=d3      !D6          
        
        !B2!B2
        aLat(4,4)= d0      !B2
        aLat(4,5)= d1      !B4
        aLat(4,6)= d1      !B5
        aLat(4,7)= d2      !C1
        aLat(4,8)= d5      !C3
        aLat(4,9)= d4      !C5
        aLat(4,10)=d2      !D1
        aLat(4,11)=d4      !D4
        aLat(4,12)=d5      !D6  
        
        !B4!B4
        aLat(5,5)= d0      !B4
        aLat(5,6)= d1      !B5
        aLat(5,7)= d5      !C1
        aLat(5,8)= d2      !C3
        aLat(5,9)= d4      !C5
        aLat(5,10)=d4      !D1
        aLat(5,11)=d3      !D4
        aLat(5,12)=d4      !D6  
        
        !B5!B5
        aLat(6,6)= d0      !B5
        aLat(6,7)= d4      !C1
        aLat(6,8)= d4      !C3
        aLat(6,9)= d3      !C5
        aLat(6,10)=d5      !D1
        aLat(6,11)=d4      !D4
        aLat(6,12)=d2      !D6  
        
        !C1!C1
        aLat(7,7)= d0      !C1
        aLat(7,8)= d1      !C3
        aLat(7,9)= d1      !C5
        aLat(7,10)=d3      !D1
        aLat(7,11)=d4      !D4
        aLat(7,12)=d4      !D6  
        
        !C3!C3
        aLat(8,8)= d0      !C3
        aLat(8,9)= d1      !C5
        aLat(8,10)=d4      !D1
        aLat(8,11)=d2      !D4
        aLat(8,12)=d5      !D6  
        
        !C5!C5
        aLat(9,9)= d0      !C5
        aLat(9,10)=d4      !D1
        aLat(9,11)=d5      !D4
        aLat(9,12)=d2      !D6  
        
        !D1!D1
        aLat(10,10)=d0     !D1
        aLat(10,11)=d1     !D4
        aLat(10,12)=d1     !D6  
                           
        !D1!D1             
        aLat(11,11)=d0     !D4
        aLat(11,12)=d1     !D6  
                           
        !D1!D1             
        aLat(12,12)=d0     !D6          
        
        do i=1,12
        do j=1,12
          if(j.lt.i) aLat(i,j) = aLat(j,i)
        enddo
        enddo   
        
        !make alatent
        do i=1,24
        do j=1,24
          alatent(i,j,iph)=alat(iDirSys(i,iph),iDirSys(j,iph))
        enddo
        enddo
        
      !CRSS latent
        aLat(1,1)= g0      !A2
        aLat(1,2)= g1      !A3
        aLat(1,3)= g1      !A6
        aLat(1,4)= g3      !B2
        aLat(1,5)= g4      !B4
        aLat(1,6)= g4      !B5
        aLat(1,7)= g2      !C1
        aLat(1,8)= g4      !C3
        aLat(1,9)= g5      !C5
        aLat(1,10)=g2      !D1
        aLat(1,11)=g5      !D4
        aLat(1,12)=g4      !D6
        
        !A3!A3
        aLat(2,2)= g0      !A3
        aLat(2,3)= g1      !A6
        aLat(2,4)= g4      !B2
        aLat(2,5)= g2      !B4
        aLat(2,6)= g5      !B5
        aLat(2,7)= g4      !C1
        aLat(2,8)= g3      !C3
        aLat(2,9)= g4      !C5
        aLat(2,10)=g5      !D1
        aLat(2,11)=g2      !D4
        aLat(2,12)=g4      !D6        
        
        !A6!A6
        aLat(3,3)= g0      !A6
        aLat(3,4)= g5      !B2
        aLat(3,5)= g5      !B4
        aLat(3,6)= g2      !B5
        aLat(3,7)= g5      !C1
        aLat(3,8)= g4      !C3
        aLat(3,9)= g2      !C5
        aLat(3,10)=g4      !D1
        aLat(3,11)=g4      !D4
        aLat(3,12)=g3      !D6          
        
        !B2!B2
        aLat(4,4)= g0      !B2
        aLat(4,5)= g1      !B4
        aLat(4,6)= g1      !B5
        aLat(4,7)= g2      !C1
        aLat(4,8)= g5      !C3
        aLat(4,9)= g4      !C5
        aLat(4,10)=g2      !D1
        aLat(4,11)=g4      !D4
        aLat(4,12)=g5      !D6  
        
        !B4!B4
        aLat(5,5)= g0      !B4
        aLat(5,6)= g1      !B5
        aLat(5,7)= g5      !C1
        aLat(5,8)= g2      !C3
        aLat(5,9)= g4      !C5
        aLat(5,10)=g4      !D1
        aLat(5,11)=g3      !D4
        aLat(5,12)=g4      !D6  
        
        !B5!B5
        aLat(6,6)= g0      !B5
        aLat(6,7)= g4      !C1
        aLat(6,8)= g4      !C3
        aLat(6,9)= g3      !C5
        aLat(6,10)=g5      !D1
        aLat(6,11)=g4      !D4
        aLat(6,12)=g2      !D6  
        
        !C1!C1
        aLat(7,7)= g0      !C1
        aLat(7,8)= g1      !C3
        aLat(7,9)= g1      !C5
        aLat(7,10)=g3      !D1
        aLat(7,11)=g4      !D4
        aLat(7,12)=g4      !D6  
        
        !C3!C3
        aLat(8,8)= g0      !C3
        aLat(8,9)= g1      !C5
        aLat(8,10)=g4      !D1
        aLat(8,11)=g2      !D4
        aLat(8,12)=g5      !D6  
        
        !C5!C5
        aLat(9,9)= g0      !C5
        aLat(9,10)=g4      !D1
        aLat(9,11)=g5      !D4
        aLat(9,12)=g2      !D6  
        
        !D1!D1
        aLat(10,10)=g0     !D1
        aLat(10,11)=g1     !D4
        aLat(10,12)=g1     !D6  
                           
        !D1!D1             
        aLat(11,11)=g0     !D4
        aLat(11,12)=g1     !D6  
                           
        !D1!D1             
        aLat(12,12)=g0     !D6          
        
        do i=1,12
        do j=1,12
          if(j.lt.i) aLat(i,j) = aLat(j,i)
        enddo
        enddo   
        
        !make alatent
        do i=1,24
        do j=1,24
          alatentDD(i,j,iph)=alat(iDirSys(i,iph),iDirSys(j,iph))
        enddo
        enddo
        
      RETURN
      END SUBROUTINE
      
      SUBROUTINE READ_PROC_PAR(fileproc,edot_macro,temp_s)
      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c      SUBROUTINE READ_PROC_PAR
c     Main purpes is to read edot_macro from process file. It's only called first 
c     time Subroutine CRSS_disloc_dens1 is called, because at this point edot_macro
c     hasn't been read from edot_macro and it's needed in CRSS_disloc_dens1.
c     Later on edot_macro is read by SUBROUTINE data_process.      
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      
      use const
      
      character(len=150) fileproc(NPROCX)
      
      CHARACTER*78 prosa
      
    1 FORMAT(a)

      OPEN(unit=301,file=fileproc(1))
      do i=1,29
          READ(301,1) PROSA
      enddo
      READ(301,*) temp_s
      do i=31,40
          READ(301,1) PROSA
      enddo
      READ(301,*) edot_macro
      CLOSE(301)
      
      RETURN
      END SUBROUTINE
      
      FUNCTION sum_latent(iopt,ns1,ng)
      
      use mphase_props
      
      if (iopt.eq.1) then !for strength use alatent
      ! ns1 - active system (sum all)
        sum_latent = 0.0
        imode1 = iSysMode(ns1,ngrnph(ng))
        do ns2=1,nslsys(ngrnph(ng))
            imode2 = iSysMode(ns2,ngrnph(ng))
            sum_latent = sum_latent + alatent(ns1,ns2,ngrnph(ng))
     #                   *rho_tot(ns2,ng) !sum all rho_tot
        enddo
c        sum_latent = sum_latent - 2.0*alatent(ns1,ns2,ngrnph(ng))
c       #             *rho_tot(ns1,ng) ! erase the ones from active sys
        sum_latent = sum_latent/2.0 ! devide by 2.0 to exclude repetition
      endif
      
      if (iopt.eq.2) then !for evolution use alatentDD
        ! ns1 - active system (sum all)
        sum_latent = 0.0
        imode1 = iSysMode(ns1,ngrnph(ng))
        do ns2=1,nslsys(ngrnph(ng))
            imode2 = iSysMode(ns2,ngrnph(ng))
            if(ns2.ne.ns1.and.ns2.ne.iOpposite(ns1)
     #          .or.irevlaw(ngrnph(ng)).eq.1)then
                sum_latent = sum_latent + alatentDD(ns1,ns2,ngrnph(ng))
     #                   *rho_tot(ns2,ng) !sum all rho_tot
            elseif (ns2.eq.ns1) then
                sum_latent = sum_latent + alatentDD(ns1,ns2,ngrnph(ng))
     #                   *(rho_forw(ns1,ng)+rho_rev(ns1,ng))
            elseif (ns2.eq.iOpposite(ns1)) then
                sum_latent = sum_latent + alatentDD(ns1,ns2,ngrnph(ng))
     #                   *(rho_forw(ns1,ng)+rho_rev(ns1,ng))
            endif
            
        enddo
        sum_latent = sum_latent/2.0 ! devide by 2.0 to exclude repetition
        !subtract the opposite to active reversible dislocations
        
      endif      
      END FUNCTION
      
      FUNCTION f_rev(n1,ng)
        
      !n1 - active system
      !ng - grain
      !rev_iso - retard the evolution of dislocatins if there are already
      !present dislocations in the opposite direction       
      if(rho_rev(n1,ng).gt.rho_rev(iOpposite(n1),ng)) then
          f_rev = 1.0
      else
          f_rev = 1.0
      endif
      
      END FUNCTION
      
      FUNCTION dtau_drho(n1,n2,ng)
      
      use mphase_props
      
      imode1=iSysMode(n1,ngrnph(ng))
      imode2=iSysMode(n2,ngrnph(ng))      
      iph1 = ngrnph(ng)
      
      dtau_drho = 0.0
      
      !find the square root

      sqrt_root = sqrt(sum_latent(1,n1,ng))
      
      !latent hardening
      if (n1.ne.n2.and.iOpposite(n1).ne.n2) then 
          dtau_drho = chi_inter(iph1)*alatent(n1,n2,ngrnph(ng))
     #         *burg(imode1,iph1)*shearmod(imode1,ngrnph(ng))/
     #           (2.*sqrt_root)   
      !self hardening
      else 
          dtau_drho = chi_inter(iph1)
     #         *burg(imode1,iph1)*shearmod(imode1,ngrnph(ng))/
     #           (2.*sqrt_root)  
      endif     
      
c      if (n1.ne.n2) dtau_drho = 0.0
          
      END FUNCTION
      
      FUNCTION HP_CG(ng,iph1,ns)      
      
      use mphase_props 
      use sample_props_v, only : ngParent
      use twinning_v
      
      !HP_CG - initial deinition
      HP_CG=0.0 
      
      if (ng.le.ngParent) then 
          if (IPTSGRC(ng).gt.0) then
              iPTS_abs=IPTSGRC(ng)
              iPTS=IPTSGRC(ng)-nslsys(ngrnph(ng))
              imod1=iSysMode(ns,ngrnph(ng))
              imod2=iSysMode(iPTS_abs,ngrnph(ng))
              imodet=iSysMode(iPTS_abs,ngrnph(ng))-nslmod(ngrnph(ng))
              d_c=grsze*TWTHRES(1,imod2,iph1)             
              HP1fac=HPKCG(imod1,imodet,iph1)*sqrt(burg(imod1,iph1))
     #           *shearmod(imod1,ngrnph(ng))
              if (ns.ne.iPTS_abs) then
              !HP_CG - for parents
              HP_CG=(TWFRSY(iPTS,ng)-PTS_THRES(imodet,iph1))/ 
     #         (TWFR_THRES(imodet,iph1)-PTS_THRES(imodet,iph1))*(HP1fac/
     #         sqrt(d_c*(1.0-TWFRSY(iPTS,ng))*sysmfp(ns,iPTS)))    
              endif
          endif
      else
          ngp=iParentGrain(ng)  
          ntp=iParentSystem(ng)-nslsys(ngrnph(ng))
          ntp_abs=iParentSystem(ng)
          imtp=iParentMode(ng)
          imod2=iSysMode(ntp_abs,ngrnph(ng))  
          imod1=iSysMode(ns,ngrnph(ng))
          imodt=imtp-nslmod(ngrnph(ng))
          d_c=grsze*TWTHRES(1,imod2,iph1)
          HP1fac=HPKCG(imod1,imodt,iph1)*sqrt(burg(imod1,iph1))
     #          *shearmod(imod1,ngrnph(ng)) 
          !HP_CG - for twins
          HP_CG= HP1fac/sqrt(TWFRSY(ntp,ngp)*d_c*sysmfp(ns,ntp)) 
      endif

      RETURN
      END FUNCTION
          
      FUNCTION HP_CG_THRES(ng,iph1,ns)      

      use mphase_props 
      use twinning_v
      use sample_props_v, only : ngParent
      
      !HP_CG_THRES - initial definition
      HP_CG_THRES=0.0 
      
      if (ng.le.ngParent) then 
          if (IPTSGRC(ng).gt.0) then
              iPTS_abs=IPTSGRC(ng)
              iPTS=IPTSGRC(ng)-nslsys(ngrnph(ng))
              imod1=iSysMode(ns,ngrnph(ng))
              imod2=iSysMode(iPTS_abs,ngrnph(ng))
              imodet=iSysMode(iPTS_abs,ngrnph(ng))-nslmod(ngrnph(ng))
              d_c=grsze*TWTHRES(1,imod2,iph1)             
              HP1fac=HPKCG(imod1,imodet,iph1)*sqrt(burg(imod1,iph1))
     #           *shearmod(imod1,ngrnph(ng))
              !HP_CG_THRES - initial definition
              HP_CG_THRES=(PTS_THRES(imodet,iph1)-PTS_THRES(imodet
     #         ,iph1))/ (TWFR_THRES(imodet,iph1)-PTS_THRES(imodet,iph1)
     #         )*(HP1fac/sqrt(d_c*(1.0-PTS_THRES(imodet,iph1))
     #         *sysmfp(ns,iPTS))) 
          endif
      else
          ngp=iParentGrain(ng)  
          ntp=iParentSystem(ng)-nslsys(ngrnph(ng))
          ntp_abs=iParentSystem(ng)
          imtp=iParentMode(ng)
          imod2=iSysMode(ntp_abs,ngrnph(ng))  
          imod1=iSysMode(ns,ngrnph(ng))
          imodt=imtp-nslmod(ngrnph(ng))
          d_c=grsze*TWTHRES(1,imod2,iph1)
          HP1fac=HPKCG(imod1,imodt,iph1)*sqrt(burg(imod1,iph1))
     #       *shearmod(imod1,ngrnph(ng))  
          
          !HP_CG_THRES - initial definition
          HP_CG_THRES= HP1fac/sqrt(PTS_THRES(imodt,iph1)*d_c
     #      *sysmfp(ns,ntp))
      endif

      RETURN
      END FUNCTION
      
      SUBROUTINE hl_twin_mfp(tau)
        
      use const
      use flags, only : iDtwMfp
      use twinning_v
      use mphase_props
      use mvoigt
      use grain_props
      use grain_state, only : stcs
      use grain_rate, only : gamd
      use sample_props_v, only : ngrain,ngParent
      
      real,intent(out) :: tau(NSLS,NGR)
      
      DIMENSION i_sl_chg(NGR),dmfp(NSLS)
      DIMENSION tauHP_dif_temp(NSLS,NGR),auxSL(NSLS),aux(NSLS)
      DIMENSION dtauHP(NSLS),rss(NSLS)
      
c*******************************************************************************      
c     This is for change of CRSS of twins due to CG
c*******************************************************************************      
      ! Find derivative dtauHP_dgamma and increment dtauHP
      do ng=ngParent+1,ngrain !+ over all twin grains
        iph=ngrnph(ng)
c        if (igr_exist(ng).eq.1) then!if grain hasn't detwinned
        if (newGR(ng).eq.1) then !when new grain is created softening is not considered    
        ! Identify parent grain and twin system.  
        ngp=iParentGrain(ng)  
        ntp=iParentSystem(ng)-nslsys(iph)
        ntp_abs=iParentSystem(ng)
        imtp=iParentMode(ng)
        
        if (TWFRSY(ntp,ngp).gt.PTS_THRES(imtp-nslmod(iph)
     #      ,iph)) then!only for twin larger than IPTS thres
        do ns1=1,nsys(iph) !++ over all slip systems
          iflag=0
          auxSL(ns1)=0.0 
          dtauHP(ns1)=0.0
          
          imode1=iSysMode(ns1,ngrnph(ng))
          imode2=iSysMode(ntp_abs,ngrnph(ng))
          imodet2=imode2-nslmod(iph)
          
          d_c=grsze*TWTHRES(1,imode2,iph)
          HP1fac=HPKCG(imode1,imodet2,iph)*sqrt(burg
     #         (imode1,iph))*shearmod(imode1,iph)  
          !derivaive
          auxSL(ns1)=
     #         (  -0.5*HP1fac/(TWFRSY(ntp,ngp)*d_c*sysmfp(ns1,ntp))
     #         **1.5*d_c*sysmfp(ns1,ntp)   )
     #         *wgt(ngp)/wgtcg_all(ngp)/stw(imtp,iph)              
            
          !HP effect on detwin sys canceled
          if (ns1.eq.ntp_abs) auxSL(ns1)=0.0
          !increment
          dtauHP(ns1)=auxSL(ns1)*gamd(ntp_abs,ngp)
          !applay increment to CRSS
          crss=tau(ns1,ng)+dtauHP(ns1)
          !find RSS
          rss(ns1)=0.0
          do i=1,6
            rss(ns1)=rss(ns1)+mcs(i,ns1,ng)*stcs(i,ng)*profac(i)
          enddo
          !compare RSS with CRSS and correct increment
          if (rss(ns1).gt.crss) then
c            dtauHP(ns1)=rss(ns1)-crss
            dtauHP(ns1)=0.0              
            if (iTwinSys(ns1).eq.0) then
              if (iOpposite(ns1).lt.ns1) then
                  tau(iOpposite(ns1),ng)=tau(iOpposite(ns1),ng)-
     #              dtauHP(iOpposite(ns1))
              else
                  iflag=1
              endif
            endif
          endif
          !cancel increment if original tauHP is larger than current tauHP_CG
          if (tauHP(imode1,iph).gt.HP_CG(ng,iph,ns1)) then
              dtauHP(ns1)=0.0
          endif
          !if opposite system was out of ys cancel increment
          if (iflag.eq.1) then
              dtauHP(ns1)=0.0
          endif
          !for detwin sys
          if (iParentSystem(ng).eq.ns1.and.iDtwMfp.eq.1) then
              dtauHP(ns1)=0.0
          endif
          !calculate tau (CRSS) with correct increment
          if (dtauHP(ns1).le.0.0) then !USED ONLY FOR SOFTENING (REST IN HARDENING MATRIX)
              tau(ns1,ng)=tau(ns1,ng)+dtauHP(ns1) !COMMENT TO EXCLUDE SOFTENING OF HP IN TWINS
          endif
        enddo !++ over all slip systems
        endif !only for twins larger than IPTSTHRES
c        if (ng.eq.115) write(71,'(10f15.4)') (auxSL(i),i=1,6),
c     #   gamd(ntp_abs,ngp),(dtauHP(i),i=1,6),(HP_CG(115,i),i=1,6)
        endif
c        endif
      enddo!+ over all twin grains
c-------------------------------------------------------------------------------      
      ! Correct CRSS of newly formed twins.
      do ng=ngParent+1,ngrain !+ over all twin grains
c        if (igr_exist(ng).eq.1) then  
        ngp=iParentGrain(ng) 
        imtp=iParentMode(ng)  
        ntp=iParentSystem(ng)-nslsys(iph)
        if (newGR(ng).eq.0) then !(TWFRSY(ntp,ngp).eq.TwinFrac(imtp)) then !change this condition for when twins are first created
          newGR(ng)=1
          do ns=1,nsys(iph)
            if (iParentSystem(ng).ne.ns) then !exclude hardening of detwin sys
              imode=iSysMode(ns,ngrnph(ng))
              if (HP_CG_THRES(ng,iph,ns).gt.tauHP(imode,iph)) then
                tau(ns,ng)=tau(ns,ng)-tauHP(imode,iph)
     #           +HP_CG_THRES(ng,iph,ns)
              endif
            endif
          enddo
        endif
c        endif
      enddo !+ over all twin grains
c*******************************************************************************      
      !this changes CRSS of parents due to detwinning
c*******************************************************************************
      if (iDtwMfp.eq.1) then      
      do ng=1,ngParent !+ over all parent grains
       dtauHP(:)=0.0   
       if (IPTSGRC(ng).ne.0.) then 
        n2=IPTSGRC(ng)  
        ngc=iChildGrain(n2,ng)   
        if (gamd(n2,ngc).ne.0.0) then
        nt2=n2-nslsys(iph)
        do ns1=1,nsys(iph) !++ over all slip systems
          auxSL(ns1)=0.0 
          dtauHP(ns1)=0.0
          
          imode1=iSysMode(ns1,ngrnph(ng))
          imode2=iSysMode(n2,ngrnph(ng))
          imodet2=imode2-nslmod(iph)
          
          d_c=grsze*TWTHRES(1,imode2,iph)
          HP1fac=HPKCG(imode1,imodet2,iph)*sqrt(burg
     #         (imode1,iph))*shearmod(imode1,iph)   
          !derivaive
          auxSL(ns1)=
     #      (  1.0/((TWFR_THRES(imodet2,iph)-PTS_THRES(imodet2,iph)))
     #      *HP1fac/((1.0-TWFRSY(nt2,ng))*d_c*sysmfp(ns1,nt2))**0.5
     #      + 0.5*(TWFRSY(nt2,ng)-PTS_THRES(imodet2,iph))
     #      /(TWFR_THRES(imodet2,iph)-PTS_THRES(imodet2,iph))
     #      *HP1fac/((1.0-TWFRSY(nt2,ng))*d_c*sysmfp(ns1,nt2))
     #      **1.5*d_c*sysmfp(ns1,nt2)   )
     #      *wgt(ngc)/wgtcg_all(ngc)/stw(imode2,iph)            
            
          !increment
          dtauHP(ns1)=-auxSL(ns1)*gamd(n2,ngc)
          !applay increment to CRSS
          crss=tau(ns1,ng)+dtauHP(ns1)
          !find RSS
          rss(ns1)=0.0
          do i=1,6
            rss(ns1)=rss(ns1)+mcs(i,ns1,ng)*stcs(i,ng)*profac(i)
          enddo
          !compare RSS with CRSS and correct increment
          if (rss(ns1).gt.crss) then
c              dtauHP(ns1)=rss(ns1)-crss
              dtauHP(ns1)=0.0
          endif
          !cancel increment if original tauHP is larger than current tauHP_CG
          if (tauHP(imode1,iph).gt.HP_CG(ng,iph,ns1)) then
              dtauHP(ns1)=0.0
          endif
c          if (ng.eq.55.and.ns1.eq.2) write(71,'(5f9.2)',advance='no')
c     #     TWFRSY(nt2,ng),HP_CG(ng,ns1),tauHP(imode1),rss(ns1),crss
c          if (ng.eq.55) write(71,'(f9.1)',advance='no') dtauHP(ns1) !auxSL(ns1)
c          if (ng.eq.55.and.ns1.eq.30) write(71,*)          
          !compare rss
          !calculate tau (CRSS) with correct increment
          if (dtauHP(ns1).lt.0.0) then
              tau(ns1,ng)=tau(ns1,ng)+dtauHP(ns1) !COMMENT TO EXCLUDE SOFTENING IN PARENTS DUE TO DETWINNING
          endif      
        enddo !++ over all slip systems
        endif      
      endif
c      if (ng.eq.55) write(71,'(30f9.1)') !dtauHP(1:30)
      enddo !+ over all parent grains
      endif   
      RETURN
      END SUBROUTINE
      
      SUBROUTINE hl_pts_damp(iproc,istep,tau)
      
      use mphase_props
      use sample_props_v
      use twinning_v
      use flags, only : iTwPh
      
      integer, intent(in) :: iproc,istep
      real, intent(out) :: tau(NSLS,NGR)
      
      !PTS IS BLOCKED IN PARENT WHEN THRES2 IS REACHED
      !ALL THE TWIN SYSTEMS ARE BLOCKED WHEN THE SUM OF THEIR VOL FRACTION REACHES 1 
      do ng=1,ngrain
          iph=ngrnph(ng)
          if (ng.le.ngParent.and.iTwPh(ngrnph(ng)).eq.1) then !exclude twins frop PTDAMP effect
              !Calculation of overall accumulated twin volume fraction in a grain
              !Calculation of accumulated volume fraction of predominant twin system            
              TWFRGRAIN=0.
              DO IT=1,ntwsys(ngrnph(ng))
                  TWFRGRAIN=TWFRGRAIN+TWFRSY(IT,ng)/2.0
              ENDDO
              
              do iPTS=1,ntwsys(ngrnph(ng)) !MTV
c                  if (KTWSMX(iPTS,ng).gt.0) then
c                      TWFRPTS=TWFRSY(KTWSMX(iPTS,ng)-nslsys,ng) !MTV
c                  endif
              enddo !MTV
              !Threshold defines limit in TWFRSY, when further volume transfere is stopped. 
              !This parameter is read from single crystal data file for each mode. 
              !Calculation of PTSDAMP based on the TWFR_THRES and TWFRGRAIN.
              do iPTS=1,ntwsys(ngrnph(ng)) !MTV      
                  imode=iSysMode(iPTS+nslsys(iph),iph)-nslmod(iph) !defines mode of twin system (relative to twin modes)
                  is=iPTS+nslsys(iph)
                  TWFRPTS=TWFRSY(iPTS,ng) !MTV - value of TVF on sys
                  PTSDAMP(is,ng)=1.                            
                  PTSDAMP(is,ng)=1.+(TWFRPTS/TWFR_THRES(imode,iph))**20 
     #                 +(TWFRGRAIN/0.95)**100 
                  tau(is,ng)=tau(is,ng) !MTV
     #                *PTSDAMP(is,ng)/PTSDAMP_OLD(is,ng) !MTV
                  PTSDAMP_OLD(is,ng)=PTSDAMP(is,ng) !MTV                            
              enddo !MTV   
          endif
      enddo ! over grains    

      RETURN
      END SUBROUTINE
      
      end module hard_law1
c      
c***********************************************************************
c *** hard_law2 ********************************************************
c***********************************************************************
      module hard_law2
      
      use const
      use hard_law1_v, only : iDirSys,iOpposite,iact_sys,itemp_act_sys
     #    ,gam_acc,alpha_reg,alatent,d_mod,g_mod
      
      real CRSS_0(NMOD,NPHM),hard(NMOD,NPHM),aLat(NSLS,NSLS)   
      real dgamma(NSLS),taus(NSLS),taul(NSLS),taut(NSLS)
      real tautd(NSLS),Q(NMOD,NPHM),b(NMOD,NPHM),d0,d1,d2,d3,d4,d5
      
      !state variables definition
      type state_hl
          real :: tau(NSLS),gamtot,gam_acc(NSLS)
      end type state_hl 
      
      private
      
      public hl_read,hl_ini_tau,hl_update_statv,hl_hd,gam_acc,state_hl
     #    ,iact_sys,iOpposite,iDirSys,set_hl,get_hl,write_statv_hl
     #    ,read_statv_hl,itemp_act_sys
      contains
      
      !hl_read: reads parameters for hardening law
      !hl_ini_tau: calculates initial CRSS and assigns it to tau
      !hl_update_statv: updates state variables of hardening law
      !hl_hd: calculates hardening matrix
      !set_hl: sets statv object
      !get_hl: gets values from statv object
      !write_statv_hl: writes statv object to statv array
      !read_statv_hl: reads statv object from statv array
      
      subroutine set_hl(ng,hl)
      
      use grain_state, only : tau,gamtot
      
      integer, intent(in) :: ng
      type(state_hl),intent(out) :: hl
      
      hl%tau=tau(:,ng)
      hl%gamtot=gamtot(ng)
      hl%gam_acc=gam_acc(:,ng)
      
      end subroutine set_hl
      
      subroutine get_hl(ng,hl)
      
      use grain_state, only : tau,gamtot
      
      integer, intent(in) :: ng
      type(state_hl),intent(in) :: hl
      
      tau(:,ng)=hl%tau
      gamtot(ng)=hl%gamtot
      gam_acc(:,ng)=hl%gam_acc
      
      end subroutine get_hl
      
      subroutine write_statv_hl(hl,ng,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_hl),intent(in) :: hl
      integer,intent(in) :: ng,NSTATV
      integer,intent(inout) :: ns
      real,intent(out) :: STATEV(NSTATV)
      
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%tau(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns)=hl%gamtot
      ns=ns+1
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=hl%gam_acc(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      
      end subroutine write_statv_hl
      
      subroutine read_statv_hl(hl,ng,iph,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_hl),intent(out) :: hl
      integer,intent(in) :: ng,NSTATV
      real,intent(in) :: STATEV(NSTATV)
      integer,intent(inout) :: ns
      
      hl%tau(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      hl%gamtot=STATEV(ns)
      ns=ns+1
      hl%gam_acc(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      
      end subroutine read_statv_hl      
      
      SUBROUTINE hl_read(iph)
      
      use const
      use mphase_props
      use flags, only : kCL, iDiag
      use twinning_v, only : TwinFrac
      
      !initialize iact arrays
      call initialize_iact(iph)
      
      read(1,*)
      do imo=1,nslmod(iph)
          read(1,*) CRSS_0(imo,iph),Q(imo,iph),b(imo,iph)
      enddo
      do imo=nslmod(iph)+1,nmodes(iph)
          read(1,*) CRSS_0(imo,iph),Q(imo,iph),b(imo,iph)
     #        ,TwinFrac(imo,iph)
      enddo
      read(1,*) iDiag,alpha_reg
      
      !latent hardening
      READ(1,*) iFCC
      if (iFCC.ne.1) then
          do is1=1,nslmod(iph)
              READ(1,*) (d_mod(is1,is2,iph),is2=1,(nslmod(iph)+1))
          enddo
      else
          READ(1,*) d0,d1,d2,d3,d4,d5
      endif
      
      !initialize latent hardening matrix
      if (iFCC.eq.0) then
          call make_alat(iph)
      else
          call make_alat_FCC(iph)
      endif
      
      do ns=1,nsys(iph)
          if(iTwinSys(ns).eq.1)then
              alatent(ns,ns,iph)=1.0
          endif
      enddo
      
      RETURN
      END SUBROUTINE     
      
      SUBROUTINE hl_ini_tau(ng,tau,iopt)
      
      use const
      use files, only : fileproc
      use mphase_props
      
      integer,intent (in) :: ng
      real, intent(out) :: tau(NSLS,NGR)
      
      is=0              
      do imo=1,nmodes(ngrnph(ng))    ! first we consider the slip modes
        do isy=1,nsm(imo,ngrnph(ng))
          is=is+1              
          tau(is,ng)=CRSS_0(imo,ngrnph(ng))
        enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE hl_update_statv(ng)
      
      use const
      use mphase_props
      use grain_state, only : iact,nact
      use grain_rate, only : gamd
      
      integer, intent(in) :: ng
      
      do ns=1,nslsys(ngrnph(ng))
          gam_acc(ns,ng)=gam_acc(ns,ng)+gamd(ns,ng)
      enddo      
      
c     Set up array of active systems.
      do is=1,nact(ng)
          n1=iact(is,ng)
          if (iTwinSys(n1).EQ.0) then
              if (mod(n1,2).EQ.1) then
                  iact_sys(iDirSys(n1,ngrnph(ng)),ng)=1
              else
                  iact_sys(iDirSys(n1,ngrnph(ng)),ng)=2
              endif
          else
              iact_sys(iDirSys(n1,ngrnph(ng)),ng)=1
          endif
      enddo  
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE hl_hd(ng,hd)
      
      use const     
      use mphase_props
      use grain_state, only : nact,iact,gamtot
      use grain_rate, only : gamd
      
      !in
      integer, intent(in) :: ng
      
      !out
      real, intent(out) :: hd(NSLS,NSLS)
      
      hd(:,:)=0. 
      !total accumulated shear strain
      aux=gamtot(ng)
      do ns1=1,nact(ng)
        n1=iact(ns1,ng)          
        aux=aux+dgamma(n1)
      enddo  
      !Define the hardening matrix for the nact(ng) active systems
      do n1=1,nsys(ngrnph(ng)) ! goes over all sys, not only active
        imode1=iSysMode(n1,ngrnph(ng))
        do ns2=1,nact(ng)
          n2=iact(ns2,ng)
          imode2=iSysMode(n2,ngrnph(ng))
          dtaus_dgamma=b(imode2,ngrnph(ng))*Q(imode2,ngrnph(ng))
     #            *exp(-b(imode2,ngrnph(ng))*(gam_acc(n2,ng))) !decoupled hard !explicite +dgamma(n2)
c          dtaus_dgamma=b*Q*exp(-b*aux) !coupled hard
          !define elements of hd
          hd(n1,n2)=dtaus_dgamma
     #            *alatent(n1,n2,ngrnph(ng))
c          if (hd(n2,n2).lt.1e-3) hd(n2,n2)=1e-3
        enddo
      enddo
c     Ensure all components of the hardening matrix are larger than 1.0
      do ns1=1,nact(ng)
        n1=iact(ns1,ng)
        do ns2=1,nact(ng)
          n2=iact(ns2,ng)
            if(hd(n1,n2).LT.1.0) then 
                hd(n1,n2)=1.0
                if (iTwinSys(n1).NE.1) then
                    hd(iOpposite(n1),n2)=1.0
                endif
            endif
c          if(ng.eq.1) write(*,'(E8.2,'' '')',advance='no') hd(n1,n2)              
        enddo
      enddo     
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE initialize_iact(iph)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
c      SUBROUTINE ACTIVE_SYS(ng,iopt)      
c      initiates arrays used for determination of currently active sys, disloc
c         dens, CRSS.
c         iDirSys(is,nph) - for given index of duble system (is) returns value of direct system
c         iOpposite(is) - for given value of double slip system (is) returns value of 
c                         opposite system (system with positive/negative shear)        
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
      use mphase_props

c     Set up direct system from double.          
      is=0
      do is1=1,nslsys(iph)/2
          do is2=1,2
              is=is+1
              iDirSys(is,iph)=is1
          enddo    
      enddo
      do it=1,ntwsys(iph)
          is=is+1
          iDirSys(is,iph)=nslsys(iph)/2+it
      enddo    
c     Set up opposite system.
      do is=1,nslsys(iph)
          if (mod(is,2).EQ.1) then
              iOpposite(is)=is+1
          else 
              iOpposite(is)=is-1
          endif    
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE make_alat_FCC(iph)
      
      use mphase_props_v, only : nslsys
      real :: aLat(12,12)
      
      !CRSS latent
        aLat(1,1)= d0      !A2
        aLat(1,2)= d1      !A3
        aLat(1,3)= d1      !A6
        aLat(1,4)= d3      !B2
        aLat(1,5)= d4      !B4
        aLat(1,6)= d4      !B5
        aLat(1,7)= d2      !C1
        aLat(1,8)= d4      !C3
        aLat(1,9)= d5      !C5
        aLat(1,10)=d2      !D1
        aLat(1,11)=d5      !D4
        aLat(1,12)=d4      !D6
        
        !A3!A3
        aLat(2,2)= d0      !A3
        aLat(2,3)= d1      !A6
        aLat(2,4)= d4      !B2
        aLat(2,5)= d2      !B4
        aLat(2,6)= d5      !B5
        aLat(2,7)= d4      !C1
        aLat(2,8)= d3      !C3
        aLat(2,9)= d4      !C5
        aLat(2,10)=d5      !D1
        aLat(2,11)=d2      !D4
        aLat(2,12)=d4      !D6        
        
        !A6!A6
        aLat(3,3)= d0      !A6
        aLat(3,4)= d5      !B2
        aLat(3,5)= d5      !B4
        aLat(3,6)= d2      !B5
        aLat(3,7)= d5      !C1
        aLat(3,8)= d4      !C3
        aLat(3,9)= d2      !C5
        aLat(3,10)=d4      !D1
        aLat(3,11)=d4      !D4
        aLat(3,12)=d3      !D6          
        
        !B2!B2
        aLat(4,4)= d0      !B2
        aLat(4,5)= d1      !B4
        aLat(4,6)= d1      !B5
        aLat(4,7)= d2      !C1
        aLat(4,8)= d5      !C3
        aLat(4,9)= d4      !C5
        aLat(4,10)=d2      !D1
        aLat(4,11)=d4      !D4
        aLat(4,12)=d5      !D6  
        
        !B4!B4
        aLat(5,5)= d0      !B4
        aLat(5,6)= d1      !B5
        aLat(5,7)= d5      !C1
        aLat(5,8)= d2      !C3
        aLat(5,9)= d4      !C5
        aLat(5,10)=d4      !D1
        aLat(5,11)=d3      !D4
        aLat(5,12)=d4      !D6  
        
        !B5!B5
        aLat(6,6)= d0      !B5
        aLat(6,7)= d4      !C1
        aLat(6,8)= d4      !C3
        aLat(6,9)= d3      !C5
        aLat(6,10)=d5      !D1
        aLat(6,11)=d4      !D4
        aLat(6,12)=d2      !D6  
        
        !C1!C1
        aLat(7,7)= d0      !C1
        aLat(7,8)= d1      !C3
        aLat(7,9)= d1      !C5
        aLat(7,10)=d3      !D1
        aLat(7,11)=d4      !D4
        aLat(7,12)=d4      !D6  
        
        !C3!C3
        aLat(8,8)= d0      !C3
        aLat(8,9)= d1      !C5
        aLat(8,10)=d4      !D1
        aLat(8,11)=d2      !D4
        aLat(8,12)=d5      !D6  
        
        !C5!C5
        aLat(9,9)= d0      !C5
        aLat(9,10)=d4      !D1
        aLat(9,11)=d5      !D4
        aLat(9,12)=d2      !D6  
        
        !D1!D1
        aLat(10,10)=d0     !D1
        aLat(10,11)=d1     !D4
        aLat(10,12)=d1     !D6  
                           
        !D1!D1             
        aLat(11,11)=d0     !D4
        aLat(11,12)=d1     !D6  
                           
        !D1!D1             
        aLat(12,12)=d0     !D6          
        
        do i=1,nslsys(iph)/2
        do j=1,nslsys(iph)/2
          if(j.lt.i) aLat(i,j) = aLat(j,i)
        enddo
        enddo   
        
        !make alatent
        do i=1,nslsys(iph)
        do j=1,nslsys(iph)
          alatent(i,j,iph)=alat(iDirSys(i,iph),iDirSys(j,iph))
        enddo
        enddo
        
      RETURN
      END SUBROUTINE
      
      SUBROUTINE make_alat(iph)
      
      use mphase_props, only : iSysMode,nslsys
      
      do is1=1,nslsys(iph)
          im1=iSysMode(is1,iph)
          do is2=1,nslsys(iph)
              im2=iSysMode(is2,iph)
              if (is1.eq.is2.or.is1.eq.iOpposite(is2)) then
                  alatent(is1,is2,iph)=d_mod(im1,1,iph)
              else
                  alatent(is1,is2,iph)=d_mod(im1,im2+1,iph)
              endif
          enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      end module hard_law2
c      
c***********************************************************************
c *** twinning *********************************************************
c***********************************************************************
      module twinning
      
      use twinning_v
      
      contains
      
      !set_hl: sets statv object
      !get_hl: gets values from statv object
      !write_statv_hl: writes statv object to statv array
      !read_statv_hl: reads statv object from statv array
      
      subroutine set_tw(ng,tw)
      
      integer, intent(in) :: ng
      type(state_tw),intent(out) :: tw
      
      tw%TwFrSy=TwFrSy(:,ng)
      tw%wgtcg_all=wgtcg_all(ng)
      tw%iParentGrain=iParentGrain(ng)
      tw%iParentMode=iParentMode(ng)
      tw%iParentSystem=iParentSystem(ng)
      tw%iPTSgr=iPTSgr(:,ng)
      tw%iPTSgrc=iPTSgrc(ng)
      tw%iChildGrain=iChildGrain(:,ng)
      tw%PTSdamp=PTSdamp(:,ng)
      tw%newGR=newGR(ng)
      tw%newGRSh=newGRSh(ng)
      tw%axis3=axis3(ng)
c      tw%tau0_tw=tau0_tw(:,ng)
      
      end subroutine set_tw
      
      subroutine get_tw(ng,tw)
      
      integer, intent(in) :: ng
      type(state_tw),intent(in) :: tw
      
      TwFrSy(:,ng)=tw%TwFrSy
      wgtcg_all(ng)=tw%wgtcg_all
      iParentGrain(ng)=tw%iParentGrain
      iParentMode(ng)=tw%iParentMode
      iParentSystem(ng)=tw%iParentSystem
      iPTSgr(:,ng)=tw%iPTSgr
      iPTSgrc(ng)=tw%iPTSgrc
      iChildGrain(:,ng)=tw%iChildGrain
      PTSdamp(:,ng)=tw%PTSdamp
      newGR(ng)=tw%newGR
      newGRSh(ng)=tw%newGRSh
      axis3(ng)=tw%axis3
c      tau0_tw(:,ng)=tw%tau0_tw
      
      end subroutine get_tw
      
      subroutine clean_tw_state(ng)
      
      integer, intent(in) :: ng
      
      TwFrSy(:,ng)=0.0
      wgtcg_all(ng)=0.0
      iParentGrain(ng)=0
      iParentMode(ng)=0
      iParentSystem(ng)=0
      iPTSgr(:,ng)=0
      iPTSgrc(ng)=0
      iChildGrain(:,ng)=0
      PTSdamp(:,ng)=0.0
      newGR(ng)=0.0
	  newGRSh(ng)=0
c      tau0_tw(:,ng)=tw%tau0_tw
      
      end subroutine clean_tw_state
      
      subroutine write_statv_tw(tw,ng,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_tw),intent(in) :: tw
      integer,intent(in) :: ng,NSTATV
      integer,intent(inout) :: ns
      real,intent(out) :: STATEV(NSTATV)
      
      STATEV(ns:ns+ntwsys(ngrnph(ng))-1)=tw%TwFrSy(1:ntwsys(ngrnph(ng)))
      ns=ns+ntwsys(ngrnph(ng))
      STATEV(ns)=tw%wgtcg_all
      ns=ns+1
      STATEV(ns)=tw%iParentGrain
      ns=ns+1
      STATEV(ns)=tw%iParentMode
      ns=ns+1
      STATEV(ns)=tw%iParentSystem
      ns=ns+1
      STATEV(ns:ns+ntwsys(ngrnph(ng))-1)=tw%iPTSgr(1:ntwsys(ngrnph(ng)))
      ns=ns+ntwsys(ngrnph(ng))
      STATEV(ns)=tw%iPTSgrc
      ns=ns+1
      STATEV(ns:ns+ntwsys(ngrnph(ng))-1)=
     #    tw%iChildGrain(1+nslsys(ngrnph(ng)):nsys(ngrnph(ng)))
      ns=ns+ntwsys(ngrnph(ng))
      STATEV(ns:ns+ntwsys(ngrnph(ng))-1)=
     #    tw%PTSdamp(1+nslsys(ngrnph(ng)):nsys(ngrnph(ng)))
      ns=ns+ntwsys(ngrnph(ng))
      STATEV(ns)=tw%newGR
      ns=ns+1
      STATEV(ns)=tw%newGRSh
      ns=ns+1
      STATEV(ns)=tw%axis3
      ns=ns+1
      
      end subroutine write_statv_tw
      
      subroutine read_statv_tw(tw,ng,iph,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_tw),intent(out) :: tw
      integer,intent(in) :: ng,NSTATV,iph
      real,intent(in) :: STATEV(NSTATV)
      integer,intent(inout) :: ns

      tw%TwFrSy(1:ntwsys(iph))=STATEV(ns:ns+ntwsys(iph)-1)
      ns=ns+ntwsys(iph)
      tw%wgtcg_all=STATEV(ns)
      ns=ns+1
      tw%iParentGrain=STATEV(ns)
      ns=ns+1
      tw%iParentMode=STATEV(ns)
      ns=ns+1
      tw%iParentSystem=STATEV(ns)
      ns=ns+1
      tw%iPTSgr(1:ntwsys(iph))=STATEV(ns:ns+ntwsys(iph)-1)
      ns=ns+ntwsys(iph)
      tw%iPTSgrc=STATEV(ns)
      ns=ns+1
      tw%iChildGrain(1+nslsys(iph):nsys(iph))=
     #    STATEV(ns:ns+ntwsys(iph)-1)
      ns=ns+ntwsys(iph)
      tw%PTSdamp(1+nslsys(iph):nsys(iph))=
     #    STATEV(ns:ns+ntwsys(iph)-1)
      ns=ns+ntwsys(iph)
      tw%newGR=STATEV(ns)
      ns=ns+1
      tw%newGRSh=STATEV(ns)
      ns=ns+1
      tw%axis3=STATEV(ns)
      ns=ns+1
      
      end subroutine read_statv_tw
      
      SUBROUTINE initialize_tw(ngrain)
      
      use const
      use grain_props_v, only : wgt
      use sample_props_v, only : ngParent
      
      do ng=1,ngrain
        !Setting the flag for child grains to zero for all deformation modes (including non-twin systems for simplicity)
        do i=1,NSLS
            iChildGrain(i,ng)=0
        enddo
        iParentGrain(ng)=0  !the array that gives the parent for each twin
        iParentSystem(ng)=0 !the array that gives the Twin system in the parent grain for each twin
        !Setting the counter for twin level (parent, 1st twin, secondary twinning, ..., etc) to zero for the parent grains
        iTwinLevel(ng)=0      
      enddo
      ngParent=ngrain
      do ng=1,ngrain
        wgtcg_all(ng)=wgt(ng)
        PTSDAMP_OLD(:,ng)=1.
        PTSDAMP(:,ng)=1.
      enddo
      
      RETURN
      END SUBROUTINE
      
      subroutine twin_fraction(step)
      !update twfrsy, twfrgr, iptsgr
      use const
      use sample_props_v, only : ngrain, ngParent
      use mphase_props
      use flags, only : iDeTwOpt,iTwPh
      use grain_state, only : tau,iact,nact
      use grain_props_v, only : wgt
      use grain_rate, only : gamd
      
      dimension aux1(NSLS),iaux1(NSLS)
      
      do ng=1,ngrain     
        ! NEW EVOLUTION
        ! only active twin modes need to be considered
        if(iTwPh(ngrnph(ng)).eq.1) then
          do ns1=1,nact(ng)
            n1=iact(ns1,ng)
            if (iTwinSys(n1).eq.1) then
              its=n1-nslsys(ngrnph(ng))
              itm=iSysMode(n1,ngrnph(ng))
              !Calculation of accumulated twin volume fraction.                    
              GSYST=gamd(n1,ng)*step
              GMODE=GMODE+GSYST
              TwFrInc=GSYST/stw(itm,ngrnph(ng))*wgt(ng)/WGTCG_all(ng)
              if(mod(n1,2).eq.0) then 
                TwFrInc=-TwFrInc
                TwFrSy(its-1,ng)=TwFrSy(its-1,ng)+TwFrInc !accumulated twin volume fraction per sys per grain
                if (TwFrSy(its-1,ng).lt.0.0) TwFrSy(its-1,ng)=0.0
                TwFrSy(its,ng)=TwFrSy(its-1,ng)
              else
                TwFrSy(its,ng)=TwFrSy(its,ng)+TwFrInc !accumulated twin volume fraction per sys per grain
                TwFrSy(its+1,ng)=TwFrSy(its,ng)
              endif
              TwFrGr(itm,ng)=TwFrGr(itm,ng)+TwFrInc !accumulated twin volume fraction per mode per grain
              !Detwinning - transfers volume fraction back to parent if pts for detwinning is active in child
              if (iDeTwOpt.eq.1) then
                if (iParentGrain(ng).gt.0) then
c                  do iPTS=1,ntwsys !MTV
                  if (n1.eq.iParentSystem(ng))
     #             then
                      TwFrSy(its,iParentGrain(ng))=
     #                TwFrSy(its,iParentGrain(ng))-TwFrInc
                     if (TwFrSy(its,iParentGrain(ng)).lt.0.) then
                          TwFrSy(its,iParentGrain(ng))=0.0
c                          tau(its_abs,ng)=10000.0
                      endif
                  endif
c                  enddo
                endif
              endif   
            endif
          enddo
        ! OLD EVOLUTION
c        TWFMAX=0.
c        its_abs=nslsys(ngrnph(ng))
c        itm_abs=nslmod(ngrnph(ng))
c        its=0
c        do itm=1,ntwmod(ngrnph(ng)) !goes over twinning modes
c          itm_abs=itm_abs+1 !absolute counter of modes
c          gmode=0.   
c          TWFRGR(ITM,ng)=0.  
c          TWSHX=stw(itm_abs) !characteristic twin shear - in EPSC its given for all modes
c          do IT=1,nsm(itm_abs,1) !goes over all systems in twin mode
c            its=its+1   ! counter for twin systems in all twin modes
c            its_abs=its_abs+1   ! absolute counter over all systems         
c            !Calculation of accumulated twin volume fraction.                    
c            GSYST=gamd(its_abs,ng)*step
c            GMODE=GMODE+GSYST
c            TwFrInc=GSYST/TWSHX*wgt(ng)/WGTCG_all(ng)
cc            if (ilog_dtw(its_abs,ng).eq.1) then
cc              TWFRINC=-TWFRINC !detwinning
cc            endif
c            TwFrSy(its,ng)=TwFrSy(its,ng)+TwFrInc !accumulated twin volume fraction per sys per grain
c            TWFRGR(ITM,ng)=TWFRGR(ITM,ng)+TwFrSy(its,ng) !accumulated twin volume fraction per mode per grain
c            rand=1.0 !to add stochasticity to selection - not implemented
c            !Detwinning - transfers volume fraction back to parent if pts for detwinning is active in child
c            if (iDeTwOpt.eq.1) then
c              if (iParentGrain(ng).gt.0) then
cc                do iPTS=1,ntwsys !MTV
c                if (its_abs.eq.iParentSystem(ng))
c     #           then
c                    TwFrSy(its,iParentGrain(ng))=
c     #              TwFrSy(its,iParentGrain(ng))-TWFRINC
c                   if (TwFrSy(its,iParentGrain(ng)).lt.0.) then
c                        TwFrSy(its,iParentGrain(ng))=0.0
c                        tau(its_abs,ng)=10000.0
c                    endif
c                endif
cc                enddo
c              endif
c            endif                    
c          enddo !over all systems in twin mode
c        enddo ! over twinning modes
        endif  
      enddo ! over all  grains
      
      !Searching for PTS & PTM       
      do ng=1,ngParent
        if(iTwPh(ngrnph(ng)).eq.1) then
          aux1=TwFrSy(:,ng)
          ntw=1
          do itws=1,ntwsys(ngrnph(ng))
              if (iChildGrain(itws+nslsys(ngrnph(ng)),ng).ne.0) then
                  aux1(itws)=100.0
                  ntw=ntw+1
              endif
          enddo
          call sorting(aux1,iaux1,ng)
          iPTSgr(ntw:NSLS,ng)=iaux1(ntw:NSLS)
        endif
      enddo
      !Finde IPTSGRC - array that defines one PTS per grain using condition
      !that PTS is the first twin to reach 5% VF.
      do ng=1,ngParent
        if(iTwPh(ngrnph(ng)).eq.1) then
          IPTS=IPTSGR(1,ng) !find current PTS in grain
          IPTM=iSysMode(IPTS,ngrnph(ng))
          if(IPTSGRC(ng).eq.0.and.TWFRSY(IPTS-nslsys(ngrnph(ng)),ng).ge.
     #        PTS_THRES(IPTM-nslmod(ngrnph(ng)),ngrnph(ng))) then
              IPTSGRC(ng)=IPTS
          endif
          
          if (IPTSGRC(ng).gt.0) then
              if (TWFRSY(IPTSGRC(ng)-nslsys(ngrnph(ng)),ng).
     #             lt.PTS_THRES(IPTM-nslmod(ngrnph(ng)),ngrnph(ng)))then
                  IPTSGRC(ng)=0
              endif
          endif    
        endif
      enddo
      
      PTVFM(:)=0.
      CTVFM(:)=0.        
      do ng=1,ngParent
        if(iTwPh(ngrnph(ng)).eq.1) then
          it1=1            
          do im=1,ntwmod(ngrnph(ng))
              do it=1,nsm(nslmod(ngrnph(ng))+im,ngrnph(ng))
                  PTVFM(im)=PTVFM(im)+TWFRSY(it1,ng)*wgtcg_all(ng)/2.0 !/2.0 twin in both directions
                  it1=it1+1                    
              enddo
          enddo   
        endif
      enddo
      do ng=ngParent+1,ngrain
        ngp=iParentGrain(ng)  
        it1=1
        do im=1,ntwmod(ngrnph(ng))
            do it=1,nsm(nslmod(ngrnph(ng))+im,ngrnph(ng))
                if (iParentSystem(ng).NE.it1+nslsys(ngrnph(ng))) then !exclude detwinning
                    CTVFM(im)=CTVFM(im)+TWFRSY(it1,ng)*wgtcg_all(ng)/2.0!/2.0 twin in both directions
                endif
                it1=it1+1
            enddo
        enddo            
      enddo
      
      end subroutine
      
      subroutine sorting(a,ind,ng)

      use const
      use mphase_props
      
      dimension a(NSLS),an(NSLS),ind(NSLS)
      
      an=-a
      nsize=ntwsys(ngrnph(ng))
      do irow = 1, nsize
          krow = minloc( an( 1:nsize ), dim=1 )
          ind(irow) = krow
          an(krow)=1000.0
      enddo
      !absolute index for iPTSgr
      ind(1:ntwsys(ngrnph(ng)))=ind(1:ntwsys(ngrnph(ng)))
     #    +nslsys(ngrnph(ng))
      
      end subroutine
        
      SUBROUTINE TWIN_BARRIERS
      !update sysmfp     
      use mphase_props
      use sample_props_v, only : nph
      
      DIMENSION sinalfa(NSLS,NSLS),bur(3),aux33(3,3),T_cr_tw(3,3)
      REAL ncc_tw(3,NSLS)
      
      !Calculates:
      !Angle between twin sys normal and slip sys normal
      !Find sinalfa, by doing a dot product between plane normals
      do ns=1,nsys(nph)
          do nt=1,ntwsys(nph)
              nt_abs=nslsys(nph)+nt
              cosalfa=ncc(1,ns,nph)*ncc(1,nt_abs,nph)+ncc(2,ns,nph)
     #             *ncc(2,nt_abs,nph)+ncc(3,ns,nph)*ncc(3,nt_abs,nph)
              if (abs(cosalfa).gt.0.9999) cosalfa=sign(0.9999,cosalfa)
              sinalfa(ns,nt)=(1.0-cosalfa**2)**0.5
              sysmfp(ns,nt)=1.0/(1.0-cosalfa**2)**0.5
          enddo
      enddo

      RETURN
      END SUBROUTINE
      
      subroutine wgtd_calc_tw(wgtd)
      
      use const
      use flags, only : iDeTwOpt,iTwPh
      use mphase_props
      use grain_props_v, only : wgt
      use grain_state, only : iact,nact,gamtot
      use grain_rate, only : gamd
      use sample_props_v, only : ngrain,ngParent
      
      real, intent(out) :: wgtd(NGR)
      
      LOGICAL nuc
      DIMENSION aux33(3,3),aux6(6)
     
      do ngP=1,ngParent
        if(iTwPh(ngrnph(ngP)).eq.1) then
          !restart parent weight
c          wgt(ngP)=wgtcg_all(ngp) !if wgt is updated directly
          !restart weight increment
          wgtd(ngP)=0.0
          if(sum(iChildGrain(1:ntwsys(ngrnph(ngP)),ngP)).ne.0)then
              do itws=1,ntwsys(ngrnph(ngP))
                  if(iChildGrain(itws,ngP).ne.0)then
                      wgtd(iChildGrain(itws,ngP))=0.0
                  endif
              enddo
          endif
          do itws=1,ntwsys(ngrnph(ngP))
              iPTS=iPTSgr(itws,ngP)
              iPTM=iSysMode(iPTS,ngrnph(ngP))
              iPTStw=iPTS-nslsys(ngrnph(ngP))
              !without children - nucleate
              if (iChildGrain(iPTS,ngP).eq.0.and.mod(iPTS,2).eq.1) then
                  !nucleation criteria
                  nuc=nucleation(TwFrSy(iPTStw,ngP))
                  if(nuc) then
                      !add new grain
                      ngrain=ngrain+1
                      iChildGrain(iPTS,ngP)=-ngrain
                      iChildGrain(iPTS+1,ngP)=ngrain
                      iParentGrain(ngrain)=ngP
                      iParentMode(ngrain)=iPTM
                      iParentSystem(ngrain)=iPTS
                      ngC=ngrain
                      wgtcg_all(ngC)=wgtcg_all(ngP)
                      nphngr(ngrnph(ngP))=nphngr(ngrnph(ngP))+1
                      ngrnph(ngC)=ngrnph(ngP)
                      call add_twin(ngP,ngC,iPTS,iPTM)
                  endif
              endif
              !with children - evolve weight
              if (iChildGrain(iPTS,ngP).ne.0.and.mod(iPTS,2).eq.1) then
                  ngC=abs(iChildGrain(iPTS,ngP))
                  wgtd(ngC)=0.0
                  !update weight through wgtd
                  wgtd(ngC)=wgtcg_all(ngp)*TwFrSy(iPTStw,ngP)-wgt(ngC)
                  !check for detwinning
                  if (wgt(ngC)+wgtd(ngC).le.0.0) then
                      wgtd(ngC)=-wgt(ngC)
                  else
                      wgtd(ngP)=wgtd(ngP)-wgtd(ngC)
                  endif
                  !update weight directly
c                  wgt(ngC)=wgtcg_all(ngp)*TwFrSy(iPTStw,ngP)
c                  wgt(ngP)=wgt(ngP)-wgt(ngC)
              endif
          enddo
        endif  
      enddo
      
      end subroutine
      
      function nucleation(TwFr)
      real, intent(in) :: TwFr
      logical nucleation
      
      TWINTHRES=0.01
      if (TwFr.ge.TWINTHRES) nucleation=.true.
      if (TwFr.lt.TWINTHRES) nucleation=.false.
      
      end function
      
      subroutine add_twin(ngP,ngC,iPTS,iPTM)
      
      use const
      use mphase_props
      use grain_props
      use grain_state
      use hard_law1_v
      use miscellaneous_sub, only : euler
      
      integer, intent(in) :: ngP,ngC,iPTS,iPTM
      
      dimension aux33(3,3),bur(3),aux(3,3)
      
      ! orientation and weight
      ! variant selection
      !calculates orientation of the child               
      do i=1,3
          bur(i)=bcc(i,iPTS,ngrnph(ngP))
          do j=1,3
              aux33(i,j)=r(j,i,ngP)
          enddo
      enddo
      call twinor(bur,aux33) !calculates transformation matrix from twinned crystal to sample
      do i=1,3
          do j=1,3
              r(i,j,ngC)=aux33(j,i)
              aux(i,j)=aux33(j,i)
          enddo
      enddo
      call euler(1,phi(ngC),the(ngC),ome(ngC),aux) !MZ calculate euler angles of twin orientation
      
      ! crystal properties
      call cr_to_sa(ngC,ngC,0)
      call cr_to_sa(ngC,ngC,1) 
      
      ! hard law state variables
      tau(:,ngC)=tau(:,ngP)
      rho_rev(:,ngC)=rho_rev(:,ngP)
      rho_forw(:,ngC)=rho_forw(:,ngP)
      rho_tot(:,ngC)=rho_tot(:,ngP)
      rho_tot_max(:,ngC)=rho_tot_max(:,ngP)
      rho_act(:,ngC)=rho_act(:,ngP)
      tau_act(:,ngC)=tau_act(:,ngP)
      rho_deb(ngC)=rho_deb(ngP)
      iact_sys(:,ngC)=iact_sys(:,ngP)
      TwFrSy(iPTS-nslsys(ngrnph(ngp)),ngp)=TwinFrac(iPTM,ngrnph(ngC))
      !correction for PTSdamp
      do ns=nslsys(ngrnph(ngP))+1,nsys(ngrnph(ngP))
          tau(ns,ngC)=tau(ns,ngP)/PTSDAMP(ns,ngp)          
      enddo
      
      !stress and strain in grains
c      stcs(:,ngC)=stcs(:,ngP)
c      etcs(:,ngC)=etcs(:,ngP)
      !finite initial fraction
      call InitializeChild(ngC) !corrects tau 
      
      end subroutine
      
      SUBROUTINE TWINOR(BUR,A)
      
      use miscellaneous_sub, only : euler
C
C *****************************************************************************
C     SUBROUTINE TWINOR      --->      VERSION OF 30/June/2005 BJCL
C
C       Changed the call to EULER so that the angles are given in degrees
C
C     GIVES THE TRANSFORMATION MATRIX 'A' (FROM CRYSTAL TO SAMPLE) OF THE
C     TWINNED RELATED RELATED CRYSTAL.
C     GIVEN THE BURGERS VECTOR OF THE TWIN SYSTEM (IN CRYSTAL AXES) MAKES
C     A ROTATION OF 180 DEG. AROUND IT TO DEFINE THE TWINNED CRYSTAL.
C *****************************************************************************
C
      DIMENSION BUR(3),HPI(3,3),HTW(3,3),A(3,3),AUX(3,3),ATW(3,3)
C
      DATA HPI/-1.,0.,0.,0.,-1.,0.,0.,0.,1./
      PI=4.0*ATAN(1.0)
C
      ANG1=ATAN2(BUR(2),BUR(1))+PI/2.0
      ANG2=SQRT(BUR(1)**2+BUR(2)**2)
      ANG2=ATAN2(ANG2,BUR(3))
      CALL EULER (2,ANG1*180.0/PI,ANG2*180.0/PI,0.,AUX)

      DO 10 I=1,3
      DO 10 J=1,3
      HTW(I,J)=0.
      DO 10 K1=1,3
      DO 10 K2=1,3
   10 HTW(I,J)=HTW(I,J)+AUX(K1,I)*HPI(K1,K2)*AUX(K2,J)

      DO 20 I=1,3
      DO 20 J=1,3
      ATW(I,J)=0.
      DO 20 K=1,3
   20 ATW(I,J)=ATW(I,J)+A(I,K)*HTW(K,J)
c
      DO I=1,3
      DO J=1,3
        A(I,J)=ATW(I,J)
      ENDDO
      ENDDO
c
      RETURN
      END SUBROUTINE
      
      SUBROUTINE InitializeChild(iChild)
c     Calculate the stress and elastic strain state in the twin child
c     using the following rules for continuity across the boundary
c     between twin and parent (axis 3 is the twin plane normal and
c     axis 1 is along the burgers vector):
c     Sig33_T=Sig33_P, Sig13_T=Sig13_P and Sig23_T=Sig23_P
c     Eps11_ET=Eps11_EP, Eps22_ET=Eps22_EP and Eps12_ET=Eps12_EP
c
c     The stress and elastic strain tensors of both parent and child grains
c     and the stiffness tensor of the child grain must be rotated to the
c     coordinate system of the twin system before the boundary conditions
c     can be applied, and then rotated back to the sample system again
c
c     The relsoved shear stresses (RSS) for all systems are calculated,
c     and if any are larger than the critical resolved shear stresses (CRSS),
c     the value of CRSS is set to the calculated RSS for that system
c      
      use mvoigt
      use miscellaneous_sub, only : ludcmpc,lubksbc
      use mphase_props
      use grain_props
      use grain_state
      use grain_rate
      
      integer, intent(in) :: iChild
      
      DIMENSION etplParent(6),aux6(6),BackStress(6),TwinStrain(6)
c

      imo=iParentMode(iChild)
      igr=iParentGrain(iChild)
      ist=iParentSystem(iChild)
      iph=ngrnph(igr)

      if(TwinFrac(imo,iph).gt.0.0) then
          Gamma0=TwinFrac(imo,iph)*stw(imo,iph)
        do i=1,6
          TwinStrain(i)=mcs(i,ist,igr)*Gamma0
        enddo
        do i=1,6
          BackStress(i)=0.0
          do j=1,6
            BackStress(i)=BackStress(i)-
     #        ccs2(i,j,iChild)*TwinStrain(j)*profac(j)
          enddo
c            stcs(i,iChild)=BackStress(i)
c  Merkel, 05/2010: stress in child is stress in parent + backstress
            stcs(i,iChild) = stcs(i,igr) + BackStress(i)
        enddo
      else
        call CalcChildStress(iChild)
      endif
c
c       write(*,*) 'Child Stress in sample'
c      do i=1,6
c       write(*,*) stcs(i,iChild)
c      enddo
c
c     Calculate the plastic strain of the parent
      iParent=iParentGrain(iChild)
c Merkel 05/2010, if large strains, we already stored the elastic strain
      if (kSM.eq.1) then
        do i=1,6
          etplParent(i)=etcs(i,iParent)-etelcs(i,iParent)
        enddo
c if not, calculate it using Hooke
      else
        do i=1,6
          aux6(i)=0.0
          do j=1,6
            aux6(i)=aux6(i)+scs2(i,j,iParent)*stcs(j,iParent)*profac(j)
          enddo
        enddo
        do i=1,6
          etplParent(i)=etcs(i,iParent)-aux6(i)
        enddo
      endif
c
c     Calculate the strain in the twin (elastic + parent plastic)
      do i=1,6
        aux6(i)=0.0
        do j=1,6
          aux6(i)=aux6(i)+scs2(i,j,iChild)*stcs(j,iChild)*profac(j)
        enddo
      enddo
      do i=1,6
        etcs(i,iChild)=aux6(i)+etplParent(i)
      enddo
c Merkel 05/2010, if large strains, we will store the elastic strain
      if (kSM.eq.1) then
        do i=1,6
          etelcs(i,iChild)=aux6(i)
        enddo
       ! We also need to set the hydrostatic strains
       ! in the newly formed twin. We copy the hydrostatic strain in
       ! the parent grain
       do i=1,6
         etelhycs(i,iChild) = etelhycs(i,iParent)
       enddo
      endif

c     Initialize the state of the new grain
      wgt(iChild)=0.0
      wgtd(iChild)=0.0
c       gamtot(iChild)=gamtot(iParent)
      gamtot(iChild)=0.0
      ns1=0
      do mo=1,nmodes(iph)
        do isys=1,nsm(mo,iph)
          ns1=ns1+1
c-------------------------------------------------------------------------------
c       Initialize the internal state of the newly form twin
c-------------------------------------------------------------------------------
          if(kCL.NE.2.and.kCL.NE.3) then !ADDEDMZ
c            tau(ns1,iChild)=tau(ns1,iParent)
c            tau(ns1,iChild)=tau0(ns1)
          else
            if (iTwinLaw.eq.1) then 
                tau(ns1,iChild)=tau(ns1,iParent) 
c                rho_for(ns1,iChild)=rho_for(ns1,iParent) 
c                rho_deb(iChild)=rho_deb(iParent) 
            endif
          endif
        enddo
      enddo
c       Ensuring that the yield surface of the Child is not exceeded by the actual stress state
      do ns1=1,nsys(iph)
        rss=0.0
        do i=1,6
          rss=rss+mcs(i,ns1,iChild)*stcs(i,iChild)*profac(i)
        enddo
        if(rss.gt.tau(ns1,iChild)) tau(ns1,iChild)=rss
      enddo
      do i=1,6
        etrcs(i,iChild)=0.0
        strcs(i,iChild)=0.0
      stcsref(i,iChild)=0.0
c        stcsref(i,iChild)=stcsref(i,iParent)
      enddo
c Merkel, 05/2010, for large strains, reference stresses have to be set
c because strains are calculated incrementally
      if (kSM.eq.1) then
        do i=1,6
          stcsref(i,iChild)=stcs(i,iChild)
        enddo
      endif

      do ns1=1,nsys(iph)
        taud(ns1,iChild)=0.0
        gamd(ns1,iChild)=0.0
      enddo
c
      END SUBROUTINE
c _________________________________________________________________________
c
      SUBROUTINE CalcChildStress(iChild)
c     Calculate the stress and elastic strain state in the twin child
c     using the following rules for continuity across the boundary
c     between twin and parent (axis 3 is the twin plane normal and
c     axis 1 is along the burgers vector):
c     Sig33_T=Sig33_P, Sig13_T=Sig13_P and Sig23_T=Sig23_P
c     Eps11_ET=Eps11_EP, Eps22_ET=Eps22_EP and Eps12_ET=Eps12_EP
c
c     The stress and elastic strain tensors of both parent and child grains
c     and the stiffness tensor of the child grain must be rotated to the
c     coordinate system of the twin system before the boundary conditions
c     can be applied, and then rotated back to the sample system again
c
c     The relsoved shear stresses (RSS) for all systems are calculated,
c     and if any are larger than the critical resolved shear stresses (CRSS),
c     the value of CRSS is set to the calculated RSS for that system
c
      use mvoigt
      use miscellaneous_sub, only : ludcmpc,lubksbc
      use mphase_props
      use grain_props
      use grain_state
      use grain_rate
      
      integer, intent(in) :: iChild
      
      DIMENSION stParent(6),etelParent(6),stChild(6),etelChild(6)
      DIMENSION aChild(6,6),aux66(6,6),aux6(6),rTW2SA(3,3),rSA2TW(3,3)
      DIMENSION istTW(6),ietTW(6),aux3333(3,3,3,3),aux33(3,3),indx(6)
c     The stress and elastic strain twin boundary conditions as given above
      data istTW/0,0,1,1,1,0/ !MZ S_F
      data ietTW/1,1,0,0,0,1/ !MZ S_F
c
      iParent=iParentGrain(iChild)
      isys=iParentSystem(iChild)


c     Define rTW2SA which is the rotation from twin system coordinates to sample coordinates
c     The twin system burgers vector is 1 axis
      vlen=0.0
      do i=1,3
          vlen=vlen+bcs(i,isys,iParent)**2
      enddo
      vlen=sqrt(vlen)
      do i=1,3
          rTW2SA(i,1)=bcs(i,isys,iParent)/vlen
      enddo
c     The twin system plane normal is the 3 axis
      vlen=0.0
      do i=1,3
          vlen=vlen+ncs(i,isys,iParent)**2
      enddo
      vlen=sqrt(vlen)
      do i=1,3
          rTW2SA(i,3)=ncs(i,isys,iParent)/vlen
      enddo
c     The 2 axis is found as the cross product of the 3 and 1 axes
      rTW2SA(1,2)=rTW2SA(2,3)*rTW2SA(3,1)-rTW2SA(3,3)*rTW2SA(2,1)
      rTW2SA(2,2)=rTW2SA(3,3)*rTW2SA(1,1)-rTW2SA(1,3)*rTW2SA(3,1)
      rTW2SA(3,2)=rTW2SA(1,3)*rTW2SA(2,1)-rTW2SA(2,3)*rTW2SA(1,1)
c     Define rSA2TW as the inverse (transpose) of rTW2SA
      do i=1,3
        do j=1,3
          rSA2TW(i,j)=rTW2SA(j,i)
        enddo
      enddo
c
c     Rotate the Parent stress and elastic strain and the Child stiffness to the twin coordiante system
      do i=1,6
        aux6(i)=stcs(i,iParent)
        enddo
      call voigt(aux6,aux33,aux66,aux3333,1)
      do ij=1,6
        i=ijv(ij,1)
        j=ijv(ij,2)
        stParent(ij)=0.0
        do i1=1,3
          do j1=1,3
            stParent(ij)=stParent(ij)
     #      +rSA2TW(i,i1)*rSA2TW(j,j1)*aux33(i1,j1)
          enddo
        enddo
      enddo
c S. Merkel, large elastic strain, we store the elastic strain
      if (kSM.eq.1) then
        do i=1,6
          aux6(i)=etelcs(i, iParent)
        enddo
      else
        do i=1,6
          aux6(i)=0.0
          do j=1,6
            aux6(i)=aux6(i)+scs2(i,j,iParent)*
     #                stcs(j,iParent)*profac(j)
          enddo
        enddo
      endif
      call voigt(aux6,aux33,aux66,aux3333,1)
      do ij=1,6
        i=ijv(ij,1)
        j=ijv(ij,2)
        etelParent(ij)=0.0
        do i1=1,3
          do j1=1,3
            etelParent(ij)=etelParent(ij)
     #      +rSA2TW(i,i1)*rSA2TW(j,j1)*aux33(i1,j1)
            enddo
          enddo
        enddo
        do i=1,6
          do j=1,6
            aux66(i,j)=acs2(i,j,iChild)
          enddo
        enddo
      call voigt(aux6,aux33,aux66,aux3333,3)
      do ij=1,6
        i=ijv(ij,1)
        j=ijv(ij,2)
        do kl=1,6
          k=ijv(kl,1)
          l=ijv(kl,2)
          aChild(ij,kl)=0.0
          do i1=1,3
            do j1=1,3
              do k1=1,3
                do l1=1,3
                   aChild(ij,kl)=aChild(ij,kl)
     #             +rSA2TW(i,i1)*rSA2TW(j,j1)*rSA2TW(k,k1)
     #             *rSA2TW(l,l1)*aux3333(i1,j1,k1,l1)
                enddo
              enddo
              enddo
            enddo
        enddo
      enddo
c
c     Calculate the stress and elastic strain state in the child grain
      do i=1,6
       aux6(i)=-1.0*istTW(i)*stParent(i)
        do j=1,6
          aux6(i)=aux6(i)+aChild(i,j)*ietTW(j)*etelParent(j)*profac(j)
          aux66(i,j)=ietTW(j)*(i/j)*(j/i)-
     #             istTW(j)*aChild(i,j)*profac(j)
        enddo
      enddo
      call ludcmpc(aux66,6,6,indx,d)
      call lubksbc(aux66,6,6,indx,aux6)
      do i=1,6
        etelChild(i)=ietTW(i)*etelParent(i)+istTW(i)*aux6(i)
        stChild(i)=istTW(i)*stParent(i)+ietTW(i)*aux6(i)
      enddo
c
c     Rotate the stresses and elastic strains back into sample coordinate system
      call voigt(stChild,aux33,aux66,aux3333,1)
      do ij=1,6
        i=ijv(ij,1)
        j=ijv(ij,2)
        stcs(ij,iChild)=0.0
        do i1=1,3
          do j1=1,3
            stcs(ij,iChild)=stcs(ij,iChild)
     #      +rTW2SA(i,i1)*rTW2SA(j,j1)*aux33(i1,j1)
          enddo
        enddo
      enddo
c Merkel, for large strains we also need to store the elastic strain in the child
      if (kSM.eq.1) then
        call voigt(etelChild,aux33,aux66,aux3333,1)
        do ij=1,6
          i=ijv(ij,1)
          j=ijv(ij,2)
          etelcs(ij,iChild)=0.0
          do i1=1,3
            do j1=1,3
              etelcs(ij,iChild)=etelcs(ij,iChild)
     #        +rTW2SA(i,i1)*rTW2SA(j,j1)*aux33(i1,j1)
            enddo
          enddo
        enddo
        ! Copy hydrostatic strain of the parent into the child
        do i=1,6
          etelhycs(i,iChild) = etelhycs(i,iParent)
        enddo
      endif

      END SUBROUTINE
      
      subroutine corr_p_st(stcs)
      
      use flags, only : iTwPh
      use mphase_props
      use sample_props_v, only : ngrain
      real,intent(inout) :: stcs(6,NGR)
      
      dimension aux6(6)
      !If new grains were added in this step the stress in the parent grain
      !is corrected with the back-stress   
      do ng=1,ngrain
        if(iTwPh(ngrnph(ng)).eq.1) then
          do i=1,6
            aux6(i)=0.0
          enddo
          dumtemp=0.0
          do ist=1,nsys(ngrnph(ng))
            j=iChildGrain(ist,ng)
            if(j.lt.0) then
              j=-j
              iChildGrain(ist,ng)=j
              imo=iParentMode(j)
              do i=1,6
                aux6(i)=aux6(i)+stcs(i,j)*TwinFrac(imo,ngrnph(ng))
              enddo
              dumtemp=dumtemp+TwinFrac(imo,ngrnph(ng))
            endif
          enddo
          if(dumtemp.ne.0.0) then
            do i=1,6
              stcs(i,ng)=(stcs(i,ng)+aux6(i))/(1-dumtemp)
            enddo
          endif
        endif  
      enddo
      
      end subroutine
      
      SUBROUTINE TWIN_SHAPE(KKK,KGX,IPTS,AXISGR) !KGX is the number of parent and KKK is the number of twin
      !Defines ellipsoid of grain as in VPSC.       
      !AXISGR is assigned and updated for twins.

      use const
      use mphase_props_v, only : twthres,nslsys,ngrnph,ncc,bcc
      use miscellaneous_sub, only : det
      use grain_props_v, only : r
c      use grain_state_v, only : AXISGR
      
      integer,intent(in) :: KKK,KGX,IPTS
      real,intent(inout) :: AXISGR(0:3,3,NGR)
      
      DIMENSION ANG(2),DUMMY(3,3),AA(3,3)
      REAL MAT_PROJ(3,3), PROJECTION(3), NORM_PROJ
      REAL RATIO
      DIMENSION W(3),BX(3,3),B(3,3),BT(3,3),TWIN_ORIENT(3,3,NGR)
           
      !Creation of ellipsoid representing the twin
      !(1) Definition of orientation of elliopsoid axes. 
      !    Using b and n vectors and their cross product, we get axes of ellipsoid.
      !    These are stored as [b , n x b , n] in TWIN_ORIENT (column vectors).
      !    n and b are defined in data_crystal (ncc and bcc).
      IPTM = iParentMode(KKK)
      IF(newGrSh(KKK).EQ.0) THEN          
          do i = 1,3
              TWIN_ORIENT(i,1,KKK) = bcc(i,IPTS,ngrnph(KKK))
              TWIN_ORIENT(i,3,KKK) = ncc(i,IPTS,ngrnph(KKK))
          enddo
          TWIN_ORIENT(1,2,KKK)=TWIN_ORIENT(2,3,KKK)*TWIN_ORIENT(3,1,KKK)
     #                    -TWIN_ORIENT(3,3,KKK)*TWIN_ORIENT(2,1,KKK)
          TWIN_ORIENT(2,2,KKK)=TWIN_ORIENT(3,3,KKK)*TWIN_ORIENT(1,1,KKK)
     #                    -TWIN_ORIENT(1,3,KKK)*TWIN_ORIENT(3,1,KKK)
          TWIN_ORIENT(3,2,KKK)=TWIN_ORIENT(1,3,KKK)*TWIN_ORIENT(2,1,KKK)
     #                    -TWIN_ORIENT(2,3,KKK)*TWIN_ORIENT(1,1,KKK)
           
      !** DEFINE TWIN ELLIPSOID ORIENTATION IN SAMPLE FRAME
      ! Takes the TWIN_ORIENT (which has in it b and n) and transforms it 
      ! to sample frame.           
           DO I=1,3
           DO J=1,3
              DUMMY(I,J)=0
              DO K=1,3
                 DUMMY(I,J)=DUMMY(I,J)+r(K,I,KGX)*TWIN_ORIENT(K,J,KKK) !transpose of r taken
              ENDDO
           ENDDO
           ENDDO

           DO I=1,3
           DO J=1,3
              TWIN_ORIENT(I,J,KKK)=DUMMY(I,J)
              AXISGR(I,J,KKK)=DUMMY(I,J)
           ENDDO
           ENDDO          
           
      !(2) Calculation of length of twin ellipsoid axes. 
      !First, the orientation of twin axes in ellipsoid of parent in spherical coordinates
      !is found (ANG(1) and ANG(2)). Then vector from origin of ellipsoid axes to surface 
      !of parent ellipsoid in direction of twin axes is found. Then it's length is 
      !calculated and used as starting ellipsoid length of twins. Afterwards third axes of
      !twin ellipsoid is: volume fraction of predominant twin system x threshold value
           
      !**   DEFINE AXES LENGTHS FOR ELLIPSOID REPRESENTING THE TWIN
      !**   WE NEED TO CALCULATE THE LENGTHS OF EACH AXIS IN THE ELLIPSOID OF THE GRAIN
          DO I=1,3
      !** WE FIRST NEED TO CALCULATE THE POLAR (ANG(2)) AND AZIMUTHAL (ANG(1)) ANGLES              
             DOT_PRODUCT=0
             DO J=1,3
                DOT_PRODUCT=DOT_PRODUCT+AXISGR(I,J,KKK)*AXISGR(3,J,KGX)
             ENDDO
             if(dot_product.gt.1.0) dot_product=1.0
             if(dot_product.lt.-1.0) dot_product=-1.0
             ANG(2)=ACOS(DOT_PRODUCT)
             IF(ANG(2).LT.0) ANG(2)=-ANG(2)
             IF(ANG(2).GT.3.1415926535898)
     #              ANG(2)=2*3.1415926535898-ANG(2)
             DO L=1,3
             DO M=1,3
                IF(L.EQ.M) THEN
                   MAT_PROJ(L,M)=1
                ELSE
                   MAT_PROJ(L,M)=0
                ENDIF
                MAT_PROJ(L,M)=MAT_PROJ(L,M)-AXISGR(3,L,KGX)
     #                                     *AXISGR(3,M,KGX)
             ENDDO
             ENDDO
             DO L=1,3
                PROJECTION(L)=0
                DO M=1,3
                   PROJECTION(L)=PROJECTION(L)
     #                +MAT_PROJ(L,M)*AXISGR(I,M,KKK)
                ENDDO
             ENDDO
             NORM_PROJ=0
             DO L=1,3
                NORM_PROJ=NORM_PROJ+PROJECTION(L)*PROJECTION(L)
             ENDDO
             NORM_PROJ=SQRT(NORM_PROJ)
             DO L=1,3
                PROJECTION(L)=PROJECTION(L)/NORM_PROJ
             ENDDO
             DOT_PRODUCT=0
             DO J=1,3
                DOT_PRODUCT=DOT_PRODUCT+PROJECTION(J)*AXISGR(1,J,KGX)
             ENDDO
             if(dot_product.gt.1.0) dot_product=1.0
             if(dot_product.lt.-1.0) dot_product=-1.0
             ANG(1)=ACOS(DOT_PRODUCT)
             X1=AXISGR(0,2,KGX)*SIN(ANG(2))*COS(ANG(1))
             X2=AXISGR(0,1,KGX)*SIN(ANG(2))*SIN(ANG(1))
             X3=AXISGR(0,3,KGX)*COS(ANG(2))
             AXISGR(0,I,KKK)=SQRT(X1*X1+X2*X2+X3*X3)
          ENDDO
          AXIS3(KKK)=AXISGR(0,3,KKK) ! NEEDS TO BE KEPT IN MEMORY FOR UPDATING TWIN THICKNESS
          twthres1=twthres(1,IPTM,ngrnph(KKK))
          AXISGR(0,3,KKK)=TWFRSY(IPTS-nslsys(ngrnph(KKK)),KGX)*TWTHRES1
     #                   *AXISGR(0,3,KKK)
          RATIO=AXISGR(0,2,KKK)/AXISGR(0,3,KKK)
          IF(RATIO.GT.20.0) AXISGR(0,3,KKK)=AXISGR(0,2,KKK)/20.0    !LIMIT THE ELLIPSOID RATIO TO BE GREATER THAN 20
                                                                !AVOID PROBLEM FOR ESHELBY TENSOR CALCULATION

          !THE LONGEST AXIS SHOULD BE THE 2ND ONE
          IF(AXISGR(0,1,KKK).GT.AXISGR(0,2,KKK)) THEN
            EXCHANGE=AXISGR(0,1,KKK)
            AXISGR(0,1,KKK)=AXISGR(0,2,KKK)
            AXISGR(0,2,KKK)=EXCHANGE
            DO I=1,3
               EXCHANGE=AXISGR(1,I,KKK)
               AXISGR(1,I,KKK)=AXISGR(2,I,KKK)
               AXISGR(2,I,KKK)=EXCHANGE
            ENDDO
          ENDIF

          !   **    MAKE THE SYSTEMS RIGHT HANDED
          DO I=1,3
             DO J=1,3
                AA(I,J)=AXISGR(I,J,KKK)
             ENDDO
          ENDDO

          IF(DET(AA).LE.0.) THEN
            DO I=1,3
               AXISGR(2,I,KKK)=-AXISGR(2,I,KKK)
            ENDDO
          ENDIF           
          newGrSh(KKK)=1.                 
      ELSE
c          if (KKK.eq.ngParent+1) write(72,'(31f15.4)') axisgr(0,3,KKK)
           THRES1=TWTHRES(1,IPTM,ngrnph(KKK))
           AXISGR(0,3,KKK)=TWFRSY(IPTS-nslsys(ngrnph(KKK)),KGX)*THRES1
     #                  *AXIS3(KKK)
           RATIO=AXISGR(0,2,KKK)/AXISGR(0,3,KKK)
           IF(RATIO.GT.20.0) AXISGR(0,3,KKK)=AXISGR(0,2,KKK)/20.0
           THRES2=TWTHRES(2,IPTM,ngrnph(KGX))
c           if(kkk.le.(2*ngr(1))) then !update shape of parent only if we have a 1st generation twin
           !MZ (comment) - following condition won't be satisfied at any point in case of THRES2 = 1, so following
           !code won't be executed. No changes were made in copying from VPSC and it's not clear if it will
           !work properly.           
             IF(TWFRSY(IPTS-nslsys(ngrnph(KGX)),KGX).GE.THRES2) THEN
               if(newGrSh(KGX).EQ.0) then
                  DO I=1,3
                  DO J=1,3
                     AXISGR(I,J,KGX)=AXISGR(I,J,KKK)  !GIVE TO PARENT SAME ELLIPSOID ORIENTATION AS TWIN
                  ENDDO
                  ENDDO
                  DO I=1,2
                    AXISGR(0,I,KGX)=AXISGR(0,I,KKK) !2 OF ELLIPSOID AXES HAVE SAME LENGTH FOR PARENT AND TWIN
                  ENDDO
                  newGrSh(kgx)=1
               endif ! if(newGrSh(KGX).EQ.0)
               AXISGR(0,3,KGX)=(1-TWFRSY(IPTS-nslsys(ngrnph(KKK)),KGX))
     #                  *THRES1*AXIS3(KKK)
               RATIO=AXISGR(0,2,KGX)/AXISGR(0,3,KGX)
               IF(RATIO.GT.20.0) AXISGR(0,3,KGX)=AXISGR(0,2,KGX)/20.0
             ELSE
                 newGrSh(KGX)=0
             ENDIF !  IF(TWFRSY(IPTS,KGX).GE.THRES)
c           endif !  if(kkk.le.(2*ngr(1)))    
      ENDIF             
      
      
      RETURN
      END SUBROUTINE
        
      end module twinning        
c      
c***********************************************************************
c *** back_stress ******************************************************
c***********************************************************************
      module back_stress
      
      use back_stress_v
      
      private 
      
      public bs_read,bs_hd,bs_update_statv,tau_bcst,set_bs,get_bs
     #    ,write_statv_bs,read_statv_bs,state_bs,ph_bcst
      
      contains
      
      !bs_read: reads parameters for backstress law
      !bs_update_statv: updates state variables of backstress law
      !bs_hd: calculates hardening matrix for bs
      !set_bs: sets statv object
      !get_bs: gets values from statv object
      !write_statv_bs: writes statv object to statv array
      !read_statv_bs: reads statv object from statv array
      
      subroutine set_bs(ng,bs)
      
      use hard_law1_v, only : gam_acc
      
      integer, intent(in) :: ng
      type(state_bs),intent(out) :: bs
      
      bs%tau_bcst=tau_bcst(:,ng,:)
      bs%gam_acc=gam_acc(:,ng)
      
      end subroutine set_bs
      
      subroutine get_bs(ng,bs)
      
      use hard_law1_v, only : gam_acc
      
      integer, intent(in) :: ng
      type(state_bs),intent(in) :: bs
      
      tau_bcst(:,ng,:)=bs%tau_bcst
      gam_acc(:,ng)=bs%gam_acc
      
      end subroutine get_bs
      
      subroutine write_statv_bs(bs,ng,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_bs),intent(in) :: bs
      integer,intent(in) :: ng,NSTATV
      integer,intent(inout) :: ns
      real,intent(out) :: STATEV(NSTATV)
      
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=bs%tau_bcst(1:nsys(ngrnph(ng)),1)
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=bs%tau_bcst(1:nsys(ngrnph(ng)),2)
      ns=ns+nsys(ngrnph(ng))
      STATEV(ns:ns+nsys(ngrnph(ng))-1)=bs%gam_acc(1:nsys(ngrnph(ng)))
      ns=ns+nsys(ngrnph(ng))
      
      end subroutine write_statv_bs
      
      subroutine read_statv_bs(bs,ng,iph,ns,NSTATV,STATEV)
      
      use mphase_props
      
      type(state_bs),intent(out) :: bs
      integer,intent(in) :: ng,NSTATV,iph
      real,intent(in) :: STATEV(NSTATV)
      integer,intent(inout) :: ns
      
      bs%tau_bcst(1:nsys(iph),1)=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      bs%tau_bcst(1:nsys(iph),2)=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      bs%gam_acc(1:nsys(iph))=STATEV(ns:ns+nsys(iph)-1)
      ns=ns+nsys(iph)
      
      end subroutine read_statv_bs
      
      subroutine bs_read(iph)
      
      !initialize (read) parameters
      
      use flags
      use mphase_props
      
      ilatBS=2
      ibs=1
      if (ibs.eq.1) then
        read(1,*)
        do imo=1,nslmod(iph)
            read(1,*) tau_sat(imo,iph),ni(imo,iph) !et_sat,
            read(1,*) gam_b(imo,iph),fact_1(imo,iph)
        enddo
        read(1,*) iphBsCtrl(iph) !0 independent, 1 controlling, -1 controlled
        read(1,*) iDiagBs,alpha_reg_bs
        
        ! find latent bs matrix    
        if (ilatBS.eq.2) then
            do is1=1,nslsys(iph)
            do is2=1,nslsys(iph)    
                aux=0.0
                do i=1,3
                    do j=1,3
                       aux = aux + 2.0*mc2(i,j,is1,iph)*mc2(i,j,is2,iph)
                    enddo
                enddo
c                if (aux.gt.0.0) then
                aL_bs(is1,is2,iph)=aux
c                endif
            enddo
            enddo
           
            do is1=1,nslsys(iph)
            do is2=1,nslsys(iph)
c                aL_bs(is1,is2,iph)=0.0
            enddo
            enddo
            
        endif
      elseif(ibs.eq.2)then
        read(1,*) 
        read(1,*) vol_frac 
      elseif(ibs.eq.3)then        
        read(1,*)
        !read only 1st mode and use for whole grain
        read(1,*) tau_sat(1,iph),ni(1,iph) !et_sat,
        read(1,*) gam_b(1,iph),fact_1(1,iph)
        read(1,*) iphBsCtrl(iph) !0 independent, 1 controlling, -1 controlled
        
        do is1=1,nslsys(iph)
        do is2=1,nslsys(iph)    
            aux=0.0
            do i=1,3
                do j=1,3
                   aux = aux + 2.0*mc2(i,j,is1,iph)*mc2(i,j,is2,iph)
                enddo
            enddo
c            if (aux.gt.0.0) then
            aL_bs(is1,is2,iph)=aux
c            endif
        enddo
        enddo
      endif
      
      end subroutine 
      
      subroutine bs_hd(ng,hd_bs)  
      
      !Calculation of hd_bs

      use const
      use mvoigt
      use mphase_props
      use grain_props_v
      use grain_rate, only : gamd
      use grain_state, only : nact, iact
      use hard_law1_v, only : itemp_act_sys,iOpposite,iDirSys
     #    ,iact_sys
      use flags, only : ilatBS
      
      integer, intent(in) :: ng
      real, intent(out) :: hd_bs(NSLS,NSLS,2)
      dimension CA(6,6)
      
      iph=ngrnph(ng)
      
      if (ibs.eq.1) then
        !find reversing and active systems   
        !Make hd_bs be zero matrix        
        do ns1=1,nsys(iph)
            do ns2=1,nsys(iph)
                hd_bs(ns1,ns2,:)=0.0
            enddo
        enddo    
      
        !update iact_sys array
        do is=1,nact(ng)
            n1=iact(is,ng)
            if (iTwinSys(n1).EQ.0) then
                if (mod(n1,2).EQ.1) then
                    iact_sys(iDirSys(n1,iph),ng)=1
                else
                    iact_sys(iDirSys(n1,iph),ng)=2
                endif
            else
                iact_sys(iDirSys(n1,iph),ng)=1
            endif
        enddo
        
        !find iDifference (not equal 0 if there is reversal of slip system)
        do is1=1,nact(ng)
          n1=iact(is1,ng)
          id=iDirSys(n1,iph)
          
          if (itemp_act_sys(id,ng).NE.0) then !system was active at some point
              iDifference=itemp_act_sys(id,ng)-iact_sys(id,ng)
          endif
             
        !define current back stress (tau_0b) and acc. shear strain (gam_0b) at point of reversal  
          if (iDifference.ne.0) then
          
          endif
        enddo                     
        !-------------------------------------define elements of hardening matrix based on evolution law
        ! there are two cases
        ! if back stress > 0 then evolution is normal (shear loop evolution law)
        ! if back stress < 0 then evolution is rapid decay since this means that reversal
        !    has happened on sys and large back stress in opposite direction needs to be 
        !    destroyed (fact times larger)
        ! in both cases opposite system has correct evolution by multiplication with fact1
        
        ! (1) old evolution
        do ns=1,nact(ng)
            n1=iact(ns,ng)
            imo=iSysMode(n1,iph)
            if (tau_sat(imo,iph).ne.0.0.and.iTwinSys(n1).eq.0) then
                if (tau_bcst(n1,ng,1).ge.0.0) then
                    if (tau_bcst(n1,ng,1).ge.tau_sat(imo,iph)) then !in case backstress overshoots gam_star is NaN
                        gam_star=0.0
                        hd_bs(n1,n1,1)=0.0
                    else
                        gam_star=-1.0/ni(imo,iph)*log( 1.0-
     #                         tau_bcst(n1,ng,1)/tau_sat(imo,iph) )
                        hd_bs(n1,n1,1)=bs_org_der(tau_sat(imo,iph)
     #                      ,ni(imo,iph),gam_star)
                        fact1=fact_1(imo,iph)
                    endif
                endif      
                
                if (tau_bcst(n1,ng,1).lt.0.0) then
                   gam_star_1 = -gam_b(imo,iph)*log(  
     #                  ( tau_sat(imo,iph) - tau_bcst(n1,ng,1)  ) / 
     #                  ( (fact_1(imo,iph)+1.0)*tau_sat(imo,iph) ) )
                   hd_bs(n1,n1,1) = bs_new_der(tau_sat(imo,iph)
     #                   ,gam_b(imo,iph),fact_1(imo,iph),gam_star_1)
                   if (tau_bcst(n1,ng,1).lt.0.0) then
                        fact1=1.0/fact_1(imo,iph)
                   else
                        fact1=fact_1(imo,iph)
                   endif
                endif    
                
                if (hd_bs(n1,n1,1).le.1e-3) then
                    hd_bs(n1,n1,1)=0.0
                endif
                hd_bs(iOpposite(n1),n1,1)=-fact1*hd_bs(n1,n1,1)
            endif
        enddo         

        !evolve tau_bcst with sum of shears i.e. sum hd_bs matrix along diagonal
c        aux=0.0
c        aux1=0.0
c        do ns=1,nact(ng)
c            n1=iact(ns,ng)
c            aux=aux+hd_bs(n1,n1,1)/real(nact(ng))
c            aux1=aux1+hd_bs(iOpposite(n1),n1,1)/real(nact(ng))
c        enddo 
c        do ns=1,nact(ng)
c            n1=iact(ns,ng)
c            hd_bs(n1,n1,1)=aux
c            hd_bs(iOpposite(n1),n1,1)=aux1
c        enddo
        
        ! include latent bs matrix
        if (ilatBS.eq.2) then
          do ns=1,nact(ng)
              n1=iact(ns,ng)
              hd_bs(n1,n1,2)=hd_bs(n1,n1,1)
              hd_bs(iOpposite(n1),n1,2)=hd_bs(iOpposite(n1),n1,1)
          enddo
          do ns1=1,nslsys(iph)
              do ns2=1,nact(ng)
                  n1=iact(ns2,ng)
                  if (iTwinSys(n1).eq.0) then
                      if (ns1.ne.n1.and.ns1.ne.iOpposite(n1)) then
                          if (tau_bcst(n1,ng,1).ge.0.0) then
                              hd_bs(ns1,n1,2) = hd_bs(ns1,n1,2) + 
     #                              aL_bs(ns1,n1,iph)*hd_bs(n1,n1,1)
                          else
                              hd_bs(ns1,n1,2) = hd_bs(ns1,n1,2) + 
     #                              aL_bs(ns1,iOpposite(n1),iph)
     #                              *hd_bs(iOpposite(n1),n1,1)
                          endif
                      endif
                  endif
              enddo
          enddo      
        endif
      
      
c     MZ write hd to file
        if (ng.eq.1) then
c        do i=1,48
c            write(71,'(48f9.1)') (hd_bs(i,j),j=1,48)
c        enddo          
c        write(71,'(24f9.1)') hd_bs(3,3:4)
c        write(71,'(24f9.1)') hd_bs(4,3:4)
c        write(71,*) '***********************************************'
        endif            
          
      elseif(ibs.eq.2) then
        call rod_esh(r(:,:,ng),CA)
        hd_bs = 0.0
        do n1=1,nsys(iph)
            !n1=iact(ns1,ng)
            do ns2=1,nact(ng)
                n2=iact(ns2,ng)
                do i=1,6
                    do j=1,6
                        hd_bs(n1,n2,1) = hd_bs(n1,n2,1) - mcs(i,n1,ng)*
     #                     mcs(j,n2,ng)*vol_frac/(1.0-vol_frac)*
     #                     CA(i,j)*profac(i)*profac(j)!minus sign since the backstress is subtracted from stress
                    enddo
                enddo
                !if (hd_bs(n1,n2).le.0.0) hd_bs(n1,n2)=1.0
            enddo
        enddo
        
        
      elseif (ibs.eq.3) then        
        !total bs per mode in function of sum of shear strains
        !Make hd_bs be zero matrix        
        do ns1=1,nsys(iph)
            do ns2=1,nsys(iph)
                hd_bs(ns1,ns2,:)=0.0
            enddo
        enddo    
        !find Gamma* from current bacstress
        if (tau_bcst_tot(ng).ge.0.0) then
          if (tau_bcst_tot(ng).ge.tau_sat(1,iph)) then !in case backstress overshoots gam_star is NaN
              gam_star=0.0
              hd_bs_tot=0.0
          else
              gam_star=-1.0/ni(1,iph)*log( 1.0-
     #               tau_bcst_tot(ng)/tau_sat(1,iph) )
              hd_bs_tot=bs_org_der(tau_sat(1,iph)
     #            ,ni(1,iph),gam_star)
              fact1=fact_1(1,iph)
          endif
        endif

        if (tau_bcst_tot(ng).lt.0.0) then
           gam_star_1 = -gam_b(1,iph)*log(  
     #          ( tau_sat(1,iph) - tau_bcst_tot(ng)  ) / 
     #          ( (fact_1(1,iph)+1.0)*tau_sat(1,iph) ) )
           hd_bs_tot = bs_new_der(tau_sat(1,iph)
     #           ,gam_b(1,iph),fact_1(1,iph),gam_star_1)
           if (tau_bcst_tot(ng).lt.0.0) then
                fact1=1.0/fact_1(1,iph)
           else
                fact1=fact_1(1,iph)
           endif
         endif  
         
         if (hd_bs_tot.le.1e-3) then
             hd_bs_tot=0.0
         endif
         hd_bs_tot_opp=-fact1*hd_bs_tot
           
           
        !construct the digonal hd_bs
        do ns=1,nact(ng)
            n1=iact(ns,ng)
            hd_bs(n1,n1,1)=hd_bs_tot
            hd_bs(iOpposite(n1),n1,1)=hd_bs_tot_opp
            hd_bs(n1,n1,2)=hd_bs_tot
            hd_bs(iOpposite(n1),n1,2)=hd_bs_tot_opp
        enddo 
        
        ! include latent bs matrix
        if (ilatBS.eq.2) then
          do ns1=1,nslsys(iph)
              do ns2=1,nact(ng)
                  n1=iact(ns2,ng)
                  if (iTwinSys(n1).eq.0) then
                      if (ns1.ne.n1.and.ns1.ne.iOpposite(n1)) then
                          if (tau_bcst(n1,ng,1).ge.0.0) then
                              hd_bs(ns1,n1,2) = hd_bs(ns1,n1,2) + 
     #                              aL_bs(ns1,n1,iph)*hd_bs(n1,n1,1)
                          else
                              hd_bs(ns1,n1,2) = hd_bs(ns1,n1,2) + 
     #                              aL_bs(ns1,iOpposite(n1),iph)
     #                              *hd_bs(iOpposite(n1),n1,1)
                          endif
                      endif
                  endif
              enddo
          enddo      
        endif      
      
      endif
      
      end subroutine
      
      subroutine bs_update_statv(ng)
      
      use const
      use mvoigt
      use mphase_props
      use grain_props_v
      use grain_rate, only : gamd,taud_bcst
      use grain_state, only : nact, iact
      use hard_law1_v, only : itemp_act_sys,iOpposite,iDirSys
     #    ,iact_sys
      use flags, only : ilatBS
      
      iph=ngrnph(ng)
      
      if (nact(ng).ne.0.and.ibs.eq.1) then !update is accurate when there is nact(ng)>0
        do ns=1,nslsys(iph)
c            gam_acc(ns,ng)=gam_acc(ns,ng)+gamd(ns,ng) !done in hd_update_statv
        enddo
        !update backstress using hd_bs for active and define inactive as:
        !   tau_bcst (active) > 0 ->  tau_bcst(opposite(active))=-A*tau_bcst(active)
        !   tau_bcst (active) < 0 ->  tau_bcst(opposite(active))=-1/A*tau_bcst(active)
        if (ilatBS.eq.1) then
          do ns=1,nact(ng)
              n1=iact(ns,ng)
              if (iTwinSys(n1).eq.0) then
                  aux=tau_bcst(n1,ng,1)
                  auxOpp=tau_bcst(iOpposite(n1),ng,1)
                  imo=iSysMode(n1,iph)
                  tau_bcst(n1,ng,ilatBS)=tau_bcst(n1,ng,ilatBS)
     #                +taud_bcst(n1,ng,ilatBS)
                  if (tau_bcst(n1,ng,ilatBS).lt.0.0) then
                       fact1=1.0/fact_1(imo,iph)
                  else
                       fact1=fact_1(imo,iph)
                  endif
                  tau_bcst(iOpposite(n1),ng,ilatBS)=-fact1
     #                *tau_bcst(n1,ng,ilatBS)
                  !find taud_bcst_sys
                  if (iphBsCtrl(ngrnph(ng)).eq.1) then
                      taud_bcst_sys(n1,ng)=tau_bcst(n1,ng,1)-aux
                      taud_bcst_sys(iOpposite(n1),ng)=
     #                    tau_bcst(iOpposite(n1),ng,1)-auxOpp
                  endif
              endif
          enddo    
          ! use function
c          do ns=1,nact(ng)
c              n1=iact(ns,ng)
c              if (tau_bcst(n1,ng,1).ge.0.0) then
c                  if(tau_bcst(n1,ng,1).ge.tau_sat(ngrnph(ng)))then !in case backstress overshoots gam_star is NaN
c                      gam_star=0.0
c                  else
c                      gam_star=-1/ni(ngrnph(ng))*log( 1.0-
c     #                       tau_bcst(n1,ng,1)/tau_sat(ngrnph(ng)) )
c                      tau_bcst(n1,ng,1)=bs_org(tau_sat(ngrnph(ng))
c     #                    ,ni(ngrnph(ng)),gam_star+gamd(n1,ng))
c                      fact1=fact_1(ngrnph(ng))
c                  endif
c              endif      
c              
c              if (tau_bcst(n1,ng,1).lt.0.0) then
c                 gam_star_1 = -gam_b(ngrnph(ng))*log(  
c     #                ( tau_sat(ngrnph(ng)) - tau_bcst(n1,ng,1)  ) /
c     #                ( (fact_1(ngrnph(ng))+1.0)*tau_sat(ngrnph(ng)) )  )
c                 tau_bcst(n1,ng,1) = bs_new(tau_sat(ngrnph(ng))
c     #                 ,gam_b(ngrnph(ng)),fact_1(ngrnph(ng))
c     #                 ,gam_star_1+gamd(n1,ng))                
c                 fact1=1.0/fact_1(ngrnph(ng))
c              endif    
c              tau_bcst(iOpposite(n1),ng,1)=-fact1
c     #              *tau_bcst(n1,ng,1)
c        enddo    
          
        endif
        !latent
        if (ilatBS.eq.2) then
          ! update latent using hardening matrix
          do ns=1,nslsys(iph)
            tau_bcst(ns,ng,2)=tau_bcst(ns,ng,2)
     #               +taud_bcst(ns,ng,2)
          enddo          
          ! update per system bs using function
          do ns=1,nact(ng)
              n1=iact(ns,ng)
              if (iTwinSys(n1).eq.0) then
                  aux=tau_bcst(n1,ng,1)
                  auxOpp=tau_bcst(iOpposite(n1),ng,1)
                  imo=iSysMode(n1,iph)
                  if (tau_bcst(n1,ng,1).ge.0.0) then
                      if(tau_bcst(n1,ng,1).ge.tau_sat(imo,iph))then !in case backstress overshoots gam_star is NaN
                          gam_star=0.0
                      else
                          gam_star=-1.0/ni(imo,iph)*log( 1.0-
     #                           tau_bcst(n1,ng,1)/tau_sat(imo,iph) )
                          tau_bcst(n1,ng,1)=bs_org(tau_sat(imo,iph)
     #                        ,ni(imo,iph),gam_star+gamd(n1,ng))
                          fact1=fact_1(imo,iph)
                      endif
                  endif      
                  
                  if (tau_bcst(n1,ng,1).lt.0.0) then
                     gam_star_1 = -gam_b(imo,iph)*log(  
     #                   ( tau_sat(imo,iph) - tau_bcst(n1,ng,1)  ) /
     #                   ( (fact_1(imo,iph)+1.0)*tau_sat(imo,iph) ) )   
                     tau_bcst(n1,ng,1) = bs_new(tau_sat(imo,iph)
     #                     ,gam_b(imo,iph),fact_1(imo,iph)
     #                     ,gam_star_1+gamd(n1,ng))                
                     fact1=1.0/fact_1(imo,iph)
                  endif    
                  tau_bcst(iOpposite(n1),ng,1)=-fact1
     #                  *tau_bcst(n1,ng,1)
                  !find taud_bcst_sys
                  if (iphBsCtrl(ngrnph(ng)).eq.1) then
                      taud_bcst_sys(n1,ng)=tau_bcst(n1,ng,1)-aux
                      taud_bcst_sys(iOpposite(n1),ng)=
     #                    tau_bcst(iOpposite(n1),ng,1)-auxOpp
                  endif
              endif
          enddo 
           
        endif
        
      elseif(ibs.eq.2)then
        do ns = 1,nsys(iph)  
            tau_bcst(ns,ng,2) = tau_bcst(ns,ng,2) + taud_bcst(ns,ng,2)
        enddo
      endif
      
      if (nact(ng).ne.0.and.ibs.eq.3) then !update is accurate when there is nact(ng)>0
          ! update using hardening matrix
          do ns=1,nslsys(iph)
            tau_bcst(ns,ng,1)=tau_bcst(ns,ng,1)
     #               +taud_bcst(ns,ng,1)
            tau_bcst(ns,ng,2)=tau_bcst(ns,ng,2)
     #               +taud_bcst(ns,ng,2)
          enddo
          ! update tot bs in grain using an average over active slip systems
          tau_bcst_tot(ng)=0.0
          do ns=1,nact(ng)
              n1=iact(ns,ng)
              tau_bcst_tot(ng)=tau_bcst_tot(ng)+tau_bcst(n1,ng,1)
     #            /real(nact(ng))
          enddo
          ! update bs for use in other martensite phase
          if (iphBsCtrl(ngrnph(ng)).eq.1) then
              do ns=1,nact(ng)
                  n1=iact(ns,ng)
                  taud_bcst_sys(n1,ng)=taud_bcst(ns,ng,1)
              enddo
          endif
        endif
      
      end subroutine 

      FUNCTION bs_org(tau_sat,ni,gam)      
      
      REAL ni
      !variables are scalars for particular grain and system with same names as arrays
      
      !original backstress evolution
      bs_org=tau_sat*(  1.0 - exp( -ni*gam) )     
      
      END FUNCTION

      FUNCTION bs_org_der(tau_sat,ni,gam)      
      
      REAL ni
      !variables are scalars for particular grain and system with same names as arrays
      
      !original backstress evolution
      bs_org_der=ni*tau_sat*exp( -ni*gam)      
      
      END FUNCTION

      FUNCTION bs_new(tau_sat,gam_b,fact,gam)    
      
      !variables are scalars for particular grain and system with same names as arrays
      
      !new backstress evolution
      bs_new = -(fact+1.0)*tau_sat
     #            *exp(-gam/gam_b)+tau_sat
      
      END FUNCTION
      
      FUNCTION bs_new_der(tau_sat,gam_b,fact,gam)    
      
      !variables are scalars for particular grain and system with same names as arrays
      
      !new backstress evolution
      bs_new_der = (1.0+fact)*tau_sat/gam_b
     #            *exp(-gam/gam_b)
      
      END FUNCTION
      
      SUBROUTINE rod_esh(r,CA)
C **********************************************************************
C *** input: crystal orientation                                     ***
C *** output: (average ???) eshelby tensor for three variants        ***
C ***         of precipitates multiplied with isotropic stiffness    ***
C ***         of Al in sample frame                                  ***
c **********************************************************************
c **********************************************************************
      
      use mvoigt
      
      DIMENSION S(3,3,3,3),S1(3,3,3,3),S2(3,3,3,3),S3(3,3,3,3),SV(6,6)
     #          ,aux6(6),aux31(3,3),aux32(3,3),aux33(3,3),r(3,3),C(6,6)
     #          ,aux66(6,6),CA(6,6)
      REAL nu,mu,lambda,r

      
      nu=0.345
      E = 70300.0
      mu = E/2.0/(1.0 + nu)
      lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu)
      C = 0.0
      
      ! elastic stiffnes of Al
      do i=1,3
          do j=1,3
              C(i,j) = lambda
          enddo
      enddo
      do i=1,6
          if(i.le.3) C(i,i) = C(i,i) + 2.0*mu
          if(i.gt.3) C(i,i) = mu          
      enddo
            
      ! eshelby tensor
      
      S = 0.0 
      S1 = 0.0
      SV1 = 0.0
      S2 = 0.0
      SV2 = 0.0      
      S3 = 0.0
      SV3 = 0.0
      
      ! find the eshelby tensor for precipitate along z axis
      S(1,1,1,1) = (5.0-4.0*nu)/8.0/(1.0-nu)
      S(2,2,2,2) = S(1,1,1,1)
      S(1,1,2,2) = (4.0*nu-1.0)/8.0/(1.0-nu)
      S(2,2,1,1) = S(1,1,2,2)
      S(1,1,3,3) = nu/2.0/(1.0-nu)
      S(2,2,3,3) = S(1,1,3,3)
      S(1,2,1,2) = (3.0-4.0*nu)/8.0/(1.0-nu)
      S(1,2,2,1) = S(1,2,1,2)
      S(2,1,2,1) = S(1,2,1,2)
      S(2,1,1,2) = S(1,2,1,2)
      S(1,3,1,3) = 0.25
      S(1,3,3,1) = 0.25
      S(3,1,3,1) = 0.25
      S(3,1,1,3) = 0.25
      S(2,3,2,3) = 0.25
      S(2,3,3,2) = 0.25
      S(3,2,3,2) = 0.25
      S(3,2,2,3) = 0.25
      
      ! (var1) find the eshelby tensor for precipitate along z axis
      ! transformation matrix from cylinder frame to sample frame
      ! transformation matrix from cylinder frame to crystal frame
      aux31 = 0.0
      aux31(1,1) = 1.0
      aux31(2,2) = 1.0
      aux31(3,3) = 1.0
      ! transformation matrix from cylinder frame to sample frame
      aux32 = matmul(transpose(r),aux31)
      
      DO I=1,3
      DO J=1,3
      DO M=1,3
      DO N=1,3
        DUMMY=0.
        DO I1=1,3
        DO J1=1,3
        DO M1=1,3
        DO N1=1,3
          DUMMY=DUMMY+aux32(I,I1)*aux32(J,J1)*aux32(M,M1)
     #         *aux32(N,N1)*S(I1,J1,M1,N1)
        END DO
        END DO
        END DO
        END DO
        S1(I,J,M,N)=DUMMY
      END DO
      END DO
      END DO
      END DO
      
      ! (var2)  transform eshelby tensor from cylinder frame to crystal
      ! transformation matrix from cylinder frame to crystal frame
      aux31 = 0.0
      aux31(1,1) = 1.0
      aux31(3,2) = -1.0
      aux31(2,3) = 1.0
      ! transformation matrix from cylinder frame to sample frame
      aux32 = matmul(transpose(r),aux31)
      
      DO I=1,3
      DO J=1,3
      DO M=1,3
      DO N=1,3
        DUMMY=0.
        DO I1=1,3
        DO J1=1,3
        DO M1=1,3
        DO N1=1,3
          DUMMY=DUMMY+aux32(I,I1)*aux32(J,J1)*aux32(M,M1)
     #         *aux32(N,N1)*S(I1,J1,M1,N1)
        END DO
        END DO
        END DO
        END DO
        S2(I,J,M,N)=DUMMY
      END DO
      END DO
      END DO
      END DO
      
      ! (var3)  transform eshelby tensor from cylinder frame to crystal
      ! transformation matrix from cylinder frame to crystal frame
      aux31 = 0.0
      aux31(3,1) = -1.0
      aux31(2,2) = 1.0
      aux31(1,3) = 1.0
      ! transformation matrix from cylinder frame to sample frame
      aux32 = matmul(transpose(r),aux31)
      
      DO I=1,3
      DO J=1,3
      DO M=1,3
      DO N=1,3
        DUMMY=0.
        DO I1=1,3
        DO J1=1,3
        DO M1=1,3
        DO N1=1,3
          DUMMY=DUMMY+aux32(I,I1)*aux32(J,J1)*aux32(M,M1)
     #         *aux32(N,N1)*S(I1,J1,M1,N1)
        END DO
        END DO
        END DO
        END DO
        S3(I,J,M,N)=DUMMY
      END DO
      END DO
      END DO
      END DO      
      
      ! find average eshelby tensor
      S1 = 1.0/3.0*(S1 + S2 + S3)
      call VOIGT(aux6,aux33,SV,S1,4)
      
      call TENS_MULT(CA,C,(SV-id2))      
      
      RETURN
      END SUBROUTINE
      
      subroutine ph_bcst
      !volume average backstress over iphBsCtrl==1 phase is used to 
      !define volume average backstress over iphBsCtrl==0 phases
      
      use const
      use mvoigt
      use flags, only : iLatBS
      use mphase_props_v
      use sample_props_v
      use grain_state_v, only : iact,nact
      use grain_props_v, only : mcs,wgt
      use grain_rate_v, only : taud_bcst
      
      real :: strss_bs(6,NPHM),strssva_bsCting(6),strss_bsCted(6)
      
      !volume average of increment bs in iphBsCtrl==1
      strssva_bsCting=0.0
      do iph=1,nph
        if (iphBsCtrl(iph).eq.1) then
          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2=SUM(nphngr(1:iph))
          strss_bs(:,iph)=0.0
          do ng=ng1,ng2
            do ns=1,nact(ng)
              n1=iact(ns,ng)
              imo=iSysMode(n1,iph)
              if (tau_bcst(n1,ng,1).gt.0.0) then
                strss_bs(:,iph)=strss_bs(:,iph)+ 2.0*mcs(:,n1,ng)
     #             *taud_bcst_sys(n1,ng)*wgt(ng)
              else
                strss_bs(:,iph)=strss_bs(:,iph)+ 2.0*mcs(:,n1,ng)
     #             *taud_bcst_sys(n1,ng)/fact_1(imo,iph)*wgt(ng)
              endif
            enddo   
          enddo
          strssva_bsCting=strssva_bsCting+strss_bs(:,iph) 
        endif
      enddo
      
      !volume average of increment bs in iphBsCtrl==0 and nact==0
      wgt_phCted=0.0
      do iph=1,nph
        if (iphBsCtrl(iph).eq.-1) then
          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2=SUM(nphngr(1:iph))
          do ng=ng1,ng2
            if (nact(ng).eq.0) then
              wgt_phCted=wgt_phCted+wgt(ng)
            endif
          enddo
        endif
      enddo
      if(wgt_phCted.ne.0.0) strss_bsCted=-strssva_bsCting/wgt_phCted
      
      !calcualte taud_bcst in iphBsCtrl==0
      do iph=1,nph
        if (iphBsCtrl(iph).eq.-1) then
          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2=SUM(nphngr(1:iph))
          do ng=ng1,ng2
            taud_bcst(1:nslsys(ngrnph(ng)),ng,iLatBS)=0.0
            if (nact(ng).eq.0) then
              do ns=1,nslsys(ngrnph(ng))  
                do i=1,6  
                  taud_bcst(ns,ng,iLatBS)=taud_bcst(ns,ng,iLatBS)+ 
     #                mcs(i,ns,ng)*strss_bsCted(i)*profac(i)
                enddo
                !update bcst
                tau_bcst(ns,ng,iLatBS)=tau_bcst(ns,ng,iLatBS)+
     #                taud_bcst(ns,ng,iLatBS)
              enddo
            endif
          enddo
        endif
      enddo      
      
c      write(194,'(100e25.5)') tau_bcst(1,1,iLatBs),strssva_bsCting
c     #    ,-strss_bsCted*wgt_phCted
c      write(194,'(100e25.5)') strssva_bsCting
c     #  ,-strss_bsCted,wgt_phCted
      
      end subroutine
      
      end module
c      
c***********************************************************************
c *** phase_transf *****************************************************
c***********************************************************************
      module phase_transf
      
      use phase_transf_v
      
      contains
      
      subroutine set_pt(ng,pt)
      
      use twinning_v, only : iParentGrain,iChildGrain
      use grain_rate_v, only : etrcs_eig
      use grain_state_v, only : gamtot
      
      integer, intent(in) :: ng
      type(state_pt),intent(out) :: pt
 
      pt%etrcs_eig=etrcs_eig(:,ng)
      pt%etcs_pt=etcs_pt(:,ng)
      pt%iParentGrain=iParentGrain(ng)
      pt%iChildGrain=iChildGrain(1,ng)
      pt%gamtot=gamtot(ng)
      
      end subroutine set_pt
      
      subroutine get_pt(ng,pt)
      
      use twinning_v, only : iParentGrain,iChildGrain
      use grain_rate_v, only : etrcs_eig
      use grain_state_v, only : gamtot
      
      integer, intent(in) :: ng
      type(state_pt),intent(in) :: pt
      
      etrcs_eig(:,ng)=pt%etrcs_eig
      etcs_pt(:,ng)=pt%etcs_pt
      iParentGrain(ng)=pt%iParentGrain
      iChildGrain(1,ng)=pt%iChildGrain
      gamtot(ng)=pt%gamtot
      
      end subroutine get_pt
      
      subroutine write_statv_pt(pt,ng,ns,NSTATV,STATEV)
      
      use twinning_v, only : iParentGrain,iChildGrain
      use grain_rate_v, only : etrcs_eig
      use grain_state_v, only : gamtot
      use mphase_props
      
      type(state_pt),intent(in) :: pt
      integer,intent(in) :: ng,NSTATV
      integer,intent(inout) :: ns
      real,intent(out) :: STATEV(NSTATV)
      
      STATEV(ns:ns+6-1)=pt%etrcs_eig
      ns=ns+6
      STATEV(ns:ns+6-1)=pt%etcs_pt
      ns=ns+6
      STATEV(ns)=pt%iParentGrain
      ns=ns+1
      STATEV(ns)=pt%iChildGrain
      ns=ns+1
      STATEV(ns)=pt%gamtot
      ns=ns+1
      
      end subroutine write_statv_pt
      
      subroutine read_statv_pt(pt,ng,iph,ns,NSTATV,STATEV)
      
      use twinning_v, only : iParentGrain,iChildGrain
      use grain_rate_v, only : etrcs_eig
      use grain_state_v, only : gamtot
      use mphase_props
      
      type(state_pt),intent(out) :: pt
      integer,intent(in) :: ng,NSTATV,iph
      real,intent(in) :: STATEV(NSTATV)
      integer,intent(inout) :: ns

      pt%etrcs_eig=STATEV(ns:ns+6-1)
      ns=ns+6
      pt%etcs_pt=STATEV(ns:ns+6-1)
      ns=ns+6
      pt%iParentGrain=STATEV(ns)
      ns=ns+1
      pt%iChildGrain=STATEV(ns)
      ns=ns+1
      pt%gamtot=STATEV(ns)
      ns=ns+1
      
      end subroutine read_statv_pt

      subroutine phtr_read(filetemp)
      
      DIMENSION aux33(3,3),aux6(6)
      character(len=40) filetemp
      
      !read parameter, initialize arrays
      open(unit=1,file=filetemp)
      read(1,*) alpha0,alphaK,beta0,betaK,anexp,sfe,wgtcr
      call variant_selection(0,idummy,idummy,aux33,aux6,dummy)
      close(1)
      
      end subroutine 
      
      subroutine wgtd_calc(wgtd)
      
      use const
      use twinning_v, only : iChildGrain,iParentGrain
      use mphase_props, only : nphngr,ngrnph
      use grain_props_v, only : wgt
      use grain_state, only : iact,nact,gamtot
      use grain_rate, only : gamd
      use sample_props_v, only : ngrain
      
      real, intent(out) :: wgtd(NGR)
      
      LOGICAL nuc
      DIMENSION aux33(3,3),aux6(6)
     
      do ng=1,nphngr(1)
          !find alpha and beta
          alpha=alpha_escgr(alpha0,alphaK,escgr(ng))    
          beta=beta_triax(beta0,betaK,triax(ng))
          !restart weight increment
          wgtd(ng)=0.0
          if(iChildGrain(1,ng).ne.0.0) wgtd(iChildGrain(1,ng))=0.0
          !without children - nucleate
          if (iChildGrain(1,ng).eq.0) then
              wgt_aust0=wgt(ng)
              !nucleation criteria
              nuc=nucleation(alpha,beta,anexp,wgtcr,wgt_aust0
     #        ,gamtot_phtr(ng))
              if(nuc) then
                  !add new grain
                  ngrain=ngrain+1
                  nphngr(2)=nphngr(2)+1
                  iChildGrain(1,ng)=ngrain
                  iParentGrain(ngrain)=ng
                  ngP=ng
                  ngC=ngrain
                  ngrnph(ngC)=2
                  call add_grain(ngP,ngC,wgtcr)
              endif
          endif
          !with children - evolve weight
          if (iChildGrain(1,ng).ne.0) then
              wgt_aust0=wgt(ng)+wgt(iChildGrain(1,ng))
              !calculate wgtd
              ! linear law
c              aux1 = 0.0
c              do i=1,6
c                  aux1=aux1+etrss(i)**2
c              enddo
c              aux1=0.4*sqrt(aux1)                     
c              wgtd(ng)=-aux1
c              wgtd(iChildGrain(1,ng))=aux1
              !OC law
              aux1=0.0
              do ns1=1,nact(ng)
                  n1=iact(ns1,ng)
                  aux1=aux1+wgtd_martOC_det(alpha,beta,anexp
     #                     ,wgt_aust0,gamtot_phtr(ng))*gamd(n1,ng)
              enddo
              wgtd(ng)=-aux1
              wgtd(iChildGrain(1,ng))=aux1
              
              !prevent negative volume fraction
              if (wgt(ng)+wgtd(ng).le.0.0) then
                  wgtd(ng)=0.0
                  wgtd(iChildGrain(1,ng))=0.0
              endif
          endif
      enddo
      
      end subroutine wgtd_calc
      
      FUNCTION nucleation(alpha,beta,anexp,wgtcr,wgt_aust0,et)
      
      LOGICAL nucleation
      
      !critical strain criteria      
c      aux = 0.0
c      do i=1,6
c          aux=aux+et(i)**2
c      enddo
c      aux=sqrt(aux)
c      
c      if(aux.gt.0.03) then
c          nucleation=.true.
c      else
c          nucleation=.false.
c      endif
      
      !critical volume fraction criteria
      if (wgt_martOC(alpha,beta,anexp,wgt_aust0,et).gt.wgtcr) then
          nucleation=.true.
      else
          nucleation=.false.
      endif
      
      RETURN
      END FUNCTION
      
      SUBROUTINE add_grain(ngP,ngC,wgtcr)
      
      use flags
      use miscellaneous_sub, only : euler
      use grain_props_v, only : phi,the,ome,r
      use grain_props, only : cr_to_sa
      use grain_state, only : stcs,etcs,tau
      use twinning_v, only : iParentGrain
      use grain_state_v, only : fijgr,fijgr_pt
      use hard_law1, only : hl_ini_tau1=>hl_ini_tau
      use hard_law2, only : hl_ini_tau2=>hl_ini_tau
      
      DIMENSION aux33(3,3),r_SV(3,3),etcs_ptSV(6),hn(3)
      DIMENSION Ram_e(3,3),U_e(3,3),e1(3),e2(3),aux3(3)
      
      !single crystal variables:
      ! crystal properties
      !     ccs2(ngrain),scs2(ngrain),bcs(ngrain),ncs(ngrain),mcs(ngrain)
      ! orientation and weight
      !     r(ngrain),phi(ngrain),the(ngrain),ome(ngrain)
      !     wgt(ngrain)
      ! stress and strain in grains
      !     stcs(ngrain),etcs(ngrain)
      ! hard law
      !     tau(ngrain),rho_rev(ngrain),rho_forw(ngrain)
      !     rho_tot(ngrain),rho_tot_max(ngrain),rho_deb(ngrain)
      !     rho_act(ngrain),tau_act(ngrain),iact_sys(ngrain)
      ! backstress law
      !     itemp_act_sys(nsys,ngrain),tau_0b(nsys,ngrain)
      !     gam_0b(nsys,ngrain),i_memory_1(nsys,ngrain)
      !     i_memory(nsys,ngrain),gam_acc(nsys,ngrain)
      !     tau_bcst(nsys,ngrain),aL_bs(nsys,nsys,nph),
      !     tau_bcst_l(nsys,ngrain)
      !     stss_mbs(6),stss_fbs(6),stss_lbs(6)
      ! elipsoid shape
      !     eulerph(nph),axisph(nph),fijph(nph),aef(nph)
      !     auxsample(nph)
      ! twinning
      !     volume fraction: TWFRSY(ns,ng),TWFRGR(ng),WGTCG_all(ng)
      !     PTS: KTWSMX(ns,ng)
      !     just_created(ng),create_twin(ng),igr_exist(ng),iChildGrain(ns,ng)
      !     iParentGrain(ng),iParentMode(ng),IPTSGRC(ng)
      !     PTSDAMP(NSLS,NGR),PTSDAMP_OLD(NSLS,NGR),tau0_tw(NSLS,NGR)      

      ! orientation and weight
      ! variant selection
      call variant_selection(1,ngP,iSelVar,r_SV,etcs_ptSV,E_max)
      r(:,:,ngC)=r_SV
      etcs_pt(:,ngC)=etcs_ptSV
      call euler(1,phi(ngC),the(ngC),ome(ngC),r_SV)
c      wgt(ngC)=wgtcr
      
      ! crystal properties
      call cr_to_sa(ngC,ngC,0)
      call cr_to_sa(ngC,ngC,1) 
      
      ! stress and strain in grains
      stcs(:,ngC)=stcs(:,iParentGrain(ngC))
      etcs(:,ngC)=etcs(:,iParentGrain(ngC))
      
      ! initialize shape
      if (iShape.eq.2) then
          !habit plane in austenite frame
          hn_a(1)=-0.1846
          hn_a(2)=-0.7821
          hn_a(3)=-0.5951
          !habit plane in sample frame for selected variant
          aux33=matmul(O_cubic_sym(:,:,iSelVar),r(:,:,ngP))
          hn=matmul(transpose(aux33),hn_a)
          !e1 and e2 axes in sample frame
          e1_a=(/-0.9733,0.2297,0.0/)
          e2_a=(/0.1367,0.5792,-0.8036/)
          
          e1=matmul(transpose(aux33),e1_a)
          e2=matmul(transpose(aux33),e2_a)
          !rotation from sphere in sample axes to mart ellipsoid
          Rsm_e(:,1,ngC)=e1
          Rsm_e(:,2,ngC)=e2
          Rsm_e(:,3,ngC)=hn
c          aux3=(/0.0,0.0,1.0/)
c          hn=matmul(Rsm_e(:,:,ngC),aux3) !check if rotation of sample z axis is to hn
          !streching tensor applied to sphere in sample frame
c          U_e(1,1)=axis(1,2)
c          U_e(2,2)=axis(2,2)
c          U_e(3,3)=axis(3,2)
c          !initial deformation gradient for phase deformation
c          fijgr_pt(:,:,ngC)=matmul(Rsm_e(:,:,ngC),U_e)
c          !initial def grad for def
          do i=1,3
              fijgr(i,i,ngC)=1.0
          enddo
      endif
      
      ! hard law - defined by 2nd phase
      if (kCL.eq.1) call hl_ini_tau1(ngC,tau,0)
      if (kCL.eq.2) call hl_ini_tau2(ngC,tau,0)
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE variant_selection(iopt,ngP,iSelVar,r_SV,etcs_ptSV
     #                ,E_max)
      
      use mvoigt
      use grain_props_v, only : r
      use grain_state, only : stcs
      
      DIMENSION r_m(3,3,24),etcs_pt33(3,3,24)
      DIMENSION etcs_ptV(6,24),E_i(24),r_SV(3,3),etcs_ptSV(6)
      DIMENSION aux3333(3,3,3,3),aux66(6,6),aux31(3,3)

      ! initialize arrays
      if (iopt.eq.0) then
          call cubic_sym_operators !get all symmetry operators
          call transformation_a_m(r_m) !get all variant transformation matrices
          call phase_trans_et(etcs_pt33)! get all phase transformation strains
      endif
      
      ! variant selection
      if (iopt.eq.1) then
          ! change frame from austenite to sample for phase trans strain
          ! (r is sample to crystal transformation)
          do i=1,24
              aux31 = matmul(etcs_pt33(:,:,i),r(:,:,ngP))
              aux31 = matmul(transpose(r(:,:,ngP))
     #                           ,aux31)
              ! change to voigt notation 
              call VOIGT(etcs_ptV(:,i),aux31,aux66,aux3333,2)
          enddo
          
          ! find max mechanical driving force to select the variant
          do i=1,24
              aux = 0.0
              do j=1,6
                  aux = aux + stcs(j,ngP)*etcs_ptV(j,i)*profac(j)
              enddo
              E_i(i) = aux
          enddo
          E_max = maxval(E_i,1)
          iSelVar = maxloc(E_i,1)
          
          ! define output: iSelVar, corresponding transformation matrix and phase transformation strain
          r_SV = matmul(r_m(:,:,iSelVar),r(:,:,ngP))
          etcs_ptSV=etcs_ptV(:,iSelVar)
      endif
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE cubic_sym_operators
      
      DIMENSION temp(216),temp1(72,3)
      
      !for an input array populates it with cubic symmetry operators
      ! temporary array containing information about symmetry operators
      temp=(/1,0,0,0,0,-1,0,0,1,-1,0,0,0,0,-1,1,0,0,1,0,0,1,0,0,0,1
     #,0,-1,0,0,0,-1,0,0,1,0,0,0,1,0,-1,0,0,0,-1,0,0,-1,0,1,0,0,-1
     #,0,0,0,1,0,1,0,-1,0,0,0,0,1,0,-1,0,-1,0,0,0,1,0,0,-1,0,0,1,0
     #,0,1,0,0,1,0,0,0,1,0,-1,0,0,0,-1,-1,0,0,0,-1,0,1,0,0,0,0,1,1
     #,0,0,0,0,1,-1,0,0,1,0,0,0,0,-1,0,0,-1,-1,0,0,1,0,0,0,0,1,0,-1
     #,0,-1,0,0,0,0,-1,0,0,1,-1,0,0,-1,0,0,0,0,-1,1,0,0,0,-1,0,0,0
     #,-1,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,0,1,0,-1,0,0,0,1,0,0,-1,0
     #,-1,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,1,0,0,0,0,-1,0,-1,0/)

      ! reshape the temp
      temp1 = reshape(temp,(/72, 3/))
      
      ! populate O_cubic_sym
      do i=1,24
          O_cubic_sym(:,:,i) = temp1(((i-1)*3+1):3*i,1:3)
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE transformation_a_m(r_m)
      
      DIMENSION r_m(3,3,24),r_m1(3,3)
      
      r_m1(1,1:3) = (/0.7239, -0.6896, -0.0183/)
      r_m1(2,1:3) = (/0.6778, 0.716, -0.167/)
      r_m1(3,1:3) = (/0.1283, 0.1085, 0.9858/)
      
      !get all martensite transformation matrices 
      do i=1,24
          r_m(:,:,i) = matmul(r_m1,O_cubic_sym(:,:,i))
      enddo      
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE phase_trans_et(etcs_pt33)
      
      DIMENSION etcs_pt33(3,3,24),etcs_pt331(3,3)
      DIMENSION aux33(3,3)
      
      etcs_pt331(1,1:3) = (/-0.0078,0.0,-0.0253/)
      etcs_pt331(2,1:3) = (/0.0,0.1408,0.0/)
      etcs_pt331(3,1:3) = (/-0.0253,0.0,-0.0815/)
      
      !get all martensite transformation  strains in austenite frame
      do i=1,24
          aux33 = matmul(etcs_pt331,O_cubic_sym(:,:,i))
          etcs_pt33(:,:,i) = matmul(transpose(O_cubic_sym(:,:,i)),aux33)
      enddo      
      
      RETURN
      END SUBROUTINE
      
      FUNCTION wgt_martOC(alpha,beta,anexp,wgt_aust0,et)
      
      wgt_martOC = 
     #wgt_aust0*( 1.0 - exp(-beta* (1.0 - exp(-alpha*et) )**anexp ) )
      
      END FUNCTION
      
      FUNCTION wgtd_martOC_det(alpha,beta,anexp,wgt_aust0,et)
      
      wgtd_martOC_det = alpha*beta*wgt_aust0*anexp
     #*exp( -beta* (1.0 - exp(-alpha*et) )**anexp - alpha*et)
     #*(1.0-exp(-alpha*et))**(anexp-1.0)
      
      END FUNCTION
      
      FUNCTION alpha_escgr(alpha0,alphaK,escgr)
      
      alpha_escgr=alpha0+alphaK*escgr
      if (alpha_escgr.lt.0.0) alpha_escgr=0.0
      
      END FUNCTION 
      
      FUNCTION beta_triax(beta0,betaK,triax)
      
      beta_triax=beta0+betaK*triax
      
      END FUNCTION
      
      SUBROUTINE eig_calc
      
      use mphase_props, only : nphngr
      use grain_props, only : wgt
      use twinning_v, only : iChildGrain,iParentGrain
      use grain_rate_v, only : etrcs_eig,wgtd
      
      
      do ngP=1,nphngr(1)
          ngC=iChildGrain(1,ngP)
          if(ngC.ne.0)then
              !split eigenstrain between parent and child
              etrcs_eig(:,ngP)=wgtd(ngC)*etcs_pt(:,ngC)
     #         /(wgt(ngP)+wgt(ngC))
              etrcs_eig(:,ngC)=wgtd(ngC)*etcs_pt(:,ngC)
     #         /(wgt(ngP)+wgt(ngC))
              !give eigenstrain to parent only
c              etrcs_eig(:,ngP)=wgtd(ngC)/wgt(ngP)*etcs_pt(:,ngC)
          endif
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE SFW_calc(ng1,ng2)
      
      use const
      use mvoigt
      use grain_props_v, only : ncs,bcs,mcs,wgt
      use grain_state_v, only : stcs,nact,iact
      use grain_rate_v, only : gamd
      use mphase_props_v, only : nslsys,ngrnph,nsm,wgt_ph
      use output_v, only : etsspleq,TEMP
      
      DIMENSION sfw(NSLS,NGR),rss(NSLS),ics(NSLS),sfwn(NSLS,NGR)
      DIMENSION esc(NSLS,NGR),aux3(3),aux33(3,3),esceig(NSLS,NGR)
      DIMENSION aux66(6,6),aux3333(3,3,3,3),aux6(6),sig(24)
      DATA (sig(i),i=1,24)/-1.,1.,       1.,-1.,        -1.,1.   
     #                   ,1.,-1.,        -1.,1.,       1.,-1.   
     #                   ,1.,-1.,        -1.,1.,        -1.,1.   
     #                    ,-1.,1.,       1.,-1.,       1.,-1./

      DIMENSION sfwngr(NGR),sbsgr(NSLS,NGR),sbgr(NGR)
      
      !define constants
c      sfe=11.8 !mJ/m^2
c      sfe=20.0 !mJ/m^2
      a=0.36469 !nm
      bp=a/sqrt(6.0)
      amu=77000.0 !MPa
      anu=0.2885
      theta=90.0 !edge dislocation
      theta1=(90.0+30.0)
      theta2=(90.0-30.0)
      pi=4.0*atan(1.0)
      theta1=theta1*pi/180.0
      theta2=theta2*pi/180.0
      fT1T2=cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)
     #  /(1.0-anu)
      c=amu*bp**2*fT1T2/pi !MPa*nm^2
      !equilibrium sfw
      sfw_eq=c/2.0/sfe
      do ng=ng1,ng2
        do ns=1,nslsys(ngrnph(ng))
          !find non shcmid tensor 
          aux3(1)=ncs(2,ns,ng)*bcs(3,ns,ng)-ncs(3,ns,ng)*bcs(2,ns,ng)
          aux3(2)=ncs(3,ns,ng)*bcs(1,ns,ng)-ncs(1,ns,ng)*bcs(3,ns,ng)
          aux3(3)=ncs(1,ns,ng)*bcs(2,ns,ng)-ncs(2,ns,ng)*bcs(1,ns,ng)
          do i=1,3
            do j=1,3
              aux33(i,j)=sig(ns)*1.0/2.0*(aux3(i)*ncs(j,ns,ng)
     #             +aux3(j)*ncs(i,ns,ng))
            enddo
          enddo
          call VOIGT(aux6,aux33,aux66,aux3333,2)
          !find tau2-tau1 (esceig stress)
          aux=0.0
          auxa=0.0
          do i=1,6
            aux=aux+stcs(i,ng)*aux6(i)*profac(i)
            auxa=auxa+stcs(i,ng)*mcs(i,ns,ng)*profac(i)
          enddo
          !store aux
          esceig(ns,ng)=aux
          !find if shear band can be created on a slip system (2*sfe+(tau2-tau1)<=0)
          sbsgr(ns,ng)=0.0
          if (2.0*sfe-aux*bp.le.0.0.and.auxa.gt.0.0) sbsgr(ns,ng)=1.0
          if (auxa.ne.0.0) then
            esc(ns,ng)=aux/auxa
          else 
            esc(ns,ng)=0.0
          endif
          !find sfw and sfwn
          if (auxa.gt.0.0) then
            sfw(ns,ng)=c/(2.0*sfe-bp*aux) 
            sfwn(ns,ng)=sfw(ns,ng)/sfw_eq
            if (sfwn(ns,ng).lt.0.0.or.sfwn(ns,ng).gt.10.0) then
              sfwn(ns,ng)=10.0
            endif
          else
            sfw(ns,ng)=0.0
            sfwn(ns,ng)=0.0
          endif
        enddo
      enddo
        
      !average over active slip systems
      asfw=0.0
      asfwn=0.0
      wgt_sfw=0.0
      nsfw=0
      do ng=ng1,ng2
        sfwngr(ng)=0.0
        sbgr(ng)=0.0
        escgr(ng)=0.0
        asfw_gr=0.0
        if (nact(ng).ne.0) then
          do ns=1,nact(ng)
            n1=iact(ns,ng)
            escgr(ng)=escgr(ng)+esc(n1,ng)/real(nact(ng)) !if(esc(n1,ng).le.0.0)
            sbgr(ng)=sbgr(ng)+sbsgr(n1,ng)/real(nact(ng))
            sfwngr(ng)=sfwngr(ng)+sfwn(n1,ng)/real(nact(ng))
            asfwn=asfwn+sfwn(n1,ng)*wgt(ng)/real(nact(ng))
            if(sfwn(n1,ng).eq.10.0) then
              nsfw=nsfw+1
            endif
          enddo
          wgt_sfw=wgt_sfw+wgt(ng)
        endif
      enddo
      if (wgt_sfw.ne.0.0) then
        asfwn=asfwn/wgt_sfw
      endif
      
      !shear band onset
      if(sum(sbgr(ng1:ng2)).ne.0.0)sbgr(ng1:ng2)=1.0
      do ng=ng1,ng2
          if (sbgr(ng).gt.0.0) then
              do ns1=1,nact(ng)
                  n1=iact(ns1,ng)
                  gamtot_phtr(ng)=gamtot_phtr(ng)+gamd(n1,ng)
              enddo
          endif
      enddo      
      
      !percent of slip systems with diverged sfw
      fnsfw=real(nsfw)/real(nsm(1,1))/real(ng2-ng1)*2.0
      
      !volume average
      escavg=0.0
      do ng=ng1,ng2
        escavg=escavg+escgr(ng)*wgt(ng)
      enddo
      escavg=escavg/wgt_ph(1)   
      
      TEMP(1)=escavg
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE lode_triax(ng1,ng2)
      
      use output_v, only : TEMP
      use const
      use mvoigt
      use grain_props_v, only : wgt
      use grain_state_v, only : stcs
      use miscellaneous_sub, only : det
      
      DIMENSION st33(3,3),stdev(3,3),sthyd(3,3),aux66(6,6),aux33(3,3)
      DIMENSION aux3333(3,3,3,3),alode(NGR)
      pi=4.0*atan(1.0)
      
      triaxavg=0.0
      alodeavg=0.0
      
      do ng=ng1,ng2
        !find hydrostatic and deviatoric stresses
        call VOIGT(stcs(:,ng),st33,aux66,aux3333,1)
        aux=( st33(1,1)+st33(2,2)+st33(3,3) )/3.0
        sthyd=0.0
        do i=1,3
          sthyd(i,i)=aux
        enddo
        stdev=st33-sthyd
        !find invariants J2 and J3 of dev stress
        aJ2=0.0
        do i=1,3
          do j=1,3
            aJ2=aJ2+0.5*stdev(i,j)*stdev(i,j)
          enddo
        enddo
        aJ3=det(stdev)
        steq=sqrt(aJ2)*sqrt(3.0)
        
        !find triaxiality and lode angle
        if (steq.ne.0.0) then
          triax(ng)=aux/steq
          alode(ng)=1.0-2.0/pi*acos(3.0*sqrt(3.0)/2.0*aJ3/sqrt(aJ2**3))
        else
          triax(ng)=0.0
          alode(ng)=0.0
        endif

        triaxavg=triaxavg+wgt(ng)*triax(ng)
        alodeavg=alodeavg+wgt(ng)*alode(ng)
        
        TEMP(2)=alodeavg
        TEMP(3)=triaxavg
        
      enddo
      
      END SUBROUTINE
      
      subroutine update_fijgr_pt
      
      use mphase_state_v, only : axis
      use mphase_props_v, only : nphngr
      use grain_state_v, only : fijgr_pt
      use twinning_v, only : iParentGrain
      use grain_props_v, only : wgt
      
      dimension U_e(3,3)
      
      !total def grad at t+dt due to phase transormation
      ng1=SUM(nphngr(1:2))-nphngr(2)+1
      ng2=SUM(nphngr(1:2))
      do ng=ng1,ng2
          ngC=ng
          ngP=iParentGrain(ngC)
          !martensite
          U_e(1,1)=axis(1,2)
          U_e(2,2)=axis(2,2)
          U_e(3,3)=axis(3,2)+(1.0-axis(3,2)*2.0)
     #        *wgt(ngC)/(wgt(ngP)+wgt(ngC))
          fijgr_pt(:,:,ngC)=matmul(Rsm_e(:,:,ngC),U_e)
          !austenite
          U_e(1,1)=axis(1,1)
          U_e(2,2)=axis(2,1)
          U_e(3,3)=1.0-U_e(3,3)
          fijgr_pt(:,:,ngP)=matmul(Rsm_e(:,:,ngC),U_e)
      enddo
          
      end subroutine
      
      subroutine initialize_fijgr_pt(iph)
      
      use flags
      use mphase_state_v, only : eulerph,axis
      use grain_state_v, only : fijgr_pt    
      use sample_props_v, only : ngrain
      
      dimension aux33(3,3)
      
      if (iph.eq.1.and.(eulerph(1,1).ne.0.0
     #    .or.eulerph(2,1).ne.0.0
     #    .or.eulerph(3,1).ne.0.0
     #    .or.axis(1,1).ne.1.0
     #    .or.axis(2,1).ne.1.0      
     #    .or.axis(3,1).ne.1.0)) then
          write(*,*) 'MUST USE SPHEREICAL SHAPE WITH IPHTR=1 AND ISHAPE=
     #2'
          STOP  
      endif      
      do ng=1,ngrain
          do i=1,3
              fijgr_pt(i,i,ng)=1.0
          enddo
      enddo
      
      end subroutine
      
      end module phase_transf
c      
c***********************************************************************
c *** sc_estimate ******************************************************
c***********************************************************************
      module sc_estimate
      
      private 
      public sc_new,av_modulus,loc_tens
      
      contains
      !sc_new       :Solves the iterative SC equation

      SUBROUTINE sc_new(iopt,liter,iflagcond,e2,interaction,istep,ng1,
     #    ng2,etrss,alfass,ass2,ilevel,iSkip_ph,
     #    axisgr,axisph,css2,einvsa,einvsagr,escr4,escr4gr,aef,
     #    aefgr,acs2,alfacs,aloc2,itmax_mod,ngrain,wgt,a_guess,
     #    error_mod,meffc,etrcs_eig,etrss_eig,xmmin_in)
c **********************************************************************
c *** For given 'acs2' and 'alfacs':
c *** If iopt=0 solves selfconsistent equation for the sample elasto-
c     plastic 'ass2' and thermal 'alfass' moduli ('itmax_mod' iterations)
c *** If iopt=1 calculates 'ass2' using self-consistent equation
c     (iterations are controlled from calling program).
c *** If iopt=2 calculates 'ass2' and 'alfass' using self-consistent
c     equations (iterations are controlled from calling program).
c **********************************************************************
c *** USES:    eshelby     m_effective     invten       voigt
c **********************************************************************
c *** VERSION: 19/may/2002                                           ***
c **********************************************************************
      use mvoigt
      use const
      use mphase_props, only : ngrnph
      use flags, only : ishape,iOutput,nCoatedPh,nCoatingPh,iUmat
      use miscellaneous_sub, only : tmismatch
      use grain_state_v, only : nact
      
      DIMENSION e2i(6,6),ass4(3,3,3,3),anew(6,6)
      DIMENSION aux11(6),aux12(6),aux6(6),aux33(3,3)
      DIMENSION aux21(6,6),aux22(6,6),aux23(6,6),aux24(6,6)
     #         ,aux25(6,6),aux26(6,6),aux27(6,6),aux28(6,6),aux29(6,6)
      DIMENSION esim4(3,3,3,3,NPHM),ESIM4TEMP(3,3,3,3)
     #         ,ESCR4TEMP(3,3,3,3)
      DIMENSION axb(3),EIGB(3,3),aux34(6)
      DIMENSION escr2(6,6),e2inv(6,6),aux31(6),aux32(6)
      DIMENSION E2INVGRTEMP(6,6),EINVSATEMP(3,3,3,3),E2GRTEMP(6,6)
      DIMENSION E2IGRTEMP(6,6)!,ass2(6,6),alfass(6),etrss(6)
      DIMENSION aux13(6),aux14(6),aux15(6),aux16(6) !iEig
      

      real, intent(in) :: axisph(0:3,3,NPHM),axisgr(0:3,3,NGR),css2(6,6)
     #    ,wgt(NGR),acs2(6,6,NGR),alfacs(6,NGR),etrss(6),a_guess(2),
     #    error_mod(2),etrcs_eig(6,NGR),xmmin_in
      real, intent(out) ::einvsa(3,3,3,3,NPHM),einvsagr(3,3,3,3,NGR)
     #         ,escr4(3,3,3,3,NPHM),escr4gr(3,3,3,3,NGR),aef(6,6,NPHM)
     #         ,aefgr(6,6,NGR),meffc(NGR),e2(6,6),aloc2(6,6,NGR)
      real, intent(inout) :: ass2(6,6),alfass(6),etrss_eig(6)
      
      integer, intent (in) :: itmax_mod,ngrain,iopt,liter
     #         ,interaction,istep,ng1,ng2,ilevel,iSkip_ph
     
      integer, intent(inout) :: iflagcond
      
c ______________________________________________________________________

    2 FORMAT(1h ,6d12.4)
    4 FORMAT(1h ,'For iteration = ',i3,' the error is ',d12.4)
c    5 FORMAT(1h+,'ITER: ',i3,'   ERROR: ',d12.4)
    5 FORMAT(1h+,'ITER: ',i3,'   ERROR: ',d12.4,'  NG: ',I5)
    6 FORMAT(1h ,'CONVERGENCE IN SUBROUTINE SC NOT ACHIEVED AFTER ',i3,
     #         ' ITERATIONS',/,
     #1h ,'ABNORMAL PROGRAM STOP......................')
c ______________________________________________________________________

c *** Does the LOOP for the iterative procedure                      ***
c *** Sets counters and number of iterations depending on IOPTION.
      iguess=1
      if (iopt.eq.0) then
        niter=itmax_mod
      else
        niter=1
      endif

      do while (iguess.le.niter)

c *** Calc. Eshelby tensor 'e2' in axis of the ellipsoid -
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++
c *** rotation routine loosly copied from vpsc6

        call voigt(aux6,aux33,ass2,ass4,3)
c        call update_shape !ADDEDMZ commented out
        
c *** ADDEDMZ (comment)    
c     Update of deformation gradient is happening at the end of the deformation step (not in the
c     loop over self-consistent iterations) so update shape will be
c     called many times with same input (deformation gradient) and will
c     give same output. 
c     It should probably be moved after call to subroutine update_fij 
c     but also before first call to subroutine sc_new, since update_shape
c     initializes array AXISPH (and AXISGR if ISHAPE = 2) which contains 
c     info about grain ellipsoids.
c *** END        
        do iph=ngrnph(ng1),ngrnph(ng2)  !MZ_pseudo go over phases and find S(iph),P(iph),Leff(iph) - ishape.lt.2
        if(iph.eq.iskip_ph) cycle !BUG    
        DO J=1,3
          AXB(J)=AXISPH(0,J,iph)
          DO I=1,3
            EIGB(I,J)=AXISPH(I,J,iph)
          ENDDO
        ENDDO

        call stiffness_rotation(ASS4,EIGB,AXB,ESIM4(:,:,:,:,iph)
     #   ,ESCR4(:,:,:,:,iph)) !MZ_pseudo ESIM4 and ESCR4 are now per phase and inverse values as well (E2,E2INV,EINVSA,e2i,aef)

        call voigt(aux6,aux33,e2,esim4(:,:,:,:,iph),4)

        call INVTEN(E2,E2INV)

        CALL VOIGT(AUX6,AUX33,E2INV,EINVSA(:,:,:,:,iph),3)

        IF (ISHAPE.GE.2) THEN
          DO IGR=ng1,ng2
            DO J=1,3
              AXB(J)=AXISGR(0,J,IGR)
              DO I=1,3
                EIGB(I,J)=AXISGR(I,J,IGR)
              END DO
            END DO
            CALL stiffness_rotation(ASS4,EIGB,AXB,ESIM4TEMP,ESCR4TEMP)
            CALL VOIGT(AUX6,AUX33,E2GRTEMP,ESIM4TEMP,4)
            CALL INVTEN(E2GRTEMP,E2INVGRTEMP)
            CALL VOIGT(AUX6,AUX33,E2INVGRTEMP,EINVSATEMP,3)
            DO I=1,3
              DO J=1,3
                DO K=1,3
                  DO L=1,3
                    EINVSAGR(I,J,K,L,IGR)=EINVSATEMP(I,J,K,L)
                    ESCR4GR(I,J,K,L,IGR)=ESCR4TEMP(I,J,K,L)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c          ENDDO  !ADDEDMZ commented out 

          call invten(e2grtemp,e2igrtemp)
c     ADDEDMZ (comment) 
c     Lines CALL INVTEN(E2GRTEMP,E2INVGRTEMP) and call invten(e2grtemp,e2igrtemp)
c     perform same operation so e2igrtemp = E2INVGRTEMP 
c     END          
          do i=1,6
            do j=1,6
              aefgr(i,j,igr)=0.0
              do k=1,6
                aefgr(i,j,igr)=aefgr(i,j,igr)+ass2(i,k)*
     #                           (e2igrtemp(k,j)-id2(k,j))*profac(k)
              enddo
            enddo
          enddo
        ENDDO !ADDEDMZ moved to here
        END IF

c *** Calc. effective stiffnes: aef = ass2 * ( e2**(-1) - I )
        call invten(e2,e2i)
        do i=1,6
          do j=1,6
            aef(i,j,iph)=0.0
            do k=1,6
              aef(i,j,iph)=aef(i,j,iph)+ass2(i,k)*(e2i(k,j)-
     #          id2(k,j))*profac(k)
            enddo
          enddo
        enddo
        enddo !MZ_pseudo end loop over phases

c *** Defines an effective compliance for each grain 'meffc(ng)' for
c     using in the interaction equation as 'meffc*aef'
c     Stiff grain --> less deformation --> 'tangenty' --> meffc < 1
c     Compl grain --> more deformation --> 'secanty'  --> meffc = 1
c *** Forces the secant for the elastic regime (MEFF=1) and widens the
c     interval up to 0.5<MEFF<1.0 as the elastoplastic stiffness decreases
c     Elastic --> ass= css --> tdif=0 --> tmis=0 --> xmmin=1
c     Plastic --> ass<<css --> tdif~tave --> tmis~1 --> xmmin~0.5

        tmis=tmismatch(ass2,css2,6,6)       ! 0<tmis<1
c        xmmin=1.0-0.5*tmis    
        xmmin=1.0                       ! WARNING  --> SECANT IS HARD WIRED !!!!
        if (sum(nact(ng1:ng2)).ne.0) xmmin=xmmin_in
        
        if(interaction.eq.1) then
          do ng=ng1,ng2
            meffc(ng)=xmmin
          enddo
        else
          call m_effective(NG1,NG2,XMMIN,etrss,ass2,wgt,acs2,MEFFC)
c          call m_effective1(NG1,NG2,XMMIN,meffc)
c          write(12,'('' xmmin'',f10.3)') xmmin
c          write(12,'('' meffc'',/,(10f7.2))') (meffc(ng),ng=ng1,ng2,50)
        endif

        !change the interaction tensor with meffc
        if (ishape.eq.0.or.ishape.eq.1) then
          do iph=ngrnph(ng1),ngrnph(ng2)
            do i=1,6
              do j=1,6
                  aef(i,j,iph)=meffc(1)*aef(i,j,iph)
              enddo
            enddo
          enddo
        elseif (ishape.eq.2) then
          do ng=ng1,ng2      ! STARTS LOOP OVER GRAINS
            if(ngrnph(ng).eq.iskip_ph) cycle !BUG  
            if(wgt(ng).ge.0)then !2site (10)
              do i=1,6
                do j=1,6
                  aefgr(i,j,ng)=meffc(ng)*aefgr(i,j,ng)
                enddo
              enddo
            endif
          enddo
        endif
        
c *** Calculates localization and elasto-plastic stiffness tensors:
c
c     aloca = (acs2+aef)**(-1)*(ass2+aef)
c     aloca = (  aux24 )**(-1)*(  aux25 )
c     aloca =       aux21     *   aux25
c     aloca =               aux23
c
c     anew  = < acs2 * aloca > * < aloca >**(-1)
c     anew  = < acs2 * aux23 > * < aux23 >**(-1)
c     anew =        aux22      *   aux26  **(-1)
c     anew =        aux22      *      aux27

!cot_inc for ilevel = 1 and phase = 4 localization tensor for layer needs
!to be used 



        do i=1,6
          aux11(i)=0.0
          aux12(i)=0.0 !iEig
          aux13(i)=0.0
          aux14(i)=0.0
          aux15(i)=0.0
          aux16(i)=0.0    
          do j=1,6
            aux22(i,j)=0.0
            aux26(i,j)=0.0
            aux27(i,j)=0.0
            aux29(i,j)=0.0
          enddo
        enddo

        do ng=ng1,ng2      ! STARTS LOOP OVER GRAINS
          if(ngrnph(ng).eq.iskip_ph) cycle !BUG  
          if(wgt(ng).ge.0)then !2site (10)
            do i=1,6
              do j=1,6
                if (ishape.eq.0.or.ishape.eq.1) then
                  aux24(i,j)=acs2(i,j,ng)+aef(i,j,ngrnph(ng))!meffc(ng)*
                  aux25(i,j)=ass2(i,j) +
     #                  aef(i,j,ngrnph(ng))!meffc(ng)*
                ELSEIF (ishape.ge.2) then
                  aux24(i,j)=acs2(i,j,ng)+aefgr(i,j,ng) !meffc(ng)*
                  aux25(i,j)=ass2(i,j)+aefgr(i,j,ng) !meffc(ng)*
                end if
              enddo
            enddo

            call invten(aux24,aux21)
                
            if (ngrnph(ng).eq.nCoatedPh.or.ngrnph(ng).eq.nCoatingPh)then
              call loc_tens(ng,ass2,axisgr,axisph,acs2,wgt,aux23)  
            else
              do i=1,6
                do j=1,6
                  aux23(i,j)=0.0
                  do k=1,6
                    aux23(i,j)=aux23(i,j)+aux21(i,k)*aux25(k,j)
     #                *profac(k)
                  enddo
                enddo
              enddo 
            endif
            
            !store localization tensor for grain
            aloc2(:,:,ng)=aux23(:,:)
            
            do i=1,6
              do j=1,6
                aux26(i,j)=aux26(i,j)+aux23(i,j)*wgt(ng)
                do k=1,6
                  aux22(i,j)=aux22(i,j)+acs2(i,k,ng)*aux23(k,j)*
     #              profac(k)*wgt(ng)
                enddo
              enddo
            enddo

c ************************************************************************
c    auxiliar tensors required for calculating the overall thermal tensor

            if(iopt.ne.1) then

              do i=1,6
                do j=1,6
                  aux28(i,j)=0.0
                  do k=1,6
                   aux28(i,j)=aux28(i,j)+aux25(i,k)*aux21(k,j)*profac(k)
                  enddo
                enddo
              enddo

              do i=1,6
                aux12(i)=0.0
                do j=1,6
                  do k=1,6
                    aux12(i)=aux12(i)+aux28(i,j)*acs2(j,k,ng)*
     #                alfacs(k,ng)*profac(j)*profac(k)
                  enddo
                enddo
              enddo

              do i=1,6
                aux11(i)=aux11(i)+aux12(i)*wgt(ng)
                do j=1,6
                  aux29(i,j)=aux29(i,j)+aux28(i,j)*wgt(ng)
                enddo
              enddo
c ************************************************************************
c     auxiliar tensors required for calculating the overall eigenstrain (etrss_eig) iEig
c
c     aloca_1 = (acs2+aef)**(-1)*(acs2*etrcs_eig-ass2*etrss_eig)
c     aloca_1 = (  aux24 )**(-1)*(     aux13    -    aux14     )
c     aloca_1 =    aux21*aux13
c     aloca_1 =    aux14
c            
              call TENS_MULT_1(aux13,acs2(:,:,ng),etrcs_eig(:,ng))
              call TENS_MULT_1(aux14,ass2,etrss_eig)
              aux13 = aux13 - aux14
              call TENS_MULT_1(aux14,aux21,aux13)
c
c     etrss_eig = < aloca_1 > - ass2**(-1)*< acs2*(aloca_1 - etrcs_eig) >   
c     etrss_eig = < aux14 >   - ass2**(-1)*< acs2*(aux14 - etrcs_eig) > 
c     etrss_eig = < aux14 >   - ass2**(-1)*< aux13 > 
c     etrss_eig =   aux15     - ass2**(-1)*< aux13 >
c     etrss_eig =   aux15     - anew**(-1)*aux16 
c      
              call TENS_MULT_1(aux13,acs2(:,:,ng),aux14-etrcs_eig(:,ng))
              aux16 = aux16 + aux13*wgt(ng)
              aux15 = aux15 + aux14*wgt(ng)
              
c ************************************************************************
            endif
c ************************************************************************
          endif
        enddo      ! ENDS LOOP OVER GRAINS

        call invten(aux26,aux27)

        do i=1,6
          do j=1,6
            anew(i,j)=0.0
            do k=1,6
              anew(i,j)=anew(i,j)+aux22(i,k)*aux27(k,j)*profac(k)
            enddo
          enddo
        enddo

c *** Enforces the symmetry of the calculated moduli
        do i=1,6
          do j=i+1,6
            anew(i,j)=0.50*(anew(i,j)+anew(j,i))
            anew(j,i)=anew(i,j)
          enddo
        enddo
c
c *** Checks for the convergence of (anew-ass2) to zero

        error=tmismatch(ass2,anew,6,6)
        error_max = error !MZ
        do i=1,6
          do j=1,6
            ass2(i,j)=a_guess(ilevel)*anew(i,j)+
     #                   (1.0-a_guess(ilevel))*ass2(i,j) !MZ
            !ass2_sc(i,j)=anew(i,j)
          enddo
        enddo

        if (iopt.eq.0) then
          if(ilevel.eq.1) then  !2site write only 1st level 
            if(iOutput.eq.1) write(11,4) iguess,error
c            write(*,5)  iguess,error
            if(iUmat.eq.0)then
              write(*,5)  iguess,error,ngrain
            else
              write(5500,5)  iguess,error,ngrain
            endif
          else
            if(iOutput.eq.1) write(72,4) iguess,error
c            write(*,5)  iguess,error
          endif  
        else
          if(ilevel.eq.1) then  !2site write only 1st level 
            if(iOutput.eq.1) write(11,4) liter,error
c            write(*,5)  liter,error
            if(iUmat.eq.0)then
              write(*,5)  liter,error,ngrain
            else
              write(5500,5)  liter,error,ngrain
            endif
          else
            if(iOutput.eq.1) write(72,4) liter,error
          endif
        endif
c       BJCL
        if (iopt.eq.0) then
          if (error.le.error_mod(ilevel)) then
            iguess=itmax_mod+1
            iflagcond=0
          else
            if (iguess.eq.itmax_mod.and.iguess.gt.1) then
              if(ilevel.eq.1) then  !2site write only 1st level 
                if(iOutput.eq.1) write(11,*)
                if(iOutput.eq.1) write(11,6) itmax_mod
                if(iUmat.eq.0)then
                  write(*,*)
                  write(*,6)  itmax_mod
                else
                  write(5500,*)
                  write(5500,6)  itmax_mod
                endif
              else
                if(iOutput.eq.1) write(72,*)
                if(iOutput.eq.1) write(72,6) itmax_mod 
              endif              
              stop
            endif
          iguess=iguess+1
          iflagcond=1
          endif
        else
          if (error.le.error_mod(ilevel).and.liter.gt.1) then       !bjcl has .ge.0   jn has .ge.1   jnflag
            iguess=itmax_mod+1
            iflagcond=0
          else
            if (iguess.eq.itmax_mod) then
              if(ilevel.eq.1) then  !2site write only 1st level                 
                if(iOutput.eq.1) write(11,*)
                if(iOutput.eq.1) write(11,6) itmax_mod
                if(iUmat.eq.0)then
                  write(*,*)
                  write(*,6)  itmax_mod
                else
                  write(5500,*)
                  write(5500,6)  itmax_mod
                endif
              else
                if(iOutput.eq.1) write(72,*)
                if(iOutput.eq.1) write(72,6) itmax_mod
              endif              
              stop
            endif
            iguess=iguess+1
            iflagcond=1
          endif
        endif
c     BJCL
      enddo      ! Closes DO WHILE (IGUESS.LE.NITER)

c ************************************************************************
c *** Evaluates the sc thermal expansion tensor for the sample

      if(iopt.ne.1) then
        call invten(aux29,aux24)
        call invten(anew ,aux22)
        do i=1,6
          alfass(i)=0.0
          do j=1,6
            do k=1,6
              alfass(i)=alfass(i)+aux22(i,j)*aux24(j,k)*aux11(k)
     #                            *profac(j)*profac(k)
            enddo
          enddo
        enddo

c *** Evaluates the sc eigenstrain (etrss_eig)   iEig      
c     etrss_eig =  aux15 - anew**(-1)*aux16 
c     etrss_eig =  aux15 - aux22*aux16 
c     etrss_eig =  aux15 -    aux14
        
        call TENS_MULT_1(aux14,aux22,aux16)
        etrss_eig = aux15 - aux14 
        
        if (iflagcond.eq.0) then
c           write(11,*)
c           write(11,*) 'SELFCONSISTENT MODULI ass2 :'
c           write(11,2) ass2
c           write(11,*) 'THERMAL COEFFICIENTS alfass :'
c           write(11,2) alfass
        endif

      endif
c ************************************************************************

      return
      END SUBROUTINE

      SUBROUTINE stiffness_rotation(STIFFNESS,ROTB,ELIPAXIS,ES4,EAS4)
      
      use miscellaneous_sub, only : ran2
      use meshelby
      
      DIMENSION STIFFNESS(3,3,3,3),ROTB(3,3),ELIPAXIS(3)
      DIMENSION ASS4GA(3,3,3,3),E4GA(3,3,3,3),AUX3333(3,3,3,3)
      DIMENSION ES4(3,3,3,3),EAS4(3,3,3,3)

C *** ROTATION OF STIFFNESS 'C4' FROM SAMPLE TO ELLIPSOID AXES
        DO I=1,3
        DO J=1,3
        DO M=1,3
        DO N=1,3
          DUMMY=0.
          DO I1=1,3
          DO J1=1,3
          DO M1=1,3
          DO N1=1,3
            DUMMY=DUMMY+ROTB(I1,I)*ROTB(J1,J)*ROTB(M1,M)
     #           *ROTB(N1,N)*STIFFNESS(I1,J1,M1,N1)
          END DO
          END DO
          END DO
          END DO
          ass4GA(I,J,M,N)=DUMMY
        END DO
        END DO
        END DO
        END DO

c        CALL ESHELBY (ELIPAXIS,ass4GA,0.,E4GA,AUX3333,1)
          call eshelbyb(ELIPAXIS,ass4GA,0.,E4GA,AUX3333,1)

C *** ROTATES THE ESHELBY TENSOR FOR THE PHASE BACK INTO SAMPLE AXES.
        DO I=1,3
        DO J=1,3
        DO M=1,3
        DO N=1,3
          DUMMYE=0.
          DO I1=1,3
          DO J1=1,3
          DO M1=1,3
          DO N1=1,3
            DUMMYE=DUMMYE+ROTB(I,I1)*ROTB(J,J1)*ROTB(M,M1)
     #           *ROTB(N,N1)*E4GA(I1,J1,M1,N1)
          END DO
          END DO
          END DO
          END DO
          ES4(I,J,M,N)=DUMMYE
        END DO
        END DO
        END DO
        END DO

        DO I=1,3
        DO J=1,3
        DO M=1,3
        DO N=1,3
          DUMMYE=0.
          DO I1=1,3
          DO J1=1,3
          DO M1=1,3
          DO N1=1,3
            DUMMYE=DUMMYE+ROTB(I,I1)*ROTB(J,J1)*ROTB(M,M1)
     #           *ROTB(N,N1)*AUX3333(I1,J1,M1,N1)
          END DO
          END DO
          END DO
          END DO
          EAS4(I,J,M,N)=DUMMYE
        END DO
        END DO
        END DO
        END DO

      END SUBROUTINE
  
      SUBROUTINE M_EFFECTIVE(NG1,NG2,XMMIN,etrss,ass2,wgt,acs2,MEFFC)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SUBROUTINE M_EFFECTIVE  --->   VERSION 26/JAN/2000
C
C     CALCULATES THE PROJECTION OF THE AVERAGE SECANT STIFFNESS, AND
C     EACH INDIVIDUAL GRAIN SECANT STIFFNESS, AGAINST THE IMPOSED
C     STRAIN RATE TENSOR.
C     USES THE RELATIVE VALUE OF THE GRAIN AND OVERALL PROJECTION TO
C     GIVE AN EFFECTIVE VALUE OF 'm' FOR THE INTERACTION EQUATION. SUCH
C     VALUE IS IN THE INTERVAL: 0 < meffc < 1
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      use const
      use mvoigt, only :  profac
      
      integer, intent(in) :: ng1,ng2
      real, intent(in) :: etrss(6),ass2(6,6),wgt(NGR),acs2(6,6,NGR)
      
      real, intent(out) :: MEFFC(NGR)
      
c      use sample_state, only : ass2
c      use grain_props_v, only : wgt
c      use grain_state, only : acs2


      DIMENSION ASS2P(6,6),ACS2P(6,6),PROJ(6,6)
      DIMENSION RC(NGR)

c     write(*,'('' etrss='',6e12.5)') etrss
c     write(*,'('' ass2 ='',6e12.5)') ass2
c     pause

C *** CALCULATES PROJECTOR OVER STRAIN-RATE DIRECTION
      ETRSSNORM=0.
      DO I=1,6
        ETRSSNORM=ETRSSNORM+ETRSS(I)*ETRSS(I)*PROFAC(I)
      ENDDO
        ETRSSNORM=SQRT(ETRSSNORM)
      DO I=1,6
      DO J=1,6
        PROJ(I,J)=ETRSS(I)*ETRSS(J)/ETRSSNORM**2
      ENDDO
      ENDDO

C *** CALCULATES PROJECTION OF OVERALL STIFFNESS AGAINST STRAIN-RATE
      DO I=1,6
      DO J=1,6
        ASS2P(I,J)=0.
        DO K=1,6
          ASS2P(I,J)=ASS2P(I,J)+PROJ(I,K)*ASS2(K,J)*PROFAC(K)
        ENDDO
      ENDDO
      ENDDO

      ASS2PNORM=0.
      DO I=1,6
      DO J=1,6
        ASS2PNORM=ASS2PNORM+ASS2P(I,J)*ASS2P(I,J)*PROFAC(I)*PROFAC(J)
      ENDDO
      ENDDO
        ASS2PNORM=SQRT(ASS2PNORM)

c     write(*,'(''   etrssnorm    ass2pnorm'')')
c     write(*,'(2e12.4)') etrssnorm,ass2pnorm
c     pause

C *** CALCULATES PROJECTION OF GRAIN'S STIFFNESS AGAINST STRAIN-RATE,
C *** PLUS RELATIVE STIFFNESS AND RELATIVE RATE OF EACH GRAIN WITH
C *** RESPECT TO THE OVERALL ONE.

      DO KKK=NG1,NG2

        DO I=1,6
          DO J=1,6
            ACS2P(I,J)=0.
            DO K=1,6
              ACS2P(I,J)=ACS2P(I,J)+PROJ(I,K)*ACS2(K,J,KKK)*PROFAC(K)
            ENDDO
          ENDDO
        ENDDO

        ACS2PNORM=0.
        DO I=1,6
        DO J=1,6
          ACS2PNORM=ACS2PNORM+ACS2P(I,J)*ACS2P(I,J)*PROFAC(I)*PROFAC(J)
        ENDDO
        ENDDO
        ACS2PNORM=SQRT(ACS2PNORM)
        RC(KKK)=ACS2PNORM/ASS2PNORM

      ENDDO

c     write(*,'(''   rc(ng)'')')
c     write(*,'(6e12.4)') (rc(ng),ng=1,ngrain,50)
c     pause

C *** MAKES STATISTIC OF RELATIVE COMPLIANCE OVER ALL THE GRAINS.
C *** ONLY RCMAX & RCMIN ARE REQUIRED BY THE RELAT DIRECT COMPL CRITERION
      RCAVE=0.
      RCDEV=0.
      RCMAX=RC(1)
      RCMIN=RC(1)
      DO KKK=NG1,NG2
        RCAVE=RCAVE+RC(KKK)        *WGT(KKK)
        RCDEV=RCDEV+RC(KKK)*RC(KKK)*WGT(KKK)
        IF(RC(KKK).LT.RCMIN) RCMIN=RC(KKK)
        IF(RC(KKK).GT.RCMAX) RCMAX=RC(KKK)
      ENDDO
      RCDEV=SQRT(RCDEV-RCAVE**2)
      RCTOP=RCAVE+RCDEV
      RCBOT=RCAVE-RCDEV
      IF(RCTOP.GT.RCMAX) RCTOP=RCMAX
      IF(RCBOT.LT.RCMIN) RCBOT=RCMIN
C
C *** (SUPERFLUOUS) CALCULATION OF VOLUME FRACTION WITHIN STANDARD DEV
      VOLFRAC=0.
      DO KKK=NG1,NG2
        IF(RC(KKK).GE.RCBOT .AND. RC(KKK).LE.RCTOP) THEN
          VOLFRAC=VOLFRAC+WGT(KKK)
        ENDIF
      ENDDO
C
C *** DOES A LINEAR INTERPOLATION BETWEEN 1 AND XMMIN
C *** WHEN RC=RCMAX --> MAXIMUM STIFFNESS --> MEFFC = XMMIN  --> TANGENT
C *** WHEN RC=RCMIN --> MINIMUM STIFFNESS --> MEFFC = 1      --> SECANT

      AVMEFF=0.
      DO KKK=NG1,NG2
        XINTERP=(RC(KKK)-RCMIN)/(RCMAX-RCMIN)
        MEFFC(KKK)=1.+(XMMIN-1.)*XINTERP
        AVMEFF=AVMEFF+MEFFC(KKK)*WGT(KKK)
      ENDDO

c      write(12,'(''avmeff, rcmin, rcave, rcmax, rcdev,volfrac'',/,
c     #         6f7.3)') avmeff,rcmin,rcave,rcmax,rcdev,volfrac

      RETURN
      END SUBROUTINE
      
      subroutine m_effective1(ng1,ng2,xmmin,meffc)
      
      use const
      use mvoigt, only :  profac,voigt
      use miscellaneous_sub, only : det
      use grain_state_v, only : stcs,nact
      use sample_state_v, only : stss
      
      integer, intent(in) :: ng1,ng2
      real, intent(in) :: xmmin
      
      real, intent(out) :: MEFFC(NGR)
      
      DIMENSION sthyd(3,3),aux66(6,6),aux3333(3,3,3,3),st33(3,3)
      DIMENSION stdev(3,3),hard_r(NGR)
      
      !find stress ratio (i.e. hardness ratio)
      
      !find hydrostatic and deviatoric stresses
      call VOIGT(stss,st33,aux66,aux3333,1)
      aux=( st33(1,1)+st33(2,2)+st33(3,3) )/3.0
      sthyd=0.0
      do i=1,3
        sthyd(i,i)=aux
      enddo
      stdev=st33-sthyd
      !find invariants J2 and J3 of dev stress
      aJ2=0.0
      do i=1,3
        do j=1,3
          aJ2=aJ2+0.5*stdev(i,j)*stdev(i,j)
        enddo
      enddo
      aJ3=det(stdev)
      steqss=sqrt(aJ2)*sqrt(3.0)
      
      hard_rmin=1.0
      hard_rmax=0.0
      do ng=ng1,ng2
        if (steqss.ne.0) then
          !find hydrostatic and deviatoric stresses
          call VOIGT(stcs(:,ng),st33,aux66,aux3333,1)
          aux=( st33(1,1)+st33(2,2)+st33(3,3) )/3.0
          sthyd=0.0
          do i=1,3
            sthyd(i,i)=aux
          enddo
          stdev=st33-sthyd
          !find invariants J2 and J3 of dev stress
          aJ2=0.0
          do i=1,3
            do j=1,3
              aJ2=aJ2+0.5*stdev(i,j)*stdev(i,j)
            enddo
          enddo
          aJ3=det(stdev)
          steq=sqrt(aJ2)*sqrt(3.0)
          hard_r(ng)=steq/steqss
        else
          hard_r(ng)=1.0
        endif     
        if(hard_r(ng).lt.hard_rmin) hard_rmin=hard_r(ng)
        if(hard_r(ng).gt.hard_rmax) hard_rmax=hard_r(ng)
      enddo
      
      xmmax=1.0
      DO ng=NG1,NG2
        if (steqss.ne.0) then
          XINTERP=(hard_r(ng)-hard_rmin)/(hard_rmax-hard_rmin)
          MEFFC(ng)=xmmax+(XMMIN-xmmax)*XINTERP
        else
          MEFFC(ng)=1.
        endif
      ENDDO
      
      end subroutine      

      SUBROUTINE av_modulus(iopt,ng1,ng2,wgt,scs2,ccs2,ass2,css2,sss2)
c **********************************************************************
c *** Calc. the VOIGT, REUSS and HILL averages for the elastic mod.  ***
c **********************************************************************
c *** USES:    invten                                                ***
c **********************************************************************
c *** VERSION: 12/OCT/94                                             ***
c **********************************************************************
c
      use const
      use mvoigt
      use flags, only : iOutput

      integer, intent(in) :: iopt,ng1,ng2
      real, intent(in) :: wgt(NGR),scs2(6,6,NGR),ccs2(6,6,NGR)
      real, intent(out) :: ass2(6,6),css2(6,6),sss2(6,6)
      
c ______________________________________________________________________
c
      DIMENSION css2vo(6,6),css2re(6,6)
      DIMENSION sss2vo(6,6),sss2re(6,6)
c ______________________________________________________________________
c
    1 FORMAT(1h ,'ELASTIC PROPERTIES AVERAGE:')
    2 FORMAT(1h ,6d12.4)
    3 FORMAT(1h ,78('*'))
c ______________________________________________________________________
c
      do i=1,6
        do j=1,6
          css2vo(i,j)=0.0
          sss2re(i,j)=0.0
        enddo
      enddo
c
      do ng=ng1,ng2
        do i=1,6
          do j=1,6
            css2vo(i,j)=css2vo(i,j)+wgt(ng)*ccs2(i,j,ng)
            sss2re(i,j)=sss2re(i,j)+wgt(ng)*scs2(i,j,ng)
          enddo
        enddo
      enddo
c
      call invten(css2vo,sss2vo)
      call invten(sss2re,css2re)      
c
      if(iopt.eq.1)then
c
        if(iOutput.eq.1)then
          write(11,1)
          write(11,*)
          write(11,*) ' VOIGT average stiffness matrix'
          write(11,2) ((css2vo(i,j),j=1,6),i=1,6)
          write(11,*)
          write(11,*) ' VOIGT average compliance matrix'
          write(11,2) ((sss2vo(i,j),j=1,6),i=1,6)
          write(11,*)
          write(11,*) ' REUSS average stiffness matrix'
          write(11,2) ((css2re(i,j),j=1,6),i=1,6)
          write(11,*)
          write(11,*) ' REUSS average compliance matrix'
          write(11,2) ((sss2re(i,j),j=1,6),i=1,6)
        endif
      endif
      do i=1,6
        do j=1,6
          css2(i,j)=(css2vo(i,j)+css2re(i,j))/2.0
          ass2(i,j)=css2(i,j)
        enddo
      enddo
c
      if(iopt.eq.1)then
        if(iOutput.eq.1)then
          write(11,*)
          write(11,*)
          write(11,*) ' HILL average stiffness matrix'
          write(11,2) ((css2(i,j),j=1,6),i=1,6)
        endif
        call invten(css2,sss2)   
        if(iOutput.eq.1)then
          write(11,*)
          write(11,*) ' HILL average compliance matrix'
          write(11,2) ((sss2(i,j),j=1,6),i=1,6)
          write(11,*)
          write(11,3)
        endif
      endif
c ______________________________________________________________________
c
      return
      END SUBROUTINE

      subroutine loc_tens(ng,ass2,axisgr,axisph,acs2,wgt,aloc)
      
      use miscellaneous_sub, only : svdcmp
      use mphase_props_v, only : ngrnph,nphngr
      use const
      use mvoigt
      use flags
      
      DIMENSION aux6(6),aux7(6),aux8(6),aux9(6),aux33(3,3)
   
      DIMENSION aux21(6,6),aux22(6,6),aux23(6,6),aux24(6,6),aux25(6,6)
     #         ,aux26(6,6),aux27(6,6),aux28(6,6),aux29(6,6)
      
      DIMENSION aux3333(3,3,3,3),ass4(3,3,3,3),acs4(3,3,3,3)
      
      DIMENSION v(6,6),W(6)
      
      DIMENSION AXB(3),EIGB(3,3),ESIM4_eff(3,3,3,3)
     #         ,ESIM4_l(3,3,3,3),e2_eff(6,6),e2_l(6,6)
     #         ,T_eff(6,6),T_2(6,6),aloc_1(6,6),aloc_2(6,6)
     #         ,w_2_1(6,6),alfa_1(6,6),alfa_2(6,6),aloc_I(6,6)
      
      integer, intent(in) :: ng
      real, intent(in) :: axisph(0:3,3,NPHM),axisgr(0:3,3,NGR)
     #    ,wgt(NGR),acs2(6,6,NGR),ass2(6,6)
      real, intent(out) :: aloc(6,6)

      if (ngrnph(ng).ne.nCoatedPh.and.ngrnph(ng).ne.nCoatingPh) then
        ! loc tensor for inclusion without coating   
        ! S(Leff); T(Leff)   
        iph=ngrnph(ng)
        call voigt(aux6,aux33,ass2,ass4,3)
        DO J=1,3
          AXB(J)=AXISPH(0,J,iph)
          DO I=1,3
            EIGB(I,J)=AXISPH(I,J,iph)
          ENDDO
        ENDDO
        
        call stiffness_rotation(ASS4,EIGB,AXB,ESIM4_eff(:,:,:,:)
     #   ,aux3333) 
        
        call voigt(aux6,aux33,e2_eff,esim4_eff(:,:,:,:),4)
        
        call INVTEN(ass2,aux29)  
        call TENS_MULT(T_eff,e2_eff,aux29)       
        ! A = [ I + S(Leff)*Leff**(-1)*(L1 - Leff) ]**(-1)
        ! A = [ I + T_eff*aux29 ]**(-1)
        ! A = [ I + aux28 ]**(-1)
        aux29 = acs2(:,:,ng) - ass2
        call TENS_MULT(aux28,T_eff,aux29)
        aux28 = id2 + aux28
        call INVTEN(aux28,aloc)
        
        ! L* = Leff*(S(Leff)**(-1) - I)
        ! L* = Leff*(aux29 - I)
        ! L* = Leff*aux29
        ! L* = aux28
        call INVTEN(e2_eff,aux29)
        aux29 = aux29 - id2
        call TENS_MULT(aux28,ass2,aux29)
        
        ! A = (Lc + L*)**(-1)*(L* + Leff) 
        ! A = (Lc + aux28)**(-1)*(aux28 + Leff) 
        ! A = aux29**(1)*aux27
        ! A = aux26*aux27
        aux29 = acs2(:,:,ng) + aux28
        aux27 = aux28 + ass2
        call INVTEN(aux29,aux26)
        call TENS_MULT(aloc,aux26,aux27)
        
        return
      endif
      
      ! particle (p) = ng_1
      ! layer (l) = ng_2      
      if (ngrnph(ng).eq.nCoatedPh) then
        ngCtds=SUM(nphngr(1:nCoatedPh-1))
        ng_1=ng
        ng_1_rel=ng_1-ngCtds
        ngCtgs=SUM(nphngr(1:nCoatingPh-1))
        ng_2=ngCtgs+ng_1_rel
      elseif (ngrnph(ng).eq.nCoatingPh) then
        ngCtgs=SUM(nphngr(1:nCoatingPh-1))
        ng_2=ng
        ng_2_rel=ng_2-ngCtgs
        ngCtds=SUM(nphngr(1:nCoatedPh-1))
        ng_1=ngCtds+ng_2_rel
      endif
      
      ! step 1    
      ! S = S(Leff); S = S(Ll)
      ! T(Leff) = S(Leff) * Leff**(-1); T(Ll) = S(Ll) * Ll**(-1)
      
      ! S(Leff); T(Leff)
      call voigt(aux6,aux33,ass2,ass4,3)
      DO J=1,3
        AXB(J)=AXISPH(0,J,nCoatingPh)
        DO I=1,3
          EIGB(I,J)=AXISPH(I,J,nCoatingPh)
        ENDDO
      ENDDO

      call stiffness_rotation(ASS4,EIGB,AXB,ESIM4_eff(:,:,:,:)
     # ,aux3333) 

      call voigt(aux6,aux33,e2_eff,esim4_eff(:,:,:,:),4)
      
      call INVTEN(ass2,aux29)  
      call TENS_MULT(T_eff,e2_eff,aux29)
      
      ! S(Ll,ng_2); T(Ll,ng_2)
      call voigt(aux6,aux33,acs2(:,:,ng_2),acs4,3)
      DO J=1,3
        AXB(J)=AXISPH(0,J,nCoatedPh)
        DO I=1,3
          EIGB(I,J)=AXISPH(I,J,nCoatedPh)
        ENDDO
      ENDDO

      call stiffness_rotation(acs4,EIGB,AXB,ESIM4_l(:,:,:,:)
     # ,aux3333) 

      call voigt(aux6,aux33,e2_l,esim4_l(:,:,:,:),4)     
      
      !regular inversion
      call INVTEN(acs2(:,:,ng_2),aux29)  
      !SVD inversion
c      aux21=acs2(:,:,ng_2)
c      call svdcmp(aux21,6,6,6,6,w,v)
c      ! condition number
c      cond_number=maxval(w(1:6))/minval(w(1:6))
c      wmax=0. !Will be the maximum singular value obtained.
c      do j=1,6
c          if(w(j).gt.wmax)wmax=w(j)
c      enddo 
c      wmin=wmax*1.0e-12 !This is where we set the threshold for singular values
c                       !allowed to be nonzero. The constant is typical,but not universal.           
c      aux29=0.0
c      do i=1,6
c          if (w(i).gt.wmin) then
c              aux29(i,i)=1.0/w(i)
c          else
c              aux29(i,i)=0.0
c          endif
c      enddo
c      aux29=matmul(aux29,transpose(aux21))
c      aux29=matmul(v,aux29)    
c      do i=1,6
c          do j=1,6
c              aux29(i,j)=aux29(i,j)*INVFAC(i,j)
c          enddo
c      enddo
      
      call TENS_MULT(T_2,e2_l,aux29)
      
      if (wgt(ng_2).eq.0.0) T_2 = 0.0
      
      ! step 2
      ! alfa_1 = (I*phi_1 + phi_2*w_2_1)**(-1)
      ! alfa_2 = w_2_1*alfa_1
      ! w_2_1 = I - T_2*(L_2-L_1)
      phi_1 = wgt(ng_1)/(wgt(ng_1)+wgt(ng_2))
      phi_2 = wgt(ng_2)/(wgt(ng_1)+wgt(ng_2))
      
      if ((wgt(ng_1)+wgt(ng_2)).eq.0.0) then
          phi_1 = 1.0
          phi_2 = 0.0
      endif
      
      aux29 = acs2(:,:,ng_2) - acs2(:,:,ng_1)
      call TENS_MULT(aux28,T_2,aux29)
      w_2_1 = id2 -aux28
      
      aux27 = (id2*phi_1 + phi_2*w_2_1)
      call INVTEN(aux27,alfa_1)
      call TENS_MULT(alfa_2,w_2_1,alfa_1)
      
      ! step 4
      ! A_I = [ I + T_eff*[phi_1*(L_1 - L_eff)*alfa_1 + phi_2*(L_2 - L_eff)*alfa_2] ]**(-1)
      ! A_I = [ I + T_eff*[aux28 + aux27] ]**(-1)
      ! A_I = [ I + T_eff*aux29 ]**(-1)
      ! A_I = [ I + aux26 ]**(-1)
      ! A_I = aux26**(-1)
      aux29 = acs2(:,:,ng_1) - ass2
      call TENS_MULT(aux28,aux29,alfa_1)
      aux28 = phi_1*aux28
      
      aux29 = acs2(:,:,ng_2) - ass2
      call TENS_MULT(aux27,aux29,alfa_2)
      aux27 = phi_2*aux27
      
      aux29 = aux28 + aux27
      call TENS_MULT(aux26,T_eff,aux29)
      aux26 = id2 + aux26
      
      call INVTEN(aux26,aloc_I)
      
      ! step 6
      ! A_1 = alfa_1*A_I
      ! A_2 = alfa_2*A_I
      call TENS_MULT(aloc_1,alfa_1,aloc_I)
      call TENS_MULT(aloc_2,alfa_2,aloc_I)
      
      ! return loc tensor
      if (ngrnph(ng).eq.nCoatedPh) aloc = aloc_1
      if (ngrnph(ng).eq.nCoatingPh) aloc = aloc_2
      
      end subroutine
      
      end module sc_estimate

c      
c***********************************************************************
c *** diagonal *********************************************************
c***********************************************************************
      module diagonal


      contains
      
      SUBROUTINE hd_diagonal(ng,gamd,nact,iact,hd,ioptbs)
      
      use const
      use mvoigt, only : id2,profac,i6
      use mphase_props, only : ngrnph,nslsys
      use hard_law1_v, only : iOpposite, alpha_reg
      use back_stress_v, only : alpha_reg_bs
      
      real, intent(in) :: gamd(NSLS,NGR)
      integer, intent(in) :: ng,nact(NGR),iact(NSLS,NGR),ioptbs
      real, intent(out) :: hd(NSLS,NSLS)
      
      DIMENSION hd_diag(NSLS),aux(NSLS),aux1(NSLS),hd_off(NSLS,NSLS)
      
      iph=ngrnph(ng)
      
      if (ioptbs.eq.1) then
          alpha_reg_tmp=alpha_reg_bs
      elseif (ioptbs.eq.0) then
          alpha_reg_tmp=alpha_reg
      endif
      
      if (sum(gamd(1:nslsys(iph),ng)).ne.0.0) then
        hd_off=hd
        do ns1=1,nslsys(iph)
          hd_off(ns1,ns1)=0.0
          hd_off(ns1,iOpposite(ns1))=0.0
        enddo
        aux = matmul( hd(1:nslsys(iph),1:nslsys(iph))
     #       ,gamd(1:nslsys(iph),ng) )
        aux1 = matmul( alpha_reg_tmp*hd_off(1:nslsys(iph),1:nslsys(iph))
     #       ,gamd(1:nslsys(iph),ng) )
        hd = 0.0
        do ns1=1,nact(ng)
            n1=iact(ns1,ng)
            hd(n1,n1) = (aux(n1)-aux1(n1))/gamd(n1,ng)
            if (hd(n1,n1).lt.1.0) hd(n1,n1)=1.0
            hd(iOpposite(n1),n1) = hd(n1,n1)
        enddo
        hd=hd+alpha_reg_tmp*hd_off
      else
        do ns1=1,nact(ng)
            n1=iact(ns1,ng)      
            aux(n1)=sum(hd(n1,:))
        enddo
        hd=0.0
        do ns1=1,nact(ng)
            n1=iact(ns1,ng)
            hd(n1,n1) = aux(n1)
            hd(iOpposite(n1),n1) = hd(n1,n1)
        enddo
      endif
      
      RETURN
      END SUBROUTINE

      SUBROUTINE hd_to_modulus(ng,hd,hd_bs,gamd,nact,iact,ccs2,mcs,nmcs,
     #    stcs,f,aux65)
      
      use const
      use mphase_props, only : ngrnph
      use miscellaneous_sub, only : lubksbc,ludcmpc,ran2
      use mvoigt, only : id2,profac,i6
      use back_stress, only : bs_hd
      use flags, only : iBackStress,inonSch,iLatBS
      
      PARAMETER (TOLER_DET=1.0e-20)
      
      real, intent(in) :: gamd(NSLS,NGR),ccs2(6,6,NGR),mcs(6,NSLS,NGR)
     #    ,nmcs(6,NSLS,NGR),stcs(6,NGR),hd(NSLS,NSLS),hd_bs(NSLS,NSLS,2)
      integer, intent(in) :: ng
      real, intent(out) :: aux65,f(6,NSLS,NGR)
      integer, intent(inout) :: iact(NSLS,NGR),nact(NGR)
      
      DIMENSION aux21(6,6),aux65(6,6),aux6(6),x(NSLS,NSLS),y(NSLS,NSLS)
      DIMENSION indx(NSLS),dgamma(NSLS),rss(NSLS)
      
      !find f matrix
      do ns1=1,nact(ng)
        n1=iact(ns1,ng)
      do ns2=1,nact(ng)
        n2=iact(ns2,ng)
        x(ns1,ns2)= hd(n1,n2)
        if(iBackStress.eq.1) x(ns1,ns2)=x(ns1,ns2)+ hd_bs(n1,n2,iLatBS)!MZ_bs added last term
        do i=1,6
          do j=1,6
            x(ns1,ns2)=x(ns1,ns2)+(mcs(i,n1,ng)+nmcs(i,n1,ng))
     #             *mcs(j,n2,ng)*ccs2(i,j,ng)*profac(i)*profac(j) !inonSch
          enddo
        enddo
      enddo
      enddo

      !force x (and in turn y) to be symmetric
c      do ns1=1,nact(ng)
c        do ns2=1,ns1
c          x(ns1,ns2)=0.5*(x(ns1,ns2)+x(ns2,ns1))
c          x(ns2,ns1)=x(ns1,ns2)
c        enddo
c      enddo
      
      call ludcmpc(x,nact(ng),NSLS,indx,d)
      do ns1=1,nact(ng)
        d=d*x(ns1,ns1)
      enddo

      if (abs(d).lt.TOLER_DET) then
        idelsys=(ran2(jran)*nact(ng)+1)
        if(idelsys.lt.nact(ng)) then     ! why not .LE. instead?
          do ns1=idelsys,nact(ng)
            iact(ns1,ng)=iact(ns1+1,ng)
          enddo
        endif
        nact(ng)=nact(ng)-1
        igverify=0
        write(*,*) ng,' DET=0 --> look for other combination'
      else
        do ns1=1,nact(ng)
          do ns2=1,nact(ng)
            y(ns1,ns2)=(ns1/ns2)*(ns2/ns1)
          enddo
        enddo
        do ns1=1,nact(ng)
          call lubksbc(x,nact(ng),NSLS,indx,y(1,ns1))
        enddo
        do ns1=1,nact(ng)
          do j=1,6
            f(j,ns1,ng)=0.0
            do ns2=1,nact(ng)
              n2=iact(ns2,ng)
              do i=1,6
                f(j,ns1,ng)=f(j,ns1,ng)+y(ns1,ns2)*(mcs(i,n2,ng)+
     #                   nmcs(i,n2,ng))*profac(i)*(ccs2(i,j,ng) !inonSch
     #                    -stcs(i,ng)*i6(j))                        !jn updated F
              enddo
            enddo
          enddo
        enddo
      endif
      
      !find modulus
      do i=1,6
        do j=1,6
          aux21(i,j)=id2(i,j)
          do ns1=1,nact(ng)
            n1=iact(ns1,ng)
            aux21(i,j)=aux21(i,j)-mcs(i,n1,ng)*f(j,ns1,ng)
          enddo
        enddo
      enddo
      do i=1,6
        do j=1,6
          aux65(i,j)=0.0
          do k=1,6
            aux65(i,j)=aux65(i,j)+ccs2(i,k,ng)*aux21(k,j)
     #                   *profac(k)
          enddo
          aux65(i,j)=aux65(i,j)-stcs(i,ng)*i6(j)              !jn added this line
        enddo
      enddo

      
      RETURN
      END SUBROUTINE
      
      end module diagonal
c      
c***********************************************************************
c *** solve_t **********************************************************
c***********************************************************************
      module solve_t
      
      use const
      use sc_estimate
      
      private
      public main_calc,inc_toler
      
      contains
      !g_actsys     :Tests the potentially active systems   
      !g_modulus    :Calc. the el-pl moduli for each grain  
      !g_state      :Calc. the stress & strain rate in grain
      !g_verify     :Checks the active loading condition    
      !s_state      :Calculates the sample stress & strain rate
      
      SUBROUTINE main_calc(iproc,istep,ng1,ng2,idiv ,iact,nact
     #                ,jran,itmax_mod,ngrain,ng_update,nout_old,ietbc
     #                ,istbc,etrss,strss,alfass,ass2,tau,tau_bcst,gamtot
     #                ,gamd,ccs2,mcs,nmcs,stcs,alfacs,deltemp,e2,aef
     #                ,aefgr,acs2,aloc2,taud,taud_bcst,etrcs,strcs
     #                ,axisph,axisgr,css2,wgt,a_guess,error_mod,einvsa
     #                ,einvsagr,escr4,escr4gr,meffc,tau_update,etrbc
     #                ,strbc,etrcs_eig,etrss_eig,xmmin_in) 
      
      use const
      use flags, only : ishape, ipileup,iSingleCry
      use mphase_props, only : nsys,ngrnph
      use bc, only : itmax_grain
      
      !local
      integer :: iSkip_ph,iflagcond,ilevel,interaction,iopt,liter
      real :: auxsample(6,NPHM)
     
      !in
      integer, intent(in) :: iproc,istep,ng1,ng2,itmax_mod,ngrain
     #    ,ietbc(6),istbc(6)
      real, intent(in) :: a_guess(2),error_mod(2),css2(6,6),etrbc(6)
     #    ,strbc(6),deltemp,axisph(0:3,3,NPHM),ccs2(6,6,NGR)
     #    ,mcs(6,NSLS,NGR),nmcs(6,NSLS,NGR),alfacs(6,NGR),wgt(NGR)
     #    ,stcs(6,NGR),gamtot(NGR),axisgr(0:3,3,NGR)
     #    ,tau_bcst(NSLS,NGR,2),etrcs_eig(6,NGR),xmmin_in
     
      !out
      integer, intent(out) :: idiv,iact(NSLS,NGR),nact(NGR)
     #    ,ng_update(NGR)
      real, intent(out) :: etrss(6),strss(6),gamd(NSLS,NGR)
     #    ,taud(NSLS,NGR),taud_bcst(NSLS,NGR,2),aloc2(6,6,NGR)
     #    ,e2(6,6),einvsagr(3,3,3,3,NGR),escr4(3,3,3,3,NPHM)
     #    ,escr4gr(3,3,3,3,NGR),meffc(NGR),einvsa(3,3,3,3,NPHM)
     #    ,etrcs(6,NGR),strcs(6,NGR),tau_update(NSLS,NGR),acs2(6,6,NGR)
      !inout
      integer, intent(inout) ::jran,nout_old
      real, intent(inout) :: alfass(6),tau(NSLS,NGR),ass2(6,6)
     #    ,aef(6,6,NPHM),aefgr(6,6,NGR),etrss_eig(6)
      
    9 FORMAT(1h ,'ABNORMAL PROGRAM STOP',/,1h ,'DOES NOT CONVERGE
     # AFTER ',i3,' ITERATIONS OVER GRAINS SYSTEMS')  
      
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c *** Does a loop to find the sets of active systems in all grains that
c *** satisfy the overall problem.
c *** Stops if convergence is not achieved before 'itmax_grain' iterations
          iSkip_ph=0    

          it_grain=0
          iverify=1
          
          !clean gamd array
          do ng=1,ngrain
            do ns=1,nsys(ngrnph(ng))
              gamd(ns,ng)=0.0
            enddo     
          enddo
		  
          do while (iverify.ne.0)
c
            it_grain=it_grain+1
c            write(11,8) it_grain
            if (it_grain.eq.itmax_grain) then
c              write(11,*)
c              write(11,9)
              write(*,*)
              write(*,9)
c              stop !split
              idiv=1
              exit ! exit the loop
            endif
c
c *** Given the grain stress 'stcs' identifies and tags the potentially
c *** active systems in each grain
            call g_actsys(nout_old,ng1,ng2,iSkip_ph,
     #    iact,nact,ng_update,tau_update,
     #    mcs,nmcs,stcs,tau,tau_bcst)

c *** Given the elasto-plastic tensor 'ass2' and the boundary conditions
c *** calculates strain rate and stress rate components for the sample
            call s_state(ietbc,etrbc,istbc,strbc,
     #    alfass,ass2,aef,deltemp,
     #    etrss,strss,auxsample,etrss_eig) !finds auxsample
            
c *** Given the potentially active systems in each grain, checks if the
c *** on-the-SCYS condition is fulfilled. Modifies active systems if not.
c *** Calculates grain magnitudes "acs2','strcs','etrcs'
            call g_modulus(ng1,ng2,etrss,strss,alfass,
     #    ass2,iSkip_ph,ccs2,mcs,nmcs,stcs,aef,
     #    aefgr,alfacs,deltemp,tau,gamtot,
     #    acs2,gamd,auxsample,etrcs,strcs,taud,taud_bcst,
     #    jran,iact,nact,etrcs_eig,etrss_eig) !uses auxsample
            
c *** Solves for the self-consistent elasto-plastic tensor 'ass2'.
            interaction=1
c            iopt=1
            iopt=2
            if(ipileup.eq.1) iopt=2   ! recalculates virtual thermal expansion

            if(iSingleCry.eq.1) then !CPFE
                ass2=acs2(:,:,1)
                iverify=0
                call s_state(ietbc,etrbc,istbc,strbc,
     #          alfass,ass2,aef,deltemp,
     #          etrss,strss,auxsample,etrss_eig) !finds auxsample
                
                !01.25.2016. strcs is found in the same way as strss
                strcs(:,1)=strss
                etrcs(:,1)=etrss
            else 
                call sc_new(iopt,it_grain,iverify,e2,interaction,istep,
     #    ng1,ng2,etrss,alfass,ass2,1,iSkip_ph,
     #    axisgr,axisph,css2,einvsa,einvsagr,escr4,escr4gr,aef,
     #    aefgr,acs2,alfacs,aloc2,itmax_mod,ngrain,wgt,a_guess,
     #    error_mod,meffc,etrcs_eig,etrss_eig,xmmin_in)
            endif
        enddo      ! end of (iverify) loop
        
        !clean gamd array
        do ng=1,ngrain
          do ns=1,nsys(ngrnph(ng))
            if (gamd(ns,ng).lt.0.0) then
              gamd(ns,ng)=0.0
            endif
          enddo     
        enddo
        
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      END SUBROUTINE

      SUBROUTINE g_modulus(ng1,ng2,etrss,strss,alfass,
     #    ass2,iSkip_ph,ccs2,mcs,nmcs,stcs,aef,
     #    aefgr,alfacs,deltemp,tau,gamtot,
     #    acs2,gamd,auxsample,etrcs,strcs,taud,taud_bcst,
     #    jran,iact,nact,etrcs_eig,etrss_eig)
c **********************************************************************
c *** Calculates the single crystal elasto-plastic incremental       ***
c *** stiffness 'acs2' for all the grains for each subset of active  ***
c *** loading systems.
c *** Checks positiveness of the shear rates and the active loading  ***
c *** condition.                                                     ***
c **********************************************************************
c *** USES:    g_verify    lubksbc    ludcmpc                        ***
c **********************************************************************
c *** VERSION: 12/OCT/94                                             ***
c **********************************************************************

      use miscellaneous_sub, only : lubksbc,ludcmpc,ran2,svdcmp
      use hard_law1, only : hl_hd1=>hl_hd
      use hard_law2, only : hl_hd2=>hl_hd
      use const
      use flags, only : idiag,kCL,iSingleCry,iBackStress,iLatBS,iDiagBs
      use mvoigt, only : id2,profac,i6
      use back_stress, only : bs_hd
      use mphase_props, only : ngrnph
      use diagonal

      PARAMETER (TOLER_DET=1.0e-20)

      DIMENSION aux21(6,6),x(NSLS,NSLS),y(NSLS,NSLS),indx(NSLS)
      DIMENSION y_t(NSLS,NSLS),f(6,NSLS,NGR),hd(NSLS,NSLS)
      DIMENSION hd_bs(NSLS,NSLS,2),w(NSLS),v(NSLS,NSLS)
      
      !in
      integer, intent(in) :: ng1,ng2,iSkip_ph
      real, intent(in) :: etrss(6),strss(6),alfass(6)
     #     ,ass2(6,6),ccs2(6,6,NGR),mcs(6,NSLS,NGR),nmcs(6,NSLS,NGR)
     #     ,stcs(6,NGR),aef(6,6,NPHM),aefgr(6,6,NGR)
     #     ,alfacs(6,NGR),deltemp,tau(NSLS,NGR)
     #     ,gamtot(NGR),etrcs_eig(6,NGR),etrss_eig(6)
      !out
      real, intent(out) :: acs2(6,6,NGR),gamd(NSLS,NGR),taud(NSLS,NGR)
     #     ,auxsample(6,NPHM),etrcs(6,NGR),strcs(6,NGR)
     #     ,taud_bcst(NSLS,NGR,2)
      
      !inout
      integer, intent(inout) :: iact(NSLS,NGR),nact(NGR),jran
c ______________________________________________________________________

      do ng=ng1,ng2
        if(ngrnph(ng).eq.iskip_ph) cycle !BUG
        igverify=0
        if (nact(ng).eq.0) then
          do i=1,6
            do j=1,6
              acs2(i,j,ng)=ccs2(i,j,ng)-stcs(i,ng)*i6(j)   !jnflag: added term
            enddo
          enddo
          igverify=1
          call g_state(ng,etrss,strss,alfass,ass2
     #                ,acs2,aef,aefgr,alfacs,etrcs,deltemp
     #                ,strcs,auxsample,etrcs_eig,etrss_eig)
        else
          do while (igverify.eq.0.and.nact(ng).gt.0)
          
            if (kCl.eq.1) call hl_hd1(ng,hd)
            if (kCl.eq.2) call hl_hd2(ng,hd)
            if (iBackStress.eq.1) call bs_hd(ng,hd_bs)

            do ns1=1,nact(ng)
              n1=iact(ns1,ng)
            do ns2=1,nact(ng)
              n2=iact(ns2,ng)
              x(ns1,ns2)= hd(n1,n2)
              if(iBackStress.eq.1) x(ns1,ns2)=x(ns1,ns2)+
     #            hd_bs(n1,n2,iLatBS)
              do i=1,6
                do j=1,6
                  x(ns1,ns2)=x(ns1,ns2)+(mcs(i,n1,ng)+nmcs(i,n1,ng))
     #                   *mcs(j,n2,ng)*ccs2(i,j,ng)*profac(i)*profac(j)!inonSch
                enddo
              enddo
            enddo
            enddo

            call ludcmpc(x,nact(ng),NSLS,indx,d)
            do ns1=1,nact(ng)
              d=d*x(ns1,ns1)
            enddo

            if (abs(d).lt.TOLER_DET) then
              idelsys=(ran2(jran)*nact(ng)+1)
              if(idelsys.lt.nact(ng)) then     ! why not .LE. instead?
                do ns1=idelsys,nact(ng)
                  iact(ns1,ng)=iact(ns1+1,ng)
                enddo
              endif
              nact(ng)=nact(ng)-1
              igverify=0
              write(*,*) ng,' DET=0 --> look for other combination'
            else
              do ns1=1,nact(ng)
                do ns2=1,nact(ng)
                  y(ns1,ns2)=(ns1/ns2)*(ns2/ns1)
                enddo
              enddo
              do ns1=1,nact(ng)
                call lubksbc(x,nact(ng),NSLS,indx,y(1,ns1))
              enddo
              
            !SVD inversion
c            if (.true.) then
c              n=nact(ng)
c              call svdcmp(x,n,n,NSLS,NSLS,w,v)
c              ! condition number
cc              cond_number=maxval(w(1:n))/minval(w(1:n))
c              
c              wmax=0. !Will be the maximum singular value obtained.
c              do j=1,n
c                  if(w(j).gt.wmax)wmax=w(j)
c              enddo 
c              wmin=wmax*1.0e-12 !This is where we set the threshold for singular values
c                               !allowed to be nonzero. The constant is typical,but not universal.           
c              y=0.0
c              do i=1,n
c                  if (w(i).gt.wmin) then
c                      y(i,i)=1.0/w(i)
c                  else
c                      y(i,i)=0.0
c                  endif
c              enddo
c              y(1:n,1:n)=matmul(y(1:n,1:n)
c     #          ,transpose(x(1:n,1:n)))
c              y(1:n,1:n)=matmul(v(1:n,1:n),y(1:n,1:n))    
              
              do ns1=1,nact(ng)
                do j=1,6
                  f(j,ns1,ng)=0.0
                  do ns2=1,nact(ng)
                    n2=iact(ns2,ng)
                    do i=1,6
                      f(j,ns1,ng)=f(j,ns1,ng)+y(ns1,ns2)*(mcs(i,n2,ng)+
     #                          nmcs(i,n2,ng))*profac(i)*(ccs2(i,j,ng)!inonSch
     #                          -stcs(i,ng)*i6(j))                     !jn updated F
                    enddo
                  enddo
                enddo
              enddo

              do i=1,6
                do j=1,6
                  aux21(i,j)=id2(i,j)
                  do ns1=1,nact(ng)
                    n1=iact(ns1,ng)
                    aux21(i,j)=aux21(i,j)-mcs(i,n1,ng)*f(j,ns1,ng)
                  enddo
                enddo
              enddo
              do i=1,6
                do j=1,6
                  acs2(i,j,ng)=0.0
                  do k=1,6
                    acs2(i,j,ng)=acs2(i,j,ng)+ccs2(i,k,ng)*aux21(k,j)
     #                           *profac(k)
                  enddo
                  acs2(i,j,ng)=acs2(i,j,ng)-stcs(i,ng)*i6(j)              !jn added this line
                enddo
              enddo

c *** Enforces the symmetry of the calc. moduli                      ***
c              do i=1,6
c                do j=i+1,6
c                  acs2(i,j,ng)=0.50*(acs2(i,j,ng)+acs2(j,i,ng))
c                  acs2(j,i,ng)=acs2(i,j,ng)
c                enddo
c              enddo

              call g_state(ng,etrss,strss,alfass,ass2
     #                    ,acs2,aef,aefgr,alfacs,etrcs,deltemp
     #                    ,strcs,auxsample,etrcs_eig,etrss_eig)
              call g_verify(ng,alfacs,deltemp,f,etrcs,strcs,
     #                tau,mcs,nmcs,gamtot,
     #                igverify,taud,taud_bcst,gamd,
     #                iact,nact,jran,etrcs_eig)
            endif
          enddo

          if (nact(ng).eq.0) then
c           write(*,*) ng,' !!! WARNING: reaches zero system state'
            do i=1,6
              do j=1,6
                acs2(i,j,ng)=ccs2(i,j,ng)-stcs(i,ng)*i6(j)          !jnflag new term on end
              enddo
            enddo
            igverify=1
            call g_state(ng,etrss,strss,alfass,ass2
     #                    ,acs2,aef,aefgr,alfacs,etrcs,deltemp
     #                    ,strcs,auxsample,etrcs_eig,etrss_eig)
          endif
        
          !iDiag - diagonalize hd and find corrected acs2
          if (iDiag.eq.1) then
            !diagonalize the hd matrix (hd is input)
            call hd_diagonal(ng,gamd,nact,iact,hd,0)
          endif
          if (iDiagBs.eq.1) then
            !diagonalize the hd_bs matrix (hd_bs is input)
            call hd_diagonal(ng,gamd,nact,iact,hd_bs(:,:,2),1)
          endif
          if (iDiag.eq.1.or.IDiagBs.eq.1) then
            !find modulus and additional stress with given hd (hd is input)
            call hd_to_modulus(ng,hd,hd_bs,gamd,nact,iact,ccs2,mcs,nmcs
     #          ,stcs,f,acs2(:,:,ng))
          endif
        endif        
      enddo

c ______________________________________________________________________
c
      return
      end subroutine

      SUBROUTINE G_VERIFY (ng,alfacs,deltemp,f,etrcs,strcs,
     #                tau,mcs,nmcs,gamtot,
     #                igverify,taud,taud_bcst,gamd,
     #                iact,nact,jran,etrcs_eig)
      
c **********************************************************************
c *** Tests if the condition of active-loading is fulfilled for the
c *** tentative set in routine G_MODULUS. If it is verified (igverify=
c *** =1), if not then (igverify=0) and reduce the NACT number.
c **********************************************************************
c *** USES:    g_state                                               ***
c **********************************************************************
c *** VERSION: 13/FEB/98                                             ***
c **********************************************************************
      use miscellaneous_sub, only : ran2
      use hard_law1, only : hl_hd1=>hl_hd
      use hard_law2, only : hl_hd2=>hl_hd
      use const
      use flags, only : idiag,kCL,iBackStress,iLatBS
      use mphase_props
      use mvoigt, only : profac  
      use back_stress, only : bs_hd

      PARAMETER (ERROR_LOAD=0.001) ! Tolerance for active loading condition
      dimension gammamod(NMOD),hd(NSLS,NSLS),hd_bs(NSLS,NSLS,2)
      
      !in
      integer, intent(in) :: ng
      real, intent(in) :: alfacs(6,NGR),deltemp,f(6,NSLS,NGR)
     #    ,etrcs(6,NGR),strcs(6,NGR),tau(NSLS,NGR),mcs(6,NSLS,NGR)
     #    ,nmcs(6,NSLS,NGR),gamtot(NGR),etrcs_eig(6,NGR)
      !out
      integer, intent(out) :: igverify
      real, intent(out) :: taud(NSLS,NGR),taud_bcst(NSLS,NGR,2)
     #    ,gamd(NSLS,NGR)
      
      !inout
      integer, intent(inout) :: iact(NSLS,NGR),nact(NGR),jran

c ______________________________________________________________________
c
      igverify=1
      if (nact(ng).ne.0) then
        do ns1=1,nsys(ngrnph(ng))
          gamd(ns1,ng)=0.0
        enddo
C        write(11,*) nact(ng),ng
        do ns1=1,nact(ng)
          n1=iact(ns1,ng)
          do i=1,6
            gamd(n1,ng)=gamd(n1,ng)+f(i,ns1,ng)*(etrcs(i,ng)-
     #                  alfacs(i,ng)*deltemp-etrcs_eig(i,ng))*profac(i)!iEig
          enddo
        enddo

c********************************************************************
        ns1=1
        do while(ns1.le.nact(ng))
          n1=iact(ns1,ng)
          if (gamd(n1,ng).le.0.0) then
c            write(*,*)  ng,' Look for other combination GAMMA<0'
c            write(12,*) ng,' Look for other combination GAMMA<0'
c            write(12,*) 'in system',n1
            igverify=0
            nact(ng)=nact(ng)-1
            do ns2=ns1,nact(ng)
              iact(ns2,ng)=iact(ns2+1,ng)
            enddo
            ns1=nact(ng)
          endif
          ns1=ns1+1
        enddo
        if (igverify.ne.0) then
          neload=0

c           do ns1=1,nsys
c             voce=thet1(ns1)
c             if(tau1(ns1).gt.0.001*tau0(ns1)) then
c               thet0x=thet0(ns1)
c               thet1x=thet1(ns1)
c               fact  =gamtot(ng)*thet0x/tau1(ns1)
c               voce  =voce+(thet0x-thet1x+thet1x*fact)*exp(-fact)
c             endif
c             taud(ns1,ng)=0.0
c             do ns2=1,nact(ng)
c               n2=iact(ns2,ng)
c               taud(ns1,ng)=taud(ns1,ng)+voce*h(ns1,n2)*gamd(n2,ng)
c             enddo
c            enddo
c         else

          if (kCl.eq.1) call hl_hd1(ng,hd)
          if (kCl.eq.2) call hl_hd2(ng,hd)
          
          do ns1=1,nsys(ngrnph(ng))
            taud(ns1,ng)=0.0
            do ns2=1,nact(ng)
              n2=iact(ns2,ng)
              taud(ns1,ng)=taud(ns1,ng)+hd(ns1,n2)*gamd(n2,ng)
            enddo
          enddo
          
          if (ibackstress.eq.1) then
            call bs_hd(ng,hd_bs)
            do ns1=1,nsys(ngrnph(ng))
              taud_bcst(ns1,ng,1:2)=0.0
              do ns2=1,nact(ng)
                n2=iact(ns2,ng)
                taud_bcst(ns1,ng,1)=taud_bcst(ns1,ng,1)+
     #           hd_bs(ns1,n2,1)*gamd(n2,ng)
                taud_bcst(ns1,ng,2)=taud_bcst(ns1,ng,2)+
     #           hd_bs(ns1,n2,2)*gamd(n2,ng)
              enddo
            enddo    
          endif
         
          do ns1=1,nact(ng)
            n1=iact(ns1,ng)
            rssd=0.0
            do i=1,6
              rssd=rssd+(mcs(i,n1,ng)+nmcs(i,n1,ng))*strcs(i,ng)
     #             *profac(i)
            enddo
            
            if (iBackstress.eq.1) rssd=rssd-taud_bcst(n1,ng,iLatBS)
            
            control_load=abs((rssd-taud(n1,ng))/taud(n1,ng))
c            cl2=ABS((rssd1-taud(n1,ng))/taud(n1,ng))
            if (control_load.gt.ERROR_LOAD) then
              neload=neload+1
c              write(12,*)
c              write(12,*) 'NG: ',ng,'     RSSD>TAUD in SYSTEM: ',n1
c              write (*,*)  'NG: ',ng,'     RSSD>TAUD in SYSTEM: ',n1
c     #                   ,control_load,cl2,taud(n1,ng),rssd,rssd1
c              read (*,*)
            endif
          enddo
          if (neload.ne.0) then
c            write(*,*) ng,' Look for other combination RSSD>TAUD'
            idelsys=(ran2(jran)*nact(ng)+1)
            try=ran2(jran)   !MZ_bs
c              idelsys=1
c              idelsys=nact(ng)
c              write(*,*) 'ng, nact, idelsys',ng,nact(ng),idelsys
            nact(ng)=nact(ng)-1
            if(idelsys.gt.0.and.idelsys.lt.nact(ng)) then
                do ns1=idelsys,nact(ng)
                  iact(ns1,ng)=iact(ns1+1,ng)
                enddo
              endif
            igverify=0
          endif
        endif
      endif
c ______________________________________________________________________
c
      return
      END SUBROUTINE

      SUBROUTINE G_STATE (ng,etrss,strss,alfass,ass2,acs2,aef,aefgr
     #                    ,alfacs,etrcs,deltemp,strcs,auxsample
     #                    ,etrcs_eig,etrss_eig)
c **********************************************************************
c *** Evaluates the stress rate and strain rate in each grain:
c       etrcs=(Lc+L~)^(-1):[(Lbar+L~):etrss-(Lbar:alfas-Lc:alfac)*delT]
c       auxsample=(Lbar+L~):etrss-Lbar:alfas*delT] --> calc ins S_STATE
c       strcs= Lc:(etrcs-alfac*delT)
c **********************************************************************
c *** USES:    invten                                                ***
c **********************************************************************
c *** VERSION: 12/OCT/94                                             ***
c **********************************************************************
      use miscellaneous_sub, only : lubksbc,ludcmpc
      use const
      use flags, only : ishape,iSingleCry,nCoatingPh,nCoatedPh
      use mvoigt
      use mphase_props, only : ngrnph
      use grain_props_v, only : wgt
      use grain_state_v, only : axisgr
      use mphase_state_v, only : axisph
      use sc_estimate, only : loc_tens

      DIMENSION aux11(6),aux21(6,6),aux22(6,6),aux6(6)
      DIMENSION aloc(6,6) !cot_inc
      
      !in
      integer, intent(in) :: ng
      real, intent(in) :: acs2(6,6,NGR),aef(6,6,NPHM),aefgr(6,6,NGR)
     #    ,alfacs(6,NGR),deltemp,etrss(6),strss(6)
     #    ,alfass(6),ass2(6,6),etrcs_eig(6,NGR),etrss_eig(6)
     
      !out
      real, intent(out) :: strcs(6,NGR),etrcs(6,NGR),auxsample(6,NPHM)
c ______________________________________________________________________
c
c
      if (iSingleCry.ne.1.and.nCoatedPh+nCoatingPh.eq.0) then
      do i=1,6
        do j=1,6
          IF (ishape.eq.0.or.ishape.eq.1) THEN
            aux21(i,j)=aef(i,j,ngrnph(ng))+acs2(i,j,ng)
          ELSEIF (ISHAPE.GE.2) THEN
            aux21(i,j)=aefgr(i,j,ng)+acs2(i,j,ng)
c *** ADDEDMZ calculation of auxsample for each grain.
c     auxsample is in function of effective stiffness and should be 
c     different for each grain when ISHAPE = 2. auxsample is calculated
c     in s_state with average effective stiffness. 
c     Here part of the code from s_state that calculates auxsample is 
c     coppied but using the effective stiffness of each grain instead.            
          do i_mz=1,6
            aux11(i_mz)=0.0
            auxsample(i_mz,ngrnph(ng))=0.0
            do j_mz=1,6
              aux11(i_mz)=aux11(i_mz)+(aefgr(i_mz,j_mz,ng)+
     #        ass2(i_mz,j_mz))*etrss(j_mz)*profac(j_mz)
              auxsample(i_mz,ngrnph(ng))=auxsample(i_mz,ngrnph(ng))+
     #        ass2(i_mz,j_mz)*alfass(j_mz)*deltemp*
     #        profac(j_mz)
     #        +ass2(i_mz,j_mz)*etrss_eig(j_mz)* !iEig
     #        profac(j_mz) !iEig
            enddo
            auxsample(i_mz,ngrnph(ng))=aux11(i_mz)-
     #        auxsample(i_mz,ngrnph(ng))
          enddo       
c     END                
          ENDIF
        enddo
      enddo

      call invten(aux21,aux22)
      do i=1,6
        aux11(i)=0.0
        do j=1,6
          aux11(i)=aux11(i)+acs2(i,j,ng)*
     #             alfacs(j,ng)*deltemp*profac(j)
     #             +acs2(i,j,ng)*etrcs_eig(j,ng)*profac(j) !iEig
        enddo
      enddo

      do i=1,6
        etrcs(i,ng)=0.0
        do j=1,6
          etrcs(i,ng)=etrcs(i,ng)+aux22(i,j)*profac(j)*
     #                (auxsample(j,ngrnph(ng))+aux11(j))
        enddo
      enddo

      do i=1,6
        strcs(i,ng)=0.0
        do j=1,6
          strcs(i,ng)=strcs(i,ng)+acs2(i,j,ng)*(etrcs(j,ng)-
     #                alfacs(j,ng)*deltemp
     #                -etrcs_eig(j,ng))*profac(j) !iEig
        enddo
      enddo
      
      elseif(nCoatedPh+nCoatingPh.ne.0)then
      ! find loc tensors
      ! et_c = A_c*E
      ! st_c = L_c*(et_c - et_eig_c) 
      call loc_tens(ng,ass2,axisgr,axisph,acs2,wgt,aloc) !in here appropriate ass2 is used depending on the phase
      call TENS_MULT_1(etrcs(:,ng),aloc,etrss) 
      call TENS_MULT_1(strcs(:,ng),acs2(:,:,ng),etrcs(:,ng))         
      
      elseif(iSingleCry.eq.1)then !CPFE
        do i=1,6
          etrcs(i,ng)=0.0
          do j=1,6
            etrcs(i,ng)=etrss(i)
          enddo
        enddo
        
        do i=1,6
          strcs(i,ng)=0.0
          do j=1,6
            strcs(i,ng)=strcs(i,ng)+acs2(i,j,ng)*(etrcs(j,ng)-
     #                  alfacs(j,ng)*deltemp)
          enddo
        enddo          
      
      endif  
      
      RETURN
      END SUBROUTINE

      SUBROUTINE g_actsys(nout_old,ng1,ng2,iSkip_ph,
     #    iact,nact,ng_update,tau_update,
     #    mcs,nmcs,stcs,tau,tau_bcst)
c **********************************************************************
c *** Identifies the active systems in each grain                    ***
c *** Modify the critical stress if the stress state is out of SCYS  ***
c **********************************************************************
c *** VERSION: 12/OCT/94                                             ***
c **********************************************************************
c
      use const
      use mvoigt, only : profac
      use mphase_props, only : nsys,ngrnph,iTwinSys,nslsys
      use flags, only : iBackStress,ilatBS,iOutput
      use twinning, only : TwFrSy
      
      ! in
      integer, intent(in) :: ng1,ng2,iSkip_ph
      real, intent(in) :: mcs(6,NSLS,NGR),stcs(6,NGR)
     #    ,tau_bcst(NSLS,NGR,2),nmcs(6,NSLS,NGR)

      !out
      integer, intent(out) :: iact(NSLS,NGR),nact(NGR),ng_update(NGR)

      !inout
      real, intent(inout) :: tau_update(NSLS,NGR),tau(NSLS,NGR)
      integer, intent(inout) :: nout_old
      
      
      DIMENSION rss(nsls)
      REAL    delt
      LOGICAL log
c ______________________________________________________________________
c  
      slack=0.98
      nout=0
      do ng=ng1,ng2 
        if(ngrnph(ng).eq.iskip_ph) cycle !BUG
        iflag=0
        delt=0.0
        nout_flag=0
        nact(ng)=0
        do ns1=1,nsys(ngrnph(ng))
          rss(ns1)=0.0
          if(iTwinSys(ns1).eq.1) then
            log=(iTwinSys(ns1).eq.1.and
     #        .TwFrSy(ns1-nslsys(ngrnph(ng)),ng).gt.0.0)
     #        .or.(iTwinSys(ns1).eq.1.and.mod(ns1,2).eq.1)
          endif
          if (iTwinSys(ns1).eq.0.or.log) then
            do i=1,6 
              rss(ns1)=rss(ns1)+(mcs(i,ns1,ng)+nmcs(i,ns1,ng))
     #                 *stcs(i,ng)*profac(i)  !inonSch
            enddo     
            if(iBackStress.eq.1) rss(ns1)=rss(ns1)-tau_bcst(ns1,ng
     #         ,ilatBS)
            rss(ns1)=rss(ns1)/tau(ns1,ng)       
            if (rss(ns1).ge.slack) then
              nact(ng)=nact(ng)+1
              iact(nact(ng),ng)=ns1
              if (rss(ns1).gt.1.00) then
                nout_flag=1
                tau(ns1,ng)=rss(ns1)*tau(ns1,ng) !DEBUG_2site
c                if(ng.eq.24.and.ns1.eq.2)
c     #           write(190,'(1000f25.8)')rss(ns1)!DEBUG_2site
                delt=delt+(rss(ns1)-1)**2
                  iflag=1
c                tau(ns1,ng)=rss(ns1)*tau(ns1,ng)
c                tau_update(ns1,ng)=rss(ns1)*tau_update(ns1,ng)
                ng_update(ng)=1
              endif
            endif
          endif
        enddo
        delt=(delt)**(0.5)
        tau_update(1,ng)=tau_update(1,ng)*(delt+1.0)

        if (nout_flag.eq.1) nout=nout+1
      enddo

      if (nout.gt.nout_old) then
        nout_old=nout
        if(iOutput.eq.1)then
          write(12,*)
          write(12,*) 'WARNING'
          write(12,*) 'STRESS IS OUT OF THE SCYS FOR ',nout,' GRAINS'
        endif
      endif
c ______________________________________________________________________
c
      return
      end subroutine
      
      SUBROUTINE s_state(ietbc,etrbc,istbc,strbc,
     #    alfass,ass2,aef,deltemp,
     #    etrss,strss,auxsample,etrss_eig)
c **********************************************************************
c *** Evaluates the stress and strain rate in sample using:
c     STRSS(i)=ASS2(i,j)*(ETRBC(j)-ALFASS(j)*DELTATEMP)
c *** Shuffles know and unknow components of stress and strain rate
c **********************************************************************
c *** USES:    ludcmpc     lubksbc                                   ***
c **********************************************************************
c *** VERSION: 12/OCT/94                                             ***
c **********************************************************************
      use miscellaneous_sub, only : lubksbc,ludcmpc,cond_numb
      use const
      use sample_props_v, only : nph
      use mvoigt, only : profac
      
      !in
      integer, intent(in) :: ietbc(6),istbc(6)
      real, intent(in) :: etrbc(6),strbc(6),ass2(6,6)
     #    ,alfass(6),aef(6,6,NPHM),deltemp,etrss_eig(6)
      
      !out
      real, intent(out) ::  etrss(6),strss(6),auxsample(6,NPHM)
      

      DIMENSION aux11(6),aux21(6,6),indx(6)

c ______________________________________________________________________
c
      do i=1,6
       aux11(i)=-1.0*istbc(i)*strbc(i)
        do j=1,6
          aux11(i)=aux11(i)+ass2(i,j)*(ietbc(j) !2site
     #                     *etrbc(j)-alfass(j)*deltemp
     #                     -etrss_eig(j))*profac(j) !iEig
          aux21(i,j)=ietbc(j)*(i/j)*(j/i)-
     #               istbc(j)*ass2(i,j)*profac(j)
        enddo
      enddo

c      call cond_numb(6,6,aux21,cond)
c      if (cond.gt.1e7) pause
      call ludcmpc(aux21,6,6,indx,d)
      call lubksbc(aux21,6,6,indx,aux11)
      do i=1,6
        etrss(i)=ietbc(i)*etrbc(i)+istbc(i)*aux11(i) !2site
        strss(i)=istbc(i)*strbc(i)+ietbc(i)*aux11(i)
      enddo

c     write(11,*)
c     write(11,*) 'Polycrystal Stress Rate State:'
c     write(11,'(1h ,6d12.4)') strss
c     write(11,*) 'Polycrystal Strain Rate State:'
c     write(11,'(1h ,6d12.4)') etrss

c *** Calculate the auxiliar vector 'AUXSAMPLE' for using in G_STATE.
      do iph=1,nph !MZ_pseudo go over phases
      do i=1,6
        aux11(i)=0.0
        auxsample(i,iph)=0.0
        do j=1,6
          aux11(i)=aux11(i)+(aef(i,j,iph)+ass2(i,j))
     #      *etrss(j)*profac(j)
          auxsample(i,iph)=auxsample(i,iph)+
     #      ass2(i,j)*alfass(j)*deltemp*profac(j)
     #      +ass2(i,j)*etrss_eig(j)*profac(j) !iEig
        enddo
        auxsample(i,iph)=aux11(i)-auxsample(i,iph)
      enddo
      enddo !MZ_pseudo end loop over phases
c______________________________________________________________________
c
      return
      END SUBROUTINE

      SUBROUTINE inc_toler(ilevel,i_converge,error_mod,error_mod_org,
     # ass2,ass2_old,a_guess)
      
      !in
      integer, intent(in) :: ilevel
      real, intent(in) :: error_mod_org(2),ass2_old(6,6)
     
      !out
      integer, intent(out) :: i_converge
      real, intent(out) :: error_mod(2),ass2(6,6)
      
      !inout
      real, intent(inout) :: a_guess(2)
        
      if (ilevel.eq.1) then
      ! use original guess
      ass2 = ass2_old  
      ! incease tolerance
        if (error_mod(ilevel).eq.error_mod_org(ilevel)) then
            error_mod(ilevel) = error_mod_org(ilevel)*10.0
        elseif (a_guess(ilevel).eq.2.0/3.0) then
            !pause
            a_guess(ilevel) = 9.0/10.0
            error_mod(ilevel) = error_mod_org(ilevel)*10.0      
        elseif (a_guess(ilevel).eq.9.0/10.0) then
            !pause
            error_mod(ilevel) = error_mod_org(ilevel)*20.0  
        elseif (error_mod(ilevel).eq.error_mod_org(ilevel)*20.0) then
            !pause
            error_mod(ilevel) = error_mod_org(ilevel)*50.0              
        endif
      endif
      if (ilevel.eq.2) then
        if (error_mod(ilevel).eq.error_mod_org(ilevel)) then
            error_mod(ilevel) = error_mod_org(ilevel)*10.0
        elseif (error_mod(ilevel).eq.error_mod_org(ilevel)*10.0) then
            !pause
            error_mod(ilevel) = error_mod_org(ilevel)*50.0
        elseif (error_mod(ilevel).ge.error_mod_org(ilevel)*50.0) then  
            error_mod(ilevel) = error_mod(ilevel) + 
     #        error_mod_org(ilevel)*20.0            
            if(error_mod(ilevel).gt.error_mod_org(ilevel)*100.0) 
     #        i_converge=1
            write(72,*) 'TOLERANCE INCREASED ABOVE *50'
        endif  
      endif
      
      RETURN
      END SUBROUTINE
      
      end module solve_t
c      
c***********************************************************************
c *** state_var ********************************************************
c***********************************************************************
      module state_var

      use flags
      use const
      use hard_law1, only : state_hl1=>state_hl
      use hard_law2, only : state_hl2=>state_hl
      use back_stress, only : state_bs
      use twinning, only : state_tw
      use phase_transf, only : state_pt
      
      type grain
          real :: st(6),et(6),fij(3,3),axor(0:3,3),r(3,3),wgt,ngrnph
          type(state_hl1) :: hl1
          type(state_hl2) :: hl2
          type(state_bs) :: bs
          type(state_tw) :: tw
          type(state_pt) :: pt
      end type grain
      
      type phase
          integer :: ngrain
          real :: st(6),et(6),fij(3,3),axor(0:3,3),euler(3),aef(6,6)
          type(grain), dimension(NGR) :: grain
      end type phase
      
      type sample
          integer :: i_prev_proc,nph,ngrain,ngParent,jran
          real :: st(6),et(6),etrav(6),ass2(6,6),alfass(6),etrss_eig(6)
          type (phase), dimension(NPHM) :: phase
      end type sample

      contains
      
      SUBROUTINE set(s)
      
      use sample_props_v
      use mphase_props
      use grain_props_v
      
      use sample_state
      use mphase_state
      use grain_state
      
      use sample_rate_v
      use mphase_rate
      use grain_rate
      
      use hard_law1, only : set_hl1=>set_hl
      use hard_law2, only : set_hl2=>set_hl
      use back_stress, only : set_bs
      use twinning, only : set_tw
      use phase_transf, only : set_pt
      
      type (sample), intent(out) :: s
      
      s%i_prev_proc = i_prev_proc
      s%nph = nph
      s%ngrain = ngrain
      s%ngParent = ngParent
      s%jran = jran
      s%st = stss
      s%et = etss
c      s%etrav = 
      s%ass2 = ass2
      s%alfass = alfass
      if (iPhTr.eq.1) s%etrss_eig = etrss_eig
      do iph=1,nph
          s%phase(iph)%ngrain = nphngr(iph)
          s%phase(iph)%fij = fijph(:,:,iph)
          s%phase(iph)%axor = axisph(:,:,iph)
          s%phase(iph)%euler = eulerph(:,iph)
          s%phase(iph)%aef = aef(:,:,iph)
          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2=SUM(nphngr(1:iph))
          ng_ph=0
          do ng=ng1,ng2
              ng_ph=ng_ph+1
              s%phase(iph)%grain(ng_ph)%st = stcs(:,ng)
              s%phase(iph)%grain(ng_ph)%et = etcs(:,ng)
              s%phase(iph)%grain(ng_ph)%r = r(:,:,ng)
              s%phase(iph)%grain(ng_ph)%wgt = wgt(ng)
              s%phase(iph)%grain(ng_ph)%ngrnph = ngrnph(ng)
              !hard law
              if(kCL.eq.1)call set_hl1(ng,s%phase(iph)%grain(ng_ph)%hl1)
              if(kCL.eq.2)call set_hl2(ng,s%phase(iph)%grain(ng_ph)%hl2)
              if(iBackStress.eq.1) then
                  call set_bs(ng,s%phase(iph)%grain(ng_ph)%bs)
              endif
              if(itwinning.eq.1) then
                  call set_tw(ng,s%phase(iph)%grain(ng_ph)%tw)
              endif
              if(iPhTr.eq.1)then
                  call set_pt(ng,s%phase(iph)%grain(ng_ph)%pt)
              endif
          enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE get(s)
      
      use sample_props_v
      use mphase_props
      use grain_props_v
      
      use sample_state
      use mphase_state
      use grain_state
      
      use sample_rate_v
      use mphase_rate
      use grain_rate
      
      use hard_law1, only : get_hl1=>get_hl
      use hard_law2, only : get_hl2=>get_hl
      use back_stress, only : get_bs
      use twinning, only : get_tw
      use phase_transf, only : get_pt
      
      type (sample), intent(in) :: s
      
      if (iUmat.eq.1) i_prev_proc = s%i_prev_proc !umat only
      nph = s%nph
      ngrain = s%ngrain
      ngParent = s%ngParent
      jran = s%jran
      stss = s%st
      etss = s%et
c     = s%etrav  
      ass2 = s%ass2
      alfass = s%alfass
      if (iPhTr.eq.1) etrss_eig = s%etrss_eig
      do iph=1,nph
          nphngr(iph) = s%phase(iph)%ngrain
          fijph(:,:,iph) = s%phase(iph)%fij
          axisph(:,:,iph) = s%phase(iph)%axor
          eulerph(:,iph) = s%phase(iph)%euler
          aef(:,:,iph) = s%phase(iph)%aef
          ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2=SUM(nphngr(1:iph))
          ng_ph=0
          do ng=ng1,ng2
              ng_ph=ng_ph+1
              stcs(:,ng) = s%phase(iph)%grain(ng_ph)%st
              etcs(:,ng) = s%phase(iph)%grain(ng_ph)%et
              r(:,:,ng) = s%phase(iph)%grain(ng_ph)%r
              wgt(ng) = s%phase(iph)%grain(ng_ph)%wgt
              ngrnph(ng) = s%phase(iph)%grain(ng_ph)%ngrnph
              !hard law
              if(kCl.eq.1)call get_hl1(ng,s%phase(iph)%grain(ng_ph)%hl1)
              if(kCl.eq.2)call get_hl2(ng,s%phase(iph)%grain(ng_ph)%hl2)
              if(iBackStress.eq.1) then
                  call get_bs(ng,s%phase(iph)%grain(ng_ph)%bs)
              endif
              if(itwinning.eq.1) then
                  call get_tw(ng,s%phase(iph)%grain(ng_ph)%tw)
              endif
              if(iPhTr.eq.1)then
                  call get_pt(ng,s%phase(iph)%grain(ng_ph)%pt)
              endif
          enddo
      enddo
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE read_statv(s,fileprev,NSTATV,STATEV)
      
      use mphase_props
      use hard_law1, only : read_statv_hl1=>read_statv_hl
      use hard_law2, only : read_statv_hl2=>read_statv_hl
      use back_stress, only : read_statv_bs
      use twinning, only : read_statv_tw
      use phase_transf, only : read_statv_pt
      
      real :: STATEV(NSTATV),field(6,6)  
      logical :: mask(6,6)
      character(len=150), intent(in) :: fileprev
      type (sample), intent(out) :: s  
    
      !read
      if (iUmat.eq.0) then
          iunit=1
          open(unit=iunit,file=fileprev,status='unknown')    
          read(iunit,*)  !title
          do i=1,NSTATV
              read(1,*,IOSTAT=IOstatus) STATEV(i)
              if (IOstatus.lt.0) exit
          enddo
          close(iunit)
      endif
      
      !restore state from statev array     
      mask = .true.
      ns=1
      s%i_prev_proc = STATEV(ns)
      ns=ns+1
      !if it is 0 STATEV array don't reset the total number
      if (iUmat.eq.1.and.s%i_prev_proc.ne.0.or.iUmat.eq.0) then 
          s%nph = STATEV(ns)
      endif
      ns=ns+1
      !if it is 0 STATEV array don't reset the total number
      if (iUmat.eq.1.and.s%i_prev_proc.ne.0.or.iUmat.eq.0) then 
          s%ngrain = STATEV(ns)
      endif
      ns=ns+1
      !if it is 0 STATEV array don't reset the total number
      if (iUmat.eq.1.and.s%i_prev_proc.ne.0.or.iUmat.eq.0) then 
          s%ngParent = STATEV(ns)
      endif
      ns=ns+1
      s%jran = STATEV(ns)
      ns=ns+1
      s%st = STATEV(ns:ns+5)
      ns=ns+6
      s%et = STATEV(ns:ns+5)
c     = s%etrav  
      ns=ns+6
      s%ass2 = unpack(STATEV(ns:ns+35),mask,field)
      ns=ns+36
      s%alfass = STATEV(ns:ns+5)
      ns=ns+6
      ng=0
      do iph=1,s%nph
          !if it is 0 STATEV array don't reset the total number
          if (iUmat.eq.1.and.s%i_prev_proc.ne.0.or.iUmat.eq.0) then 
              s%phase(iph)%ngrain = STATEV(ns)
          endif
          ns=ns+1
          s%phase(iph)%fij = unpack(STATEV(ns:ns+8),mask(1:3,1:3)
     #        ,field(1:3,1:3))
          ns=ns+9
          s%phase(iph)%axor = unpack(STATEV(ns:ns+11),mask(1:4,1:3)
     #        ,field(1:4,1:3))
          ns=ns+12
          s%phase(iph)%euler = STATEV(ns:ns+2)
          ns=ns+3
          s%phase(iph)%aef = unpack(STATEV(ns:ns+35),mask,field)
          ns=ns+36
          ns=ns+6
          do ng_ph=1,s%phase(iph)%ngrain
              ng=ng+1
              s%phase(iph)%grain(ng_ph)%st = STATEV(ns:ns+5)
              ns=ns+6
              s%phase(iph)%grain(ng_ph)%et = STATEV(ns:ns+5)
              ns=ns+6
              s%phase(iph)%grain(ng_ph)%r = unpack(STATEV(ns:ns+8)
     #            ,mask(1:3,1:3),field(1:3,1:3))
              ns=ns+9
              s%phase(iph)%grain(ng_ph)%wgt = STATEV(ns)
              ns=ns+1
              !if it is 0 STATEV array don't reset the total number
              if (iUmat.eq.1.and.s%i_prev_proc.ne.0.or.iUmat.eq.0) then 
                  s%phase(iph)%grain(ng_ph)%ngrnph = STATEV(ns)
              endif
              ns=ns+1
              if(kCL.eq.1)call read_statv_hl1(
     #            s%phase(iph)%grain(ng_ph)%hl1,ng,iph,ns,NSTATV,STATEV)
              if(kCL.eq.2)call read_statv_hl2(
     #            s%phase(iph)%grain(ng_ph)%hl2,ng,iph,ns,NSTATV,STATEV)
              if (iBackStress.eq.1) then
                  call read_statv_bs(s%phase(iph)%grain(ng_ph)%bs,ng,iph
     #            ,ns,NSTATV,STATEV)
              endif
              if (itwinning.eq.1) then
                  call read_statv_tw(s%phase(iph)%grain(ng_ph)%tw,ng,iph
     #            ,ns,NSTATV,STATEV)
              endif
              if (iPhTr.eq.1) then
                  call read_statv_pt(s%phase(iph)%grain(ng_ph)%pt,ng,iph
     #            ,ns,NSTATV,STATEV)
              endif
          enddo
      enddo
      
      end subroutine read_statv
      
      SUBROUTINE write_statv(s,fileprev,NSTATV,STATEV)
      
      use mphase_props
      use hard_law1, only : write_statv_hl1=>write_statv_hl
      use hard_law2, only : write_statv_hl2=>write_statv_hl
      use back_stress, only : write_statv_bs
      use twinning, only : write_statv_tw
      use phase_transf, only : write_statv_pt
      
      real :: STATEV(NSTATV)
      character(len=150), intent(in) :: fileprev
      type (sample), intent(in) :: s
      
    1 FORMAT(1h ,1d37.30)
      
      !store to statv array      
      ns=1
      STATEV(ns) = s%i_prev_proc
      ns=ns+1
      STATEV(ns) = s%nph
      ns=ns+1
      STATEV(ns) = s%ngrain
      ns=ns+1
      STATEV(ns) = s%ngParent
      ns=ns+1
      STATEV(ns) = s%jran
      ns=ns+1
      STATEV(ns:ns+5) = s%st
      ns=ns+6
      STATEV(ns:ns+5) = s%et
c     = s%etrav  
      ns=ns+6
      STATEV(ns:ns+35) = pack(s%ass2,.true.)
      ns=ns+36
      STATEV(ns:ns+5) = s%alfass
      ns=ns+6
      ng=0
      do iph=1,s%nph
          STATEV(ns) = s%phase(iph)%ngrain
          ns=ns+1
          STATEV(ns:ns+8) = pack(s%phase(iph)%fij,.true.)
          ns=ns+9
          STATEV(ns:ns+11) = pack(s%phase(iph)%axor,.true.)
          ns=ns+12
          STATEV(ns:ns+2) = s%phase(iph)%euler
          ns=ns+3
          STATEV(ns:ns+35) = pack(s%phase(iph)%aef,.true.)
          ns=ns+36
          ns=ns+6
          do ng_ph=1,s%phase(iph)%ngrain
              ng=ng+1
              STATEV(ns:ns+5) = s%phase(iph)%grain(ng_ph)%st
              ns=ns+6
              STATEV(ns:ns+5) = s%phase(iph)%grain(ng_ph)%et
              ns=ns+6
              STATEV(ns:ns+8) = pack(s%phase(iph)%grain(ng_ph)%r,.true.)
              ns=ns+9
              STATEV(ns) = s%phase(iph)%grain(ng_ph)%wgt
              ns=ns+1
              STATEV(ns) = s%phase(iph)%grain(ng_ph)%ngrnph
              ns=ns+1
              if(kCL.eq.1)call write_statv_hl1(
     #            s%phase(iph)%grain(ng_ph)%hl1,ng,ns,NSTATV,STATEV)
              if(kCL.eq.2)call write_statv_hl2(
     #            s%phase(iph)%grain(ng_ph)%hl2,ng,ns,NSTATV,STATEV)
              if (iBackStress.eq.1) then
                  call write_statv_bs(s%phase(iph)%grain(ng_ph)%bs,ng
     #            ,ns,NSTATV,STATEV)
              endif
              if (itwinning.eq.1) then
                  call write_statv_tw(s%phase(iph)%grain(ng_ph)%tw,ng
     #            ,ns,NSTATV,STATEV)
              endif
              if (iPhTr.eq.1) then
                  call write_statv_pt(s%phase(iph)%grain(ng_ph)%pt,ng
     #            ,ns,NSTATV,STATEV)
              endif
          enddo
      enddo
      
      !write
      if (iUmat.eq.0) then
          iunit=1
          open(unit=iunit,file=fileprev,status='unknown')    
          write(iunit,*) 'STATE VARIABLES' !title
          do i=1,ns
              write(1,1) STATEV(i)
          enddo
          close(iunit)
      endif
      
      end subroutine write_statv
      
      end module state_var   
c      
c***********************************************************************
c *** output ***********************************************************
c***********************************************************************
      module output
      
      CHARACTER*78 prosa
      
      contains 
      
      subroutine open_file
      
      use flags
      use sample_props_v, only : nph
      
      character(len=40) :: line

  200 FORMAT('epsc7_ph',i2.2,'.out') 
  201 FORMAT('epsc3_ph',i2.2,'.out')  
  202 FORMAT('epsc8_ph',i2.2,'.out')   
  203 FORMAT('epsc10_ph',i2.2,'.out')
  204 FORMAT('epsc9_ph',i2.2,'.out')
      

      OPEN(unit=11,file='epsc1.out',status='unknown') 
      OPEN(unit=12,file='epsc2.out',status='unknown')
      OPEN(unit=13,file='epsc3.out',status='unknown')
      OPEN(unit=14,file='epsc4.out',status='unknown')
      OPEN(unit=15,file='epsc5.out',status='unknown')
      OPEN(unit=16,file='epsc6.out',status='unknown')
      OPEN(unit=18,file='epsc8.out',status='unknown')
     

      do iph=1,nph
        write(line,200) iph
        iunit=170+iph
        OPEN(unit=iunit,file=line,status='unknown')   
        write(line,201) iph        
        iunit=130+iph
        OPEN(unit=iunit,file=line,status='unknown')    
        write(line,202) iph        
        iunit=180+iph
        OPEN(unit=iunit,file=line,status='unknown')    
        !epsc9 & epsc10
        if (i_diff_dir.eq.1) then
          write(line,204) iph        
          iunit=190+iph
          OPEN(unit=iunit,file=line,status='unknown')           
          write(line,203) iph        
          iunit=200+iph
          OPEN(unit=iunit,file=line,status='unknown')  
        endif
      enddo     
      
      !open temporary files
      OPEN(300,file='temp.out',status='unknown')
      
      end subroutine open_file
      
      subroutine close_file
      
      do iunit=11,20
        CLOSE(unit=iunit)
      enddo

      do iph=1,nph
          iunit=130+iph
          close(iunit)
          iunit=170+iph
          close(iunit)
          iunit=180+iph
          close(iunit)     
          iunit=190+iph
          close(iunit) 
          iunit=200+iph
          close(iunit) 
      enddo
      
      !open temporary files
      close(300)
      
      end subroutine close_file
      
      subroutine write_header(label)
      
      use flags
      use sample_props_v, only : nph
      use diffract
      use files
      use mphase_props
      
      character(len=150),intent(in) :: label
      
    2 FORMAT(1h ,a)
    3 FORMAT(1h ,78('*'))     
    4 FORMAT(1h ,8('*'),' SAMPLE data - File: ',a,9('*'))
    
      
      !Writes the title in output files
      do iunit=11,19
        if (iunit.ne.17.and.iunit.ne.19.and.iunit.ne.20) then
          write(iunit,3)
          write(iunit,2) label
          write(iunit,3)
        endif
      enddo
      
      !copy input files to epsc1.out
      call write_epsc_1_ini
    
      !epsc3
      write(13,2) 'COMPONENTS 11 22 33 OF SAMPLE STRAIN,
     # STRESS, ELASTIC STRAIN and AVACS' 
      write(13,3)
      write(13,'(8x,''et11'',8x,''et22'',8x,''et33'',8x,''et23'',8x,
     #           ''et13'',8x,''et12'',11x,
     #           ''st11'',8x,''st22'',8x,''st33'',8x,
     #           ''st23'',8x,''st13'',8x,''st12'',9x,
     #           ''etel11'',6x,''etel22'',6x,''etel33'',8x,
     #           ''etth11'',6x,''etth22'',6x,''etth33'')') !,8x''
c     #           ''avacs'',8x,
c     #           ''rho_avg'',8x,''rho_for'',8x,''rho_deb'')')
      do iph=1,nph
        iunit=130+iph
        write(iunit,2) 'COMPONENTS 11 22 33 OF SAMPLE STRAIN,
     # STRESS, ELASTIC STRAIN and AVACS' 
        write(iunit,3)
        write(iunit,'(8x,''et11'',8x,''et22'',8x,''et33'',8x,''et23'',8x
     #           ,''et13'',8x,''et12'',11x,
     #           ''st11'',8x,''st22'',8x,''st33'',8x,
     #           ''st23'',8x,''st13'',8x,''st12'',9x,
     #           ''etel11'',6x,''etel22'',6x,''etel33'',8x,
     #           ''wgt_ph'')') 
      enddo
     
      !epsc7
      do iph=1,nph
        iunit=170+iph     
        write(iunit,3)
        write(iunit,2) label
        write(iunit,3)
        if(iTotStep.eq.0) then !MZ_2ph output write identification for each file
          write(iunit,'(''           Ctrl'')', advance='no')
          do i=1,3
            write(iunit,'(''            EP'',i1)', advance='no') i
          enddo
          do i=1,3
            write(iunit,'(''           SIG'',i1)', advance='no') i
          enddo
          do i=1,nmodes(iph)
            write(iunit,'(''          MODE'',i1)', advance='no') i
          enddo
          write(iunit,'(''          ActAv'')', advance='no')
          do i=1,nmodes(iph)
            write(iunit,'(''         PMODE'',i1)', advance='no') i
          enddo
          write(iunit,'(''         PActAv'')', advance='no')
          do i=1,nmodes(iph)
            write(iunit,'(''         CMODE'',i1)', advance='no') i
          enddo
          write(iunit,'(''         CActAv'')', advance='no')
          do i=nslmod(iph)+1,nmodes(iph)
            write(iunit,'(''       PTVFMOD'',i1)', advance='no') i
          enddo
          do i=nslmod(iph)+1,nmodes(iph)
            write(iunit,'(''       CTVFMOD'',i1)', advance='no') i
          enddo
c          write(iunit,'(''        PTWINVF'')', advance='no')
c          write(iunit,'(''        CTWINVF'')', advance='no')
          write(iunit,'(''         nGrain'')', advance='no')
          write(iunit,'(''       MaxTwins'')', advance='no')
          write(iunit,*)
        endif
      enddo      
      
      !epsc8
      write(18,2) 'EQUIVALENT STATES'
      write(18,2) 'EQ ET - EQ PL ET - EQ EIG ET - EQ ST - EQ ETR
     #- EQ PL ETR'
     #           ,' - EQ STR - VOLUME - PRESSURE - WTOT - WPLTOT'
     #           ,' - STET - STETAV'
      
      !epsc9 & epsc10
      if (i_diff_dir.eq.1) then
      do iph=1,nph
        iunit=190+iph
        write(iunit,3)
        write(iunit,2) label
        write(iunit,3)
        write(iunit,2) 'EVOLUTION OF INTERNAL STRAINS'
        write(iunit,3)

        iunit=200+iph
        write(iunit,3)
        write(iunit,2) label
        write(iunit,3)
        write(iunit,'(''  Ctrl       '')', advance='no')
        do i=1,3
          write(iunit,'(''  EP'',I1,''        '')', advance='no') i
        enddo
        do i=1,3
          write(iunit,'(''  SIG'',I1,''       '')', advance='no') i
        enddo
        do i=1,3
          write(iunit,'(''  EPAV'',I1,''      '')', advance='no') i
        enddo
        do i=1,3
          write(iunit,'(''  SIGAV'',I1,''     '')', advance='no') i
        enddo
        do i=1,ndiff(iph)
          write(iunit,'(''  DIF'',I3.3,''     '')', advance='no') i
        enddo
        do i=1,ndiff(iph)
          write(iunit,'(''  STDEV'',I3.3,''   '')', advance='no') i
        enddo
        do i=1,ndiff(iph)
          write(iunit,'(''  DEV_DIF'',I3.3,'' '')', advance='no') i
        enddo
        write(iunit,'('','')',advance='no')
        do i=1,ndiff(iph)
         write(iunit,'('' DEV_STDEV'',I3.3,'''')', advance='no') i
        enddo
        do i=1,ndiff(iph)
          write(iunit,'(''  WGT'',I3.3,''     '')', advance='no') i
          wgtsetini(i)=wgtset(i)
        enddo
        do i=1,ndiff(iph)
         write(iunit,'('' DEC11'',I3.3,''    '')', advance='no') i
        enddo
        do i=1,ndiff(iph)
         write(iunit,'('' DEC22'',I3.3,''    '')', advance='no') i
        enddo
        write(iunit,*)
      enddo
      endif
      
      end subroutine write_header
      
      subroutine write_file(istep,iproc,step,temp)
      
      use diffract, only : dif_planes
      
      use const
      use flags
      use mphase_props 
      use grain_props_v
      use sample_props_v 
      use mphase_state
      use grain_state
      use sample_state 
      use mphase_rate
      use grain_rate
      use sample_rate_v
      use bc
      use solve_t
      use hard_law1_v
      
      real :: aux6(6)
      character(len=150) :: cdummy
      
      !epsc3
      aux6=0.0
      do ng=1,ngrain
          aux6=aux6+etthcs(:,ng)*wgt(ng)
      enddo
      write(13,'(2(6e12.4,3x),2(3e12.4,3x),f10.4,3x
     #   ,3(e12.4,3x))')(etss(i)-etssref(i),i=1,6),
     #   (stss(j),j=1,6),(etelss(k),k=1,3),(aux6(k),k=1,3)
      do iph=1,nph
          iunit=130+iph
          write(iunit,'(2(6e12.4,3x),2(3e12.4,3x),f10.4,3x
     #        ,3(e12.4,3x))')(etss_ph(i,iph)-etssref(i),i=1,6),
     #        (stss_ph(j,iph),j=1,6),(etelss_ph(k,iph),k=1,3)
     #        ,wgt_ph(iph)
      enddo
      
      !epsc8
      call EFFECTIVE_MAGNITUDES(iproc,istep,step)
      
      !epsc9 & epsc10
      if (i_diff_dir.eq.1) then
          do iph=1,nph
              call dif_planes(cdummy,0.0
     #           ,istep,1,iph)
          enddo
      endif
      
      !write texture files
      if(ivarBC.eq.0.or.iproc.eq.nproc) then
          if(itexskip.ne.0.and.istep.ne.0) then
              nskip=(istep/itexskip)*itexskip
              if(istep.eq.1.or.istep.eq.nsteps.or.istep
     #              .eq.nskip) then
                  nFile=nFile+1
                  call WriteTexFile(nFile,iproc,istep)
              endif
          else
              if(istep.eq.nsteps) then
                  nFile=nFile+1
                  call WriteTexFile(nFile,iproc,istep)
              endif
          endif 
      endif
      !epsc7_ph
      do iph=1,nph
          call plasticity(iproc,istep,step,temp,1,iph)
      enddo
      
      end subroutine write_file
      
      subroutine write_epsc_1_ini
      
      use files
      
    2 FORMAT(1h ,78('*'))
    3 FORMAT(1h ,8('*'),' CRYSTAL data - File: ',a,9('*'))
    4 FORMAT(1h ,8('*'),' SAMPLE data - File: ',a,9('*'))
      
      !write initial
      !Open and read main control file: "EPSC4.IN"   -> unit=1
      OPEN(unit=1,file='epsc4.in',status='old')     
      !writes simulation master file into 'epsc1.out'
      WRITE(11,*)
      WRITE(11,'('' ****** SIMULATION MASTER FILE *******'')')
      DO IDUM=1,200
        READ(UNIT=1,END=99,FMT='(A)') PROSA
        WRITE(11,'(A)') PROSA
      ENDDO
   99 REWIND 1
      WRITE(11,'('' ****** END OF SIMULATION MASTER FILE '')')
      WRITE(11,*)
      close(1)
	  
c *** Opens the file with the CRYSTAL data "filecrys"
      OPEN(unit=1,file=filecrys(1),status='old')
      write(11,3) filecrys(1)
      write(11,2)
C     WRITES CRYSTAL DATA FILE INTO 'EPSC1.OUT'
      WRITE(11,*)
      WRITE(11,'('' ****** CRYSTAL DATA FILE *******'')')
      DO IDUM=1,200
        READ(UNIT=1,END=98,FMT='(A)') PROSA
        WRITE(11,'(A)') PROSA
      ENDDO
   98 REWIND 1
      WRITE(11,'('' ****** END OF CRYSTAL DATA FILE '')')
      WRITE(11,*)
      close(1)
	  
c *** Opens the file with the SAMPLE data "filesamp"
      OPEN(unit=1,file=filesamp,status='old')
      write(11,4) filesamp(1)
      write(11,2)
C     WRITES TEXTURE DATA FILE INTO 'EPSC1.OUT'
      OPEN(unit=1,file=filesamp(1),status='old')
      WRITE(11,*)
      WRITE(11,'('' ****** TEXTURE DATA FILE (partial) *******'')')
      DO IDUM=1,20
        READ(UNIT=1,END=97,FMT='(A)') PROSA
        WRITE(11,'(A)') PROSA
      ENDDO
   97 REWIND 1
      WRITE(11,'('' ****** END OF TEXTURE DATA FILE '')')
      WRITE(11,*)
      close(1)
	  
      OPEN(unit=1,file=fileproc,status='old')
C     WRITES PROCESS DATA FILE INTO 'EPSC1.OUT'
      WRITE(11,*)
      WRITE(11,'('' ****** PROCESS DATA FILE *******'')')
      DO IDUM=1,200
        READ(UNIT=1,END=96,FMT='(A)') PROSA
        WRITE(11,'(A)') PROSA
      ENDDO
   96 REWIND 1
      WRITE(11,'('' ****** END OF PROCESS DATA FILE '')')
      WRITE(11,*)
      close(1)
      
      if(i_diff_dir.eq.1) then
          OPEN(unit=1,file=filediff(1),status='old')
C -----------------------------------------------------------------
C     WRITES DIFFRACTION DATA FILE INTO 'EPSC1.OUT'
          WRITE(11,*)
          WRITE(11,'('' ****** DIFFRACTION DATA FILE *******'')')
          DO IDUM=1,200
            READ(UNIT=1,END=95,FMT='(A)') PROSA
            WRITE(11,'(A)') PROSA
          ENDDO
   95     REWIND 1
          WRITE(11,'('' ****** END OF DIFFRACTION DATA FILE '')')
          WRITE(11,*)
          close(1)
      endif
      
      end subroutine
      
      SUBROUTINE WriteTexFile(nFile,iProc,iStep)
      
      use miscellaneous_sub, only : euler
      use sample_props_v, only : ngrain
      use sample_state, only : etss,stss
      use mphase_state, only : axisph,eulerph 
      use grain_props_v, only : r,phi,the,ome,wgt
      use bc_v, only : iTotStep
      
      DIMENSION RG(3,3)
      CHARACTER*20 filename
        CHARACTER Number(10)
        data Number/'0','1','2','3','4','5','6','7','8','9'/
10      FORMAT(' EPSC Step= ',I5,I5,I5)
15      FORMAT(' EPSC Step= ',I5,I5,I5,E13.4,E13.4,E13.4,E13.4
     #         ,E13.4,E13.4)
20      FORMAT(' B ',I10,' 0')
30      FORMAT(F12.2,F12.2,F12.2,F12.8)
c       Creating the filename using the number passed
c       to the subroutine (nFile)
        if(nFile.gt.9999) then
                print*,'Maximum number of texture files is 9999'
                stop
        else
                itho=nFile/1000
                ihun=(nFile-itho*1000)/100
                iten=(nFile-itho*1000-ihun*100)/10
                ione=nFile-itho*1000-ihun*100-iten*10
                filename='tex'//Number(itho+1)//Number(ihun+1)//
     #            Number(iten+1)//Number(ione+1)//'.out'
        endif
      OPEN(unit=22,file=filename,status='unknown')
        write(22,'(A)')'TEXTURE AT STRAIN' !MZ
        write(22,'(6f9.3)') (axisph(0,I,1),I=1,3),(eulerph(i,1),i=1,3)
c        write(22,10) iTotStep,iProc,iStep
        write(22,15) iTotStep,iProc,iStep
     #    ,etss(1),etss(2),etss(3),stss(1),stss(2),stss(3)
        write(22,20) ngrain

        DO NG=1,NGRAIN
          DO I=1,3
            DO J=1,3
              rg(i,j)=r(I,J,NG)
            END DO
          END DO
          CALL EULER(1,PHI(NG),THE(NG),OME(NG),RG)
          write(22,30) phi(ng),the(ng),ome(ng),wgt(ng)
        enddo

c        do ng=1,ngrain
c                write(22,30) phi(ng),the(ng),ome(ng),wgt(ng)
c        enddo
      CLOSE(unit=22)
      END SUBROUTINE
      
      SUBROUTINE plasticity (iproc,istep,step,temp,ioption,iph) !MZ_2ph output over phases, each write unit depends on input phase num  
c
c **********************************************************************
c
c     SUBROUTINE plasticity   --->   version 11/jul/02
c
c *** Evaluates the plastic activity in each step of the process     ***
c *** Average active systems, shears, twinning volume fractions      ***
c *** Relative activity of each plastic mechanism                    ***
c *** Plastic activity in diffracting sets (for Tom Holden)          ***
c **********************************************************************
      
      use const
      use flags
      use mphase_props
      use mphase_state
      use mphase_rate
      use grain_props
      use grain_state
      use grain_rate
      use sample_props_v
      use sample_state
      use sample_rate_v
      use twinning
      use diffract
      use bc
      
      DIMENSION WP(NGR),WC(NGR),aux(NMOD),aux_pa(NMOD),aux_ch(NMOD)
      DIMENSION shear_mod_acum(NMOD),shear_mod(NMOD),
     #                shear_dif_acum(NDIFFX)
c       BJCL include activity in child grains
      DIMENSION shear_mod_acum_ch(NMOD),
     #          shear_mod_acum_pa(NMOD),
     #          shear_mod_ch(NMOD),shear_dif_acum_ch(NDIFFX),
     #          shear_mod_pa(NMOD),shear_dif_acum_pa(NDIFFX)
      DIMENSION jact(0:NSLS)
c       BJCL include activity in child grains

    1 FORMAT(1h ,'STATISTIC OF SLIP ACTIVITY:')
    2 FORMAT(1h ,'Average number of load_active systems:',2x,f8.4)
    3 FORMAT(1h ,'No plastic activity')
    4 FORMAT(1h ,'Fraction of total shear:')
    5 FORMAT(1h ,'Mode #:',i3,2x,'Activity:',1x,f8.4,1x,'%')
    6 Format(1h ,'Twinning volume fraction:')
    7 FORMAT(1h ,'No twinning activation.')
    8 FORMAT(1h ,'Twin mode:',i3,2x,'Volume fraction:',1x,f8.4)
    9 FORMAT(1h ,'Total number of grains that plastify:',i6)
c   10 FORMAT(1h , f10.5,3x,15f7.3)
   10 FORMAT(1h , 100F15.5)
   18 FORMAT(1h ,'NO PLASTIC ACTIVITY')
   19 FORMAT(1h ,'Percent respect to total shear:')
   20 FORMAT(1h ,'Mode #:',i3,2x,'Activity:',1x,f8.4,1x,'%')
   21 FORMAT(1h ,'NO TWIN ACTIVITY')
   22 FORMAT(1h ,'Twin mode:',i3,2x,'Volume fraction:',1x,f10.6)
   23 FORMAT(1h ,'               ','Percent:',1x,f8.4)
c       BJCL
   30 FORMAT(1h , f11.6,3x,30f11.6)     ! BJCL
     
      !MZ_2ph output define start and end grain number
      ng1=SUM(nphngr(1:iph))-nphngr(iph)+1
      ng2=SUM(nphngr(1:iph))
      !END
c ______________________________________________________________________

      IF(IOPTION.EQ.0) THEN
        if(ndiff(iph).gt.0) then
          do nd=1,ndiff(iph)
            shear_dif_acum(nd)=0.0
            shear_dif_acum_ch(nd)=0.0
            shear_dif_acum_pa(nd)=0.0
          enddo
        endif
        do mo=1,nmodes(iph)  !MZ_2ph output no need to go over phases
          shear_mod_acum(mo)=0.0
          shear_mod_acum_ch(mo)=0.0
          shear_mod_acum_pa(mo)=0.0
          vfrac_mod_acum(mo)=0.0
        enddo
      ENDIF

c ______________________________________________________________________

      IF(IOPTION.EQ.1) THEN

      write(12,*) 
      ngtotact=0
      actav=0.0
      do mo=1,nmodes(iph) 
        shear_mod(mo)=0.0
      enddo

      do i=0,10
        jact(i)=0
        actwgt(i)=0.0                         !added for james by jn
      end do
      jminact=25
      jmaxact=0

c     First the total average
      wgtph=0.0 !MZ_2ph output
      do ng=ng1,ng2
        wgtph=wgtph+wgt(ng) !MZ_2ph output
        if (nact(ng).gt.jmaxact) jmaxact=nact(ng)                  !added for james by jn
        if (nact(ng).lt.jminact) jminact=nact(ng)
        jact(nact(ng))=jact(nact(ng))+1
        actwgt(nact(ng))=actwgt(nact(ng))+wgt(ng)

        if (nact(ng).ne.0) then
          ngtotact=ngtotact+1
          actav=actav+nact(ng)*wgt(ng)
          nst=0
          do mo=1,nmodes(iph) 
            do isys=1,nsm(mo,iph)
              nst=nst+1
              shear_mod(mo)=shear_mod(mo)+gamd(nst,ng)*wgt(ng)
              shear_mod_acum(mo)=shear_mod_acum(mo)
     #                           +gamd(nst,ng)*wgt(ng)
            enddo
          enddo
        endif
      enddo
      actav=actav/wgtph !MZ_2ph output

c      IF (istep.eq.1.and.iproc.eq.1) THEN
c        WRITE(21,*)                                                     !added for james by jn
c        WRITE(21,*)
c        WRITE(21,*) 'SYSTEM ACTIVITY STATISTICS FOR JAMES'
c        WRITE(21,*) ' STEP      STRAIN     STRESS     AVACT ',
c     #     'ACTMIN ACTMAX   N0   N1   N2   N3   ',
c     #     'N4   N5   N6   N7   N8   N9  N10     ',
c     #     'W0     W1     W2     W3     W4     W5     W6    ',
c     #     'W7     W8     W9     W10'
c      endif
c
c      if (istep.ge.1) then
c        WRITE(21,'(I6,2E12.4,F9.5,2I7,11I5,11F7.3)') istep,                  !added for james by jn
c     #  etss(i_control_var)-etssref(i_control_var),stss(i_control_var),
c     #  actav,jminact,jmaxact,(jact(i),i=0,10),(actwgt(i),i=0,10)
c      endif



c     Then the parent grains
      ngtotact_pa=0
      actav_pa=0.0
      do mo=1,nmodes(iph)
        shear_mod_pa(mo)=0.0
      enddo
      TWP=0.0
      do ng=1,ngParent
        TWP=TWP+wgt(ng)
        enddo
      do ng=1,ngParent
        WP(ng)=wgt(ng)/TWP
        enddo
      do ng=1,ngParent
        if (nact(ng).ne.0) then
          ngtotact_pa=ngtotact_pa+1
          actav_pa=actav_pa+nact(ng)*WP(ng)
          nst=0
          do mo=1,nmodes(iph) 
            do isys=1,nsm(mo,iph)
              nst=nst+1
              shear_mod_pa(mo)=shear_mod_pa(mo)+gamd(nst,ng)*WP(ng)
              shear_mod_acum_pa(mo)=shear_mod_acum_pa(mo)
     #                           +gamd(nst,ng)*WP(ng)
            enddo
          enddo
        endif
      enddo
c     Then the child grains
      ngtotact_ch=0
      actav_ch=0.0
      do mo=1,nmodes(iph) 
        shear_mod_ch(mo)=0.0
      enddo
      TWC=0.0
      do ng=ngParent+1,ngrain
          TWC=TWC+wgt(ng)
        enddo
      do ng=ngParent+1,ngrain
          WC(ng)=wgt(ng)/TWC
        enddo
      do ng=ngParent+1,ngrain
        if (nact(ng).ne.0) then
          ngtotact_ch=ngtotact_ch+1
          actav_ch=actav_ch+nact(ng)*WC(ng)
          nst=0
          do mo=1,nmodes(iph) 
            do isys=1,nsm(mo,iph)
              nst=nst+1
              shear_mod_ch(mo)=shear_mod_ch(mo)+gamd(nst,ng)*WC(ng)
              shear_mod_acum_ch(mo)=shear_mod_acum_ch(mo)
     #                           +gamd(nst,ng)*wgt(ng)
            enddo
          enddo
        endif
      enddo
      shear_tot=0.0
      shear_tot_ch=0.0
      shear_tot_pa=0.0
      do mo=1,nmodes(iph) 
        shear_tot=shear_tot+shear_mod(mo)
        shear_tot_ch=shear_tot_ch+shear_mod_ch(mo)
        shear_tot_pa=shear_tot_pa+shear_mod_pa(mo)
      enddo

c *** Write the results in "EPSC?.OUT"
      if (icvx.eq.0) xref=temp
      if (icvx.ge.1.and.icvx.le.6) xref=etss(icvx)-etssref(icvx)
      if (icvx.ge.7) xref=stss(icvx-6)

      write(12,1)
      write(12,2) actav
      write(12,3)
        if(shear_tot.ne.0.0) then
          do mo=1,nmodes(iph) 
            aux(mo)=shear_mod(mo)/shear_tot
          enddo
        else
          do mo=1,nmodes(iph) 
            aux(mo)=shear_mod(mo)
          enddo
        endif
        if(shear_tot_pa.ne.0.0) then
          do mo=1,nmodes(iph) 
            aux_pa(mo)=shear_mod_pa(mo)/shear_tot_pa
          enddo
        else
          do mo=1,nmodes(iph)
            aux_pa(mo)=shear_mod_pa(mo)
          enddo
        endif
        if(shear_tot_ch.ne.0.0) then
          do mo=1,nmodes(1) 
            aux_ch(mo)=shear_mod_ch(mo)/shear_tot_ch
          enddo
        else
          do mo=1,nmodes(1) 
            aux_ch(mo)=shear_mod_ch(mo)
          enddo
        endif
        iunit=170+iph !MZ_2ph output
        write(iunit,'(1400F15.5)') xref, !MZ_2ph output depends on the phase 
     #    ((etss(i)-etssref(i)),i=1,3),(stss(i),i=1,3),
c     #    (etav(i),i=1,3),(stav(i),i=1,3),
     #    (aux(mo),mo=1,nmodes(iph)),actav 
     #    ,(aux_pa(mo),mo=1,nmodes(iph)),actav_pa, !MZ_2ph output commented out
     #    (aux_ch(mo),mo=1,nmodes(iph)),actav_ch, !MZ_2ph output commented out
     #    (PTVFM(mo),mo=1,ntwmod(iph)),
     #    (CTVFM(mo),mo=1,ntwmod(iph)),real(ngrain),real(MaxTwins) !MZ_2ph output commented out

c     ngrprn=1
c     do ng=1,ngrprn
c       if (nact(ng).ne.0) then
c         write(12,*)
c         write(12,'(''NG:'',i3,''  NACT:',i3,''  IACT:'',24i3)')
c    #                 ng,nact(ng),(iact(ns1,ng),ns1=1,nact(ng))
c       endif
c     enddo

c *** Calculates accum. plastic shear for each set of diffracting grains.
c *** Does an average per grain instead of a weighted average. The idea is
c *** to make the result independent of the volume fraction of grains
c *** contained in each diffracting set.

c      if (ndiff.gt.0) then
c        do nd=1,ndiff
c          dummy=0.0
c          do ng=1,ngrset(nd)
c            ngset=igrset(nd,ng)
c            nst=0
c            do mo=1,nmodes(1)
c              do isys=1,nsm(mo)
c                nst=nst+1
c                dummy=dummy+abs(gamd(nst,ngset))
c              enddo
c            enddo
c          enddo
c          shear_dif_acum(nd)=shear_dif_acum(nd)+dummy/ngrset(nd)
c        enddo
c        write(20,'(i5,f10.4,/,(5x,6f9.5))') ns,xref,
c     #                           (shear_dif_acum(nd),nd=1,ndiff)
c      endif

      ENDIF
c ______________________________________________________________________

      IF(IOPTION.EQ.2) THEN

        write(12,*)
        shear_tot_acum=0.0
        vfrac_tot_acum=0.0
        do mo=1,nmodes(iph) 
          shear_tot_acum=shear_tot_acum+shear_mod_acum(mo)
          vfrac_tot_acum=vfrac_tot_acum+vfrac_mod_acum(mo)
        enddo
        if (shear_tot_acum.eq.0.0) then
          write(12,18)
        else
          write(12,19)
          do mo=1,nmodes(iph) 
            write(12,20) mo,shear_mod_acum(mo)/shear_tot_acum*100.0
          enddo
        endif
        if (vfrac_tot_acum.eq.0.0) then
          write(12,21)
        else
          do mo=1,nmodes(iph) 
            if (itw(mo).eq.1) then
              write(12,22) mo,vfrac_mod_acum(mo)
              write(12,23) vfrac_mod_acum(mo)/vfrac_tot_acum*100.0
            endif
          enddo
        endif

      ENDIF
c ______________________________________________________________________
      return
      end subroutine
      
      subroutine write_temp
      
      use mvoigt
      use meshelby
      use sc_estimate
      use miscellaneous_sub
      use diffract
      use state_var
      use back_stress
      use back_stress_v
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
      use solve_t
      use hard_law1_v
      use phase_transf
      use twinning
      use output_v
      
      DIMENSION CRSS_avg_all(nsls),rho_avg_all(nsls)
      DIMENSION TVF_p(2)
      real :: aux6(6),aux66(6,6),aux5(6),rss(NSLS)
      
      !find plastic strain
c      aux5=0.0
c      do ng=1,ngrain
c          call invten(ccs2(:,:,ng),aux66)
c          call tens_mult_1(aux6,aux66,stcs(:,ng)) !el strain
c          aux6=etcs(:,ng)-aux6 !pl strain
c          aux5=aux5+aux6*wgt(ng)
c      enddo
c      
c      WRITE(300,'(50f24.12)') aux5
      
c      WRITE(300,'(50f24.12)') TEMP(1:3)
      
      tau_avg=0.0
      do ng=1,nphngr(1)
          tau_avg=tau_avg+sum(tau(1:nslsys(1),ng))*wgt(ng)/48.0
      enddo
      
      rho_avg=0.0
      do ng=1,nphngr(1)
          rho_avg=rho_avg+sum(rho_tot(1:nslsys(1),ng))*wgt(ng)/2.0
      enddo
      
      rho_avg_act=0.0
      rho_avg_inact=0.0
      rho_tot_0_avg=0.0
      do ng=1,nphngr(1)
          do ns=1,nslsys(1)/2
              if (iact_sys(ns,ng).eq.1) then
                  rho_avg_act=rho_avg_act
     #                +rho_rev(ns*2-1,ng)*wgt(ng)
                  rho_avg_inact=rho_avg_inact
     #                +rho_rev(ns*2,ng)*wgt(ng)
                  
              elseif (iact_sys(ns,ng).eq.2) then
                  rho_avg_act=rho_avg_act
     #                +rho_rev(ns*2,ng)*wgt(ng)
                  rho_avg_inact=rho_avg_inact
     #                +rho_rev(ns*2-1,ng)*wgt(ng)
              endif
          enddo  
      enddo
      
      ng=1
      do is=1,nslsys(ngrnph(ng))
          imo=iSysMode(is,ngrnph(ng))
          !find rss
          rss(is)=0.0
          do i=1,6 
              rss(is)=rss(is)+(mcs(i,is,ng)+nmcs(i,is,ng))
     #               *stcs(i,ng)*profac(i)  !inonSch
          enddo
          if(iBackStress.eq.1) rss(is)=rss(is)-tau_bcst(is,ng
     #       ,ilatBS)
      enddo
      
      !average out
c      WRITE(300,'(50e25.12)') etss(1),stss(1),rho_avg,rho_avg_act
c     #    ,rho_avg_inact,tau_avg
      
      !grain out
c      WRITE(300,'(1000e25.12)')rssmin(1:nslsys(1),1),rss(1:nslsys(1))
c     #    ,gam_acc(1:nslsys(1),1),rho_tot(1:nslsys(1),1)
c     #    ,rho_rev(1:nslsys(1),1),rho_forw(1:nslsys(1),1)
      
      WRITE(300,'(1000e25.12)')etss(1),tau_bcst_tot(1)
     #    ,tau_bcst(5:6,1,1),tau_bcst(5:6,1,2)
     #    ,gam_acc(5:6,1)

      end subroutine
      
      subroutine write_umat(NOEL,NPT,KSPT,TIME,DTIME,STATEV,NSTATV
     #    ,DSTRAN,DFGRD1,DROT,STRESS,NTENS,error_avg_SDV,error_max_SDV
     #    ,lenoutdir,outdir,istart) 
      
      use mvoigt
      use meshelby
      use sc_estimate
      use miscellaneous_sub
      use diffract
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
      use solve_t
      use hard_law1_v
      use phase_transf
      use twinning
      
      dimension TIME(2),STATEV(NSTATV),STRESS(NTENS),DSTRAN(NTENS)
     #    ,DROT(3,3),DFGRD1(3,3),DFGRDPL(3,3),DFGRDEL(3,3),aux33(3,3)
     #    ,aux34(3,3),indx(3),CAUGREEN(3,3),ELLAG(3,3),ELSTIFF(3,3,3,3)
     #    ,CONJSTRSS(3,3),CAUSTRSS(3,3)
      real :: aux6(6),aux66(6,6),aux5(6)
      character(len=150) :: outdir
      
      if (istart.eq.1) then
          write(5500,*) 'NOEL: ',NOEL, 'KSPT: ',KSPT
          write(5500,*) 'TIME(2): ',TIME(2)
          write(5500,'(21e26.16)') (DSTRAN(i),i=1,NTENS),((DROT(i,j)
     #        ,i=1,3),j=1,3),(STRESS(i),i=1,NTENS)

      elseif(istart.eq.2) then
          
          !field output 
          STATEV(NSTATV-20-3)=iout
          STATEV(NSTATV-20-2)=error_avg_SDV
          STATEV(NSTATV-20-1)=error_max_SDV
          !field output 
c          do ng=1,ngrain
c            CALL EULER(1,PHI(ng),THE(ng),OME(ng),r(:,:,ng))
c          enddo
c          ns=NSTATV-27-ngrain*4
c          do ng=1,ngrain
c            STATEV(ns)=phi(ng)
c            ns=ns+1
c            STATEV(ns)=the(ng)
c            ns=ns+1
c            STATEV(ns)=ome(ng)
c            ns=ns+1
c            STATEV(ns)=wgt(ng)
c            ns=ns+1
c          enddo
          
          !perform volume average of stcs and etcs
          aux5=0.0
          aux6=0.0
          do ng=1,ngrain
              aux5=aux5+etcs(:,ng)*wgt(ng) !pl strain
              aux6=aux6+stcs(:,ng)*wgt(ng)
          enddo
          STATEV(NSTATV-200:NSTATV-200+5)=aux5
          STATEV(NSTATV-200+6:NSTATV-200+11)=aux6
          
          !find plastic strain
          aux5=0.0
          do ng=1,ngrain
              call invten(ccs2(:,:,ng),aux66)
              call tens_mult_1(aux6,aux66,stcs(:,ng)) !el strain
              aux6=etcs(:,ng)-aux6 !pl strain
              aux5=aux5+aux6*wgt(ng)
          enddo
c          STATEV(NSTATV-6:NSTATV-1)=aux5
c          
c          STATEV(NSTATV-50:NSTATV-50+5)=stcs(:,1)
c          STATEV(NSTATV-44:NSTATV-44+5)=etcs(:,1)
c          STATEV(NSTATV-38:NSTATV-38+8)=reshape(drotcs(:,:,1),(/9/))
          
          !state variable output
c          nstv=NSTATV-200
c          STATEV(nstv:nstv+23)=gam_acc(1:24,1)
c          nstv=nstv+24
cc          STATEV(nstv:nstv+23)=tau_bcst(1:24,1,1)
c          STATEV(nstv:nstv+23)=tau(1:24,1)
c          nstv=nstv+24
          
      elseif(istart.eq.3) then
          !write state variables
          write(5500,*) 'Error in'
          write(5500,*) 'NOEL: ',NOEL, 'KSPT: ',KSPT
          write(5500,*) 'TIME(2): ',TIME(2)
          write(5500,'(21e26.16)') (DSTRAN(i),i=1,NTENS),((DROT(i,j)
     #        ,i=1,3),j=1,3),(STRESS(i),i=1,NTENS)
          write(5500,*) 'Write state variables'
          do i=1,NSTATV
              write(5500,'(e36.16)') STATEV(i)
          enddo
          STOP
      endif
      
      end subroutine
	  
      SUBROUTINE EFFECTIVE_MAGNITUDES (iproc,istep,step)

c
c **********************************************************************
c     SUBROUTINE EFFECTIVE_MAGNITUDES   --->   version 14/sep/01 (CNT)
c
c     Calculates and writes in unit #18 equivalent states and energies.
c **********************************************************************
c
      use mvoigt
      use meshelby
      use sc_estimate
      use miscellaneous_sub
      use diffract
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
      use solve_t
      use hard_law1_v
      use phase_transf
      use twinning
      use output_v
      
      real :: etsspl(6),etrsspl(6),stsshyd(6),stssdev(6),aux6(6)
     #    ,et_inc(6),aux7(6),etrsspl1(6)
      
      if(iproc.eq.1 .and. istep.eq.0) then
        wtotal=0.0
        wplastic=0.0
      endif
      
      !equivalent stress calculation
      strsseq=0.0
      pressure=(1.0/3.0)*(stss(1)+stss(2)+stss(3))
      stsshyd(1:3)=pressure
      stssdev=stss-stsshyd
      aJ2=0.0
      do i=1,6
          aJ2=aJ2+0.5*stssdev(i)*stssdev(i)*profac(i)
      enddo
      stsseq=sqrt(3.0*aJ2)
      
      !find plastic strain from gamd
      etrsspl1=0.0
      do ng=1,ngrain
          aux7=0.0
          do ns1=1,nact(ng)
              n1=iact(ns1,ng)
              aux7=aux7+mcs(:,n1,ng)*gamd(n1,ng)
          enddo
          etrsspl1=etrsspl1+aux7*wgt(ng)
      enddo
      !find plastic strain from crystal stress
      etrsspl1=0.0
      do ng=1,ngrain
          aux7=etrcs(:,ng)-etrcs_eig(:,ng)
          do i=1,6
              do j=1,6
                  aux7(i)=aux7(i)-scs2(i,j,ng)*strcs(j,ng)*profac(j)
              enddo
          enddo
          etrsspl1=etrsspl1+aux7*wgt(ng)
      enddo
      
          
      !equivalent strain measures
      etrsseq=0.0
      etsseigeq=0.0
      et_inc=etss-aux6
      do i=1,6
        !use stiffness for pl strain calculation
c        etrsspl(i)=et_inc(i)-etrss_eig(i)
c        do j=1,6
c          etrsspl(i)=etrsspl(i)-sss2(i,j)*strss(j)*profac(j)
c        enddo
        !use volume average over crystals
        etrsspl(i)=etrsspl1(i)
        etsseigeq=etsseigeq+etss_eig(i)*etss_eig(i)*profac(i)
        etrsseq=etrsseq+et_inc(i)*et_inc(i)*profac(i)
      enddo
      etsseigeq=sqrt(2.0)/3.0*sqrt(etsseigeq)
      etrsseq=sqrt(2.0)/3.0*sqrt(etrsseq)
      
      etsspl=etsspl+etrsspl
      
      !change of volume
      volume=(etss(1)+etss(2)+etss(3))
      
      !work calculation
      stet=0.0
      wrplastic=0.0
      wrtotal=0.0
      do i=1,6
        wrtotal=wrtotal+stss(i)*et_inc(i)*profac(i)*step
        wrplastic=wrplastic+stss(i)*etrsspl(i)*profac(i)*step
        stet=stet+stss(i)*etss(i)*profac(i)
      enddo
      wtotal=wtotal+wrtotal
      stetav=0.0
      do ng=1,ngrain
        do i=1,6
          stetav=stetav+stcs(i,ng)*etcs(i,ng)*profac(i)*wgt(ng)
        enddo
      enddo
      
      !redefine so that wplastic=etsspleq*stsseq
      if (stsseq.eq.0.0) then
          etrsspleq=0.0
      else
          etrsspleq=wrplastic/stsseq
          etrsseq=wrtotal/stsseq
      endif
      etsspleq=etsspleq+etrsspleq
      etsseq=etsseq+etrsseq
      
      !store etss of current step
      aux6=etss
      
      write(18,'(1x,13d12.4)') etsseq,etsspleq,etsseigeq,stsseq,etrsseq
     #           ,etrsspleq,strsseq,volume,pressure,wtotal,wplastic
     #           ,stet,stetav

      return
      end SUBROUTINE
	  
      end module output
c      
c***********************************************************************
c *** mumat ************************************************************
c***********************************************************************
      module mumat
      
      contains
      
      SUBROUTINE EPSCTOABQ_J(t1,t2)

      DIMENSION t1(6,6),t2(6,6),temp(6,6)
      
      temp = t1
      temp(:,4) = t1(:,6)
      temp(:,6) = t1(:,4)
      t2 = temp
      t2(4,:) = temp(6,:)
      t2(6,:) = temp(4,:)
      
      
      RETURN
      END SUBROUTINE
      
      SUBROUTINE J_shell(ass2,avgass2,nsteps)
      
      use miscellaneous_sub, only : ludcmpc,lubksbc
      
      DIMENSION avgass2(6,6),ass2(6,6)
      DIMENSION aux13(3,3),aux23(3,3),aux33(3,3),aux43(3,3)
      DIMENSION aux53(3,3),index1(3),index2(3),indx(3),y(3,3)
      DATA index1 /1,2,6/, index2 /3,4,5/
      
      ! define submatrices
      do i=1,3
          do j=1,3
              i1=index1(i) !1,2,6
              j1=index1(j) !1,2,6
              i2=index2(i) !3,4,5
              j2=index2(j) !3,4,5              
              aux13(i,j)=ass2(i1,j1)
              aux23(i,j)=ass2(i1,j2)
              aux33(i,j)=ass2(i2,j2)
              aux43(i,j)=ass2(i2,j1)
          enddo
      enddo
      
      ! define J without averaging
      ! y = aux33^-1
      do i=1,3 !Set up identity matrix.
          do j=1,3
              y(i,j)=0.
          enddo 
          y(i,i)=1.
      enddo 
      call ludcmpc(aux33,3,3,indx,d) !Decompose the matrix just once.
      do  j=1,3 !Find inverse by columns.
          call lubksbc(aux33,3,3,indx,y(1,j))
          !Note that FORTRAN stores two-dimensional matrices by column, so y(1,j) is the
          !address of the jth column of y.
      enddo      
      
      ! aux33 = aux13 - aux23*y*aux43
      aux53=matmul(y,aux43)
      aux33=matmul(aux23,aux53)
      aux33=aux13-aux33
      
      ! store the averaged matrix in avgass2
      avgass2(1:3,1:3)=avgass2(1:3,1:3)+aux33/real(nsteps)
      
      RETURN
      END SUBROUTINE
      
      subroutine add_bc_umat(NTENS,DSTRAN,STRAN,DROT,DTIME,TEMP,DTEMP
     #    ,i_temp_cij,i_ref_et,i_ref_st,i_bc_mode,i2d,etssref)
      
      use mvoigt
      use bc
      use sample_state_v, only : etss,stss
      
      !in
      integer, intent(in) :: NTENS
      real, intent(in)  :: DSTRAN(NTENS),DROT(3,3),DTIME,STRAN(NTENS)
      !out
      integer, intent(out) :: i_temp_cij,i_ref_et,i_ref_st,i_bc_mode,i2d
      real, intent(out) :: etssref(6)
      
      real :: aux61(6),det_full(3,3),aux66(6,6),aux3333(3,3,3,3)
     #    ,aux31(3,3),drot_full(3,3)
      
      !defines data from process file
      i_control_var=0
      i_bc_mode=0
      if(ntens.eq.6) i2d=0
      if(ntens.eq.3) i2d=1
      !define bc
      if (i2d.eq.0) then  
          ietbc=1
          aux61(1) = STRAN(1) + DSTRAN(1) - etss(1)
          aux61(2) = STRAN(2) + DSTRAN(2) - etss(2)
          aux61(3) = STRAN(3) + DSTRAN(3) - etss(3)
          aux61(4) = (STRAN(6) + DSTRAN(6))/2.0 - etss(4)
          aux61(5) = (STRAN(5) + DSTRAN(5))/2.0 - etss(5)
          aux61(6) = (STRAN(4) + DSTRAN(4))/2.0 - etss(6)
          call VOIGT(aux61,fulletbc,aux66,aux3333,1)
          istbc=0
          stbc=0.0
      endif
      !change to mixed BC for shell elements 2DUMAT        
      if (i2d.eq.1) then
          ! define strain tensor enfoced
          ietbc=0
          ietbc(1:2)=1
          ietbc(6)=1
          aux61(1) = DSTRAN(1) + STRAN(1)- etss(1)
          aux61(2) = DSTRAN(2) + STRAN(2)- etss(2)
          aux61(3) = 0.0
          aux61(4) = 0.0
          aux61(5) = 0.0
          aux61(6) = DSTRAN(3)/2.0 + STRAN(3)/2.0- etss(6)
          call VOIGT(aux61,fulletbc,aux66,aux3333,1)         
          ! define stress tensor enfoced
          istbc=0
          istbc(3:5)=1
          stbc(3)=0.0 !-stss(3)
          stbc(4)=0.0 !-stss(4)
          stbc(5)=0.0 !-stss(5)
      endif
      ! temperature (start = final, deltemp = 0)
      temp_s = TEMP
      temp_f = TEMP+DTEMP
      deltemp = DTEMP
      ! temp and referance st and et are not present
      i_temp_cij=1
      i_ref_et=0
      i_ref_st=0
      
      ! find strain rate
c      edot_macro = tnorm(det_full/DTIME,3,3)
      edot_macro = 0.0005 !fixed strain rate
      
      !find spin tensor from DROT and time inc
c      call drot2spin(DROT,DTIME,aux31)
c      drot_full = aux31*DTIME ! EPSC takes increment in rotation angle as input, not the spin
      
      ! divide strain inc from ABAQUS to smaller steps
      nproc = 1
c      nsteps= 1
      amax = MAXVAL(ABS(fulletbc))
c      thres_et=5.0e-4 !1e-4
      nsteps=NINT(amax/thres_et)
      if(nsteps.lt.nsteps_min) nsteps=nsteps_min

      call load_conditions
      
      !Calculates symmetric and antisymmetric portions of strain tensor
      call define_bc
      call define_bcr(dummy1,dummy1,i_bc_mode,0)

      !reference state for macro and grain state prior to beginning of process
      call et_st_ref(i_ref_et,i_ref_st,etssref)
      
      end subroutine add_bc_umat
      
      SUBROUTINE drot2spin(drot,dtime,W_app)
      !====================================
      !Calculates spin based on incremental rotation.
      !input:
      !  drot - incremental rotation (3x3)
      !  dtime - time increment
      !output
      !  W_app - applied spin
      !====================================
      REAL*8 drot(3,3),dtime,W_app(3,3)
      REAL*8 n_dual(3),theta_dot,N(3,3),n_dual_norm
      
      n_dual(1)=drot(3,2)-drot(2,3)!dual vector
      n_dual(2)=-drot(3,1)+drot(1,3)!corr
      n_dual(3)=drot(2,1)-drot(1,2)
      n_dual_norm=sqrt(n_dual(1)**2+n_dual(2)**2+n_dual(3)**2)
      
      if(n_dual_norm.eq.0.0)then
          n_dual=0.0
      else
          n_dual=n_dual/n_dual_norm
      endif
          
      !problem, solution according to Rollet
      if(abs(0.5*(drot(1,1)+drot(2,2)+drot(3,3)-1)).ge.1.0)then
          theta_dot=0.0
      else
          theta_dot=acos(0.5*(drot(1,1)+drot(2,2)+drot(3,3)-1.0))/dtime!angular velocity 
      endif
          
      N=0.0!skew tensor
      N(1,2)=-n_dual(3)
      N(1,3)=n_dual(2)
      N(2,3)=-n_dual(1)
      N=N-transpose(N)
      
      W_app=N*theta_dot!spin
      
      END SUBROUTINE
      
      SUBROUTINE update_config(drot,ng1,ng2)
      
      use flags
      use mvoigt
      use grain_props_v, only : r,wgt
      use grain_state_v, only : etcs,stcs
      use sample_state_v, only : etss,stss,etss_eig
      use mphase_props, only : ngrnph,nphngr,wgt_ph
      use grain_props, only : cr_to_sa
      use phase_transf, only : etcs_pt
      use sample_props_v, only : nph
      
      real, intent(in) :: drot(3,3)
      integer, intent(in) :: ng1,ng2
      
      ! updates variables to current config at t+dt
      !  - update crystal variables
      !  - update sample stress
      
      DIMENSION aux33(3,3),aux66(6,6),aux3333(3,3,3,3)
      
      !update rotation matrix 
      do ng=ng1,ng2
        r(:,:,ng)=MATMUL(r(:,:,ng),transpose(drot(:,:)))
      enddo
          
      !update all variables (or include phase division !???)
      do iph=1,nph
          ng1ph=SUM(nphngr(1:iph))-nphngr(iph)+1
          ng2ph=SUM(nphngr(1:iph))
          wgt_ph(iph)=SUM(wgt(ng1ph:ng2ph))          
          if (ng1ph.le.ng2ph) then
              call cr_to_sa(ng1ph,ng2ph,0)
              call cr_to_sa(ng1ph,ng2ph,1)
          endif
      enddo
      
      !update total stress (crystal)
      do ng=ng1,ng2
        call VOIGT(stcs(:,ng),aux33,aux66,aux3333,1)
        aux33=MATMUL(aux33,transpose(drot(:,:)))
        aux33=MATMUL(drot(:,:),aux33)
        call VOIGT(stcs(:,ng),aux33,aux66,aux3333,2)
      enddo
      
      !update total strain (crystal)
      do ng=ng1,ng2
        call VOIGT(etcs(:,ng),aux33,aux66,aux3333,1)
        aux33=MATMUL(aux33,transpose(drot(:,:)))
        aux33=MATMUL(drot(:,:),aux33)
        call VOIGT(etcs(:,ng),aux33,aux66,aux3333,2)
      enddo        
      
      !rotate phase transformation strain (crystal)
      if (iPhTr.eq.1) then
        do ng=ng1,ng2
          call VOIGT(etcs_pt(:,ng),aux33,aux66,aux3333,1)
          aux33=MATMUL(aux33,transpose(drot(:,:)))
          aux33=MATMUL(drot(:,:),aux33)
          call VOIGT(etcs_pt(:,ng),aux33,aux66,aux3333,2)
        enddo  
      endif
      
      !update total stress (macroscopic)      
      call VOIGT(stss,aux33,aux66,aux3333,1)
      aux33=MATMUL(aux33,transpose(drot))
      aux33=MATMUL(drot,aux33)
      call VOIGT(stss,aux33,aux66,aux3333,2)
      
      !update total strain (macroscopic)
      call VOIGT(etss,aux33,aux66,aux3333,1)
      aux33=MATMUL(aux33,transpose(drot))
      aux33=MATMUL(drot,aux33)
      call VOIGT(etss,aux33,aux66,aux3333,2)   
      
      !update eigenstrain (macroscopic_
      call VOIGT(etss_eig,aux33,aux66,aux3333,1)
      aux33=MATMUL(aux33,transpose(drot))
      aux33=MATMUL(drot,aux33)
      call VOIGT(etss_eig,aux33,aux66,aux3333,2)  
      
      RETURN
      END SUBROUTINE
      
      subroutine read_main_umat(nph,axis,eulerph,nproc,nproc_cycle
     #    ,ncycle,filecrys,filesamp,fileproc,fileprev,filediff,filetemp
     #    ,error_mod,itmax_mod,itmax_grain,label,xmmin_in,outdir
     #    ,lenoutdir,thres_et,nsteps_min)
      use const
      use flags
      
    1 FORMAT(a)      
    2 FORMAT(1h ,a)
    3 FORMAT(1h ,78('*'))
  100 FORMAT(1h ,14('*'),' SELF-CONSISTENT THERMO-ELASTOPLASTIC CODE'
     #,' "EPSC" ',14('*'))
      
      !in
      character(len=150),intent(in) :: outdir
      integer, intent(in) :: lenoutdir
      !out
      integer, intent(out) :: nph,nproc,nsteps_min
      real, intent(out) :: axis(3,NPHM),eulerph(3,NPHM),error_mod(2)
     #    ,xmmin_in,thres_et
      character(len=150), intent(out) :: filecrys(NPHM),filesamp(NPHM)
     #             ,fileproc(NPROCX),fileprev,filediff(NPHM),filetemp(5)
     #             ,label
     
      CHARACTER*78 prosa
      character*150 filename

      !Open and read main control file: "EPSC4.IN"   -> unit=1
      filename = OUTDIR(1:LENOUTDIR)//'epsc4.in'
      OPEN(unit=1,file=filename,status='old')    
      
      read(1,1) label          ! simulation label
      if (iUmat.eq.0) then
          write(*,3)
          write(*,100)
          write(*,2) label
          write(*,3)
      endif
      READ(1,*) nph           
      READ(1,1) prosa
      READ(1,*) ishape          
      do iph=1,nph !MZ_2ph          
          READ(1,*) (axis(i,iph),i=1,3)  !MZ_pseudo read axis and eulerph for each phase
          READ(1,*) (eulerph(i,iph),i=1,3) !initial ellipsoid orientation angles
      enddo
      !Name of TEXTURE file
      read(1,1) prosa
      do iph=1,nph !MZ_2ph       
          read(1,1) filesamp(iph)
      enddo !MZ_2ph 
      READ(1,*) irot            
      !Name of MATERIAL file
      read(1,1) prosa
      do iph=1,nph !MZ_2ph      
          read(1,1) filecrys(iph) 
      enddo !MZ_2ph    
      !Hardening law, backstress MZ_bs      
      read(1,1) prosa !MZ_bs
      read(1,*) kCL,iTwinLaw,iBackStress,iPhTr,iOutput,nCoatedPh
     #    ,nCoatingPh,ivarBC,inonSch,xmmin_in
      !Precision Settings
      read(1,1) prosa     
      read(1,*) itmax_mod,nsteps_min,thres_et
      READ(1,*) error_mod(1)
      read(1,*) itmax_grain
      !Flag for previous procedure and the file name 
      read(1,1) prosa
      read(1,*) i_prev_proc
      read(1,1) fileprev
      READ(1,*) itexskip
      !Flag for difraction calculation and the file with diff directions
      read(1,*) i_diff_dir
      do iph=1,nph   
          read(1,1) filediff(iph)
      enddo
      !Flag for strain pole figure calculation
      read(1,*) i_strpf
      !Number of thermomechanical process in this simulation 
      read(1,1) prosa
      read(1,*) nproc,nproc_cycle,ncycle
      if (nproc.gt.NPROCX) then
        write(*,'(1h ,''ERROR: Number of processes greater''
     #  ,'' than code dimension !!!'',/,1h ,''DIMENSION in code = ''
     #  ,i3)') NPROCX
        write(*,*)
        write(*,'(1h ,''STOP IN MAIN '')')
        stop
      endif
      !Names of PROCESS files 
      read(1,1) prosa
      do n=1,nproc
        read(1,1) fileproc(n)
      enddo
      read(1,1) prosa    
      read(1,1) filetemp(1)
      
      CLOSE(unit=1)         
      
      !umat change file location to include full path
      do iph=1,nph
          filesamp(iph)=outdir(1:lenoutdir)//filesamp(iph)
          filecrys(iph)=outdir(1:lenoutdir)//filecrys(iph)
          filediff(iph)=outdir(1:lenoutdir)//filediff(iph)
      enddo
      do n=1,nproc
          fileproc(n)=outdir(1:lenoutdir)//fileproc(n)
      enddo
      do n=1,1
          filetemp(n) = outdir(1:lenoutdir)//filetemp(n)
      enddo
      
      end subroutine read_main_umat
      
      end module mumat