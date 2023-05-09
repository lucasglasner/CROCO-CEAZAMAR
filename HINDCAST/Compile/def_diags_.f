      subroutine def_diags(ncid,total_rec,ierr)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=150,  MMm0=250,  N=50)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=3,  NP_ETA=20,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=60)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Ntides
      parameter (Ntides=8)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=8*NT)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real ustar(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_ustar/ustar
      integer*4 kbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 kbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_kbl/ kbl
      common /lmd_kpp_hbbl/ hbbl
      common /lmd_kpp_kbbl/ kbbl
      real hbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /lmd_kpp_hbl/ hbls
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV
      parameter (indxU=6, indxV=7)
      integer*4 indxT
      parameter (indxT=indxV+1)
      integer*4 indxS
      parameter (indxS=indxV+ntrc_temp+1)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxTXadv,indxTYadv,indxTVadv,
     &        indxTHmix,indxTVmix,indxTForc,indxTrate
      parameter (indxTXadv=indxV+ntrc_temp+ntrc_salt+ntrc_pas+
     &           ntrc_bio+ntrc_sed+1,
     &           indxTYadv=indxTXadv+NT,
     &           indxTVadv=indxTYadv+NT,
     &           indxTHmix=indxTVadv+NT,
     &           indxTVmix=indxTHmix+NT,
     &           indxTForc=indxTVmix+NT,
     &           indxTrate=indxTForc+NT
     &                                         )
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio
     &                      +ntrc_sed+ntrc_substot
     &           +ntrc_diats+ntrc_diauv+ntrc_diavrt+ntrc_diaek
     &           +ntrc_diapv+ntrc_diaeddy+ntrc_surf+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxAks
      parameter (indxAks=indxAkv+ntrc_temp+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkv+ntrc_temp+5)
      integer*4 indxHbbl
      parameter (indxHbbl=indxAkv+ntrc_temp+6)
      integer*4 indxSSH
      parameter (indxSSH=indxAkv+ntrc_temp+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSUSTR+5)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWSPD,indxTAIR,indxRHUM,indxRADLW,indxRADSW,
     &        indxPRATE,indxUWND,indxVWND,indxPATM
      parameter (indxWSPD=indxSST+3,  indxTAIR=indxSST+4,
     &           indxRHUM=indxSST+5,  indxRADLW=indxSST+6,
     &           indxRADSW=indxSST+7, indxPRATE=indxSST+8,
     &           indxUWND=indxSST+9,  indxVWND=indxSST+10,
     &           indxPATM=indxSST+11)
      integer*4 indxShflx_rlw,indxShflx_lat,indxShflx_sen
      parameter (indxShflx_rlw=indxSST+12,
     &           indxShflx_lat=indxSST+13, indxShflx_sen=indxSST+14)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
     &       , ntbtf(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4 rstHbl
      integer*4 rstHbbl
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
      integer*4 hisT(NT)
      integer*4 nciddia, nrecdia, nrpfdia
     &      , diaTime, diaTime2, diaTstep
     &      , diaTXadv(NT), diaTYadv(NT), diaTVadv(NT)
     &      , diaTHmix(NT), diaTVmix(NT)
     &      , diaTForc(NT), diaTrate(NT)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw, avgBhflx, avgBwflx
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
     &      , avgbvf
      integer*4 avgT(NT)
      integer*4 avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen
      integer*4 nciddia_avg, nrecdia_avg, nrpfdia_avg
     &      , diaTime_avg, diaTime2_avg, diaTstep_avg
     &      , diaTXadv_avg(NT), diaTYadv_avg(NT), diaTVadv_avg(NT)
     &      , diaTHmix_avg(NT), diaTVmix_avg(NT)
     &      , diaTForc_avg(NT), diaTrate_avg(NT)
      logical wrthis(500+NT)
     &      , wrtavg(500+NT)
     &      , wrtdia3D(NT+1)
     &      , wrtdia2D(NT+1)
     &      , wrtdia3D_avg(NT+1)
     &      , wrtdia2D_avg(NT+1)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &     ,  nttclm, ntstf, nttsrc, ntbtf
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
     & ,   rstT
     &      , rstHbl
     &      , rstHbbl
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisShflx_rlw
     &      , hisShflx_lat, hisShflx_sen
     &      , nciddia, nrecdia, nrpfdia
     &      , diaTime, diaTime2, diaTstep
     &      , diaTXadv, diaTYadv, diaTVadv, diaTHmix
     &      , diaTVmix, diaTForc, diaTrate
     &      , nciddia_avg, nrecdia_avg, nrpfdia_avg
     &      , diaTime_avg, diaTime2_avg, diaTstep_avg
     &      , diaTXadv_avg, diaTYadv_avg, diaTVadv_avg
     &      , diaTHmix_avg, diaTVmix_avg, diaTForc_avg
     &      , diaTrate_avg
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgBhflx, avgBwflx
     &      , avgU,    avgV
     &      ,     avgT
     &      ,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , avgbvf
     &      , avgShflx_rlw
     &      , avgShflx_lat, avgShflx_sen
     &      , wrthis
     &      , wrtavg
     &      , wrtdia3D
     &      , wrtdia2D
     &      , wrtdia3D_avg
     &      , wrtdia2D_avg
      character*80 date_str, title, start_date
      character*80 origin_date, start_date_run
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8)             :: origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,  avgname
     &                                ,  dianame
     &                                ,  dianame_avg
     &                                ,   bry_file
      character*75  vname(20, 500)
      common /cncscrum/   date_str,   title,  start_date
     &         ,   origin_date, start_date_run
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,  avgname
     &                                ,  dianame
     &                                ,  dianame_avg
     &                                ,   bry_file
     &                                ,   vname
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      integer*4 nwrtdia
      integer*4 ntsdia_avg, nwrtdia_avg
      logical ldefhis
      logical got_tini(NT)
      logical ldefdia
      logical ldefdia_avg
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , got_tini
     &                      , ldefdia, nwrtdia
     &                      , ldefdia_avg
     &                      , nwrtdia_avg
     &                      , ntsdia_avg
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer*4 mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, 
     &                                                     SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
      real TXadv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TYadv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVadv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real THmix(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVmix(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TForc(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real Trate(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real timedia_avg
      real TXadv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TYadv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVadv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real THmix_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVmix_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TForc_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real Trate_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      common /diag_TXadv/TXadv
     &       /diag_TYadv/TYadv
     &       /diag_TVadv/TVadv
     &       /diag_THmix/THmix
     &       /diag_TVmix/TVmix
     &       /diag_TForc/TForc
     &       /diag_Trate/Trate
      common /diag_timedia_avg/timedia_avg
      common /diag_TXadv_avg/TXadv_avg
     &       /diag_TYadv_avg/TYadv_avg
     &       /diag_TVadv_avg/TVadv_avg
     &       /diag_THmix_avg/THmix_avg
     &       /diag_TVmix_avg/TVmix_avg
     &       /diag_TForc_avg/TForc_avg
     &       /diag_Trate_avg/Trate_avg
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      logical create_new_file, res
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),u2dgrd(3),v2dgrd(3),auxil(2),checkdims
     &      , r3dgrd(4),u3dgrd(4),v3dgrd(4), w3dgrd(4), itrc
      if (may_day_flag.ne.0) return
      ierr=0
      lstr=lenstr(dianame)
      if (nrpfdia.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfdia))
        call insert_time_index (dianame, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefdia
      if (ncid.ne.-1) create_new_file=.false.
      if (mynode.gt.0) create_new_file=.false.
 10   if (create_new_file) then
        lstr=lenstr(dianame)
        ierr=nf_create(dianame(1:lstr),nf_64bit_offset, ncid)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) dianame(1:lstr)
          may_day_flag=3
          return
        endif
        call put_global_atts (ncid, ierr)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) dianame(1:lstr)
          may_day_flag=3
          return
        endif
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
      if (total_rec.le.1) then
         call def_grid_3d(ncid, r2dgrd, u2dgrd, v2dgrd
     &        ,r3dgrd, w3dgrd)
      endif
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 diaTstep)
        ierr=nf_put_att_text (ncid, diaTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, diaTime)
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, diaTime, 'long_name',
     &       lvar, vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, diaTime, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, diaTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, diaTime, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, diaTime, indxTime, 5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, diaTime2)
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2, 'long_name',
     &       lvar, vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, diaTime2, indxTime2, 5,
     &       NF_REAL, ierr)
        do itrc=1,NT
          if (wrtdia3D(itrc)) then
          lvar=lenstr(vname(1,indxTXadv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxTXadv+itrc-1)(1:lvar),
     &          NF_REAL, 4, r3dgrd, diaTXadv(itrc))
          lvar=lenstr(vname(2,indxTXadv+itrc-1))
          ierr=nf_put_att_text (ncid, diaTXadv(itrc), 'long_name',
     &          lvar, vname(2,indxTXadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(3,indxTXadv+itrc-1))
          ierr=nf_put_att_text (ncid, diaTXadv(itrc), 'units', lvar,
     &          vname(3,indxTXadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxTXadv+itrc-1))
          ierr=nf_put_att_text (ncid,diaTXadv(itrc), 'field',
     &          lvar, vname(4,indxTXadv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTXadv(itrc),indxTXadv+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTYadv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTYadv+itrc-1)(1:lvar),
     &                           NF_REAL, 4, r3dgrd, diaTYadv(itrc))
           lvar=lenstr(vname(2,indxTYadv+itrc-1))
           ierr=nf_put_att_text (ncid, diaTYadv(itrc), 'long_name',
     &                     lvar, vname(2,indxTYadv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxTYadv+itrc-1))
           ierr=nf_put_att_text (ncid, diaTYadv(itrc), 'units', lvar,
     &                             vname(3,indxTYadv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTYadv+itrc-1))
           ierr=nf_put_att_text (ncid,diaTYadv(itrc), 'field',
     &                       lvar, vname(4,indxTYadv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTYadv(itrc),indxTYadv+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTVadv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTVadv+itrc-1)(1:lvar),
     &                           NF_REAL, 4, r3dgrd, diaTVadv(itrc))
           lvar=lenstr(vname(2,indxTVadv+itrc-1))
           ierr=nf_put_att_text (ncid, diaTVadv(itrc), 'long_name',
     &                     lvar, vname(2,indxTVadv+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxTVadv+itrc-1))
           ierr=nf_put_att_text (ncid, diaTVadv(itrc), 'units', lvar,
     &                             vname(3,indxTVadv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTVadv+itrc-1))
           ierr=nf_put_att_text (ncid,diaTVadv(itrc), 'field',
     &                       lvar, vname(4,indxTVadv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTVadv(itrc),indxTVadv+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTHmix+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTHmix+itrc-1)(1:lvar),
     &                            NF_REAL, 4, r3dgrd, diaTHmix(itrc))
           lvar=lenstr(vname(2,indxTHmix+itrc-1))
           ierr=nf_put_att_text (ncid, diaTHmix(itrc), 'long_name',
     &                       lvar, vname(2,indxTHmix+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxTHmix+itrc-1))
           ierr=nf_put_att_text (ncid, diaTHmix(itrc), 'units', lvar,
     &                             vname(3,indxTHmix+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTHmix+itrc-1))
           ierr=nf_put_att_text (ncid,diaTHmix(itrc), 'field',
     &                       lvar, vname(4,indxTHmix+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTHmix(itrc),indxTHmix+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTVmix+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTVmix+itrc-1)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, diaTVmix(itrc))
           lvar=lenstr(vname(2,indxTVmix+itrc-1))
           ierr=nf_put_att_text (ncid, diaTVmix(itrc), 'long_name',
     &                     lvar, vname(2,indxTVmix+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxTVmix+itrc-1))
           ierr=nf_put_att_text (ncid, diaTVmix(itrc), 'units', lvar,
     &                             vname(3,indxTVmix+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTVmix+itrc-1))
           ierr=nf_put_att_text (ncid,diaTVmix(itrc), 'field',
     &                       lvar, vname(4,indxTVmix+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTVmix(itrc),indxTVmix+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTForc+itrc-1))
           ierr=nf_def_var (ncid,vname(1,indxTForc+itrc-1)(1:lvar),
     &                        NF_REAL, 4, r3dgrd, diaTForc(itrc))
           lvar=lenstr(vname(2,indxTForc+itrc-1))
           ierr=nf_put_att_text (ncid, diaTForc(itrc), 'long_name',
     &                      lvar, vname(2,indxTForc+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxTForc+itrc-1))
           ierr=nf_put_att_text (ncid, diaTForc(itrc),'units',
     &                     lvar, vname(3,indxTForc+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTForc+itrc-1))
           ierr=nf_put_att_text (ncid, diaTForc(itrc),'field',
     &                     lvar, vname(4,indxTForc+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTForc(itrc),indxTForc+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTrate+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTrate+itrc-1)(1:lvar),
     &                           NF_REAL, 4, r3dgrd, diaTrate(itrc))
           lvar=lenstr(vname(2,indxTrate+itrc-1))
           ierr=nf_put_att_text (ncid, diaTrate(itrc), 'long_name',
     &                     lvar, vname(2,indxTrate+itrc-1)(1:lvar))
           lvar=lenstr(vname(3,indxTrate+itrc-1))
           ierr=nf_put_att_text (ncid, diaTrate(itrc),'units',
     &                     lvar, vname(3,indxTrate+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTrate+itrc-1))
           ierr=nf_put_att_text (ncid, diaTrate(itrc),'field',
     &                     lvar, vname(4,indxTrate+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTrate(itrc),indxTrate+itrc-1,5,
     &       NF_REAL, ierr)
         endif
       enddo
        ierr=nf_enddef(ncid)
        if (mynode.eq.0) write(stdout,'(6x,4A,1x,A,i4)')
     &        'DEF_DIAG - Created ',
     &                'new netCDF file ''',
     &                 dianame(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open (dianame(1:lstr), nf_write, ncid)
        if (ierr. ne. nf_noerr) then
          ierr=checkdims (ncid, dianame, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfdia.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfdia))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &     ) 'DEF_DIAGS/DIAGS_AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  dianame(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfdia.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
      if (mynode.eq.0)  write(stdout,'(/1x,4A,2x,A,I4/)')
     &         'DEF_DIAGS_HIS/AVG ERROR: ',
     &         'Cannot open file ''', dianame(1:lstr), '''.'
     &                   ,' mynode =', mynode
            goto 99
          endif
        endif
        ierr=nf_inq_varid (ncid, 'time_step', diaTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', dianame(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),diaTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), dianame(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),diaTime2)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), dianame(1:lstr)
          goto 99
        endif
        do itrc=1,NT
          if (wrtdia3D(itrc)) then
          lvar=lenstr(vname(1,indxTXadv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTXadv+itrc-1)(1:lvar),
     &                                               diaTXadv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTXadv+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTYadv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTYadv+itrc-1)(1:lvar),
     &                                               diaTYadv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTYadv+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTVadv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTVadv+itrc-1)(1:lvar),
     &                                               diaTVadv(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTVadv+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTHmix+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTHmix+itrc-1)(1:lvar),
     &                                               diaTHmix(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTHmix+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTVmix+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTVmix+itrc-1)(1:lvar),
     &                                               diaTVmix(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTVmix+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTForc+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTForc+itrc-1)(1:lvar),
     &                                               diaTForc(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTForc+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTrate+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTrate+itrc-1)(1:lvar),
     &                                               diaTrate(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTrate+itrc-1)(1:lvar),
     &                                         dianame(1:lstr)
            goto 99
          endif
        endif
      enddo
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_DIAG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_DIAG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      else
        ierr=nf_open (dianame(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
        if (mynode.eq.0)  write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'DEF_DIAG ERROR: ',
     &                'Cannot open file ''', dianame(1:lstr), '''.'
     &                 ,' mynode =', mynode
          goto 99
        endif
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
       write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_DIAG ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'DEF_DIAG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
      if (total_rec.le.1) call wrt_grid (ncid, dianame, lstr)
  11  format(/' DEF_DIAGS - unable to create diag file: ',a)
  20  format(/' DEF_DIAGS - error while writing variable: ',a,
     &        /,15x,'into diag  file: ',a)
  99  return
      end
      subroutine def_diags_avg(ncid,total_rec,ierr)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=150,  MMm0=250,  N=50)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=3,  NP_ETA=20,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=60)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Ntides
      parameter (Ntides=8)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=8*NT)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real ustar(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_ustar/ustar
      integer*4 kbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 kbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_kbl/ kbl
      common /lmd_kpp_hbbl/ hbbl
      common /lmd_kpp_kbbl/ kbbl
      real hbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /lmd_kpp_hbl/ hbls
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV
      parameter (indxU=6, indxV=7)
      integer*4 indxT
      parameter (indxT=indxV+1)
      integer*4 indxS
      parameter (indxS=indxV+ntrc_temp+1)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxTXadv,indxTYadv,indxTVadv,
     &        indxTHmix,indxTVmix,indxTForc,indxTrate
      parameter (indxTXadv=indxV+ntrc_temp+ntrc_salt+ntrc_pas+
     &           ntrc_bio+ntrc_sed+1,
     &           indxTYadv=indxTXadv+NT,
     &           indxTVadv=indxTYadv+NT,
     &           indxTHmix=indxTVadv+NT,
     &           indxTVmix=indxTHmix+NT,
     &           indxTForc=indxTVmix+NT,
     &           indxTrate=indxTForc+NT
     &                                         )
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio
     &                      +ntrc_sed+ntrc_substot
     &           +ntrc_diats+ntrc_diauv+ntrc_diavrt+ntrc_diaek
     &           +ntrc_diapv+ntrc_diaeddy+ntrc_surf+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxAks
      parameter (indxAks=indxAkv+ntrc_temp+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkv+ntrc_temp+5)
      integer*4 indxHbbl
      parameter (indxHbbl=indxAkv+ntrc_temp+6)
      integer*4 indxSSH
      parameter (indxSSH=indxAkv+ntrc_temp+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSUSTR+5)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWSPD,indxTAIR,indxRHUM,indxRADLW,indxRADSW,
     &        indxPRATE,indxUWND,indxVWND,indxPATM
      parameter (indxWSPD=indxSST+3,  indxTAIR=indxSST+4,
     &           indxRHUM=indxSST+5,  indxRADLW=indxSST+6,
     &           indxRADSW=indxSST+7, indxPRATE=indxSST+8,
     &           indxUWND=indxSST+9,  indxVWND=indxSST+10,
     &           indxPATM=indxSST+11)
      integer*4 indxShflx_rlw,indxShflx_lat,indxShflx_sen
      parameter (indxShflx_rlw=indxSST+12,
     &           indxShflx_lat=indxSST+13, indxShflx_sen=indxSST+14)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
     &       , ntbtf(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4 rstHbl
      integer*4 rstHbbl
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
      integer*4 hisT(NT)
      integer*4 nciddia, nrecdia, nrpfdia
     &      , diaTime, diaTime2, diaTstep
     &      , diaTXadv(NT), diaTYadv(NT), diaTVadv(NT)
     &      , diaTHmix(NT), diaTVmix(NT)
     &      , diaTForc(NT), diaTrate(NT)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw, avgBhflx, avgBwflx
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
     &      , avgbvf
      integer*4 avgT(NT)
      integer*4 avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen
      integer*4 nciddia_avg, nrecdia_avg, nrpfdia_avg
     &      , diaTime_avg, diaTime2_avg, diaTstep_avg
     &      , diaTXadv_avg(NT), diaTYadv_avg(NT), diaTVadv_avg(NT)
     &      , diaTHmix_avg(NT), diaTVmix_avg(NT)
     &      , diaTForc_avg(NT), diaTrate_avg(NT)
      logical wrthis(500+NT)
     &      , wrtavg(500+NT)
     &      , wrtdia3D(NT+1)
     &      , wrtdia2D(NT+1)
     &      , wrtdia3D_avg(NT+1)
     &      , wrtdia2D_avg(NT+1)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &     ,  nttclm, ntstf, nttsrc, ntbtf
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
     & ,   rstT
     &      , rstHbl
     &      , rstHbbl
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisShflx_rlw
     &      , hisShflx_lat, hisShflx_sen
     &      , nciddia, nrecdia, nrpfdia
     &      , diaTime, diaTime2, diaTstep
     &      , diaTXadv, diaTYadv, diaTVadv, diaTHmix
     &      , diaTVmix, diaTForc, diaTrate
     &      , nciddia_avg, nrecdia_avg, nrpfdia_avg
     &      , diaTime_avg, diaTime2_avg, diaTstep_avg
     &      , diaTXadv_avg, diaTYadv_avg, diaTVadv_avg
     &      , diaTHmix_avg, diaTVmix_avg, diaTForc_avg
     &      , diaTrate_avg
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgBhflx, avgBwflx
     &      , avgU,    avgV
     &      ,     avgT
     &      ,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , avgbvf
     &      , avgShflx_rlw
     &      , avgShflx_lat, avgShflx_sen
     &      , wrthis
     &      , wrtavg
     &      , wrtdia3D
     &      , wrtdia2D
     &      , wrtdia3D_avg
     &      , wrtdia2D_avg
      character*80 date_str, title, start_date
      character*80 origin_date, start_date_run
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8)             :: origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,  avgname
     &                                ,  dianame
     &                                ,  dianame_avg
     &                                ,   bry_file
      character*75  vname(20, 500)
      common /cncscrum/   date_str,   title,  start_date
     &         ,   origin_date, start_date_run
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,  avgname
     &                                ,  dianame
     &                                ,  dianame_avg
     &                                ,   bry_file
     &                                ,   vname
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      integer*4 nwrtdia
      integer*4 ntsdia_avg, nwrtdia_avg
      logical ldefhis
      logical got_tini(NT)
      logical ldefdia
      logical ldefdia_avg
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , got_tini
     &                      , ldefdia, nwrtdia
     &                      , ldefdia_avg
     &                      , nwrtdia_avg
     &                      , ntsdia_avg
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer*4 mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, 
     &                                                     SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      integer*4 max_opt_size
      parameter (max_opt_size=4400)
      character*4400 Coptions,srcs
      common /strings/ Coptions,srcs
      real TXadv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TYadv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVadv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real THmix(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVmix(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TForc(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real Trate(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real timedia_avg
      real TXadv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TYadv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVadv_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real THmix_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TVmix_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real TForc_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      real Trate_avg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,NT)
      common /diag_TXadv/TXadv
     &       /diag_TYadv/TYadv
     &       /diag_TVadv/TVadv
     &       /diag_THmix/THmix
     &       /diag_TVmix/TVmix
     &       /diag_TForc/TForc
     &       /diag_Trate/Trate
      common /diag_timedia_avg/timedia_avg
      common /diag_TXadv_avg/TXadv_avg
     &       /diag_TYadv_avg/TYadv_avg
     &       /diag_TVadv_avg/TVadv_avg
     &       /diag_THmix_avg/THmix_avg
     &       /diag_TVmix_avg/TVmix_avg
     &       /diag_TForc_avg/TForc_avg
     &       /diag_Trate_avg/Trate_avg
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      logical create_new_file, res
      integer*4 ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),u2dgrd(3),v2dgrd(3),auxil(2),checkdims
     &      , r3dgrd(4),u3dgrd(4),v3dgrd(4), w3dgrd(4), itrc
        character*60 text
      if (may_day_flag.ne.0) return
      ierr=0
      lstr=lenstr(dianame_avg)
      if (nrpfdia_avg.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, nrpfdia_avg))
        call insert_time_index (dianame_avg, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif
      create_new_file=ldefdia_avg
      if (ncid.ne.-1) create_new_file=.false.
      if (mynode.gt.0) create_new_file=.false.
 10   if (create_new_file) then
        lstr=lenstr(dianame_avg)
        ierr=nf_create(dianame_avg(1:lstr),nf_64bit_offset, ncid)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) dianame_avg(1:lstr)
          may_day_flag=3
          return
        endif
        call put_global_atts (ncid, ierr)
        if (ierr.ne.nf_noerr) then
          write(stdout,11) dianame_avg(1:lstr)
          may_day_flag=3
          return
        endif
        ierr=nf_def_dim (ncid, 'xi_rho',   xi_rho,   r2dgrd(1))
        ierr=nf_def_dim (ncid, 'xi_u',     xi_u,     u2dgrd(1))
        ierr=nf_def_dim (ncid, 'eta_rho',  eta_rho,  r2dgrd(2))
        ierr=nf_def_dim (ncid, 'eta_v',    eta_v,    v2dgrd(2))
        ierr=nf_def_dim (ncid, 's_rho',    N,        r3dgrd(3))
        ierr=nf_def_dim (ncid, 's_w',      N+1,      w3dgrd(3))
        ierr=nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=nf_def_dim (ncid, 'auxil',    4,        auxil(1))
        auxil(2)=timedim
        r2dgrd(3)=timedim
        u2dgrd(2)=r2dgrd(2)
        u2dgrd(3)=timedim
        v2dgrd(1)=r2dgrd(1)
        v2dgrd(3)=timedim
        r3dgrd(1)=r2dgrd(1)
        r3dgrd(2)=r2dgrd(2)
        r3dgrd(4)=timedim
        u3dgrd(1)=u2dgrd(1)
        u3dgrd(2)=r2dgrd(2)
        u3dgrd(3)=r3dgrd(3)
        u3dgrd(4)=timedim
        v3dgrd(1)=r2dgrd(1)
        v3dgrd(2)=v2dgrd(2)
        v3dgrd(3)=r3dgrd(3)
        v3dgrd(4)=timedim
      if (total_rec.le.1) then
         call def_grid_3d(ncid, r2dgrd, u2dgrd, v2dgrd
     &        ,r3dgrd, w3dgrd)
      endif
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 diaTstep_avg)
        ierr=nf_put_att_text (ncid, diaTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, diaTime_avg)
        text='avg'/ /vname(2,indxTime)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, diaTime_avg, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(2,indxTime))
        ierr=nf_put_att_text (ncid, diaTime_avg, 'long_name', lvar,
     &                                vname(2,indxTime)(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, diaTime_avg, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))
        lvar=lenstr(vname(4,indxTime))
        ierr=nf_put_att_text (ncid, diaTime_avg, 'field',  lvar,
     &                                vname(4,indxTime)(1:lvar))
        call nf_add_attribute(ncid, diaTime_avg, indxTime, 5,
     &       NF_REAL, ierr)
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_def_var (ncid, vname(1,indxTime2)(1:lvar),
     &                            NF_DOUBLE, 1, timedim, diaTime2_avg)
        text='avg'/ /vname(2,indxTime2)
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, diaTime2_avg, 'long_name',
     &       lvar, text(1:lvar))
        lvar=lenstr(vname(2,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2_avg, 'long_name', lvar,
     &                                vname(2,indxTime2)(1:lvar))
        lvar=lenstr(vname(3,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2_avg, 'units',  lvar,
     &                                vname(3,indxTime2)(1:lvar))
        lvar=lenstr(vname(4,indxTime2))
        ierr=nf_put_att_text (ncid, diaTime2_avg, 'field',  lvar,
     &                                vname(4,indxTime2)(1:lvar))
        call nf_add_attribute(ncid, diaTime2_avg, indxTime2, 5,
     &       NF_REAL, ierr)
        do itrc=1,NT
          if (wrtdia3D_avg(itrc)) then
          lvar=lenstr(vname(1,indxTXadv+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxTXadv+itrc-1)(1:lvar),
     &          NF_REAL, 4, r3dgrd, diaTXadv_avg(itrc))
          text='averaged '/ /vname(2,indxTXadv+itrc-1)
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, diaTXadv_avg(itrc), 'long_name',
     &          lvar, text(1:lvar))
          lvar=lenstr(vname(3,indxTXadv+itrc-1))
          ierr=nf_put_att_text (ncid, diaTXadv_avg(itrc), 'units', lvar,
     &          vname(3,indxTXadv+itrc-1)(1:lvar))
          lvar=lenstr(vname(4,indxTXadv+itrc-1))
          ierr=nf_put_att_text (ncid,diaTXadv_avg(itrc), 'field',
     &          lvar, vname(4,indxTXadv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTXadv_avg(itrc),indxTXadv+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTYadv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTYadv+itrc-1)(1:lvar),
     &                           NF_REAL, 4, r3dgrd, diaTYadv_avg(itrc))
           text='averaged '/ /vname(2,indxTYadv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diaTYadv_avg(itrc), 'long_name',
     &       lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxTYadv+itrc-1))
           ierr=nf_put_att_text (ncid, diaTYadv_avg(itrc), 'units', 
     &                                                             lvar,
     &                             vname(3,indxTYadv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTYadv+itrc-1))
           ierr=nf_put_att_text (ncid,diaTYadv_avg(itrc), 'field',
     &                       lvar, vname(4,indxTYadv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTYadv_avg(itrc),indxTYadv+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTVadv+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTVadv+itrc-1)(1:lvar),
     &                           NF_REAL, 4, r3dgrd, diaTVadv_avg(itrc))
           text='averaged '/ /vname(2,indxTVadv+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diaTVadv_avg(itrc), 'long_name',
     &       lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxTVadv+itrc-1))
           ierr=nf_put_att_text (ncid, diaTVadv_avg(itrc), 'units', 
     &                                                             lvar,
     &                             vname(3,indxTVadv+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTVadv+itrc-1))
           ierr=nf_put_att_text (ncid,diaTVadv_avg(itrc), 'field',
     &                       lvar, vname(4,indxTVadv+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTVadv_avg(itrc),indxTVadv+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTHmix+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTHmix+itrc-1)(1:lvar),
     &                            NF_REAL, 4, r3dgrd, 
     &                                               diaTHmix_avg(itrc))
           text='averaged '/ /vname(2,indxTHmix+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diaTHmix_avg(itrc), 'long_name',
     &                                            lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxTHmix+itrc-1))
           ierr=nf_put_att_text (ncid, diaTHmix_avg(itrc), 'units', 
     &                                                             lvar,
     &                             vname(3,indxTHmix+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTHmix+itrc-1))
           ierr=nf_put_att_text (ncid,diaTHmix_avg(itrc), 'field',
     &                       lvar, vname(4,indxTHmix+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTHmix_avg(itrc),indxTHmix+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTVmix+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTVmix+itrc-1)(1:lvar),
     &                             NF_REAL, 4, r3dgrd, 
     &                                               diaTVmix_avg(itrc))
           text='averaged '/ /vname(2,indxTVmix+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diaTVmix_avg(itrc), 'long_name',
     &       lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxTVmix+itrc-1))
           ierr=nf_put_att_text (ncid, diaTVmix_avg(itrc), 'units', 
     &                                                             lvar,
     &                             vname(3,indxTVmix+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTVmix+itrc-1))
           ierr=nf_put_att_text (ncid,diaTVmix_avg(itrc), 'field',
     &                       lvar, vname(4,indxTVmix+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTVmix_avg(itrc),indxTVmix+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTForc+itrc-1))
           ierr=nf_def_var (ncid,vname(1,indxTForc+itrc-1)(1:lvar),
     &                        NF_REAL, 4, r3dgrd, diaTForc_avg(itrc))
           text='averaged '/ /vname(2,indxTForc+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diaTForc_avg(itrc), 'long_name',
     &          lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxTForc+itrc-1))
           ierr=nf_put_att_text (ncid, diaTForc_avg(itrc),'units',
     &                     lvar, vname(3,indxTForc+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTForc+itrc-1))
           ierr=nf_put_att_text (ncid, diaTForc_avg(itrc),'field',
     &                     lvar, vname(4,indxTForc+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTForc_avg(itrc),indxTForc+itrc-1,5,
     &       NF_REAL, ierr)
           lvar=lenstr(vname(1,indxTrate+itrc-1))
           ierr=nf_def_var (ncid, vname(1,indxTrate+itrc-1)(1:lvar),
     &                           NF_REAL, 4, r3dgrd, diaTrate_avg(itrc))
           text='averaged '/ /vname(2,indxTrate+itrc-1)
           lvar=lenstr(text)
           ierr=nf_put_att_text (ncid, diaTrate_avg(itrc), 'long_name',
     &                                          lvar, text(1:lvar))
           lvar=lenstr(vname(3,indxTrate+itrc-1))
           ierr=nf_put_att_text (ncid, diaTrate_avg(itrc),'units',
     &                     lvar, vname(3,indxTrate+itrc-1)(1:lvar))
           lvar=lenstr(vname(4,indxTrate+itrc-1))
           ierr=nf_put_att_text (ncid, diaTrate_avg(itrc),'field',
     &                     lvar, vname(4,indxTrate+itrc-1)(1:lvar))
       call nf_add_attribute(ncid,diaTrate_avg(itrc),indxTrate+itrc-1,5,
     &       NF_REAL, ierr)
         endif
       enddo
        ierr=nf_enddef(ncid)
        if (mynode.eq.0) write(stdout,'(6x,4A,1x,A,i4)')
     &        'DEF_DIAG_AVG - Created ',
     &                'new netCDF file ''',
     &                 dianame_avg(1:lstr), '''.'
     &                 ,' mynode =', mynode
      elseif (ncid.eq.-1) then
        ierr=nf_open (dianame_avg(1:lstr), nf_write, ncid)
        if (ierr. ne. nf_noerr) then
          ierr=checkdims (ncid, dianame_avg, lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (nrpfdia_avg.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, nrpfdia_avg))
            endif
            if (ierr.gt.0) then
              if (mynode.eq.0) write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &     ) 'DEF_DIAGS/DIAGS_AVG WARNING: Actual number of records',
     &               rec,  'in netCDF file',  '''',  
     &                                              dianame_avg(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (nrpfdia_avg.eq.0) then
              total_rec=rec+1
              if (mynode.gt.0) total_rec=total_rec-1
            endif
            ierr=nf_noerr
          endif
        endif
        if (ierr. ne. nf_noerr) then
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
      if (mynode.eq.0)  write(stdout,'(/1x,4A,2x,A,I4/)')
     &         'DEF_DIAGS_HIS/AVG ERROR: ',
     &         'Cannot open file ''', dianame_avg(1:lstr), '''.'
     &                   ,' mynode =', mynode
            goto 99
          endif
        endif
        ierr=nf_inq_varid (ncid, 'time_step', diaTstep_avg)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', dianame_avg(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),diaTime_avg)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), dianame_avg(1:lstr)
          goto 99
        endif
        lvar=lenstr(vname(1,indxTime2))
        ierr=nf_inq_varid (ncid,vname(1,indxTime2)(1:lvar),diaTime2_avg)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime2)(1:lvar), 
     &                                               dianame_avg(1:lstr)
          goto 99
        endif
        do itrc=1,NT
          if (wrtdia3D_avg(itrc)) then
          lvar=lenstr(vname(1,indxTXadv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTXadv+itrc-1)(1:lvar),
     &                                               diaTXadv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTXadv+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTYadv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTYadv+itrc-1)(1:lvar),
     &                                               diaTYadv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTYadv+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTVadv+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTVadv+itrc-1)(1:lvar),
     &                                               diaTVadv_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTVadv+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTHmix+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTHmix+itrc-1)(1:lvar),
     &                                               diaTHmix_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTHmix+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTVmix+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTVmix+itrc-1)(1:lvar),
     &                                               diaTVmix_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTVmix+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTForc+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTForc+itrc-1)(1:lvar),
     &                                               diaTForc_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTForc+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
          lvar=lenstr(vname(1,indxTrate+itrc-1))
          ierr=nf_inq_varid (ncid, vname(1,indxTrate+itrc-1)(1:lvar),
     &                                               diaTrate_avg(itrc))
          if (ierr .ne. nf_noerr) then
            write(stdout,1) vname(1,indxTrate+itrc-1)(1:lvar),
     &                                         dianame_avg(1:lstr)
            goto 99
          endif
        endif
      enddo
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_DIAG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
        if (mynode.eq.0) write(*,'(6x,2A,i4,1x,A,i4)')
     &                     'DEF_DIAG -- Opened ',
     &                     'existing file  from record =', rec
     &                      ,' mynode =', mynode
      else
        ierr=nf_open (dianame_avg(1:lstr), nf_write, ncid)
        if (ierr .ne. nf_noerr) then
        if (mynode.eq.0)  write(stdout,'(/1x,4A,2x,A,I4/)')
     &                'DEF_DIAG ERROR: ',
     &                'Cannot open file ''', dianame_avg(1:lstr), '''.'
     &                 ,' mynode =', mynode
          goto 99
        endif
      endif
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
       write(*,'(6x,2A,i4,1x,A,i4)') 'DEF_DIAG ERROR: Cannot ',
     &    'switch to ''nf_nofill'' more; netCDF error code =', ierr
      endif
   1  format(/1x,'DEF_DIAG ERROR: Cannot find variable ''',
     &                   A, ''' in netCDF file ''', A, '''.'/)
      if (total_rec.le.1) call wrt_grid (ncid, dianame_avg, lstr)
  11  format(/' DEF_DIAGS - unable to create diag file: ',a)
  20  format(/' DEF_DIAGS - error while writing variable: ',a,
     &        /,15x,'into diag  file: ',a)
  99  return
      end
