program main

   use omp_lib  !OpenMP使用

   !use linear_equation_of_state
   !use nonlinear_equation_of_state

   implicit none

   character(8)  :: date
   character(10) :: time
   character(5)  :: zone
   integer, dimension(8) :: values

   character(15):: &
      expid='ex4-case1' ! 実験名

   integer,parameter:: &
      im = 100, & ! x方向格子数
      jm =  40, & ! y方向格子数
      km =  12   ! z方向格子数
       
   real(8),parameter:: &
      t0_npm2  = -0.1d0, &                               ! 東西方向風応力
      gr_mps2  = 9.80d0, &                               ! 重力加速度
      r0_kgpm3 = 1.d3, &                                 ! 基準海水密度
      f0_psec  = 0.d-4, &                                ! コリオリ係数
      bt_pmsec = 1.d-11, &                               ! ベータ
      er_m     = 6370.d3,&                               ! 地球半径
      om_psec  = 7.29d-5,&                               ! 地球自転速度
      ah_m2ps  = 1.d5, &                                 ! 水平渦粘性係数
      av_m2ps  = 1.d-3, &                                ! 鉛直渦粘性係数
      kh_m2ps  = 0.d0, &                                 ! 水位の水平拡散係数
      kv_m2ps  = 1.d-4,  &                               ! 鉛直渦拡散係数
      gm_psec  = 1.d0/50.d0/86400.d0, &                  ! 温位緩和係数
      alphat_pk  = 2.0d-4, &                             ! 熱膨張係数
      alphas_ppsu = 7.6d-4, &                            ! 塩分収縮係数
      t0_c = 30.0d0, &                                   ! 基準温度
      dx_m     = 100.d3, &                               ! x方向格子幅
      dy_m     = 100.d3, &                               ! y方向格子幅
      dz_m     = 50.d0, &                                ! z方向格子幅
      dt_sec   = 300.d0, &                               ! 時間刻み幅
      time_to_start_sec   = 0.d0, &                      ! 実験開始時刻
      time_to_end_sec     = 86400.d0*30.d0*12.d0*5.d0, & ! 実験終了時刻（５年）
      !time_to_end_sec     = dt_sec * 10.d0 * 5.d0, &    ! 実験終了時刻（１年）
      output_interval_sec = 86400.d0*30.d0, &            ! 出力時間間隔（３０日）
      !output_interval_sec = dt_sec , &                  ! 出力時間間隔（１０日）
      asf      = 0.5d0                                   ! アセリンフィルター係数
   real(8),parameter :: fric_non = 0.0d0 / dz_m  !- linear friction coefficient (for Stommel 1948)
       
   real(8),dimension(:,:,:):: &
      ua_mps(0:im+1,0:jm+1,0:km+1), & ! 東西流速
      ub_mps(0:im+1,0:jm+1,0:km+1), &
      uc_mps(0:im+1,0:jm+1,0:km+1), &
      va_mps(0:im+1,0:jm+1,0:km+1), & ! 南北流速
      vb_mps(0:im+1,0:jm+1,0:km+1), &
      vc_mps(0:im+1,0:jm+1,0:km+1), &
      wa_mps(0:im+1,0:jm+1,0:km+1), & ! 鉛直流速
      wb_mps(0:im+1,0:jm+1,0:km+1), &
      wc_mps(0:im+1,0:jm+1,0:km+1), &
      ta_c(0:im+1,0:jm+1,0:km+1),   & ! 温位
      tb_c(0:im+1,0:jm+1,0:km+1),   &
      tc_c(0:im+1,0:jm+1,0:km+1),   &
      sa_psu(0:im+1,0:jm+1,0:km+1), & ! 塩分 （34.0からの偏差）
      sb_psu(0:im+1,0:jm+1,0:km+1), &
      sc_psu(0:im+1,0:jm+1,0:km+1)

   real(8),dimension(:,:,:):: &
      ww_mps(0:im+1,0:jm+1,0:km+1), &     ! 鉛直流速
      pp_npm2(0:im+1,0:jm+1,0:km+1), &    ! 圧力
      rr_kgpm3(0:im+1,0:jm+1,0:km+1)      ! 密度

   real(8),dimension(:,:):: &
      ea_m(0:im+1,0:jm+1), &       ! 水位
      eb_m(0:im+1,0:jm+1), &
      ec_m(0:im+1,0:jm+1), &
      tx_npm2(0:im+1,0:jm+1), &    ! 風応力（東西）
      ty_npm2(0:im+1,0:jm+1), &    ! 風応力（南北）
      at_c(0:im+1,0:jm+1),      &  ! 緩和温位（気温）
      sf_psumps(0:im+1,0:jm+1)     ! 海面塩分フラックス [psu m/s]

   real(8),dimension(:):: &
      fs_psec(0:jm+1) ! コリオリ係数

   real(8),dimension(:):: &
      ys(0:jm), & ! 計算用ｙ
      zs(0:km)    ! 計算用z

   real(8):: &
      time_sec, & ! 時間
      time_to_output_sec, & ! 次期出力時刻
      dpi, & ! π
      gu, gv, ge, gt, gs, &
      pre, adx, ady, adz, adv_x, adv_y, adv_z, cor, dfx, dfy, dfz, frc, rr_f, &
      ax, ay, az, div, kx, ky, &
      tav, sav ! 作業用

   ! 強制力（貿易風＋偏西風：NECC/EUC向け）
   real(8) :: lat_deg
   real(8), parameter :: pi = 4.d0*datan(1.d0)

   ! 風応力パラメータ
   real(8), parameter :: tau_tr = 0.10d0   ! 貿易風振幅
   real(8), parameter :: tau_w  = 0.06d0   ! 偏西風振幅
   real(8), parameter :: lat_tr = 8.d0     ! 貿易風中心
   real(8), parameter :: wid_tr = 6.d0     ! 貿易風幅
   real(8), parameter :: lat_w  = 30.d0    ! 偏西風中心
   real(8), parameter :: wid_w  = 6.d0     ! 偏西風幅

   real(8):: &
      er_rec, pi_rec, wid_tr_rec, wid_w_rec, &
      dx_rec, dy_rec, dz_rec, &
      dx_2_rec, dy_2_rec, dz_2_rec, &
      r0_rec, av_rec, kv_rec, &
      lat_max_deg, rec_600

   integer:: &
      i,j,k,n ! 作業用

   character(80):: &
      buff ! 作業用

   ! プログラム開始時の時刻を表示
   call date_and_time(date, time, zone, values)
   print *, '========================================='
   print *, 'プログラム開始'
   print *, '日付: ', date(1:4), '/', date(5:6), '/', date(7:8)
   print *, '時刻: ', time(1:2), ':', time(3:4), ':', time(5:10)
   print *, '========================================='
   print *

   !<--- ★ OpenMP情報表示（追加）
   print *, 'OpenMP max threads: ', omp_get_max_threads()
   print *

   er_rec = 1.d0/er_m
   pi_rec = 1.d0/pi
   wid_tr_rec = 1.d0/wid_tr
   wid_w_rec = 1.d0/wid_w
   dx_rec = 1.d0/dx_m
   dy_rec = 1.d0/dy_m
   dz_rec = 1.d0/dz_m
   dx_2_rec = 1.d0/dx_m**2
   dy_2_rec = 1.d0/dy_m**2
   dz_2_rec = 1.d0/dz_m**2
   r0_rec = 1.d0/r0_kgpm3
   av_rec = 1.d0/av_m2ps
   kv_rec = 1.d0/kv_m2ps

   ys(0) = 0.d0
   do j = 1, jm
      ys(j) = (dble(j) - 0.5d0)*dy_m
   end do

   zs(0) = 0.d0
   do k = 1, km
      zs(k) = -(dble(km) - dble(k) + 0.5d0)*dz_m
   end do

   lat_max_deg = 1.d0/(ys(jm)*er_rec*180.d0*pi_rec)
   rec_600 = 1.d0/600.d0

   ! コリオリ係数
   do j = 0, jm+1
      fs_psec(j) = 0
      !fs_psec(j) = 2.d0 * om_psec * sin( (dble(j)-0.5d0) * dy_m * er_rec )
   end do

   ! 強制力
   !dpi=4.d0*datan(1.d0)
   !do j=1,jm
      !do i=1,im
         !tx_npm2(i,j) = t0_npm2 * cos( dpi*( (dble(j)-0.5d0) / dble(jm) ) ) !亜熱帯循環
         !ty_npm2(i,j) = 0.d0
         !at_c(i,j)      = -20.d0 *( dble(j)-0.5d0 ) / dble(jm) + 30.d0 ! 東西一様水温
         !sf_psumps(i,j) = 0.d0
      !end do
   !end do

   do j = 1, jm
      ! 赤道から北への緯度 [deg]
      lat_deg = (ys(j) * er_rec) * 180.d0 * pi_rec

      do i = 1, im
         ! 東西風応力：低緯度で西向き，高緯度で東向き
         tx_npm2(i,j)   = -tau_tr * exp( -((lat_deg - lat_tr)*wid_tr_rec)**2 ) &
                          +tau_w  * exp( -((lat_deg - lat_w )*wid_w_rec )**2 )
         ty_npm2(i,j)   = 0.d0
         at_c(i,j)      = 30.d0 - 20.d0 * lat_deg * lat_max_deg  ! 東西一様水温
         sf_psumps(i,j) = 0.d0
      end do
   end do


   ! 初期化
   uc_mps(:,:,:) = 0.d0
   vc_mps(:,:,:) = 0.d0
   ec_m(:,:)     = 0.d0

   ww_mps(:,:,:)  = 0.d0
   pp_npm2(:,:,:) = 0.d0

   adv_x = 0.d0
   adv_y = 0.d0
   adv_z = 0.d0
   rr_f = 0.d0
   cor = 0.d0
   pre = 0.d0
   ax = 0.d0
   ay = 0.d0
   az = 0.d0
   div = 0.d0
   kx = 0.d0
   ky = 0.d0

   ! 初期設定
   if( time_to_start_sec == 0.d0 ) then ! 初期化
      ua_mps(:,:,:) = 0
      ub_mps(:,:,:) = 0
      va_mps(:,:,:) = 0
      vb_mps(:,:,:) = 0
      ea_m(:,:)     = 0
      eb_m(:,:)     = 0
      time_sec      = 0

      do k = 1, km
         ! zs(k) は負。深いほど小さい（例：-575）
         ! ここでは簡単に 0〜600m で線形に
         tb_c(:,:,k) = 5.d0 + (30.d0-5.d0) * ( (zs(k)+600.d0) * rec_600 )
         ta_c(:,:,k) = tb_c(:,:,k)
      end do

      ! 境界用のダミーも埋めておく（安心）
      tb_c(:,:,0)    = tb_c(:,:,1)
      tb_c(:,:,km+1) = tb_c(:,:,km)
      ta_c(:,:,0)    = ta_c(:,:,1)
      ta_c(:,:,km+1) = ta_c(:,:,km)

      sa_psu(:,:,:) = 0.d0 !- 34からの偏差
      sb_psu(:,:,:) = 0.d0
     
   else ! 継続
      open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
      read(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m
      close(10)
   endif

   time_to_output_sec = time_sec + output_interval_sec

   do
      if( time_sec > time_to_end_sec ) exit

      !======================================================================
      ! OpenMP並列領域開始
      !======================================================================
      !$omp parallel default(none) &
      !$omp shared(ua_mps, ub_mps, uc_mps, va_mps, vb_mps, vc_mps, &
      !$omp        wa_mps, wb_mps, wc_mps, ta_c, tb_c, tc_c, &
      !$omp        sa_psu, sb_psu, sc_psu, ww_mps, pp_npm2, rr_kgpm3, &
      !$omp        ea_m, eb_m, ec_m, tx_npm2, ty_npm2, at_c, sf_psumps, &
      !$omp        fs_psec, ys, zs, &
      !$omp        dx_rec, dy_rec, dz_rec, dx_2_rec, dy_2_rec, dz_2_rec, &
      !$omp        r0_rec, av_rec, kv_rec) &
      !$omp private(i, j, k, n, adv_x, adv_y, adv_z, cor, pre, &
      !$omp         ax, ay, az, div, kx, ky, adx, ady, adz, &
      !$omp         dfx, dfy, dfz, frc, rr_f, gt, gs, tav, sav)

      ! 診断量の計算

     ! w
     !$omp do schedule(static)
     do j = 1, jm
        do i = 1, im
           ww_mps(i,j,0) = 0.d0
           do k = 1, km
              ww_mps(i,j,k) = ww_mps(i,j,k-1) - (ub_mps(i,j,k) - ub_mps(i-1,j,k))*dz_m*dx_rec - (vb_mps(i,j,k) - vb_mps(i,j-1,k))*dz_m*dy_rec
           end do
           ww_mps(i,j,km) = 0.d0
           wb_mps(i,j,0:km) = ww_mps(i,j,0:km)
        end do
     end do
     !$omp end do

      ! r
      !$omp do schedule(static) collapse(2)
      do k = km, 1, -1
         !pre_npm2 = pre_npm2 + gr_mps2 * r0_kgpm3 * dz_m !- 圧力はrrで計算すべきだがr0で近似
         do j = 1, jm
            do i = 1, im
               rr_kgpm3(i,j,k) = r0_kgpm3*(1 - alphat_pk*(tb_c(i,j,k) - t0_c) + alphas_ppsu*sb_psu(i,j,k))
               ! rr_kgpm3(i,j,k) = r0_kgpm3 !温度・塩分効果考慮なしver確認
            end do
         end do
      end do
      !$omp end do

      ! p
      !$omp do schedule(static)
      do j = 1, jm
         do i = 1, im
            rr_f = 0.d0
            pp_npm2(i,j,km) = gr_mps2*r0_kgpm3*(dz_m*0.5d0 + eb_m(i,j))
            do k = km-1, 1, -1
               rr_f = rr_f + gr_mps2*(rr_kgpm3(i,j,k) - r0_kgpm3)*dz_m
               pp_npm2(i,j,k) = gr_mps2*r0_kgpm3*(eb_m(i,j) - zs(k)) + rr_f
               ! pp_npm2(i,j,k) = gr_mps2*r0_kgpm3*(eb_m(i,j) - zs(k)) !温度・塩分効果考慮なしver確認
            end do
         end do
      end do
      !$omp end do

      ! 予報量の計算
      ! u
      !$omp do schedule(static) collapse(2)
      do k = 1, km
         do j = 1, jm
            do i = 1, im - 1
               adv_x = -(((ub_mps(i+1,j,k) + ub_mps(i,j,k))**2)*0.25d0 - ((ub_mps(i,j,k) + ub_mps(i-1,j,k))**2)*0.25d0)*dx_rec
               adv_y = -(((ub_mps(i,j+1,k) + ub_mps(i,j,k))*0.5d0*(vb_mps(i+1,j,k) + vb_mps(i,j,k)))*0.5d0 - ((ub_mps(i,j,k) + ub_mps(i,j-1,k))*0.5d0*(vb_mps(i+1,j-1,k) + vb_mps(i,j-1,k)))*0.5d0)*dy_rec
               adv_z = -(((ub_mps(i,j,k+1) + ub_mps(i,j,k))*0.5d0*(wb_mps(i+1,j,k) + wb_mps(i,j,k)))*0.5d0 - ((ub_mps(i,j,k) + ub_mps(i,j,k-1))*0.5d0*(wb_mps(i+1,j,k-1) + wb_mps(i,j,k-1)))*0.5d0)*dz_rec
               cor = (f0_psec + bt_pmsec*ys(j))*(vb_mps(i,j-1,k) + vb_mps(i,j,k) + vb_mps(i+1,j-1,k) + vb_mps(i+1,j,k))*0.25d0
               pre = -(pp_npm2(i+1,j,k) - pp_npm2(i,j,k))*r0_rec*dx_rec
               ax = ah_m2ps*(ua_mps(i+1,j,k) - 2*ua_mps(i,j,k) + ua_mps(i-1,j,k))*dx_2_rec
               ay = ah_m2ps*(ua_mps(i,j+1,k) - 2*ua_mps(i,j,k) + ua_mps(i,j-1,k))*dy_2_rec
               az = av_m2ps*(ua_mps(i,j,k+1) - 2*ua_mps(i,j,k) + ua_mps(i,j,k-1))*dz_2_rec
               uc_mps(i,j,k) = ua_mps(i,j,k) + 2*dt_sec*(adv_x + adv_y + adv_z + cor + pre + ax + ay + az)
            end do
         end do
      end do
      !$omp end do
     
      ! v
      !$omp do schedule(static) collapse(2)
      do k = 1, km
         do j = 1, jm - 1
            do i = 1, im
               adv_x = -(((vb_mps(i+1,j,k) + vb_mps(i,j,k))*0.5d0*(ub_mps(i,j+1,k) + ub_mps(i,j,k)))*0.5d0 - ((vb_mps(i,j,k) + vb_mps(i-1,j,k))*0.5d0*(ub_mps(i-1,j+1,k) + ub_mps(i-1,j,k)))*0.5d0)*dx_rec
               adv_y = -(((vb_mps(i,j+1,k) + vb_mps(i,j,k))**2)*0.25d0 - ((vb_mps(i,j,k) + vb_mps(i,j-1,k))**2)*0.25d0)*dy_rec
               adv_z = -(((vb_mps(i,j,k+1) + vb_mps(i,j,k))*0.5d0*(wb_mps(i,j+1,k) + wb_mps(i,j,k)))*0.5d0 - ((vb_mps(i,j,k) + vb_mps(i,j,k-1))*0.5d0*(wb_mps(i,j+1,k-1) + wb_mps(i,j,k-1)))*0.5d0)*dz_rec
               cor = -(f0_psec + bt_pmsec*ys(j))*(ub_mps(i-1,j+1,k) + ub_mps(i,j+1,k) + ub_mps(i-1,j,k) + ub_mps(i,j,k))*0.25d0
               pre = -(pp_npm2(i,j+1,k) - pp_npm2(i,j,k))*r0_rec*dy_rec
               ax = ah_m2ps*(va_mps(i+1,j,k) - 2*va_mps(i,j,k) + va_mps(i-1,j,k))*dx_2_rec
               ay = ah_m2ps*(va_mps(i,j+1,k) - 2*va_mps(i,j,k) + va_mps(i,j-1,k))*dy_2_rec
               az = av_m2ps*(va_mps(i,j,k+1) - 2*va_mps(i,j,k) + va_mps(i,j,k-1))*dz_2_rec
               vc_mps(i,j,k) = va_mps(i,j,k) + 2*dt_sec*(adv_x + adv_y + adv_z + cor + pre + ax + ay + az)
            end do
         end do
      end do
      !$omp end do

      ! e
      !$omp do schedule(static)
      do j = 1, jm
         do i = 1, im
            div = 0.d0
            do k = 1, km
                div = div - (ub_mps(i,j,k) - ub_mps(i-1,j,k))*dz_m*dx_rec - (vb_mps(i,j,k) - vb_mps(i,j-1,k))*dz_m*dy_rec
            end do
            kx = kh_m2ps*(ea_m(i+1,j)-2*ea_m(i,j)+ea_m(i-1,j))*dx_2_rec
            ky = kh_m2ps*(ea_m(i,j+1)-2*ea_m(i,j)+ea_m(i,j-1))*dy_2_rec
            ec_m(i,j) = ea_m(i,j) + 2.d0*dt_sec*(div + kx + ky)
         end do
      end do
      !$omp end do

      ! t
      !$omp do schedule(static) collapse(2)
      do k = 1, km
         do j = 1, jm
            do i = 1, im
               adx= -((tb_c(i+1,j,k)+tb_c(i,j,k))*ub_mps(i,j,k)*0.5d0 - (tb_c(i,j,k)+tb_c(i-1,j,k))*ub_mps(i-1,j,k)*0.5d0)*dx_rec
               ady= -((tb_c(i,j+1,k)+tb_c(i,j,k))*vb_mps(i,j,k)*0.5d0 - (tb_c(i,j,k)+tb_c(i,j-1,k))*vb_mps(i,j-1,k)*0.5d0)*dy_rec
               adz= -((tb_c(i,j,k+1)+tb_c(i,j,k))*wb_mps(i,j,k)*0.5d0 - (tb_c(i,j,k)+tb_c(i,j,k-1))*wb_mps(i,j,k-1)*0.5d0)*dz_rec
               dfx= kh_m2ps*(tb_c(i+1,j,k) - 2*tb_c(i,j,k) + tb_c(i-1,j,k))*dx_2_rec
               dfy= kh_m2ps*(tb_c(i,j+1,k) - 2*tb_c(i,j,k) + tb_c(i,j-1,k))*dy_2_rec
               dfz= kv_m2ps*(tb_c(i,j,k+1) - 2*tb_c(i,j,k) + tb_c(i,j,k-1))*dz_2_rec

               if ( k == km ) then
                  frc = gm_psec * ( at_c(i,j) - tb_c(i,j,km) )
               else
                  frc = 0.d0
               endif

               gt = adx+ady+adz+dfx+dfy+dfz+frc
               tc_c(i,j,k) = ta_c(i,j,k) + 2.d0 * dt_sec * gt
            end do
         end do
      end do
      !$omp end do

      ! s
      !$omp do schedule(static) collapse(2)
      do k = 1, km
         do j = 1, jm
            do i = 1, im
               adx= -((sb_psu(i+1,j,k)+sb_psu(i,j,k))*ub_mps(i,j,k)*0.5d0 - (sb_psu(i,j,k)+sb_psu(i-1,j,k))*ub_mps(i-1,j,k)*0.5d0)*dx_rec
               ady= -((sb_psu(i,j+1,k)+sb_psu(i,j,k))*vb_mps(i,j,k)*0.5d0 - (sb_psu(i,j,k)+sb_psu(i,j-1,k))*vb_mps(i,j-1,k)*0.5d0)*dy_rec
               adz= -((sb_psu(i,j,k+1)+sb_psu(i,j,k))*wb_mps(i,j,k)*0.5d0 - (sb_psu(i,j,k)+sb_psu(i,j,k-1))*wb_mps(i,j,k-1)*0.5d0)*dz_rec
               dfx= kh_m2ps*(sb_psu(i+1,j,k) - 2*sb_psu(i,j,k) + sb_psu(i-1,j,k))*dx_2_rec
               dfy= kh_m2ps*(sb_psu(i,j+1,k) - 2*sb_psu(i,j,k) + sb_psu(i,j-1,k))*dy_2_rec
               dfz= kv_m2ps*(sb_psu(i,j,k+1) - 2*sb_psu(i,j,k) + sb_psu(i,j,k-1))*dz_2_rec
              
               gs = adx+ady+adz+dfx+dfy+dfz
               sc_psu(i,j,k) = sa_psu(i,j,k) + 2.d0 * dt_sec * gs
            end do
         end do
      end do
      !$omp end do


     ! 対流調節（単純化。2回繰り返す）
     do n = 1, 2
        ! 予報したTSで密度を計算する
        !$omp do schedule(static) collapse(2)
        do k = km, 1, -1
           do j = 1, jm
              do i = 1, im
                 rr_kgpm3(i,j,k) = r0_kgpm3 * (1.d0 - alphat_pk *   ( tc_c(i,j,k)   - t0_c ) &
                                &                   + alphas_ppsu * sc_psu(i,j,k) )
              end do
           end do
        end do
        !$omp end do

        ! 当該格子の密度と上の格子の密度を比べ、上の方が重ければ混ぜる
        !$omp do schedule(static)
        do j = 1, jm
           do i = 1, im
              do k = km-1, 1, -1
                 if( rr_kgpm3(i,j,k+1) <= rr_kgpm3(i,j,k) ) cycle

                 tav = 0.5d0 * ( tc_c(i,j,k) + tc_c(i,j,k+1) )
                 tc_c(i,j,k  ) = tav
                 tc_c(i,j,k+1) = tav

                 sav = 0.5d0 * ( sc_psu(i,j,k) + sc_psu(i,j,k+1) )
                 sc_psu(i,j,k  ) = sav
                 sc_psu(i,j,k+1) = sav
              end do
           end do
        end do
        !$omp end do
     end do


     ! 東西境界条件
     !$omp single
     uc_mps(0   ,:,:) = 0.d0
     uc_mps(im  ,:,:) = 0.d0
     uc_mps(im+1,:,:) = 0.d0  !- dummy

     vc_mps(0   ,:,:) = vc_mps(1   ,:,:)
     vc_mps(im+1,:,:) = vc_mps(im   ,:,:)

     ec_m(0   ,:)     = ec_m(1   ,:) 
     ec_m(im+1,:)     = ec_m(im   ,:)

     tc_c(0   ,:,:) = tc_c(1 ,:,:)
     tc_c(im+1,:,:) = tc_c(im,:,:)

     sc_psu(0   ,:,:) = sc_psu(1 ,:,:)
     sc_psu(im+1,:,:) = sc_psu(im,:,:)

     ! 南北境界条件
     uc_mps(:,0   ,:) = uc_mps(:,1   ,:) 
     uc_mps(:,jm+1,:) = uc_mps(:,jm   ,:)

     vc_mps(:,0   ,:) = 0.d0
     vc_mps(:,jm  ,:) = 0.d0
     vc_mps(:,jm+1,:) = 0.d0

     ec_m(:,0   )     = ec_m(:,1   ) 
     ec_m(:,jm+1)     = ec_m(:,jm   )

     tc_c(:,0   ,:) = tc_c(:,1 ,:)
     tc_c(:,jm+1,:) = tc_c(:,jm,:)

     sc_psu(:,0   ,:) = sc_psu(:,1 ,:)
     sc_psu(:,jm+1,:) = sc_psu(:,jm,:)

     ! 上下境界条件
     uc_mps(:,:,km+1) = uc_mps(:,:,km) + tx_npm2*dz_m*r0_rec*av_rec
     uc_mps(:,:,0   ) = uc_mps(:,:,1)

     vc_mps(:,:,km+1) = vc_mps(:,:,km) + ty_npm2*dz_m*r0_rec*av_rec
     vc_mps(:,:,0   ) = vc_mps(:,:,1)

     tc_c(:,:,km+1) = tc_c(:,:,km)
     tc_c(:,:,0   ) = tc_c(:,:,1 )

     sc_psu(:,:,km+1) = sc_psu(:,:,km) + dz_m * kv_rec * sf_psumps(:,:)
     sc_psu(:,:,0   ) = sc_psu(:,:,1 )
     !$omp end single

     ! n+1ステップが求まったので配列をシフトさせる
     ! ua_mps=ub_mps, ub_mps=uc_mps 
     !$omp workshare
     ua_mps(:,:,:)    = ub_mps(:,:,:) +0.5d0*asf* ( ua_mps(:,:,:) - 2.d0*ub_mps(:,:,:) + uc_mps(:,:,:) )
     va_mps(:,:,:)    = vb_mps(:,:,:) +0.5d0*asf* ( va_mps(:,:,:) - 2.d0*vb_mps(:,:,:) + vc_mps(:,:,:) )
     ea_m(:,:)        = eb_m(:,:)     +0.5d0*asf* ( ea_m(:,:)     - 2.d0*eb_m(:,:)     + ec_m(:,:)     )
     ta_c(:,:,:)      = tb_c(:,:,:)   +0.5d0*asf* ( ta_c(:,:,:)   - 2.d0*tb_c(:,:,:)   + tc_c(:,:,:)   )
     sa_psu(:,:,:)    = sa_psu(:,:,:) +0.5d0*asf* ( sa_psu(:,:,:) - 2.d0*sb_psu(:,:,:) + sc_psu(:,:,:) )
     
     ! n+1 step -> ub_mps,vb_mps,eb_m,tb,sb
     ub_mps(:,:,:)    = uc_mps(:,:,:)
     vb_mps(:,:,:)    = vc_mps(:,:,:)
     eb_m(:,:)        = ec_m(:,:)
     tb_c(:,:,:)      = tc_c(:,:,:)
     sb_psu(:,:,:)    = sc_psu(:,:,:)
     !$omp end workshare

     !======================================================================
     ! OpenMP並列領域終了
     !======================================================================
     !$omp end parallel

     time_sec = time_sec + dt_sec

     if( time_sec < time_to_output_sec ) cycle

     n = idnint( time_sec / output_interval_sec )

     !======================================================================
     ! OpenMP並列領域開始（診断量再計算用）
     !======================================================================
     !$omp parallel default(none) &
     !$omp shared(ub_mps, vb_mps, wb_mps, ww_mps, &
     !$omp        tb_c, sb_psu, rr_kgpm3, pp_npm2, eb_m, &
     !$omp        dx_rec, dy_rec, dz_rec, zs) &
     !$omp private(i, j, k, rr_f)

     ! 診断量の再計算（記録用）
     ! w
     !$omp do schedule(static)
     do j = 1, jm
        do i = 1, im
           ww_mps(i,j,0) = 0.d0
           do k = 1, km
              ww_mps(i,j,k) = ww_mps(i,j,k-1) - (ub_mps(i,j,k) - ub_mps(i-1,j,k))*dz_m*dx_rec - (vb_mps(i,j,k) - vb_mps(i,j-1,k))*dz_m*dy_rec
           end do
           ww_mps(i,j,km) = 0.d0
           wb_mps(i,j,0:km) = ww_mps(i,j,0:km)
        end do
     end do
     !$omp end do

     ! r
     !$omp do schedule(static) collapse(2)
     do k = km, 1, -1
        do j = 1, jm
           do i = 1, im
              rr_kgpm3(i,j,k) = r0_kgpm3*(1 - alphat_pk*(tb_c(i,j,k) - t0_c) + alphas_ppsu*sb_psu(i,j,k))
              ! rr_kgpm3(i,j,k) = r0_kgpm3 !温度・塩分効果考慮なしver確認
           end do
        end do
     end do
     !$omp end do

     ! p
     !$omp do schedule(static)
     do j = 1, jm
        do i = 1, im
           rr_f = 0.d0
           pp_npm2(i,j,km) = gr_mps2*r0_kgpm3*(dz_m*0.5d0 + eb_m(i,j))
           do k = km-1, 1, -1
              rr_f = rr_f + gr_mps2*(rr_kgpm3(i,j,k) - r0_kgpm3)*dz_m
              pp_npm2(i,j,k) = gr_mps2*r0_kgpm3*(eb_m(i,j) - zs(k)) + rr_f
              ! pp_npm2(i,j,k) = gr_mps2*r0_kgpm3*(eb_m(i,j) - zs(k)) !温度・塩分効果考慮なしver確認
           end do
        end do
     end do
     !$omp end do

     !======================================================================
     ! OpenMP並列領域終了
     !======================================================================
     !$omp end parallel

     ! 記録用（解析・作図用）
     write(buff,'(a,i6.6)') trim(expid)//'.n',n
     open(10,file=trim(buff),form='unformatted',access='stream')
     write(10) time_sec,ub_mps,vb_mps,eb_m,ww_mps,pp_npm2,tb_c,sb_psu,rr_kgpm3
     close(10)


     ! 継続用
     open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
     write(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m,ta_c,tb_c,sa_psu,sb_psu
     close(10)

     ! 画面で監視
     write(6,*) 'time [s] = ',time_sec
     write(6,*) 'umin: ',minval(ub_mps(1:im,1:jm,1:km)),'(',minloc(ub_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'umax: ',maxval(ub_mps(1:im,1:jm,1:km)),'(',maxloc(ub_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'vmin: ',minval(vb_mps(1:im,1:jm,1:km)),'(',minloc(vb_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'vmax: ',maxval(vb_mps(1:im,1:jm,1:km)),'(',maxloc(vb_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'wmin: ',minval(ww_mps(1:im,1:jm,1:km)),'(',minloc(ww_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'wmax: ',maxval(ww_mps(1:im,1:jm,1:km)),'(',maxloc(ww_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'emin: ',minval(eb_m(1:im,1:jm)),'(',minloc(eb_m(1:im,1:jm)),')'
     write(6,*) 'emax: ',maxval(eb_m(1:im,1:jm)),'(',maxloc(eb_m(1:im,1:jm)),')'

     time_to_output_sec = time_to_output_sec + output_interval_sec

   end do

   ! プログラム終了時の時刻を表示
   call date_and_time(date, time, zone, values)
   print *, '========================================='
   print *, 'プログラム終了'
   print *, '日付: ', date(1:4), '/', date(5:6), '/', date(7:8)
   print *, '時刻: ', time(1:2), ':', time(3:4), ':', time(5:10)
   print *, '========================================='
    

end program main