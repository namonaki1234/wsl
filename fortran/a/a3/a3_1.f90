program a3
    implicit none

    ! ==========================
    ! パラメータ設定
    ! 格子条件
    double precision, parameter :: H = 10.0d0
    integer, parameter :: IM = 30, JM = 30
    double precision :: dx, r_y, a ! dx は格子幅, r_y は公比, a はy格子生成用の初項

    ! 流れ場の初期条件 (スカラー)
    double precision, parameter :: M0 = 2.9d0, T0 = 293.0d0, R0 = 287.1d0
    double precision, parameter :: gam = 1.4d0, rho0 = 1.2d0, dt = 1.0d-6

    ! 変数宣言
    integer :: i, j, k_dim, n_step ! k_dim は QEF の4成分用
    double precision, dimension(0:IM, 0:JM) :: x, y
    double precision, dimension(0:IM, 0:JM) :: rho, u, v, p, energy, T_gas, E_bar

    ! 保存量ベクトルと流束ベクトル
    double precision, dimension(0:IM, 0:JM, 4) :: Q, E, F

    ! メトリクスとヤコビアン、反変速度
    double precision, dimension(0:IM, 0:JM) :: xe, ye, xh, yh, Jn
    double precision, dimension(0:IM, 0:JM) :: xix, xiy, ex, ey
    double precision, dimension(0:IM, 0:JM) :: uu, vv

    double precision :: u_init, p_init, energy_init, c_sound_init

    ! ==========================
    ! パラメータ初期化
    dx = 3.0d0*H/dble(IM)
    r_y = 1.1d0

    ! === 格子生成 ===
    if (abs(r_y - 1.0d0) < 1.0d-12) then
        print *, "Error: r_y (公比) が1.0です。この格子生成式では使用できません。"
        stop
    end if
    do i = 0, IM
        do j = 0, JM
            x(i,j) = dx * dble(i)
        end do
    end do
    do j = 0, JM
        a = H * (r_y - 1.0d0) / (r_y**dble(JM) - 1.0d0)
        do i = 0, IM
            if (j == 0) then
                y(i,j) = 0.0d0
            else
                y(i,j) = a * (r_y**dble(j) - 1.0d0) / (r_y - 1.0d0)
            end if
        end do
    end do

    ! ==========================
    ! 初期条件設定
    c_sound_init = dsqrt(gam*R0*T0)
    u_init = M0*c_sound_init
    p_init = rho0*R0*T0
    energy_init = p_init/(gam-1.0d0) + 0.5d0*rho0*u_init**2

    do j = 0, JM
        do i = 0, IM
            rho(i,j)    = rho0
            u(i,j)      = u_init
            v(i,j)      = 0.0d0
            p(i,j)      = p_init
            energy(i,j) = energy_init
            T_gas(i,j)  = T0
        end do
    end do
    ! 流入条件 (i=0) は初期値で固定 (上記で設定済み)

    ! ==========================
    ! メトリクスとヤコビアンの計算 (先輩のコードの定義に従う)
    ! xe, ye, xh, yh (∂x/∂ξ, ∂y/∂ξ, ∂x/∂η, ∂y/∂η に対応する差分)
    do i=0,IM
        do j=1,JM-1
            xh(i,j)=(x(i,j+1)-x(i,j-1))/2.0d0
            yh(i,j)=(y(i,j+1)-y(i,j-1))/2.0d0
        end do
    end do
    do j=0,JM
        do i=1,IM-1
            xe(i,j)=(x(i+1,j)-x(i-1,j))/2.0d0
            ye(i,j)=(y(i+1,j)-y(i-1,j))/2.0d0
        end do
    end do
    ! 境界
    do j=0,JM
        xe(0,j)=(-3.0d0*x(0,j)+4.0d0*x(1,j)-x(2,j))/2.0d0
        xe(IM,j)=(3.0d0*x(IM,j)-4.0d0*x(IM-1,j)+x(IM-2,j))/2.0d0
        ye(0,j)=(-3.0d0*y(0,j)+4.0d0*y(1,j)-y(2,j))/2.0d0
        ye(IM,j)=(3.0d0*y(IM,j)-4.0d0*y(IM-1,j)+y(IM-2,j))/2.0d0
    end do
    do i=0,IM
        xh(i,0)=(-3.0d0*x(i,0)+4.0d0*x(i,1)-x(i,2))/2.0d0
        xh(i,JM)=(3.0d0*x(i,JM)-4.0d0*x(i,JM-1)+x(i,JM-2))/2.0d0
        yh(i,0)=(-3.0d0*y(i,0)+4.0d0*y(i,1)-y(i,2))/2.0d0
        yh(i,JM)=-(-3.0d0*y(i,JM)+4.0d0*y(i,JM-1)-y(i,JM-2))/2.0d0
    end do

    ! Jn (ヤコビアンの逆数 1/J_geom), xix, xiy, ex, ey (メトリック係数)
    do i=0,IM
        do j=0,JM
            Jn(i,j)=1.0d0/(xe(i,j)*yh(i,j)-ye(i,j)*xh(i,j))
            if (abs(xe(i,j)*yh(i,j)-ye(i,j)*xh(i,j)) < 1.0d-12) then
                 Jn(i,j) = 1.0d0 / 1.0d-12
            end if
            xix(i,j)=Jn(i,j)*yh(i,j)
            xiy(i,j)=Jn(i,j)*xh(i,j)
            ex(i,j)=-Jn(i,j)*ye(i,j)
            ey(i,j)=Jn(i,j)*xe(i,j)
        end do
    end do

    ! ==========================
    ! メインループ
    do n_step = 1, 10000
        call calc_QEF(IM, JM, rho, u, v, energy, p, &
                                  Jn, xix, xiy, ex, ey, uu, vv, &
                                  Q, E, F)

        ! 時間発展 (内部点)
        do k_dim = 1, 4
            do j = 1, JM - 1
                do i = 1, IM - 1
                    Q(i,j,k_dim) = Q(i,j,k_dim) - dt * ( &
                        (E(i+1,j,k_dim) - E(i-1,j,k_dim))/(2.0d0) + &
                        (F(i,j+1,k_dim) - F(i,j-1,k_dim))/(2.0d0) )
                end do
            end do
        end do

        ! 物理量へ変換 (内部点)
        do j = 1, JM-1
            do i = 1, IM-1
                rho(i,j) = Q(i,j,1)*Jn(i,j)
                if (rho(i,j) < 1.0d-6) rho(i,j) = 1.0d-6
                u(i,j)   = (Q(i,j,2)*Jn(i,j))/rho(i,j)
                v(i,j)   = (Q(i,j,3)*Jn(i,j))/rho(i,j)
                energy(i,j) = Q(i,j,4)*Jn(i,j)
            end do
        end do

        ! 境界条件の適用 (Q に対して外挿)
        ! i=0 (流入) は物理量が固定なので、Q(0,j,k) も初期状態から計算した固定値とする
        ! (ここでは簡単のため、メインループの最初の calc_QEF で計算された値を維持すると解釈)

        ! 上壁 (j=JM)
        do k_dim = 1, 4
            do i = 1, IM - 1
                Q(i,JM,k_dim) = 2.0d0*Q(i,JM-1,k_dim) - Q(i,JM-2,k_dim)
            end do
        end do
        ! 下壁 (j=0)
        do k_dim = 1, 4
            do i = 1, IM - 1
                Q(i,0,k_dim) = 2.0d0*Q(i,1,k_dim) - Q(i,2,k_dim)
            end do
        end do
        ! 流出境界 (i=IM)
        do k_dim = 1, 4
            do j = 0, JM
                Q(IM,j,k_dim) = 2.0d0*Q(IM-1,j,k_dim) - Q(IM-2,j,k_dim)
            end do
        end do

        ! 境界の物理量を更新 (Q から rho, u, v, energy へ)
        ! 上壁 (j=JM)
        do i = 1, IM -1
            rho(i,JM) = Q(i,JM,1)*Jn(i,JM)
            if (rho(i,JM) < 1.0d-6) rho(i,JM) = 1.0d-6
            u(i,JM)   = (Q(i,JM,2)*Jn(i,JM))/rho(i,JM)
            v(i,JM)   = 0.0d0 ! 壁条件
            ! energyの更新: u,vが求まったので、全エネルギーから運動エネルギーを引いたものが内部エネルギー。
            ! 外挿されたQ4から求めたenergyが妥当か、あるいは圧力を外挿してエネルギーを再計算するか検討の余地あり。
            ! ここではQ4の外挿値を使い、v=0を反映させて再計算。
            energy(i,JM) = (Q(i,JM,4)*Jn(i,JM)) ! まずQ4から仮のe_total
            E_bar(i,JM) = energy(i,JM)/rho(i,JM) - 0.5d0 * (u(i,JM)**2 + v(i,JM)**2) ! 比内部エネルギー
            if (E_bar(i,JM) < 1.0d-6 * (R0 / (gam - 1.0d0)) ) then
                 E_bar(i,JM) = 1.0d-6 * (R0 / (gam - 1.0d0))
            end if
            p(i,JM) = (gam-1.0d0)*rho(i,JM)*E_bar(i,JM)
            if (p(i,JM) < 1.0d0) p(i,JM) = 1.0d0
            energy(i,JM) = rho(i,JM) * (E_bar(i,JM) + 0.5d0 * (u(i,JM)**2 + v(i,JM)**2) )
        end do
        ! 下壁 (j=0)
        do i = 1, IM -1
            rho(i,0) = Q(i,0,1)*Jn(i,0)
            if (rho(i,0) < 1.0d-6) rho(i,0) = 1.0d-6
            u(i,0)   = (Q(i,0,2)*Jn(i,0))/rho(i,0)
            v(i,0)   = 0.0d0 ! 壁条件
            energy(i,0) = (Q(i,0,4)*Jn(i,0))
            E_bar(i,0) = energy(i,0)/rho(i,0) - 0.5d0 * (u(i,0)**2 + v(i,0)**2)
            if (E_bar(i,0) < 1.0d-6 * (R0 / (gam - 1.0d0)) ) then
                 E_bar(i,0) = 1.0d-6 * (R0 / (gam - 1.0d0))
            end if
            p(i,0) = (gam-1.0d0)*rho(i,0)*E_bar(i,0)
            if (p(i,0) < 1.0d0) p(i,0) = 1.0d0
            energy(i,0) = rho(i,0) * (E_bar(i,0) + 0.5d0 * (u(i,0)**2 + v(i,0)**2) )
        end do
        ! 流出境界 (i=IM)
        do j = 0, JM
            rho(IM,j) = Q(IM,j,1)*Jn(IM,j)
            if (rho(IM,j) < 1.0d-6) rho(IM,j) = 1.0d-6
            u(IM,j)   = (Q(IM,j,2)*Jn(IM,j))/rho(IM,j)
            v(IM,j)   = (Q(IM,j,3)*Jn(IM,j))/rho(IM,j)
            energy(IM,j) = Q(IM,j,4)*Jn(IM,j)
        end do

        ! 熱力学量の再計算 (全領域、ただし流入i=0は初期値固定を維持)
        do j = 0, JM
            do i = 0, IM
                if (i == 0) then ! 流入境界は初期値固定
                    p(i,j) = p_init
                    T_gas(i,j) = T0
                    ! E_bar(i,j) は使用しないので計算不要でも可
                else
                    E_bar(i,j) = energy(i,j)/rho(i,j) - (u(i,j)**2 + v(i,j)**2)/2.0d0
                    if (E_bar(i,j) < 1.0d-6 * (R0 / (gam - 1.0d0)) ) then
                        E_bar(i,j) = 1.0d-6 * (R0 / (gam - 1.0d0))
                    end if
                    p(i,j) = (gam-1.0d0)*rho(i,j)*E_bar(i,j)
                    if (p(i,j) < 1.0d0) p(i,j) = 1.0d0
                    T_gas(i,j) = p(i,j)/(rho(i,j)*R0)
                    if (T_gas(i,j) < 1.0d0) T_gas(i,j) = 1.0d0
                end if
            end do
        end do

        if (mod(n_step, 100) == 0) then
            print*, 'n_step = ', n_step, ', time = ', dble(n_step)*dt
        end if
    end do ! メインループ終わり

! === データ出力（MicroAVS用 DAT）===
    open (10, file='a2.dat', status='replace')
    do j = 0, JM
      do i = 0, IM
        write (10, *) x(i, j), y(i, j), u(i, j), v(i, j)
      end do
    end do
    close (10)

    ! === MicroAVSのFLDヘッダ出力 ===
    open (11, file='a2.fld', status='replace')
    write (11, '(A)') '# AVS field file'
    write (11, '(A)') 'ndim = 2'
    write (11, '(A,I5)') 'dim1 =', IM + 1
    write (11, '(A,I5)') 'dim2 =', JM + 1
    write (11, '(A)') 'nspace = 2'
    write (11, '(A)') 'veclen = 2'
    write (11, '(A)') 'data = double'
    write (11, '(A)') 'field = irregular'
    write (11, '(A)') 'label = u v'
    write (11, '(A)') 'variable 1 file=a2.dat filetype=ascii skip=0 offset=2 stride=4'
    write (11, '(A)') 'variable 2 file=a2.dat filetype=ascii skip=0 offset=3 stride=4'
    write (11, '(A)') 'coord 1 file=a2.dat filetype=ascii skip=0 offset=0 stride=4'
    write (11, '(A)') 'coord 2 file=a2.dat filetype=ascii skip=0 offset=1 stride=4'
    close (11)

    print *, "→ MicroAVS用の .dat および .fld を出力しました。"

contains

    subroutine calc_QEF(IM_s, JM_s, rho_in, u_in, v_in, energy_in, p_in, &
                                    Jn_in, xix_in, xiy_in, ex_in, ey_in, &
                                    uu_out, vv_out, Q_out, E_out, F_out)
        implicit none
        integer, intent(in) :: IM_s, JM_s
        double precision, intent(in) :: rho_in(0:IM_s,0:JM_s), u_in(0:IM_s,0:JM_s), v_in(0:IM_s,0:JM_s)
        double precision, intent(in) :: energy_in(0:IM_s,0:JM_s), p_in(0:IM_s,0:JM_s)
        double precision, intent(in) :: Jn_in(0:IM_s,0:JM_s), xix_in(0:IM_s,0:JM_s), xiy_in(0:IM_s,0:JM_s)
        double precision, intent(in) :: ex_in(0:IM_s,0:JM_s), ey_in(0:IM_s,0:JM_s)

        double precision, intent(out) :: uu_out(0:IM_s,0:JM_s), vv_out(0:IM_s,0:JM_s)
        double precision, intent(out) :: Q_out(0:IM_s,0:JM_s,4)
        double precision, intent(out) :: E_out(0:IM_s,0:JM_s,4)
        double precision, intent(out) :: F_out(0:IM_s,0:JM_s,4)
        integer :: i_s, j_s

        do j_s = 0, JM_s
            do i_s = 0, IM_s
                ! 反変速度
                uu_out(i_s,j_s) = xix_in(i_s,j_s)*u_in(i_s,j_s) + xiy_in(i_s,j_s)*v_in(i_s,j_s)
                vv_out(i_s,j_s) = ex_in(i_s,j_s)*u_in(i_s,j_s)  + ey_in(i_s,j_s)*v_in(i_s,j_s)

                ! 保存量ベクトル Q = Q_tilde * J_geom (Q_tilde = Q / J_geom)
                ! Jn_in は 1/J_geom なので、rho_in/Jn_in は rho_in*J_geom
                Q_out(i_s,j_s,1) = rho_in(i_s,j_s) / Jn_in(i_s,j_s)
                Q_out(i_s,j_s,2) = (rho_in(i_s,j_s)*u_in(i_s,j_s)) / Jn_in(i_s,j_s)
                Q_out(i_s,j_s,3) = (rho_in(i_s,j_s)*v_in(i_s,j_s)) / Jn_in(i_s,j_s)
                Q_out(i_s,j_s,4) = energy_in(i_s,j_s) / Jn_in(i_s,j_s)

                ! ξ方向流束ベクトル E = E_tilde * J_geom
                E_out(i_s,j_s,1) = (rho_in(i_s,j_s)*uu_out(i_s,j_s)) / Jn_in(i_s,j_s)
                E_out(i_s,j_s,2) = (rho_in(i_s,j_s)*u_in(i_s,j_s)*uu_out(i_s,j_s) + xix_in(i_s,j_s)*p_in(i_s,j_s)) / Jn_in(i_s,j_s)
                E_out(i_s,j_s,3) = (rho_in(i_s,j_s)*v_in(i_s,j_s)*uu_out(i_s,j_s) + xiy_in(i_s,j_s)*p_in(i_s,j_s)) / Jn_in(i_s,j_s)
                E_out(i_s,j_s,4) = (energy_in(i_s,j_s) + p_in(i_s,j_s)) * uu_out(i_s,j_s) / Jn_in(i_s,j_s)

                ! η方向流束ベクトル F = F_tilde * J_geom
                F_out(i_s,j_s,1) = (rho_in(i_s,j_s)*vv_out(i_s,j_s)) / Jn_in(i_s,j_s)
                F_out(i_s,j_s,2) = (rho_in(i_s,j_s)*u_in(i_s,j_s)*vv_out(i_s,j_s) + ex_in(i_s,j_s)*p_in(i_s,j_s)) / Jn_in(i_s,j_s)
                F_out(i_s,j_s,3) = (rho_in(i_s,j_s)*v_in(i_s,j_s)*vv_out(i_s,j_s) + ey_in(i_s,j_s)*p_in(i_s,j_s)) / Jn_in(i_s,j_s)
                F_out(i_s,j_s,4) = (energy_in(i_s,j_s) + p_in(i_s,j_s)) * vv_out(i_s,j_s) / Jn_in(i_s,j_s)
            end do
        end do
    end subroutine calc_QEF

end program a3