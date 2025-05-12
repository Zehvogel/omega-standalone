! File generated automatically by O'Mega 3.1.2 release Mar 21 2023
!
!   omega_SM.opt -scatter "e- e+ -> e- nuebar u dbar" -target:parameter_module parameters_sm
!
! with all scattering amplitudes for the process(es)
!
!   flavor combinations:
!
!       1: e- e+ -> e- nuebar u dbar
!
!   color flows:
!
!       1: (  0,  0) (  0,  0) -> (  0,  0) (  0,  0) (  1,  0) (  0, -1)
!
!     NB: i.g. not all color flows contribute to all flavor
!     combinations.  Consult the array FLV_COL_IS_ALLOWED
!     below for the allowed combinations.
!
!   Color Factors:
!
!     (  1,  1): + N
!
!   vanishing or redundant flavor combinations:
!
!
!
module omega_amplitude
  use kinds
  use omega95
  use omega_color, OCF => omega_color_factor
  use parameters_sm
  implicit none
  private
  public :: number_particles_in, number_particles_out, number_color_indices, &
    reset_helicity_selection, new_event, is_allowed, get_amplitude, &
    color_sum, external_masses, openmp_supported, number_spin_states, &
    spin_states, number_flavor_states, flavor_states, number_color_flows, &
    color_flows, number_color_factors, color_factors

  ! DON'T EVEN THINK of removing the following!
  ! If the compiler complains about undeclared
  ! or undefined variables, you are compiling
  ! against an incompatible omega95 module!
  integer, dimension(7), parameter, private :: require = &
    (/ omega_spinors_2010_01_A, omega_spinor_cpls_2010_01_A, &
       omega_vectors_2010_01_A, omega_polarizations_2010_01_A, &
       omega_couplings_2010_01_A, omega_color_2010_01_A, &
       omega_utils_2010_01_A /)

  integer, parameter :: n_prt = 6
  integer, parameter :: n_in = 2
  integer, parameter :: n_out = 4
  integer, parameter :: n_cflow = 1
  integer, parameter :: n_cindex = 2
  integer, parameter :: n_flv = 1
  integer, parameter :: n_hel = 64

  ! NB: you MUST NOT change the value of N_ here!!!
  !     It is defined here for convenience only and must be
  !     compatible with hardcoded values in the amplitude!
  real(kind=default), parameter :: N_ = 3
  logical, parameter :: F = .false.
  logical, parameter :: T = .true.

  integer, dimension(n_prt,n_hel), save, protected :: table_spin_states
  data table_spin_states(:,   1) / -1, -1, -1, -1, -1, -1 /
  data table_spin_states(:,   2) / -1, -1, -1, -1, -1,  1 /
  data table_spin_states(:,   3) / -1, -1, -1, -1,  1, -1 /
  data table_spin_states(:,   4) / -1, -1, -1, -1,  1,  1 /
  data table_spin_states(:,   5) / -1, -1, -1,  1, -1, -1 /
  data table_spin_states(:,   6) / -1, -1, -1,  1, -1,  1 /
  data table_spin_states(:,   7) / -1, -1, -1,  1,  1, -1 /
  data table_spin_states(:,   8) / -1, -1, -1,  1,  1,  1 /
  data table_spin_states(:,   9) / -1, -1,  1, -1, -1, -1 /
  data table_spin_states(:,  10) / -1, -1,  1, -1, -1,  1 /
  data table_spin_states(:,  11) / -1, -1,  1, -1,  1, -1 /
  data table_spin_states(:,  12) / -1, -1,  1, -1,  1,  1 /
  data table_spin_states(:,  13) / -1, -1,  1,  1, -1, -1 /
  data table_spin_states(:,  14) / -1, -1,  1,  1, -1,  1 /
  data table_spin_states(:,  15) / -1, -1,  1,  1,  1, -1 /
  data table_spin_states(:,  16) / -1, -1,  1,  1,  1,  1 /
  data table_spin_states(:,  17) / -1,  1, -1, -1, -1, -1 /
  data table_spin_states(:,  18) / -1,  1, -1, -1, -1,  1 /
  data table_spin_states(:,  19) / -1,  1, -1, -1,  1, -1 /
  data table_spin_states(:,  20) / -1,  1, -1, -1,  1,  1 /
  data table_spin_states(:,  21) / -1,  1, -1,  1, -1, -1 /
  data table_spin_states(:,  22) / -1,  1, -1,  1, -1,  1 /
  data table_spin_states(:,  23) / -1,  1, -1,  1,  1, -1 /
  data table_spin_states(:,  24) / -1,  1, -1,  1,  1,  1 /
  data table_spin_states(:,  25) / -1,  1,  1, -1, -1, -1 /
  data table_spin_states(:,  26) / -1,  1,  1, -1, -1,  1 /
  data table_spin_states(:,  27) / -1,  1,  1, -1,  1, -1 /
  data table_spin_states(:,  28) / -1,  1,  1, -1,  1,  1 /
  data table_spin_states(:,  29) / -1,  1,  1,  1, -1, -1 /
  data table_spin_states(:,  30) / -1,  1,  1,  1, -1,  1 /
  data table_spin_states(:,  31) / -1,  1,  1,  1,  1, -1 /
  data table_spin_states(:,  32) / -1,  1,  1,  1,  1,  1 /
  data table_spin_states(:,  33) /  1, -1, -1, -1, -1, -1 /
  data table_spin_states(:,  34) /  1, -1, -1, -1, -1,  1 /
  data table_spin_states(:,  35) /  1, -1, -1, -1,  1, -1 /
  data table_spin_states(:,  36) /  1, -1, -1, -1,  1,  1 /
  data table_spin_states(:,  37) /  1, -1, -1,  1, -1, -1 /
  data table_spin_states(:,  38) /  1, -1, -1,  1, -1,  1 /
  data table_spin_states(:,  39) /  1, -1, -1,  1,  1, -1 /
  data table_spin_states(:,  40) /  1, -1, -1,  1,  1,  1 /
  data table_spin_states(:,  41) /  1, -1,  1, -1, -1, -1 /
  data table_spin_states(:,  42) /  1, -1,  1, -1, -1,  1 /
  data table_spin_states(:,  43) /  1, -1,  1, -1,  1, -1 /
  data table_spin_states(:,  44) /  1, -1,  1, -1,  1,  1 /
  data table_spin_states(:,  45) /  1, -1,  1,  1, -1, -1 /
  data table_spin_states(:,  46) /  1, -1,  1,  1, -1,  1 /
  data table_spin_states(:,  47) /  1, -1,  1,  1,  1, -1 /
  data table_spin_states(:,  48) /  1, -1,  1,  1,  1,  1 /
  data table_spin_states(:,  49) /  1,  1, -1, -1, -1, -1 /
  data table_spin_states(:,  50) /  1,  1, -1, -1, -1,  1 /
  data table_spin_states(:,  51) /  1,  1, -1, -1,  1, -1 /
  data table_spin_states(:,  52) /  1,  1, -1, -1,  1,  1 /
  data table_spin_states(:,  53) /  1,  1, -1,  1, -1, -1 /
  data table_spin_states(:,  54) /  1,  1, -1,  1, -1,  1 /
  data table_spin_states(:,  55) /  1,  1, -1,  1,  1, -1 /
  data table_spin_states(:,  56) /  1,  1, -1,  1,  1,  1 /
  data table_spin_states(:,  57) /  1,  1,  1, -1, -1, -1 /
  data table_spin_states(:,  58) /  1,  1,  1, -1, -1,  1 /
  data table_spin_states(:,  59) /  1,  1,  1, -1,  1, -1 /
  data table_spin_states(:,  60) /  1,  1,  1, -1,  1,  1 /
  data table_spin_states(:,  61) /  1,  1,  1,  1, -1, -1 /
  data table_spin_states(:,  62) /  1,  1,  1,  1, -1,  1 /
  data table_spin_states(:,  63) /  1,  1,  1,  1,  1, -1 /
  data table_spin_states(:,  64) /  1,  1,  1,  1,  1,  1 /

  integer, dimension(n_prt,n_flv), save, protected :: table_flavor_states
  data table_flavor_states(:,   1) /  11, -11,  11, -12,   2,  -1 / ! e- e+ e- nuebar u dbar

  integer, dimension(n_cindex,n_prt,n_cflow), save, protected :: table_color_flows
  data table_color_flows(:,:,   1) / 0,0,  0,0,  0,0,  0,0,  1,0,  0,-1 /

  logical, dimension(n_prt,n_cflow), save, protected :: table_ghost_flags
  data table_ghost_flags(:,   1) / F,  F,  F,  F,  F,  F /

  integer, parameter :: n_cfactors = 1
  type(OCF), dimension(n_cfactors), save, protected :: table_color_factors
  real(kind=default), parameter, private :: color_factor_000001 = +N_
  data table_color_factors(     1) / OCF(1,1,color_factor_000001) /

  logical, dimension(n_flv, n_cflow), save, protected ::  flv_col_is_allowed
  data flv_col_is_allowed(:,   1) / T /

  complex(kind=default), dimension(n_flv, n_cflow, n_hel), save :: amp

  logical, dimension(n_hel), save :: hel_is_allowed = T
  real(kind=default), dimension(n_hel), save :: hel_max_abs = 0
  real(kind=default), save :: hel_sum_abs = 0, hel_threshold = 1E10_default
  integer, save :: hel_count = 0, hel_cutoff = 100
  integer :: i
  integer, save, dimension(n_hel) :: hel_map = (/(i, i = 1, n_hel)/)
  integer, save :: hel_finite = n_hel

    type(momentum) :: p1, p2, p3, p4, p5, p6
    type(momentum) :: p12, p123, p124, p125, p126, p13, p134, p135, p136, &
      p24, p34, p56
    type(spinor) :: owf_d1_1__6, owf_n1_4, owf_l1_1
    type(conjspinor) :: owf_u1b__1_5, owf_l1b_3, owf_l1b_2
    type(spinor) :: owf_d1_1__136, owf_d1_1__126, owf_n1_134, owf_n1_124
    type(conjspinor) :: owf_u1b__1_135, owf_u1b__1_125, owf_l1b_123
    type(vector) :: owf_a_13, owf_a_12, owf_z_13, owf_z_12, owf_wm_56, &
      owf_wp_34, owf_wp_24
    complex(kind=default) :: oks_l1l1bl1n1bu1_1_d1b__1

contains

  pure function number_particles_in () result (n)
    integer :: n
    n = n_in
  end function number_particles_in

  pure function number_particles_out () result (n)
    integer :: n
    n = n_out
  end function number_particles_out

  pure function number_spin_states () result (n)
    integer :: n
    n = size (table_spin_states, dim=2)
  end function number_spin_states

  pure subroutine spin_states (a)
    integer, dimension(:,:), intent(out) :: a
    a = table_spin_states
  end subroutine spin_states

  pure function number_flavor_states () result (n)
    integer :: n
    n = size (table_flavor_states, dim=2)
  end function number_flavor_states

  pure subroutine flavor_states (a)
    integer, dimension(:,:), intent(out) :: a
    a = table_flavor_states
  end subroutine flavor_states

  pure subroutine external_masses (m, flv)
    real(kind=default), dimension(:), intent(out) :: m
    integer, intent(in) :: flv
    select case (flv)
    case (  1)
      m( 1) = mass(11)
      m( 2) = mass(11)
      m( 3) = mass(11)
      m( 4) = mass(12)
      m( 5) = mass(2)
      m( 6) = mass(1)
    end select
  end subroutine external_masses

  pure function openmp_supported () result (status)
    logical :: status
    status = .false.
  end function openmp_supported

  pure function number_color_indices () result (n)
    integer :: n
    n = size (table_color_flows, dim=1)
  end function number_color_indices

  pure function number_color_flows () result (n)
    integer :: n
    n = size (table_color_flows, dim=3)
  end function number_color_flows

  pure subroutine color_flows (a, g)
    integer, dimension(:,:,:), intent(out) :: a
    logical, dimension(:,:), intent(out) :: g
    a = table_color_flows
    g = table_ghost_flags
  end subroutine color_flows

  pure function number_color_factors () result (n)
    integer :: n
    n = size (table_color_factors)
  end function number_color_factors

  pure subroutine color_factors (cf)
    type(OCF), dimension(:), intent(out) :: cf
    cf = table_color_factors
  end subroutine color_factors

  function color_sum (flv, hel) result (amp2)
    integer, intent(in) :: flv, hel
    real(kind=default) :: amp2
    amp2 = real (omega_color_sum (flv, hel, amp, table_color_factors))
  end function color_sum

  subroutine new_event (p)
    real(kind=default), dimension(0:3,*), intent(in) :: p
    logical :: mask_dirty
    integer :: hel
    call calculate_amplitudes (amp, p, hel_is_allowed)
    if ((hel_threshold .gt. 0) .and. (hel_count .le. hel_cutoff)) then
      call omega_update_helicity_selection (hel_count, amp, hel_max_abs, &
              hel_sum_abs, hel_is_allowed, hel_threshold, hel_cutoff, &
              mask_dirty)
      if (mask_dirty) then
        hel_finite = 0
        do hel = 1, n_hel
          if (hel_is_allowed(hel)) then
            hel_finite = hel_finite + 1
            hel_map(hel_finite) = hel
          end if
        end do
      end if
    end if
  end subroutine new_event

  subroutine reset_helicity_selection (threshold, cutoff)
    real(kind=default), intent(in) :: threshold
    integer, intent(in) :: cutoff
    integer :: i
    hel_is_allowed = T
    hel_max_abs = 0
    hel_sum_abs = 0
    hel_count = 0
    hel_threshold = threshold
    hel_cutoff = cutoff
    hel_map = (/(i, i = 1, n_hel)/)
    hel_finite = n_hel
  end subroutine reset_helicity_selection

  pure function is_allowed (flv, hel, col) result (yorn)
    logical :: yorn
    integer, intent(in) :: flv, hel, col
    yorn = hel_is_allowed(hel) .and. flv_col_is_allowed(flv,col)
  end function is_allowed

  pure function get_amplitude (flv, hel, col) result (amp_result)
    complex(kind=default) :: amp_result
    integer, intent(in) :: flv, hel, col
    amp_result = amp(flv, col, hel)
  end function get_amplitude



  subroutine calculate_amplitudes (amp, k, mask)
    complex(kind=default), dimension(:,:,:), intent(out) :: amp
    real(kind=default), dimension(0:3,*), intent(in) :: k
    logical, dimension(:), intent(in) :: mask
    integer, dimension(n_prt) :: s
    integer :: h, hi
    p1 = - k(:,1) ! incoming
    p2 = - k(:,2) ! incoming
    p3 =   k(:,3) ! outgoing
    p4 =   k(:,4) ! outgoing
    p5 =   k(:,5) ! outgoing
    p6 =   k(:,6) ! outgoing
    p12 = p1 + p2
    p13 = p1 + p3
    p24 = p2 + p4
    p34 = p3 + p4
    p56 = p5 + p6
    p123 = p2 + p13
    p124 = p1 + p24
    p134 = p1 + p34
    p125 = p5 + p12
    p135 = p5 + p13
    p126 = p6 + p12
    p136 = p6 + p13
    amp = 0
    if (hel_finite == 0) return
    do hi = 1, hel_finite
      h = hel_map(hi)
      s = table_spin_states(:,h)
      owf_l1_1 = u (mass(11), - p1, s(1))
      owf_l1b_2 = vbar (mass(11), - p2, s(2))
      owf_l1b_3 = ubar (mass(11), p3, s(3))
      owf_n1_4 = v (mass(12), p4, s(4))
      owf_u1b__1_5 = ubar (mass(2), p5, s(5))
      owf_d1_1__6 = v (mass(1), p6, s(6))
      call compute_fusions_0001 ()
      call compute_fusions_0002 ()
      call compute_brakets_0001 ()
      amp(1,1,h) = oks_l1l1bl1n1bu1_1_d1b__1
    end do
  end subroutine calculate_amplitudes
  subroutine compute_fusions_0001 ()
      owf_a_12 = pr_feynman(p12, + v_ff(qlep,owf_l1b_2,owf_l1_1))
      owf_z_12 = pr_unitarity(p12,mass(23),wd_tl(p12,width(23)),.false., &
         + va_ff(gnclep(1),gnclep(2),owf_l1b_2,owf_l1_1))
      owf_a_13 = pr_feynman(p13, + v_ff(qlep,owf_l1b_3,owf_l1_1))
      owf_z_13 = pr_unitarity(p13,mass(23),wd_tl(p13,width(23)),.false., &
         + va_ff(gnclep(1),gnclep(2),owf_l1b_3,owf_l1_1))
      owf_wp_24 = pr_unitarity(p24,mass(24),wd_tl(p24,width(24)),.false., &
         + vl_ff(gcc,owf_l1b_2,owf_n1_4))
      owf_wp_34 = pr_unitarity(p34,mass(24),wd_tl(p34,width(24)),.false., &
         + vl_ff(gcc,owf_l1b_3,owf_n1_4))
      owf_wm_56 = pr_unitarity(p56,mass(24),wd_tl(p56,width(24)),.false., &
         + vl_ff(gcc,owf_u1b__1_5,owf_d1_1__6))
      owf_l1b_123 = pr_psibar(p123,mass(11),wd_tl(p123,width(11)),.false., &
         + f_fv(qlep,owf_l1b_2,owf_a_13) &
         + f_fva(gnclep(1),gnclep(2),owf_l1b_2,owf_z_13) &
         - f_fv(qlep,owf_l1b_3,owf_a_12) &
         - f_fva(gnclep(1),gnclep(2),owf_l1b_3,owf_z_12))
      owf_n1_124 = pr_psi(p124,mass(12),wd_tl(p124,width(12)),.false., &
         + f_vlf(gcc,owf_wp_24,owf_l1_1) &
         - f_vaf(gncneu(1),gncneu(2),owf_z_12,owf_n1_4))
      owf_n1_134 = pr_psi(p134,mass(12),wd_tl(p134,width(12)),.false., &
         + f_vlf(gcc,owf_wp_34,owf_l1_1) &
         - f_vaf(gncneu(1),gncneu(2),owf_z_13,owf_n1_4))
  end subroutine compute_fusions_0001
  subroutine compute_fusions_0002 ()
      owf_u1b__1_125 = pr_psibar(p125,mass(2),wd_tl(p125,width(2)),.false., &
         - f_fv(qup,owf_u1b__1_5,owf_a_12) &
         - f_fva(gncup(1),gncup(2),owf_u1b__1_5,owf_z_12))
      owf_u1b__1_135 = pr_psibar(p135,mass(2),wd_tl(p135,width(2)),.false., &
         - f_fv(qup,owf_u1b__1_5,owf_a_13) &
         - f_fva(gncup(1),gncup(2),owf_u1b__1_5,owf_z_13))
      owf_d1_1__126 = pr_psi(p126,mass(1),wd_tl(p126,width(1)),.false., &
         - f_vf(qdwn,owf_a_12,owf_d1_1__6) &
         - f_vaf(gncdwn(1),gncdwn(2),owf_z_12,owf_d1_1__6))
      owf_d1_1__136 = pr_psi(p136,mass(1),wd_tl(p136,width(1)),.false., &
         - f_vf(qdwn,owf_a_13,owf_d1_1__6) &
         - f_vaf(gncdwn(1),gncdwn(2),owf_z_13,owf_d1_1__6))
  end subroutine compute_fusions_0002
  subroutine compute_brakets_0001 ()
      oks_l1l1bl1n1bu1_1_d1b__1 = 0
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_l1b_123*( &
         - f_vlf(gcc,owf_wm_56,owf_n1_4))
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + ( &
         + f_fvl(gcc,owf_l1b_2,owf_wm_56))*owf_n1_134
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + ( &
         - f_fvl(gcc,owf_l1b_3,owf_wm_56))*owf_n1_124
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_u1b__1_135* &
        ( + f_vlf(gcc,owf_wp_24,owf_d1_1__6))
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_u1b__1_125* &
        ( - f_vlf(gcc,owf_wp_34,owf_d1_1__6))
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + ( &
         + f_fvl(gcc,owf_u1b__1_5,owf_wp_24))*owf_d1_1__136
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + ( &
         - f_fvl(gcc,owf_u1b__1_5,owf_wp_34))*owf_d1_1__126
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_a_12*( &
         + g_gg(iqw,owf_wm_56,p56,owf_wp_34,p34))
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_a_13*( &
         - g_gg(iqw,owf_wm_56,p56,owf_wp_24,p24))
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_z_12*( &
         + g_gg(igzww,owf_wm_56,p56,owf_wp_34,p34))
      oks_l1l1bl1n1bu1_1_d1b__1 = oks_l1l1bl1n1bu1_1_d1b__1 + owf_z_13*( &
         - g_gg(igzww,owf_wm_56,p56,owf_wp_24,p24))
      oks_l1l1bl1n1bu1_1_d1b__1 = &
         - oks_l1l1bl1n1bu1_1_d1b__1 ! 4 vertices, 3 propagators
      ! unit symmetry factor
  end subroutine compute_brakets_0001

end module omega_amplitude
