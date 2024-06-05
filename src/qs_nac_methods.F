!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_nac_methods
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              nac_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                              cp_fm_scale_and_add
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_init_random,&
                                              cp_fm_release,&
                                              cp_fm_set_element,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE input_constants,                 ONLY: cholesky_restore
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              duplicate_mo_set,&
                                              init_mo_set,&
                                              mo_set_type,&
                                              reassign_allocated_mos
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_nac_methods'

   PUBLIC :: nac

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE nac(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac'

      INTEGER                                            :: handle
      LOGICAL                                            :: initial
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: nac_matrix
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: current_mo
      TYPE(nac_control_type), POINTER                    :: nac_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: input, nac_section, print_section

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (input, dft_control, scf_env, current_mo)
      CALL get_qs_env(qs_env, &
                      input=input, &
                      dft_control=dft_control, &
                      scf_env=scf_env)

      NULLIFY (nac_control)
      nac_control => dft_control%nac_control

      NULLIFY (nac_section, print_section)
      nac_section => section_vals_get_subs_vals(input, "PROPERTIES%NAC")
      print_section => section_vals_get_subs_vals(nac_section, "PRINT")

      initial = .FALSE.

      IF (.NOT. ASSOCIATED(qs_env%last_mo)) THEN
         initial = .TRUE.
         CALL nac_init(qs_env, scf_env, nac_control)
      END IF

      IF (.NOT. initial .AND. nac_control%dt .LT. 0) THEN
         nac_control%dt = cp_unit_from_cp2k(qs_env%sim_time/qs_env%sim_step, "fs")
      END IF

      IF (nac_control%do_eigen) CALL nac_compute_mo(qs_env, current_mo)

      IF (.NOT. initial) THEN
         CALL nac_compute(qs_env, nac_control, current_mo, nac_matrix)
         CALL nac_output(print_section, nac_matrix)
      END IF

      IF (nac_control%do_eigen) CALL nac_save_last_mo(qs_env, current_mo)

      IF (.NOT. initial) DEALLOCATE (nac_matrix)

      CALL timestop(handle)

   END SUBROUTINE nac

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param scf_env ...
!> \param nac_control ...
! **************************************************************************************************
   SUBROUTINE nac_init(qs_env, scf_env, nac_control)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(nac_control_type), POINTER                    :: nac_control

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_init'

      INTEGER                                            :: handle
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: last_mo, mos

      CALL timeset(routineN, handle)

      NULLIFY (last_mo)
      CALL get_qs_env(qs_env, &
                      mos=mos, &
                      last_mo=last_mo)
      IF (scf_env%method .EQ. ot_method_nr .AND. nac_control%vir_active_space .GT. 0) THEN
         nac_control%do_eigen = .TRUE.
      ELSE IF (scf_env%method .NE. ot_method_nr) THEN
         CPABORT("NAC calculations only implemented for OT method now.")
      ELSE IF (nac_control%vir_active_space .LE. 0) THEN
         CPABORT("No virtial orbital No NAC calculation.")
      END IF

      CALL nac_init_last_mo(qs_env, last_mo)

      IF (nac_control%occ_active_space .LE. 0) nac_control%occ_active_space = mos(1)%homo
      nac_control%occ_active_space = MIN(nac_control%occ_active_space, mos(1)%homo)
      IF (nac_control%vir_active_space .LE. 0) nac_control%vir_active_space = 0
      nac_control%vir_active_space = MIN(nac_control%vir_active_space, mos(1)%nao - mos(1)%homo)

      nac_control%dt = -1

      CALL timestop(handle)

   END SUBROUTINE nac_init

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param last_mo ...
! **************************************************************************************************
   SUBROUTINE nac_init_last_mo(qs_env, last_mo)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: last_mo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_init_last_mo'

      INTEGER                                            :: handle, ispin
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(nac_control_type), POINTER                    :: nac_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, mos)
      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      para_env=para_env, &
                      dft_control=dft_control, &
                      mos=mos)

      NULLIFY (nac_control)
      nac_control => dft_control%nac_control

      ALLOCATE (last_mo(SIZE(mos)))

      DO ispin = 1, SIZE(mos)
         CALL cp_fm_struct_create(fm_struct_tmp, &
                                  context=blacs_env, &
                                  para_env=para_env, &
                                  nrow_global=mos(ispin)%nao, &
                                  ncol_global=mos(ispin)%homo + nac_control%vir_active_space)
         CALL allocate_mo_set(last_mo(ispin), &
                              mos(ispin)%nao, &
                              mos(ispin)%homo + nac_control%vir_active_space, &
                              mos(ispin)%nelectron, &
                              mos(ispin)%n_el_f, &
                              mos(ispin)%maxocc, &
                              mos(ispin)%flexible_electron_count)
         CALL init_mo_set(mo_set=last_mo(ispin), &
                          fm_struct=fm_struct_tmp, &
                          name="For saving last MO")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_init_random(last_mo(ispin)%mo_coeff)
      END DO

      CALL set_qs_env(qs_env, last_mo=last_mo)

      CALL timestop(handle)

   END SUBROUTINE nac_init_last_mo

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param current_mo ...
! **************************************************************************************************
   SUBROUTINE nac_compute_mo(qs_env, current_mo)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: current_mo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_compute_mo'

      INTEGER                                            :: handle, ispin, cholesky_method
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fm_struct
      TYPE(cp_fm_type)                                   :: matrix_ks_fm, ortho_fm, work_fm
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: last_mo
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, matrix_ks, matrix_s, current_mo)
      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      scf_env=scf_env, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      last_mo=last_mo)

      ALLOCATE (current_mo(SIZE(last_mo)))

      CALL cp_fm_struct_create(ao_ao_fm_struct, nrow_global=last_mo(1)%nao, ncol_global=last_mo(1)%nao, context=blacs_env)

      DO ispin = 1, SIZE(last_mo)
         CALL duplicate_mo_set(current_mo(ispin), last_mo(ispin))

         CALL cp_fm_create(matrix_ks_fm, ao_ao_fm_struct)
         CALL cp_fm_create(ortho_fm, ao_ao_fm_struct)
         CALL cp_fm_create(work_fm, ao_ao_fm_struct)

         CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, ortho_fm)
         CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix, matrix_ks_fm)
         CALL cp_fm_cholesky_decompose(ortho_fm)
         cholesky_method = cholesky_restore
         CALL eigensolver(matrix_ks_fm=matrix_ks_fm, mo_set=current_mo(ispin), ortho=ortho_fm, &
                          work=work_fm, cholesky_method=cholesky_method, &
                          do_level_shift=.FALSE., level_shift=0.0_dp, use_jacobi=.FALSE.)

         CALL cp_fm_release(work_fm)
         CALL cp_fm_release(ortho_fm)
         CALL cp_fm_release(matrix_ks_fm)
      END DO

      CALL timestop(handle)

   END SUBROUTINE nac_compute_mo

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nac_control ...
!> \param current_mo ...
!> \param nac_matrix ...
!> \return ...
! **************************************************************************************************
   SUBROUTINE nac_compute(qs_env, nac_control, current_mo, nac_matrix)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(nac_control_type), POINTER                    :: nac_control
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: current_mo
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: nac_matrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_compute'

      INTEGER                                            :: handle, imo, ispin, nspin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: current_eigenvalues, last_eigenvalues
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: current_mo_coeff, last_mo_coeff, &
                                                            nac_matrix_back
      TYPE(cp_fm_type), POINTER                          :: current_smo, last_smo
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: last_mo

      CALL timeset(routineN, handle)

      NULLIFY (last_mo)
      CALL get_qs_env(qs_env, &
                      last_mo=last_mo, &
                      matrix_s=matrix_s)

      nspin = SIZE(last_mo)
      ALLOCATE (last_mo_coeff(nspin), current_mo_coeff(nspin), nac_matrix_back(nspin), nac_matrix(nspin))

      CALL nac_get_active_space(qs_env, &
                                last_mo, &
                                nac_control%occ_active_space, &
                                nac_control%vir_active_space, &
                                last_mo_coeff, &
                                last_eigenvalues)
      CALL nac_get_active_space(qs_env, &
                                current_mo, &
                                nac_control%occ_active_space, &
                                nac_control%vir_active_space, &
                                current_mo_coeff, &
                                current_eigenvalues)

      ALLOCATE (current_smo, last_smo)

      DO ispin = 1, nspin
         CALL nac_get_smo(matrix_s(1)%matrix, last_mo_coeff(ispin), last_smo)
         CALL nac_get_smo(matrix_s(1)%matrix, current_mo_coeff(ispin), current_smo)
         CALL nac_compute_nac(last_smo, current_mo_coeff(ispin), nac_matrix(ispin))
         CALL nac_compute_nac(current_smo, last_mo_coeff(ispin), nac_matrix_back(ispin))
         CALL cp_fm_scale_and_add(0.658218_dp/2.0_dp/nac_control%dt/13.60569253_dp, &
                                  nac_matrix(ispin), &
                                  -0.658218_dp/2.0_dp/nac_control%dt/13.60569253_dp, &
                                  nac_matrix_back(ispin))
         DO imo = 1, nac_control%occ_active_space + nac_control%vir_active_space
            CALL cp_fm_set_element(nac_matrix(ispin), &
                                   imo, &
                                   imo, &
                                   (current_eigenvalues(ispin, imo) + last_eigenvalues(ispin, imo)) &
                                   *27.2114_dp/2.0_dp/13.60569253_dp)
         END DO
      END DO

      DEALLOCATE (last_mo_coeff, current_mo_coeff, nac_matrix_back, &
                  current_smo, last_smo, current_eigenvalues, last_eigenvalues)

      CALL timestop(handle)

   END SUBROUTINE nac_compute

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mo ...
!> \param occ_active_space ...
!> \param vir_active_space ...
!> \param mo_coeff ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE nac_get_active_space(qs_env, mo, occ_active_space, vir_active_space, mo_coeff, eigenvalues)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo
      INTEGER                                            :: occ_active_space, vir_active_space
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mo_coeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenvalues

      CHARACTER(LEN=*), PARAMETER :: routineN = 'nac_get_active_space'

      INTEGER                                            :: handle, ispin, nao, nmo, nspin
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      blacs_env=blacs_env, &
                      para_env=para_env)

      nao = mo(1)%nao
      nmo = occ_active_space + vir_active_space
      nspin = SIZE(mo)

      CALL cp_fm_struct_create(fm_struct_tmp, &
                               context=blacs_env, &
                               para_env=para_env, &
                               nrow_global=nao, &
                               ncol_global=nmo)

      ALLOCATE (eigenvalues(nspin, nmo))

      DO ispin = 1, nspin
         CALL cp_fm_create(mo_coeff(ispin), fm_struct_tmp)
         CALL cp_fm_to_fm(mo(ispin)%mo_coeff, &
                          mo_coeff(ispin), &
                          occ_active_space + vir_active_space, &
                          mo(ispin)%nmo - vir_active_space - occ_active_space + 1, &
                          1)
         eigenvalues(ispin, :) = mo(ispin)%eigenvalues(mo(ispin)%nmo - vir_active_space - occ_active_space + 1:)
      END DO

      CALL timestop(handle)

   END SUBROUTINE nac_get_active_space

! **************************************************************************************************
!> \brief ...
!> \param matrix_s ...
!> \param mo_coeff ...
!> \param smo ...
! **************************************************************************************************
   SUBROUTINE nac_get_smo(matrix_s, mo_coeff, smo)
      TYPE(dbcsr_type), POINTER                          :: matrix_s
      TYPE(cp_fm_type)                                   :: mo_coeff, smo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_get_smo'

      INTEGER                                            :: handle, nao, nmo
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(mo_coeff, &
                          nrow_global=nao, &
                          ncol_global=nmo, &
                          context=blacs_env, &
                          para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_tmp, &
                               context=blacs_env, &
                               para_env=para_env, &
                               nrow_global=nao, &
                               ncol_global=nmo)
      CALL cp_fm_create(smo, fm_struct_tmp)
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s, &
                                   mo_coeff, &
                                   smo, &
                                   ncol=nmo, &
                                   alpha=1.0_dp, &
                                   beta=0.0_dp)
      CALL timestop(handle)

   END SUBROUTINE nac_get_smo

! **************************************************************************************************
!> \brief ...
!> \param smo ...
!> \param mo ...
!> \param nac_matrix ...
! **************************************************************************************************
   SUBROUTINE nac_compute_nac(smo, mo, nac_matrix)
      TYPE(cp_fm_type)                                   :: smo, mo, nac_matrix

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_compute_nac'

      INTEGER                                            :: handle, nao, nmo
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      CALL cp_fm_get_info(mo, &
                          nrow_global=nao, &
                          ncol_global=nmo, &
                          context=blacs_env, &
                          para_env=para_env)
      CALL cp_fm_struct_create(fm_struct_tmp, &
                               context=blacs_env, &
                               para_env=para_env, &
                               nrow_global=nmo, &
                               ncol_global=nmo)
      CALL cp_fm_create(nac_matrix, fm_struct_tmp)
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_fm_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, smo, mo, 0.0_dp, nac_matrix)

      CALL timestop(handle)
   END SUBROUTINE nac_compute_nac

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param current_mo ...
! **************************************************************************************************
   SUBROUTINE nac_save_last_mo(qs_env, current_mo)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: current_mo

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_save_last_mo'

      INTEGER                                            :: handle, ispin
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: last_mo

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, last_mo=last_mo)
      DO ispin = 1, SIZE(current_mo)
         CALL reassign_allocated_mos(last_mo(ispin), current_mo(ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE nac_save_last_mo

! **************************************************************************************************
!> \brief ...
!> \param print_section ...
!> \param mos ...
! **************************************************************************************************
   SUBROUTINE nac_output(print_section, mos)
      TYPE(section_vals_type), POINTER                   :: print_section
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'nac_output'

      INTEGER                                            :: handle, ispin, unit
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()

      IF (cp_print_key_should_output(logger%iter_info, print_section) /= 0) THEN
         unit = cp_print_key_unit_nr(logger, print_section, extension=".nacLog", file_status="REPLACE")
         DO ispin = 1, SIZE(mos)
            CALL cp_fm_write_formatted(mos(ispin), unit)
         END DO
         CALL cp_print_key_finished_output(unit, logger, print_section)
      END IF

      CALL timestop(handle)

   END SUBROUTINE nac_output

END MODULE qs_nac_methods
