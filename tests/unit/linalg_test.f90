program linalg_test
    implicit none
    integer :: ret


    ret = gemm_test()
    if (ret /= 0) then
        write(*,*) "Error in gemm_test."
        stop ret
    end if

  ! ret = gesvd_test()
  ! if (ret /= 0) then
  !     write(*,*) "Error in gesvd_test."
  !     stop ret
  ! end if

  ! ret = syev_test()
  ! if (ret /= 0) then
  !     write(*,*) "Error in syev_test."
  !     stop ret
  ! end if

    contains
    
    

    integer function gemm_test() result(ret)
        use global_defs, only : dp
        use linalg_wrapper_mod, only : gemm
        real(dp) :: mat1(2, 2)
        real(dp) :: mat2(2, 2)
        real(dp) :: tiny = 1.0e-15

        mat1(1, :) = [1.0_dp, 2.0_dp]
        mat1(2, :) = [3.0_dp, 4.0_dp]

        ret = 1
        call gemm(mat1, mat1, mat2)
        if (abs(mat2(2,2)-22.0_dp) > tiny) return
        if (abs(mat2(1,2)-10.0_dp) > tiny) return
        ret = 0
    end function gemm_test


    
end program linalg_test
