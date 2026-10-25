module pot_test_mod
    use global_defs, only : dp
    use evaluator_diabatic_mod, only : diabatic_evaluator
    implicit none

    !> Evaluator without a potential, used to test the state group setup.
    type, extends(diabatic_evaluator) :: group_evaluator
    contains
        procedure :: update_geometry => group_update_geometry
        procedure :: eval_diab => group_eval_diab
    end type group_evaluator

contains

    subroutine group_update_geometry(self, geometry, coordinates)
        class(group_evaluator), intent(inout) :: self
        real(dp), intent(in) :: geometry(:, :)
        character(len=*), intent(in) :: coordinates
    end subroutine group_update_geometry

    subroutine group_eval_diab(self)
        class(group_evaluator), intent(inout) :: self
    end subroutine group_eval_diab

end module pot_test_mod


program pot_test
    implicit none
    integer :: ret


    ret = group_test()
    if (ret /= 0) then
        write(*,*) "Error in group_test."
        stop ret
    end if

    ret = vc_group_test()
    if (ret /= 0) then
        write(*,*) "Error in vc_group_test."
        stop ret
    end if

    contains


    integer function group_test() result(ret)
        use pot_test_mod, only : group_evaluator
        use evaluator_base_mod, only : state_location
        type(group_evaluator) :: pot
        type(state_location) :: loc

        ret = 1
        ! 4 singlets and 3 triplets.
        call pot%set_group_nstate([4, 3], [1, 3])
        if (pot%n_state /= 13) return
        if (pot%n_group /= 4) return
        if (any(pot%group_nstate /= [4, 3, 3, 3])) return
        if (any(pot%group_i0 /= [1, 5, 8, 11])) return
        if (any(pot%group_mult /= [1, 3, 3, 3])) return
        if (any(pot%group_block /= [1, 2, 2, 2])) return

        ret = 2
        loc = pot%locate_state(4)
        if (loc%group /= 1 .or. loc%block /= 1 .or. loc%mult /= 1) return
        if (loc%component /= 1 .or. loc%index /= 4) return
        loc = pot%locate_state(5)
        if (loc%group /= 2 .or. loc%block /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 1 .or. loc%index /= 1) return
        loc = pot%locate_state(10)
        if (loc%group /= 3 .or. loc%block /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 2 .or. loc%index /= 3) return
        loc = pot%locate_state(13)
        if (loc%group /= 4 .or. loc%block /= 2 .or. loc%mult /= 3) return
        if (loc%component /= 3 .or. loc%index /= 3) return
        ret = 0
    end function group_test


    integer function vc_group_test() result(ret)
        use global_defs, only : dp
        use vc_model_mod, only : vc_model
        use vc_evaluator_mod, only : vc_evaluator
        use evaluator_base_mod, only : state_location
        type(vc_model), target :: model
        type(vc_evaluator) :: pot
        type(vc_evaluator) :: pot_bad
        type(state_location) :: loc
        real(dp) :: soc(8, 8)

        ! 2 singlets and 2 triplets. Blocks are numbered by multiplicity, also when there are no
        ! doublets.
        model%nstate = [2, 0, 2]
        model%nmode = 1
        soc = 0.0_dp
        ! Singlet-triplet couplings only.
        soc(1, 3) = 1.0e-4_dp
        soc(2, 6) = 2.0e-4_dp
        soc(1, 8) = 3.0e-4_dp
        soc = soc + transpose(soc)
        call model%soc%add_zero(model%nmode, soc)

        ret = 1
        call pot%initialize(model)
        if (pot%n_state /= 8) return
        if (any(pot%group_nstate /= [2, 2, 2, 2])) return
        if (any(pot%group_i0 /= [1, 3, 5, 7])) return
        if (any(pot%group_mult /= [1, 3, 3, 3])) return
        if (any(pot%group_block /= [1, 3, 3, 3])) return
        loc = pot%locate_state(6)
        if (loc%block /= 3 .or. loc%component /= 2 .or. loc%index /= 2) return

        ret = 2
        if (pot%soc_in_group) return

        ! Coupling between the two states of the second triplet component.
        ret = 3
        soc(5, 6) = 1.0e-4_dp
        soc(6, 5) = 1.0e-4_dp
        call model%soc%add_zero(model%nmode, soc)
        call pot_bad%initialize(model)
        if (.not. pot_bad%soc_in_group) return
        ret = 0
    end function vc_group_test


end program pot_test
