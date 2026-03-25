module potential_evaluator_mod
    use global_defs
    implicit none

    type, abstract :: potential_evaluator
    contains
        procedure(update_geometry), deferred :: update_geometry
    end type potential_evaluator

    abstract interface
        subroutine update_geometry(self, geometry)
            use global_defs, only : dp
            import potential_evaluator
            class(potential_evaluator) :: self
            real(dp), intent(in) :: geometry(:, :)
        end subroutine update_geometry
    end interface

end module potential_evaluator_mod