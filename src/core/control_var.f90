!--------------------------------------------------------------------------------------------------
! MODULE: control_var
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date November, 2016
!
! DESCRIPTION: 
!> @brief Module containing program control variables.
!--------------------------------------------------------------------------------------------------
module control_var
    use global_defs
    use rattle_mod, only : constraint
    implicit none

    private
    public :: ctrl

    !----------------------------------------------------------------------------------------------
    ! TYPE: control
    !> @brief Holds program control variables.
    !----------------------------------------------------------------------------------------------
    type control

        !------------------------------------------------------------------------------------------
        ! Input/Output options.
        !------------------------------------------------------------------------------------------
        character(len=:), allocatable :: maindir !< Program working directory.
        character(len=:), allocatable :: output_dir !< Main output directory name.

        logical :: qext !< Flag to call external program for QM calculations.
        character(len=:), allocatable :: qprog !< Program for QM calculations.
        character(len=:), allocatable :: mprog !< Program for MM calculations.
        character(len=:), allocatable :: oprog !< Program for overlap calculations.

        character(len=:), allocatable :: qmdir !< Work directory for QM calculation.
        character(len=:), allocatable :: bufile !< Backup file name.
        integer :: buinterval !< Number of steps between backup files.

        logical :: print(50) = .false. !< Output options (for writestep subroutine).

        !------------------------------------------------------------------------------------------
        ! Program flow control.
        !------------------------------------------------------------------------------------------
        integer :: seed(1) = [-1] !< Seed for the random number generator. A new value can be given 
        !! in the input file. If not, random_seed is called to generate a new seed.
        logical :: restart = .false. !< Restart from backup of previous run.
        real(dp) :: min01gap = -10000.0_dp !< Threshold to stop program in case of S0/S1 CI.
        integer :: target_state = -1 !< Stop dynamics after reaching target state.
        real(dp) :: target_state_time = 0.0_dp !< Additional time to propagate after reaching the
        !! target state.
        real(dp) :: max_tot_en_change_step = 100.0_dp !< Maximum change in total energy allowed
        !! during a single step.
        real(dp) :: max_tot_en_change = 100.0_dp !< Maximum change in total energy allowed over the
        !! course of a trajectory.
        real(dp) :: t0_tot_en !< Total energy at the beginning of the calculation.

        !------------------------------------------------------------------------------------------
        ! Nuclear dynamics options.
        !------------------------------------------------------------------------------------------
        real(dp) :: max_time !< Total time for the simulation.
        real(dp) :: dt !< Time step for classical dynamics.
        real(dp) :: dt_0 !< Initial dt.
        integer :: orientlvl !< Orientation of the molecule during dynamics:
        !! - 0 - None.
        !! - 1 - Project translation.
        !! - 2 - Project translation and rotation.
        logical :: constrain_dyn = .false. !< Signal to enable constraints in dynamics. (RATTLE)
        integer :: n_cns = 0 !< Number of constraints.
        type(constraint), allocatable :: cns(:) !< Constraints.

        !------------------------------------------------------------------------------------------
        ! QMMM options.
        !------------------------------------------------------------------------------------------
        logical :: qm !< Signal to enable QM calculation.
        logical :: mm !< Signal to enable MM calculation.
        real(dp) :: mmcut !< Cutoff distance for calculating interactions between the QM and MM
        !! systems. If the distance between two atoms is greater than this value, their interaction
        !! will not affect the energy of the system.
        logical :: pcharge !< Signal to use partial charges in the calculations.
        logical :: pbc !< Signal to use periodic boundary conditions in the MM calculation.

        !------------------------------------------------------------------------------------------
        ! Surface hopping options.
        !------------------------------------------------------------------------------------------
        logical :: hop !< Signal that a hop occured in the current step.
        integer :: sh !< Type of surface hopping calculation:
                      !! - 0 - No hopping.
                      !! - 1 - Landau-Zener formula.
                      !! - 2 - Fewest switches surface hopping in the adiabatic representation.
                      !! - 3 - FSSH in the locally diabatic representation.
        integer :: shnstep !< Number of steps in the integration of the TDSE.

        logical, allocatable :: couple(:) !< Signals that a state is included in the coupling
        !! calculation. When .false., the state is skipped regardless of other options.

        integer :: decohlvl !< Method of decoherence correction during dynamics:
        !!  - 0 - No decoherence correction.
        !!  - 1 - Non linear decay of mixing (Granucci, Persico; doi: 10.1063/1.2715585).

        integer :: couplvl !< Method of reduction of the coupling matrix calculation. At each step,
        !! the cmask logical array is created based on the value of this variable:
        !! - 0 - All matrix elements are calculated.
        !! - 1 - Only ctrl%coupndiff states around the current state are included.
        !! - 2 - Only states in a ctrl%coupediff energy window of the current state are included.
        !! - 3 - Couplings between each state and ctrl%coupndiff states around it are calculated.
        !! - 4 - Couplings between each states and the states in a ctrl%coupediff energy window
        !!       around it are calculated.
        integer :: coupndiff !< Number of states for couplvl 1 and 3.
        real(dp) :: coupediff !< Energy difference for couplvl 2 and 4.

        integer :: tdc_type !< Method for calculating the time-derivative coupling matrix:
        !! - 1 - Use wave function overlaps.
        !! - 2 - Use nonadiabatic coupling vectors.

        integer :: tdc_interpolate !< Method for interpolating tdcs during time step.
        !! - 1 - Finite differences method to calculate tdc at t0 + dt/2.
        !! - 2 - Linear interpolation of finite differences tdcs at t0 - dt/2 and t0 + dt/2.
        !! - 4 - Norm-preserving interpolation of Meek and Levine. (10.1021/jz5009449).

        integer :: ene_interpolate !< Method for interpolating energies during time step.
        !! - 0 - Constant energies equal to E(t + dt).
        !! - 1 - Energies equal to E(t) until t + dt/2, and to E(t + dt) afterwards.
        !! - 2 - Linear interpolation between E(t) and E(t + dt)

        integer :: vrescale !< Method for rescaling velocities after hop.
        !! - 0 - No velocity rescaling.
        !! - 1 - Rescale along velocity vector.
        !! - 2 - Rescale along gradient difference vector.
        !! - 3 - Rescale along nonadiabatic coupling vector.

        integer :: fhop !< Method for treating frustrated hops.
        !! - 1 - Return to previous state without changing velocity.
        !! - 2 - Return to previous state and invert the velocity along the relevant direction.

        integer :: phaselvl !< Method of handling phase of wave functions during dynamics.
        !! - 0 - None.
        !! - 1 - Match phase of adiabatic states.
        !! - 2 - Match phase of diabatic states (using assignment problem solution).

        real(dp) :: lz_prob_conv = 0.05_dp !< Convergence for estimated error in LZ probability.
        real(dp) :: lz_min_dt = 0.0_dp !< Minimum time step when converging LZ probability.
        real(dp) :: qm_en_err = 1.0e-6 !< Estimated error in calculated energies.

        logical :: oscill !< Signal to read oscillator strengths from the QM calculation.

    end type control

    type(control) :: ctrl
    
end module control_var
