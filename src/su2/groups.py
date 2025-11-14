# -*- coding: utf-8 -*-
from dataclasses import dataclass, field
from .enums import (
    SolverType,
    TurbulenceModel,
    ShearStressTransportModel,
    SpalartAllmarasModel,
    TransitionModel,
    LmTransitionModelOptions,
    SgsModel,
    SolutionVerification,
    MathProblem,
    UnitSystem,
)

@dataclass
class ProblemDefinition:
    """ Direct, adjoint, and linearized problem definition.

    Attributes
    ----------
    solver : SolverType
        Solver type.
    kind_turb_model : TurbulenceModel
        Turbulence model.
    sst_options : ShearStressTransportModel
        SST model versions/corrections.
    sa_options : SpalartAllmarasModel
        SA model versions/corrections.
    kind_trans_model : TransitionModel
        Transition model.
    hroughness : float
        RMS roughness for transition model.
    lm_options : LmTransitionModelOptions
        LM model versions/corrections.
    kind_sgs_model : SgsModel
        Subgrid scale model.
    kind_verification_solution : SolutionVerification
        Verification solution type.
    math_problem : MathProblem
        Mathematical problem type.
    axisymmetric : bool
        Axisymmetric simulation for 2D problems.
    gravity_force : bool
        Enable gravity force.
    restart_sol : bool
        Restart solution.
    wrt_restart_compact : bool
        Save only minimum required variables for restart.
    discard_infiles : bool
        Discard data stored in solution and geometry files.
    system_measurements : UnitSystem
        System of measurements.
    config_list : list[str]
        List of config files for multizone setup.
    """
    solver: SolverType
    kind_turb_model: TurbulenceModel       = TurbulenceModel.NONE
    sst_options: ShearStressTransportModel = ShearStressTransportModel.NONE
    sa_options: SpalartAllmarasModel       = SpalartAllmarasModel.NONE
    kind_trans_model: TransitionModel      = TransitionModel.NONE
    hroughness: float                      = 1.0e-6
    lm_options: LmTransitionModelOptions   = LmTransitionModelOptions.NONE
    kind_sgs_model: SgsModel               = SgsModel.NONE
    kind_verification_solution: SolutionVerification = SolutionVerification.NO_VERIFICATION_SOLUTION
    math_problem: MathProblem              = MathProblem.DIRECT
    axisymmetric: bool                     = False
    gravity_force: bool                    = False
    restart_sol: bool                      = False
    wrt_restart_compact: bool              = True
    discard_infiles: bool                  = False
    system_measurements: UnitSystem        = UnitSystem.SI
    config_list: list[str]                 = field(default_factory=lambda: [])
