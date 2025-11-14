# -*- coding: utf-8 -*-
from dataclasses import dataclass, field
from .enums import (
    YesNoEnum,
    SolverType,
    TurbulenceModel,
    ShearStressTransportModel,
    SpalartAllmarasModel,
    TransitionModel,
    LmTransitionModelOptions,
    SgsModel,
    Verification,
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
    turbulence_model : TurbulenceModel
        Turbulence model.
    sst_options : ShearStressTransportModel
        SST model versions/corrections.
    sa_options : SpalartAllmarasModel
        SA model versions/corrections.
    transition_model : TransitionModel
        Transition model.
    hroughness : float
        RMS roughness for transition model.
    lm_options : LmTransitionModelOptions
        LM model versions/corrections.
    sgs_model : SgsModel
        Subgrid scale model.
    solution_verification : Verification
        Verification solution type.
    math_problem : MathProblem
        Mathematical problem type.
    axisymmetric : YesNoEnum
        Axisymmetric simulation for 2D problems.
    gravity_force : YesNoEnum
        Enable gravity force.
    restart_sol : YesNoEnum
        Restart solution.
    wrt_restart_compact : YesNoEnum
        Save only minimum required variables for restart.
    discard_infiles : YesNoEnum
        Discard data stored in solution and geometry files.
    system_measurements : UnitSystem
        System of measurements.
    config_list : list[str]
        List of config files for multizone setup.
    """
    solver: SolverType                     = SolverType.EULER
    turbulence_model: TurbulenceModel      = TurbulenceModel.NONE
    sst_options: ShearStressTransportModel = ShearStressTransportModel.NONE
    sa_options: SpalartAllmarasModel       = SpalartAllmarasModel.NONE
    transition_model: TransitionModel      = TransitionModel.NONE
    hroughness: float                      = 1.0e-6
    lm_options: LmTransitionModelOptions   = LmTransitionModelOptions.NONE
    sgs_model: SgsModel                    = SgsModel.NONE
    solution_verification: Verification    = Verification.NO_VERIFICATION_SOLUTION
    math_problem: MathProblem              = MathProblem.NONE
    axisymmetric: YesNoEnum                = YesNoEnum.NO
    gravity_force: YesNoEnum               = YesNoEnum.NO
    restart_sol: YesNoEnum                 = YesNoEnum.NO
    wrt_restart_compact: YesNoEnum         = YesNoEnum.NO
    discard_infiles: YesNoEnum             = YesNoEnum.NO
    system_measurements: UnitSystem        = UnitSystem.NONE
    config_list: list[str]                 = field(default_factory=lambda: [])

    def _handle_turbulence_model(self, cfg: list[str]) -> None:
        """ Handle turbulence model related settings. """
        cfg.append(f"KIND_TURB_MODEL= {self.turbulence_model.value}")

        match self.turbulence_model:
            case TurbulenceModel.SST:
                if self.sst_options != ShearStressTransportModel.NONE:
                    cfg.append(f"SST_OPTIONS= {self.sst_options.value}")
            case TurbulenceModel.SA:
                if self.sa_options != SpalartAllmarasModel.NONE:
                    cfg.append(f"SA_OPTIONS= {self.sa_options.value}")
            case _:
                pass

    def _handle_transition_model(self, cfg: list[str]) -> None:
        """ Handle transition model related settings. """
        cfg.append(f"KIND_TRANS_MODEL= {self.transition_model.value}")
        cfg.append(f"HROUGHNESS= {self.hroughness}")

        match self.transition_model:
            case TransitionModel.LM:
                if self.lm_options != LmTransitionModelOptions.NONE:
                    cfg.append(f"LM_OPTIONS= {self.lm_options.value}")
            case _:
                pass

    def to_cfg(self) -> str:
        """ Generate configuration file entries for problem definition.

        Returns
        -------
        str
            Configuration file entries.
        """
        cfg = []

        header = "DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION"
        cfg.append(f"{80 * '%'}\n%% {header}\n{80 * '%'}")

        cfg.append(f"SOLVER= {self.solver.value}")

        if self.turbulence_model != TurbulenceModel.NONE:
            self._handle_turbulence_model(cfg)

        if self.transition_model != TransitionModel.NONE:
            self._handle_transition_model(cfg)

        if self.sgs_model != SgsModel.NONE:
            cfg.append(f"KIND_SGS_MODEL= {self.sgs_model.value}")

        if self.solution_verification != Verification.NO_VERIFICATION_SOLUTION:
            cfg.append(f"KIND_VERIFICATION_SOLUTION= {self.solution_verification.value}")

        if self.math_problem != MathProblem.NONE:
            cfg.append(f"MATH_PROBLEM= {self.math_problem.value}")

        if self.axisymmetric != YesNoEnum.NO:
            cfg.append(f"AXISYMMETRIC= {self.axisymmetric.value}")

        if self.gravity_force != YesNoEnum.NO:
            cfg.append(f"GRAVITY_FORCE= {self.gravity_force.value}")

        if self.restart_sol != YesNoEnum.NO:
            cfg.append(f"RESTART_SOL= {self.restart_sol.value}")

        if self.wrt_restart_compact != YesNoEnum.NO:
            cfg.append(f"WRT_RESTART_COMPACT= {self.wrt_restart_compact.value}")

        if self.discard_infiles != YesNoEnum.NO:
            cfg.append(f"DISCARD_INFILES= {self.discard_infiles.value}")

        if self.system_measurements != UnitSystem.NONE:
            cfg.append(f"SYSTEM_MEASUREMENTS= {self.system_measurements.value}")

        if self.config_list:
            cfg.append(f"CONFIG_LIST= {','.join(self.config_list)}")

        return "\n\n".join(cfg)
