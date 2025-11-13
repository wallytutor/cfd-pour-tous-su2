# -*- coding: utf-8 -*-
from dataclasses import dataclass, field
from .enums import (
    SolverType,
    TurbulenceModel,
    ShearStressTransportModel,
    SpalartAllmarasModel
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
    """
    solver: SolverType
    kind_turb_model: TurbulenceModel       = TurbulenceModel.NONE
    sst_options: ShearStressTransportModel = ShearStressTransportModel.NONE
    sa_options: SpalartAllmarasModel       = SpalartAllmarasModel.NONE

#     kind_trans_model: str = "NONE"
#     """Transition model (NONE, LM)"""
#     hroughness: float = 1.0e-6
#     """RMS roughness for transition model"""
#     lm_options: str = "NONE"
#     """LM model versions/correlations"""
#     kind_sgs_model: str = "NONE"
#     """Subgrid scale model"""
#     kind_verification_solution: str = "NO_VERIFICATION_SOLUTION"
#     """Verification solution type"""
#     math_problem: MathProblem = MathProblem.DIRECT
#     """Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)"""
#     axisymmetric: bool = False
#     """Axisymmetric simulation for 2D problems"""
#     gravity_force: bool = False
#     """Enable gravity force"""
#     restart_sol: bool = False
#     """Restart solution"""
#     wrt_restart_compact: bool = True
#     """Save only minimum required variables for restart"""
#     discard_infiles: bool = False
#     """Discard data stored in solution and geometry files"""
#     system_measurements: SystemMeasurements = SystemMeasurements.SI
#     """System of measurements (SI, US)"""
#     config_list: List[str] = field(default_factory=lambda: ["configA.cfg", "configB.cfg"])
#     """List of config files for multizone setup"""
