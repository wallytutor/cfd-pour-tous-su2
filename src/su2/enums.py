# -*- coding: utf-8 -*-
from enum import Enum


class SolverType(Enum):
    """ Solver type options.

    Documented at:
    https://su2code.github.io/docs_v7/Solver-Setup/#defining-the-problem
    """
    EULER              = "EULER"
    NAVIER_STOKES      = "NAVIER_STOKES"
    RANS               = "RANS"
    INC_EULER          = "INC_EULER"
    INC_NAVIER_STOKES  = "INC_NAVIER_STOKES"
    INC_RANS           = "INC_RANS"
    NEMO_EULER         = "NEMO_EULER"
    NEMO_NAVIER_STOKES = "NEMO_NAVIER_STOKES"
    FEM_EULER          = "FEM_EULER"
    FEM_NAVIER_STOKES  = "FEM_NAVIER_STOKES"
    FEM_RANS           = "FEM_RANS"
    FEM_LES            = "FEM_LES"
    HEAT_EQUATION_FVM  = "HEAT_EQUATION_FVM"
    ELASTICITY         = "ELASTICITY"


class MathProblem(Enum):
    """ Mathematical problem types.

    """
    DIRECT             = "DIRECT"
    CONTINUOUS_ADJOINT = "CONTINUOUS_ADJOINT"
    DISCRETE_ADJOINT   = "DISCRETE_ADJOINT"
