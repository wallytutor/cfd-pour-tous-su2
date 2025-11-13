# -*- coding: utf-8 -*-

SU2_VERSION = "8.3.0"

__all__ = []

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# enums
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from .enums import (
    SolverType,
    InletType,
    ConvectiveScheme,
    LinearSolver,
    Preconditioner,
    MathProblem,
    NumMethodGrad,
    MgCycle,
    TimeDiscretization
)

__all__ += [
    "SolverType",
    "InletType",
    "ConvectiveScheme",
    "LinearSolver",
    "Preconditioner",
    "MathProblem",
    "NumMethodGrad",
    "MgCycle",
    "TimeDiscretization"
]

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# EOF
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+