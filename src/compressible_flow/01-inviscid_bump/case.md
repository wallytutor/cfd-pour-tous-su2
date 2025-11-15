```python
%load_ext autoreload
%autoreload 2
```

```python
from majordome.su2 import (
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
from majordome.su2 import (
    ProblemDefinition,
    CompressibleFreeStreamDefinition,
)
```

```python
inlet_type          = InletType.TOTAL_CONDITIONS
convective_scheme   = ConvectiveScheme.JST
linear_solver       = LinearSolver.FGMRES
preconditioner      = Preconditioner.ILU
math_problem        = MathProblem.DIRECT
num_method_grad     = NumMethodGrad.GREEN_GAUSS
mg_cycle            = MgCycle.W_CYCLE
time_discretization = TimeDiscretization.EULER_IMPLICIT

# inlet_type.validate_solver(SolverType.EULER)
```

```python
text = ""

defs = ProblemDefinition(
    solver = SolverType.EULER,
)

text += defs.to_cfg()

defs = CompressibleFreeStreamDefinition(
    mach             = 0.5,
    angle_of_attack  = 0.0,
    sideslip_angle   = 0.0,
    pressure         = 101325.0,
    temperature      = 288.0,
)

text += defs.to_cfg()

print(text)
```

```python

```
