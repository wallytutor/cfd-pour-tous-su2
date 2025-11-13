```python
%load_ext autoreload
%autoreload 2
```

```python
import sys
sys.path.insert(0, "../../")

from su2 import (
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
```

```python
solver              = SolverType.EULER
inlet_type          = InletType.TOTAL_CONDITIONS
convective_scheme   = ConvectiveScheme.JST
linear_solver       = LinearSolver.FGMRES
preconditioner      = Preconditioner.ILU
math_problem        = MathProblem.DIRECT
num_method_grad     = NumMethodGrad.GREEN_GAUSS
mg_cycle            = MgCycle.W_CYCLE
time_discretization = TimeDiscretization.EULER_IMPLICIT

inlet_type.validate_solver(solver)
```

```python
str(mg_cycle.name)
```

```python

```
