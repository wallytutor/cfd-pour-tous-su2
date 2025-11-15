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
problem_definition = ProblemDefinition(
    solver = SolverType.EULER,
)
print(problem_definition.to_cfg())
```

```python

```
