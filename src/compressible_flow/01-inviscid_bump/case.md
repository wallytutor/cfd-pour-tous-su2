```python
%load_ext autoreload
%autoreload 2
```

```python
from majordome.su2 import (
    YesNoEnum,
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
    ReferenceValues,
    BoundaryConditions,
    SurfacesIdentification,
    SU2Configuration,
)
```

```python
inlet_type          = InletType.TOTAL_CONDITIONS
convective_scheme   = ConvectiveScheme.JST
linear_solver       = LinearSolver.FGMRES
preconditioner      = Preconditioner.ILU
math_problem        = MathProblem.NONE
num_method_grad     = NumMethodGrad.GREEN_GAUSS
mg_cycle            = MgCycle.W_CYCLE
time_discretization = TimeDiscretization.EULER_IMPLICIT

# inlet_type.validate_solver(SolverType.EULER)
```

```python
conf = SU2Configuration(
    problem = ProblemDefinition(
        solver = SolverType.EULER,
    )
)

conf.compressible_freestream = CompressibleFreeStreamDefinition(
    mach             = 0.5,
    angle_of_attack  = 0.0,
    sideslip_angle   = 0.0,
    pressure         = 101325.0,
    temperature      = 288.0
)

conf.reference_values = ReferenceValues(
    ref_origin_moment_x = 0.25,
    ref_origin_moment_y = 0.00,
    ref_origin_moment_z = 0.00,
    ref_length          = 1.0,
    ref_area            = 1.0
)

conf.boundary_conditions = BoundaryConditions(
    marker_euler  = ["upper_wall", "lower_wall"],
    marker_inlet  = ["inlet", 288.6, 102010.0, 1.0, 0.0, 0.0],
    marker_outlet = ["outlet", 101300.0],
    inlet_type    = InletType.TOTAL_CONDITIONS,
)

conf.surfaces_identification = SurfacesIdentification(
    marker_plotting   = ["lower_wall"],
    marker_monitoring = ["upper_wall", "lower_wall"],
)

print(conf.to_cfg())
```
