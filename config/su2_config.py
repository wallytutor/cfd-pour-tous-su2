"""
SU2 Configuration Data Classes

This module provides Python dataclasses that represent the SU2 configuration file structure.
Based on SU2 version 8.3.0 "Harrier".

Author: Auto-generated from config_template.cfg
"""

from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Union
from enum import Enum


# ============================================================================
# Enumerations
# ============================================================================

class Solver(Enum):
    """Solver type options"""
    EULER = "EULER"
    NAVIER_STOKES = "NAVIER_STOKES"
    RANS = "RANS"
    INC_EULER = "INC_EULER"
    INC_NAVIER_STOKES = "INC_NAVIER_STOKES"
    INC_RANS = "INC_RANS"
    NEMO_EULER = "NEMO_EULER"
    NEMO_NAVIER_STOKES = "NEMO_NAVIER_STOKES"
    FEM_EULER = "FEM_EULER"
    FEM_NAVIER_STOKES = "FEM_NAVIER_STOKES"
    FEM_RANS = "FEM_RANS"
    FEM_LES = "FEM_LES"
    HEAT_EQUATION_FVM = "HEAT_EQUATION_FVM"
    ELASTICITY = "ELASTICITY"


class TurbulenceModel(Enum):
    """Turbulence model options"""
    NONE = "NONE"
    SA = "SA"
    SST = "SST"


class MathProblem(Enum):
    """Mathematical problem type"""
    DIRECT = "DIRECT"
    CONTINUOUS_ADJOINT = "CONTINUOUS_ADJOINT"
    DISCRETE_ADJOINT = "DISCRETE_ADJOINT"


class SystemMeasurements(Enum):
    """System of measurements"""
    SI = "SI"  # International System
    US = "US"  # United States customary units


class TimeMarching(Enum):
    """Time marching scheme"""
    NO = "NO"
    TIME_STEPPING = "TIME_STEPPING"
    DUAL_TIME_STEPPING_1ST_ORDER = "DUAL_TIME_STEPPING-1ST_ORDER"
    DUAL_TIME_STEPPING_2ND_ORDER = "DUAL_TIME_STEPPING-2ND_ORDER"
    HARMONIC_BALANCE = "HARMONIC_BALANCE"


class FluidModel(Enum):
    """Fluid model options"""
    STANDARD_AIR = "STANDARD_AIR"
    IDEAL_GAS = "IDEAL_GAS"
    VW_GAS = "VW_GAS"
    PR_GAS = "PR_GAS"
    CONSTANT_DENSITY = "CONSTANT_DENSITY"
    INC_IDEAL_GAS = "INC_IDEAL_GAS"
    INC_IDEAL_GAS_POLY = "INC_IDEAL_GAS_POLY"
    MUTATIONPP = "MUTATIONPP"
    SU2_NONEQ = "SU2_NONEQ"
    FLUID_MIXTURE = "FLUID_MIXTURE"
    COOLPROP = "COOLPROP"
    FLUID_FLAMELET = "FLUID_FLAMELET"
    DATADRIVEN_FLUID = "DATADRIVEN_FLUID"


class ViscosityModel(Enum):
    """Viscosity model options"""
    SUTHERLAND = "SUTHERLAND"
    CONSTANT_VISCOSITY = "CONSTANT_VISCOSITY"
    POLYNOMIAL_VISCOSITY = "POLYNOMIAL_VISCOSITY"
    FLAMELET = "FLAMELET"


class ConductivityModel(Enum):
    """Thermal conductivity model options"""
    CONSTANT_CONDUCTIVITY = "CONSTANT_CONDUCTIVITY"
    CONSTANT_PRANDTL = "CONSTANT_PRANDTL"
    POLYNOMIAL_CONDUCTIVITY = "POLYNOMIAL_CONDUCTIVITY"
    FLAMELET = "FLAMELET"


class ConvectiveScheme(Enum):
    """Convective numerical scheme"""
    JST = "JST"
    JST_KE = "JST_KE"
    JST_MAT = "JST_MAT"
    LAX_FRIEDRICH = "LAX-FRIEDRICH"
    ROE = "ROE"
    AUSM = "AUSM"
    AUSMPLUSUP = "AUSMPLUSUP"
    AUSMPLUSUP2 = "AUSMPLUSUP2"
    AUSMPLUSM = "AUSMPLUSM"
    HLLC = "HLLC"
    TURKEL_PREC = "TURKEL_PREC"
    SW = "SW"
    MSW = "MSW"
    FDS = "FDS"
    SLAU = "SLAU"
    SLAU2 = "SLAU2"
    L2ROE = "L2ROE"
    LMROE = "LMROE"


class SlopeLimiter(Enum):
    """Slope limiter options"""
    NONE = "NONE"
    VENKATAKRISHNAN = "VENKATAKRISHNAN"
    VENKATAKRISHNAN_WANG = "VENKATAKRISHNAN_WANG"
    BARTH_JESPERSEN = "BARTH_JESPERSEN"
    VAN_ALBADA_EDGE = "VAN_ALBADA_EDGE"
    NISHIKAWA_R3 = "NISHIKAWA_R3"
    NISHIKAWA_R4 = "NISHIKAWA_R4"
    NISHIKAWA_R5 = "NISHIKAWA_R5"


class TimeDiscretization(Enum):
    """Time discretization scheme"""
    RUNGE_KUTTA_EXPLICIT = "RUNGE-KUTTA_EXPLICIT"
    EULER_IMPLICIT = "EULER_IMPLICIT"
    EULER_EXPLICIT = "EULER_EXPLICIT"


class LinearSolver(Enum):
    """Linear solver options"""
    BCGSTAB = "BCGSTAB"
    FGMRES = "FGMRES"
    RESTARTED_FGMRES = "RESTARTED_FGMRES"
    CONJUGATE_GRADIENT = "CONJUGATE_GRADIENT"
    SMOOTHER = "SMOOTHER"


class Preconditioner(Enum):
    """Preconditioner options"""
    ILU = "ILU"
    LU_SGS = "LU_SGS"
    LINELET = "LINELET"
    JACOBI = "JACOBI"


class OutputFormat(Enum):
    """Output file format"""
    TECPLOT_ASCII = "TECPLOT_ASCII"
    TECPLOT = "TECPLOT"
    SURFACE_TECPLOT_ASCII = "SURFACE_TECPLOT_ASCII"
    SURFACE_TECPLOT = "SURFACE_TECPLOT"
    CSV = "CSV"
    SURFACE_CSV = "SURFACE_CSV"
    PARAVIEW_ASCII = "PARAVIEW_ASCII"
    PARAVIEW_LEGACY = "PARAVIEW_LEGACY"
    SURFACE_PARAVIEW_ASCII = "SURFACE_PARAVIEW_ASCII"
    SURFACE_PARAVIEW_LEGACY = "SURFACE_PARAVIEW_LEGACY"
    PARAVIEW = "PARAVIEW"
    SURFACE_PARAVIEW = "SURFACE_PARAVIEW"
    RESTART_ASCII = "RESTART_ASCII"
    RESTART = "RESTART"
    CGNS = "CGNS"
    SURFACE_CGNS = "SURFACE_CGNS"
    STL_ASCII = "STL_ASCII"
    STL_BINARY = "STL_BINARY"


# ============================================================================
# Problem Definition
# ============================================================================

@dataclass
class ProblemDefinition:
    """Direct, adjoint, and linearized problem definition"""
    
    solver: Solver = Solver.EULER
    """Solver type"""
    
    kind_turb_model: TurbulenceModel = TurbulenceModel.NONE
    """Turbulence model (NONE, SA, SST)"""
    
    sst_options: str = "NONE"
    """SST model versions/corrections"""
    
    sa_options: str = "NONE"
    """SA model versions/corrections"""
    
    kind_trans_model: str = "NONE"
    """Transition model (NONE, LM)"""
    
    hroughness: float = 1.0e-6
    """RMS roughness for transition model"""
    
    lm_options: str = "NONE"
    """LM model versions/correlations"""
    
    kind_sgs_model: str = "NONE"
    """Subgrid scale model"""
    
    kind_verification_solution: str = "NO_VERIFICATION_SOLUTION"
    """Verification solution type"""
    
    math_problem: MathProblem = MathProblem.DIRECT
    """Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)"""
    
    axisymmetric: bool = False
    """Axisymmetric simulation for 2D problems"""
    
    gravity_force: bool = False
    """Enable gravity force"""
    
    restart_sol: bool = False
    """Restart solution"""
    
    wrt_restart_compact: bool = True
    """Save only minimum required variables for restart"""
    
    discard_infiles: bool = False
    """Discard data stored in solution and geometry files"""
    
    system_measurements: SystemMeasurements = SystemMeasurements.SI
    """System of measurements (SI, US)"""
    
    config_list: List[str] = field(default_factory=lambda: ["configA.cfg", "configB.cfg"])
    """List of config files for multizone setup"""


@dataclass
class SolverControl:
    """Solver control parameters"""
    
    iter: int = 1
    """Number of iterations for single-zone problems"""
    
    inner_iter: int = 9999
    """Maximum number of inner iterations"""
    
    outer_iter: int = 1
    """Maximum number of outer iterations (multizone)"""
    
    time_iter: int = 1
    """Maximum number of time iterations"""
    
    conv_field: str = "DRAG"
    """Convergence field"""
    
    conv_residual_minval: float = -8
    """Min value of the residual (log10)"""
    
    conv_startiter: int = 10
    """Start convergence criteria at iteration number"""
    
    conv_cauchy_elems: int = 100
    """Number of elements to apply the criteria"""
    
    conv_cauchy_eps: float = 1E-10
    """Epsilon to control series convergence"""
    
    restart_iter: int = 0
    """Iteration number to begin unsteady restarts"""
    
    window_cauchy_crit: bool = True
    """Time convergence monitoring"""
    
    conv_window_field: List[str] = field(default_factory=lambda: ["TAVG_DRAG", "TAVG_LIFT"])
    """List of time convergence fields"""
    
    conv_window_startiter: int = 0
    """Window convergence start iteration offset"""
    
    conv_window_cauchy_eps: float = 1E-3
    """Window convergence epsilon"""
    
    conv_window_cauchy_elems: int = 10
    """Number of elements for window convergence"""
    
    max_update_flow: float = 0.2
    """Maximum ratio for updating density and energy variables"""
    
    max_update_sa: float = 0.99
    """Maximum ratio for updating nu_tilde in SA model"""
    
    max_update_sst: float = 1.0
    """Maximum ratio for updating TKE and Omega in SST model"""


@dataclass
class TimeDependent:
    """Time-dependent simulation parameters"""
    
    time_domain: bool = False
    """Enable time domain simulation"""
    
    time_marching: TimeMarching = TimeMarching.NO
    """Time marching scheme"""
    
    time_step: float = 0.0
    """Time step for dual time stepping (s)"""
    
    max_time: float = 50.0
    """Total physical time for dual time stepping (s)"""
    
    unst_cfl_number: float = 0.0
    """Unsteady Courant-Friedrichs-Lewy number"""
    
    window_start_iter: int = 500
    """Time iteration to start windowed time average"""
    
    window_function: str = "SQUARE"
    """Window function (SQUARE, HANN, HANN_SQUARE, BUMP)"""
    
    unst_adjoint_iter: int = 0
    """Starting direct solver iteration for unsteady adjoint"""


@dataclass
class DESParameters:
    """DES (Detached Eddy Simulation) parameters"""
    
    hybrid_ransles: str = "SA_DDES"
    """Hybrid RANS/LES model"""
    
    des_const: float = 0.65
    """DES constant"""


# ============================================================================
# Flow Conditions
# ============================================================================

@dataclass
class CompressibleFreestream:
    """Compressible free-stream definition"""
    
    mach_number: float = 0.8
    """Mach number (non-dimensional)"""
    
    aoa: float = 1.25
    """Angle of attack (degrees)"""
    
    sideslip_angle: float = 0.0
    """Side-slip angle (degrees)"""
    
    init_option: str = "REYNOLDS"
    """Initialization option (REYNOLDS, TD_CONDITIONS)"""
    
    freestream_option: str = "TEMPERATURE_FS"
    """Free-stream option (TEMPERATURE_FS, DENSITY_FS)"""
    
    freestream_pressure: float = 101325.0
    """Free-stream pressure (Pa)"""
    
    freestream_temperature: float = 288.15
    """Free-stream temperature (K)"""
    
    freestream_temperature_ve: float = 288.15
    """Free-stream vibrational temperature (K)"""
    
    reynolds_number: float = 6.5E6
    """Reynolds number (non-dimensional)"""
    
    reynolds_length: float = 1.0
    """Reynolds length (m)"""
    
    freestream_density: float = 1.2886
    """Free-stream density (kg/m³)"""
    
    freestream_velocity: Tuple[float, float, float] = (1.0, 0.0, 0.0)
    """Free-stream velocity vector (m/s)"""
    
    freestream_viscosity: float = 1.853E-5
    """Free-stream viscosity (N·s/m²)"""
    
    freestream_turbulenceintensity: float = 0.05
    """Free-stream turbulence intensity"""
    
    freestream_intermittency: float = 1.0
    """Free-stream intermittency value"""
    
    turb_fixed_values: bool = False
    """Fix turbulence quantities to far-field values"""
    
    turb_fixed_values_domain: float = -1.0
    """Half-space shift for fixed turbulence values"""
    
    freestream_turb2lamviscratio: float = 10.0
    """Ratio between turbulent and laminar viscosity"""
    
    ref_dimensionalization: str = "DIMENSIONAL"
    """Compressible flow non-dimensionalization"""


@dataclass
class IncompressibleFlow:
    """Incompressible flow condition definition"""
    
    inc_density_model: str = "CONSTANT"
    """Density model (CONSTANT, BOUSSINESQ, VARIABLE, FLAMELET)"""
    
    inc_energy_equation: bool = False
    """Solve energy equation"""
    
    inc_density_init: float = 1.2886
    """Initial density (kg/m³)"""
    
    inc_velocity_init: Tuple[float, float, float] = (1.0, 0.0, 0.0)
    """Initial velocity (m/s)"""
    
    inc_temperature_init: float = 288.15
    """Initial temperature (K)"""
    
    inc_nondim: str = "INITIAL_VALUES"
    """Non-dimensionalization scheme"""
    
    inc_density_ref: float = 1.0
    """Reference density (kg/m³)"""
    
    inc_velocity_ref: float = 1.0
    """Reference velocity (m/s)"""
    
    inc_temperature_ref: float = 1.0
    """Reference temperature (K)"""
    
    inc_inlet_type: str = "VELOCITY_INLET"
    """Inlet type"""
    
    inc_inlet_damping: float = 0.1
    """Damping coefficient for pressure inlets"""
    
    inlet_use_normal: bool = False
    """Impose inlet velocity in normal direction"""
    
    inc_outlet_type: str = "PRESSURE_OUTLET"
    """Outlet type"""
    
    inc_outlet_damping: float = 0.1
    """Damping coefficient for mass flow outlets"""
    
    bulk_modulus: float = 1.42E5
    """Bulk modulus for computing Mach number"""
    
    beta_factor: float = 4.1
    """Epsilon² multiplier in Beta calculation"""


@dataclass
class SolidZoneHeat:
    """Solid zone heat variables"""
    
    thermal_conductivity_constant: float = 0.0
    """Thermal conductivity for heat equation"""
    
    freestream_temperature: float = 288.15
    """Solids temperature at freestream conditions (K)"""
    
    material_density: float = 2710.0
    """Density in solids (kg/m³)"""


@dataclass
class CLDriver:
    """CL driver definition for fixed lift mode"""
    
    fixed_cl_mode: bool = False
    """Activate fixed lift mode"""
    
    target_cl: float = 0.80
    """Target coefficient of lift"""
    
    dcl_dalpha: float = 0.2
    """Estimation of dCL/dAlpha (per degree)"""
    
    update_aoa_iter_limit: int = 100
    """Max iterations between AoA updates"""
    
    update_ih: int = 5
    """Number of times Alpha is updated"""
    
    iter_dcl_dalpha: int = 500
    """Iterations to evaluate dCL_dAlpha by finite differences"""
    
    eval_dof_dcx: bool = False
    """Evaluate dOF_dCL or dOF_dCMy during runtime"""
    
    netthrust_dbcthrust: float = 1.0
    """Damping factor for thrust BC (actuator disk)"""
    
    dcd_dcl_value: float = 0.0
    """Parameter for complex objective function"""
    
    dcmx_dcl_value: float = 0.0
    """Parameter for complex objective function"""
    
    dcmy_dcl_value: float = 0.0
    """Parameter for complex objective function"""
    
    dcmz_dcl_value: float = 0.0
    """Parameter for complex objective function"""


@dataclass
class ReferenceValues:
    """Reference value definition"""
    
    ref_origin_moment_x: float = 0.25
    """Reference origin X for moment computation (m)"""
    
    ref_origin_moment_y: float = 0.00
    """Reference origin Y for moment computation (m)"""
    
    ref_origin_moment_z: float = 0.00
    """Reference origin Z for moment computation (m)"""
    
    ref_length: float = 1.0
    """Reference length (m)"""
    
    ref_velocity: float = 1.0
    """Reference velocity (m/s)"""
    
    ref_viscosity: float = 1.0
    """Reference viscosity (N·s/m²)"""
    
    ref_area: float = 1.0
    """Reference area (m²)"""
    
    semi_span: float = 0.0
    """Aircraft semi-span (m)"""


# ============================================================================
# Fluid Properties
# ============================================================================

@dataclass
class FluidProperties:
    """Fluid model and properties"""
    
    fluid_model: FluidModel = FluidModel.STANDARD_AIR
    """Fluid model"""
    
    fluid_name: str = "nitrogen"
    """Fluid name for CoolProp library"""
    
    gamma_value: float = 1.4
    """Ratio of specific heats"""
    
    gas_constant: float = 287.058
    """Specific gas constant (J/kg·K)"""
    
    critical_temperature: float = 131.00
    """Critical temperature (K)"""
    
    critical_pressure: float = 3588550.0
    """Critical pressure (Pa)"""
    
    critical_density: float = 263.0
    """Critical density (kg/m³)"""
    
    acentric_factor: float = 0.035
    """Acentric factor"""
    
    thermodynamic_pressure: float = 101325.0
    """Operating pressure (Pa)"""
    
    specific_heat_cp: float = 1004.703
    """Specific heat at constant pressure (J/kg·K)"""
    
    thermal_expansion_coeff: float = 0.00347
    """Thermal expansion coefficient (1/K)"""
    
    molecular_weight: List[float] = field(default_factory=lambda: [28.96, 16.043])
    """Molecular weights of species (g/mol)"""
    
    cp_polycoeffs: Tuple[float, ...] = (0.0, 0.0, 0.0, 0.0, 0.0)
    """Temperature polynomial coefficients for Cp"""
    
    gas_model: str = "AIR-5"
    """Gas model for mixtures"""
    
    gas_composition: Tuple[float, ...] = (0.77, 0.23, 0.0, 0.0, 0.0)
    """Initial gas composition in mass fractions"""
    
    frozen_mixture: bool = False
    """Freeze chemical reactions"""
    
    interpolation_method: str = "MLP"
    """Interpolation method (MLP, LUT)"""
    
    filenames_interpolator: List[str] = field(default_factory=lambda: ["MLP_1.mlp", "MLP_2.mlp", "MLP_3.mlp"])
    """Input files for interpolator"""
    
    datadriven_newton_relaxation: float = 0.8
    """Relaxation factor for Newton solvers"""
    
    use_pinn: bool = True
    """Use physics-informed neural network"""
    
    datadriven_initial_density: float = -1
    """Initial guess for fluid density"""
    
    datadriven_initial_energy: float = -1
    """Initial guess for fluid energy"""
    
    ionization: bool = False
    """Specify if there is ionization"""
    
    vt_residual_limiting: bool = False
    """VT transfer residual limiting"""
    
    inlet_temperature_ve: float = 288.15
    """NEMO inlet vibrational temperature (K)"""
    
    inlet_gas_composition: Tuple[float, ...] = (0.77, 0.23, 0.0, 0.0, 0.0)
    """NEMO inlet gas composition"""


@dataclass
class ViscosityProperties:
    """Viscosity model definition"""
    
    viscosity_model: ViscosityModel = ViscosityModel.SUTHERLAND
    """Viscosity model"""
    
    mu_constant: float = 1.716E-5
    """Constant molecular viscosity (N·s/m²)"""
    
    mu_ref: float = 1.716E-5
    """Sutherland viscosity reference (N·s/m²)"""
    
    mu_t_ref: float = 273.15
    """Sutherland temperature reference (K)"""
    
    sutherland_constant: float = 110.4
    """Sutherland constant (K)"""
    
    mu_polycoeffs: Tuple[float, ...] = (0.0, 0.0, 0.0, 0.0, 0.0)
    """Temperature polynomial coefficients for viscosity"""


@dataclass
class ConductivityProperties:
    """Thermal conductivity model definition"""
    
    conductivity_model: ConductivityModel = ConductivityModel.CONSTANT_PRANDTL
    """Laminar conductivity model"""
    
    thermal_conductivity_constant: float = 0.0257
    """Molecular thermal conductivity (W/m·K)"""
    
    prandtl_lam: float = 0.72
    """Laminar Prandtl number"""
    
    kt_polycoeffs: Tuple[float, ...] = (0.0, 0.0, 0.0, 0.0, 0.0)
    """Temperature polynomial coefficients for conductivity"""
    
    turbulent_conductivity_model: str = "CONSTANT_PRANDTL_TURB"
    """Turbulent thermal conductivity model"""
    
    prandtl_turb: float = 0.90
    """Turbulent Prandtl number"""


# ============================================================================
# Dynamic Mesh and Motion
# ============================================================================

@dataclass
class DynamicMesh:
    """Dynamic mesh definition"""
    
    grid_movement: str = "NONE"
    """Type of dynamic mesh"""
    
    mach_motion: float = 0.8
    """Motion Mach number"""
    
    motion_origin: Tuple[float, float, float] = (0.25, 0.0, 0.0)
    """Coordinates of motion origin"""
    
    rotation_rate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Angular velocity vector (rad/s)"""
    
    pitching_omega: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Pitching angular frequency (rad/s)"""
    
    pitching_ampl: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Pitching amplitude (degrees)"""
    
    pitching_phase: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Pitching phase offset (degrees)"""
    
    translation_rate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Translational velocity (m/s)"""
    
    plunging_omega: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Plunging angular frequency (rad/s)"""
    
    plunging_ampl: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Plunging amplitude (m)"""
    
    surface_movement: str = "NONE"
    """Type of dynamic surface movement"""
    
    marker_moving: List[str] = field(default_factory=lambda: ["NONE"])
    """Moving wall boundary markers"""
    
    surface_motion_origin: float = 0.25
    """Surface motion origin coordinate"""
    
    surface_rotation_rate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface angular velocity (rad/s)"""
    
    surface_pitching_omega: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface pitching frequency (rad/s)"""
    
    surface_pitching_ampl: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface pitching amplitude (degrees)"""
    
    surface_pitching_phase: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface pitching phase (degrees)"""
    
    surface_translation_rate: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface translational velocity (m/s)"""
    
    surface_plunging_omega: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface plunging frequency (rad/s)"""
    
    surface_plunging_ampl: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Surface plunging amplitude (m)"""
    
    move_motion_origin: bool = False
    """Move motion origin for marker moving"""


@dataclass
class BuffetSensor:
    """Buffet sensor definition"""
    
    buffet_monitoring: bool = False
    """Evaluate buffet sensor"""
    
    buffet_k: float = 10.0
    """Sharpness coefficient for Heaviside function"""
    
    buffet_lambda: float = 0.0
    """Offset parameter for Heaviside function"""


@dataclass
class Aeroelastic:
    """Aeroelastic simulation (Typical Section Model)"""
    
    flutter_speed_index: float = 0.6
    """Flutter speed index"""
    
    plunge_natural_frequency: float = 100
    """Natural frequency in plunging direction (rad/s)"""
    
    pitch_natural_frequency: float = 100
    """Natural frequency in pitching direction (rad/s)"""
    
    airfoil_mass_ratio: float = 60
    """Airfoil mass ratio"""
    
    cg_location: float = 1.8
    """CG distance behind elastic axis (semichords)"""
    
    radius_gyration_squared: float = 3.48
    """Radius of gyration squared (semichords²)"""
    
    aeroelastic_iter: int = 3
    """Solve aeroelastic equations every N iterations"""


@dataclass
class GustSimulation:
    """Gust simulation parameters"""
    
    wind_gust: bool = False
    """Apply a wind gust"""
    
    gust_type: str = "NONE"
    """Type of gust"""
    
    gust_dir: str = "Y_DIR"
    """Direction of the gust"""
    
    gust_wavelength: float = 10.0
    """Gust wavelength (m)"""
    
    gust_periods: float = 1.0
    """Number of gust periods"""
    
    gust_ampl: float = 10.0
    """Gust amplitude (m/s)"""
    
    gust_begin_time: float = 0.0
    """Time to begin the gust (s)"""
    
    gust_begin_loc: float = 0.0
    """Location where gust begins (m)"""


# ============================================================================
# Boundary Conditions
# ============================================================================

@dataclass
class BoundaryConditions:
    """Boundary condition definitions"""
    
    marker_euler: List[str] = field(default_factory=lambda: ["airfoil"])
    """Euler wall boundary markers"""
    
    marker_heatflux: List[str] = field(default_factory=lambda: ["NONE"])
    """Constant heat flux wall markers"""
    
    marker_heattransfer: List[str] = field(default_factory=lambda: ["NONE"])
    """Heat transfer/convection wall markers"""
    
    marker_isothermal: List[str] = field(default_factory=lambda: ["NONE"])
    """Isothermal wall markers"""
    
    marker_far: List[str] = field(default_factory=lambda: ["farfield"])
    """Far-field boundary markers"""
    
    marker_sym: List[str] = field(default_factory=lambda: ["NONE"])
    """Symmetry boundary markers"""
    
    marker_internal: List[str] = field(default_factory=lambda: ["NONE"])
    """Internal boundary markers"""
    
    marker_nearfield: List[str] = field(default_factory=lambda: ["NONE"])
    """Near-field boundary markers"""
    
    inlet_type: str = "TOTAL_CONDITIONS"
    """Inlet boundary type"""
    
    specified_inlet_profile: bool = False
    """Read inlet profile from file"""
    
    inlet_filename: str = "inlet.dat"
    """File specifying inlet profile"""
    
    inlet_interpolation_function: str = "NONE"
    """Type of spanwise interpolation for inlet"""
    
    inlet_interpolation_data_type: str = "VRVTHETA"
    """Type of radial spanwise interpolation"""
    
    print_inlet_interpolated_data: bool = False
    """Write interpolated inlet data to file"""
    
    marker_inlet: List[str] = field(default_factory=lambda: ["NONE"])
    """Inlet boundary markers with conditions"""
    
    marker_outlet: List[str] = field(default_factory=lambda: ["NONE"])
    """Outlet boundary markers"""
    
    actdisk_double_surface: bool = False
    """Actuator disk double surface"""
    
    actdisk_type: str = "VARIABLES_JUMP"
    """Actuator disk boundary type"""
    
    marker_actdisk: List[str] = field(default_factory=lambda: ["NONE"])
    """Actuator disk boundary markers"""
    
    marker_actdisk_bem_cg: List[str] = field(default_factory=lambda: ["NONE"])
    """Blade element CG markers"""
    
    marker_actdisk_bem_axis: List[str] = field(default_factory=lambda: ["NONE"])
    """Blade element axis markers"""
    
    actdisk_filename: str = "actuatordisk.dat"
    """Actuator disk data input file"""
    
    bem_prop_filename: str = "prop_geom_alfclcd_data.txt"
    """Propeller blade element data file"""
    
    bem_prop_blade_angle: float = 25.0
    """Propeller blade angle at 0.75*radius (degrees)"""
    
    bem_freq: int = 40
    """BEM calculation frequency"""
    
    marker_supersonic_inlet: List[str] = field(default_factory=lambda: ["NONE"])
    """Supersonic inlet markers"""
    
    marker_supersonic_outlet: List[str] = field(default_factory=lambda: ["NONE"])
    """Supersonic outlet markers"""
    
    marker_periodic: List[str] = field(default_factory=lambda: ["NONE"])
    """Periodic boundary markers"""
    
    engine_inflow_type: str = "FAN_FACE_MACH"
    """Engine inflow boundary type"""
    
    marker_engine_inflow: List[str] = field(default_factory=lambda: ["NONE"])
    """Engine inflow markers"""
    
    marker_engine_exhaust: List[str] = field(default_factory=lambda: ["NONE"])
    """Engine exhaust markers"""
    
    marker_normal_displ: List[str] = field(default_factory=lambda: ["NONE"])
    """Displacement boundary markers"""
    
    marker_pressure: List[str] = field(default_factory=lambda: ["NONE"])
    """Pressure boundary markers"""
    
    marker_riemann: List[str] = field(default_factory=lambda: ["NONE"])
    """Riemann boundary markers"""
    
    marker_shroud: List[str] = field(default_factory=lambda: ["NONE"])
    """Shroud boundary markers"""
    
    marker_zone_interface: List[str] = field(default_factory=lambda: ["NONE"])
    """Zone interface markers"""
    
    marker_cht_interface: List[str] = field(default_factory=lambda: ["NONE"])
    """CHT interface markers"""
    
    marker_fluid_interface: List[str] = field(default_factory=lambda: ["NONE"])
    """Fluid interface markers"""
    
    marker_fluid_load: List[str] = field(default_factory=lambda: ["NONE"])
    """Fluid load markers"""
    
    kind_interpolation: str = "NEAREST_NEIGHBOR"
    """Interface interpolation kind"""
    
    conservative_interpolation: bool = True
    """Use conservative interpolation"""
    
    kind_radial_basis_function: str = "WENDLAND_C2"
    """Type of radial basis function"""
    
    radial_basis_function_parameter: float = 0.015
    """Radius for radial basis function"""
    
    radial_basis_function_polynomial_term: bool = True
    """Use polynomial term in RBF interpolation"""
    
    radial_basis_function_prune_tolerance: float = 0
    """Tolerance to prune RBF coefficients"""


# ============================================================================
# Numerical Methods
# ============================================================================

@dataclass
class NumericalMethods:
    """Common parameters defining numerical methods"""
    
    num_method_grad: str = "GREEN_GAUSS"
    """Numerical method for spatial gradients"""
    
    num_method_grad_recon: str = "LEAST_SQUARES"
    """Gradients for MUSCL reconstruction"""
    
    cfl_number: float = 15.0
    """CFL number"""
    
    cfl_adapt: bool = False
    """Adaptive CFL number"""
    
    cfl_adapt_param: Tuple[float, ...] = (0.1, 2.0, 10.0, 1e10, 0.001, 0)
    """Adaptive CFL parameters"""
    
    max_delta_time: float = 1E6
    """Maximum delta time in local time stepping"""
    
    ext_iter_offset: int = 0
    """External iteration offset due to restart"""
    
    rk_alpha_coeff: Tuple[float, ...] = (0.66667, 0.66667, 1.000000)
    """Runge-Kutta alpha coefficients"""
    
    objective_function: str = "DRAG"
    """Objective function in gradient evaluation"""
    
    objective_weight: float = 1.0
    """Weighting value for objective function"""
    
    custom_objfunc: str = "'DRAG + 10 * pow(fmax(0.4-LIFT, 0), 2)'"
    """Custom objective function expression"""


@dataclass
class SlopeLimiters:
    """Slope limiter and dissipation sensor definition"""
    
    muscl_flow: bool = True
    """MUSCL in flow equations"""
    
    slope_limiter_flow: SlopeLimiter = SlopeLimiter.VENKATAKRISHNAN
    """Slope limiter for flow"""
    
    muscl_turb: bool = False
    """MUSCL in turbulence equations"""
    
    slope_limiter_turb: SlopeLimiter = SlopeLimiter.VENKATAKRISHNAN
    """Slope limiter for turbulence"""
    
    muscl_adjflow: bool = True
    """MUSCL in adjoint flow equations"""
    
    slope_limiter_adjflow: SlopeLimiter = SlopeLimiter.VENKATAKRISHNAN
    """Slope limiter for adjoint flow"""
    
    muscl_adjturb: bool = False
    """MUSCL in adjoint turbulence equations"""
    
    slope_limiter_adjturb: SlopeLimiter = SlopeLimiter.VENKATAKRISHNAN
    """Slope limiter for adjoint turbulence"""
    
    venkat_limiter_coeff: float = 0.05
    """Venkatakrishnan limiter coefficient"""
    
    ref_sharp_edges: float = 3.0
    """Reference coefficient for sharp edges"""
    
    adj_sharp_limiter_coeff: float = 3.0
    """Adjoint sharp edges limiter coefficient"""
    
    sens_remove_sharp: bool = False
    """Remove sharp edges from sensitivity"""
    
    sens_smoothing: str = "NONE"
    """Sensitivity smoothing"""
    
    limiter_iter: int = 999999
    """Freeze limiter after N iterations"""
    
    lax_sensor_coeff: float = 0.15
    """Lax-Friedrichs sensor coefficient"""
    
    jst_sensor_coeff: Tuple[float, float] = (0.5, 0.02)
    """JST sensor coefficients (2nd, 4th order)"""
    
    adj_lax_sensor_coeff: float = 0.15
    """Adjoint Lax-Friedrichs coefficient"""
    
    adj_jst_sensor_coeff: Tuple[float, float] = (0.5, 0.02)
    """Adjoint JST sensor coefficients"""


@dataclass
class LinearSolverDef:
    """Linear solver definition"""
    
    linear_solver: LinearSolver = LinearSolver.FGMRES
    """Linear solver or smoother"""
    
    discadj_lin_solver: LinearSolver = LinearSolver.FGMRES
    """Linear solver for discrete adjoint"""
    
    adjturb_lin_solver: LinearSolver = LinearSolver.FGMRES
    """Linear solver for turbulent adjoint"""
    
    enable_cuda: bool = False
    """Use CUDA GPU acceleration"""
    
    adjturb_lin_prec: Preconditioner = Preconditioner.ILU
    """Preconditioner for turbulent adjoint"""
    
    adjturb_lin_iter: int = 10
    """Max iterations for turbulent adjoint solver"""
    
    linear_solver_prec: Preconditioner = Preconditioner.ILU
    """Preconditioner of Krylov solver"""
    
    discadj_lin_prec: Preconditioner = Preconditioner.ILU
    """Preconditioner for discrete adjoint"""
    
    linear_solver_ilu_fill_in: int = 0
    """ILU preconditioner fill-in level"""
    
    linear_solver_error: float = 1E-6
    """Minimum error of linear solver"""
    
    linear_solver_iter: int = 5
    """Max iterations of linear solver"""
    
    linear_solver_restart_frequency: int = 10
    """Restart frequency for RESTARTED_FGMRES"""
    
    linear_solver_smoother_relaxation: float = 1.0
    """Relaxation factor for smoother-type solvers"""


@dataclass
class Multigrid:
    """Multigrid parameters"""
    
    mglevel: int = 0
    """Multi-grid levels (0 = no multi-grid)"""
    
    mgcycle: str = "V_CYCLE"
    """Multi-grid cycle"""
    
    mg_pre_smooth: Tuple[int, ...] = (1, 2, 3, 3)
    """Pre-smoothing level"""
    
    mg_post_smooth: Tuple[int, ...] = (0, 0, 0, 0)
    """Post-smoothing level"""
    
    mg_correction_smooth: Tuple[int, ...] = (0, 0, 0, 0)
    """Jacobi implicit smoothing"""
    
    mg_damp_restriction: float = 0.75
    """Damping factor for residual restriction"""
    
    mg_damp_prolongation: float = 0.75
    """Damping factor for correction prolongation"""
    
    smooth_geometry: int = 0
    """Implicitly smooth nodal coordinates"""


@dataclass
class FlowNumericalMethod:
    """Flow numerical method definition"""
    
    conv_num_method_flow: ConvectiveScheme = ConvectiveScheme.ROE
    """Convective numerical method"""
    
    roe_low_dissipation: str = "FD"
    """Roe low dissipation function"""
    
    roe_kappa: float = 0.5
    """Roe dissipation coefficient"""
    
    min_roe_turkel_prec: float = 0.01
    """Minimum beta for Roe-Turkel preconditioner"""
    
    max_roe_turkel_prec: float = 0.2
    """Maximum beta for Roe-Turkel preconditioner"""
    
    low_mach_corr: bool = False
    """Post-reconstruction correction for low Mach"""
    
    low_mach_prec: bool = False
    """Roe-Turkel preconditioning for low Mach"""
    
    use_accurate_flux_jacobians: bool = False
    """Use numerically computed Jacobians"""
    
    use_vectorization: bool = True
    """Use vectorized version of scheme"""
    
    entropy_fix_coeff: float = 0.0
    """Entropy fix coefficient"""
    
    central_jacobian_fix_factor: float = 4.0
    """Central schemes Jacobian fix factor"""
    
    central_inc_jacobian_fix_factor: float = 1.0
    """Incompressible central Jacobian fix factor"""
    
    time_discre_flow: TimeDiscretization = TimeDiscretization.EULER_IMPLICIT
    """Time discretization"""
    
    newton_krylov: bool = False
    """Use Newton-Krylov method"""
    
    newton_krylov_iparam: Tuple[int, ...] = (10, 3, 2)
    """Newton-Krylov integer parameters"""
    
    newton_krylov_dparam: Tuple[float, ...] = (1.0, 0.1, -6.0, 1e-5)
    """Newton-Krylov double parameters"""


@dataclass
class TurbulenceNumericalMethod:
    """Turbulent numerical method definition"""
    
    conv_num_method_turb: str = "SCALAR_UPWIND"
    """Convective numerical method"""
    
    time_discre_turb: TimeDiscretization = TimeDiscretization.EULER_IMPLICIT
    """Time discretization"""
    
    cfl_reduction_turb: float = 1.0
    """CFL reduction factor"""
    
    lower_limit_k_factor: float = 1e-15
    """Lower limit constant for k (SST)"""
    
    lower_limit_omega_factor: float = 1e-5
    """Lower limit constant for omega (SST)"""
    
    use_accurate_turb_jacobians: bool = False
    """Use numerically computed exact Jacobians"""


# ============================================================================
# Output and IO
# ============================================================================

@dataclass
class ScreenHistoryVolume:
    """Screen/History/Volume output configuration"""
    
    screen_output: List[str] = field(default_factory=lambda: 
        ["INNER_ITER", "RMS_DENSITY", "RMS_MOMENTUM-X", "RMS_MOMENTUM-Y", "RMS_ENERGY"])
    """Screen output fields"""
    
    history_output: List[str] = field(default_factory=lambda: ["ITER", "RMS_RES"])
    """History output groups"""
    
    custom_outputs: str = "''"
    """User defined functions for output"""
    
    volume_output: List[str] = field(default_factory=lambda: 
        ["COORDINATES", "SOLUTION", "PRIMITIVE"])
    """Volume output fields/groups"""
    
    screen_wrt_freq_inner: int = 1
    """Writing frequency for screen output (inner)"""
    
    screen_wrt_freq_outer: int = 1
    """Writing frequency for screen output (outer)"""
    
    screen_wrt_freq_time: int = 1
    """Writing frequency for screen output (time)"""
    
    history_wrt_freq_inner: int = 1
    """Writing frequency for history output (inner)"""
    
    history_wrt_freq_outer: int = 1
    """Writing frequency for history output (outer)"""
    
    history_wrt_freq_time: int = 1
    """Writing frequency for history output (time)"""
    
    output_wrt_freq: List[int] = field(default_factory=lambda: [10, 250, 42])
    """List of writing frequencies for output files"""
    
    wrt_performance: bool = False
    """Output performance summary to console"""
    
    wrt_ad_statistics: bool = False
    """Output tape statistics (discrete adjoint)"""
    
    wrt_restart_overwrite: bool = True
    """Overwrite or append iteration to restart files"""
    
    wrt_surface_overwrite: bool = True
    """Overwrite or append iteration to surface files"""
    
    wrt_volume_overwrite: bool = True
    """Overwrite or append iteration to volume files"""
    
    wrt_forces_breakdown: bool = False
    """Write forces breakdown"""
    
    comm_level: str = "FULL"
    """MPI communication level"""
    
    visualize_cv: int = -1
    """Node number for CV visualization"""
    
    extra_output: bool = False
    """Write extra output (experimental)"""
    
    extra_heat_zone_output: int = -1
    """Write extra heat output for zone"""


@dataclass
class SurfacesIdentification:
    """Surfaces identification"""
    
    marker_plotting: List[str] = field(default_factory=lambda: ["airfoil"])
    """Surface flow solution markers"""
    
    marker_monitoring: List[str] = field(default_factory=lambda: ["airfoil"])
    """Surface where coefficients are evaluated"""
    
    marker_wall_functions: List[str] = field(default_factory=lambda: 
        ["airfoil", "NO_WALL_FUNCTION"])
    """Viscous wall markers for wall functions"""
    
    marker_python_custom: List[str] = field(default_factory=lambda: ["NONE"])
    """Surface with custom thermal BCs"""
    
    marker_designing: List[str] = field(default_factory=lambda: ["airfoil"])
    """Surface for design problem evaluation"""
    
    marker_analyze: List[str] = field(default_factory=lambda: ["airfoil"])
    """Surface for detailed analysis"""
    
    marker_analyze_average: str = "MASSFLUX"
    """Method to compute average in analysis"""


@dataclass
class InputOutputFiles:
    """Input/Output file information"""
    
    mesh_filename: str = "mesh_NACA0012_inv"
    """Mesh input file"""
    
    mesh_format: str = "SU2"
    """Mesh input file format"""
    
    mesh_box_size: Tuple[int, int, int] = (33, 33, 33)
    """Grid points in RECTANGLE/BOX grid"""
    
    mesh_box_length: Tuple[float, float, float] = (1.0, 1.0, 1.0)
    """Length of RECTANGLE/BOX grid"""
    
    mesh_box_offset: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    """Offset of RECTANGLE/BOX grid"""
    
    mesh_out_filename: str = "mesh_out"
    """Mesh output file"""
    
    solution_filename: str = "solution_flow"
    """Restart flow input file"""
    
    solution_adj_filename: str = "solution_adj"
    """Restart adjoint input file"""
    
    tabular_format: str = "CSV"
    """Output tabular file format"""
    
    output_precision: int = 10
    """Set output precision for SU2_DOT and HISTORY"""
    
    multizone_adapt_filename: bool = True
    """Extend filenames by zone number"""
    
    output_files: List[OutputFormat] = field(default_factory=lambda: 
        [OutputFormat.RESTART, OutputFormat.PARAVIEW, OutputFormat.SURFACE_PARAVIEW])
    """Files to output"""
    
    conv_filename: str = "history"
    """Convergence history output file"""
    
    breakdown_filename: str = "forces_breakdown.dat"
    """Forces breakdown output file"""
    
    restart_filename: str = "restart_flow"
    """Restart flow output file"""
    
    restart_adj_filename: str = "restart_adj"
    """Restart adjoint output file"""
    
    volume_filename: str = "flow"
    """Volume flow output file"""
    
    volume_adj_filename: str = "adjoint"
    """Volume adjoint output file"""
    
    value_objfunc_filename: str = "of_eval"
    """Objective function output file"""
    
    grad_objfunc_filename: str = "of_grad"
    """Objective function gradient file"""
    
    surface_filename: str = "surface_flow"
    """Surface flow coefficient file"""
    
    surface_adj_filename: str = "surface_adjoint"
    """Surface adjoint coefficient file"""
    
    surface_sens_filename: str = "surface_sens"
    """Surface sensitivity file (discrete adjoint)"""
    
    volume_sens_filename: str = "volume_sens"
    """Volume sensitivity file"""
    
    read_binary_restart: bool = True
    """Read binary restart files"""
    
    reorient_elements: bool = True
    """Reorient elements based on negative volumes"""


# ============================================================================
# Main Configuration Class
# ============================================================================

@dataclass
class SU2Config:
    """
    Complete SU2 configuration file representation.
    
    This dataclass encompasses all configuration parameters for SU2 version 8.3.0 "Harrier".
    Each section corresponds to a logical grouping in the configuration file.
    
    Example:
        >>> config = SU2Config()
        >>> config.problem.solver = Solver.NAVIER_STOKES
        >>> config.compressible.mach_number = 0.85
        >>> config.compressible.aoa = 2.5
    """
    
    # Problem definition and control
    problem: ProblemDefinition = field(default_factory=ProblemDefinition)
    """Direct, adjoint, and linearized problem definition"""
    
    solver_control: SolverControl = field(default_factory=SolverControl)
    """Solver control parameters"""
    
    time_dependent: TimeDependent = field(default_factory=TimeDependent)
    """Time-dependent simulation parameters"""
    
    des_params: DESParameters = field(default_factory=DESParameters)
    """DES parameters"""
    
    # Flow conditions
    compressible: CompressibleFreestream = field(default_factory=CompressibleFreestream)
    """Compressible free-stream definition"""
    
    incompressible: IncompressibleFlow = field(default_factory=IncompressibleFlow)
    """Incompressible flow condition definition"""
    
    solid_heat: SolidZoneHeat = field(default_factory=SolidZoneHeat)
    """Solid zone heat variables"""
    
    cl_driver: CLDriver = field(default_factory=CLDriver)
    """CL driver definition"""
    
    reference: ReferenceValues = field(default_factory=ReferenceValues)
    """Reference values"""
    
    # Fluid properties
    fluid: FluidProperties = field(default_factory=FluidProperties)
    """Fluid model and properties"""
    
    viscosity: ViscosityProperties = field(default_factory=ViscosityProperties)
    """Viscosity model"""
    
    conductivity: ConductivityProperties = field(default_factory=ConductivityProperties)
    """Thermal conductivity model"""
    
    # Dynamic mesh and motion
    dynamic_mesh: DynamicMesh = field(default_factory=DynamicMesh)
    """Dynamic mesh definition"""
    
    buffet: BuffetSensor = field(default_factory=BuffetSensor)
    """Buffet sensor"""
    
    aeroelastic: Aeroelastic = field(default_factory=Aeroelastic)
    """Aeroelastic simulation"""
    
    gust: GustSimulation = field(default_factory=GustSimulation)
    """Gust simulation"""
    
    # Boundary conditions
    boundaries: BoundaryConditions = field(default_factory=BoundaryConditions)
    """Boundary condition definitions"""
    
    # Numerical methods
    numerical: NumericalMethods = field(default_factory=NumericalMethods)
    """Common numerical method parameters"""
    
    limiters: SlopeLimiters = field(default_factory=SlopeLimiters)
    """Slope limiters and dissipation sensors"""
    
    linear_solver: LinearSolverDef = field(default_factory=LinearSolverDef)
    """Linear solver definition"""
    
    multigrid: Multigrid = field(default_factory=Multigrid)
    """Multigrid parameters"""
    
    flow_method: FlowNumericalMethod = field(default_factory=FlowNumericalMethod)
    """Flow numerical method"""
    
    turb_method: TurbulenceNumericalMethod = field(default_factory=TurbulenceNumericalMethod)
    """Turbulence numerical method"""
    
    # Output and IO
    output: ScreenHistoryVolume = field(default_factory=ScreenHistoryVolume)
    """Screen/History/Volume output"""
    
    surfaces: SurfacesIdentification = field(default_factory=SurfacesIdentification)
    """Surfaces identification"""
    
    io_files: InputOutputFiles = field(default_factory=InputOutputFiles)
    """Input/Output files"""
    
    def to_cfg(self, filepath: str) -> None:
        """
        Write configuration to SU2 .cfg file.
        
        Args:
            filepath: Path to output configuration file
        """
        # TODO: Implement conversion back to .cfg format
        raise NotImplementedError("Conversion to .cfg format not yet implemented")
    
    @classmethod
    def from_cfg(cls, filepath: str) -> 'SU2Config':
        """
        Load configuration from SU2 .cfg file.
        
        Args:
            filepath: Path to configuration file
            
        Returns:
            SU2Config instance populated from file
        """
        # TODO: Implement parser for .cfg format
        raise NotImplementedError("Parsing from .cfg format not yet implemented")
