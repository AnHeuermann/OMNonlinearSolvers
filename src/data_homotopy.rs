// (C) 2026, Andreas Heuermann

use crate::config::CONFIGURATION;
use crate::data_homotopy::TearingKind::Unknown;
use crate::nls_user_data::NonlinearSystemUserData;

enum TearingKind {
    Unknown,
    Strict,
    Casual,
}

enum Direction {
    Forward,
    Backward,
}

impl Direction {
    fn sign(self) -> f64 {
        match self {
            Direction::Forward => 1.0,
            Direction::Backward => -1.0,
        }
    }
}

/// Memory structure for nonlinear system solver
pub struct DataHomotopy {
    initialized: bool,

    /// Dimension (n = size)
    n: u32,
    /// Dimension (n = size + 1)
    m: u32,

    /// Tolerance for updating solution vector
    xtol_sqrd: f64,
    /// Tolerance for accepting accuracy
    ftol_sqrd: f64,

    error_f_sqrd: f64,

    /// Residual scaling
    res_scaling: Box<[f64]>,
    /// Scaled function values
    fvec_scaled: Box<[f64]>,
    /// Scaled function values
    hvec_scaled: Box<[f64]>,
    // Scaled solution
    dx_scaled: Box<[f64]>,

    /// `min` attribute of variable
    // min_values: &'a [f64],
    /// `max` attribute of variable
    // max_values: &'a [f64],
    /// `nominal` attribute of variable [x.nominal,lambda.nominal] with lambda.nominal=1.0
    // x_scaling_values: &'a [f64],

    /// Used in wrapper_
    f1: Box<[f64]>,
    /// Used in wrapper_
    f2: Box<[f64]>,
    /// Used for steepest descent method
    grad_fx: Box<[f64]>,

    // /* return value, if success info == 1 */
    // int info;
    /// Number of iterations over the whole simulation time
    number_of_iterations: u32,
    /// Number of function evaluations over the whole simulation time
    number_of_function_evaluations: u32,
    /// Maximum allowed number of Newton steps
    max_number_of_iteartions: u32,

    // /* strict tearing set or casual tearing set */
    causal_tearing_set: TearingKind, // 0 ==> Strict, 1 ==> Casual
    // int casualTearingSet;

    // Newton algorithm
    x: Box<[f64]>,
    x0: Box<[f64]>,
    x_start: Box<[f64]>,
    x1: Box<[f64]>,
    finit: Box<[f64]>,
    fx0: Box<[f64]>,
    /// n times n Jacobian matrix with additional scaling row at the end in
    /// column-major format
    f_jac: Box<[f64]>,
    /// n times n Jacobian matrix with additional scaling row at the end in
    /// column-major format
    f_jacx0: Box<[f64]>,

    // /* debug arrays */
    // double* debug_fJac;
    // double* debug_dx;

    // Homotopy parameters
    /// Homotopy method used for the initialization with lambda from the homotopy()-operator
    init_homotopy: bool,
    start_direction: Direction,
    tau: f64,
    y0: Box<[f64]>,
    y1: Box<[f64]>,
    y2: Box<[f64]>,
    yt: Box<[f64]>,
    dy0: Box<[f64]>,
    dy1: Box<[f64]>,
    dy2: Box<[f64]>,
    hvec: Box<[f64]>,
    h_jac: Box<[f64]>,
    h_jac2: Box<[f64]>,
    h_jac_init: Box<[f64]>,
    ones: Box<[f64]>,

    // Linear system
    ind_row: Box<[usize]>,
    ind_col: Box<[usize]>,

    // int (*f)         (struct DATA_HOMOTOPY*, double*, double*);
    // int (*f_con)     (struct DATA_HOMOTOPY*, double*, double*);
    // int (*fJac_f)    (struct DATA_HOMOTOPY* solverData, double* x, double* fJac);
    // int (*h_function)(struct DATA_HOMOTOPY*, double*, double*);
    // int (*hJac_dh)   (struct DATA_HOMOTOPY*, double*, double*);
    user_data: Option<NonlinearSystemUserData>,
    eq_system_number: Option<u32>,
    time_value: Option<f64>,
    mixed_system: Option<bool>,
    // DATA_HYBRD* dataHybrid;
}

pub fn allocate_homotopy_data(
    size: u32,
    user_data: Option<NonlinearSystemUserData>,
) -> DataHomotopy {
    let jac_cols = user_data
        .as_ref()
        .map_or(0, |d| d.analytic_jacobian.size_cols)
        .max(size + 1);

    DataHomotopy {
        initialized: false,
        n: size,
        m: size + 1,
        xtol_sqrd: CONFIGURATION.newton_x_tol * CONFIGURATION.newton_x_tol,
        ftol_sqrd: CONFIGURATION.newton_f_tol * CONFIGURATION.newton_f_tol,
        error_f_sqrd: 0.0,

        res_scaling: vec![0.0; size as usize].into_boxed_slice(),
        fvec_scaled: vec![0.0; size as usize].into_boxed_slice(),
        hvec_scaled: vec![0.0; size as usize].into_boxed_slice(),
        dx_scaled: vec![0.0; size as usize].into_boxed_slice(),

        f1: vec![0.0; size as usize].into_boxed_slice(),
        f2: vec![0.0; size as usize].into_boxed_slice(),
        grad_fx: vec![0.0; size as usize].into_boxed_slice(),

        number_of_iterations: 0,
        number_of_function_evaluations: 0,
        max_number_of_iteartions: size * 100,

        causal_tearing_set: Unknown,

        x: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        x0: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        x_start: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        x1: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        finit: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        fx0: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        f_jac: vec![0.0; (size * jac_cols) as usize].into_boxed_slice(),
        f_jacx0: vec![0.0; (size * jac_cols) as usize].into_boxed_slice(),

        init_homotopy: false,
        start_direction: Direction::Forward,
        tau: CONFIGURATION.tau_start,

        y0: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        y1: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        y2: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        yt: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        dy0: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        dy1: vec![0.0; (size + CONFIGURATION.backtrace_strategy) as usize].into_boxed_slice(),
        dy2: vec![0.0; (size + 1) as usize].into_boxed_slice(),
        hvec: vec![0.0; size as usize].into_boxed_slice(),
        h_jac: vec![0.0; (size * jac_cols) as usize].into_boxed_slice(),
        h_jac2: vec![0.0; ((size + 1) * (size * jac_cols)) as usize].into_boxed_slice(),
        h_jac_init: vec![0.0; (size * jac_cols) as usize].into_boxed_slice(),
        ones: vec![0.0; (size + 1) as usize].into_boxed_slice(),

        ind_row: vec![0; (size + CONFIGURATION.backtrace_strategy - 1) as usize].into_boxed_slice(),
        ind_col: vec![0; (size + CONFIGURATION.backtrace_strategy) as usize].into_boxed_slice(),

        user_data,
        eq_system_number: None,
        time_value: None,
        mixed_system: None,
    }
}
