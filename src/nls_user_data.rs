// (C) 2026, Andreas Heuermann

use std::ffi::c_void;

pub struct Jacobian {
    /// Pointer to Jacobian of type `JACOBIAN*`
    jacobian: *mut c_void,

    /// Number of Rows
    pub size_rows: u32,
    /// Number of Columns
    pub size_cols: u32,
}

pub struct NonlinearSystemUserData {
    /// Pointer to data of C type `DATA*`
    data: *mut c_void,
    /// Pointer to thread data of C type `threadData_t*`
    thread_data: *mut c_void,

    /// System index
    sys_number: u32,
    /// Pointer to nonlinear system data of type `NONLINEAR_SYSTEM_DATA*`
    nls_data: *mut c_void,
    /// Analytic Jacobian
    pub analytic_jacobian: Jacobian,

    /// Optional pointer to ODE solver data of C type `void*`.
    /// Used in NLS solving of ODE integrator step.
    solver_data: *mut c_void,
}
