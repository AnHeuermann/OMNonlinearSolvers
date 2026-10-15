/*
 * This file belongs to the OpenModelica Run-Time System
 *
 * Copyright (c) 1998-2026, Open Source Modelica Consortium (OSMC), c/o Linköpings
 * universitet, Department of Computer and Information Science, SE-58183 Linköping, Sweden. All rights
 * reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * AGPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.8. ANY
 * USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES RECIPIENT'S
 * ACCEPTANCE OF THE BSD NEW LICENSE OR THE OSMC PUBLIC LICENSE OR THE AGPL
 * VERSION 3, ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium) Public License
 * (OSMC-PL) are obtained from OSMC, either from the above address, from the URLs:
 * http://www.openmodelica.org or https://github.com/OpenModelica/ or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica distribution. GNU
 * AGPL version 3 is obtained from: https://www.gnu.org/licenses/licenses.html#GPL. The BSD NEW
 * License is obtained from: http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY
 * SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF
 * OSMC-PL.
 *
 */

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
