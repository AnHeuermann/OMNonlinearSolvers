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

use log::trace;

/// Print unknowns of the system
///
/// From `void printUnknowns(int logName, DATA_HOMOTOPY *solverData)`
pub(crate) fn print_unknowns(
    names: &[&str],
    x: &[f64],
    x_scaling: &[f64],
    min_value: &[f64],
    max_value: &[f64],
) {
    // Early return
    if !log::log_enabled!(log::Level::Trace) {
        return;
    }

    assert_eq!(
        names.len(),
        x.len(),
        "Array lengths of `names` and `x` aren't equal!"
    );
    assert_eq!(
        names.len(),
        x_scaling.len(),
        "Array lengths of `names` and `x_scaling` aren't equal!"
    );
    assert_eq!(
        names.len(),
        min_value.len(),
        "Array lengths of `names` and `min_value` aren't equal!"
    );
    assert_eq!(
        names.len(),
        max_value.len(),
        "Array lengths of `names` and `max_value` aren't equal!"
    );

    for i in 0..names.len() {
        trace!(
            "[{0}] {1} = {2} nom = {3} min = {4} max = {5}",
            i + 1,
            names[i],
            x[i],
            x_scaling[i],
            min_value[i],
            max_value[i]
        );
    }
}

// Unit tests
#[test]
fn print_unknowns_test() {
    let _ = env_logger::builder()
        .filter_level(log::LevelFilter::Trace)
        .is_test(true)
        .try_init();

    let names = ["var1", "var2"];
    let x = [1.0, 200.0];
    let x_scaling = [1.0, 1.0 / 100.0];
    let min_value = [-100.0, -100.0];
    let max_value = [100.0, 200.0];

    print_unknowns(&names, &x, &x_scaling, &min_value, &max_value);
}
