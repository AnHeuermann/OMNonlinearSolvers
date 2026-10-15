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

//! Runtime state the solver reads and writes, apart from the model functions.
//!
//! NONLINEAR_SYSTEM_DATA
//! - read: `size`, `equationIndex`, `homotopySupport`, `initHomotopy`,
//!   `mixedSystem`, `min`, `max`, `nominal`, `nlsxExtrapolation`, `nlsx`
//! - written: `nlsx` (solution), `numberOfFEval`, `numberOfIterations`,
//!   `numberOfJEval`, `numberOfCall`, `jacobianTime`
//!
//! DATA
//! - read: `localData[0]->timeValue`, `simulationInfo->discreteCall`,
//!   `simulationInfo->initial`, `simulationInfo->relations`,
//!   `simulationInfo->relationsPre`, `simulationInfo->homotopySteps`,
//!   `simulationInfo->nlsLinearSolver`, `modelData->nRelations`,
//!   `modelData->modelFilePrefix`, `callback->homotopyMethod`
//! - written: `simulationInfo->solveContinuous`

// Fallback solver, from nonlinearSolverHybrd.h:
// DATA_HYBRD* allocateHybrdData(size_t size, NLS_USERDATA* userData)
// void freeHybrdData(DATA_HYBRD* hybrdData)
// NLS_SOLVER_STATUS solveHybrd(DATA *data, threadData_t *threadData, NONLINEAR_SYSTEM_DATA* nlsData)

// Logging, from util/omc_error.h:
// void infoStreamPrint(int stream, int indentNext, const char *format, ...)
// void warningStreamPrint(int stream, int indentNext, const char *format, ...)
// void (*messageClose)(int stream)

// Timers, from util/rtclock.h:
// void rt_ext_tp_tick(rtclock_t* tick_tp)
// double rt_ext_tp_tock(rtclock_t* tick_tp)
