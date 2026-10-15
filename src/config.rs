// (C) Andreas Heuermann, 2026

pub struct Config {
    pub newton_x_tol: f64,
    pub newton_f_tol: f64,
    pub tau_start: f64,
    pub backtrace_strategy: u32,
}

pub static CONFIGURATION: Config = Config {
    newton_x_tol: 1e-12,
    newton_f_tol: 1e-12,
    tau_start: 0.2,
    backtrace_strategy: 1,
};
