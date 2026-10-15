# TODO

* Use a package for linear algebra, e.g. [nalgebra][nalgebra] (or similar),
  instead of flat vectors
* Use boolean or enum for `causal_tearing_set`
* Split DataHomotopy into multiple structs
  * Newton data
  * Homotopy parameters
  * Scaling vectors
  * Tolerances and settings
  * Statistics
* Check how to include hybrj_ from Minpack
* Move `max_number_of_iteartions` to configuration
* Replace our implementation of a newton algorithm with
 [newton_rootfinder][newton_rootfinder] or [nrfind][nrfind] (or similar)

[nalgebra]: https://nalgebra.rs
[newton_rootfinder]: https://github.com/Nateckert/newton_rootfinder
[nrfind]: https://github.com/NoraCodes/nrfind
