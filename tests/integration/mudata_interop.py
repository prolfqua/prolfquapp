"""Exercise R -> Python -> R MuData, including singleton axes and boolean masks.

Run from the package root with:
    uv run --no-project --with mudata==0.4.1 tests/integration/mudata_interop.py
"""

from pathlib import Path
import subprocess
import tempfile

import mudata
import numpy as np


def run_r(code: str, *paths: Path) -> None:
    subprocess.run(
        ["Rscript", "-e", code, *map(str, paths)], check=True,
    )


def main() -> None:
    """Verify complete modalities and maps using both native implementations."""
    with tempfile.TemporaryDirectory() as directory:
        original = Path(directory) / "r.h5mu"
        python_copy = Path(directory) / "python.h5mu"
        run_r('''
          devtools::load_all(quiet = TRUE)
          paths <- commandArgs(TRUE)
          obs <- data.frame(group = "A", row.names = "sample1")
          a <- anndataR::AnnData(X = matrix(c(2, NA), 1, 2), obs = obs,
            var = data.frame(protein = c("p1", "p1"), row.names = c("site1", "site2")),
            varm = list(present = matrix(c(TRUE, FALSE), 2, 1)))
          b <- anndataR::AnnData(X = matrix(1, 1, 1), obs = obs,
            var = data.frame(protein = "p1", row.names = "p1"))
          cf <- a[, 1]$clone(deep = TRUE)
          write_h5mu(list(enriched = a, total = b, cf = cf), paths[1], obs,
            list(stage = "test", version = "2.0.0"))
        ''', original)
        result = mudata.read_h5mu(original)
        assert list(result.mod) == ["enriched", "total", "cf"]
        assert list(result.obs_names) == ["sample1"]
        assert list(result.var_names) == ["site1", "site2", "p1"]
        np.testing.assert_array_equal(result.mod["enriched"].varm["present"], [[True], [False]])
        assert result.mod["cf"].varm["present"].shape == (1, 1)
        assert result.mod["cf"].varm["present"].dtype == np.bool_
        assert np.isnan(result.mod["enriched"].X[0, 1])
        for name in result.mod:
            assert result.obsm[name].shape == (1, 1)
        np.testing.assert_array_equal(result.varmap["cf"].ravel(), [1, 0, 0])
        result.write_h5mu(python_copy)
        run_r('''
          devtools::load_all(quiet = TRUE)
          paths <- commandArgs(TRUE)
          original <- read_h5mu(paths[1])
          restored <- read_h5mu(paths[2])
          stopifnot(identical(names(original$modalities), names(restored$modalities)),
            identical(original$obs, restored$obs), identical(original$uns, restored$uns))
          for (key in names(original$modalities)) {
            for (slot in c("X", "obs", "var", "varm")) {
              comparison <- .compare_mudata_slot(original$modalities[[key]], restored$modalities[[key]], slot)
              if (!isTRUE(comparison)) stop(key, "/", slot, ": ", paste(comparison, collapse = "; "))
            }
          }
        ''', original, python_copy)
        print("R -> Python -> R MuData compatibility passed")


if __name__ == "__main__":
    main()
