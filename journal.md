# Development journal

## 2026-11-13 — Is scikit-bio 0.5.9 outdated? (to revisit)

**Short answer:** yes, 0.5.9 is outdated (released Aug 2023; current line is 0.7.x), but
scikit-bio is still the right library — it is actively developed and the standard for this
kind of work. The question is *when* to migrate, not *what* to replace it with.

### What changed upstream since 0.5.9 that matters for this code
- **`otu_ids` → `taxa`** (0.6.0): the argument used in the two `beta_diversity(...)` UniFrac
  calls in `microbiomeutils.py`. `otu_ids` is still accepted as a deprecated alias, but
  upstream says it will be removed in a future release.
- **numpy ≥ 2.0 required** from 0.7.1.
- **Python:** 3.9 dropped in 0.7.1, 3.10 dropped in the upcoming 0.7.5 → target Python 3.11+.
- **`pcoa`:** `number_of_dimensions` → `dimensions` and `distance_matrix` → `distmat`
  renamed in 0.7.0. We do not use those arguments, so no impact.
- **UniFrac:** 0.7.4 added a faster `engine="fast"` / parallel UniFrac — a nice bonus.
- **File formats:** no changes found to lsmat distance-matrix or ordination-results formats.

### Availability on the Alliance clusters
- The wheelhouse has `scikit-bio 0.7.1.post1` for cp311, so an upgrade does not require
  internet pip installs.

### The real blocker: emperor
- The wheelhouse only has `emperor 1.0.1` (Python 2 only); the `meco/python/3.11.6` module
  has `emperor 1.0.4`.
- emperor is barely maintained; its compatibility with scikit-bio 0.7 + numpy 2 is
  **unverified**.

### Recommended plan
1. **Stay on 0.5.9 for now** — it is the version verified against the reference outputs in
   `data/` (`meco/python/3.11.6`: Python 3.11.6, scikit-bio 0.5.9, numpy 1.26.4,
   scipy 1.10.1, pandas 2.2.0, emperor 1.0.4). Pinned in `requirements.txt`.
2. **Test 0.7.1 separately on `dev`:** throwaway venv from the wheelhouse (scikit-bio 0.7.1,
   numpy 2, emperor 1.0.4), rerun the regression comparisons and the unit tests, check
   whether `emperor` still works. Do this in a short `salloc`, not on the login node.
3. **Then:**
   - if emperor works → switch to `taxa=` and require scikit-bio ≥ 0.6;
   - if not → keep `emperor` as an optional subcommand with its own pinned environment, or
     replace it with a plain plotting library for the 3D PCoA.

Possible stopgap: make `betadiv` pass `otu_ids=` or `taxa=` depending on the installed
scikit-bio version, so it works on both.

### Sources
- [scikit-bio CHANGELOG.md](https://github.com/scikit-bio/scikit-bio/blob/main/CHANGELOG.md)
- [scikit-bio on PyPI](https://pypi.org/project/scikit-bio)
- [scikit-bio releases](https://github.com/scikit-bio/scikit-bio/releases)
- [scikit-bio 0.5.9 on Zenodo (Aug 2, 2023)](https://zenodo.org/records/8209901)
