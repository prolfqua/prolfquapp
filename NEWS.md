# prolfquapp 2.3.3

- Protein abundance QC report (`QC_ProteinAbundances.Rmd`): contaminant and decoy proteins are now drawn on top of the regular points instead of being hidden behind them, and are rendered with a separate, higher opacity so they stand out. `plot_abundance_vs_percent()` gains a `highlight_alpha` argument (default `1`, so existing callers are unchanged); the report sets it to `0.8`.
- DiffExpQC report (`DiffExpQC_R6.Rmd`): reverted the four combined figures from interactive `plotly` subplots back to static `gridExtra::grid.arrange()` panels.
