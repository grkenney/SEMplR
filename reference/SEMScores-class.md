# Class for storing SEM motif binding calculations for multiple genomic ranges or variants

Class for storing SEM motif binding calculations for multiple genomic
ranges or variants

## Slots

- `ranges`:

  A `GRanges` or `VRanges` object to hold one or more genomic ranges

- `semData`:

  A `data.table` object of metadata for the SEMs with one row for each
  SEM

- `scores`:

  A `data.table` object for motif information and binding scores
