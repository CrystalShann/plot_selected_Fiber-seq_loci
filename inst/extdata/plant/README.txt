Synthetic fixtures for native package tests (not biological plant data).
Chr5 and scaffold_42 each have three reads spanning BED [0,10000).
Two reads have m6A at one-based positions 2001 and 2101; the third has only
sentinel blocks and therefore no m6A calls. FiberHMM nucleosomes are 150 bp.
The 20-bp TF call has score 80 and the 50-bp call score 20 (filtered at 50).
Annotations exercise both strands, shared and alternative TSSs, and missing
canonical tags. These fixtures test file/coordinate compatibility, not upstream
FiberHMM model applicability to plants.
