# Adamson 2016 inputs

These sample sheets use the submitted Cell Ranger 1.x BAMs. The v1 read layout
is reconstructed from the BAM tags: `CR+UR` provides the cell barcode and UMI,
and the long `SEQ` read contains the captured guide barcode at bases 41--60.
The pipeline trims that 20-base feature before KITE counting.

The pilot feature table contains the eight observed guide barcodes coupled to
the pilot H5AD labels. The UPR table contains the high-confidence barcode
couplings recovered from the ten guide gemgroups; labels below the pipeline's
minimum-cell threshold are not included. Ten legacy UPR target symbols absent
from the reference GTF use explicit stable ENSG suffixes in the feature labels.
`test_adamson_2016_targets.py` checks all 92 non-control feature labels across
82 targets against the same GTF resolver used by the pipeline. Run it in the
pipeline image with `--gtf <GRCh38.115.gtf.gz>`.

The shared GEX BAM is streamed once and partitioned by the corrected-barcode
CB/BX gemgroup suffix. Each standard count task receives its staged per-group
BAM alongside the corresponding single-group guide BAM; alignment tags and
UMIs are retained.
