# Adamson 2016 inputs

These sample sheets use the submitted Cell Ranger 1.x BAMs. The v1 read layout
is reconstructed from the BAM tags: `CR+UR` provides the cell barcode and UMI,
and the long `SEQ` read contains the captured guide barcode at bases 41--60.
The pipeline trims that 20-base feature before KITE counting.

The pilot feature table contains the eight observed guide barcodes coupled to
the pilot H5AD labels. The UPR table contains the high-confidence barcode
couplings recovered from the ten guide gemgroups; labels below the pipeline's
minimum-cell threshold are not included.
