# Retained GMM Stan sources

These six files are byte-for-byte copies of the statistical model sources in
`ot_perturbation_catalog_mave_experiments/src/` at commit
`da4b24f2f2ad6a8d86d10831e34699199c748300`. Each model includes
`functions.stan`; the include relationship is unchanged.

| File | SHA-256 |
| --- | --- |
| `functions.stan` | `c9f4083ac3130cc45a4d6e7e457b634912456feba0728f0adace0dca373d4ef9` |
| `fixed_gmm_flexible.stan` | `e8414ab496e441635b865c9fe18426a3eb16ed2047de00464779c1f10b94f81d` |
| `fixed_gmm_ko.stan` | `3c061e41c1bd855e4654f36f962de44704412028dba536ebbc0c8df8da14dc17` |
| `fixed_gmm_del.stan` | `bc750165746546c4db440aec41bb0a5c5896ffa74364ed2ec929af7488f73011` |
| `fixed_gmm_neut.stan` | `843eb7db6039e10ce07d56f5e1a68e385ae4ffc8ea7c2dbb6c1ab0b31e689357` |
| `fixed_gmm_pos.stan` | `beff4c7536f1695673365ed7dcdd2706d65470ddaa733ab9d41f1cae429c0861` |

The source revision has no Stan-file license headers, and the sibling checkout
contains no top-level `LICENSE` file or declared license. Source comments have
been retained unchanged; this note records provenance without adding a license
to the upstream code.
