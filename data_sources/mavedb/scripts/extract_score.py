"""Deprecated command alias for process_mavedb.py's consolidated exports.

--output now names a directory containing scores.tsv and metadata.tsv.
"""

import logging

if __package__:
    from .process_mavedb import main
else:
    from process_mavedb import main


if __name__ == "__main__":
    logging.getLogger(__name__).warning(
        "This extractor is deprecated; use process_mavedb.py. "
        "--output is now a directory for both consolidated tables."
    )
    main()
