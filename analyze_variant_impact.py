#!/usr/bin/env python
import sys
from pathlib import Path

# Add src to path and delegate to src.analyze_variant_impact
src_dir = Path(__file__).parent / "source"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from analyze_variant_impact import parse_args, run

if __name__ == "__main__":
    args = parse_args()
    run(args.data_dir, args.output_dir, args.show)
