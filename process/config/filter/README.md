# Filter Configurations

`filter` copies only the configured trees and branches from a trigger-like ROOT
file into a new ROOT file. A configuration contains tree sections followed by
one branch name per line:

```text
[frames]
spill

[timing]
nhits
time
```

Blank lines and comments beginning with `#` are ignored. Every requested tree
and branch must exist in the input; missing items are reported as errors. The
output retains all entries of each selected tree and preserves the original
ROOT branch types and leaflists. Unlisted trees and branches are omitted.

For variable-length arrays, include their count branch, such as `nhits`, when
selecting the array. `template.conf` lists the current trigger output trees and
branches, including timing-estimator branches.

Example:

```bash
process/bin/filter \
  --input triggered.root \
  --output filtered.root \
  --config process/config/filter/template.conf
```
