# Version identities

| Object | Current identity | Meaning |
| --- | --- | --- |
| Scientific study | MembraneCal v0.2 | Frozen cohort, methods and results |
| Repository/software package | 0.2.0, tag v0.2.0 | First standalone repository/package release |
| Repository metadata schema | 1.1.0 | Validation/derivation metadata contract |
| Numerical protocol | Existing frozen schema and SHA256 | Never renamed to match repository SemVer |
| Archive | External SHA256 | Exact distributed bytes |

Packaging/documentation fixes may advance 0.2.0 to 0.2.1 without changing study
identity, with a `scientific_changes: false` declaration. Additive tooling and
breaking interfaces follow SemVer, but a software version bump never silently
authorizes a changed scientific study. Changed cohorts/methods require a new study
identity and an explicit comparability statement.

Published tags are immutable. Citation
metadata records package version 0.2.0 while retaining study v0.2 in the title.
