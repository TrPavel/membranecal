# Release policy

1. Select a reviewed commit and preserve every frozen-study file hash.
2. Validate the allowlist, licences, provenance, citation and DOI metadata.
3. Pass repository CI and full extracted-archive replay on the declared platform.
4. Build deterministic ZIP/checksum assets. Record commit, tree, inventory,
   manifest/archive hashes, CI URLs and independent clean-clone evidence in an
   external receipt to avoid circular self-hashes.
5. Obtain owner authorization outside the scientific release payload. Visibility,
   tag creation, GitHub Release and Zenodo publication each require authorization.
6. Publish the approved bytes and verify downloaded assets. Do not rebuild or
   substitute them during publication. Reuse the assigned Zenodo record.

No workflow creates tags, GitHub Releases, Zenodo records or journal submissions.
Automatic GitHub-to-Zenodo deposition is not part of this release process.

Study v0.2 remains immutable. Corrections are additive and versioned, with their
effect on claims and comparability documented. Routine engineering PRs may use
squash merges; provenance-bearing scientific histories must not be rewritten.
Use required CI checks and protected release tags when configuring publication.

Current metadata has no invented publication date or article DOI. The primary
citation is the exact reproducibility package; any future article is a separate
related object. DOI states in the metadata schema include `reserved_draft` and
`published`. The v0.2.0 publication object requires `published`; historical RC
records retain their original state. These are content states, not instructions
to change an external service. Deployment evidence and approvals remain external.
