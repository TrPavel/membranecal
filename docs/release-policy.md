# Release policy

1. Select a clean, reviewed commit and preserve all frozen-study file hashes.
2. Validate the exact allowlist, licences, provenance, citation and reserved DOI.
3. Pass repository CI and full extracted-archive replay on the declared platform.
4. Build deterministic ZIP/checksum assets. Record commit, tree, inventory,
   manifest/archive hashes, CI URLs and independent clean-clone evidence in an
   external receipt to avoid circular self-hashes.
5. Complete `PUBLIC_RELEASE_OWNER_REVIEW.md`. Visibility change, tag creation,
   GitHub Release and Zenodo publication each require explicit owner authorization.
6. Only after that authorization, publish the approved bytes and verify downloaded
   assets. Do not rebuild or substitute them during publication.

This work package stops before step 5 approvals. No workflow creates tags,
GitHub Releases, Zenodo records or journal submissions. The existing reserved
Zenodo draft is reused; automatic GitHub-to-Zenodo integration is not enabled.

Study v0.2 remains immutable. A correction is additive and versioned, with its
effect on claims and comparability documented. Routine engineering PRs may use
squash merges; provenance-bearing scientific histories must not be rewritten.
Use required CI checks and protected release tags when configuring publication.

Current metadata has no invented publication date or article DOI. The primary
citation is the exact reproducibility package; any future article remains a
separate related object.
