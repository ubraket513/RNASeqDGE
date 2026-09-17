# Pinned binary toolchain provenance

`manifest.json` identifies the five core alignment packages (including HTSlib),
upstream source URLs and SHA-256 values declared in their recipes, binary package
SHA-256 values, licenses, runtime dependencies, and hashes of retained evidence.
The complete 48-package runtime closure is pinned in
`config/alignment-linux-64.explicit.txt` and its provenance JSON.

The package archives were SHA-256 checked against installed conda records before
copying these recipes, patches, build configurations and license texts. Version
outputs and hashes of the actual installed entry points are captured separately
by `tools/toolchain.py preflight`. Prefix relocation can change executable/wrapper
bytes; package archive hashes and installed-file hashes serve different purposes.

These are upstream build records, **not scripts invoked by this project**. The
selected deployment route stages the pinned binary packages. STAR's recipe builds
several SIMD variants and its dispatcher selects a supported variant at runtime.
The HISAT2 recipe contains an unpinned `simde-no-tests` clone; its full historical
source-build inputs cannot be reconstructed from the recipe alone. We therefore
claim reproduction of the locked binary environment, not a bit-identical source
rebuild. Upstream source archive hashes are recorded from recipes, not independently
download-verified source archives. No upstream build runs during compute.

The R/legacy staging environment has a separate 344-package lock. Its recorded
Bioconductor post-link downloads mean it is not an offline installation bundle.
Offline **execution** is tested after staging, with networking disabled.

To regenerate this evidence into an absent `vendor/toolchain` directory, use
`.deps/p0/bin/python tools/record_tool_recipes.py .deps/alignment`; it requires the
original package cache and PyYAML. Do not run the retained vendor build scripts.
