# Vendored dependencies

Builds never download dependencies. `make verify-vendor` verifies SHA-256 hashes
in `SHA256SUMS`. Headers and all embedded copyright/license notices are intact.

| Component | Version | Pinned upstream |
| --- | --- | --- |
| csv-parser | 5.3.0 | `vincentlaucsb/csv-parser`, tag commit `32e99be14236b0585f33fb96f37e4fefc363f448` |
| doctest | 2.4.12 | `doctest/doctest`, tag commit `1da23a3e8119ec5cce4f9388e91b065e20bf06f5` |

csv-parser's source-tree `single_include/csv.hpp` is a compatibility shim.
The vendored file is the [5.3.0 release asset](https://github.com/vincentlaucsb/csv-parser/releases/download/5.3.0/csv.hpp),
verified against its [published checksum](https://github.com/vincentlaucsb/csv-parser/releases/download/5.3.0/csv.hpp.sha256).
`csv-parser/LICENSE` comes from that tag commit. The release workflow generates
the asset from `single_header.json` using upstream `single_header` 0.1.0.
The parser embeds mio (MIT) and string-view-lite (Boost Software License 1.0);
their notices remain inside `csv.hpp`. No extra linked runtime is introduced.
`CSV_ENABLE_THREADS=0` and `CSVFormat::threading(false)` ensure serial parsing.

doctest's [header](https://raw.githubusercontent.com/doctest/doctest/1da23a3e8119ec5cce4f9388e91b065e20bf06f5/doctest/doctest.h)
and [MIT license](https://raw.githubusercontent.com/doctest/doctest/1da23a3e8119ec5cce4f9388e91b065e20bf06f5/LICENSE.txt)
come directly from the pinned commit. The test header is used only in the test
binary. 2.4.12 was selected for identifiable release provenance; the earlier
AntRepCLA copy's 2.5.0 label was not used as release evidence.

Context7 returned unrelated JavaScript/Python projects for these C++ libraries;
API decisions were checked against the pinned upstream header/source instead.
