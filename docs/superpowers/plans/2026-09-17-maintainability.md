# Maintainability Implementation Plan

> **For agentic workers:** Use superpowers:subagent-driven-development to implement and review each task.

**Goal:** Improve maintainability without changing validated interfaces or scientific behavior.
**Architecture:** Private workflow modules; explicit shared Make rules; readable self-contained R and helper scripts; indexed current/history documentation.
**Tech Stack:** C++20, GNU Make, R, Bash; Python-free compute.
**Spec:** docs/CODEBASE_REFACTOR_DESIGN.md (approved).

## Global constraints

Work on main, preserve existing README/license/evidence edits, no commits. Keep CLI, TSVs, numerical tolerances, process supervision, immutable snapshots and public script paths. Do not edit third-party code, raw data or scientific fixtures. Do not run concurrent heavy R checks. Production R stays self-contained so existing source hashing remains complete.

## Task 1: Native workflow boundaries

Own src/workflow.cpp, new src/workflow/*; optionally focused native tests. Root owns Makefile. Extract existing code into private config/state/execution modules with narrow internal declarations, leaving command dispatch in src/workflow.cpp. Preserve validation order, exception text, source inventory, lock lifecycle and signal semantics. Expand compressed control flow. Tell root exact new source list before building. Verify with full native check and isolated sanitizer build; no new implementation-mirroring tests.

## Task 2: Readable R and helper code

Own tools/*.R, run_deg_analysis_offline.R, optional .R formatter policy. Expand dense statements with consistent indentation and descriptive intermediate names only where needed. Preserve exact parsed expressions for mechanical formatting; retain original expression snapshots outside tracked source for comparison. Do not split production R into sourced modules. Run parse-equivalence checks, fixture/staging/toolchain helper checks. Coordinate expensive statistical check with root; independent oracle remains untouched.

## Task 3: Build, repository navigation and integration (root)

Replace repeated Make compile recipes with explicit source lists and shared rules, separate test objects under BUILD/tests; preserve BUILD/compiler/linker overrides/dependency includes/check/sanitize/workflow. Add .editorconfig, .clang-format, CONTRIBUTING.md, docs/README.md and minimal native GitHub CI. Move historical plans to docs/history and update links; retain current evidence and reports. Correct README C++20 requirement. Keep entry points stable.

## Verification and completion

- [ ] Existing native baseline passes before refactor.
- [ ] Native extraction passes full check and targeted sanitizer coverage.
- [ ] R/helper expressions and focused checks preserve behavior.
- [ ] Alternate BUILD and documented public launcher verified.
- [ ] Native CI/config/docs reviewed; links and diff whitespace checked.
- [ ] Small real installed-tool workflow and independent R oracle pass after integration.
- [ ] Independent final review; update status with measured evidence and limits.

Ruling: user already approved design and instructed execution; continue without additional planning approval. No automated commits. Use existing relevant tests instead of inventing failing tests for a behavior-preserving extraction.
