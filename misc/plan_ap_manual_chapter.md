# Plan: Manual chapter for acyclic partitions

Status: not started. API reference coverage (`docs/src/apicore/ap.md`,
`plot_morse_strata` in `apicore/plots.md`) is already done — this plan
covers only the narrative "Manual" part of the docs, which was
deliberately deferred.

This includes the top-layer functions `is_top_layer`, `stratum_top_layer`,
and `construct_ap_stratum_top` (added to `src/ap/stratum_top_layer.jl`):
their `apicore/ap.md` entries (under a "Top Layer" heading) are already in
place, so only the narrative chapter coverage described in step 7a below
remains outstanding for them, same as for everything else on this page.

## Goal

Add a `man/ap.md` chapter, analogous in depth and style to the existing
`man/conley.md`, that explains acyclic partitions `AP(X)`/`AP^c(X)`,
the Morse-vector stratification, and the atomic-refinement coverage
theory, with worked code examples for every exported function.

## Source material

`construct_apx/doc/morse-vector-strata.md` (793 lines, sections 0-8 plus
appendix) is the primary source — it's the paper note the `ap/` code
implements. Section map, for picking what to condense:

| Section | Content | Relevance to chapter |
|---|---|---|
| 0 | What is proved | Framing/motivation paragraph |
| 1.1-1.3 | Lefschetz complexes, Conley index recap, `AP(X)`/Morse vector definitions | Core definitions section |
| 1.4 | The question (does refinement always drop the Morse vector uniformly?) | Motivates the whole chapter |
| 2 | Two structural reductions | Skip or one paragraph — proof machinery, not needed for API users |
| 3.1-3.2 | Rank dictionary / deficiency `rho(V)`, stratum invariant | `block_rho` section |
| 3.3 | Splits | `block_split_deltas`/`split_spectrum`/`connected_split_deltas` section |
| 4-6 | Obstruction, realization theorem, counterexamples | Background for Theorem 7.1 classification; condense to a few sentences, skip proofs |
| 7 | Classification (uniform vs. non-uniform jumps) | `stratum_adjacency` section — this is the theorem the coloring in `plot_morse_strata` visualizes |
| 8 | What remains open | Optional "Further reading" note |
| Appendix | Exhaustive machine verification | Not needed — mention `test/test_ap.jl` instead as the code-side analog |

The top layer (`is_top_layer`, `stratum_top_layer`,
`construct_ap_stratum_top`) is **not** covered by
`morse-vector-strata.md` and has no other citable source at this time —
it comes from separate, unpublished work-in-progress notes that must not
be named or cited anywhere in the docs. Write step 7a's content directly
from the functions' own docstrings and the definitions already
established earlier in the chapter (atomic refinement, Morse vector,
`block_split_deltas`/`connected_split_deltas`), phrased independently
rather than summarized from any outside note. See "Citation" below.

Do **not** try to port the proofs (sections 2, 4-6) — the manual chapters
elsewhere in this package explain definitions/theorems and how to use the
code, not proofs. Match that register.

## Proposed chapter outline (`docs/src/man/ap.md`)

Mirror `man/conley.md`'s structure: `#` title, then `##` sections, each
with a `!!! tip "Definition: ..."` or `!!! tip "Theorem: ..."` admonition
for formal statements (existing convention throughout `man/*.md`),
followed by a runnable code example using `create_simplicial_complex` or
similar small examples (the package favors small, fully-printed examples
over abstract descriptions — see `man/tutorial.md`, `man/conley.md`).

1. `# Acyclic Partitions`
   - One-paragraph motivation (from source section 0/1.4): given a
     Lefschetz complex, its multivector fields form a poset under atomic
     refinement; does refining a multivector field always change its
     Morse vector uniformly?
2. `## Acyclic Partitions and the Morse Vector`
   - Definitions of `AP(X)` (full unrestricted set of MVF partitions),
     `AP^c(X)` (connected sub-poset), Morse vector `M(W) = sum beta(V)`
     over critical blocks (from source 1.3).
   - `construct_ap_space` example (small triangle/tetrahedron complex,
     `connected=true` vs `connected=false`, matching the caveat already
     written up in `misc/changelog` / the migration note about
     `connected=false` blowing up).
   - `construct_ap_stratum` example, including the documented caveat
     about no pruning benefit when targeting the global-minimum Morse
     vector.
3. `## Converting Between Representations`
   - `convert_partition_mvf` / `convert_mvf_partition`, brief, since this
     is plumbing rather than theory.
4. `## Atomic Refinement`
   - `is_atomic_refinement`, `atomic_distances`, with a small worked
     example showing a covering relation.
5. `## Connectivity`
   - `is_connected_block` / `is_connected_partition` / `is_closed_in_block`
     — ties back to why `AP^c(X)` is the more commonly used variant
     (tractability, from source 1.3 and the migration note's caveats).
6. `## The Rank Vector and Block Splits`
   - `block_rho` (Lemma 3.1, source 3.1-3.2).
   - `block_split_deltas` / `split_spectrum` / `connected_split_deltas`
     (source 3.3) — brute-force split enumeration, example on a single
     block.
7. `## Stratification and Refinement Coverage`
   - `stratum_partition` grouping by Morse vector.
   - `stratum_adjacency` and `stratum_jump_weight`, stating the
     Theorem 7.1 classification (source section 7) informally: when is
     coverage between adjacent strata uniform vs. not. This is the
     mathematical payoff of the chapter — spend the most space here.
7a. `## The Top Layer of a Stratum`
   - New subsection, placed right after "Stratification and Refinement
     Coverage" since it refines the same idea (a stratum's *internal*
     structure) rather than coverage *between* strata.
   - Definition admonition, phrased from scratch (no outside citation —
     see the note above the section-map table): within a fixed Morse
     stratum, the top layer is the set of partitions that admit no
     further atomic refinement preserving that Morse vector — i.e. every
     explicit block's split spectrum (section 6) avoids the zero delta.
   - `is_top_layer` worked example: take one partition from a stratum and
     check it directly; a good place to show a partition that *fails*
     the test alongside one that passes, so the "can still be refined
     conservatively" failure mode is visible, not just asserted.
   - `stratum_top_layer` example: filter a full `stratum_partition` result
     down to each stratum's top layer, on the same running example used
     in section 7 (reuse the complex/`ap` already built there rather than
     introducing a new one).
   - `construct_ap_stratum_top` example: the direct one-target-vector
     shortcut, contrasted with `construct_ap_stratum` (section 2) the
     same way `construct_ap_stratum` itself was contrasted with
     `construct_ap_space` there.
   - One or two sentences noting a Forman vector field is always in the
     top layer of its stratum, tying back to the Forman-pair language
     used in `man/tutorial.md`/`man/conley.md`.
   - Fold in (or forward-reference from) a caveat: this inherits
     `block_split_deltas`'s brute-force `length(block) <= 20` limit, and
     `stratum_top_layer` memoizes per-block rigidity across the whole
     call (worth one sentence if section 9 below ends up thin).
8. `## Visualizing Strata`
   - `plot_morse_strata` worked example with an embedded image, following
     the pattern in `man/plotting.md`'s "Plotting Morse Sets" subsection
     (code block + rendered image + brief explanation of the color coding
     tied back to Theorem 7.1).
9. `## Caveats` (or fold into relevant subsections instead of a separate
   section — check what reads better once drafted)
   - Performance notes from the migration doc's Section 8: pruning isn't
     automatic, `connected=false` can be intractable, `atomic_distances`
     emits `@info` progress lines.
   - Top layer's own caveat from step 7a (`length(block) <= 20`), if not
     already folded in there.
10. Optional closing `## [References](@id refap)` section only if/when
    the source paper gets a real citation (see below).

## Wiring

- `docs/make.jl`: add `"man/ap.md"` to the `"Manual"` list, after
  `"man/conley.md"` and before `"man/flowexamples.md"` (mirrors the
  `apicore/ap.md` placement already done).
- `docs/src/index.md`: add `"man/ap.md"` to the `@contents` `Pages` list
  in the same position, and extend the "Manual Outline" prose paragraph
  with one sentence introducing the new chapter.

## Citation

`construct_apx/doc/morse-vector-strata.md`/`.pdf` is not yet in
`docs/src/refs.bib`. If/when it's posted as a preprint or published,
add a bibtex entry (following the existing key convention,
e.g. `wanner:2Xa`) and cite it from the new chapter's intro paragraph,
matching how `man/conley.md` and others cite `wanner:25a` etc. Until
then, skip the References subsection rather than inventing a citation.

Step 7a (the top layer) is a separate case: it is **not** part of
`morse-vector-strata.md` at all, and its own source is unpublished
work-in-progress that must not be mentioned, named, filenamed, or cited
anywhere in the docs, in code comments, or in commit messages — not even
once the rest of the chapter gets a real citation above. If that source
is ever published, revisit this note before adding any reference to it.
Until then, step 7a stands entirely on its own definitions and the
package's own docstrings/tests, with no References entry of its own.

## Suggested effort / sequencing

This is a genuine prose-writing task (the existing `man/*.md` chapters
run 400-1300 lines each with full math typesetting and worked examples),
not a mechanical wiring job like the API reference page was. Recommend
drafting section-by-section against a running `make.jl --local-html`
build (catches broken `@ref`/`@docs`/math-macro issues immediately) and
reviewing the rendered HTML rather than raw markdown, same as was done
for the `apicore/ap.md` page.
