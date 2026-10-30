# T05 expanded composition pilot — prescoring design

Question: with ESM-2 8M fixed, how does combining two complete proteins change
retrieval of **each** component relative to that component alone? Examine length
fraction, component identity and order. Independently implemented, conceptually
inspired by Nishit Kalani's supplied presentation (methods pp. 11–14), without
his code and without reproducing his reported numerical results.

## Selection, before scoring

Use the named baker's-yeast S288c genes in `config/t05_expanded_selection.json`.
These ordinary chromatin, ribosomal, metabolic and protein-maintenance proteins
are a deliberately related convenience sample, not a random proteome sample.
Resolve gene names to reviewed UniProtKB/Swiss-Prot accessions using official
EMBL-EBI Proteins API records. Retain every eligible resolved accession, deduplicate
accessions, retain and explicitly group identical sequences. A single accession
can describe multiple genes. Do not fill a missing or excluded slot after scoring.
Reject unreviewed, incorrect-organism, ambiguous, noncanonical, over-1,022-aa,
hazard-associated or permission-unclear candidates. Review names and functional
annotations individually before approving. Record failures and exclusions.

Select the 16 named targets in the selection file for varied natural lengths and
functions. If one fails inclusion, report its exclusion; do not choose a replacement.
Enumerate **all unordered pairs** of distinct admitted target accessions. Exclude
pairs longer than 1,022 residues and record their lengths and reason. Test both
orders of every remaining pair. This avoids picking favorable examples or forcing
length ratios. Retain all clean targets even if some pairs are excluded.

T02's provenance, local-storage, licensing and record review gates apply. The
user's explicit instruction authorizes this protein-only digital T05 extension;
it does not amend T02 or approve a DNA construct benchmark. Codex's record review
under that instruction is not advisor or institutional certification. Attribute
UniProt data to the UniProt Consortium under CC BY 4.0, retain the source's
other-rights disclaimer and original annotation bytes locally. No sequences,
source annotation copies, weights or numeric output files enter Git.

## Fixed construction and scoring

Use unchanged `config/t05_model.json`, original encoder and last-layer mean pooling
over actual residues. CPU float32, four threads, seed zero, eval/inference mode;
no truncation. Concatenate full canonical amino-acid strings directly. No linker,
mutation, repeat, trimming, translation or claim of a real/functional fusion.
Known intervals come from construction. Independently re-encode **each** extracted
component for each ordered query, rather than reuse its clean/catalogue embedding.

Encode catalogue entries singly, caching with model revision/file hashes, settings,
encoder code hashes, numerical-library versions and sequence hashes. Validate cache
metadata and bytes before reuse. Clean and oracle encodes bypass the cache. Store
complete accession rankings, exact-score ties, near ties within 1e-6, rank interval,
cosines, reciprocal ranks, and changes from the component's own clean query.
Exact ties use accession ascending (the original ranking rule). Also report a
duplicate-sequence-aware hit at each depth; this does not replace accession ranks.

For combined queries, recall@k = number of its two expected accessions in top k / 2;
thus maximum recall@1 is **0.5**. Report any-component hit@k and both-component hit@k
separately. Component recall/hit@k uses each component as one observation. Clean
and oracle controls have one expected component. Partner hits are valid, not errors.

Prespecified fraction bins: (0,.25), [.25,.40), [.40,.60), [.60,.75), [.75,1).
Report by bin, N/C-terminal position, and target identity; compare reversal within
each protein/pair. Summarize per-target means as well as all component instances.
No hypothesis tests or independent-observation claims: target reuse, pair reversals,
homologous catalogue entries and identical sequences create dependencies.

## Checks and evidence

Before expanded scoring, freeze reviewed manifest, pairs, queries, subset IDs,
selection, licence and model hashes in a local Git commit and a freeze descriptor.
The first timing/memory check uses only original approved inputs, without ranking
or influencing new selection. Preserve original experiment files and output hashes.

Run original tests and meaningful new checks: duplicate/tie handling, recall with
two valid components, length exclusions, reversals, oracle slicing, cache key
invalidation/corruption, manifest tampering and numerical cosine agreement.
Retain real padded-versus-single pooling check. Repeat, in a separate process with
catalogue re-encoding, a subset chosen before scoring: shortest-total, longest-total,
most unequal and most balanced allowed pair (ties accession order), both orders,
plus all 16 clean targets. Require rank agreement and score/vector abs error <=1e-6.
Report clean failures, numerical checks, exceptions and exclusions, not silently
drop cases. Save peak resident memory and separated loading, catalogue/query
encoding and ranking timings. Source acquisition and plotting are separate.

Self-retrieval and oracle success are wiring checks, not generalization. Order
comparisons are controlled within this digital model experiment; length/identity
associations are confounded. Extracting and re-encoding changes context as well
as which residue vectors enter the mean. No function prediction, automatic
segmentation, biological construct validation or hazard screening is established.
