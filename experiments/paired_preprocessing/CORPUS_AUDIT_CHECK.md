# Corpus audit metadata check

The acquisition and audit stages in PR #19 record source identity and reproduce prior labels. The saved metadata needs a separate check before it can support a new evaluation. The checker works with the two private report files and prints aggregate findings. It reads JSON rather than sequence files and uses the Python standard library.

```sh
run=results/paired-corpus-batch-v1
python -m experiments.paired_preprocessing.check_corpus_audit \
  --acquisition "$run/acquisition.private.json" --audit "$run/audit.private.json"
```

The check compares the reported file hashes and version lists, reconciles the expected records with the saved rows, and recounts the audit flags. It flags absent or inconsistent region metadata, including spans outside their own CDS or across a gap in a joined CDS. Exit status 2 indicates findings or unreadable metadata. These are checks of the saved reports. Independent verification of the provider bytes requires the original source files.

## What the findings mean

- **Keep location evidence per occurrence.** The PR #19 audit uses a protein-level cache and stores its family dictionary directly in each output row. It then writes construct-specific coordinates into that same dictionary. Repeated proteins can make an earlier row inherit later coordinates. The original flag counts can still look correct because the check was performed before the later write. Separate cached evidence from per-record location evidence before trusting exported region rows. Containment detects some corruption even when flags pass. It cannot establish the exact location of a contained span. [Audit commit](https://github.com/khgng026/iri-construct-screening-2026/commit/6b55078f2c1eeda6e0c4aae4129c521a62f22964)
- **Separate anchor evidence and deposited evidence.** A region that translates consistently with an anchor hit supports an exploratory annotation. Counts of matching named source features can guide independent review of boundaries. A matching name alone leaves functional identity and independent verification open.
- **Check dependence before splitting.** The group diagnostic connects records that share any identical exported protein hash, including through another record. It reports how many groups there are and the size of the largest. Shared common components can connect otherwise distinct constructs. Use independently established lineage and component evidence to choose the final evaluation groups.

## Frozen metadata cross-check

Using the `parents` and exported `units` in the original and extension panel freezes, the shared-hash diagnostic gives 193 groups across 1,000 records, with 773 records in the largest group. Among units labelled `full`, `embedded`, or `multi_family`, 93 protein hashes recur in different constructs. Those hashes account for 1,143 labelled occurrences across 749 records. These counts show how often repeated proteins occur. They do not measure how many saved loci were overwritten, and the group diagnostic is not a proposed split. [Original freeze](../../docs/results/component-study-v1/panel.freeze.json) · [Extension freeze](../../docs/results/component-study-v1/extension-1/panel.freeze.json)

The private audit files from the reported run were unavailable in this workspace. The tests use fictional metadata to check that overwritten coordinates, gaps, record mismatches, and snapshot drift are detected.
