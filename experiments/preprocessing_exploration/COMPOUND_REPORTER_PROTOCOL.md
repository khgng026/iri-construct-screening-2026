# Compound-annotation reporter follow-up

The archive coverage run found a documented, multi-part EGFP coding region
in the public LAMPRA plasmid map. No one contiguous translation from the
existing candidate methods contained its full peptide. This follow-up
asks whether existing protein-to-DNA annotation can propose a
usable peptide from that same benign record.

Use the already pinned archive record and the same seven distinct natural
GFP references as the whole-DNA translated-search pilot. The independent
reference proteins exclude the exact EGFP peptide. Run miniprot 0.18 in
GFF mode, supplying only the archived DNA and reference proteins. Rank
predicted transcript models by their reported score, without annotations.
Check annotated CDS overlap, exon interval agreement and recovery of the
documented translated peptide in the top model and across all returned
models. Report whether candidate translation matches exactly, allowing a
terminal stop to be excluded for comparison. The map contains unresolved
source bases away from the annotated EGFP coding region.

The three annotated intervals are adjacent. Evaluate them as a compound
annotation and report the intervening base count separately. This one
reporter is a feasibility example on an imperfect archived map.
Reference-guided prediction requires a suitable reference. The comparison
does not estimate generalization to deeply divergent proteins.

Also compare a translation variant that keeps fully resolved peptide starts
and ends when another part of the archived translation contains an unknown
residue. Leave the input DNA and the evaluation label intact. Evaluate
coverage and presence of the exact documented reporter at caps 1, 2, 4,
8 and 16 and across the full pool. Recheck the yeast archive to measure how
the variant behaves on the remaining records.

[Miniprot](https://doi.org/10.1093/bioinformatics/btad014) ·
[LAMPRA plasmid maps](https://github.com/shendurelab/LAMPRA/tree/2aac5500935a8b9e6452288945f611789d839108/plasmid_maps)
