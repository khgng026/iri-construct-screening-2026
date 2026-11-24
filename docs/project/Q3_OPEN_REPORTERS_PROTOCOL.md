# Open Reporters source and preprocessing feasibility

The Open Reporters Collection releases circular cloning plasmid records in
GenBank and documents the reporter inserted in each. This creates a fresh
whole-DNA source for the project. It concentrates on one cloning architecture,
so it cannot supply diversity across construct layouts. The release dates of
these plasmid records are insufficient to qualify their encoded proteins for
a temporal model holdout. The project uses the source in place, keeps sequence
data out of this repository, and attributes the collection and protein sources.
[Open Reporters Collection](https://github.com/Reclone-org/Open-DNA-Collections/tree/55d4b762194055c8027fd6c780068d8a78a03282/Open%20Reporters%20Collection)

The preregistered inventory covers the 23 natural-derived GFP-like reporters
and five ordinary luminescent or colorimetric reporters from the source
catalogue. Exclude the engineered EiraCFP entry and the resistance-enzyme
entry from the inventory. Use the first unique annotated reporter CDS
identified from the source's single reporter and vector-marker feature
design, and check that its exactly annotated coding region translates with a
single terminal stop and no internal stops. A close match alone is insufficient
to promote a record to independently corroborated status.

Look for exact amino-acid matches of the GFP-like reporter to the independent
natural protein set published by Somermeyer and colleagues. Look for exact
matches even if the source record name uses an alternative spelling. Report
the number of corroborated proteins and their source genera. Non-GFP reporters
are catalogue-labeled leads for future out-of-family verification. Source
dates and release hashes are available, while exhaustive pretraining exposure
cannot be established. [Protein publication and data](https://doi.org/10.7554/eLife.75842)

Use the same positive-reference policy and identity metrics as in the
existing distant-reporter code. Before alignment, fix three natural reference
groups from the authors' release by taxon. Exclude all anthozoan proteins from
these positive reference groups. Report the reference groups individually and
together. Count independently corroborated reporters below 25% global
identity whose closest same-label local alignment also fails to cover at least
70% of both proteins. The individual reference groups are sensitivity
conditions and cannot be combined as independent sets.

Run the existing minimal and Pyrodigal candidate methods on raw archived DNA
only. Compare whether a candidate covers the documented CDS on the correct
strand at budgets of one, two and four, with two primary. Record exact peptide
matches to the validated translation separately. This tests reachability on
fresh archived records. Family rank, abstention, and model-assisted accuracy
require a later study with an independently frozen reference and comparison
panel and fresh inference.
