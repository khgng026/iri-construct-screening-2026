# Offline recount of the optical and yeast controls

From this folder, run the following command with Python 3.10 or later.

```sh
python3 replay.py
```

The command recalculates all 45 training cutoffs and checks all reported counts against the two archived experiments. It prints nine rows, one per method and training allowance. To see the full paired tables, calibration results, returned-call counts and equal-series estimates, run `python3 replay.py --json`. The package only requires the Python standard library.

## How to read the output

| Column | Interpretation |
| --- | --- |
| Allowance | Maximum mean uncredited calls per complete-key *training* stock, with one as the primary setting. Held-out call burdens are measured separately. |
| All keys /11 | Number of the eleven optical stocks with every documented occurrence recovered. |
| All & no extras /11 | Number of those eleven with every documented occurrence recovered and every returned call credited. |
| Reporter /7 | Number of the seven optical stocks whose named reporter was recovered. Two have partial source keys. Five also appear in the eleven-stock inventory analysis. |
| Colonies /4 | Number of the four separate physical colony assemblies with both documented products recovered. They come from two designs. |
| Yeast /4 | Number of four related yeast parent plasmids whose single documented marker was recovered. Other calls have no adjudication in this endpoint. |

`blastx` is translated nucleotide–protein alignment. `minimal` and `gene_caller` supply different protein candidates to the same ESM-2 model. The experiment evaluates these conventional positive controls. Each reader can verify the original numerical result without fitting a model or downloading a database.

## What the data contain

[Evidence](evidence.json.gz) contains 21 pseudonymous records, four anonymous calibration series and three arrays of scored-call evidence. The 2,129 candidate calls include 335 from BLASTx, 1,731 from minimal translation, and 63 from gene calling. Within each record and arm, calls are ordered by their original score, with the original tie order preserved. Each candidate is encoded as three items.

| Item | Meaning |
| --- | --- |
| Score rank | Larger integers indicate stronger scores. Equal native scores have the same integer. Rank order is defined separately within each method. The lookup of previously reported numeric cutoffs is used only after rank-based calibration selects its operating point. |
| Eligible components | A bit in this integer identifies one anonymous component this call could cover. The source comparison required matching family and strand, at least 80% coverage of the documented occurrence, and frame agreement where both locations have a single contiguous segment. Component numbering is local to each record. |
| Competing calls | Indices of earlier calls of the same family that satisfied the study’s reciprocal 80% location-overlap and applicable frame rule. A duplicate is suppressed only if a listed earlier call was retained. |

The portable program starts from these score ranks and relations. It recomputes duplicate suppression and enumerates possible one-to-one assignments to source components. For each held-out series it constructs the calibration grid from the other complete-key series and selects the lowest qualifying cutoff. The two partial keys and the four external markers use the cutoff from all eleven training stocks. The colony checks use the training fold that excludes their design series. All paired comparisons use the same records across methods.

The projection omits biological sequences, stock identifiers, reference family names and genomic coordinates. Confirming the projected relations against physical coordinates or repeating the original model inference requires the original files and model runtime. The original numerical and prediction file digests, plus those of the packaged evidence and expected results, are recorded in the [provenance file](provenance.json). [Expected aggregate results](expected.json) are opened after reconstruction for comparison.

Both command forms were tested from a clean copy with the project files and network denied to the replay process. The [test receipt](clean-test.json) lists the artifact digests, checks and runtime.
