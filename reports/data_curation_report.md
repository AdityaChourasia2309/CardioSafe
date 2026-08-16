# hERG Dataset Curation Report

## Sources
- ChEMBL (target CHEMBL240): 20073 raw bioactivity records with valid pIC50
- TDC hERG: attempted, unavailable (Harvard Dataverse under maintenance at time of pull)

## Curation Steps
1. Canonicalized all SMILES using RDKit (salt stripping, largest-fragment selection)
2. Dropped 0 unparseable SMILES
3. Converted ChEMBL IC50 (nM) to pIC50 = 9 - log10(IC50_nM)
4. Derived binary blocker label: pIC50 >= 5.0 (IC50 <= 10 uM) = blocker
5. De-duplicated on canonical SMILES; conflicting values resolved via median (pIC50) / majority vote (label)

## Final Dataset
- Unique compounds: 16159
- Blockers (label=1): 9000
- Non-blockers (label=0): 7159

## Schema
compound_id, smiles, pIC50, label, n_records, sources
