cp files/manual-overrides.yaml .
cp files/gather_manual_templates.py .
cp files/generate_pdbs.py .

ensembler init
ensembler gather_targets --gather_from uniprot --query 'accession:Q9NQR1' --uniprot_domain_regex SET
ensembler gather_templates --gather_from uniprot --query 'accession:Q9NQR1'

python generate_pdbs.py
python gather_manual_templates.py

ensembler align

cp files/alignments/1ZKK_apo.pir models/KMT5A_HUMAN_D0/KMT5A_HUMAN_1ZKK_apo/alignment.pir
cp files/alignments/4IJ8_apo.pir models/KMT5A_HUMAN_D0/KMT5A_HUMAN_4IJ8_apo/alignment.pir
cp files/alignments/Inhibitor_apo.pir models/KMT5A_HUMAN_D0/KMT5A_HUMAN_Inhibitor_apo/alignment.pir
cp files/alignments/Inhibitor_4IJ8.pir models/KMT5A_HUMAN_D0/KMT5A_HUMAN_Inhibitor_4IJ8/alignment.pir

ensembler build_models
ensembler cluster --cutoff 0
ensembler refine_implicit
ensembler solvate

cp files/nwaters-use.txt models/KMT5A_HUMAN_D0
ensembler refine_explicit --simlength 5*nanoseconds
ensembler package_models --package_for FAH
