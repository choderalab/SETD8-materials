cp files/manual-overrides.yaml .
cp files/gather_manual_templates_TDIY.py .
cp files/split_chains_TDIY.py .

ensembler init
ensembler gather_targets --gather_from uniprot --query 'accession:Q9NQR1' --uniprot_domain_regex SET
ensembler gather_templates --gather_from uniprot --query 'accession:Q9NQR1'

python split_chains_TDIY.py
python gather_manual_templates_TDIY.py

ensembler align
ensembler build_models
ensembler cluster --cutoff 0
ensembler refine_implicit
ensembler solvate

cp files/nwaters-use.txt models/KMT5A_HUMAN_D0
ensembler refine_explicit --simlength 5*nanoseconds
ensembler package_models --package_for FAH
