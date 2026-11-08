#!/bin/bash
# Pull segment info from ENA for the 19 studies with NA segments.
# LOGIN NODE ONLY (needs internet).  bash fetch_segments.sh
cd "$(dirname "$0")"
printf "study\trun\talias\ttitle\tdescription\n" > segments_lookup.tsv
echo ">> Rebollada_Merino_2020  (181 NA runs)  PRJNA643396"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA643396&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Rebollada_Merino_2020\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Ma_2024  (155 NA runs)  PRJNA1086343"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1086343&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Ma_2024\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Zhang_2026a  (144 NA runs)  PRJNA1223836"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1223836&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Zhang_2026a\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Bernard_2024  (57 NA runs)  PRJEB62837"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB62837&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Bernard_2024\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Xu_2026  (50 NA runs)  CRA019535"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=CRA019535&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Xu_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Chen_2024  (40 NA runs)  PRJNA1062743"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1062743&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Chen_2024\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Jin_2026  (36 NA runs)  PRJNA997628"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA997628&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Jin_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Cao_2020  (32 NA runs)  PRJNA682703"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA682703&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Cao_2020\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Obianwuna_2026  (32 NA runs)  PRJNA1202053"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1202053&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Obianwuna_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Tang_2026  (30 NA runs)  PRJNA1270426"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1270426&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Tang_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Xu_2023  (24 NA runs)  PRJNA830845"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA830845&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Xu_2023\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Ji_2026  (24 NA runs)  NMDC40084412"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=NMDC40084412&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Ji_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Zhang_2023  (18 NA runs)  PRJNA904521"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA904521&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Zhang_2023\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Zhou_2026  (18 NA runs)  PRJNA1413382"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1413382&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Zhou_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Nong_2026  (12 NA runs)  PRJNA1442398"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1442398&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Nong_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Zhang_2026b  (12 NA runs)  PRJNA1290854"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1290854&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Zhang_2026b\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Du_2026  (12 NA runs)  CRA024363"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=CRA024363&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Du_2026\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Biasato_2018  (10 NA runs)  SRP125741"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP125741&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Biasato_2018\t/" >> segments_lookup.tsv
sleep 0.3
echo ">> Wassie_2021  (6 NA runs)  PRJNA765947"
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA765947&result=read_run&fields=run_accession,sample_alias,sample_title,description" | tail -n +2 | sed "s/^/Wassie_2021\t/" >> segments_lookup.tsv
sleep 0.3
echo; echo "rows: $(wc -l < segments_lookup.tsv)"
echo "--- distinct descriptions per study ---"
awk -F'\t' 'NR>1{print $1"\t"$5}' segments_lookup.tsv | sort -u
