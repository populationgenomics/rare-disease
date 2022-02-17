PAP_DATE="2021-09-03"

analysis-runner \
  --dataset acute-care \
  --description "run reanalysis draft" \
  -o "reanalysis/${PAP_DATE}" \
  --access-level test \
  reanalysis/reanalysis_wrapper.py \
    --conf gs://cpg-acute-care-test/reanalysis/reanalysis_conf.json \
    --matrix gs://cpg-acute-care-main/mt/acute-care.mt \
    --pap_date ${PAP_DATE} \
    --ped gs://cpg-acute-care-test/reanalysis/cpg_id_pedigree.json
