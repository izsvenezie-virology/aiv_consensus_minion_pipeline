#USAGE: bash estrai_best_hit_influenza_v2.0.sh title num
#ESEMPIO: bash estrai_best_hit_influenza_v2.0.sh 17VIR4455-2_subset0_mod.fil.tab.entry_title 10
#title: file .entry_title prodotto con scegli_ref_vX.0.sh
#num: numero di righe da stampare; se settato su 0, stampa solo la prima riga e formatta in modo diverso l'output

if [ ${2} -gt 0 ]
then
  echo "segment 1 PB2:"
  cat ${1} | grep "segment_1" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 2 PB1:"
  cat ${1} | grep "segment_2" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 3 PA:"
  cat ${1} | grep "segment_3" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 4 HA:"
  cat ${1} | grep "segment_4" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 5 NP:"
  cat ${1} | grep "segment_5" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 6 NA:"
  cat ${1} | grep "segment_6" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 7 MP:"
  cat ${1} | grep "segment_7" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -e "\nsegment 8 NS:"
  cat ${1} | grep "segment_8" | grep "_complete_cds" | grep -v "_partial_cds" | head -n ${2}
  echo -ne "\n"
else
  special="NON-SENSE_WORD_LIKE_ASH_NAZG_DURBATULUK,_ASH_NAZG_GIMBATUL,_ASH_NAZG_THRAKATULUK,_AGH_BURZUM-ISHI_KRIMPATUL"
  echo -ne "PB2\t"
  cat ${1} | grep "segment_1" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "PB1\t"
  cat ${1} | grep "segment_2" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "PA\t"
  cat ${1} | grep "segment_3" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "HA\t"
  cat ${1} | grep "segment_4" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "NP\t"
  cat ${1} | grep "segment_5" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "NA\t"
  cat ${1} | grep "segment_6" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "MP\t"
  cat ${1} | grep "segment_7" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
  echo -ne "NS\t"
  cat ${1} | grep "segment_8" | grep "_complete_cds" | grep -v "_partial_cds" | head -n 1 | cut -f 2 | awk -v var="${special}" 'BEGIN {getline; if ($0) {print $0} else {print var}}'
fi
