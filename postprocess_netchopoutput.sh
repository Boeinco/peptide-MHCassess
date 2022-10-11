sed '/NA/d' netchop_out.txt > netchop_out_parsed.txt
awk -F'\t' 'NR>1{$0=$0"\t"NR-1} 1' netchop_out_parsed.txt > netchop_numbered.txt
