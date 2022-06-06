Tools to parse results from `samtools mpileup`

`scripts/baseparser.py` originally developed by Noushin Niknafs  

Example usage: 

```
cat example-data/SPP.test.mpileup.txt | python2 scripts/baseparser.py > example-data/SPP.test.mpileup.parsed.txt
Rscript parse-mpileup-results.R -o example-data/SPP_results.txt -s example-data/SPP.test.mpileup.parsed.txt -a SPP.test.txt
```


