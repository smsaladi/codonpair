## Regression test datasets

By using genes in the *E. coli* K-12 proteome that were identified to have
their C-terminus in the cytoplasm yields a good set to calculate tAI for.

```bash
cat test/data/Daley_gfp.fna | parallel  -N1 --recstart '>' --pipe -k './cps.pl -q -p data/ec_de3_ref_cps.tbd - ' > test/data/Daley_gfp.fna.ec_cpb
cat test/data/ecolik12.ffn | parallel  -N1 --recstart '>' --pipe -k './cps.pl -q -p data/ec_de3_cps.tbd - ' > test/data/ecolik12.ffn.ec_cpb
cat test/data/ecolik12.ffn | parallel  -N1 --recstart '>' --pipe -k './cps.pl -q -p data/sp_ref_cps.tbd - ' > test/data/ecolik12.ffn.sp_cpb
cat test/data/Daley_gfp.fna | parallel  -N1 --recstart '>' --pipe -k './cps.pl -q -p data/sp_ref_cps.tbd - ' > test/data/Daley_gfp.fna.sp_cpb
```
