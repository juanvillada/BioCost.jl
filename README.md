# BioCost.jl

https://en.wikipedia.org/wiki/Biomolecule

- With the test proteome:

```{julia}
julia --threads 2 BioCost.jl data/test_genome/faa_folder output_test_2022827000.csv AA
julia --threads 2 BioCost.jl data/test_genome/fna_folder output_test_2022827000_dna.csv DNA
```

- With a MAG folder from our collection:

```{julia}
julia --threads 2 BioCost.jl /global/homes/j/jvillada/pdata/00_genome_files/02_mags_105k_300_per_chunk_genomes_faa/chunk_1 output_test.csv AA
```

```{julia}
julia --threads 2 BioCost.jl /global/homes/j/jvillada/pdata/00_genome_files/02_mags_105k_300_per_chunk_genomes_faa/chunk_301 output_test.csv AA
```
