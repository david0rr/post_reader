### Example run of post_reader
---

Using default parameters and specifying a single species to label

```Shell
python post_reader.py --path codeml_output/F6935_Domain1/ --AA_fasta original_fasta/F6935_Domain1.fasta.vespa.mft --species DANRE05082
```

Using default parameters and specifying multiple species

```Shell
python post_reader.py --path codeml_output/F6935_Domain1/ --AA_fasta original_fasta/F6935_Domain1.fasta.vespa.mft --species_list DANRE05082 --species_list POEFO05164 --species_list TETNG00360
```

Using manual parameters and specifying multiple species

```Shell
python post_reader.py --path codeml_output/F6935_Domain1/ --AA_fasta original_fasta/F6935_Domain1.fasta.vespa.mft --species_list DANRE05082 --species_list POEFO05164 --species_list TETNG00360 --b_range 10 --dist 0.1 --gaps_cutoff 0.75 --c_range 5
```