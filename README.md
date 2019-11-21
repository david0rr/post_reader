### POST READER README
---

Post_reader is a tool to be utilised after the analysis of postive site detection using ~~VESPA or VESPA-slim~~ [Vespasian](https://github.com/bede/vespasian).

Post_reader will look to sift throught the outputs of any ~~VESPA~~ [Vespasian](https://github.com/bede/vespasian) analysis and will remove any controversial detected sites through three analyses: 
  1. unique site substitutions in desired species
  2. sequence conservation at surrounding regions 
  3. removal of sites flanked by regions which consist primarily of gaps.


### ROADMAP
---

Commands:

unique_sub ✅

blanks ✅

conservation ✅

report [create-report] :white_check_mark:

species ✅

multiple species :white_check_mark:

### UPDATES
---

* Code updated for vespasian outputs:
  1. Load original aligned AA fasta file into memory as a dictionary
  2. Parse vespasian summary.tsv file to find positively selected sites

* Ability to run filter tests on multiple species at once

* Use argparse as the arguments system to run the script - input commands are now different

* Output sites that passed and failed to separate tsv files

### IN PROGRESS
---

Run `unique_selected_sites` function only - right now it runs all filter functions (use `cli.py` file to do this)

Update the output pass function

Allow script to load `branches.yaml` file to get the list of species rather than using --species_list multiple times

Create alignment viewer image with sites labelled in command line using [Jalview](http://www.jalview.org/help/html/features/commandline.html)

### USAGE
---

The directory `example/` above shows an example on how to run post_reader

