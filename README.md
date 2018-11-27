### POST READER README
---

POST_READER is a tool to be utilised after the analysis of postive site dectection using VESPA or VESPA-slim.

Post-reader will look to sift throught the outputs of any VESPA analysis and will remove any controversial detected sites through the analysis of unique site substitutions in desired species as well as sequence conservation at surrounding regions in addition to removal of sites flanked by regions which consist primarily of gaps.


### ROADMAP
---

Commands:

unique_sub [u] ✅

blanks [b] ✅

conservation [c] ✅

report [create-report] ❌

species ✅

multiple species ❌

### UPDATES
---

Abilty to output dataframe report detailing which sites were filtered/passed at each step

Function to apply commands to particular species or lineage

