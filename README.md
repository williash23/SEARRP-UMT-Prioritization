# SEARRP-UMT-Prioritization

This repository holds scripts used to prioritize multiple objectives in a two step process across Sabah, Malaysia. All work was done using privately shared data generously provided by the Sabah Forestry Department and their collaborators, including Forest Trends and the Carnegie Airborn Observatory (CAO). 

The objectives were largely grouped to promote the persistence of biodiversity and cover seven input categories:
1. Vertebrate species ranges (binary occurence/non-occurence spatial range maps downloaded from the IUCN Spatial Data Download website)
2. Invertebrate species ranges (restricted to butterflies; outputs from collaboratively generated processed Species Distribution Models; probability of occurence continuous value from 0-1)
3. Plant species ranges (outputs from collaboratively generated processed Species Distribution Models; probability of occurence continuous value from 0-1; also includes many extremely rare species representated by <10 localities)
4. Forest types (i.e., forest foramtions; outputs from Forest Trends; binary occurence/non-occurence of forest type)
5. Aboveground Carbon Density (from CAO; pixel values indicate amount of carbon stored)
6. Elevational conectivity (output from Condatis connectivity analysis; values are continuous)
7. Dispersal corridors (output from a dispersal movement simulation and corridor prioritization anlaysis; binary prioritized corridor/not)


Processing of data inputs occured before the prioritization process and was done to (briefly), make sure all raw feature inputs aligned in resolution and extent, so that they could be processed together in a sinlge conservation problem. Study area boundaries and locations of existing forest and protected areas was also done before the prioritization. 


The scripts here cover the main two-step prioritization analysis (which occured in two steps), which was accomplished using the 'prioritizr' package (Hanson et al., https://prioritizr.net/index.html), and the process used to assess coverage and achievement of the solution (i.e., the prioritized area) for each raw input feature.
