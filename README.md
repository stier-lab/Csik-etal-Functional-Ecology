# Csik-etal-Functional-Ecology

This repository contains the code and data for reproducing Csik *et al.* **in prep** (to be submitted to Functional Ecology): **"TITLE HERE"** (**URL HERE**)

The DOI for this code and data repository is managed through Dryad with DOI number **NUMBER HERE** (**LINK HERE**).

> **Title:** 

> **Authors:** Samantha R. Csik, Bartholomew P. DiFiore, Krista Kraskura, Emily A. Hardison, Joseph S. Curtis, Erika J. Eliason, Christopher L. Jerde, Adrian C. Stier

> **Abstract:**  

### Repository Structure

```
Csik-etal-Functional-Ecology
  |_ code
  |_ data
    |_ correlations
      |_ outputs
    |_ foraging
      |_ outputs
      |_ raw
    |_ heart_rate
      |_ inventory_files
      |_ outputs
      |_ raw
    |_ metabolism
      |_ outputs
    |_ metadata
  |_ figures
    |_ main_text
    |_ other
    |_ supplemental
  |_ media
```

### Software

These analyses were performed in R (version 3.6.3).

### Code 

Code is meant to be run in the order designated by the file name numbers (e.g. `0_libraries.R` should be run first, `15_in_text_summary_stats.Rmd` should be run last). The only exception to this is the `FR_JAGS` series (number 11) -- only `11c_FR_JAGS_plotting.R` needs to be run (both `11a_FR_JAGS_functions.R` and `11b_FR_JAGS_model.R` are sourced into this file and do not need to be opened separately to run analyses).

At the start of each script and RMarkdown file, you will find a brief **summary** of the analyses to follow, an **outline** of code chunks/subsections, any **required packages**, and **required data** (along with their file paths for ease of locating).

### Data

The following three files are the primary data used for all analyses. Their metadata are detailed below: 

#### * `data/heart_rate/outputs/processed_heart_rates.csv`
* `date`: date/time of data point collection
* `temp`: acclimation temperature
* `lobster_id`: unique lobster identifier
* `hr`: heart rate (beats min<sup>-1</sup>) 
* `QI`: Quality index score ranging from 0 (good) to 3 (poor); NOTE: these data have been visually assessed for quality and therefore any hr with a QI = 1-3 was determined to still be an accurate measure.

#### * `data/metabolism/metabolic_traits.csv`
* `temp`: acclimation temperature
* `ID`: unique lobster identifier
* `BW`: body weight of lobster (g)
* `SMR`: standard metabolic rate (mg O<sub>2</sub> kg<sup>-1</sup> min<sup>-1</sup>), calculated as the 15th percentile MO<sub>2</sub> measurement
* `MMR`: maximum metabolic rate (mg O<sub>2</sub> kg<sup>-1</sup> min<sup>-1</sup>), calculated as the fastest rate of linear O<sub>2</sub> decline over a 60s inverval
* `AAS`: absolute aerobic scope (mg O<sub>2</sub> kg<sup>-1</sup> min<sup>-1</sup>); calculated as MMR-SMR
* `FAS`: factorial aerobic scope (unitless); calculated as MMR/SMR

#### * `data/metabolism/epoc.csv`
* `temp`: acclimation temperature
* `ID` : unique lobster identifier
* `spar` : the smoothing level applied to the recovery profile (see function [`smooth.spline`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/smooth.spline) parameters)
* `EPOC_full` : total EPOC (i.e. recovery costs; mg O<sub>2</sub> kg<sup>-1</sup>); represents the integrated area under the smoothed line at a particular spar
* `end_EPOC_min` : total recovery time (min); the time it takes for for the oxygen consumption rate (MO<sub>2</sub>; the smoothed line) to reach SMR
* `SMR` : standard metabolic rate (mg O<sub>2</sub> kg<sup>-1</sup> min<sup>-1</sup>), calculated as the 15th percentile MO<sub>2</sub> measurement
* `SMR_integrated_full` : cost of SMR (area integrated from SMR to 0 mg O<sub>2</sub> kg<sup>-1</sup> min<sup>-1</sup>)
* `SMR_threshold` : the value used to calculate EPOC; here, we consider an individual recovered when oxygen consumption rates fall within 20% of SMR
* `EPOC_1hr` : EPOC from the start of the measurement (i.e. MMR) until 1 hr post exercise (up to EPOC_5hr; units = mg O<sub>2</sub> kg<sup>-1</sup>)
* `MO2_1hr` : oxygen consumption rate (mg O<sub>2</sub> min<sup>-1</sup>) at 1 hour post exercise (i.e. MMR; up to MO2_5hr)
* `end_EPOC_mmr` : recovery time (min) to reach 50% of MMR
* `EPOC_mmr` : EPOC at 50% of MMR
* `MO2_mmr` : the oxygen consumption rate at 50% of MMR
* `MMR` : maximum metabolic rate (mg O<sub>2</sub> kg<sup>-1</sup> min<sup>-1</sup>), calculated as the fastest rate of linear O<sub>2</sub> decline over a 60s inverval

#### * `data/foraging/raw/foraging_assay_data.csv`
* `date`: date the foraging asssay began on (each assay ran for 24h)
* `trial`: trial number; feeding assays were replicated three times per lobster at each of the five prey densities
* `lobster_id`: unique lobster identifier 
* `Initial`: number of mussels available at the start of the feeding assay
* `Killed`: number of mussels consumed by the end of the feeding assay
* `temp`: acclimation temperature

# Figures

* `/main_text`: contains pdf versions of Figures 1, 2, 3, 4, & 5 found in the main text of the manuscript
* `/supplementary`: contains pdf version of Figures S1, S2, S3, & S4 found in the supplementary text
* `/other`: contains plots of raw heart rate data (bpm & ecgs) used for manual data quality assessments

### Media

Experimental photos used in this README.

### Curious about what some of our experimental trials look like?
![Alt text](/media/respirometry.png?raw=true "A lobster inside an intermittent-flow respirometry chamber" )
![Alt text](/media/foraging.png?raw=true "Munching on mussels during a foraging trial")
![Alt text](/media/heart.png?raw=true "Heart rate loggers are implanted under the carapace, then the incision is sealed using dental wax")
