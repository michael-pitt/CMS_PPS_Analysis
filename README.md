# CMS-PPS-Analysis

A comprehensive framework for analysis in CMS that use PPS, covering Run 2 and Run 3 data. 

## Analysis Tools
1. **Proton Asymmetry**: Correlation studies between protons and the central system (SD events)
2. **CEP Jets**: Central Exclusive Production of Jets ($M_{jj}$ vs $M_{pp}$ matching).
3. **CEP Muons**: Central Exclusive Production of Muons in High-PU runs (Run2 and Run3 data).

## Repository Structure
- `PostProcessing/`: NanoAODTools modules for physics object enrichment.
- `Common/`: Era-dependent constants and calibration managers.
- `Analysis/`: RDataFrame-based plotting and Poisson fitting logic.

## Setup
```bash
cmsrel CMSSW_16_0_4
cd CMSSW_16_0_4/src
cmsenv
git clone [https://github.com/michael-pitt/CMS-PPS-Analysis.git](https://github.com/michael-pitt/CMS-PPS-Analysis.git)
