# R Codes for Yurui Jiang's MSc Statistics Dissertation

This repository contains the R code used in Yurui Jiang's MSc Statistics dissertation:

**The Geometry of Political Ideology: Comparing Euclidean and Non-Euclidean Dimension Reduction and Latent Space Models**

The study investigates the geometric structure of political ideology using U.S. congressional roll-call voting data. It compares Euclidean and non-Euclidean dimension-reduction techniques with latent space models and evaluates how well different methods preserve the structure of voting behaviour.

## Repository Structure
The analysis is organised into five main R scripts.

### 1. `voteview dataset description.R`
This script reads, cleans, and describes the **119th U.S. House Voteview dataset**. It performs the main data-preparation and descriptive-analysis steps, including:
* filtering legislators and votes according to data-completeness criteria;
* examining party composition;
* summarising vote distributions; and
* calculating pairwise voting disagreement between members.

The cleaned voting data generated through this process provide the basis for the main empirical analyses.

### 2. `voteview unified method comparison.R`
This script contains the **main empirical analysis** of the dissertation.
It applies the full set of dimension-reduction and latent-space methods to a common cleaned sample of Members' votes, allowing the methods to be compared under the same data conditions.
The script produces the principal method-comparison results used in the dissertation.

### 3. `voteview graph.R`
This script generates graphical representations of the congressional voting data and the fitted low-dimensional configurations obtained from the different methods.
These visualisations are used to interpret the estimated ideological geometry and compare the structures recovered by alternative modelling approaches.

### 4. `voteview policy profile comparison.R`
This script conducts a supplementary **robustness analysis** based on policy-area voting profiles.
Individual votes are aggregated according to **Congressional Research Service (CRS) policy areas**, and dimension-reduction methods are then compared in terms of their ability to preserve distances between legislators in the resulting policy-profile space.
This analysis examines whether the main conclusions remain robust when political behaviour is represented at the policy-area level rather than through individual roll-call votes.

### 5. `simulation.R`
This script implements the simulation study used to evaluate method performance under alternative latent ideological structures.
It contains two main simulation settings:
* a **cyclic-policy structure**; and
* a **multi-bloc structure**.

These simulations are used to investigate how multidimensional scaling (MDS) and latent space models perform when the underlying political geometry departs from a simple Euclidean ideological continuum.

## Data
The empirical analysis uses roll-call voting data from the **119th U.S. House of Representatives**, obtained from Voteview.
The scripts include the data-cleaning and preprocessing procedures required for the analyses described in the dissertation.

## Reproducibility
The scripts in this repository correspond to the empirical analyses, robustness checks, simulations, and graphical results reported in the dissertation.
For reproducibility, users should run the relevant data-preparation procedures before executing scripts that depend on the cleaned Voteview data.

## Author

**Yurui Jiang**

MSc Statistics Dissertation
