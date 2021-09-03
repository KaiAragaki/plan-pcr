# plan-pcr <img src='plan-pcr.png' align="right" height="138" />

# Statement of Need

Planning a PCR experiment is typically a multi-context affair:
1. Think about an experiment, then fiddle around in Excel for a plate layout
2. Calculate the amount of reagents you'll need (either 📝 or on Excel, which may or may not be the same Excel sheet)
3. Figure out how to dilute your RNA (again, this could be Excel or 📝)

This sucks for a few reasons:
1. Context switching requires cognitive overhead and invites error between switches
2. Manual calculations invite error
3. Multiple files are harder to manage
4. Pen and paper is difficult to manage - it's unsearchable, and it's difficult to link to digital documents (if it even gets saved)

This process almost always requires no creativity, and therefore can be solved algorithmically, freeing our brains to work on more stimulating things.

# Description

`plan-pcr` is a `shiny` app that requires only RNA concentration levels of samples and the number of primers you'd like to use. From there, it:
- Performs dilution calculations (including if an interim dilution, or 'dilution factor', needs to be made)
- Calculates the amount of mastermix needed for each primer
- Guides plate layout

Additionally, it also prevents the user from performing misguided experiments, by: 
- Automaticaly adding a non-targeting control for each primer
- Warning if only one primer is used (implying no endogenous control is being used)
- Throwing an error if the experiment requires more wells than the current plating configuration allows for

Together, this app forms a flexible yet robust method for PCR experiment planning quickly and reproducibly.

# Usage

First, go to the [shiny app webpage](https://kai-a.shinyapps.io/plan-pcr/) for plan-pcr.

Next, upload a dataset: 
* The dataset should be a `.tsv`, `.csv`, or `.xls(x)` file. 
* If there is one column of data, it is assumed to be RNA concentrations. 
* If there are 2+ columns, the first is assumed to be the sample names, the second is assumed to be RNA concentrations.
* The file is assumed to have column names

Adjust the parameters as your experiment requires. Typically, excecution of the experiment should follow a left-to-right progression through the tabs.

# TODO

- [x] Support more file types (.csv, .xls, .xlsx)
- [ ] Report generation
- [ ] Deal with long sample/primer names more gracefully (likely trucation or excision using `stringr`)
- [x] Better aesthetics
- [x] Code refactoring
- [ ] More distinct colors, particularly when the number of primers get very high (low priority)

