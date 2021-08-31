# plan-pcr

# Usage

First, go to the [shiny app webpage](https://kai-a.shinyapps.io/plan-pcr/) for plan-pcr.

Next, upload a dataset: 
* The dataset should be a `.tsv` file. 
* If there is one column of data, it is assumed to be RNA concentrations. 
* If there are 2+ columns, the first is assumed to be the sample names, the second is assumed to be RNA concentrations.
* The file is assumed to have column names

Adjust the parameters as your experiment requires. Typically, excecution of the experiment should follow a left-to-right progression through the tabs.

# TODO

- [ ] Support more file types (.csv, .xls, .xlsx)
- [ ] Report generation
- [ ] Deal with long sample/primer names more gracefully (likely trucation or excision using `stringr`)
- [ ] Better aesthetics
- [ ] Code refactoring
- [ ] More distinct colors, particularly when the number of primers get very high (low priority)

