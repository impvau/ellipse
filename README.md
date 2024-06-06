
# Overview

Code repository for the publication: A new approximation for the perimeter of the ellipse.

| File                         | Details                                                                                          |
|------------------------------|--------------------------------------------------------------------------------------------------|
| [Main Notebook](ellipse.nb)   | Entry point for the functionality that re-generates all data files, tables, and plots            |
| [Training Dataset](data.m)   | Generates the training data points for use in the main notebook                                  |
| [Plot functionality](plot.m) | Generate plots according to scientific reports functionality                                     |
| [Ellipse functions](funcs.m) | Define a series of functions that are evaluated on the training data or plotting data            |
| [Supplementary Material](SupplementaryMaterial.pdf) | Complete PDF of supplementary material   |
| [Table Helper](MathaToLatex.xlsx) | Conversion tool between Mathematica output table and latex displayed table |

Please cite this work as

```
@inproceedings{Moscato:2023,
    author      = {Moscato, Pablo and Ciezak, Andrew},
    title       = {A new approximation for the perimeter of the ellipse},
    booktitle   = {},
    series      = {},
    year        = {2024},
    location    = {},
    pages       = {},
    numpages    = {},
    url         = {},
    doi         = {},
    isbn        = {},
    publisher   = {},
    address     = {},
}
```

## Setup

1. Open and run data.m, plot.m, funcs.m
2. Proceed to running ellipse.nb

Results will be output to the out/ directory

## Notes
- To produce the figure for plotting a = 0 to a = 80, open plot.m and comment/uncomment the appropriate sections, run the code in plots.m, then re-run the code to produce the specific plot
- Note that saving files after opening them in excel will cause truncation of the high precision to <20 decimla places
- We avoid precision issues in our table that uses scientific notation by copying plain text from Mathematica table into MathaToLatex.xlsx which generates our latex table

# Sykora Ellipse Approximations

Given the methods put forward in [Sykoras analysis](http://www.ebyte.it/library/docs/math05a/EllipsePerimeterApprox05.html) we re-define these functions within the within the [supplementary material](SupplementaryMaterial.pdf) and produce the comparsion plots that are define in this repository.