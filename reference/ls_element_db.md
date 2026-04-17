# NIST Emission Line Database

Returns a curated subset of the NIST Atomic Spectra Database containing
analytically useful emission lines for LIBS analysis, with a focus on
biomedically relevant elements. Data compiled from the NIST Atomic
Spectra Database (<https://physics.nist.gov/asd>); values are for
reference use only and users requiring high-precision spectroscopy
should consult NIST ASD directly.

## Usage

``` r
ls_element_db(
  elements = NULL,
  range_nm = c(190, 900),
  min_aki = 0,
  ionization = c(1L, 2L),
  persistent_only = FALSE
)
```

## Source

NIST Atomic Spectra Database: <https://physics.nist.gov/asd>

## Arguments

- elements:

  Character vector. Restrict to specific element symbols (e.g.
  `c("Ca", "Fe")`). Default `NULL` (all).

- range_nm:

  Numeric vector of length 2. Wavelength range filter in nm. Default
  `c(190, 900)`.

- min_aki:

  Numeric. Minimum transition probability (10^8 s^-1). Default 0 (no
  filter).

- ionization:

  Integer vector of ionization states to include: 1 = neutral, 2 =
  singly ionized. Default `c(1L, 2L)`.

- persistent_only:

  Logical. Keep only persistent lines. Default `FALSE`.

## Value

A
[`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with columns `element`, `ionization`, `wavelength_nm`, `aki`, `ei_ev`,
`ek_ev`, `persistent`.

## Examples

``` r
head(ls_element_db())
#> # A tibble: 6 × 7
#>   element ionization wavelength_nm   aki ei_ev ek_ev persistent
#>   <chr>        <int>         <dbl> <dbl> <dbl> <dbl> <lgl>     
#> 1 Ca               1          423.  2.18  0     2.93 TRUE      
#> 2 Ca               1          240.  1.87  0     5.17 FALSE     
#> 3 Ca               1          272.  0.11  0     4.56 FALSE     
#> 4 Ca               1          430.  1.36  1.88  4.76 TRUE      
#> 5 Ca               1          444.  0.63  1.88  4.67 FALSE     
#> 6 Ca               1          445.  0.87  1.88  4.66 FALSE     
ls_element_db(elements = "Ca", persistent_only = TRUE)
#> # A tibble: 6 × 7
#>   element ionization wavelength_nm   aki ei_ev ek_ev persistent
#>   <chr>        <int>         <dbl> <dbl> <dbl> <dbl> <lgl>     
#> 1 Ca               1          423. 2.18   0     2.93 TRUE      
#> 2 Ca               1          430. 1.36   1.88  4.76 TRUE      
#> 3 Ca               2          393. 1.47   0     3.15 TRUE      
#> 4 Ca               2          397. 1.4    0     3.12 TRUE      
#> 5 Ca               2          854. 0.099  1.69  3.15 TRUE      
#> 6 Ca               2          866. 0.089  1.7   3.12 TRUE      
ls_element_db(elements = c("Na", "K"), range_nm = c(580, 780))
#> # A tibble: 6 × 7
#>   element ionization wavelength_nm   aki ei_ev ek_ev persistent
#>   <chr>        <int>         <dbl> <dbl> <dbl> <dbl> <lgl>     
#> 1 Na               1          589  0.616  0     2.10 TRUE      
#> 2 Na               1          590. 0.614  0     2.10 TRUE      
#> 3 K                1          766. 0.387  0     1.62 TRUE      
#> 4 K                1          770. 0.384  0     1.61 TRUE      
#> 5 K                1          694. 0.095  1.61  3.40 FALSE     
#> 6 K                1          691. 0.058  1.62  3.41 FALSE     
```
