# Combine Multiple LIBS Spectra by Stacking Shots

Concatenates shots from multiple spectra. All spectra must share the
same wavelength axis.

## Usage

``` r
# S3 method for class 'libs_spectrum'
c(...)
```

## Arguments

- ...:

  Two or more `libs_spectrum` objects.

## Value

A new `libs_spectrum` with combined shots.
