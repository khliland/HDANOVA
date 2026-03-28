# Update a Model without Factor

Perform a model update while removing a chosen factor. Hierarchical
corresponds to type "II" sum-of-squares, i.e., obeying marginality,
while non-hierarchical corresponds to type "III" sum-of-squares.

## Usage

``` r
update_without_factor(model, fac, hierarchical = TRUE)
```

## Arguments

- model:

  `model` object to update.

- fac:

  `character` factor to remove.

- hierarchical:

  `logical` obey hierarchy when removing factor (default = TRUE).

## Value

An updated model object is returned. If the supplied model is of type
`lmerMod` and no random effects are left, the model is automatically
converted to a linear model before updating.
