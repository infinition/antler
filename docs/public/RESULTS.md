# Validated public result

The digital transport protocol implements an Abelian logical phase primitive
in the correlated-hopping ladder.

The numerical deep-limit closure uses `D=4,6,8,10,12`, all at `dt=0.125`.
The fitted law is

```text
odd_slope(D) = -1 + 0.66493 D^(-2.28624)
```

with `R^2 = 0.999983`. Leakage is below `5e-5` at every retained point.
The minimum singular value exceeds `0.999975`.

The gate family commutes. This release establishes an Abelian phase primitive,
not non-Abelian braiding or topological quantum computation.
