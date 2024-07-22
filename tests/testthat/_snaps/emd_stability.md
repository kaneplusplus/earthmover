# `emd_stability()` works for matrices

    Code
      emd_stability(X, Y)
    Output
      An object of class "earthmover_stability"
      Slot "jack_dists":
      # A tibble: 8 x 6
        var     row dist_omnibus dist_within p_omnibus p_within
        <chr> <int>        <dbl>       <dbl>     <dbl>    <dbl>
      1 x         1         3.73       1.16      0.260    0.290
      2 x         2         3.73       0.877     0.260    0.290
      3 x         3         3.07       1.17      0.260    0.290
      4 y         1         3.38       0.717     0.478    0.292
      5 y         2         3.45       0.660     0.489    0.292
      6 y         3         3.65       1.16      0.350    0.292
      7 y         4         3.47       1.05      0.489    0.292
      8 y         5         3.18       0.788     0.260    0.324
      
      Slot "p":
      [1] 2
      

# `emd_stability()` works for matrices 2

    Code
      emd_stability(X, X)
    Output
      An object of class "earthmover_stability"
      Slot "jack_dists":
      # A tibble: 6 x 6
        var     row dist_omnibus dist_within p_omnibus p_within
        <chr> <int>        <dbl>       <dbl>     <dbl>    <dbl>
      1 x         1        1.16        1.16      0.268    0.290
      2 x         2        0.877       0.877     0.268    0.290
      3 x         3        1.17        1.17      0.268    0.290
      4 y         1        1.16        1.16      0.268    0.290
      5 y         2        0.877       0.877     0.268    0.290
      6 y         3        1.17        1.17      0.268    0.290
      
      Slot "p":
      [1] 2
      

# `emd_stability()` works for data.frames

    Code
      emd_stability(iris[1:10, ], iris[11:20, ])
    Output
      An object of class "earthmover_stability"
      Slot "jack_dists":
      # A tibble: 20 x 6
         var     row dist_omnibus dist_within p_omnibus p_within
         <chr> <int>        <dbl>       <dbl>     <dbl>    <dbl>
       1 x         1        0.678       0.174     0.369   0.387 
       2 x         2        0.603       0.176     0.369   0.387 
       3 x         3        0.612       0.163     0.369   0.387 
       4 x         4        0.615       0.180     0.369   0.387 
       5 x         5        0.673       0.178     0.369   0.387 
       6 x         6        0.731       0.318     0.323   0.0544
       7 x         7        0.641       0.171     0.469   0.387 
       8 x         8        0.659       0.152     0.451   0.387 
       9 x         9        0.619       0.249     0.369   0.387 
      10 x        10        0.610       0.163     0.369   0.387 
      11 y         1        0.642       0.228     0.469   0.385 
      12 y         2        0.694       0.276     0.369   0.428 
      13 y         3        0.710       0.339     0.323   0.428 
      14 y         4        0.709       0.447     0.323   0.133 
      15 y         5        0.570       0.324     0.323   0.428 
      16 y         6        0.564       0.373     0.323   0.385 
      17 y         7        0.619       0.249     0.369   0.422 
      18 y         8        0.679       0.226     0.369   0.385 
      19 y         9        0.596       0.292     0.369   0.428 
      20 y        10        0.650       0.227     0.469   0.385 
      
      Slot "p":
      [1] 2
      

# `emd_stability()` works for data.frames 2

    Code
      emd_stability(iris[1:10, ], iris[1:10, ])
    Output
      An object of class "earthmover_stability"
      Slot "jack_dists":
      # A tibble: 20 x 6
         var     row dist_omnibus dist_within p_omnibus p_within
         <chr> <int>        <dbl>       <dbl>     <dbl>    <dbl>
       1 x         1        0.174       0.174    0.384    0.387 
       2 x         2        0.176       0.176    0.384    0.387 
       3 x         3        0.163       0.163    0.384    0.387 
       4 x         4        0.180       0.180    0.384    0.387 
       5 x         5        0.178       0.178    0.384    0.387 
       6 x         6        0.318       0.318    0.0445   0.0544
       7 x         7        0.171       0.171    0.384    0.387 
       8 x         8        0.152       0.152    0.384    0.387 
       9 x         9        0.249       0.249    0.384    0.387 
      10 x        10        0.163       0.163    0.384    0.387 
      11 y         1        0.174       0.174    0.384    0.387 
      12 y         2        0.176       0.176    0.384    0.387 
      13 y         3        0.163       0.163    0.384    0.387 
      14 y         4        0.180       0.180    0.384    0.387 
      15 y         5        0.178       0.178    0.384    0.387 
      16 y         6        0.318       0.318    0.0445   0.0544
      17 y         7        0.171       0.171    0.384    0.387 
      18 y         8        0.152       0.152    0.384    0.387 
      19 y         9        0.249       0.249    0.384    0.387 
      20 y        10        0.163       0.163    0.384    0.387 
      
      Slot "p":
      [1] 2
      

