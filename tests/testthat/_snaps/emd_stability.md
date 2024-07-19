# `emd_stability()` works for matrices

    Code
      emd_stability(X, Y)
    Output
      An object of class "earthmover_stability"
      Slot "emds":
      # A tibble: 8 x 6
        var     row dist_across dist_within p_across p_within
        <chr> <int>       <dbl>       <dbl>    <dbl>    <dbl>
      1 x         1        3.73       1.16    0.108     0.852
      2 x         2        3.73       0.877   0.102     0.839
      3 x         3        3.07       1.17    0.0538    0.372
      4 y         4        3.38       0.717   1         1    
      5 y         5        3.45       0.660   1         1    
      6 y         6        3.65       1.16    1         0.470
      7 y         7        3.47       1.05    1         1    
      8 y         8        3.18       0.788   0.968     0.372
      
      Slot "normality":
      [1] TRUE
      

# `emd_stability()` works for matrices 2

    Code
      emd_stability(X, X)
    Output
      An object of class "earthmover_stability"
      Slot "emds":
      # A tibble: 6 x 6
        var     row dist_across dist_within p_across p_within
        <chr> <int>       <dbl>       <dbl>    <dbl>    <dbl>
      1 x         1       1.16        1.16     0.868    0.868
      2 x         2       0.877       0.877    0.372    0.372
      3 x         3       1.17        1.17     0.824    0.824
      4 y         4       1.16        1.16     0.868    0.868
      5 y         5       0.877       0.877    0.372    0.372
      6 y         6       1.17        1.17     0.824    0.824
      
      Slot "normality":
      [1] TRUE
      

# `emd_stability()` works for data.frames

    Code
      emd_stability(iris[1:10, ], iris[11:20, ])
    Output
      An object of class "earthmover_stability"
      Slot "emds":
      # A tibble: 20 x 6
         var     row dist_across dist_within p_across p_within
         <chr> <int>       <dbl>       <dbl>    <dbl>    <dbl>
       1 x         1       0.678       0.174    1        1    
       2 x         2       0.603       0.176    1        1    
       3 x         3       0.612       0.163    1        1    
       4 x         4       0.615       0.180    1        1    
       5 x         5       0.673       0.178    1        1    
       6 x         6       0.731       0.318    0.541    0.169
       7 x         7       0.641       0.171    1        1    
       8 x         8       0.659       0.152    1        1    
       9 x         9       0.619       0.249    1        1    
      10 x        10       0.610       0.163    1        1    
      11 y        11       0.642       0.228    1        1    
      12 y        12       0.694       0.276    1        1    
      13 y        13       0.710       0.339    0.536    1    
      14 y        14       0.709       0.447    0.588    1    
      15 y        15       0.570       0.324    0.360    0.917
      16 y        16       0.564       0.373    0.252    0.738
      17 y        17       0.619       0.249    1        1    
      18 y        18       0.679       0.226    1        1    
      19 y        19       0.596       0.292    1        1    
      20 y        20       0.650       0.227    1        1    
      
      Slot "normality":
      [1] TRUE
      

# `emd_stability()` works for data.frames 2

    Code
      emd_stability(iris[1:10, ], iris[1:10, ])
    Output
      An object of class "earthmover_stability"
      Slot "emds":
      # A tibble: 20 x 6
         var     row dist_across dist_within p_across p_within
         <chr> <int>       <dbl>       <dbl>    <dbl>    <dbl>
       1 x         1       0.174       0.174   1        1     
       2 x         2       0.176       0.176   1        1     
       3 x         3       0.163       0.163   1        1     
       4 x         4       0.180       0.180   1        1     
       5 x         5       0.178       0.178   1        1     
       6 x         6       0.318       0.318   0.0730   0.0730
       7 x         7       0.171       0.171   1        1     
       8 x         8       0.152       0.152   1        1     
       9 x         9       0.249       0.249   1        1     
      10 x        10       0.163       0.163   1        1     
      11 y        11       0.174       0.174   1        1     
      12 y        12       0.176       0.176   1        1     
      13 y        13       0.163       0.163   1        1     
      14 y        14       0.180       0.180   1        1     
      15 y        15       0.178       0.178   1        1     
      16 y        16       0.318       0.318   0.0730   0.0730
      17 y        17       0.171       0.171   1        1     
      18 y        18       0.152       0.152   1        1     
      19 y        19       0.249       0.249   1        1     
      20 y        20       0.163       0.163   1        1     
      
      Slot "normality":
      [1] TRUE
      

