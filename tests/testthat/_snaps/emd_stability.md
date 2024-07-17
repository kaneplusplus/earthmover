# `emd_stability()` works for matrices

    Code
      emd_stability(X, Y)
    Output
      # A tibble: 8 x 6
        var     row    diff p_overall p_within p_across
        <chr> <int>   <dbl>     <dbl>    <dbl>    <dbl>
      1 x         1 -0.326      0.143     0.5     0    
      2 x         2 -0.331      0         0       0    
      3 x         3  0.328      0         0       0    
      4 y         1 -0.0243     0.429     0.25    0.333
      5 y         2  0.0480     0.429     0.5     0.333
      6 y         3  0.246      0.143     0       0.333
      7 y         4  0.0707     0.286     0.25    0.333
      8 y         5 -0.218      0.286     0       0.333

# `emd_stability()` works for matrices 2

    Code
      emd_stability(X, X)
    Output
      # A tibble: 6 x 6
        var     row   diff p_overall p_within p_across
        <chr> <int>  <dbl>     <dbl>    <dbl>    <dbl>
      1 x         1 -1.16        0.2      0.5        0
      2 x         2 -0.877       0.4      0          0
      3 x         3 -1.17        0        0          0
      4 y         1  1.16        0.2      0.5        0
      5 y         2  0.877       0.4      0          0
      6 y         3  1.17        0        0          0

# `emd_stability()` works for data.frames

    Code
      emd_stability(iris[1:10, ], iris[11:20, ])
    Output
      # A tibble: 20 x 6
         var     row     diff p_overall p_within p_across
         <chr> <int>    <dbl>     <dbl>    <dbl>    <dbl>
       1 x         1 -0.0540     0.105     0.111      0.1
       2 x         2  0.0207     0.263     0          0.5
       3 x         3  0.0117     0.421     0.222      0.4
       4 x         4  0.00825    0.474     0.333      0.4
       5 x         5 -0.0495     0.211     0.222      0.2
       6 x         6 -0.108      0         0          0  
       7 x         7 -0.0177     0.368     0.444      0.3
       8 x         8 -0.0355     0.263     0.333      0.2
       9 x         9  0.00519    0.474     0.444      0.4
      10 x        10  0.0141     0.368     0.111      0.4
      11 y         1  0.0180     0.316     0.444      0.1
      12 y         2  0.0699     0.105     0.222      0  
      13 y         3  0.0867     0         0          0  
      14 y         4  0.0848     0.0526    0.111      0  
      15 y         5 -0.0535     0.158     0.111      0.2
      16 y         6 -0.0600     0.0526    0          0.1
      17 y         7 -0.00438    0.421     0.333      0.5
      18 y         8  0.0554     0.158     0.333      0  
      19 y         9 -0.0280     0.316     0.222      0.4
      20 y        10  0.0264     0.211     0.444      0  

# `emd_stability()` works for data.frames 2

    Code
      emd_stability(iris[1:10, ], iris[1:10, ])
    Output
      # A tibble: 20 x 6
         var     row   diff p_overall p_within p_across
         <chr> <int>  <dbl>     <dbl>    <dbl>    <dbl>
       1 x         1 -0.174    0.263     0.444        0
       2 x         2 -0.176    0.211     0.444        0
       3 x         3 -0.163    0.421     0.111        0
       4 x         4 -0.180    0.105     0.222        0
       5 x         5 -0.178    0.158     0.333        0
       6 x         6 -0.318    0         0            0
       7 x         7 -0.171    0.316     0.333        0
       8 x         8 -0.152    0.474     0            0
       9 x         9 -0.249    0.0526    0.111        0
      10 x        10 -0.163    0.368     0.222        0
      11 y         1  0.174    0.263     0.444        0
      12 y         2  0.176    0.211     0.444        0
      13 y         3  0.163    0.421     0.111        0
      14 y         4  0.180    0.105     0.222        0
      15 y         5  0.178    0.158     0.333        0
      16 y         6  0.318    0         0            0
      17 y         7  0.171    0.316     0.333        0
      18 y         8  0.152    0.474     0            0
      19 y         9  0.249    0.0526    0.111        0
      20 y        10  0.163    0.368     0.222        0

