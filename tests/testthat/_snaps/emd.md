# `emd()` works for matrices

    Code
      emd(X, Y)
    Output
      $dist
      [1] 3.400067
      
      $pairs
      # A tibble: 7 x 3
         from    to   mass
        <dbl> <dbl>  <dbl>
      1     1     1 0.133 
      2     1     4 0.2   
      3     2     1 0.0667
      4     2     2 0.0667
      5     2     5 0.200 
      6     3     2 0.133 
      7     3     3 0.2   
      
      attr(,"class")
      [1] "earthmover_dist"

---

    Code
      emd(X, X)
    Output
      $dist
      [1] 0
      
      $pairs
      # A tibble: 3 x 3
         from    to  mass
        <dbl> <dbl> <dbl>
      1     1     1 0.333
      2     2     2 0.333
      3     3     3 0.333
      
      attr(,"class")
      [1] "earthmover_dist"

# `emd()` works for data.frames

    Code
      emd(iris[1:75, ], iris[76:150, ])
    Output
      $dist
      [1] 3.853397
      
      $pairs
      # A tibble: 75 x 3
          from    to   mass
         <dbl> <dbl>  <dbl>
       1     1    36 0.0133
       2     2     5 0.0133
       3     3    22 0.0133
       4     4    68 0.0133
       5     5    30 0.0133
       6     6    46 0.0133
       7     7    40 0.0133
       8     8    53 0.0133
       9     9    32 0.0133
      10    10    59 0.0133
      # i 65 more rows
      
      attr(,"class")
      [1] "earthmover_dist"

# `emd()` works for data.frames 2

    Code
      emd(iris[1:75, ], iris[1:75, ])
    Output
      $dist
      [1] 0
      
      $pairs
      # A tibble: 75 x 3
          from    to   mass
         <dbl> <dbl>  <dbl>
       1     1     1 0.0133
       2     2     2 0.0133
       3     3     3 0.0133
       4     4     4 0.0133
       5     5     5 0.0133
       6     6     6 0.0133
       7     7     7 0.0133
       8     8     8 0.0133
       9     9     9 0.0133
      10    10    10 0.0133
      # i 65 more rows
      
      attr(,"class")
      [1] "earthmover_dist"

