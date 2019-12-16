##  HW5: Lasso implementation in C++

# Introduction
In this homework, you will practice C++ by implementing the LASSO coordinate-descent algorithm from HW4. We will only focus on fitLASSOstandardized functions and corresponding helpers as you will see the biggest speed improvement there. For simplicity, we will also avoid doing any compatibility checks in C++. In practice, you will often have an R code wrapper with all compatibility checks that then calls the corresponding function in C++.

# Functions Instructions

**LassoInC.cpp** contains the starter for C++ code. You will have to modify this starter to create the following functions:

  - **soft_c** - this should be the C++ version of the soft function from HW4. It should return the soft-thresholded value of a. 
  - **lasso_c** - this should be the C++ version of the lasso function from HW4. It should return the value of lasso objective function at given arguments.
  - **fitLASSOstandardized_c** - this should be the C++ version of the fitLASSOstandardixed function from HW4 with the following modification: it only returns $\beta$ vector.
  - **fitLASSOstandardized_seq_c** - this should be the C++ version of the fitLASSOstandardixed function from HW4 with the following modifications: it only takes used-supplied lambda sequence and it only returns matrix of $\beta$ values.
  
All these functions should be directly available from R after executing **sourceCpp** command with corresponding source file. 
  
**Keep in mind:** for this assignment, we will use Armadillo library classes for matrices and vectors. You should check the examples in class and references for Armadillo library to find how to translate your R code in C++ appropriately with that library.  
  
**TestsCpp.R**  contains starter file for testing your C++ functions against your R functions. You should do the following

  - upload your LASSOfunctions.R script from HW4 to this project repository
  - develop 2 tests for each function **to check equality of the corresponding outputs between your R and C versions**
  - do 2 microbenchmarks comparisons as per comments
  - do speed test on riboflavin data as coded in the end of that file
  
# Grading criteria

Your assingment will be graded based on 

 * correctness (50%)

Take advantage of the fact that you have a correct R function to extensively test your C code. You also have to develop tests yourself for this particular assignment.
 
 * speed (30%) 
 
You should expect to see an amazing speed improvement when moving your code to C++. On riboflavin data, my C++ fitLASSOstandardized_seq code with 30 lambdas runs in 77 milliseconds (compared to several seconds for corresponding R version). As always, you will get +5 points if your correct C++ code is considerably faster than mine.

 * code style/documentation (10%)

You need to comment different parts of the code so it's clear what they do, have good identation, readable code with names that make sense.
 
 * version control/commit practices (10%)
 
 I expect you to start early on this assignment, and work gradually. You want to commit often, have logically organized commits with short description that makes sense.  
  
  

