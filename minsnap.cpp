#include <iostream>
#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <list>
#include <math.h>
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif
// program and solution types
typedef CGAL::Quadratic_program_from_iterators
<int**,                                                // for A
 int*,                                                 // for b
 CGAL::Const_oneset_iterator<CGAL::Comparison_result>, // for r
 bool*,                                                // for fl
 int*,                                                 // for l
 bool*,                                                // for fu
 int*,                                                 // for u
 int**,                                                // for D
 int*>                                                 // for c 
Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

const int m = 12;
const int t = 10;

void print(float matrix[t-1][m+1][m+1])
{
    int i, j, z;
    for (z = 0; z < t-1; z++)
    {
      for (i = 0; i < m+1; ++i)
      {
          for (j = 0; j < m+1; ++j)
              printf("%f ", matrix[z][i][j]);
          printf("\n");
      }
  }
}

// int * concat(int list1[], int list2[]){
//   static int list[sizeof(list1)/sizeof(*list1) + sizeof(list2)/sizeof(*list2)];
//   for (int counter = 0; counter < (sizeof(list)/sizeof(*list)); counter++)
//   {
//     if (counter < sizeof(list1)/sizeof(*list1))
//     {
//       list[counter] = list1[counter];
//     }
//     else
//     {
//       list[counter] = list2[counter-(sizeof(list1)/sizeof(*list1))];
//     }
//   }
//   return list;
// }

int main() {
  
  int C[m+1];
  float pH[t-1][m+1][m+1];
  float nH[t-1][m+1][m+1];
  float H[t-1][m+1][m+1];
  float T[t];
  for (int i=0; i < t; i++)
  {
    T[i] = (float)i/(float)t;
  }

  for (int q=0; q < t-1; q++)
  {
    for (int i=0; i <= m ; i++)
    {
      C[i] = i*(i-1)*(i-2)*(i-3);
    }
    for (int i=0; i <= m ; i++)
    {
      for (int z=0; z <= m; z++)
      {
        if (z > 3 && i > 3)
        {
          pH[q][i][z] = (1./(z-3+i-4))*(float)C[i]*(float)C[z]*pow(T[q+1],(z-3+i-4));
          nH[q][i][z] = (1./(z-3+i-4))*(float)C[i]*(float)C[z]*pow(T[q],(z-3+i-4));
        }
        else
        {
          pH[q][i][z] = 0;
          nH[q][i][z] = 0;
        }
        H[q][i][z] = pH[q][i][z] - nH[q][i][z];             
      }
    }
  }
  // print(pH);

  // int *p = concat(list1,list2);
  // int list[sizeof(list1)/sizeof(*list1) + sizeof(list2)/sizeof(*list2)];
  // for (int counter = 0; counter < (sizeof(list)/sizeof(*list)); counter++)
  // {
  //   list[counter] = *(p+counter);
  //   //std::cout << *(p+counter) <<std::endl;
  // }
  int  Ax[] = {1, -1};                        // column for x
  int  Ay[] = {1,  2};                        // column for y
  int*  A[] = {Ax, Ay};                       // A comes columnwise
  int   b[] = {7, 4};                         // right-hand side
  CGAL::Const_oneset_iterator<CGAL::Comparison_result> 
        r(    CGAL::SMALLER);                 // constraints are "<="
  bool fl[] = {true, true};                   // both x, y are lower-bounded
  int   l[] = {0, 0};
  bool fu[] = {false, true};                  // only y is upper-bounded
  int   u[] = {0, 4};                         // x's u-entry is ignored
  int  D1[] = {2};                            // 2D_{1,1}
  int  D2[] = {0, 8};                         // 2D_{2,1}, 2D_{2,2}
  int*  D[] = {D1, D2};                       // D-entries on/below diagonal
  int   c[] = {0, -32};
  int  c0   = 64;                             // constant term
  // now construct the quadratic program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program qp (2, 2, A, b, r, fl, l, fu, u, D, c, c0);
  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_quadratic_program(qp, ET());
  // output solution
  std::cout << s;
  return 0;
}