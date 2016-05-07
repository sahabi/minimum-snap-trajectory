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

const int m = 6;
const int t = 10;

void print(float matrix[(m+1)*(t-1)][(m+1)*(t-1)])
{
    int i, j;
      for (i = 0; i < (m+1)*(t-1); ++i)
      {
          for (j = 0; j < (m+1)*(t-1); ++j)
              printf("%.0f ", matrix[i][j]);
          printf("\n");
      
  }
}


float** create2DArray(unsigned height, unsigned width, float array[][m+1][m+1], int length)
{
  float** array2D = 0;
  array2D = new float*[height];
  for (int i = 0; i < length; i++)
  {
   for (int h = 0; h < m+1; h++)
    {
          array2D[i*(m+1)+h] = new float[width];

          for (int w = 0; w < m+1; w++)
          {
                array2D[i*(m+1)+h][i*(m+1)+w] = array[i][h][w];
                //printf("%1f ",array[i][h][w]);
          }
    }
  }

  return array2D;
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
  float Q[(m+1)*(t-1)][(m+1)*(t-1)];
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
  float** my2DArray = create2DArray((m+1)*(t-1),(m+1)*(t-1),H,t-1);
  //print(my2DArray);
  for (int h = 0; h < (m+1)*(t-1); h++)
      {
            for (int w = 0; w < (m+1)*(t-1); w++)
            {
                  Q[h][w] = 2*my2DArray[h][w];
            }
      }
  print(Q);

  // int *p = concat(list1,list2);
  // int list[sizeof(list1)/sizeof(*list1) + sizeof(list2)/sizeof(*list2)];
  // for (int counter = 0; counter < (sizeof(list)/sizeof(*list)); counter++)
  // {
  //   list[counter] = *(p+counter);
  //   //std::cout << *(p+counter) <<std::endl;
  // }
  int*  A[] = {};                       // A comes columnwise
  int   b[] = {};                         // right-hand side
  CGAL::Const_oneset_iterator<CGAL::Comparison_result> 
        r(    CGAL::SMALLER);                 // constraints are "<="
  bool fl[(m+1)*(t-1)];                   // both x, y are lower-bounded
  for (int h = 0; h < (m+1)*(t-1); h++)
  {
    fl[h] = true;
  } 
  int   l[(m+1)*(t-1)];
  for (int h = 0; h < (m+1)*(t-1); h++)
  {
    l[h] = 0;
  } 
  bool fu[(m+1)*(t-1)];                  // only y is upper-bounded
  for (int h = 0; h < (m+1)*(t-1); h++)
  {
    fu[h] = true;
  } 
  int   u[(m+1)*(t-1)];                         // x's u-entry is ignored
  for (int h = 0; h < (m+1)*(t-1); h++)
  {
    u[h] = 1000000;
  }   
  int*  D[(m+1)*(t-1)];                       // D-entries on/below diagonal
  for (int h = 0; h < (m+1)*(t-1); h++)
  {
    D[h] = (int*)Q[h];
  }  
  int   c[] = {};
  int  c0   = 0;                             // constant term
  // now construct the quadratic program; the first two parameters are
  // the number of variables and the number of constraints (rows of A)
  Program qp ((m+1)*(t-1), 0, A, b, r, fl, l, fu, u, D, c, c0);
  // solve the program, using ET as the exact type
  Solution s = CGAL::solve_quadratic_program(qp, ET());
  // output solution
  std::cout << s;
  return 0;
}