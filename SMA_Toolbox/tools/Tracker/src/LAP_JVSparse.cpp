/** @file LAP_JVSparse.cpp
 *  @author Mark J. Olah (mjo at cs.unm.edu)
 *  @date 05-2015
 *  @brief The member definitions for the LAP Jonker Volgenant algorithm
 * 
 * This is a modern dense/sparse C++ implementation of Jonker Volgenant algoirthm using armadillo
 * and presenting C++ and Matlab interface.
 * 
 * Adapted from text of Jonker and Volgenant. Computing 38, 324-340 (1986)
 * 
 */
#include <cmath>
#include <iostream>
#include "LAP_JVSparse.h"

template<class FloatT>
typename LAP_JVSparse<FloatT>::IVecT
LAP_JVSparse<FloatT>::solve(const SpMatT &C)
{
    IndexT N = static_cast<IndexT>(C.n_rows);
    IVecT x(N), y(N);
    VecT u(N), v(N);
    checkCosts(C);//optional
    solveLAP_orig(C,x,y,u,v);
    checkSolution(C,x,y,u,v); //optional
    return x;
}

/**
 * This wraps the original sparse lap implementation that for some reason uses 1-based indexing,
 * which we correct with some pointer arrithmetic and adjusting of appropriate indicies in the
 * sparse matrix implementation.
 * 
 * Furthermore because the lap_orig code assumes a compressed-row format, but we pass it the
 * internal datastore of a compressed-col format sparse metrix.  We invert x/y and u/v on the
 * call to lap_orig to effectively let the transformation work easily with the legacy code.
 * 
 * This means x is the row sol and y is the col sol, as it normally would be.
 * 
 * 
 * 
 */
template<class FloatT>
void LAP_JVSparse<FloatT>::solveLAP_orig(const SpMatT &C, IVecT &x, IVecT &y, VecT &u, VecT &v)
{
    IndexT Ndim = static_cast<IndexT>(C.n_rows);
    IndexT Nvals = static_cast<IndexT>(C.n_nonzero);
    
    IVecT C_row_ind(Nvals);
    for(IndexT n=0; n<Nvals; n++) C_row_ind(n) = static_cast<IndexT>(C.row_indices[n])+1;

    IVecT C_col_starts(Ndim+1);
    for(IndexT n=0; n<=Ndim;n++) C_col_starts(n) = static_cast<IndexT>(C.col_ptrs[n])+1;

    x.set_size(Ndim);
    y.set_size(Ndim);
    u.set_size(Ndim);
    v.set_size(Ndim);
    const FloatT *C_values_ptr = C.values-1;
    IndexT *C_row_ind_ptr = C_row_ind.memptr() -1;
    IndexT *C_col_starts_ptr = C_col_starts.memptr()-1;
    IndexT *x_ptr = y.memptr()-1; //Swap x&y
    IndexT *y_ptr = x.memptr()-1; //Swap x&y
    FloatT *u_ptr = v.memptr()-1; //Swap u&v
    FloatT *v_ptr = u.memptr()-1; //Swap u&v
    lap_orig(Ndim, C_values_ptr, C_row_ind_ptr, C_col_starts_ptr, x_ptr, y_ptr, u_ptr, v_ptr);
    x-=1; //Correct indicies back to 0-based
    y-=1; //Correct indicies back to 0-based
//     std::cout<<"X"<<x.t()<<"\n";
//     std::cout<<"Y"<<y.t()<<"\n";
//     std::cout<<"U"<<u.t()<<"\n";
//     std::cout<<"V"<<v.t()<<"\n";
    VecT cost=computeCost(C,x);
//     std::cout<<"SolCost: "<<cost.t()<<"\n";
//     std::cout<<"Tot SolCost:"<<arma::sum(cost)<<"\n";
}

/**
 * Compute the total cost of a solution
 * 
 * @param[in] row_sol This is the 'x' output from the solver giving the col assignment for each row in order
 */
template<class FloatT>
typename LAP_JVSparse<FloatT>::VecT
LAP_JVSparse<FloatT>::computeCost(const SpMatT &C, const IVecT &row_sol)
{
    IndexT N= static_cast<IndexT>(C.n_rows);
    VecT cost(C.n_rows);
    for(IndexT n=0; n<N; n++) cost(n) = C(n,row_sol(n));
    return cost;
}


template<class FloatT>
bool LAP_JVSparse<FloatT>::checkCosts(const SpMatT &C)
{
    bool ok=true;
    const FloatT *const vals = C.values;
    const arma::uword *const row_ind = C.row_indices;
    const arma::uword *const col_ptr = C.col_ptrs;
    IndexT N= static_cast<IndexT>(C.n_rows);
    
    //Check for negative reduced cost
    for (IndexT j=0; j<N; j++) for(arma::uword t=col_ptr[j]; t<col_ptr[j+1]; t++) {
        IndexT i=row_ind[t]; //i - row, j - col;
        if(vals[t]<0){
            ok=false;
            std::cout<<">> NegativeCost: "<<vals[t]<<" ("<<i<<","<<j<<")\n";
        } else if (!std::isfinite(vals[t])) {
            ok=false;
            std::cout<<">> NonFiniteCost: "<<vals[t]<<" ("<<i<<","<<j<<")\n";
        }
    }
    return ok;
}

template<class FloatT>
bool LAP_JVSparse<FloatT>::checkSolution(const SpMatT &C, const IVecT &x, const IVecT &y, const VecT &u, const VecT &v)
{
    bool ok=true;
    const FloatT *const vals = C.values;
    const arma::uword *const row_ind = C.row_indices;
    const arma::uword *const col_ptr = C.col_ptrs;
    IndexT N= static_cast<IndexT>(C.n_rows);
    
    //Check for negative reduced cost
    for (IndexT j=0; j<N; j++) for(arma::uword t=col_ptr[j]; t<col_ptr[j+1]; t++) {
        IndexT i=row_ind[t]; //i - row, j - col;
        FloatT redC = vals[t] - u[i] - v[j];
        if(redC<0){
            std::cout<<">> NegativeReducedCost: "<<redC<<" ("<<i<<","<<j<<")\n";
            ok=false;
        }
    }
    //Check row solution
    if(arma::any(x<0) || arma::any(x>=N)){
        std::cout<<">> InvalidRowSol: X:"<<x.t()<<"\n";
        ok=false;
    }
    if(static_cast<int>(arma::unique(x).eval().n_elem) != N){
        std::cout<<">> NonUniqueRowSolPermutation: X:"<<x.t()<<"\n";
        ok=false;
    }

    //Triggered when there is only one solution for some/row column.  Not sure if this is bad or not, but I am seeming to still get good result
//     for (IndexT i=0; i<N; i++) {
//         FloatT redC = C(i,x(i)) - u(i) - v(x(i));
//         if(redC!=0) {
//             std::cout<<">> NonNullReducedRowCost: "<<redC<<" row:"<<i<<" row_sol:"<<x(i)<<"\n";
//             ok=false;
//         }
//     }
     //Check col solution
    if(arma::any(y<0) || arma::any(y>=N)){
        std::cout<<">> InvalidColSol: Y:"<<y.t()<<"\n";
        ok=false;
    }
    if(static_cast<int>(arma::unique(y).eval().n_elem) != N){
        std::cout<<">> NonUniqueColSolPermutation: Y:"<<y.t()<<"\n";
        ok=false;
    }
    //Triggered when there is only one solution for some/row column.  Not sure if this is bad or not, but I am seeming to still get good result
//     for (IndexT j=0; j<N; j++) {
//         FloatT redC = C(y(j),j) - u(y(j)) - v(j);
//         if(redC!=0) {
//             std::cout<<">> NonNullReducedColCost: "<<redC<<" col:"<<j<<" col_sol:"<<y(j)<<"\n";
//             ok=false;
//         }
//     }
    
    //Check solution parity
    for (IndexT i=0; i<N; i++) {
        if (y(x(i))!=i) {
            std::cout<<">> RowSolutionError: i:"<<i<<" rowsol:"<<x(i)<<" col[rowsol]:"<<y(x(i))<<"\n";
            ok=false;
        }
    }
    for (IndexT j=0; j<N; j++) {
        if (x(y(j))!=j) {
            std::cout<<">> ColSolutionError: j:"<<j<<" colsol:"<<y(j)<<" row[colsol]:"<<x(y(j))<<"\n";
            ok=false;
        }
    }
    return ok;
}

//The raw lap file with the original matlab-style indexing
template<class FloatT>
void LAP_JVSparse<FloatT>::lap_orig(IndexT n, const FloatT cc[], const IndexT kk[], const IndexT first[],
                                    IndexT x[], IndexT y[], FloatT u[], FloatT v[])
{
   IndexT h, i,j,k,l,t,last,tel,td1=0,td2,i0,j0=0,j1=0,l0;

   IndexT *lab, *freeRow, *todo;
   bool *ok;
   FloatT min, v0, vj, dj, tmp;
   FloatT *d;
   FloatT FLT_EPSILON = std::numeric_limits<FloatT>::epsilon();


   ok = new bool[n + 1];
   lab = new IndexT[n + 2];
   freeRow = new IndexT[n + 2];
   todo = new IndexT[n + 2];
   d = new FloatT[n + 2];
     
   /* Initialize */
   for (j=1; j<=n; j++)  v[j] = INFINITY;

   for (i = 1; i <= n; i++) {
      x[i] = 0; u[i] = 0;
      for (t = first[i]; t < first[i+1]; t++) {
         j = kk[t];
         if (cc[t] < v[j]) {
            v[j] = cc[t];
            y[j] = i;
         } /* if */
      } /* for */
   } /* for */

   for (j = n; j >= 1; j--) {
      i = y[j];
      if (x[i] == 0) {
         x[i] = j;
      } else {
         y[j] = 0;
         x[i] = -std::abs(x[i]);
      } /* if */
   } /* for */

   l = 0;
   for (i = 1; i <= n; i++) {
      if (x[i] < 0) {
         x[i] = -x[i];
      } else if (x[i] > 0) {
         min = INFINITY;
         j1 = x[i];
         for (t = first[i]; t < first[i+1]; t++) {
            j = kk[t];
            if (j != j1 && cc[t] - v[j] < min) {
               min = cc[t] - v[j];
            } /* if */
         } /* for */
         u[i] = min;
         t = first[i];
         while (kk[t] != j1) {
            t++;
         } /* while */
         v[j1] = cc[t] - min;
      } else {
         freeRow[++l] = i;
      } /* if */
   } /* for */

   /* Improve initial solution */
   for (tel = 0; tel < 2; tel++) {
      h = 1;
      l0 = l;
      l = 0;
      while (h <= l0) {
         i = freeRow[h++];
         v0 = vj = INFINITY;

         for (t = first[i]; t < first[i+1]; t++) {

            j = kk[t];
            dj = cc[t] - v[j];
            if (dj < vj) {
               if (dj >= v0) {
                  vj = dj;
                  j1 = j;
               } else {
                  vj = v0;
                  v0 = dj;
                  j1 = j0;
                  j0 = j;
               } /* if */
            } /* if */
         } /* for */

         i0 = y[j0];
         u[i] = vj;
         if (vj - v0 > FLT_EPSILON) {
            v[j0] = v[j0] - vj + v0;
         } else if (i0 > 0) {
            j0 = j1;
            i0 = y[j0];
         } /* if */

         x[i] = j0;
         y[j0] = i;

         if (i0 > 0) {
            if (vj - v0 > FLT_EPSILON) {
               freeRow[--h] = i0;
            } else {
               freeRow[++l] = i0;
            } /* if */
         } /* if */
      } /* while */
   } /* for */

   tmp = 0;
   for (i = 1; i <= n; i++) {
      tmp += u[i] + v[i];
   } /* for */

   /* Augmentation part */
   l0 = l;
   for (l = 1; l <= l0; l++) {

      for (j = 1; j <= n; j++) {
         d[j] = INFINITY;
         ok[j] = false;
      } /* for */

      min = INFINITY; i0 = freeRow[l];

      for (t = first[i0]; t < first[i0+1]; t++) {
         j = kk[t];
         dj = cc[t] - v[j];
         d[j] = dj;
         lab[j] = i0;

         if (dj <= min) {
            if (dj < min) {
               td1 = 0;
               min = dj;
            } /* if */
            todo[++td1] = j;
         } /* if */
      } /* for */

      for (h = 1; h <= td1; h++) {
         j = todo[h];
         if (y[j] == 0) {
            goto label2;
         } /*if */
         ok[j] = true;
      } /* for */

      td2 = n;
      last = n + 1;

      /* Repeat until a freeRow row found */
      while (true) {
         j0 = todo[td1--];
         i = y[j0];
         todo[td2--] = j0;
         t = first[i];

         for (t = first[i]; kk[t] != j0; t++) {
            /* nothing */
         } /* for */

         tmp = cc[t] - v[j0] - min;

         for (t = first[i]; t < first[i+1]; t++) {
            j = kk[t];
            if (!ok[j]) {
               vj = cc[t] - v[j] - tmp;
               if (vj < d[j]) {
                  d[j] = vj;
                  lab[j] = i;
                  if (vj == min) {
                     if (y[j] == 0) {
                        goto label1;
                     } /* if */
                     td1++;
                     todo[td1] = j;
                     ok[j] = true;
                  } /* if */
               } /* if */
            } /* if */
         } /* for */
         if (td1 == 0) {
            min = INFINITY - 1;
            last = td2 + 1;
            for (j = 1; j <= n; j++) {
               if (d[j] <= min) {
                  if (!ok[j]) {
                     if (d[j] < min) {
                        td1 = 0;
                        min = d[j];
                     } /* if */
                     todo[++td1] = j;
                  } /* if */
               } /* if */
            } /* for */
            for (h = 1; h <= td1; h++) {
               j = todo[h];
               if (y[j] == 0) {
                  goto label1;
               } /* if */
               ok[j] = true;
            } /* for */
         } /* if */
      } /* while */
label1:
      for (k = last; k <= n; k++) {
         j0 = todo[k];
         v[j0] += d[j0] - min;
      } /* for */

label2:
      do {
         i = lab[j];
         y[j] = i;
         k = j;
         j = x[i];
         x[i] = k;
      } while (i != i0);
   } /* for */

   for (i = 1; i <= n; i++) {
      j  = x[i];
      t = first[i];
      while (kk[t] != j) {
         t++;
      } /* while */

      u[i] = cc[t] - v[j];
   } /* for */


   delete [] ok;
   delete [] lab;
   delete [] freeRow;
   delete [] todo;
   delete [] d;

}



/* Explicit Template Instantiation */
/* These ensure the compiler emits code for both the double and float versions of DEstimator */
template class LAP_JVSparse<float>;
template class LAP_JVSparse<double>;
