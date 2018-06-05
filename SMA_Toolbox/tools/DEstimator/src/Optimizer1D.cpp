/**
 * @file Optimizer1D.cpp
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 08-21-2014
 * @brief Class definitions for Optimizer1D
 */
#include <utility>

#include "Optimizer1D.h"
#include "DCompLib.h"
using std::swap;

template<class FloatT>
const FloatT Optimizer1D<FloatT>::phi=arma::Datum<FloatT>::gratio;

template<class FloatT>
const FloatT Optimizer1D<FloatT>::phi_inv=1/arma::Datum<FloatT>::gratio;

template<class FloatT>
const FloatT Optimizer1D<FloatT>::phi_conj=1-1/arma::Datum<FloatT>::gratio;

template<class FloatT>
const FloatT Optimizer1D<FloatT>::x_tolerance=sqrt(std::numeric_limits<FloatT>::epsilon());

template<class FloatT>
const FloatT Optimizer1D<FloatT>::eps=std::numeric_limits<FloatT>::epsilon();

template<class FloatT>
const FloatT Optimizer1D<FloatT>::max_search_ratio=100;

// template<class FloatT>
// const FloatT Optimizer1D<FloatT>::max_iter=1000;
template<class FloatT>
Optimizer1D<FloatT>::Optimizer1D(const FuncT &_func, int _max_eval) 
    : func(_func), max_eval(_max_eval), N(0), 
      X(max_eval), F(max_eval)
{
}

template<class FloatT>
void Optimizer1D<FloatT>::getStats(VecT &_X, VecT &_F) const
{
    _X=X.subvec(0,N-1);
    _F=F.subvec(0,N-1);
}

template<class FloatT>
int Optimizer1D<FloatT>::getNFcalls() const
{
    return N;
}

/**
 * 
 * 
 * 
 */
template<class FloatT>
void Optimizer1D<FloatT>::maximize(FloatT xA, FloatT xB, FloatT &xmax, FloatT &Fmax)
{
    N=0;
    maximize_mode=true;
    auto bracket=bracket_min(xA,xB);
    int x=brent_min(bracket);
    xmax=X(x);
    Fmax=-F(x);
}

template<class FloatT>
void Optimizer1D<FloatT>::minimize(FloatT xA, FloatT xB, FloatT &xmin, FloatT &Fmin)
{
    N=0;
    maximize_mode=false;
    auto bracket=bracket_min(xA,xB);
    int x=brent_min(bracket);
    xmin=X(x);
    Fmin=F(x);
}


template<class FloatT>
int Optimizer1D<FloatT>::eval(FloatT x)
{
    FloatT Fval=func(x);
    if(maximize_mode) Fval=-Fval;
    X(N)=x;
    F(N)=Fval;
    return N++; //index to recorded eval
}

template<class FloatT>
inline
FloatT Optimizer1D<FloatT>::golden_step(int alpha, int beta) const
{
    return X(alpha)+phi*(X(alpha)-X(beta));
}

template<class FloatT>
FloatT Optimizer1D<FloatT>::parabolic_min(int a, int b, int c) const
{
    FloatT d1=X(b)-X(a);
    FloatT d2=X(b)-X(c);
    FloatT Q1=d1*(F(b)-F(c));
    FloatT Q2=d2*(F(b)-F(a));
    FloatT numer=(d1*Q1-d2*Q2);
    FloatT denom=-2*(Q1-Q2);
    denom= std::max(static_cast<FloatT>(fabs(denom)),eps)*sgn(denom);
    return numer/denom;
}



template<class FloatT>
typename Optimizer1D<FloatT>::IVecT 
Optimizer1D<FloatT>::bracket_min(FloatT xA, FloatT xB)
{
    if(xA==xB) throw std::logic_error("Initial points are equal");
    int a=eval(xA);
    int b=eval(xB);
    if(F(b)==F(a)) throw std::logic_error("Initial points have equal Function values");
    if(F(b)>F(a)) swap(a,b); //Now Fa=>Fb and we move downhill
    int c=eval(golden_step(b,a));
    int u;
    while(F(b)>F(c)) { //Looking for Fa > Fb < Fc
        FloatT ux=parabolic_min(b,a,c);
        FloatT ux_limit=X(b)+max_search_ratio*(X(c)-X(b));
        if((X(b)-ux)*(ux-X(c)) > 0) { // ux between b and c
            u=eval(ux);
            printf("ux:%.16g between b:%.16g and c:%.16g\n",ux,X(b),X(c));
            if(F(u)<F(c)) {//found a minimum between b and c
                printf("u is minimum between b and c\n");
                shift(a,b,u);
                break;
            } else if (F(u) > F(b)) { //b is the minium bracketed between a and u
                printf("b is minimum between a and u\n");
                c=u;
                break;
            }
            ux=golden_step(c,b); //Parabolic u no good. Try golden step instead
            u=eval(ux);
            printf("Parabolic no good.  Trying u=F(%.16g)=%.16g \n",X(u),F(u));
        } else if ((X(c)-ux)*(ux-ux_limit) > 0) {
            u=eval(ux);
            printf("ux:%.16g between c:%.16g and ulim:%.16g\n",ux,X(b),ux_limit);
        } else if ((ux-ux_limit)*(ux_limit-X(c))>=0) {
            ux=ux_limit;
            u=eval(ux);
            printf("Past ux limit, trying to use ux_limit=F(%.16g)=%.16g\n",X(u),F(u));
        } else {
            u=eval(golden_step(c,b));
            printf("Fail, trying a golden step u=F(%.16g)=%.16g\n",X(u),F(u));
        }
        printf("PreShift: a=F(%.16g)=%.16g] b=F(%.16g)=%.16g] c=F(%.16g)=%.16g]\n",X(a),F(a),X(b),F(b),X(c),F(c));
        shift(a,b,c,u);
        printf("PostShift: a=F(%.16g)=%.16g] b=F(%.16g)=%.16g] c=F(%.16g)=%.16g]\n",X(a),F(a),X(b),F(b),X(c),F(c));
    }
    for(int i=0; i<N-1; i++){
        printf("%i: F(%.16g)=%.16g\n",i,X(i),F(i));
    }
    IVecT bracket={a,b,c};
    printf("A[F(%.16g)=%.16g] B[F(%.16g)=%.16g] C[F(%.16g)=%.16g]\n",X(bracket(0)),F(bracket(0)),X(bracket(1)),F(bracket(1)),X(bracket(2)),F(bracket(2)));
    if ((X(a)>=X(b) && X(c)>=X(b)) || (X(a)<=X(b) && X(c)<=X(b))) throw std::logic_error("Bracket X values out of order");
    if (F(a)<=F(b) || F(c)<=F(b)) throw std::logic_error("Bracket F values out of order");
    return bracket;
}

template<class FloatT>
int Optimizer1D<FloatT>::golden_min(const IVecT &bracket)
{
    int a=bracket(0);
    int b=bracket(1);
    int c;
    int d=bracket(2);
    if(fabs(X(a)-X(b)) < fabs(X(b)-X(d))) {
        c=eval(X(b) + phi_conj*(X(d)-X(b)));
    } else {
        c=bracket(1);
        b=eval(X(c) + phi_conj*(X(a)-X(c)));
    }
    while( fabs(X(d)-X(a)) > x_tolerance*(fabs(X(b))+fabs(X(c)))) {
//         printf("Golden Delta(a-d):%.15f a=F(%.16g)=%.16g] b=F(%.16g)=%.16g] c=F(%.16g)=%.16g]\n",fabs(X(d)-X(a)),X(a),F(a),X(b),F(b),X(c),F(c));
        if(N==max_eval)      break;
        else if(F(c)==F(b))  break;
        else if(F(c) < F(b)) shift(a,b,c,eval(phi_inv*X(c) + phi_conj*X(d)));
        else                 shift(d,c,b,eval(phi_inv*X(b) + phi_conj*X(a)));
    }
    return F(b)<=F(c) ? b : c;
}

template<class FloatT>
int Optimizer1D<FloatT>::brent_min(const IVecT &bracket)
{
    int a=std::min(bracket(0),bracket(2)); //Lower bracket for minimum
    int b=std::max(bracket(0),bracket(2)); //upper bracket for minimum
    int x,v,w;
    x=v=w=bracket(1);
//     FloatT ux;
    int N0=N;
    VecT steps(max_eval-N0,arma::fill::zeros);
    for(int i=0; i<max_eval-N0; i++){
        FloatT delta=X(b)-X(a);
        FloatT xm=0.5*(X(a)+X(b));
        FloatT tol=x_tolerance*fabs(X(x))+eps*1e-3;
        FloatT tol2=tol*2;
        if( fabs(X(x)-xm) <= tol2-0.5*(X(b)-X(a))) return x;//Found minimum
        if(F(a)==F(b) && F(a)==F(x)) return x; //Found minimum
        printf("[I:%i]\n",i);
        FloatT pstep=parabolic_min(w,x,v);
        if(!std::isfinite(pstep) || fabs(pstep)<tol) pstep= std::signbit(pstep) ? -tol : tol;
        FloatT gstep=phi_conj*( X(x)>=xm ? X(a)-X(x) : X(b)-X(x) );
        FloatT igstep=phi*(X(x)>=xm ? X(x)-X(b) : X(x)-X(a));
        FloatT pmax_size=max_interval_size(a,b,x,pstep);
        FloatT gmax_size=max_interval_size(a,b,x,gstep);
        FloatT igmax_size=max_interval_size(a,b,x,igstep);
        printf("Parabolic step=%.16g %.16g->%.16g: max_size=%.16g\n",pstep,X(x),X(x)+pstep, pmax_size);
        printf("Golden step=%.16g %.16g->%.16g: max_size=%.16g\n",gstep,X(x),X(x)+gstep, gmax_size);
        printf("InvGolden step=%.16g %.16g->%.16g: max_size=%.16g\n",igstep,X(x),X(x)+igstep, igmax_size);
        FloatT ratio= (X(x)-X(a))/(X(b)-X(x));
        if(pstep>0) ratio=1/ratio;
        printf("ratio:%.16g\n",ratio);
        if (ratio<0.1*phi){
            printf("[%i][InvGolden Step Chosen]\n",i);
            steps(i)= igstep;
        } else if (std::isfinite(pstep) && ratio>0.1*phi_inv && pmax_size<=gmax_size){
            printf("[%i][Parabolic Step Chosen]\n",i);
            steps(i)= pstep;
        } else {
            printf("[%i][Golden Step Chosen]\n",i);
            steps(i)= gstep;
        }
        /*
            if(!std::isfinite(step) || ux<=X(a) || ux>=X(b) || fabs(step)>0.5*fabs(step3) ){
                //Rejecting ineffectual step.  Going to golden step
                step2= X(x)>=xm ? X(a)-X(x) : X(b)-X(x);
                step=phi_conj*step2;
                printf("Rejecting parabolic, new golden step: %.16g %.16g->%.16g\n",step,X(x),X(x)+step);
            }
            if( ux-X(a) < tol2 || X(b)-ux < tol2){
                //Changed the logic from NR I think they have a bug
                step=tol*sgn(step); //Don't evaluate too close to the endpoints
                printf("Adjusted step to: %.16g %.16g->%.16g\n",step,X(x),X(x)+step);
            }
        } else { //golden step
            //last step was poor.  Going to golden step.
            step2= X(x)>=xm ? X(a)-X(x) : X(b)-X(x);
            step=phi_conj*step2;
        }*/
        if(fabs(steps(i))<tol) steps(i)+=tol*sgn(steps(i));
        int u=eval(X(x)+steps(i));
        printf("I=%i Final step=%.16g to : u=F(%.16g)=%.16g\n",i,steps(i),X(u),F(u));
        if(F(u)<F(x)){
            printf("Found new minimum!  u=F(%.16g)=%.16g\n",X(u),F(u));
            if(X(u)>=X(x)) a=x;
            else           b=x;
            shift(v,w,x,u);
        } else {
            
            if(X(u)<X(x)) a=u;
            else          b=u;
            if(F(u)<=F(w) || X(w)==X(x)) shift(v,w,u);
            else if(F(u)<=F(v) || v==x || v==w) v=u;
        }
        printf("a=F(%.16g)=%.16g, b=F(%.16g)=%.16g, x=F(%.16g)=%.16g, w=F(%.16g)=%.16g, v=F(%.16g)=%.16g\n",
                X(a),F(a),X(b),F(b),X(x),F(x),X(w),F(w),X(v),F(v));
        printf("[%i] Delta=%.16g Delta'=%.16g, relative change:%.16g\n",i,delta,X(b)-X(a), delta/(X(b)-X(a)));
    }
    throw std::logic_error("Too many iterations");
    return x;
}

template<class FloatT>
inline
FloatT Optimizer1D<FloatT>::max_interval_size(int a, int b, int x, FloatT step) const
{
    using std::max;
    FloatT ux=X(x)+step;
    return max(max(fabs(X(a)-X(x)), fabs(X(b)-X(x))), max(fabs(ux-X(a)), fabs(ux-X(b))));
}
/*
template<class FloatT>
int Optimizer1D<FloatT>::brent_min(const IVecT &bracket)
{
    int a=std::min(bracket(0),bracket(2)); //Lower bracket for minimum
    int b=std::max(bracket(0),bracket(2)); //upper bracket for minimum
    int x,v,w;
    x=v=w=bracket(1);
    FloatT step=0.0,step2=0.0;
    FloatT ux;
    int N0=N;
    while(N<max_eval) {
        FloatT xm=0.5*(X(a)+X(b));
        FloatT tol=x_tolerance*fabs(X(x))+eps*1e-3;
        FloatT tol2=tol*2;
        if( fabs(X(x)-xm) <= tol2-0.5*(X(b)-X(a)) ) return x; //Found minimum
        //Save old step sizes
        FloatT step3=step2;
        step2=step;
        printf("[I:%i] Step2: %.16g Step3: %.16g\n",N-N0,step2,step3);
        if(fabs(step3)>tol && x!=v) { //parabolic step
            step=parabolic_min(w,x,v);
            ux=X(x)+step;
            printf("Parabolic step=%.16g %.16g->%.16g\n",step,X(x),ux);
            if(!std::isfinite(step) || ux<=X(a) || ux>=X(b) || fabs(step)>0.5*fabs(step3) ){
                //Rejecting ineffectual step.  Going to golden step
                step2= X(x)>=xm ? X(a)-X(x) : X(b)-X(x);
                step=phi_conj*step2;
                printf("Rejecting parabolic, new golden step: %.16g %.16g->%.16g\n",step,X(x),X(x)+step);
            }
            if( ux-X(a) < tol2 || X(b)-ux < tol2){
                //Changed the logic from NR I think they have a bug
                step=tol*sgn(step); //Don't evaluate too close to the endpoints
                printf("Adjusted step to: %.16g %.16g->%.16g\n",step,X(x),X(x)+step);
            }
        } else { //golden step
            //last step was poor.  Going to golden step.
            step2= X(x)>=xm ? X(a)-X(x) : X(b)-X(x);
            step=phi_conj*step2;
            printf("Golden step=%.16g %.16g->%.16g\n",step,X(x),X(x)+step);
        }
        if(fabs(step)<tol) step+=tol*sgn(step);
        int u=eval(X(x)+step);
        printf("I=%i Final step=%.16g to : u=F(%.16g)=%.16g\n",N-N0,step,X(u),F(u));
        if(F(u)<F(x)){
            printf("Found new minimum!  u=F(%.16g)=%.16g\n",X(u),F(u));
            if(X(u)>=X(x)) a=x;
            else           b=x;
            shift(v,w,x,u);
        } else {
            
            if(X(u)<X(x)) a=u;
            else          b=u;
            if(F(u)<=F(w) || X(w)==X(x)) shift(v,w,u);
            else if(F(u)<=F(v) || v==x || v==w) v=u;
        }
        printf("a=F(%.16g)=%.16g, b=F(%.16g)=%.16g, x=F(%.16g)=%.16g, w=F(%.16g)=%.16g, v=F(%.16g)=%.16g\n",
                X(a),F(a),X(b),F(b),X(x),F(x),X(w),F(w),X(v),F(v));
    }
    throw std::logic_error("Too many iterations");
    return x;
}*/

// 
// 
//     
//     min_x=bracket[1]; //Minimum found so far
//     FloatT min_x2=min_x; //Second best minimum found so far
//     FloatT min_x3=min_x; //Previous x2
//     min_F=Fbracket[1];
//     FloatT min_F2=min_F, min_F3=min_F;
//     FloatT step=0; 
//     FloatT step2=0; //Previous step
//     FloatT cand_x;
//     for(int i=0; i<max_iter; i++){
//         FloatT xmid=0.5*(A+B);
//         FloatT tol=x_tolerance*fabs(min_x)+eps*1e-3;
//         FloatT tol2=tol*2;
//         printf("I:%i [A:%.9g B:%.9g] [minx:%.9g minx2:%.9g minx3:%.9g] [minF:%.9g minF2:%.9g minF3:%.9g]\n",i,A,B,min_x,min_x2,min_x3,min_F, min_F2, min_F3);
//         printf("Xmid:%.9g tol:%.9g tol2:%.9g\n",xmid,tol,tol2);
//         if( fabs(min_x-xmid) <= tol2-0.5*(B-A) ) return //converged
//         FloatT step=parabolic_max(min_x2, min_x, min_x3);
//         if(cand_x-A < tol2 || B-cand_x < tol2) step=tol1*sgn(xmid-min_x); //Don't evaluate too close to the endpoints
//         FloatT cand_x=min_x+step;
//         if(fabs(step2) <= tol1 || //Previous step was small
//            fabs(step) >= 0.5*fabs(step2) //step not significantly smaller
//            cand_x<=A || cand_x>=B ) { //OOB
//             //Parabolic failed Golden step
//             printf("Parabolic step failed, falling back to golden section\n");
//             if(fabs(step2) <= tol1) printf("^^^^&&&& Weird! \n");
//             step2= min_x >= xmid ? A-min_x : B-min_x;
//             step=phi*step2;
//             cand_x=min_x+step;
//         } else {
//             printf("Parabolic step SUCCESS\n");
//             step2=step; //Accept step and record it
//         }
//         cand_F = F(cand_x);
//         printf("step:%.9g cand_x:%.9x candF:%.9g\n",step,cand_x,candF);
//         if (cand_F <= min_F) { //Found new min
//             if (cand_x >= min_x) A=min_x; else B=min_x;
//             shift(min_x3, min_x2, min_x, cand_x);
//             shift(min_F3, min_F2, min_F, cand_F);
//         } else { // min_x is still the best.  Narrow brackets
//             if (cand_x < min_x) A=cand_x; else B=cand_x;
//             if (cand_F <= min_F2 || min_x2 == min_x){
//                 swap(min_x2, min_x3);
//                 swap(min_F2, min_F3);
//             } else if (cand_F <= min_F3 || min_x3 == min_x || min_x3==min_x2 ) {
//                 min_x3=cand_x;
//                 min_F3=cand_F;
//             }
//         }
//     }
//     
// }


/* Explicit Template Intatiations */
template class Optimizer1D<float>;
template class Optimizer1D<double>;
