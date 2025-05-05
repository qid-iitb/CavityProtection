//EDITED FUNCTIONS FOR OPTIMIZATION

#include "FF.h"

#include <iostream>
#include <armadillo>
#include <math.h>
#include <algorithm>
#include <numeric>



//------------------------------------------------------------------------------------------------------------

inline arma::mat eyed(int n)
{
    arma::mat Identity(n,n,arma::fill::eye);
    return Identity;
}

inline arma::sp_cx_mat speyed(int n)
{
    arma::sp_cx_mat Sident = arma::speye<arma::sp_cx_mat>(n,n);
    return Sident;
}


inline arma::mat nullM(int n)
{
    arma::mat NullM(n,n,arma::fill::zeros);
    return NullM;
}

inline arma::mat Aop(int n)
    {
       arma::mat aop(n,n,arma::fill::zeros);
       for(int i = 0; i < n-1; i++){;
       aop(i,i+1)= sqrt(i+1);
       }
     return aop;
    }

inline arma::mat Adag(int n)
    {
       arma::mat adag(n,n,arma::fill::zeros);
       for(int i = 0; i < n; i++){;
       adag(i,i)= i;
       }
     return adag;
    }


    
    
arma::cx_mat UTT(arma::cx_mat A, double tstep){
       arma::cx_vec eigval;
       arma::cx_mat eigvec;
       eig_gen(eigval,eigvec,A);
       arma::cx_mat ut = eigvec*arma::diagmat(exp(eigval*tstep))*inv(eigvec);
       return(ut);
       }

arma::cx_cube UT(const arma::cx_cube &Le, double tstep){
        int Nspin = Le.n_slices;
        arma::cx_cube Uti;
        Uti.copy_size(Le);
        for(int i=0; i < Nspin; i++){;
          arma::cx_vec eigval;
          arma::cx_mat eigvec;
          eig_gen(eigval,eigvec,Le.slice(i));
          Uti.slice(i) = eigvec*arma::diagmat(exp(eigval*tstep))*inv(eigvec);
        }
        return Uti;
}

/*arma::vec SwapParts(arma::vec A, int k,int l){
       arma::vec B = A.swap_rows(k,l);
       return B;
}*/

arma::vec SwapParts(arma::vec A, int k,int l){
       arma::vec SP;
       SP = A;
       SP.swap_rows(k,l);
       return(SP);
       }



double prodd(const arma::vec A, int k, int l){
        double pr = 1.0;
        for(int i = k-1; i < l; i++) {
           pr = pr * A(i);
        }
        return pr;
}
   

arma::vec Breakbase(int num, const arma::vec &A){
       int lnt = A.n_rows;
       arma::vec Fi(lnt);
       for(int i = 1; i < lnt+1; i++){;
          div_t divresult;
          divresult = div(num,prodd(A,i+1,lnt));
          int uni = A(i-1);
          int pp = std::min(divresult.quot, uni-1);
          Fi(i-1) = pp;num=num-pp*prodd(A,i+1,lnt);
          }  
       return Fi;
}


int Base2Digit2(const arma::vec &B, const arma::vec &dim){
          int lnt = dim.n_rows;
          double num = 0.0;
          for(int i = 1; i < lnt+1; i++){;
            num = num + B(i-1)*prodd(dim,i+1,lnt);
            }
          return((int)num);
}       


arma::vec Breakbase2(int num, double S1, double S2){
       return {(double)(div(num,S2).quot),(double)(num - div(num,S2).quot*S2)};
}



/*int Base2Digit(arma::vec B, const arma::vec &dim){
        int lnt = dim.n_rows;
        double num = B(0)*dim(1)+B(1);
        for(int i = 2; i < lnt; i++){
            num = num*dim(i) + B(i);
            }
        return (int)num;
} */      



int Base2Digit(const arma::vec &B, const arma::vec &dim){
          int lnt = dim.n_rows;
          double num = 0.0;
          for(int i = 1; i < lnt+1; i++){;
            num += B(i-1)*prodd(dim,i+1,lnt);
            }
          return((int)num);
}       


arma::mat Col2rho(arma::vec A){
           arma::mat B;
           int lnt = A.n_rows;
           B = reshape(A,std::sqrt(lnt),std::sqrt(lnt));
           return (B.t());
}


arma::cx_mat cCol2rho(const arma::cx_vec &A){
             int lnt = std::sqrt(A.n_rows);
             return reshape(A,lnt,lnt).st();
}

arma::vec Rho2col(arma::mat A){
           arma::vec B;
           int lnt = A.n_rows;
           B = reshape(A.t(),pow(lnt,2),1);
           return (B);
}

arma::cx_vec cRho2col(const arma::cx_mat &A){
             int lnt = A.n_rows;
             return reshape(A.st(),pow(lnt,2),1);
}


arma::vec Mapleft(arma::vec A, arma::mat U, arma::vec dim){
          arma::vec mapl;
          mapl = arma::kron(U,eyed(dim(1)))*A;
          return(mapl);
}

arma::cx_vec cMapleft(arma::cx_vec A, arma::cx_mat U, arma::vec dim){
             return arma::kron(U,eyed(dim(1)))*A;
}


arma::cx_vec ctMapleft(const arma::cx_vec &A, const arma::cx_mat &U, double dim0, double dim1){
             arma::cx_mat Z = U*reshape(A,dim1,dim0).st();
             return reshape(Z.st(),Z.n_rows*Z.n_cols,1);
             //return(maple);
             //arma::cx_vec maple = reshape(Z.st(),Z.n_rows*Z.n_cols,1);
}


arma::vec Mapright(arma::vec A, arma::mat U, arma::vec dim){
          arma::vec mapr;
          mapr = arma::kron(eyed(dim(0)),U)*A;
          return(mapr);
}


arma::cx_vec cMapright(arma::cx_vec A, arma::cx_mat U, arma::vec dim){
             arma::cx_vec mapr;
             mapr = arma::kron(eyed(dim(0)),U)*A;
             return(mapr);
}

arma::cx_vec ctMapright(const arma::cx_vec &A, const arma::cx_mat &U, double dim0, double dim1){
             arma::cx_mat Z = reshape(A,dim1,dim0).st() * U.st();
             return reshape(Z.st(),Z.n_rows*Z.n_cols,1);
             //return(mapri);
             //arma::cx_vec mapri = reshape(Z.st(),Z.n_rows*Z.n_cols,1);
}


arma::vec Traceop(double si){
          double sq = std::sqrt(si);
          arma::vec B(si,arma::fill::zeros);
          for(int i=0;i<si;i+= sq+1){;
              B(i)= 1.0;
          }
          return(B);
}


double Ttr(arma::vec A, arma::vec dim){
          double len = dim.n_rows;
          arma::vec trop;
          trop = 1;
          for(int i=0;i<len;i++){
            trop = arma::kron(trop,Traceop(dim(i)));
          }
          double tr = std::inner_product(std::begin(A), std::end(A), std::begin(trop), 0.0);
          return(tr);
}


double cTtr(const arma::cx_vec &A, const arma::vec &dim){
          double len = dim.n_rows;
          arma::vec trop;
          trop = 1;
          for(int i=0;i<len;i++){
            trop = arma::kron(trop,Traceop(dim(i)));
          }
          return abs(std::inner_product(std::begin(A), std::end(A), std::begin(trop), std::complex<double>()));
}




int BLay(int Num, arma::vec dim){
     int len = dim.n_rows;
     arma::vec Dim(2);Dim.fill(std::sqrt(prodd(dim,1,len)));
     arma::vec q = Breakbase(Num,Dim);
     arma::mat ek(len,2);
     ek.col(0) = Breakbase(q(0),arma::sqrt(dim));
     ek.col(1) = Breakbase(q(1),arma::sqrt(dim));
     arma::vec New(len);arma::mat p = ek.t();
     arma::mat dsq(len,2,arma::fill::zeros);
     dsq.each_col() = arma::sqrt(dim);
     arma::mat Dsq = dsq.t();
     for(int i=0;i<len;i++){
        New(i)=Base2Digit(p.col(i), Dsq.col(i));
        }
     int bl = (int)Base2Digit(New,dim);
     return(bl);
}





arma::uvec DM2Reg(arma::vec dim){
           int len = dim.n_rows;
           int ran = (int)prodd(dim,1,len);
           arma::uvec ff(ran);
           for(int i=0;i<ran;i++){
             ff(i) = BLay(i,dim);
            }
           return ff;
}


int FLay(int Num, arma::vec dim){
     arma::vec q = Breakbase(Num,dim);
     int len = dim.n_rows;
     int leng = q.n_rows;
     arma::vec Dim(2);Dim.fill(std::sqrt(prodd(dim,1,len)));
     arma::mat dsq(len,2,arma::fill::zeros);
     dsq.each_col() = arma::sqrt(dim);
     arma::mat Dsq = dsq.t();
     arma::mat P1(2,leng);
     for(int i=0;i<leng;i++){
        P1.col(i) = Breakbase((int)q(i),Dsq.col(i));
        }
     arma::vec p(2);arma::mat PQ = P1.t();
     p.row(0) = Base2Digit(PQ.col(0), arma::sqrt(dim));
     p.row(1) = Base2Digit(PQ.col(1), arma::sqrt(dim));
     int fl = (int)Base2Digit(p,Dim);
     return(fl);
}
  

arma::uvec Reg2DM(arma::vec dim){
           int len = dim.n_rows;
           int ran = (int)prodd(dim,1,len);
           arma::uvec ff(ran);
           for(int i=0;i<ran;i++){
               ff(i) = FLay(i,dim);
           }
           return(ff);
}


arma::cx_vec DiReg2DM3(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); 
     for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                       int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                       int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                       int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                       int e1 = g30 + sqrt(S3)*(g20 +sqrt(S2)*g10);
                       int e2 = g31 + sqrt(S3)*(g21 +sqrt(S2)*g11);
                       ABCD(e1*sqrt(S1*S2*S3) + e2) = AA(i3 + S3*(i2 +S2*i1));
                       }
                }
            }
       return ABCD;
}


arma::cx_vec DiReg2DM4(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); double S4 = dim(3);
     for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                       for(int i4 = 0;i4<S4;i4++){
                             int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                             int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                             int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                             int g40 = div(i4,(int)sqrt(S4)).quot; int g41 = i4 - g40*sqrt(S4);
                             int e1 = g40 + sqrt(S4)*(g30 + sqrt(S3)*(g20 +sqrt(S2)*g10));
                             int e2 = g41 + sqrt(S4)*(g31 + sqrt(S3)*(g21 +sqrt(S2)*g11));
                             ABCD(e1*sqrt(S1*S2*S3*S4) + e2) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                       }
                }
            }
       }
      return ABCD;
}


arma::cx_vec DiReg2DM5(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); double S4 = dim(3); double S5 = dim(4);
     for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                     for(int i4 = 0;i4<S4;i4++){
                          for(int i5 = 0;i5<S5;i5++){
                             int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                             int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                             int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                             int g40 = div(i4,(int)sqrt(S4)).quot; int g41 = i4 - g40*sqrt(S4);
                             int g50 = div(i5,sqrt(S5)).quot; int g51 = i5 - g50*sqrt(S5);
                             int e1 = g50 + sqrt(S5)*(g40 + sqrt(S4)*(g30 + sqrt(S3)*(g20 +sqrt(S2)*g10)));
                             int e2 = g51 + sqrt(S5)*(g41 + sqrt(S4)*(g31 + sqrt(S3)*(g21 +sqrt(S2)*g11)));
                             ABCD(e1*sqrt(S1*S2*S3*S4*S5) + e2) = AA(i5 + S5*(i4 + S4*(i3 + S3*(i2 +S2*i1))));
                          }
                     }
                }
            }
       }
      return ABCD;
}


arma::cx_vec DiReg2DM6(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); double S4 = dim(3); double S5 = dim(4);double S6 = dim(5);
     for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                     for(int i4 = 0;i4<S4;i4++){
                         for(int i5 = 0;i5<S5;i5++){
                           for(int i6 = 0;i6<S6;i6++){
                             int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                             int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                             int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                             int g40 = div(i4,(int)sqrt(S4)).quot; int g41 = i4 - g40*sqrt(S4);
                             int g50 = div(i5,(int)sqrt(S5)).quot; int g51 = i5 - g50*sqrt(S5);
                             int g60 = div(i6,(int)sqrt(S6)).quot; int g61 = i6 - g60*sqrt(S6);
                             int e1 = g60 + sqrt(S6)*(g50 + sqrt(S5)*(g40 + sqrt(S4)*(g30 + sqrt(S3)*(g20 +sqrt(S2)*g10))));
                             int e2 = g61 + sqrt(S6)*(g51 + sqrt(S5)*(g41 + sqrt(S4)*(g31 + sqrt(S3)*(g21 +sqrt(S2)*g11))));
                             ABCD(e1*sqrt(S1*S2*S3*S4*S5*S6) + e2) = AA(i6 + S6*(i5 + S5*(i4 + S4*(i3 + S3*(i2 +S2*i1)))));
                           }
                        }
                     }
                }
            }
       }
      return ABCD;
}


arma::cx_vec DiReg2DM(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     int check = (int)dim.n_rows;
     if(check == 3){
        double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); 
        for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                       int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                       int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                       int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                       int e1 = g30 + sqrt(S3)*(g20 +sqrt(S2)*g10);
                       int e2 = g31 + sqrt(S3)*(g21 +sqrt(S2)*g11);
                       ABCD(e1*sqrt(S1*S2*S3) + e2) = AA(i3 + S3*(i2 +S2*i1));
                }
           }
        }
     }  
     else{
       if(check == 4){
           double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); double S4 = dim(3);
           for(int i1 = 0;i1<S1;i1++){
              for(int i2 = 0;i2<S2;i2++){
                  for(int i3 = 0;i3<S3;i3++){
                       for(int i4 = 0;i4<S4;i4++){
                             int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                             int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                             int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                             int g40 = div(i4,(int)sqrt(S4)).quot; int g41 = i4 - g40*sqrt(S4);
                             int e1 = g40 + sqrt(S4)*(g30 + sqrt(S3)*(g20 +sqrt(S2)*g10));
                             int e2 = g41 + sqrt(S4)*(g31 + sqrt(S3)*(g21 +sqrt(S2)*g11));
                             ABCD(e1*sqrt(S1*S2*S3*S4) + e2) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                       }
                  }
              }
           }
       }
       else{
         if(check == 5){  
            double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); double S4 = dim(3); double S5 = dim(4);
            for(int i1 = 0;i1<S1;i1++){
              for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                   for(int i4 = 0;i4<S4;i4++){
                     for(int i5 = 0;i5<S5;i5++){
                             int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                             int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                             int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                             int g40 = div(i4,(int)sqrt(S4)).quot; int g41 = i4 - g40*sqrt(S4);
                             int g50 = div(i5,(int)sqrt(S5)).quot; int g51 = i5 - g50*sqrt(S5);
                             int e1 = g50 + sqrt(S5)*(g40 + sqrt(S4)*(g30 + sqrt(S3)*(g20 +sqrt(S2)*g10)));
                             int e2 = g51 + sqrt(S5)*(g41 + sqrt(S4)*(g31 + sqrt(S3)*(g21 +sqrt(S2)*g11)));
                             ABCD(e1*sqrt(S1*S2*S3*S4*S5) + e2) = AA(i5 + S5*(i4 + S4*(i3 + S3*(i2 +S2*i1))));
                     }
                   }
                }
              }
            }
         }
         else{
           if(check == 6){
              double S1 = dim(0); double S2 = dim(1); double S3 = dim(2); double S4 = dim(3); double S5 = dim(4);double S6 = dim(5);
              for(int i1 = 0;i1<S1;i1++){
                 for(int i2 = 0;i2<S2;i2++){
                    for(int i3 = 0;i3<S3;i3++){
                       for(int i4 = 0;i4<S4;i4++){
                          for(int i5 = 0;i5<S5;i5++){
                             for(int i6 = 0;i6<S6;i6++){
                             int g10 = div(i1,(int)sqrt(S1)).quot; int g11 = i1 - g10*sqrt(S1);
                             int g20 = div(i2,(int)sqrt(S2)).quot; int g21 = i2 - g20*sqrt(S2);
                             int g30 = div(i3,(int)sqrt(S3)).quot; int g31 = i3 - g30*sqrt(S3);
                             int g40 = div(i4,(int)sqrt(S4)).quot; int g41 = i4 - g40*sqrt(S4);
                             int g50 = div(i5,(int)sqrt(S5)).quot; int g51 = i5 - g50*sqrt(S5);
                             int g60 = div(i6,(int)sqrt(S6)).quot; int g61 = i6 - g60*sqrt(S6);
                             int e1 = g60 + sqrt(S6)*(g50 + sqrt(S5)*(g40 + sqrt(S4)*(g30 + sqrt(S3)*(g20 +sqrt(S2)*g10))));
                             int e2 = g61 + sqrt(S6)*(g51 + sqrt(S5)*(g41 + sqrt(S4)*(g31 + sqrt(S3)*(g21 +sqrt(S2)*g11))));
                             ABCD(e1*sqrt(S1*S2*S3*S4*S5*S6) + e2) = AA(i6 + S6*(i5 + S5*(i4 + S4*(i3 + S3*(i2 +S2*i1)))));
                            }
                          }
                       }
                    }
                 }
              }
           }
           else{
           ABCD.zeros();}}}}
      return ABCD;
}


arma::cx_vec DiDM2Reg(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     int check = (int)dim.n_rows;
     if(check == 3){
           double S1 = dim(0); double S2 = dim(1); double S3 = dim(2);
           for(int i1 = 0;i1<S1;i1++){
              for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                    int num = i3 + S3*(i2 +S2*i1);
                    int e1 = div(num,(int)sqrt(S1*S2*S3)).quot; int e2 = num - e1*sqrt(S1*S2*S3);
                    int g10 = div(e1,(int)sqrt(S2*S3)).quot; int g20 = div(e1-g10*sqrt(S2*S3),(int)sqrt(S3)).quot; int g30 = e1-g10*sqrt(S2*S3)-g20*sqrt(S3);
                    int g11 = div(e2,(int)sqrt(S2*S3)).quot; int g21 = div(e2-g11*sqrt(S2*S3),(int)sqrt(S3)).quot; int g31 = e2-g11*sqrt(S2*S3)-g21*sqrt(S3);
                    int f1 = g10*sqrt(S1)+g11;int f2 = g20*sqrt(S2)+g21;int f3 = g30*sqrt(S3)+g31;
                    ABCD(f3 + S3*(f2 +S2*f1)) = AA(i3 + S3*(i2 +S2*i1));
                }
              }
           }
     }
     else{
       if(check == 2){
          double S1 = dim(0); double S2 = dim(1); 
          for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
                    int num = i2 +S2*i1;
                    int e1 = div(num,(int)sqrt(S1*S2)).quot; int e2 = num - e1*sqrt(S1*S2);
                    int g10 = div(e1,(int)sqrt(S2)).quot; int g20 = e1-g10*sqrt(S2);
                    int g11 = div(e2,(int)sqrt(S2)).quot; int g21 = e2-g11*sqrt(S2);
                    int f1 = g10*sqrt(S1)+g11;int f2 = g20*sqrt(S2)+g21;
                    ABCD(f2 +S2*f1) = AA(i2 +S2*i1);
            }
          }
       }
       else{
           ABCD.zeros();}}
     return ABCD;
}

arma::cx_vec DiDM2Reg3(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     double S1 = dim(0); double S2 = dim(1); double S3 = dim(2);
     for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                    int num = i3 + S3*(i2 +S2*i1);
                    int e1 = div(num,(int)sqrt(S1*S2*S3)).quot; int e2 = num - e1*sqrt(S1*S2*S3);
                    int g10 = div(e1,(int)sqrt(S2*S3)).quot; int g20 = div(e1-g10*sqrt(S2*S3),(int)sqrt(S3)).quot; int g30 = e1-g10*sqrt(S2*S3)-g20*sqrt(S3);
                    int g11 = div(e2,(int)sqrt(S2*S3)).quot; int g21 = div(e2-g11*sqrt(S2*S3),(int)sqrt(S3)).quot; int g31 = e2-g11*sqrt(S2*S3)-g21*sqrt(S3);
                    int f1 = g10*sqrt(S1)+g11;int f2 = g20*sqrt(S2)+g21;int f3 = g30*sqrt(S3)+g31;
                    ABCD(f3 + S3*(f2 +S2*f1)) = AA(i3 + S3*(i2 +S2*i1));
                }
           }
     }
     return ABCD;
}

arma::cx_vec DiDM2Reg2(const arma::cx_vec &AA, const arma::vec &dim){
     arma::cx_vec ABCD(AA.n_rows);
     double S1 = dim(0); double S2 = dim(1); 
     for(int i1 = 0;i1<S1;i1++){
           for(int i2 = 0;i2<S2;i2++){
                    int num = i2 +S2*i1;
                    int e1 = div(num,(int)sqrt(S1*S2)).quot; int e2 = num - e1*sqrt(S1*S2);
                    int g10 = div(e1,(int)sqrt(S2)).quot; int g20 = e1-g10*sqrt(S2);
                    int g11 = div(e2,(int)sqrt(S2)).quot; int g21 = e2-g11*sqrt(S2);
                    int f1 = g10*sqrt(S1)+g11;int f2 = g20*sqrt(S2)+g21;
                    ABCD(f2 +S2*f1) = AA(i2 +S2*i1);
           }
     }
     return ABCD;
}

arma::uvec cSwapref(arma::vec dim, int k, int l){
          int len = dim.n_rows;
          int ran = (int)prodd(dim,1,len);
          arma::vec dimswp = SwapParts(dim,k-1,l-1);
          arma::uvec ff(ran);
          for(int i=0;i<ran;i++){
             arma::vec po = SwapParts(Breakbase(i,dim),k-1,l-1);
             ff(i) = Base2Digit(po,dimswp);
             }             
          return(ff);
}



arma::vec RefSort(arma::vec State, arma::uvec Ref){
          arma::uvec indic = sort_index(Ref);
          arma::vec ref = State.elem(indic);
          return(ref);
} 

arma::cx_vec cRefSort(const arma::cx_vec &AA, const arma::uvec &Ref){
            arma::uvec indic = sort_index(Ref);
            return AA.elem(indic);
}    


arma::mat ULay(arma::mat Uni, arma::vec dim, arma::vec truncdim){
     arma::uvec R1 = sort_index(DM2Reg(dim));
     arma::uvec R2 = sort_index(DM2Reg(truncdim));
     arma::mat uni = Uni.cols(R1);
     uni = uni.rows(R2);
     return(uni);
}

arma::cx_mat cULay(const arma::cx_mat &Uni, const arma::vec &dim, const arma::vec &truncdim){
     arma::uvec R1 = sort_index(DM2Reg(dim));
     arma::uvec R2 = sort_index(DM2Reg(truncdim));
     arma::cx_mat uni = Uni.cols(R1);
     return uni.rows(R2);
}


arma::cx_mat cULay2(const arma::cx_mat &Uni, const arma::vec &dim, const arma::vec &truncdim){
       arma::cx_mat tuni(size(Uni));
       int check = (int)dim.n_rows;
       if(check == 3){
           double S1 = dim(0); double S2 = dim(1); double S3 = dim(2);
           for(int i1 = 0;i1<S1;i1++){
              for(int i2 = 0;i2<S2;i2++){
                for(int i3 = 0;i3<S3;i3++){
                    int num = i3 + S3*(i2 +S2*i1);
                    int e1 = div(num,(int)sqrt(S1*S2*S3)).quot; int e2 = num - e1*sqrt(S1*S2*S3);
                    int g10 = div(e1,(int)sqrt(S2*S3)).quot; int g20 = div(e1-g10*sqrt(S2*S3),(int)sqrt(S3)).quot; int g30 = e1-g10*sqrt(S2*S3)-g20*sqrt(S3);
                    int g11 = div(e2,(int)sqrt(S2*S3)).quot; int g21 = div(e2-g11*sqrt(S2*S3),(int)sqrt(S3)).quot; int g31 = e2-g11*sqrt(S2*S3)-g21*sqrt(S3);
                    int f1 = g10*sqrt(S1)+g11;int f2 = g20*sqrt(S2)+g21;int f3 = g30*sqrt(S3)+g31;
                    tuni.col(f3 + S3*(f2 +S2*f1)) = Uni.col(i3 + S3*(i2 +S2*i1));
                }
              }
           }
       }
       else{
        if(check == 2){
          double S1 = dim(0); double S2 = dim(1); 
          for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
                    int num = i2 +S2*i1;
                    int e1 = div(num,(int)sqrt(S1*S2)).quot; int e2 = num - e1*sqrt(S1*S2);
                    int g10 = div(e1,(int)sqrt(S2)).quot; int g20 = e1-g10*sqrt(S2);
                    int g11 = div(e2,(int)sqrt(S2)).quot; int g21 = e2-g11*sqrt(S2);
                    int f1 = g10*sqrt(S1)+g11;int f2 = g20*sqrt(S2)+g21;
                    tuni.col(f2 +S2*f1) = Uni.col(i2 +S2*i1);
            }
          }
        }
        else{
           tuni.zeros(arma::size(Uni));}}
     return tuni;
}


std::complex<double> InnProd(arma::cx_vec A,arma::cx_vec B){
                     arma::cx_vec Bt = conj(B);
                     std::complex<double> ip = std::inner_product(std::begin(A), std::end(A), std::begin(Bt), std::complex<double>());
                     return (ip);
}

arma::vec SwapSite(arma::vec A, arma::vec dim, int k, int l){
          int len = A.n_rows;
          arma::vec Sw(len,arma::fill::zeros);
          for(int i=0;i<len;i++){
             arma::vec po = SwapParts(Breakbase(i,dim),k,l);
             int co = Base2Digit(po,SwapParts(dim,k,l));
             Sw(co) = A(i);
             }             
          return(Sw);
}

arma::cx_vec cSwapSite(arma::cx_vec A, arma::vec dim, int k, int l){
          int len = A.n_rows;
          arma::cx_vec Sw(len,arma::fill::zeros);
          for(int i=0;i< len;i++){
             arma::vec po = SwapParts(Breakbase(i,dim),k-1,l-1);
             int co = Base2Digit(po,SwapParts(dim,k-1,l-1));
             Sw(co) = A(i);
             }             
          return(Sw);
}


arma::mat PartTrace(const arma::vec &State, arma::vec dim){
             arma::mat Int = Col2rho(State);
             int kd = (int)sqrt(dim(0));
             int che = (int)sqrt(dim(1));
             arma::mat TT(kd,kd);
             for(int i=0;i<kd;i++){
               for(int j=0;j<kd;j++){
                 TT(i,j) = arma::trace(Int.submat((i*che),(j*che),(i*che)+che-1,(j*che)+che-1));
               }
             }  
             return TT;
}


arma::cx_mat cPartTrace(const arma::cx_vec &State, double dim0, double dim1){
             arma::cx_mat Int = cCol2rho(State);
             int kd = (int)sqrt(dim0);
             int che = (int)sqrt(dim1);
             arma::cx_mat TT(kd,kd);
             for(int i=0;i<kd;i++){
               for(int j=0;j<kd;j++){
                 TT(i,j) = arma::trace(Int.submat((i*che),(j*che),(i*che)+che-1,(j*che)+che-1));
               }
             }  
             return TT;
}
         
arma::cx_mat cPartTrace2(const arma::cx_vec &A, double dim0, double dim1){
             int kd = (int)sqrt(dim0);
             int che = (int)sqrt(dim1);
             arma::cx_mat Int = arma::reshape(A,kd*che,kd*che).st();
             arma::cx_mat TT(kd,kd);
             for(int i=0;i<kd;i++){
               for(int j=0;j<kd;j++){
                 TT(i,j) = arma::trace(Int.submat((i*che),(j*che),(i*che)+che-1,(j*che)+che-1));
               }
             }  
             return TT;
}

arma::cx_mat cPartTraceD(const arma::cx_vec &A, double dim0, double dim1){
             int kd = (int)dim0;
             int che = (int)dim1;
             arma::cx_mat Int = arma::reshape(A,che,kd).st();
             arma::cx_mat TT(kd,kd);
             for(int i=0;i<kd;i++){
               for(int j=0;j<kd;j++){
                 TT(i,j) = cdot(Int.row(j),Int.row(i));
               }
             }  
             return TT;
}

arma::cx_mat cPartTraceE(const arma::cx_vec &A, double dim0, double dim1){
             int kd = (int)dim0;
             int che = (int)dim1;
             arma::cx_mat Int = arma::reshape(A,che,kd).st();
             return Int*Int.t();
}

arma::cx_mat cPartTraceX(const arma::cx_vec &A, double dim0, double dim1){
             int kd = (int)dim0;
             int che = (int)dim1;
             arma::cx_mat TT(kd,kd);
             for(int i=0;i<kd;i++){
               arma::cx_vec Ape =  A.subvec(i*che,i*che+che-1);
               TT(i,i) = cdot(Ape,Ape);
               for(int j=i+1;j<kd;j++){
                 TT(i,j) = cdot(A.subvec(j*che,j*che+che-1),Ape);
                 TT(j,i) = std::conj(TT(i,j));
               }
             }  
             return TT;
}

arma::vec Join(const std::vector< arma::vec > &A){
          int len = A.size();
          arma::vec joi = A.at(0);
          for(int i=1;i<len;i++){
             joi = join_cols(joi,A.at(i));
             }
          return(joi);
}   


arma::mat KP(const std::vector<arma::mat> &A){
             int len = A.size();
             arma::mat kp = A.at(0);
             for(int i=1;i<len;i++){
             kp = arma::kron(kp,A.at(i));
                }
             return(kp);
}

arma::cx_mat cKP(const std::vector<arma::cx_mat> &A){
             int len = A.size();
             arma::cx_mat kp = A.at(0);
             for(int i=1;i<len;i++){
             kp = arma::kron(kp,A.at(i));
                }
             return(kp);
}


arma::cx_vec cFoursiteSwap(const arma::cx_vec &AA, const arma::vec dim, const arma::uvec res){
         arma::cx_vec ABC(AA.n_rows);
             for(int i1 = 0;i1<dim(0);i1++){
               for(int i2 = 0;i2<dim(1);i2++){
                for(int i3 = 0;i3<dim(2);i3++){
                 for(int i4 = 0;i4<dim(3);i4++){
                    int ex[5] = {0,i1,i2,i3,i4};
                    ABC(ex[res(3)] + dim(res(3)-1)*(ex[res(2)] + dim(res(2)-1)*(ex[res(1)] + dim(res(1)-1)*ex[res(0)]))) = AA(i4 + dim(3)*(i3 + dim(2)*(i2 + dim(1)*i1)));
                 }
                }
               }
             }
         return ABC;
}


arma::cx_vec cFoursiteSwapexact13(const arma::cx_vec &AA, double S1, double S2, double S3, double S4){
         arma::cx_vec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  for(int i4 = 0;i4<S4;i4++){
                      ABC(i4 + S4*(i1 + S1*(i2 +S2*i3))) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                 }
                }
               }
             }
         return ABC;
}

arma::cx_vec cFoursiteSwapexact23(const arma::cx_vec &AA, double S1, double S2, double S3, double S4){
         arma::cx_vec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  for(int i4 = 0;i4<S4;i4++){
                      ABC(i4 + S4*(i2 + S2*(i3 +S3*i1))) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                 }
                }
               }
             }
         return ABC;
}

arma::uvec uvFoursiteSwapexact23(const arma::uvec &AA, double S1, double S2, double S3, double S4){
         arma::uvec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  for(int i4 = 0;i4<S4;i4++){
                      ABC(i4 + S4*(i2 + S2*(i3 +S3*i1))) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                 }
                }
               }
             }
         return ABC;
}


arma::uvec uvFoursiteSwapexact13(const arma::uvec &AA, double S1, double S2, double S3, double S4){
         arma::uvec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  for(int i4 = 0;i4<S4;i4++){
                      ABC(i4 + S4*(i1 + S1*(i2 +S2*i3))) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                 }
                }
               }
             }
         return ABC;
}



arma::cx_vec cFoursiteSwapexact12(const arma::cx_vec &AA, double S1, double S2, double S3, double S4){
         arma::cx_vec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  for(int i4 = 0;i4<S4;i4++){
                      ABC(i4 + S4*(i3 + S3*(i1 +S1*i2))) = AA(i4 + S4*(i3 + S3*(i2 +S2*i1)));
                 }
                }
               }
             }
         return ABC;
}




arma::cx_vec cThreesiteSwapexact13(const arma::cx_vec &AA, double S1, double S2, double S3){
         arma::cx_vec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  ABC(i1 + S1*(i2 +S2*i3)) = AA(i3 + S3*(i2 +S2*i1));
                 }
                }
               }
         return ABC;
}


arma::cx_vec cThreesiteSwapexact12(const arma::cx_vec &AA, double S1, double S2, double S3){
         arma::cx_vec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  ABC(i3 + S3*(i1 +S1*i2)) = AA(i3 + S3*(i2 +S2*i1));
                 }
                }
               }
         return ABC;
}


arma::uvec uvThreesiteSwapexact12(const arma::uvec &AA, double S1, double S2, double S3){
         arma::uvec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  ABC(i3 + S3*(i1 +S1*i2)) = AA(i3 + S3*(i2 +S2*i1));
                 }
                }
               }
         return ABC;
}


arma::cx_vec cThreesiteSwapexact23(const arma::cx_vec &AA, double S1, double S2, double S3){
         arma::cx_vec ABC(AA.n_rows);
           for(int i1 = 0;i1<S1;i1++){
            for(int i2 = 0;i2<S2;i2++){
               for(int i3 = 0;i3<S3;i3++){
                  ABC(i2 + S2*(i3 +S3*i1)) = AA(i3 + S3*(i2 +S2*i1));
                 }
                }
               }
         return ABC;
}


double Concurrence(const arma::cx_mat Rho){
       arma::cx_mat csy0(2,2,arma::fill::zeros);
       csy0(0,1) = {0.0,-1.0};
       csy0(1,0) = {0.0,1.0};
       arma::cx_mat RhoP = Rho*arma::kron(csy0,csy0)*arma::conj(Rho)*arma::kron(csy0,csy0);
       arma::vec eigval1;
       arma::cx_mat eigvec1;
       eig_sym(eigval1,eigvec1,RhoP);
       arma::vec Srt = arma::sort(arma::sqrt(eigval1),"descend");
       double Con = std::max(0.0,Srt(0)-Srt(1)-Srt(2)-Srt(3));
       return Con;
}

//-----------------------------------------------------------------------------------------------------------------------------------

/*arma::mat ULay(arma::mat Uni, arma::vec dim, arma::vec truncdim){
     arma::uvec R1 = DM2Reg(dim);
     arma::uvec R2 = DM2Reg(truncdim);
     arma::mat uni = Uni;
     int len = dim.n_rows;
     int led = truncdim.n_rows;
     for(int i=0;i<len ;i++){
            uni.col(i) = RefSort(uni.col(i),R1);
            }
     arma::mat tuni = uni.st();
     for(int i=0;i<led;i++){
           tuni.col(i) = RefSort(tuni.col(i),R2);
           }
     return(tuni.st());
}

    
arma::cx_mat cULay(arma::cx_mat Uni, arma::vec dim, arma::vec truncdim){
     arma::uvec R1 = DM2Reg(dim);
     arma::uvec R2 = DM2Reg(truncdim);
     arma::cx_mat uni = Uni;
     int len = R1.n_rows;
     int led = R2.n_rows;
            for(int i=0;i<len ;i++){
            uni.col(i) = cRefSort(uni.col(i),R1);
            }
     arma::cx_mat tuni = uni.st();
           for(int i=0;i<led;i++){
           tuni.col(i) = cRefSort(tuni.col(i),R2);
           }
     return(tuni.st());


arma::mat PartTrace(arma::vec State, arma::vec kpdim, arma::vec tracedim){
             int len = State.n_rows;
             int keepdim = (int)prodd(kpdim,1,kpdim.n_rows);
             int t1 = (int)tracedim(0);
             int t2 = (int)tracedim(1);
             arma::cube Ageis(len,1,1);
             Ageis.slice(0) = State;
             arma::cube Re = reshape(Ageis,t2,t1,keepdim);
             Ageis.clear();
             arma::vec Ro(keepdim);
             int i = 0;
             while(i<keepdim){
                double summ = 0;int j =0;
                  while(j<t1){
                   int k =0;
                     while(k<t2){
                     summ = summ + Re(k,j,i);k=k+std::sqrt(t2)+1;
                     }
                  j=j+std::sqrt(t1)+1;
                  }
             Ro(i) = summ;i++;
             }
             Re.clear();
             arma::uvec R1 = Reg2DM(kpdim);
             arma::vec Rope = RefSort(Ro,R1);
             arma::mat PT = Col2rho(Rope);
             R1.clear();
             Rope.clear();
             return(PT);
}

arma::cx_mat cPartTrace(arma::cx_vec State, arma::vec kpdim, arma::vec tracedim){
             int len = State.n_rows;
             int keepdim = (int)prodd(kpdim,1,kpdim.n_rows);
             int t1 = (int)tracedim(0);
             int t2 = (int)tracedim(1);
             arma::cx_cube Ageis(len,1,1);
             Ageis.slice(0) = State;
             arma::cx_cube Re = reshape(Ageis,t2,t1,keepdim);
             Ageis.clear();
             arma::cx_vec Ro(keepdim);
             int i = 0;
             while(i<keepdim){
                std::complex<double> summ = {0,0};int j =0;
                  while(j<t1){
                   int k =0;
                     while(k<t2){
                     summ = summ + Re(k,j,i);k=k+std::sqrt(t2)+1;
                     }
                  j=j+std::sqrt(t1)+1;
                  }
             Ro(i) = summ;i++;
             }
             Re.clear();
             arma::uvec R1 = Reg2DM(kpdim);
             arma::cx_vec Rope = cRefSort(Ro,R1);
             arma::cx_mat PT = cCol2rho(Rope);
             R1.clear();
             Rope.clear();
             return(PT);
}
             
}*/
