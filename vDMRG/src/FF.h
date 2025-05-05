//EDITED FUNCTIONS FOR OPTIMIZATION

#ifndef ADD_HEADER
#define ADD_HEADER
 
#include <iostream>
#include <armadillo>
#include <math.h>
#include <algorithm>
#include <numeric>

// Header file for functions_file.cc


arma::cx_mat UTT(arma::cx_mat A, double tstep);

arma::cx_cube UT(const arma::cx_cube &Le, double tstep);

arma::vec SwapParts(arma::vec A, int k,int l);

double prodd(arma::vec A, int k, int l);
   

arma::vec Breakbase(int num, const arma::vec &A);

arma::vec Breakbase2(int num, double S1, double S2);


int Base2Digit(const arma::vec &B, const arma::vec &dim);

int Base2Digit2(const arma::vec &B, const arma::vec &dim);

arma::mat Col2rho(arma::vec A);


//arma::cx_mat cCol2rho2(arma::cx_vec A);

arma::cx_mat cCol2rho(const arma::cx_vec &A);

arma::vec Rho2col(arma::mat A);

arma::cx_vec cRho2col(const arma::cx_mat &A);



arma::vec Mapleft(arma::vec A, arma::mat U, arma::vec dim);

arma::cx_vec cMapleft(arma::cx_vec A, arma::cx_mat U, arma::vec dim);

arma::cx_vec ctMapleft(const arma::cx_vec &A, const arma::cx_mat &U, double dim0, double dim1);

arma::vec Mapright(arma::vec A, arma::mat U, arma::vec dim);

arma::cx_vec cMapright(arma::cx_vec A, arma::cx_mat U, arma::vec dim);

arma::cx_vec ctMapright(const arma::cx_vec &A, const arma::cx_mat &U, double dim0, double dim1);

arma::vec Traceop(double si);

arma::vec Join(const std::vector< arma::vec > &A);

arma::cx_mat KP(const std::vector<arma::cx_mat> &A);

arma::cx_mat cKP(const std::vector<arma::cx_mat> &A);

double Ttr(arma::vec A, arma::vec dim);

double cTtr(const arma::cx_vec &A, const arma::vec &dim);


int BLay(int Num, arma::vec dim);

int FLay(int Num, arma::vec dim);
  

arma::uvec Reg2DM(arma::vec dim);

arma::uvec DM2Reg(arma::vec dim);

arma::cx_vec DiReg2DM3(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiReg2DM4(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiReg2DM5(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiReg2DM6(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiReg2DM(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiDM2Reg(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiDM2Reg2(const arma::cx_vec &AA, const arma::vec &dim);

arma::cx_vec DiDM2Reg3(const arma::cx_vec &AA, const arma::vec &dim);



arma::vec RefSort(arma::vec State, arma::uvec Ref);
   
arma::cx_vec cRefSort(const arma::cx_vec &AA, const arma::uvec &Ref);


arma::mat ULay(arma::mat Uni, arma::vec dim, arma::vec truncdim);
    
arma::cx_mat cULay(const arma::cx_mat &Uni, const arma::vec &dim, const arma::vec &truncdim);

arma::cx_mat cULay2(const arma::cx_mat &Uni, const arma::vec &dim, const arma::vec &truncdim);

std::complex<double> InnProd(arma::cx_vec A,arma::cx_vec B);

arma::vec SwapSite(arma::vec A, arma::vec dim, int k, int l);

arma::cx_vec cSwapSite(arma::cx_vec A, arma::vec dim, int k, int l);

arma::cx_vec cFoursiteSwap(const arma::cx_vec &A, arma::vec dim, arma::uvec res);

arma::cx_vec cFoursiteSwapexact13(const arma::cx_vec &AA, double S1, double S2, double S3, double S4);

arma::cx_vec cFoursiteSwapexact23(const arma::cx_vec &AA, double S1, double S2, double S3, double S4);

arma::uvec uvFoursiteSwapexact13(const arma::uvec &AA, double S1, double S2, double S3, double S4);

arma::uvec uvFoursiteSwapexact23(const arma::uvec &AA, double S1, double S2, double S3, double S4);

arma::cx_vec cFoursiteSwapexact12(const arma::cx_vec &AA, double S1, double S2, double S3, double S4);

arma::cx_vec cThreesiteSwapexact13(const arma::cx_vec &AA, double S1, double S2, double S3);

arma::cx_vec cThreesiteSwapexact12(const arma::cx_vec &AA, double S1, double S2, double S3);

arma::uvec uvThreesiteSwapexact12(const arma::uvec &AA, double S1, double S2, double S3);

arma::cx_vec cThreesiteSwapexact23(const arma::cx_vec &AA, double S1, double S2, double S3);

arma::uvec cSwapref(arma::vec dim, int k, int l);

arma::mat PartTrace(arma::vec State, arma::vec dim);

arma::cx_mat cPartTrace(const arma::cx_vec &State, double dim0, double dim1);

arma::cx_mat cPartTrace2(const arma::cx_vec &A, double dim0, double dim1);

arma::cx_mat cPartTraceD(const arma::cx_vec &A, double dim0, double dim1);

arma::cx_mat cPartTraceE(const arma::cx_vec &A, double dim0, double dim1);

arma::cx_mat cPartTraceX(const arma::cx_vec &A, double dim0, double dim1);

double Concurrence(const arma::cx_mat Rho);

#endif //ADD_HEADER

