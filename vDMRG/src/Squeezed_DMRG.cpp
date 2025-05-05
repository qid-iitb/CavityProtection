/* Version dated: 2025/04/30 
Original code written by: Himadri S. Dhar
Minor modifications done by: Siddharth Tiwary and Harsh Sharma
--------------------------------------------------------------------------------|
compile: g++ --std=c++11 Squeezed_DMRG.cpp FF.cpp -O3 -larmadillo
    run: ./a.out < ParameterFile.txt
--------------------------------------------------------------------------------|
Purpose: Calculate the full dynamics of a mesoscopic number of spins (Nspin)
    with frequency distribution given by ww(Nspin) strongly coupled to a single
    cavity mode. We consider strong but not ultra-strong coupling using a rota-
    ting wave approximation - counter rotating terms are neglected. In addition
    we assume a Markovian environment resulting in a Lindblad-type evolution.
    For details of the model see ref.[1]
         The cavity mode might be driven by some constant external pulse
    (amplitude: eta0 or eta1, frequency: wp, duration: Tpulse). Typically the
    spins are unexcited and cavity is empty or in some desired initial state.

Files:
    - "***_DMRG.cpp"          ... main program
    - "ParameterFile.txt"     ... parameters loaded to main program
    - "FF.cpp"                ... most important self-defined functions
    - "FF.h"                  ... headers

Basic nomenclature and procedure of DMRG:
    the Nspin spins are subdivided into 2 blocks:
       - 'A' (sometimes also called 'sys' or old notation 'L')
       - and 'B' ('env' or old notation 'R')
    The cavity ('cav') is always treated exactly. An additional free spin which
    is not included in 'A' or 'B' is denoted by 'free'.
    Density operator is vectorized and denoted as Rh***
    RhABC -> total density matrix containing 'A'+'B'+'cav'
    sometimes in the documentation we will explicitly mark the position of the
    free spin by 'f'
    RhABfC -> density matrix containing 'A'+'B'+'free'+'cav'
    RhAfBC -> density matrix containing 'A'+'free'+'B'+'cav'

    Time evolution matrix UT_**
    UT_fC -> time evolution is always applied to 1 spin ("free spin") + cavity;
    the free spin was chosen to be on rightmost position for this matrix multi-
    plication to be efficient -> spin released from 'A' or 'B' has to be swapp-
    ed with mostlikely 'B' to be next to 'cav'.

Notes:
for details of the methods see:
[1] main and Suppl of ... https://doi.org/10.1103/PhysRevLett.121.133601
[2] armadillo library ... http://arma.sourceforge.net/docs.html

------------------------------------------------------------------------------*/
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>  
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <string>
#include "FF.h"
using arma::kron; // -> to use kron instead of arma::kron since arma::kron is
                  // needed more than 50 times this simplfies reading

//------------------------------------------------------------------------------

inline arma::mat eyed(int n)           // identity matrix(nxn)
{
    arma::mat Identity(n,n,arma::fill::eye);
    return Identity;
}

inline arma::mat nullM(int n)          // matrix(nxn) containing zeros only
{
    arma::mat NullM(n,n,arma::fill::zeros);
    return NullM;
}

inline arma::cx_mat cxnullM(int n)     // same as nullM but complex valued
{
    arma::cx_mat cxNullM(nullM(n),nullM(n));
    return cxNullM;
}

inline arma::cx_mat ceyed(int n)       // same as eyed but complex valued
{
    arma::cx_mat cId(eyed(n),nullM(n));
    return cId;
}

inline arma::mat Aop(int n)            // annhiliation operator a
{
    arma::mat aop(n,n,arma::fill::zeros);
    for(int i = 0; i < n-1; i++){;
        aop(i,i+1)= sqrt(i+1);
    }
    return aop;
}

inline arma::cx_mat cAop(int n)        // same as Aop but complex valued
{
    arma::cx_mat aop(n,n,arma::fill::zeros);
    for(int i = 0; i < n-1; i++){;
        aop(i,i+1)= {sqrt(i+1),0};
    }
    return aop;
}

inline arma::mat Adag(int n)           // number operator a^dag a
{
    arma::mat adaga(n,n,arma::fill::zeros);
    for(int i = 0; i < n; i++){;
        adaga(i,i)= i;
    }
    return adaga;
}

inline arma::cx_mat cAdag(int n)      // same as Adag but complex valued
{
    arma::cx_mat adaga(n,n,arma::fill::zeros);
    for(int i = 0; i < n; i++){;
        adaga(i,i)= i;
    }
    return adaga;
}

inline arma::cx_mat cDisplace(int n,std::complex<double> alpha)
{
    /* Single-mode displacement operator.   exp(alpha a^dag-alpha^* a)
    ------------------------------------------------------------------
    n ... number of fock states
    alpha ... displacement amplitude
    usage:
    arma::cx_mat displace1=cDisplace(Ncav,{Re(alpha),Im(alpha)});
    psi_cav1 = displace1*psi_cav;
    //... where psi_cav is the cav.wavefct. (vec(Ncav)) to be displaced
    -----------------------------------------------------------------*/
    arma::cx_mat D = arma::expmat( alpha*arma::trans(cAop(n))
                                  - std::conj(alpha)*cAop(n) ) ;

    return D;
}

void write2tofile(std::ofstream &outputfile,double nr1, double nr2){
    // in principle here I should include a check if this file exists
    outputfile.precision(10);
    outputfile<<std::scientific;
    outputfile.width(20);
    outputfile << nr1;
    outputfile.width(20);
    outputfile << nr2 <<std::endl;
}

void write4tofile(std::ofstream &outputfile,double nr1, double nr2,
                   double nr3, double nr4){
    // in principle here I should include a check if this file exists
    outputfile.precision(10);
    outputfile<<std::scientific;
    outputfile.width(20);
    outputfile << nr1;
    outputfile.width(20);
    outputfile << nr2;
    outputfile.width(20);
    outputfile << nr3;
    outputfile.width(20);
    outputfile << nr4 <<std::endl;
}

void save2file(std::ofstream &outputfile, int M, arma::vec outvec){
    // write the whole vector of dim M into one line of outputfile
    for (int j=0;j<M;j++){;
        outputfile.precision(7);
        outputfile<<std::scientific;
        outputfile.width(17);
        outputfile << outvec(j);
    }
    outputfile <<std::endl;
}
//------------------------------------------------------------------------------
//---------------------MAIN CODE------------------------------------------------

int main(int argc, char *argv[])
{
    // Loading data from Parameter File to an array
    std::vector<double> data;
    data.reserve(11);    // avoid reallocations

    //-----------Read line‑by‑line from stdin -----------------------------------
    std::string line;
    while (std::getline(std::cin, line) && data.size() < 11) {
        // Strip inline or whole‑line comments beginning with ‘#’ 
        const auto hash = line.find('#');
        if (hash != std::string::npos)
            line.resize(hash);        // keep only what’s before '#'

        if (line.empty()) continue;

        std::stringstream ss(line);
        double value;
        while (ss >> value) {         
            data.push_back(value);
            if (data.size() == 11) break;
        }
    }
    
    for(int start_val = 0; start_val < 1;start_val++){   // Start (L999)

    //--------Input number of spins and cavity: --------------------------------
    
    int Nspin = data[0];                 // Total number of spins (1st col. of *)
    int Ncav = data[1];                  // Fock states in the cavity
    int bondDim = data[2];               // Par_Mat(start_val,3); // "Bond dimension" gives the max number
                                         // of eigenvec. kept in renorm. subspace
    double Nspi = (double)Nspin;
    int st_val = data[10];

    //--------------Defining the parameter space: ------------------------------

    double rval = data[9];        // Squeezing parameter
    double wc = data[3];          // Cavity frequency in rad. per second
    double ws = 1/cosh(2*rval);   // freq. of central spin in units of wc
    double wp = 0.0;              // Driving frequency in units of wc
    double ka = data[4]*wc/wc;    // Cavity decay (L(a))
    double gp = data[5]*wc/wc;    // Non-radiative dephasing (L(sz))
    double gh = data[6]*wc/wc;    // Spin decay, radiative dephasing (L(sm))
    double std = data[8]/wc;      // std.dev. of coupling strengths

    double Om_in = data[7]/(sqrt(Nspin)*wc);       // Coupling strength for each spin
    
    std::complex<double> eta0(tanh(2*rval)*(wc/wc-wp/2),0.00);  // squeezed driving
    double Et = eta0.real();

    //--------------Setting up the spin distribution: --------------------------
    
     arma::vec ww(Nspin);                // gaussian distribution just used for individual spins

     std::default_random_engine gen;
     std::normal_distribution<double> dist(ws,std);
    
     for (int i=0; i<Nspin; ++i) {
         ww[i] = dist(gen);
     }

    std::ofstream spin_distr_file;           // store spin distr. to file
    spin_distr_file.open("Spin_distribution_"+std::to_string(start_val)+".dat");
    spin_distr_file<<(ww)<<std::endl<<std::endl;
    spin_distr_file.close();

    //--------------Defining the coupling strengths: ---------------------------
    
    arma::vec g_distr(Nspin);                // distribution of coupling
    g_distr.ones();                          // uniform coupling
    arma::vec g = Om_in*g_distr;             // distr. of coupling strengths g(i)

    double Om = sqrt(arma::sum(arma::square(g))); // collective coupl. strength

    //--------------Defining the program parameters: ---------------------------
    double Trabi = 2.0*M_PI/(Om);                        // Rabiperiod given by coll.coupl.
    double Tfinal = 6.0*Trabi*std::exp(-rval)*2;         // final time in multipl. of Trabi
    int number_of_Tsteps = 10000*std::exp(-rval)*2;      // nr. of total timesteps (noTs)
    double disc_step = number_of_Tsteps;                 // ( nr. of gridpoints = noTs)
    double tstep = Tfinal/disc_step;                     // size of one timestep




    // define vectors for time-gridpoints and pulsesequence, the latter is used
    // to determine at each time-gridpoint the eta value for the next timestep
    // via int values: 0 (nodrive), 1 (+eta), and -1 (-eta)
    arma::vec time_gridpoints = arma::regspace(0,tstep,Tfinal);

//--Write used parameters to "Parameter_file**.dat" **start_val: ---------------
    std::ofstream pfile;
    pfile.open("Parameter_file"+std::to_string(start_val)+".dat");
    pfile<<"Number of spins:  "<<Nspin<<std::endl;
    pfile<<"Number of Fock states in the cavity:  "<<Ncav<<std::endl;
    pfile<<"Driving frequency offset: "<<(ws-wp)<<std::endl;
    pfile<<"Squeezing parameter: "<<rval<<std::endl;
    pfile<<"Squeezed driving strength: "<<eta0<<std::endl<<std::endl;
    pfile<<"Non-radiative dephasing (sig_z): "<<gp<<std::endl;
    pfile<<"Spin decay (sig_min): "<<gh<<std::endl;
    pfile<<"Cavity decay: "<<ka<<std::endl<<std::endl;
    pfile<<"Standard deviation of the distribution: "<<std<<std::endl;
    pfile<<"Spin coupling for central frequency: "<<Om_in<<std::endl;
    pfile<<"Spin coupling per spin: "<<std::endl<<g<<std::endl<<std::endl;
    pfile<<"Total coupling:  "<<Om<<std::endl<<std::endl;
    pfile<<"Bond dimension chosen:  "<<bondDim<<std::endl;
    pfile<<"Total time steps:  "<<disc_step<<std::endl;
    pfile<<"Superposition of "<<st_val<<" and "<<st_val+1<<std::endl;
    pfile<<"---------------Allocation done------------------------"<<std::endl;
//------------------------------------------------------------------------------
    std::cout<<"Trabi: "<<Trabi<<std::endl;

    //--------------Defining the operators: ------------------------------------
    // c... labels complex matrices or vectors
    // nullM(N) is a selfdefined NxN matrix containing only zeros

    const std::complex<double> ic(0.0,1.0);   // c-number "i"
    const std::complex<double> re(1.0,0.0);
    arma::cx_mat csy0(2,2,arma::fill::zeros); // sigma_y matrix
    csy0(0,1) = {0.0,-1.0};
    csy0(1,0) = {0.0,1.0};

    arma::mat sx0(2,2,arma::fill::zeros);     // sigma_x matrix
    sx0(0,1) = 1.0;
    sx0(1,0) = 1.0;

    arma::mat sz0(2,2,arma::fill::zeros);     // sigma_z matrix
    sz0(0,0) = 1.0;
    sz0(1,1) = -1.0;

    arma::mat sp0(2,2, arma::fill::zeros);    // sigma_+ matrix
    sp0(0,1) = 1.0;

    arma::cx_mat csx0(sx0,nullM(2));          // create complex valued versions
    arma::cx_mat csz0(sz0,nullM(2));
    arma::cx_mat csp0(sp0,nullM(2));
    arma::cx_mat csm0(sp0.t(),nullM(2));      // sigma_- matrix


    //--------------Defining the Hamiltonian: ----------------------------------

    arma::cx_mat Hcav0 = (1.0/Nspi)*((1.0/cosh(2*rval))*kron(eyed(2),Adag(Ncav)))*re;

    arma::cube Hspin(2*Ncav, 2*Ncav, Nspin); // spin and interaction Hamiltonian
    for(int i=0; i < Nspin; i++){;           // for each spin arranged as 'cube'
        Hspin.slice(i) = (1.0/2.0)*(ww[i])*kron(sz0,eyed(Ncav))
                        + g[i]*std::exp(rval)*(kron(sp0.t()+sp0,Aop(Ncav).t()+Aop(Ncav)))/2
            + g[i]*std::exp(-rval)*(kron(sp0.t()-sp0,Aop(Ncav).t()-Aop(Ncav)))/2;
    }

    arma::cx_cube HH0(2*Ncav, 2*Ncav, Nspin);  // tot. Hamiltonian for eta
    for(int i=0; i < Nspin; i++){;
       HH0.slice(i) = Hcav0+Hspin.slice(i);
    }

    //--------------Defining the Lindbladian operators -------------------------
    // the Lindbladian will act 'locally' only at one spin+cavity =>
    // dim of Hilbertsp. (2*Ncav)^2 in vectorized space this becomes (2*Ncav)^4
    // note that (C^dag)^T = C^*; kron(A,B)^T=kron(A^T,B^T)

    // first only dissipators are defined:

    arma::mat Lcav(pow(2*Ncav,2),pow(2*Ncav,2));
    Lcav = (1.0*ka/Nspi)*(  kron( kron(eyed(2),Aop(Ncav)) , kron(eyed(2),conj(Aop(Ncav))) )
          -0.5*kron( kron(eyed(2),Adag(Ncav)) , kron(eyed(2),eyed(Ncav)) )
          -0.5*kron( kron(eyed(2),eyed(Ncav)) , kron(eyed(2),Adag(Ncav)) )  );

    arma::mat Lspin(pow(2*Ncav,2),pow(2*Ncav,2));
    Lspin = (gh)*(  kron( kron(conj(sp0.t()),eyed(Ncav)) , kron(sp0.t(),eyed(Ncav)) )
           -0.5*kron( kron(sp0*conj(sp0.t()),eyed(Ncav)) , kron(eyed(2),eyed(Ncav)) )
           -0.5*kron( kron(eyed(2),eyed(Ncav)) , kron(conj(sp0)*sp0.t(),eyed(Ncav)) )  )
           +(gp)*(  kron( kron(sz0,eyed(Ncav)) , kron(conj(sz0),eyed(Ncav)) )
                  - kron( kron(eyed(2),eyed(Ncav)) , kron(eyed(2),eyed(Ncav)) )  );

    // now construct whole Lindbladian:      

    arma::cx_cube LL0(pow(2*Ncav,2),pow(2*Ncav,2), Nspin);
    for(int i=0; i < Nspin; i++){;
        LL0.slice(i) = -ic*(  kron(HH0.slice(i),eyed(2*Ncav))
                            - kron(eyed(2*Ncav),conj(HH0.slice(i)))  ) + Lcav + Lspin;
    }

    // Note that the Lindbladians are vectorized in the combined basis vec(rho)
    // => also UT (and ch1) will be in combined basis;  However, we want to
    // perform all our manipulations on individual parts "A","B","free",or "cav"
    // => we have to change from the combined vec. basis to a individually vec-
    // torized basis vec("A"), vec("B"), ... this is done by cULay (see FF.cpp)

    arma::cx_cube ch1 = UT(LL0,tstep/2.0); // ch1: UT ... time evolution matrix
    arma::cx_cube UT_fC0;
    UT_fC0.copy_size(ch1);
    for(int i = 0; i < ch1.n_slices;i++){
        UT_fC0.slice(i) = cULay(ch1.slice(i),{4,pow(Ncav,2)},{4,pow(Ncav,2)});
        // -> UT_fC ... time evolution matrix now in indiv. vectorized basis
        //              acting on "free" and "cav" (=>_fC)
    }
    ch1.clear();

//------ Infinite DMRG to build the initial system -----------------------------
// This is in principle not needed as long as the initial state is trivial but
// let's keep it this way: start with 1 spin in "A", 1 spin in "B" and the
// the cavity "cav" in some state psi_cav (here cat state)

    arma::cx_vec psi_f = {{0.0,0.0},{1.0,0.0}};  // unexc. spin (0,1)   [c-val.]
    arma::cx_vec Rhin = cRho2col(psi_f * psi_f.t()); // vectorize(rho=ket*bra)
    arma::cx_vec Rsys = Rhin;                        // start with 1 spin in
    arma::cx_vec Renv = Rhin;                        // in "sys" and "env"
    arma::cx_vec Rh_f = Rhin;                        // Rh_f ... is "free"

    arma::cx_vec psi_cav(Ncav,arma::fill::zeros);
    psi_cav(st_val)={1.0/std::sqrt(2),0.0};                        // cav. state: |1>
    psi_cav(st_val+1)={1.0/std::sqrt(2),0.0};                      // cav. state: |2>
    
    arma::cx_vec Rcav = cRho2col(psi_cav * psi_cav.t());  // vect.(rho=ket*bra)
    psi_cav.clear();
    psi_f.clear();
    Rhin.clear();

    // std::cout<<Rcav<<std::endl;

    // dimensions of "sys" (A), "env" (B), "cav" (c), and "free" (f)
    double Ss = Rsys.n_rows;
    double Ee = Renv.n_rows;
    double Cc = Rcav.n_rows;
    double Si = Rh_f.n_rows;

    arma::cx_vec RhABC = kron(kron(Rsys,Renv),Rcav); // RhABC total density op-
    // erator 'A'+'B'+'cav' in individually vectorized form, now filled with one
    // spin in 'A'('sys') and 'B'('env') => '1'+'1'+'cav'

    // Renormalization matrices: 1 free spin is added to the renormalized sub-
    // space by the matrix Renorm* (maps a (d+2)x(d+2) matrix to a dxd matrix)
    // RenormA ... put free spin from right to "A";  Rh(Af)BC -> Rh(A)BC
    // RenormB ... put free spin from left to "B";   RhA(fB)C -> RhA(B)C
    // transpose .t() ... releases spin from the named subspace "A" or "B"
    arma::cx_mat* RenormA = new arma::cx_mat[Nspin];
    arma::cx_mat* RenormB = new arma::cx_mat[Nspin];

    std::complex<double> Tra;      // value of the trace of RhABC

    // computing the trace in the vectorized space is like a vec multiplication
    arma::cx_vec TTA[Nspin+1];    // trace vec for "A"; start indexation with 1
    TTA[1] = cRho2col(ceyed(2));  // first only one spin will be in "A" ("sys")
    arma::cx_vec TTf = cRho2col(ceyed(2));    // trace vec for "free"
    arma::cx_vec TTB[Nspin+1];    // trace vec for "B"; start indexation with 1
    TTB[Nspin] = cRho2col(ceyed(2));          // last spin will be in "B"("env")
    arma::cx_vec TTC = cRho2col(ceyed(Ncav)); // trace vec for "cav"

    int exactspin = 1; // number of exact spins you start with (keep 1 or change
    int startspin = 1+exactspin;              // the code for the sweep)
    int endspin = Nspin-exactspin;

    // just one sweep to set up initial state: goal is to have Nspin-1 spins in
    // "A", 1 spin in "B", and the "cav"
    for(int i = startspin;i < endspin+1; i++){  // not best to write in c++
        arma::cx_vec RhABC1 = kron(RhABC,Rh_f); // add one spin to right of the
        arma::cx_vec TT = kron(TTA[i-1],TTf);   // total system -> RhABCf and
                                                // increase TTA by one spin

        // now swap 2nd and 3rd site of 3 listed sites with listed dimensions,
        // so that the added free spin is on right place (between "A" and "B"):
        RhABC1 = cThreesiteSwapexact23(RhABC1,Ss,Ee*Cc,Si); // RhABCf -> RhAfBC

        // get reduced density vec of "A"+"free" (in order to renormalize it):
        arma::cx_mat red = cPartTraceE(RhABC1,Ss*Si,Ee*Cc);

        // get renormalization matrix:
        arma::vec eigval;
        arma::cx_mat eigvec;
        eig_sym(eigval,eigvec,red);
        arma::uvec ind = arma::sort_index(eigval,"descend"); // sorted indices
        arma::cx_mat Uni = eigvec.cols(ind);                 // of max eigenvals
        arma::vec ref = eigval.elem(ind);
        int tru = bondDim;
        if(bondDim > Uni.n_rows){
            tru = Uni.n_rows;}
        Uni = Uni.cols(0,tru-1); // keep only relevant rows that correspond to
                                 // max eigenvals

        double error = arma::sum(real(ref)) - arma::sum(real(ref.rows(0,tru-1)));

        RenormA[i] = Uni.t(); // renorm. matrix for putting i-th spin to "A"
        // inverse will release the i-th spin from "A"

        // Renormalize: get from Rh(Af)BC to RhABC where f is now in A
        RhABC = ctMapleft(RhABC1,RenormA[i],Ss*Si,Ee*Cc); // RhABC1 contains the
              // free spin -> ctMapleft multiplies RenormA[i] to left part of
              // RhABC1 which is of dim Ss*Si to map it to the renormalized dim
              // (this way is more efficient than creating large matrix and
              //  multiply whole RhABC1 because right part remains unchanged)

        TTA[i] = RenormA[i]*TT; // update trace operator! because in renormal-
                                // ized space trace is not the trivial vec(diag)

        Ss = (double)tru;       // update dim of "A" because of extra spin
                                // once bonddim is reached this will not change
        Tra = arma::cdot(cKP({TTA[i],TTB[Nspin],TTC}),RhABC);
                                // check if trace of rho = 1
    }
    // now spin_1 to spin_N-1 are in left part "A" ("sys")

//------------------------------ Output files ----------------------------------

    //1) "Error_File_*.dat" will contain: t, truncation error, trace of rho,
    std::ofstream error_file;           // occupation of highest Fockstate
    error_file.open("Error_File_"+std::to_string(start_val)+".dat");

    //2) "AdagA_*.dat"     will contain: t, <AdagA> (, Re(<A>), Im(<A>))
    std::ofstream adaga_file;
    adaga_file.open("AdagA_"+std::to_string(start_val)+".dat");

    //3) "Cav_State_Re_*.dat" will contain: vectorized Real(rho_cav)
    std::ofstream Rcav_Re_file;
    Rcav_Re_file.open("Cav_State_Re_"+std::to_string(start_val)+".dat");

    //4) "Cav_State_Im_*.dat" will contain: vectorized Imag(rho_cav)
    std::ofstream Rcav_Im_file;
    Rcav_Im_file.open("Cav_State_Im_"+std::to_string(start_val)+".dat");

//------------------------------ Time-adaptive DMRG  ---------------------------
// the time evolution is performed in a Suzuki-Totter way acting always on the
// cavity and one free spin only; to do this efficient the density vector has
// to be of the type RhABfC; apply time evolution for timestep/2 -> sweep twice
// through whole ensemble to perform one time step (left sweep; right sweep):

    // indexvector to keep track of the sweep [N,N-1,...,2,1,1,2,3,...,N]
    arma::uvec Ispin = arma::linspace<arma::uvec>(Nspin,1,Nspin);
    Ispin = join_cols(Ispin,arma::linspace<arma::uvec>(1,Nspin,Nspin));

    arma::cx_mat cav(Ncav,Ncav);       // matrix to store the cav density matrix
    arma::cx_vec sz = cRho2col(csz0.t()); // bras to construct the expect. value
    arma::cx_vec sx = cRho2col(csx0.t());
    arma::cx_vec sy = cRho2col(csy0.t());
    arma::cx_vec sp = cRho2col(csp0.t());
    arma::cx_vec sm = cRho2col(csm0.t());
    arma::cx_vec adag = cRho2col(cAdag(Ncav).st());

    arma::cx_cube UT_fC = UT_fC0;      // start with eta
    std::cout<<"initialization complete for cycle...."<<start_val<<std::endl;



    // time evolution is applied only to fC of RhABfC => we always have to swap
    arma::uvec Base0[Ispin.n_rows];  // the swapping is performed by swapping
    arma::uvec Base1[Ispin.n_rows];  // the indices of the rho vec. this is done
    arma::uvec Base2[Ispin.n_rows];  // by base (since this permutations do not
    arma::uvec MBase0[Ispin.n_rows]; // change we calc. base only for 1. time
    arma::uvec MBase1[Ispin.n_rows]; // step and save the result)

    // calculate observables for t=0: - - - - - - - - - - CALC OBSERVABLES - - -
    arma::cx_vec Tt1 = cKP({TTA[Nspin-1],TTB[Nspin]});

    arma::cx_vec par = cRho2col(cxnullM(Ncav));
    for(int i=0;i<Ncav;i++){
        for(int j=0;j<Ncav;j++){
            par(i*Ncav+j) = {1.0,0.0};
            cav(i,j) = arma::cdot(cKP({Tt1,par}),RhABC);
            par(i*Ncav+j) = {0.0,0.0};
        }
    }

    std::complex<double> Navg2(0.0,0.0);

    for(int i = 0;i<Ncav;i++){
        std::complex<double> ne((double)i,0.0);
        std::complex<double> ne2(std::sqrt((double)i),0.0);  // for <a>
        if(i>0){Navg2 = Navg2 + ne*cav(i,i);}
    }

    double Fock_err = std::abs(cav(Ncav-1,Ncav-1));  // occ. final fockstate
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // write values for t=0 to files before we do first time step:
    // (1)
    write4tofile(error_file,0.0,0.0,0.0,Fock_err);      // "Error_File_..."
    // (2)
    write2tofile(adaga_file,0.0,std::abs(Navg2));       // "AdagA_*.dat"
    // (3)
    arma::vec Rcav_Re_vec=real(cRho2col(cav));
    save2file(Rcav_Re_file,Ncav*Ncav,Rcav_Re_vec);      // "Cav_State_Re..."
    // (4)
    arma::vec Rcav_Im_vec=imag(cRho2col(cav));
    save2file(Rcav_Im_file,Ncav*Ncav,Rcav_Im_vec);      // "Cav_State_Im..."

    int count = 0;                     // set counter for start
//--------------------------------------------------------- Start of (TIME LOOP)
    for(int t_ind = 0; t_ind < number_of_Tsteps; t_ind++){
        double t = time_gridpoints(t_ind+1);          // time at next gridpoint
        std::cout<< count <<std::endl;       // print progress to std-output
        count = count + 1;
        // Et = eta0.real();

        double error = 0;           // to calculate truncation error
        double le = Ispin.n_rows;

        // each timestep does double sweep through ensemble: -------------------
        for(int i = 0;i<le;i++){                             // Start of (SWEEP)

            unsigned ispin = Ispin(i); // contains index of the current spin
                                       // start with Nth spin
            bool mark = true;   // defines if left (true) or right (false) sweep
            if(i >= (int)le/2){ // left sweep is over
                mark = false;
            } // (this might only work for exactspin=1)

            if(ispin == Nspin){ // start with Nth spin -> this one is in "B" and
                                // is already "free" since at the start of each
                                // timestep there is just 1 spin in "B"
                RhABC = ctMapright(RhABC, UT_fC.slice(ispin-1),Ss,Si*Cc);
                                // note that index of UT starts with 0  => -1
                Tra = arma::cdot(cKP({TTA[ispin-1],TTB[ispin],TTC}),RhABC);
                // cKP self defined complex kronecker product to compute trace
                RhABC = RhABC/Tra;
             }
            else{               // if it is not the Nth spin:
                if(ispin > exactspin){ // but also not the 1st spin:

                    if(mark == true){  // LEFT SWEEP: (from right to left)
                      // release next spin from "A": Rh(A)BC -> Rh(Af)BC
                      RhABC = ctMapleft(RhABC,RenormA[ispin].t(),Ss,Ee*Cc);
                      Tra = arma::cdot(cKP({TTA[ispin-1],TTf,TTB[ispin+1],TTC}),RhABC);
                      RhABC = RhABC/Tra;
                      Ss = RenormA[ispin].n_cols/4; // update dim of "A"

                      if(count == 1){ // Base does the index swap but is the
                                      // same for all timesteps => needs to be
                                      // calculated only at first timestep
                          Base0[i] = uvFoursiteSwapexact23(arma::linspace<arma::uvec>(0,Ss*Si*Ee*Cc-1,Ss*Si*Ee*Cc),Ss,Si,Ee,Cc);
                      }

                      // rearange density vec indices such that RhAfBC -> RhABfC
                      RhABC = RhABC.elem(Base0[i]);
                      // apply time evolution to "free" and "cav"
                      RhABC = ctMapright(RhABC, UT_fC.slice(ispin-1),Ee*Ss,Si*Cc);

                      if(count == 1){ // get the permutation for the swap "free"
                          // with "A" (1,3); needs to be calculated only at firt
                          // timestep
                          Base1[i] = uvFoursiteSwapexact13(arma::linspace<arma::uvec>(0,Ss*Si*Ee*Cc-1,Ss*Si*Ee*Cc),Ss,Ee,Si,Cc);
                      }
                      RhABC = RhABC.elem(Base1[i]); // swap: RhABfC -> RhfBAC
                      // [this is done so that we can use mapLeft to renorm "B"
                      //  mapRight would act also on cavity (not what we want)]

                      // get reduced density vec of "free"+"B" (for renorm.):
                      arma::cx_mat red = cPartTraceE(RhABC,Si*Ee,Ss*Cc);

                      // get renormalization matrix:
                      arma::vec eigval;
                      arma::cx_mat eigvec;
                      eig_sym(eigval,eigvec,red);
                      arma::uvec ind = sort_index(eigval,"descend"); // max eig.
                      arma::cx_mat Uni = eigvec.cols(ind);           // vals.
                      arma::vec ref = eigval.elem(ind);
                      int tru = Uni.n_rows;
                      if(bondDim < Uni.n_rows){
                         tru = bondDim;
                      }
                      Uni = Uni.cols(0,tru-1); // keep only relevant rows that
                      // correspond to max eigenvals

                      error = error + (arma::sum(ref) - arma::sum(ref.rows(0,tru-1)));

                      RenormB[ispin] = Uni.t(); // renorm. matrix for putting
                      // i-th spin to "B"; inverse will release the i-th spin

                      // Renorm.: get from Rh(fB)AC to Rh(B)AC where f is now in "B"
                      RhABC = ctMapleft(RhABC, RenormB[ispin],Ee*Si,Ss*Cc);

                      // update trace operations
                      TTB[ispin]= RenormB[ispin]*kron(TTf,TTB[ispin+1]);
                      Ee = (double)tru;    // update dimensions

                      if(count == 1){ // get the permutation to swap back to
                                      // "normal" order;
                          Base2[i] = uvThreesiteSwapexact12(arma::linspace<arma::uvec>(0,Ss*Ee*Cc-1,Ss*Ee*Cc),Ee,Ss,Cc);
                      }

                      // swap back to "normal" order "A","B","C":
                      RhABC = RhABC.elem(Base2[i]); // RhBAC -> RhABC
                      Tra = arma::cdot(cKP({TTA[ispin-1],TTB[ispin],TTC}),RhABC);
                      RhABC = RhABC/Tra;
                  }
                    else{   // RIGHT SWEEP: (from left to right)
                            // after left sweep all spins but 1 is in "B"

                      if(count == 1){ // get perm matrix for all following times
                          Base0[i] = uvThreesiteSwapexact12(arma::linspace<arma::uvec>(0,Ss*Ee*Cc-1,Ss*Ee*Cc),Ss,Ee,Cc);
                      }

                      RhABC = RhABC.elem(Base0[i]); // sweep RhABC -> RhBAC
                      // [this is done so that we can use mapLeft to release "B"
                      //  mapRight would act also on cavity (not what we want)]

                      // Release next spin from "B": Rh(B)AC -> Rh(fB)AC
                      RhABC = ctMapleft(RhABC,RenormB[ispin].t(),Ee,Ss*Cc);

                      Ee = RenormB[ispin].n_cols/4;
                      Tra = arma::cdot(cKP({TTf,TTB[ispin+1],TTA[ispin-1],TTC}),RhABC);
                      RhABC = RhABC/Tra;

                      if(count == 1){ // get perm matrix for all following times
                          Base1[i] = uvFoursiteSwapexact13(arma::linspace<arma::uvec>(0,Ss*Si*Ee*Cc-1,Ss*Si*Ee*Cc),Si,Ee,Ss,Cc);
                      }

                      // rearange density vec indices such that RhfBAC -> RhABfC
                      RhABC = RhABC.elem(Base1[i]);
                      // apply time evolution to "free" and "cav"
                      RhABC = ctMapright(RhABC,UT_fC.slice(ispin-1),Ee*Ss,Si*Cc);

                      if(count == 1){ // get perm matrix for all following times
                          Base2[i] = uvFoursiteSwapexact23(arma::linspace<arma::uvec>(0,Ss*Si*Ee*Cc-1,Ss*Si*Ee*Cc),Ss,Ee,Si,Cc);
                      }
                      // rearange density vec indices such that RhABfC -> RhAfBC
                      RhABC = RhABC.elem(Base2[i]);

                      // now f can be put back into "A"
                      // get renormalization matrix:
                      arma::cx_mat red = cPartTraceE(RhABC,Ss*Si,Ee*Cc);
                      arma::vec eigval;
                      arma::cx_mat eigvec;
                      eig_sym(eigval,eigvec,red);
                      arma::uvec ind = sort_index(eigval,"descend");
                      arma::cx_mat Uni = eigvec.cols(ind);
                      arma::vec ref = eigval.elem(ind);
                      int tru = Uni.n_rows;
                      if(bondDim < Uni.n_rows){
                         tru = bondDim;
                      }
                      Uni = Uni.cols(0,tru-1); // keep only relevant rows of
                                               // max eigenvals

                      error = error + (arma::sum(ref)-arma::sum(ref.rows(0,tru-1)));

                      RenormA[ispin] = Uni.t(); // renorm. matrix for putting
                      //  i-th spin to "A"; inverse to release the i-th spin

                      // Renorm.: get from Rh(Af)BC to RhABC where f is now in A
                      RhABC = ctMapleft(RhABC, RenormA[ispin], Ss*Si,Ee*Cc);

                      TTA[ispin]=RenormA[ispin]*kron(TTA[ispin-1],TTf);
                      Ss = (double)tru;
                      Tra = arma::cdot(cKP({TTA[ispin],TTB[ispin+1],TTC}),RhABC);
                      RhABC = RhABC/Tra;
                  }
                }
                else{   // for 1st spin -> this is the last part of the left
                        // sweep and first part of right sweep; in both cases
                        // now there is only 1 spin in "A" and N-1 spins in "B"

                      if(count == 1){ // get perm matrix for all following times
                          Base0[i] = uvThreesiteSwapexact12(arma::linspace<arma::uvec>(0,Si*Ee*Cc-1,Si*Ee*Cc),Si,Ee,Cc);
                      }

                      // rearange density vec indices such that RhABC -> RhBAC
                      RhABC = RhABC.elem(Base0[i]);

                      // since there is only 1 spin in "A" this can be treated
                      // exactly as "free" -> apply time evolution to "A"+"cav"
                      RhABC = ctMapright(RhABC,UT_fC.slice(ispin-1),Ee,Si*Cc);

                      if(count == 1){ // get perm matrix for all following times
                          Base1[i] = uvThreesiteSwapexact12(arma::linspace<arma::uvec>(0,Si*Ee*Cc-1,Si*Ee*Cc),Ee,Si,Cc);
                      }
                      // rearange density vec indices back to "normal" order
                      RhABC = RhABC.elem(Base1[i]);        // RhBAC -> RhABC
                      Tra = arma::cdot(cKP({TTA[ispin],TTB[ispin+1],TTC}),RhABC);RhABC = RhABC/Tra;

                }
            }
        } //----------------------------------------------------- End of (SWEEP)
          // now one timestep is complete so let's calculate
          // some observables: - - - - - - - - - - - - - CALC. OBSERVABLES - - -

        arma::cx_vec Tt1 = cKP({TTA[Nspin-1],TTB[Nspin]});

        arma::cx_vec par = cRho2col(cxnullM(Ncav));
        for(int i=0;i<Ncav;i++){
          for(int j=0;j<Ncav;j++){
              par(i*Ncav+j) = {1.0,0.0};
              cav(i,j) = arma::cdot(cKP({Tt1,par}),RhABC);
              par(i*Ncav+j) = {0.0,0.0};
          }
        }

        std::complex<double> Navg2(0.0,0.0);

        for(int i = 0;i<Ncav;i++){
            std::complex<double> ne((double)i,0.0);
            std::complex<double> ne2(std::sqrt((double)i),0.0);  // for <a>
            if(i>0){Navg2 = Navg2 + ne*cav(i,i);}

        }

        double Fock_err = std::abs(cav(Ncav-1,Ncav-1));  // occ. final fockstate
        double TraRe = std::real(Tra);                   // real part of trace
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // at the end of each time step write results to files:
        // (1)
        write4tofile(error_file,t/wc*1e7,error,TraRe,Fock_err); // "Error_File_..."
        // (2)
        write2tofile(adaga_file,t/wc*1e7,std::abs(Navg2));      // "AdagA_*.dat"
        // (3)
        arma::vec Rcav_Re_vec=real(cRho2col(cav));
        save2file(Rcav_Re_file,Ncav*Ncav,Rcav_Re_vec);      // "Cav_State_Re..."
        // (4)
        arma::vec Rcav_Im_vec=imag(cRho2col(cav));
        save2file(Rcav_Im_file,Ncav*Ncav,Rcav_Im_vec);      // "Cav_State_Im..."

    }  //---------------------------------------------------- End of (TIME LOOP)

    error_file.close();     // (1)
    adaga_file.close();     // (2)
    Rcav_Re_file.close();   // (3)
    Rcav_Im_file.close();   // (4)

    delete []RenormA;
    delete []RenormB;
    }                                                      // End of (L999)

return 0;
}
//------------------------------------------------------------------------------
//----------------------------- End of progam ----------------------------------
