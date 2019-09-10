/*
#include "defines.hpp"
#include "sptr.hpp"
#include "exception.hpp"
#include "utils.hpp"
#include "math/stats/normaldistribution.hpp"
#include "math/matrix.hpp"
#include "math/interpol/piecewisepolynomial.hpp"
#include "market/yieldcurve.hpp"
#include "math/random/rng.hpp"
#include "methods/pde/tridiagonalops1d.hpp"
#include "methods/pde/pdegrid.hpp"
#include "methods/pde/pdeparams.hpp"
*/
#include "pricer/simplepricers.hpp"
#include "market/market.hpp"
#include "math/optim/polyfunc.hpp"
#include "math/optim/roots.hpp"
#include "math/stats/meanvarcalculator.hpp"
#include "products/europeancallput.hpp"
#include "products/asianbasketcallput.hpp"
#include "math/linalg/linalg.hpp"
#include "methods/montecarlo/eulerpathgenerator.hpp"
#include "methods/montecarlo/antitheticpathgenerator.hpp"
#include "pricer/bsmcpricer.hpp"
#include "pricer/multiassetbsmcpricer.hpp"
#include "methods/pde/pdebase.hpp"
#include "methods/pde/pde1dsolver.hpp"
#include "products/americancallput.hpp"

#include <vector>

#include <iostream>

int main() {
  /** Test for PDE and americancallput */
  try {
    Vector tMat{0.08, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5};
    double tExp = 5;
    double r = 0.02;
    // Vector spotRate{0.01, 0.012, 0.013, 0.014, 0.015, 0.02, 0.025, 0.03};
    Vector fwdVol{0.1, 0.15, 0.22, 0.18, 0.2, 0.19, 0.21, 0.23, 0.25};

    auto SptrYsp =  make_shared<YieldCurve>(&tExp, &tExp + 1, &r, &r + 1,
                                            YieldCurve::InputType::SPOTRATE);
    auto SptrVsp =  make_shared<VolatilityTermStructure>(tMat.cbegin(), tMat.cend(),
                                                         fwdVol.cbegin(), fwdVol.cend(),
                                                         VolatilityTermStructure::VolType::SPOTVOL);

    int payoffType = -1;
    double strike = 100, spot = 100, divyield = 0, T = 1, vol = 0.2;

    cout << "EuroCallPut:" << endl;
    auto prodEuro = SPtrProduct(new EuropeanCallPut(payoffType, strike, T));

    // cout << prodEuro->nAssets() << endl;
    Pde1DResults results;
    bool storeAllReulsts = false;
    PdeParams params{50, 50, 4.0, 1.0};
   
    Pde1DSolver pdeEuro(prodEuro, SptrYsp, spot, divyield, SptrVsp, results, storeAllReulsts);
    pdeEuro.solve(params);

    cout << results.prices[0] << endl;
    cout << prodEuro->payTimes() << endl;
    cout << prodEuro->payAmounts() << endl;
    
    cout << "AmericanCallPut:" << endl;
    Pde1DResults resultsAm;
    auto prodAm = SPtrProduct(new AmericanCallPut(payoffType, strike, T));
    // cout << prodAm->payTimes() << endl;
    // cout << prodAm->payAmounts() << endl;
    Pde1DSolver pdeAm(prodAm, SptrYsp, spot, divyield, SptrVsp, resultsAm, storeAllReulsts);
    
    pdeAm.solve(params);

    cout << resultsAm.prices[0] << endl;
    // cout << prodAm->payTimes() << endl;
    // cout << prodAm->payAmounts() << endl;
    
    
  }
  catch(Exception& e) {
    cout << e.what() << endl;
  }

  /** Test for pdegrid 
  try {
    PdeParams pdeparam;
    GridAxis grax;
    grax.coordinateChange->init(pdeparam);
    double S1 = exp(2);
    double X1 = grax.coordinateChange->fromRealToDiffused(S1);
    cout << X1 << endl;
    double S2 = grax.coordinateChange->fromDiffusedToReal(X1);
    cout << S2 << endl;

    double fwd = 10, vol = 0.4, T = 1;
    grax.coordinateChange->forwardAndVariance(fwd, vol, T);
    cout << fwd << endl;

    double Xmin, Xmax, drift, variance, finalVol;
    grax.coordinateChange->bounds(X1,
                                  fwd,
                                  vol,
                                  T,
                                  pdeparam.nStdDevs[0],
                                  Xmin,
                                  Xmax);   
    cout << Xmin << " to " << Xmax << endl;

    grax.coordinateChange->driftAndVariance(S1,
                                            fwd,
                                            grax.coordinateChange->theta,
                                            0.1,
                                            0.3,
                                            0.3,
                                            0.3,
                                            drift,
                                            variance,
                                            finalVol);
                                  
    cout << "drift = " << drift << endl;
    cout << "variance = " << variance << endl;
    cout << "finalVol = " << finalVol << endl;
  }
  catch(Exception& e) {
    cout << e.what() << endl;
  }
  */
  
  /** Test for tridiagonalops1d 
  try {
    size_t N = 10;
    double DT = 3, DX = 2, theta = 1;
    Vector drifts(N), variances(N);
    drifts.fill(1);
    variances.fill(1);
    
    IdentityOp1D<> identity(N);

    DeltaOp1D<> delta(drifts, DT, DX, theta);

    GammaOp1D<> gamma(variances, DT, DX, theta);

    Vector x(N + 2), y(N + 2), z(N + 2);
    for (size_t i = 0; i < N + 2; ++i) {
      x(i) = i;
      y(i) = 0;
      z(i) = -1;
    }

    
    // identity *= 2;
    // identity += delta;

    auto Op1D = delta - identity;
   
    // Op1D.adjustStandardBoundaryConditions(DX);

    cout << x << endl;
    Op1D.apply(x, y);
    cout << y << endl;
    Op1D.applyInverse(y, z);
    cout << z << endl;
   
  }
  catch(Exception& e) {
    cout << e.what() << endl;
  }
  */

  /** Test for multiassetbsmcpricer with antitheticpathgenerator  
      try {
      int payoffType = 1;
      double strike = 700;
      Vector fixingTimes = {1,2,3,4,5};
      size_t nassets = 6;
      Vector assetQuantities = {2, 3, 2, 1, 1, 1};
      auto prodAsian = SPtrProduct(new AsianBasketCallPut(payoffType, strike, fixingTimes, assetQuantities));

      McParams mcparams {McParams::UrngType::MT19937,
      McParams::PathGenType::EULER, 
      McParams::ControlVarType::NONE};

      Vector spots = {800, 600, 300, 600, 700, 800};
      Vector vols = {0.4, 0.3, 0.2, 0.1, 0.1, 0.25};
      Vector divylds(nassets, arma::fill::zeros);
      vector<double> tMat{1, 2, 3, 4, 5};
      vector<double> spotRate{0.01, 0.01, 0.01, 0.01, 0.01};
      auto SptrYsp =  make_shared<YieldCurve>(tMat.cbegin(), tMat.cend(),
      spotRate.cbegin(), spotRate.cend(),
      YieldCurve::InputType::SPOTRATE);

      MeanVarCalculator<double*> mvcal(1);
      arma::arma_rng::set_seed_random();
      Matrix A(nassets, nassets, arma::fill::randu);
      Matrix corrMat = A.t()*A;
      for (size_t i = 0; i < nassets; ++i)
      corrMat(i, i) =1;
      cout << corrMat << endl;

      unsigned long npaths;
      cout << "npaths = " << endl;
      cin >> npaths;
      for (int ifan = 0; ifan < 2; ++ifan) {
      // cout << "If antithetic? 1: yes; 0: no" << endl;
      // cin >> ifan;
      if (ifan) {
      mcparams.controlVarType = McParams::ControlVarType::ANTITHETIC;
      cout << "Antithetic results:" << endl;
      }
      else
      cout  << "Normal results:" << endl;
        
      MultiAssetBsMcPricer bsmcpricerAsian(prodAsian, SptrYsp, vols, spots,
      divylds, corrMat, mcparams);
      bsmcpricerAsian.simulate(mvcal, npaths);

      auto res = mvcal.results();
      cout << "Price = " << res(0, 0) << endl;
      cout << "Stderr = " << sqrt(res(1, 0) / npaths) << endl;
      }
        
      }
      catch(Exception& e) {
      cout << e.what() << endl;
      }
  */
     
  

  /** Test for linalg
      try {
      arma::arma_rng::set_seed_random();
      Matrix A(5, 5, arma::fill::randu);
      Matrix B = A.t()*A;
      for (size_t i = 0; i < 5; ++i)
      B(i, i) =1;
      cout << B << endl;
        Vector eigval;
        Matrix eigvec;
        eigensym(B, eigval, eigvec);
        cout << eigval << endl;
        cout << eigvec << endl;
        

        spectrunc(B);
        cout << B << endl;
        Matrix L;
        eigensym(B, eigval, eigvec);
        choldcmp(B, L);

        cout << eigval << endl;
        cout << eigvec << endl;
        cout << L << endl;
        
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */

    /** Test for bsmcpricer and with antitheticpathgenerator 
    try {
        int payoffType = 1;
        double strike = 100, T = 1;
        auto prodEuro = SPtrProduct(new EuropeanCallPut(payoffType, strike, T));
        
        McParams mcparams {McParams::UrngType::MT19937,
            McParams::PathGenType::EULER,
            // McParams::ControlVarType::ANTITHETIC};  
            McParams::ControlVarType::NONE};  
        double divYield = 0;
        double spot = 100;
        
        Vector tMat{0.08, 0.25, 0.5, 0.75, 1, 2, 3, 5};
        double rate = 0.2;
        double vol = 0.4;
        Vector spotRate{0.01, 0.012, 0.013, 0.014, 0.015, 0.02, 0.025, 0.03};
        Vector fwdVol(8);
        // Vector spotRate(5), fwdVol(5);
        // spotRate.fill(rate);
        fwdVol.fill(vol);

        auto SptrYsp =  make_shared<YieldCurve>(tMat.cbegin(), tMat.cend(),
                                                spotRate.cbegin(), spotRate.cend(),
                                                YieldCurve::InputType::SPOTRATE);
        auto SptrVsp =  make_shared<VolatilityTermStructure>(tMat.cbegin(), tMat.cend(),
                                                             fwdVol.cbegin(), fwdVol.cend(),
                                                             VolatilityTermStructure::VolType::SPOTVOL);
        MeanVarCalculator<double*> mvcal(1);
        int ifan;
        cout << "If antithetic? 1: yes; 0: no" << endl;
        cin >> ifan;
        if (ifan)
          mcparams.controlVarType = McParams::ControlVarType::ANTITHETIC;
        BsMcPricer bsmcpricerEuro(prodEuro, SptrYsp, divYield, SptrVsp, spot, mcparams);
        
        unsigned long npaths;
        cout << "npaths = " << endl;
        cin >> npaths;
        bsmcpricerEuro.simulate(mvcal, npaths);

        auto res = mvcal.results();
        cout << "Price = " << res(0, 0) << endl;
        cout << "Stderr = " << sqrt(res(1, 0) / npaths) << endl;
        // cout << mvcal.results() << endl;
        // cout << europeanOptionBS(payoffType, spot, strike, T, rate, divYield, vol) << endl;

    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    
    

    /** // Test AsianBasketCallPut
    try {
        size_t ntimesteps = 5, nfactors = 3;
        EulerPathGenerator<NormalRngMinStdRand> pathgen(ntimesteps, nfactors);
        Matrix A(ntimesteps, nfactors);
        pathgen.next(A);
        A += 10;
        cout << A << endl;
        
        int payoffType = 1;
        double strike = 2;
        Vector fixingTimes = {1,2,3,4,5};
        Vector assetQuantities = {2, 3, 2, 1, 1};
        AsianBasketCallPut ac(payoffType, strike, fixingTimes, assetQuantities);

        ac.eval(A);

        cout << ac.nAssets() << endl;
        cout << ac.payTimes() << endl;
        cout << ac.payAmounts() << endl;
        
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */

    
    /** Test EuropeanCallPut  
    try {
        size_t nfactors = 3;
        Vector timesteps = {1, 2, 3};
        arma::arma_rng::set_seed_random();
        Matrix A(nfactors, nfactors, arma::fill::randu);
        Matrix B = A.t()*A;
        for (size_t i = 0; i < nfactors; ++i)
            B(i, i) =1;

        // cout << B << endl;

        spectrunc(B);
        cout << B << endl;
        Matrix L;
        choldcmp(B, L);
        cout << L << endl;
        
        // EulerPathGenerator<NormalRngMinStdRand> pathgen(timesteps.begin(), timesteps.end(), nfactors, B);
        // EulerPathGenerator<NormalRngMinStdRand> pathgen(timesteps.begin(), timesteps.end(), nfactors);
        
        EulerPathGenerator<NormalRngSobol> pathgen(timesteps.begin(), timesteps.end(), nfactors, B);
        // EulerPathGenerator<NormalRngSobol> pathgen(timesteps.begin(), timesteps.end(), nfactors);      
        
        Matrix C(timesteps.size(), nfactors);
        pathgen.next(C);
        // C += 10;
        cout << C << endl;
        
        int payoffType = 1;
        double strike = 2, T = 10;
        EuropeanCallPut ec(payoffType, strike, T);
     
        ec.eval(A);

        cout << ec.nAssets() << endl;
        cout << ec.payTimes() << endl;
        cout << ec.payAmounts() << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    
    
    
    /** Test EulerPathGenerator
    try {
        size_t ntimesteps = 10, nfactors = 2;
        EulerPathGenerator<NormalRngMinStdRand> pathgen(ntimesteps, nfactors);
        Matrix A(ntimesteps, nfactors);
        pathgen.next(A);
        cout << A << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    
    /** Test rng
     try {
         NormalRngMinStdRand rng(2, 10, 1);
         // NormalRng<minstd_rand> rng(2, 10, 1);
         cout << rng.dim() << endl;
         Vector res(2);
         // rng.next<Vector::iterator>(res.begin(), res.end());
         rng.next(res.begin(), res.end());   // template argument deduction also works here
         cout << res << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
     
    /** Test for meanvarcalculator 
    try {
        MeanVarCalculator<Vector::iterator> mvcal(2);
        Vector v1 = {2,2}, v2 = {2,3}, v3 = {2,10};
        mvcal.addSample(v1.begin(), v1.end());
        mvcal.addSample(v2.begin(), v2.end());
        mvcal.addSample(v3.begin(), v3.end());

        Matrix res = mvcal.results();
        cout << res << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    

    /** Test Polynomial 
    try {
        Vector coeffs = {-1, -1, 6};
        Polynomial poly(coeffs);
        // cout << poly(2) << endl;
        
        double x1 = -1, x2 = 1;
        int n = 20, nroot;
        Vector xb1(10), xb2(10);
        zbrak(poly, x1, x2, n, xb1, xb2, nroot);
        cout << nroot << endl;
        cout << rtsec(poly, x1, x2, 0.01) << endl;
        for (size_t i = 0; i < nroot; ++i) 
            cout << "(" << xb1(i) << ", " << xb2(i) << ")"<< endl;
                                                             
        
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    

    /** Test market
    try {
        auto yc = market().yieldCurves();
        auto vt = market().volatilities();
        // cout << yc.version() << endl;
        vector<double> tMat{1, 2, 3, 4, 5};
        vector<double> spotRate{0.1, 0.1, 0.2, 0.4, 0.1};
        vector<double> fwdVol{0.2, 0.4, 0.3, 0.2, 0.3};
        auto SptrYsp =  make_shared<YieldCurve>(tMat.cbegin(), tMat.cend(),
                                                spotRate.cbegin(), spotRate.cend(),
                                                YieldCurve::InputType::SPOTRATE);
        auto SptrVsp =  make_shared<VolatilityTermStructure>(tMat.cbegin(), tMat.cend(),
                                                             fwdVol.cbegin(), fwdVol.cend(),
                                                             VolatilityTermStructure::VolType::FWDVOL);
        // yc.set("USD", SptrYsp);

        yc.set("USD", SptrYsp);
        yc.set("USD", SptrYsp);
        yc.set("CNY", SptrYsp);


        vt.set("USD", SptrVsp);
        vt.set("USD", SptrVsp);
        vt.set("CNY", SptrVsp);
        
        //  yc.set("USD", RptrYsp);
        // yc.set("CNY", RptrYsp);

        // delete RptrYsp;
        // yc.set("USD", shared_ptr<YieldCurve>(& ysp));
        vector<string> names = vt.list();
        for (string c: names)
            cout << c << endl;
        cout << vt.version("USD") << endl;
        
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    
    /**
    try {
        string s{"  not white  "};
        cout << trim(s) << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    
    /** Test YieldCurve class
    try {
        vector<double> tMat{1, 2, 3, 4, 5};
        vector<double> spotRate{0.1, 0.1, 0.2, 0.4, 0.1};
        vector<double> fwdRate{0.1, 0.1, 0.2, 0.4, 0.1};
        vector<double> zeroBond{0.9, 0.8, 0.7, 0.6, 0.5};
        YieldCurve ysp(tMat.cbegin(), tMat.cend(), spotRate.cbegin(), spotRate.cend(),
                      YieldCurve::InputType::SPOTRATE);
        YieldCurve yfwd(tMat.cbegin(), tMat.cend(), fwdRate.cbegin(), fwdRate.cend(),
                      YieldCurve::InputType::FWDRATE);
        YieldCurve ybond(tMat.cbegin(), tMat.cend(), zeroBond.cbegin(), zeroBond.cend(),
                         YieldCurve::InputType::ZEROBOND);

        cout << ysp.fwdRate(2, 3) << endl;
        cout << yfwd.fwdRate(2, 3) << endl;
        cout << ybond.fwdRate(2, 3) << endl;

        cout << ysp.spotRate(3) << endl;
        cout << yfwd.spotRate(3) << endl;
        cout << ybond.spotRate(3) << endl;

        cout << ysp.discount(3) << endl;
        cout << yfwd.discount(3) << endl;
        cout << ybond.discount(3) << endl;

        cout << ysp.fwdDiscount(2, 3) << endl;
        cout << yfwd.fwdDiscount(2, 3) << endl;
        cout << ybond.fwdDiscount(2, 3) << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */
    
    /** Example to use ASSERT 
    try {
        ASSERT(a > 0, "a must be positive");
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */

    
    /** Test Normal Distribution

    NormalDistribution nb;
    cout << nb.pdf(1) << endl;
    cout << nb.cdf(1) << endl;
    cout << nb.invcdf(0.6) << endl;
    */

    /** Test simplepricers

    try {
        double spot = 10;
        double strike = 10;
        double timeToExp = 1;
        double intRate = 0.1;
        double divYield = 0.05;
        double volatility = 0.1;
        cout << digitalOptionBS(1, spot, strike, timeToExp,
                                 intRate, divYield, volatility) << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }
    */

    /** Test Matrix

    Matrix A(1, 5, arma::fill::randu);
    cout << A << endl;
    cout << A.t() << endl;
    */

    /** Test PiecewisePolynomial 

    try {
        vector<double> x{0.1, 0.2, 0.3, 0.5};
        vector<double> y{1, 2, 3, 10};
        PiecewisePolynomial pp(x.begin(), x.end(), y.begin(), 1);
        cout << "size: " << pp.size() << endl;
        cout << "order: " << pp.order() << endl;

        cout << pp.eval(0.4, 2) << endl;

        x = {1, 2, 3, 5, 10};
        y = {1, 2, 3, 10, -2};
        PiecewisePolynomial p2(x.begin(), x.end(), y.begin(), 1);

        PiecewisePolynomial psum = pp + p2;
        PiecewisePolynomial pprod = pp * p2;
        
        vector<double> xnew{-1, -0.5, 0, 0.15, 0.25, 0.35, 0.4, 0.6, 1};
        vector<double> ynew(9);
        vector<double> integ(9);
        
        double xStart = -1;
        
        pprod.eval(xnew.cbegin(), xnew.cend(), ynew.begin(), 1);
        cout << "Integral = " << pprod.integral(0.35, 0.4) << endl;
        for (auto v: ynew)
            cout << v << endl;

        cout << "Integral = " << endl;
        psum.integral(xStart, xnew.cbegin(), xnew.cend(), integ.begin(), true);
        for (auto v: integ)
            cout << v << endl;
    }
    catch(Exception& e) {    
        cout << e.what() << endl;
    }
    */

    
    
    return 0;
}
