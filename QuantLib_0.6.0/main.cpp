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
*/
#include "pricer/simplepricers.hpp"
#include "market/market.hpp"
#include "math/optim/polyfunc.hpp"
#include "math/optim/roots.hpp"
#include "math/stats/meanvarcalculator.hpp"
#include "methods/montecarlo/eulerpathgenerator.hpp"
#include "products/europeancallput.hpp"
#include "products/asianbasketcallput.hpp"
#include "pricer/bsmcpricer.hpp"

#include <vector>

#include <iostream>

int main() {
    
    try {
        int payoffType = -1;
        double strike = 2, T = 5;
        auto prod = SPtrProduct(new EuropeanCallPut(payoffType, strike, T));
        McParams mcparams;   // default values
        double divYield = 0.2;
        double spot = 8;
        vector<double> tMat{1, 2, 3, 4, 5};
        vector<double> spotRate{0.1, 0.1, 0.1, 0.1, 0.1};
        vector<double> fwdVol{0.2, 0.2, 0.2, 0.2, 0.2};
        auto SptrYsp =  make_shared<YieldCurve>(tMat.cbegin(), tMat.cend(),
                                                spotRate.cbegin(), spotRate.cend(),
                                                YieldCurve::InputType::SPOTRATE);
        auto SptrVsp =  make_shared<VolatilityTermStructure>(tMat.cbegin(), tMat.cend(),
                                                             fwdVol.cbegin(), fwdVol.cend(),
                                                             VolatilityTermStructure::VolType::SPOTVOL);

        MeanVarCalculator<double*> mvcal(1);
        BsMcPricer bsmcpricer(prod, SptrYsp, divYield, SptrVsp, spot, mcparams);
        unsigned long npaths = 10000;
        bsmcpricer.simulate(mvcal, npaths);

        cout << mvcal.results() << endl;

        cout << bsmcpricer.payamts_ << endl;
        cout << bsmcpricer.discfactors_ << endl;
        cout << bsmcpricer.drifts_ << endl;
        cout << bsmcpricer.stdevs_ << endl;
        cout << europeanOptionBS(payoffType, spot, strike, T, spotRate[0], divYield, fwdVol[0]) << endl;
    }
    catch(Exception& e) {
        cout << e.what() << endl;
    }

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
        size_t ntimesteps = 1, nfactors = 1;
        EulerPathGenerator<NormalRngMinStdRand> pathgen(ntimesteps, nfactors);
        Matrix A(ntimesteps, nfactors);
        pathgen.next(A);
        A += 10;
        cout << A << endl;
        
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
        Vector v1 = {1,2}, v2 = {2,3}, v3 = {4,10};
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
