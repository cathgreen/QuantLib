#ifndef PRODUCT_HPP
#define PRODUCT_HPP

#include "../defines.hpp"
#include "../exception.hpp"
#include "../sptr.hpp"
#include "../math/matrix.hpp"

/** The abstract base class for all financial products.
    It must be inherited by specific product payoffs.
*/
class Product {
public:
    /** Default Ctor */
    explicit Product(const string& payccy = "USD")  
        : payccy_(payccy) {}
    
    virtual ~Product() {}

    /** Returns the fixing (observation) times */
    Vector const& fixTimes() const {
        return fixTimes_;
    }

    /** Returns the payment times */
    Vector const& payTimes() const {
        return payTimes_;
    }

    /** Returns the payment amounts */
    Vector const& payAmounts() const {
        return payAmounts_;
    }

    /** The number of assets this product depends on */
    virtual size_t nAssets() const = 0;
    
     /** Evaluates the product given the passed-in path
         The "pricePath" matrix must have nrows = the number of fixing times
         have ncols = the number of products
         write the results onto payAmounts_
     */
    virtual void eval(const Matrix& pricePath) = 0;

    /** Useful for PDE pricing of products with early exercise features
        Evalutes the product at fixing time index idx, for a vector of current spots,
        and a given continuation value.
    */
    // virtual void eval(size_t idx, Vector const& spot, double contValue) = 0; 

protected:
    string payccy_;
    Vector fixTimes_;        // the fixing (observation) times
    Vector payTimes_;        // the payment times
    
    /** the payment amounts, same size as payTimes_
        note, one product only have one vector of payAmounts
        if a basket of stock, then add all the payAmounts at one time
    */
    Vector payAmounts_;      
    
};

using SPtrProduct = shared_ptr<Product>;
             

#endif // PRODUCT_HPP
