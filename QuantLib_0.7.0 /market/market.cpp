#include "market.hpp"

Market& Market::instance() {
    static Market theMarket;
    return theMarket;
}

Market& market() {
    return Market::instance();
}

void Market::clear() {
    ycmap_.clear();
    volmap_.clear();
}
