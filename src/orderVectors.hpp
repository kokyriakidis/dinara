#ifndef DINARA_ORDER_VECTORS_HPP
#define DINARA_ORDER_VECTORS_HPP

#include "vector.hpp"

// Classes to sort vectors by size.

namespace dinara {

    template<class T> class OrderVectorsByIncreasingSize;
    template<class T> class OrderVectorsByDecreasingSize;

}



template<class T> class dinara::OrderVectorsByIncreasingSize {
public:
     bool operator()(const vector<T>& x, const vector<T>& y) const
    {
         return x.size() < y.size();
    }
};



template<class T> class dinara::OrderVectorsByDecreasingSize {
public:
     bool operator()(const vector<T>& x, const vector<T>& y) const
    {
         return x.size() > y.size();
    }
};

#endif

