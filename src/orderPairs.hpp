#ifndef DINARA_ORDER_PAIRS_HPP
#define DINARA_ORDER_PAIRS_HPP


// Classes to sort pairs using various criteria.

namespace dinara {

    template<class First, class Second> class OrderPairsByFirstOnly;
    template<class First, class Second> class OrderPairsByFirstOnlyGreater;

    template<class First, class Second> class OrderPairsBySecondOnly;
    template<class First, class Second> class OrderPairsBySecondOnlyGreater;

    template<class First, class Second> class OrderPairsBySecondThenByFirst;
    template<class First, class Second> class OrderPairsBySecondGreaterThenByFirstLess;

    template<class First, class Second> class OrderPairsBySizeOfSecondGreater;

}



template<class First, class Second> class dinara::OrderPairsByFirstOnly {
public:
    using Pair = pair<First, Second>;
    bool operator()(const Pair& x, const Pair& y) const
    {
         return x.first < y.first;
    }
};



template<class First, class Second> class dinara::OrderPairsByFirstOnlyGreater {
public:
    using Pair = pair<First, Second>;
    bool operator()(const Pair& x, const Pair& y) const
    {
         return x.first > y.first;
    }
};



template<class First, class Second> class dinara::OrderPairsBySecondOnly {
public:
    using Pair = pair<First, Second>;
    bool operator()(const Pair& x, const Pair& y) const
    {
         return x.second < y.second;
    }
};



template<class First, class Second> class dinara::OrderPairsBySecondOnlyGreater {
public:
    using Pair = pair<First, Second>;
    bool operator()(const Pair& x, const Pair& y) const
    {
         return x.second > y.second;
    }
};



template<class First, class Second> class dinara::OrderPairsBySecondThenByFirst {
public:
    bool operator()(const pair<First, Second>& x, const pair<First, Second>& y) const
    {
        if(x.second < y.second) return true;
        if(y.second < x.second) return false;
        return x.first < y.first;
    }
};



template<class First, class Second> class dinara::OrderPairsBySecondGreaterThenByFirstLess {
public:
    using Pair = pair<First, Second>;
    bool operator()(const Pair& x, const Pair& y) const
    {
        if(x.second > y.second) return true;
        if(y.second > x.second) return false;
        return x.first < y.first;
    }
};



template<class First, class Second> class dinara::OrderPairsBySizeOfSecondGreater {
public:
    using Pair = pair<First, Second>;
    bool operator()(const Pair& x, const Pair& y) const
    {
        return x.second.size() > y.second.size();
    }
};



#endif

