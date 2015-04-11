#ifndef __NUMBER_H__INCLUDED__
#define __NUMBER_H__INCLUDED__

# include <cstddef>


// Shift bits such that the given bit is free.
template <class Int>
Int skip_bit(Int pool, size_t bit_index)
{
    Int bit = Int(1) << bit_index;
    Int left = (pool & ~(bit-1)) << 1;
    Int right = pool & (bit-1);
    return left | right;
}


// binomial coefficient [n choose r]
template <class Int>
Int nCr(Int n, Int r)
{
    if (r > n)
        return 0;
    if (r * 2 > n)
        r = n - r;
    if (r == 0)
        return 1;
    Int result = n;
    for (Int i = 2; i <= r; ++i) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


// return sign of argument (-1 or +1)
template <class Int>
Int sign(Int a)
{
    return a < 0 ? -1 : 1;
}


// Return greatest common divisor using Euclid's Algorithm.
template <class Int>
Int gcd(Int a, Int b)
{
    while (b) {
        Int c = a % b;
        a = b;
        b = c;
    }
    return a;
}

#endif  // include guard
