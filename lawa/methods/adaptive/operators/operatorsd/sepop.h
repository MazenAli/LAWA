#ifndef LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_SEPOP_H
#define LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_SEPOP_H 1

#include <vector>

namespace lawa
{

template <typename Index, typename Optype>
class Sepop
{
public:
    typedef typename std::vector<Optype*>               Opvec;
    typedef typename std::vector<Optype*>::size_type    size_type;

private:
    Opvec               ops_;
    const size_type     rank;
    const size_type     dim;

public;
    Sepop() = delete;

};

} // namespace lawa

#endif // LAWA_METHODS_ADAPTIVE_OPERATORS_LOCALOPERATORS_SEPOP_H
