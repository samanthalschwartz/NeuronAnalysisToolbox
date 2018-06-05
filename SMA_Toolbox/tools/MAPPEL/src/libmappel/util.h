
#ifndef _UTIL_H
#define _UTIL_H
#include <cmath>
#include <cassert>
#include <memory>
#include <utility>
#include <armadillo>
#include <stdexcept>

typedef arma::Col<double> VecT;
typedef arma::Col<int> IVecT;
typedef arma::Mat<double> MatT;
typedef arma::Mat<int> IMatT;
// typedef arma::umat UMatT;
typedef arma::Cube<double> CubeT;
typedef arma::field<VecT> VecFieldT;

void enable_all_cpus();

bool istarts_with(const char* s, const char* pattern);
const char * icontains(const char* s, const char* pattern);

inline double restrict_value_range(double val, double minval, double maxval)
{
//     assert(std::isfinite(val));
    return (val<minval) ? minval : ((val>maxval) ? maxval : val);
}

class MappelBadInputException : public std::exception
{
    std::string msg;
    public:
        MappelBadInputException(const std::string &message)
        {
            std::ostringstream stream;
            stream<<"Mappel:BadInput:"<<message;
            msg=stream.str();
        }

        const char* what() const throw()
        {
            return msg.c_str();
        }
};

namespace std {
    template<typename T, typename ...Args>
    std::unique_ptr<T> make_unique( Args&& ...args )
    {
        return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
    }
}/* namespace std */


namespace arma{
    int maxidx(const VecT &v);
}/* namespace arma */

/* Statistics Functions */



#endif /* _UTIL_H */
