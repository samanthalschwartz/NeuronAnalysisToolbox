/**
 * @file Maxima.h
 * @author Mark J. Olah (mjo\@cs.unm.edu)
 * @date 07-28-2014
 * @brief The class declaration for Maxima2D and Maxima3D, local maxima finders.
 */
#ifndef _MAXIMA_H
#define _MAXIMA_H
#include <armadillo>


template <class FloatT>
class Maxima2D
{
    public:
        typedef int IndexT;
        typedef arma::Col<IndexT> IVecT;
        typedef arma::Mat<IndexT> IMatT;
        typedef arma::Col<FloatT> VecT;
        typedef arma::Mat<FloatT> ImageT;
        const int dim=2;
        IVecT size;
        int boxsize;
        Maxima2D(const IVecT &sizeX, int boxsize=3);
        int find_maxima(const ImageT &im);
        int find_maxima(const ImageT &im, IMatT &maxima_out, VecT &max_vals_out);
        void read_maxima(int Nmaxima, IMatT &maxima_out, VecT &max_vals_out) const;
        void test_maxima(const ImageT &im);
private:
        int max_maxima;//size of maxima and max_vals array
        IMatT maxima;// 2xN.
        VecT max_vals; //Nx1
        IMatT skip_buf;

        void detect_maxima(int &Nmaxima, int x, int y, FloatT val);
        int maxima_3x3(const ImageT &im);
        int maxima_3x3_edges(const ImageT &im);
        int maxima_3x3_slow(const ImageT &im);
        int maxima_5x5(const ImageT &im);
        int maxima_nxn(const ImageT &im, int filter_size);
};


template <class FloatT>
class Maxima3D
{
public:
    typedef int IndexT;
    typedef arma::Col<int> IVecT;
    typedef arma::Mat<int> IMatT;
    typedef arma::Cube<int> ICubeT;
    typedef arma::Col<FloatT> VecT;
    typedef arma::Cube<FloatT> ImageT;
    const int dim=3;
    IVecT size;
    int boxsize;
    Maxima3D(const IVecT &size, int boxsize=3);
    int find_maxima(const ImageT &im);
    int find_maxima(const ImageT &im, IMatT &maxima_out, VecT &max_vals_out);
    void read_maxima(int Nmaxima, IMatT &maxima_out, VecT &max_vals_out) const;
    void test_maxima(const ImageT &im);
private:
    int max_maxima;//size of maxima and max_vals array
    IMatT maxima;// 2xN.
    VecT max_vals; //Nx1
    IMatT skip_buf;
    ICubeT skip_plane_buf;

    void detect_maxima(int &Nmaxima, int x, int y, int z, FloatT val);
    int maxima_3x3(const ImageT &im);
    int maxima_3x3_edges(const ImageT &im);
    int maxima_3x3_slow(const ImageT &im);
    int maxima_5x5(const ImageT &im);
    int maxima_nxn(const ImageT &im, int filter_size);
};


#endif /* _MAXIMA_H */
