/** @file filter_kernels.h
 * @author Mark J. Olah (mjo@@cs.unm.edu)
 * @date 07-25-2014
 * @brief Low level filters templated on their floating point type
 *
 */
#ifndef _FILTER_KERNELS_H
#define _FILTER_KERNELS_H

#include <armadillo>

/* 1D Gauss FIR Filters */
template <class FloatT>
void gaussFIR_1D(const arma::Col<FloatT> &data, arma::Col<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_1D(int size, const FloatT data[], FloatT fdata[], int hw, const FloatT kernel[]);

template <class FloatT>
void gaussFIR_1D_small(int size, const FloatT data[], FloatT fdata[], int hw, const FloatT kernel[]);

template <class FloatT>
void gaussFIR_1D_arma(const arma::Col<FloatT> &data, arma::Col<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_1D_inplace_arma(arma::Col<FloatT> &data, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_1D_inplace(int size, FloatT data[], int hw, const FloatT kernel[]);

/* 2D Gauss FIR Filters */
template <class FloatT>
void gaussFIR_2Dx(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_2Dx_small(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_2Dx_arma(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_2Dy(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_2Dy_rowmajor(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_2Dy_colmajor(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_2Dy_small(const arma::Mat<FloatT> &data, arma::Mat<FloatT> &fdata, const arma::Col<FloatT> &kernel);

/* 3D Gauss FIR Filters */
template <class FloatT>
void gaussFIR_3Dx(const arma::Cube<FloatT> &data, arma::Cube<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_3Dx_small(const arma::Cube<FloatT> &data, arma::Cube<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_3Dy(const arma::Cube<FloatT> &data, arma::Cube<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_3Dy_small(const arma::Cube<FloatT> &data, arma::Cube<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_3Dz(const arma::Cube<FloatT> &data, arma::Cube<FloatT> &fdata, const arma::Col<FloatT> &kernel);

template <class FloatT>
void gaussFIR_3Dz_small(const arma::Cube<FloatT> &data, arma::Cube<FloatT> &fdata, const arma::Col<FloatT> &kernel);



#endif /* _FILTER_KERNELS_H */
