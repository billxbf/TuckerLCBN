/*
 * Tucker_Vector.hpp
 *
 *  Created on: Mar 15, 2017
 *      Author: amklinv
 */

#ifndef SERIAL_TUCKER_VECTOR_HPP_
#define SERIAL_TUCKER_VECTOR_HPP_

#include "Tucker_Matrix.hpp"

namespace Tucker {

class Vector : public Matrix {
public:
  /** \brief Constructor
   * \param[in] nrows Number of rows
   * \param[in] ncols Number of columns
   */
  Vector(const int nrows) :
    Matrix(nrows)
  {

  }

  /** \brief Element access
   *
   * \param[in] i Index of entry to be returned
   *
   * \exception std::out_of_range \a i is not in the range [0, nsz_)
   */
  double& operator[](const int i) {
    if(i < 0 || i >= nrows())
      throw std::out_of_range("invalid index");
    return data_[i];
  }

  /** \brief Const element access
   *
   * \param[in] i Index of entry to be returned
   *
   * \exception std::out_of_range \a i is not in the range [0, nsz_)
   */
  const double& operator[](const int i) const {
    if(i < 0 || i >= nrows())
      throw std::out_of_range("invalid index");
    return data_[i];
  }

  /// Returns the number of rows
  int nrows() const
  {
    return I_[0];
  }

  /// Returns the number of columns
  int ncols() const
  {
    return 1;
  }

private:
  /// @cond EXCLUDE
  // Disables the copy constructor
  Vector(const Vector& m);
  /// @endcond
};

} // end of namespace Tucker



#endif /* SERIAL_TUCKER_VECTOR_HPP_ */
