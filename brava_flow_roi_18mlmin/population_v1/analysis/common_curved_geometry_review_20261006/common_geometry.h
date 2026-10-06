// Validation-only prototype. Included only when SVMP_VALIDATION_COMMON_GEOMETRY
// is defined in the isolated build. Production sources/binaries are untouched.
#pragma once
#include "mat_fun.h"
#include <cmath>
#include <stdexcept>

namespace validation_common_geometry {
struct Geometry {
  Vector<double> x = Vector<double>(3);
  Array<double> J = Array<double>(3,3);
  Array<double> inverse = Array<double>(3,3);
  Array<double> metric = Array<double>(3,3);
  double determinant = 0.0;
};

// Only the geometry basis is allowed to construct physical coordinates/J.
// Field derivatives enter transform() separately and cannot overwrite geometry.
inline Geometry evaluate(const Vector<double>& Ngeom,
                         const Array<double>& Dgeom, const Array<double>& X) {
  Geometry geo;
  double h = 0.0;
  for (int a=0; a<X.ncols(); ++a)
    for (int i=0; i<3; ++i) h=std::max(h,std::abs(X(i,a)-X(i,3)));
  if (!(h>0.0)) throw std::runtime_error("Degenerate validation geometry");
  for (int i=0; i<3; ++i) {
    geo.x(i)=X(i,3);
    for (int a=0; a<X.ncols(); ++a) {
      const double dx=X(i,a)-X(i,3);
      geo.x(i)+=Ngeom(a)*dx;
      for (int j=0; j<3; ++j) geo.J(i,j)+=dx*Dgeom(j,a);
    }
  }
  // Scale only to avoid a dimensional small determinant being called singular.
  Array<double> scaled(3,3);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) scaled(i,j)=geo.J(i,j)/h;
  const double det_scaled=mat_fun::mat_det(scaled,3);
  geo.determinant=det_scaled*h*h*h;
  if (!(geo.determinant>0.0) || !std::isfinite(geo.determinant))
    throw std::runtime_error("Nonpositive/nonfinite common curved Jacobian");
  geo.inverse=mat_fun::mat_inv(scaled,3);
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) geo.inverse(i,j)/=h;
  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
    for (int k=0; k<3; ++k) geo.metric(i,j)+=geo.inverse(k,i)*geo.inverse(k,j);
  return geo;
}

inline void transform(const Array<double>& Dreference, const Geometry& geo,
                      Array<double>& Dphysical) {
  Dphysical=0.0;
  for (int a=0; a<Dreference.ncols(); ++a) for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j) Dphysical(i,a)+=Dreference(j,a)*geo.inverse(j,i);
}

// Hessian chain rule with the same inverse, including the geometry-curvature
// term. Packed ordering is the existing native xx,yy,zz,xy,yz,xz convention.
inline void transform_hessian(const Array<double>& Hreference,
                              const Array<double>& X, const Geometry& geo,
                              const Array<double>& Dphysical,
                              Array<double>& Hphysical) {
  constexpr int pair[6][2]={{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};
  Array<double> curvature(3,6);
  for (int i=0; i<3; ++i) for (int k=0; k<6; ++k)
    for (int a=0; a<X.ncols(); ++a)
      curvature(i,k)+=(X(i,a)-X(i,3))*Hreference(k,a);
  Hphysical=0.0;
  for (int a=0; a<X.ncols(); ++a) {
    Array<double> B(3,3);
    for (int k=0; k<6; ++k) {
      double v=Hreference(k,a);
      for (int i=0; i<3; ++i) v-=Dphysical(i,a)*curvature(i,k);
      B(pair[k][0],pair[k][1])=v;B(pair[k][1],pair[k][0])=v;
    }
    for (int k=0; k<6; ++k) for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
      Hphysical(k,a)+=geo.inverse(i,pair[k][0])*B(i,j)*geo.inverse(j,pair[k][1]);
  }
}
} // namespace validation_common_geometry
