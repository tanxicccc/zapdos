#ifndef TRADITIONALIONADVECTIONBC_H
#define TRADITIONALIONADVECTIONBC_H

#include "IntegratedBC.h"

class TraditionalIonAdvectionBC;

template<>
InputParameters validParams<TraditionalIonAdvectionBC>();

class TraditionalIonAdvectionBC : public IntegratedBC
{
public:

	TraditionalIonAdvectionBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real _r_units;
  Real _r;

  // Coupled variables

  const VariableGradient & _grad_potential;
  unsigned int _potential_id;

  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _e;
  const MaterialProperty<Real> & _sgn;

  Real _a;
};

#endif //TRADITIONALIONADVECTIONBC_H
