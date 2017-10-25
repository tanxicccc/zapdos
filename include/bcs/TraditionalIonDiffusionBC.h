#ifndef TRADITIONALIONDIFFUSIONBC_H
#define TRADITIONALIONDIFFUSIONBC_H

#include "IntegratedBC.h"

class TraditionalIonDiffusionBC;

template<>
InputParameters validParams<TraditionalIonDiffusionBC>();

class TraditionalIonDiffusionBC : public IntegratedBC
{
public:

	TraditionalIonDiffusionBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  Real _r_units;
  Real _r;

  const MaterialProperty<Real> & _kb;
  const MaterialProperty<Real> & _T;
  const MaterialProperty<Real> & _mass;

  Real _v_thermal;
  Real _user_velocity;
};

#endif //TRADITIONALIONDIFFUSIONBC_H
