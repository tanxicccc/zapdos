#ifndef DRIFTDIFFUSIONDONOTHINGBC_H
#define DRIFTDIFFUSIONDONOTHINGBC_H

#include "IntegratedBC.h"

class DriftDiffusionDoNothingBC;

template <>
InputParameters validParams<DriftDiffusionDoNothingBC>();

// This diffusion kernel should only be used with species whose values are in the logarithmic form.

class DriftDiffusionDoNothingBC : public IntegratedBC
{
public:
  DriftDiffusionDoNothingBC(const InputParameters & parameters);
  virtual ~DriftDiffusionDoNothingBC();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real _r_units;

  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _sign;
  MaterialProperty<Real> _user_mu;
  MaterialProperty<Real> _user_sign;

  const MaterialProperty<Real> & _diffusivity;
  MaterialProperty<Real> _user_diff;

  unsigned int _potential_id;
  const VariableGradient & _grad_potential;
  VariableGradient _minus_e_field;

  Real _d_diffusivity_d_u;
};

#endif /* DRIFTDIFFUSIONDONOTHINGBC_H */
