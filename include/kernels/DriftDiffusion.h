#ifndef DRIFTDIFFUSION_H
#define DRIFTDIFFUSION_H

#include "Kernel.h"

class DriftDiffusion;

template <>
InputParameters validParams<DriftDiffusion>();

// This diffusion kernel should only be used with species whose values are in the logarithmic form.

class DriftDiffusion : public Kernel
{
public:
  DriftDiffusion(const InputParameters & parameters);
  virtual ~DriftDiffusion();

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

#endif /* DRIFTDIFFUSION_H */
