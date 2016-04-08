#ifndef ELECTRONENERGYRATEIONIZATION_H
#define ELECTRONENERGYRATEIONIZATION_H

#include "Kernel.h"

class ElectronEnergyRateIonization;

template<>
InputParameters validParams<ElectronEnergyRateIonization>();

class ElectronEnergyRateIonization : public Kernel
{
public:
  ElectronEnergyRateIonization(const InputParameters & parameters);
  virtual ~ElectronEnergyRateIonization();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _kiz;
  const MaterialProperty<Real> & _d_kiz_d_actual_mean_en;
  const MaterialProperty<Real> & _Eiz;

  const VariableValue & em;
  unsigned int _em_id;
};


#endif /* ELECTRONENERGYRATEIONIZATION_H */