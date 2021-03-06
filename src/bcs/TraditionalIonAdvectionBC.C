#include "TraditionalIonAdvectionBC.h"

// MOOSE includes
#include "MooseVariable.h"

template<>
InputParameters validParams<TraditionalIonAdvectionBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredCoupledVar("potential","The electric potential");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  return params;
}

TraditionalIonAdvectionBC::TraditionalIonAdvectionBC(const InputParameters & parameters) :
  IntegratedBC(parameters),

  _r_units(1. / getParam<Real>("position_units")),
  _r(getParam<Real>("r")),

// Coupled Variables
  _grad_potential(coupledGradient("potential")),
  _potential_id(coupled("potential")),

  _mu(getMaterialProperty<Real>("mu" + _var.name())),
  _e(getMaterialProperty<Real>("e")),
  _sgn(getMaterialProperty<Real>("sgn" + _var.name())),
  _a(0.5)
{}

Real
TraditionalIonAdvectionBC::computeQpResidual()
{
  if ( _normals[_qp] * _sgn[_qp] * -_grad_potential[_qp] > 0.0) {
    _a = 1.0;
  }
  else {
    _a = 0.0;
  }

  return _test[_i][_qp] * _r_units * _a * _sgn[_qp] * _mu[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_u[_qp]) * _normals[_qp];
}

Real
TraditionalIonAdvectionBC::computeQpJacobian()
{
  if ( _normals[_qp] * _sgn[_qp] * -_grad_potential[_qp] > 0.0) {
    _a = 1.0;
  }
  else {
    _a = 0.0;
  }

  return _test[_i][_qp] * _r_units * _a * _sgn[_qp] * _mu[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_u[_qp]) * _phi[_j][_qp] * _normals[_qp];
}

Real
TraditionalIonAdvectionBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_id)
  {
    if ( _normals[_qp] * _sgn[_qp] * -_grad_potential[_qp] > 0.0) {
      _a = 1.0;
    }
    else {
      _a = 0.0;
    }

    return _test[_i][_qp] * _r_units * _a * _sgn[_qp] * _mu[_qp] * -_grad_phi[_j][_qp] * _r_units * std::exp(_u[_qp]) * _normals[_qp];
  }

  else
    return 0.0;
}
