#include "TraditionalEnergyBC.h"

// MOOSE includes
#include "MooseVariable.h"

template<>
InputParameters validParams<TraditionalEnergyBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<Real>("r", "The reflection coefficient");
    params.addRequiredCoupledVar("potential","The electric potential");
    params.addRequiredCoupledVar("em", "The electron density.");
    params.addRequiredCoupledVar("ip", "The ion density.");
    params.addRequiredParam<Real>("position_units", "Units of position.");
    return params;
}

TraditionalEnergyBC::TraditionalEnergyBC(const InputParameters & parameters) :
  IntegratedBC(parameters),

  _r_units(1. / getParam<Real>("position_units")),
  _r(getParam<Real>("r")),

// Coupled Variables
  _grad_potential(coupledGradient("potential")),
  _potential_id(coupled("potential")),
  _em(coupledValue("em")),
  _em_id(coupled("em")),
  _ip_var(*getVar("ip",0)),
  _ip(coupledValue("ip")),
  _grad_ip(coupledGradient("ip")),
  _ip_id(coupled("ip")),

  _muem(getMaterialProperty<Real>("muem")),
  _d_muem_d_actual_mean_en(getMaterialProperty<Real>("d_muem_d_actual_mean_en")),
  _massem(getMaterialProperty<Real>("massem")),
  _e(getMaterialProperty<Real>("e")),
  _sgnip(getMaterialProperty<Real>("sgn" + _ip_var.name())),
  _muip(getMaterialProperty<Real>("mu" + _ip_var.name())),
  _Dip(getMaterialProperty<Real>("diff" + _ip_var.name())),
  _se_coeff(getMaterialProperty<Real>("se_coeff")),
  _se_energy(getMaterialProperty<Real>("se_energy")),
  _mumean_en(getMaterialProperty<Real>("mumean_en")),
  _d_mumean_en_d_actual_mean_en(getMaterialProperty<Real>("d_mumean_en_d_actual_mean_en")),

  _a(0.5),
  _v_thermal(0),
  _ion_flux(0, 0, 0),
  _n_gamma(0),
  _d_v_thermal_d_u(0),
  _d_v_thermal_d_em(0),
  _d_ion_flux_d_potential(0, 0, 0),
  _d_ion_flux_d_ip(0, 0, 0),
  _d_n_gamma_d_potential(0),
  _d_n_gamma_d_ip(0),
  _d_n_gamma_d_u(0),
  _d_n_gamma_d_em(0),
  _actual_mean_en(0)
  {}

Real
TraditionalEnergyBC::computeQpResidual()
{

  
  if ( _normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0) {
    _a = 1.0;
  }
  else {
    _a = 0.0;
  }
  
  _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
  _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));

   return 5. / 3. * _test[_i][_qp] * _r_units * (0.25 * _v_thermal * std::exp(_u[_qp]));
//   return 5. / 3. * _test[_i][_qp] * _r_units * (-_a * _mumean_en[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_u[_qp]) * _normals[_qp] + 0.25 * _v_thermal * std::exp(_u[_qp]));
}

Real
TraditionalEnergyBC::computeQpJacobian()
{
  if ( _normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0) {
    _a = 1.0;
  }
  else {
    _a = 0.0;
  }

  _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
  _actual_mean_en = std::exp(_u[_qp] - _em[_qp]);
  // _d_n_gamma_d_u = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) + std::numeric_limits<double>::epsilon()) * -_grad_potential[_qp] * _r_units * _normals[_qp] * _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * _phi[_j][_qp];
  _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
  _d_v_thermal_d_u = 0.5 * _v_thermal * _phi[_j][_qp];
  
      

  return  5. / 3. * _test[_i][_qp] * _r_units * (0.25 * _d_v_thermal_d_u * std::exp(_u[_qp]) + 0.25 * _v_thermal * std::exp(_u[_qp]) * _phi[_j][_qp]);
//  return  5. / 3. * _test[_i][_qp] * _r_units * (-_a * _mumean_en[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_u[_qp]) * _phi[_j][_qp] * _normals[_qp] - _a * _d_mumean_en_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp] * -_grad_potential[_qp] * _r_units * std::exp(_u[_qp]) * _normals[_qp] + 0.25 * _d_v_thermal_d_u * std::exp(_u[_qp]) + 0.25 * _v_thermal * std::exp(_u[_qp]) * _phi[_j][_qp]);
}

Real
TraditionalEnergyBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  
  if (jvar == _potential_id)
  {
    if ( _normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
      _a = 1.0;
    else
      _a = 0.0;

    _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
    _d_ion_flux_d_potential = _sgnip[_qp] * _muip[_qp] * -_grad_phi[_j][_qp] * _r_units * std::exp(_ip[_qp]);
    // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
    // _d_n_gamma_d_potential = (1. - _a) * _se_coeff[_qp] / _muem[_qp] * (_d_ion_flux_d_potential * _normals[_qp] / (-_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon()) - _ion_flux * _normals[_qp] / (std::pow(-_grad_potential[_qp] * _r_units * _normals[_qp], 2.) + std::numeric_limits<double>::epsilon()) * -_grad_phi[_j][_qp] * _r_units * _normals[_qp]);
    _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));

	
    return 0.0;
	//    return 5. / 3. * _test[_i][_qp] * _r_units * (-_a * _mumean_en[_qp] * -_grad_phi[_j][_qp] * _r_units * std::exp(_u[_qp]) * _normals[_qp]);
  }

  else if (jvar == _em_id)
  {
    if ( _normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0) {
      _a = 1.0;
    }
    else {
      _a = 0.0;
    }
    _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
    _d_v_thermal_d_em = 0.5 *_v_thermal * -_phi[_j][_qp];
    // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
    _actual_mean_en = std::exp(_u[_qp] - _em[_qp]);
    // _d_n_gamma_d_em = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) + std::numeric_limits<double>::epsilon()) * -_grad_potential[_qp] * _r_units * _normals[_qp] * _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp];

    return  5. / 3. * _test[_i][_qp] * _r_units * (0.25 * _d_v_thermal_d_em * std::exp(_u[_qp]));
//    return  5. / 3. * _test[_i][_qp] * _r_units * (-_a * _d_mumean_en_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp] * -_grad_potential[_qp] * _r_units * std::exp(_u[_qp]) * _normals[_qp] + 0.25 * _d_v_thermal_d_em * std::exp(_u[_qp]));
  }

  else if (jvar == _ip_id)
  {
    if ( _normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0) {
      _a = 1.0;
    }
    else {
      _a = 0.0;
    }
    _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
    _d_ion_flux_d_ip = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) * _phi[_j][_qp] - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_phi[_j][_qp] * _r_units - _Dip[_qp] * std::exp(_ip[_qp]) * _phi[_j][_qp] * _grad_ip[_qp] * _r_units;
    // _d_n_gamma_d_ip = (1. - _a) * _se_coeff[_qp] * _d_ion_flux_d_ip * _normals[_qp] / (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
    return  0.0;
  }

  else
    return 0.0;
}
