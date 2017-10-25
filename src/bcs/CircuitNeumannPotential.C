#include "CircuitNeumannPotential.h"

// MOOSE includes
#include "MooseVariable.h"
#include "Function.h"

template<>
InputParameters validParams<CircuitNeumannPotential>()
{
  InputParameters p = validParams<IntegratedBC>();
  p.addRequiredParam<UserObjectName>("data_provider","The name of the UserObject that can provide some data to materials, bcs, etc.");
  p.addRequiredCoupledVar("ip","The ion density.");
  p.addRequiredCoupledVar("em","The log of the electron density.");
  p.addRequiredCoupledVar("mean_en","The log of the product of the mean energy and the electron density.");
  p.addRequiredParam<std::string>("potential_units", "The potential units.");
  p.addRequiredParam<Real>("r", "The reflection coefficient applied to both electrons and ions");
  p.addRequiredParam<Real>("position_units", "Units of position.");
  p.addRequiredParam<FunctionName>("surface_potential", "The electrical potential applied to the surface if no current was flowing in the circuit.");
  p.addRequiredParam<std::string>("surface", "Whether you are specifying the potential on the anode or the cathode with the requirement that the other metal surface be grounded.");
  p.addRequiredParam<Real>("resist", "The ballast resistance in Ohms");
  p.addRequiredParam<Real>("position_units", "Units of position");
  p.addRequiredParam<std::string>("potential_units", "The potential units.");
  p.addRequiredParam<bool>("use_moles", "Whether to convert from units of moles to #.");
  p.addParam<Real>("A", 1., "For 1D calculations, an area has to be passed. This area also must match the units convention of position_units.");
  return p;
}

CircuitNeumannPotential::CircuitNeumannPotential(const InputParameters & parameters) :
    IntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _V_bat(getFunction("surface_potential")),
	
	
    _data(getUserObject<ProvideMobility>("data_provider")),
    _ip_var(*getVar("ip",0)),
    _ip(coupledValue("ip")),
    _grad_ip(coupledGradient("ip")),
    _ip_id(coupled("ip")),
    _mean_en(coupledValue("mean_en")),
    _mean_en_id(coupled("mean_en")),
    _em(coupledValue("em")),
    _em_id(coupled("em")),
	
    _se_coeff(getMaterialProperty<Real>("se_coeff")),
    _muip(getMaterialProperty<Real>("mu"+_ip_var.name())),
    _eps(getMaterialProperty<Real>("eps")),
    _N_A(getMaterialProperty<Real>("N_A")),
    _sgnip(getMaterialProperty<Real>("sgn" + _ip_var.name())),
    _Dip(getMaterialProperty<Real>("diff" + _ip_var.name())),
    _muem(getMaterialProperty<Real>("muem")),
    _d_muem_d_actual_mean_en(getMaterialProperty<Real>("d_muem_d_actual_mean_en")),
    _e(getMaterialProperty<Real>("e")),
    _massem(getMaterialProperty<Real>("massem")),
	_T_heavy(getMaterialProperty<Real>("T" + _ip_var.name())),
    _kb(getMaterialProperty<Real>("k_boltz")),
    _mass(getMaterialProperty<Real>("mass" + _ip_var.name())),
	
	_potential_units(getParam<std::string>("potential_units")),	
	_surface(getParam<std::string>("surface")),
    _resist(getParam<Real>("resist")),
	_convert_moles(getParam<bool>("use_moles")),
	_A(getParam<Real>("A")),


    _ion_flux(0, 0, 0),
    _n_gamma(0),
    _actual_mean_en(0),
    _v_e_th(0),
    _d_v_e_th_d_em(0),
    _d_v_e_th_d_mean_en(0),
    _v_i_th(0),
    _a(0),
    _b(0),
    _d_ion_flux_d_u(0, 0, 0),
    _d_ion_flux_d_ip(0, 0, 0),
    _d_n_gamma_d_u(0),
    _d_n_gamma_d_ip(0),
    _d_n_gamma_d_em(0),
    _d_n_gamma_d_mean_en(0),
    _numerator(0),
    _denominator(0),
    _d_numerator_d_u(0),
    _d_denominator_d_u(0),
    _d_numerator_d_ip(0),
    _d_denominator_d_ip(0),
    _d_numerator_d_em(0),
    _d_denominator_d_em(0),
    _d_numerator_d_mean_en(0),
    _d_denominator_d_mean_en(0)
{
  if (_surface.compare("anode") == 0)
    _current_sign = -1.;
  else if (_surface.compare("cathode") == 0)
    _current_sign = 1.;
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;
}

Real
CircuitNeumannPotential::computeQpResidual()
{
  if ( _normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  if ( _normals[_qp] * _sgnip[_qp] * -_grad_u[_qp] > 0.0) {
    _b = 1.0;
  }
  else {
    _b = 0.0;
  }


  _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_u[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);
  _v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
  _v_i_th = std::sqrt(8 * _kb[_qp] * _T_heavy[_qp] / (M_PI * _mass[_qp]));

  return _test[_i][_qp] * _r_units * _eps[_qp] * ((_u[_qp] - _V_bat.value(_t, _q_point[_qp]))/_e[_qp]/_current_sign * std::pow(_r_units, 2.) / _resist * _voltage_scaling / _A / _N_A[_qp] - (1 - _a) * _se_coeff[_qp]  * _ion_flux* _normals[_qp] + 0.25* _v_e_th * std::exp(_em[_qp]) + 0.25* _v_i_th * std::exp(_ip[_qp])) / (-_a * _muem[_qp] * std::exp(_em[_qp])*_r_units +_b * _sgnip[_qp] * _muip[_qp] * _r_units * std::exp(_ip[_qp]));
}

Real
CircuitNeumannPotential::computeQpJacobian()
{
  if ( _normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  if ( _normals[_qp] * _sgnip[_qp] * -_grad_u[_qp] > 0.0) {
    _b = 1.0;
  }
  else {
    _b = 0.0;
  }

  _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_u[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  _d_ion_flux_d_u = _sgnip[_qp] * _muip[_qp] * -_grad_phi[_j][_qp] * _r_units * std::exp(_ip[_qp]);
  _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);
  _d_n_gamma_d_u = (1. - _a) * _se_coeff[_qp] / _muem[_qp] * (_d_ion_flux_d_u * _normals[_qp] / (-_grad_u[_qp] * _r_units * _normals[_qp]) - _ion_flux * _normals[_qp] / (std::pow(-_grad_u[_qp] * _r_units * _normals[_qp], 2.)) * -_grad_phi[_j][_qp] * _r_units * _normals[_qp]);
  _v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
  _v_i_th = std::sqrt(8 * _kb[_qp] * _T_heavy[_qp] / (M_PI * _mass[_qp]));

  _numerator = (_u[_qp] - _V_bat.value(_t, _q_point[_qp]))/_e[_qp]/_current_sign * std::pow(_r_units, 2.) / _resist * _voltage_scaling / _A / _N_A[_qp] - (1 - _a) * _se_coeff[_qp]  * _ion_flux* _normals[_qp] + 0.25* _v_e_th * std::exp(_em[_qp]) + 0.25* _v_i_th * std::exp(_ip[_qp]);
  _denominator = -_a * _muem[_qp] * std::exp(_em[_qp])*_r_units +_b * _sgnip[_qp] * _muip[_qp] * _r_units * std::exp(_ip[_qp]);

  _d_numerator_d_u = _phi[_j][_qp]/_e[_qp] / _current_sign * std::pow(_r_units, 2.) / _resist * _voltage_scaling / _A / _N_A[_qp] - (1 - _a) * _se_coeff[_qp]  * _d_ion_flux_d_u * _normals[_qp];
  _d_denominator_d_u = 0;

  return _test[_i][_qp] * _r_units * _eps[_qp] * (_d_numerator_d_u * _denominator - _d_denominator_d_u * _numerator) / (_denominator * _denominator);
}

Real
CircuitNeumannPotential::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ip_id)
  {
  if ( _normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  if ( _normals[_qp] * _sgnip[_qp] * -_grad_u[_qp] > 0.0) {
    _b = 1.0;
  }
  else {
    _b = 0.0;
  }

    _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_u[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
    _d_ion_flux_d_ip = _sgnip[_qp] * _muip[_qp] * -_grad_u[_qp] * _r_units * std::exp(_ip[_qp]) * _phi[_j][_qp] - _Dip[_qp] * (std::exp(_ip[_qp]) * _phi[_j][_qp] * _grad_ip[_qp] * _r_units + std::exp(_ip[_qp]) * _grad_phi[_j][_qp] * _r_units);
    _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);
    _d_n_gamma_d_ip  = (1. - _a) * _se_coeff[_qp] * _d_ion_flux_d_ip * _normals[_qp] / (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);
    _v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
    _v_i_th = std::sqrt(8 * _kb[_qp] * _T_heavy[_qp] / (M_PI * _mass[_qp]));

	_numerator = (_u[_qp] - _V_bat.value(_t, _q_point[_qp]))/_e[_qp]/_current_sign * std::pow(_r_units, 2.) / _resist * _voltage_scaling / _A / _N_A[_qp] - (1 - _a) * _se_coeff[_qp]  * _ion_flux* _normals[_qp] + 0.25* _v_e_th * std::exp(_em[_qp]) + 0.25* _v_i_th * std::exp(_ip[_qp]);
	_denominator = -_a * _muem[_qp] * std::exp(_em[_qp])*_r_units +_b * _sgnip[_qp] * _muip[_qp] * _r_units * std::exp(_ip[_qp]);

    _d_numerator_d_ip = (-1. + _a) * _se_coeff[_qp] * _d_ion_flux_d_ip *  _normals[_qp] + 0.25 * std::exp(_ip[_qp]) * _phi[_j][_qp] * _v_i_th;
    _d_denominator_d_ip = _b * _sgnip[_qp] * _muip[_qp] * _r_units * std::exp(_ip[_qp]) * _phi[_j][_qp];

    return _test[_i][_qp] * _r_units * _eps[_qp] * (_d_numerator_d_ip * _denominator - _d_denominator_d_ip * _numerator) / (_denominator * _denominator);
  }

  else if (jvar == _em_id)
  {
  if ( _normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  if ( _normals[_qp] * _sgnip[_qp] * -_grad_u[_qp] > 0.0) {
    _b = 1.0;
  }
  else {
    _b = 0.0;
  }

    _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_u[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
    _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);
    _actual_mean_en = std::exp(_mean_en[_qp] - _em[_qp]);
    _d_n_gamma_d_em = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (std::pow(_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp], 2.)) * -_grad_u[_qp] * _r_units * _normals[_qp] * _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp];
    _v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
    _d_v_e_th_d_em = 0.5 / _v_e_th * 8 * _data.coulomb_charge() * 2.0 / 3 * _actual_mean_en / (M_PI * _massem[_qp]) * -_phi[_j][_qp];
    _v_i_th = std::sqrt(8 * _kb[_qp] * _T_heavy[_qp] / (M_PI * _mass[_qp]));

	_numerator = (_u[_qp] - _V_bat.value(_t, _q_point[_qp]))/_e[_qp]/_current_sign * std::pow(_r_units, 2.) / _resist * _voltage_scaling / _A / _N_A[_qp] - (1 - _a) * _se_coeff[_qp]  * _ion_flux* _normals[_qp] + 0.25* _v_e_th * std::exp(_em[_qp]) + 0.25* _v_i_th * std::exp(_ip[_qp]);
	_denominator = -_a * _muem[_qp] * std::exp(_em[_qp])*_r_units +_b * _sgnip[_qp] * _muip[_qp] * _r_units * std::exp(_ip[_qp]);

    _d_numerator_d_em = 0.25* _v_e_th * std::exp(_em[_qp])* _phi[_j][_qp] + 0.25* _v_i_th * std::exp(_ip[_qp]);
    _d_denominator_d_em = -_a * _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp] * std::exp(_em[_qp])*_r_units -_a * _muem[_qp] * std::exp(_em[_qp]) * _phi[_j][_qp]*_r_units;

    return _test[_i][_qp] * _r_units * _eps[_qp] * (_d_numerator_d_em * _denominator - _d_denominator_d_em * _numerator) / (_denominator * _denominator);
  }

  else if (jvar == _mean_en_id)
  {
  if ( _normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  if ( _normals[_qp] * _sgnip[_qp] * -_grad_u[_qp] > 0.0) {
    _b = 1.0;
  }
  else {
    _b = 0.0;
  }

    _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_u[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
    _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);
    _actual_mean_en = std::exp(_mean_en[_qp] - _em[_qp]);
    _d_n_gamma_d_mean_en = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (std::pow(_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp], 2.)) * -_grad_u[_qp] * _r_units * _normals[_qp] * _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * _phi[_j][_qp];
    _v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
    _d_v_e_th_d_mean_en = 0.5 / _v_e_th * 8 * _data.coulomb_charge() * 2.0 / 3 * _actual_mean_en / (M_PI * _massem[_qp]) * _phi[_j][_qp];
    _v_i_th = std::sqrt(8 * _kb[_qp] * _T_heavy[_qp] / (M_PI * _mass[_qp]));

	_numerator = (_u[_qp] - _V_bat.value(_t, _q_point[_qp]))/_e[_qp]/_current_sign * std::pow(_r_units, 2.) / _resist * _voltage_scaling / _A / _N_A[_qp] - (1 - _a) * _se_coeff[_qp]  * _ion_flux* _normals[_qp] + 0.25* _v_e_th * std::exp(_em[_qp]) + 0.25* _v_i_th * std::exp(_ip[_qp]);
	_denominator = -_a * _muem[_qp] * std::exp(_em[_qp])*_r_units +_b * _sgnip[_qp] * _muip[_qp] * _r_units * std::exp(_ip[_qp]);

    _d_numerator_d_mean_en = 0.25 * _d_v_e_th_d_mean_en * std::exp(_em[_qp]);
    _d_denominator_d_mean_en = -_a * _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * _phi[_j][_qp] * std::exp(_em[_qp]);

    return _test[_i][_qp] * _r_units * _eps[_qp] * (_d_numerator_d_mean_en * _denominator - _d_denominator_d_mean_en * _numerator) / (_denominator * _denominator);
  }

  else
    return 0.0;
}
