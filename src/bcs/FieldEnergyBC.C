#include "FieldEnergyBC.h"

// MOOSE includes
#include "MooseVariable.h"

template<>
InputParameters validParams<FieldEnergyBC>()
{
    InputParameters params = validParams<IntegratedBC>();
    params.addRequiredParam<Real>("r", "The reflection coefficient");
    params.addRequiredCoupledVar("potential","The electric potential");
    params.addRequiredCoupledVar("em", "The electron density.");
    params.addRequiredCoupledVar("ip", "The ion density.");
    params.addRequiredParam<Real>("position_units", "Units of position.");
	params.addRequiredParam<std::string>("potential_units", "The potential units.");
	params.addParam<Real>("tau", 1e-9, "The time constant for ramping the boundary condition."); //FE
    params.addParam<bool>("relax", false, "Use relaxation for emission."); //FE
    return params;
}

FieldEnergyBC::FieldEnergyBC(const InputParameters & parameters) :
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
  _work_function(getMaterialProperty<Real>("work_function")),    //FE
  _field_enhancement(getMaterialProperty<Real>("field_enhancement")),  //FE
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
  _actual_mean_en(0),
  _tau(getParam<Real>("tau")),   //FE
  _relax(getParam<bool>("relax")),  //FE
  _potential_units(getParam<std::string>("potential_units"))  //FE
{ //FE
  if (_potential_units.compare("V") == 0) {
    _voltage_scaling = 1.;
  } else if (_potential_units.compare("kV") == 0) {
    _voltage_scaling = 1000;
  }

  FE_a = 1.541434E-6;        // A eV/V^2
  FE_b = 6.830890E9;                                        // V/m-eV^1.5
  FE_c = 1.439964E-9;                                       // eV^2*m/V

}

Real
FieldEnergyBC::computeQpResidual()
{
  Real v;
  Real f;
  Real jFE;
  Real F;
  Real _relaxation_Expr;
  
  if ( _normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0) {
    _a = 1.0;
  }
  else {
    _a = 0.0;
  }
  
      F = -(1 - _a) * _field_enhancement[_qp] * _normals[_qp] * _grad_potential[_qp] * _r_units * _voltage_scaling;

    f = FE_c * F / std::pow(_work_function[_qp], 2) ;
    v = 1 - f + (f/6)*std::log(f);

    jFE = (FE_a / (_work_function[_qp])) * std::pow( F , 2) * std::exp(-v * FE_b * std::pow(_work_function[_qp], 1.5) / F);

    if ( _relax == true )
      _relaxation_Expr = std::tanh(_t / _tau) ;
    else
      _relaxation_Expr = 1.0 ;
  
  _ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) - _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
  _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));

   return - 5. / 3. * 2. / (1. + _r) * (1 - _a) * _relaxation_Expr * _test[_i][_qp] * _r_units * jFE/ (_e[_qp] * 6.02E23) * _se_energy[_qp];
}

Real
FieldEnergyBC::computeQpJacobian()
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
  
      

  return  0.0;
  }

Real
FieldEnergyBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real v;
  Real f;
  Real F;
  Real jFE;
  Real _d_jFE_d_potential;
  Real _relaxation_Expr;
  
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

	F = -(1 - _a) * _field_enhancement[_qp] * _normals[_qp] * _grad_potential[_qp] * _r_units * _voltage_scaling;

      f = FE_c * F / std::pow(_work_function[_qp], 2) ;
      v = 1 - f + (f/6)*std::log(f);

      jFE = (FE_a / (_work_function[_qp])) * std::pow( F , 2) * std::exp(-v * FE_b * std::pow(_work_function[_qp], 1.5) / F);
      _d_jFE_d_potential = jFE * ( 2 - ( FE_b * FE_c ) / ( 6 * sqrt( _work_function[_qp] ) )
                                   + ( FE_b * std::pow(_work_function[_qp], 1.5) / F ) ) * (_grad_phi[_j][_qp] * _normals[_qp] ) / (_grad_potential[_qp] * _normals[_qp] ) ;

      if ( _relax == true )
        _relaxation_Expr = std::tanh(_t / _tau) ;
      else
        _relaxation_Expr = 1.0 ;
	
    return -_r * 5. / 3. * _relaxation_Expr * _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1 - _a) * _d_jFE_d_potential / ( _e[_qp] * 6.02E23 ) * _se_energy[_qp];
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

    return  0.0;
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
