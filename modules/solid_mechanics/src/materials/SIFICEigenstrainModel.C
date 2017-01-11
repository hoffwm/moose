/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "SIFICEigenstrainModel.h"
#include "SymmTensor.h"

template<>
InputParameters validParams<SIFICEigenstrainModel>()
{
  return validParams<Material>();
}

SIFICEigenstrainModel::SIFICEigenstrainModel(const InputParameters & parameters ):
  VolumetricModel( parameters )
{}

SIFICEigenstrainModel::~SIFICEigenstrainModel() {}

void
SIFICEigenstrainModel::modifyStrain(const unsigned int qp,
                                    const Real scale_factor,
                                    SymmTensor & strain_increment,
                                    SymmTensor & /*dstrain_increment_dT*/)
{
  Real ypos = _q_point[qp](1);
  Real xpos = _q_point[qp](0);
  Real t = 0.21971;
  Real R = 2.1971;
  Real E = 1e4;
  Real nu = 0.3;

  Real sigma_zz = 1.0; //uniform
  // Real sigma_zz = ((1/t)*(sqrt(((xpos-R)*(xpos-R))+(ypos*ypos))-R)); //linear
  // Real sigma_zz = (((1/t)*(sqrt(((xpos-R)*(xpos-R))+(ypos*ypos))-R))*((1/t)*(sqrt(((xpos-R)*(xpos-R))+(ypos*ypos))-R))); //quadratic
  // Real sigma_zz = (((1/t)*(sqrt(((xpos-R)*(xpos-R))+(ypos*ypos))-R))*((1/t)*(sqrt(((xpos-R)*(xpos-R))+(ypos*ypos))-R))*((1/t)*(sqrt(((xpos-R)*(xpos-R))+(ypos*ypos))-R))); //cubic


  strain_increment.xx() -=  nu * sigma_zz / E;
  strain_increment.yy() -=  nu * sigma_zz / E;
  strain_increment.zz() +=  sigma_zz / E;

}
