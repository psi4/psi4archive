#include "mints.h"

using namespace boost;
using namespace psi;

/////////
// Normal two-electron repulsion integrals
/////////

ERI::ERI(const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new Taylor_Fjt(basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1, 1e-15);
}

ERI::~ERI()
{
    delete fjt_;
}

/////////
// ErfERI 
/////////

ErfERI::ErfERI(double omega, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfFundamental(omega, 
                          basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1);
}

ErfERI::~ErfERI()
{
    delete fjt_;
}

void ErfERI::setOmega(double omega)
{
    (static_cast<ErfFundamental*>(fjt_))->setOmega(omega);
}

/////////
// ErfComplementERI 
/////////

ErfComplementERI::ErfComplementERI(double omega, const IntegralFactory *integral, int deriv, double schwarz)
    : TwoElectronInt(integral, deriv, schwarz)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfComplementFundamental(omega, 
                          basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1);
}

ErfComplementERI::~ErfComplementERI()
{
    delete fjt_;
}

void ErfComplementERI::setOmega(double omega)
{
    (static_cast<ErfComplementFundamental*>(fjt_))->setOmega(omega);
}

