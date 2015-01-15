/* gbmodels.cpp -- contains the functionality defining the GB models used in
 * Amber for use with OpenMM
 */

#include <string>
#include <sstream>

#include "amber/gbmodels.h"
#include "amber/exceptions.h"

using namespace std;
using namespace Amber;

namespace Amber {
// Some useful 
void _createEnergyTerms(OpenMM::CustomGBForce *force, double solventDielectric,
                        double soluteDielectric, double offset, double cutoff,
                        double kappa, bool useSASA) {
    stringstream params;
    params << "; solventDielectric=" << solventDielectric << "; soluteDielectric="
           << soluteDielectric << "; offset=" << offset;
    if (cutoff > 0)
        params << "; cutoff=" << cutoff;
    // The main energy term (depends on salt concentration)
    if (kappa > 0) {
        params << "; kappa=" << kappa;
        string energy = "-0.5*138.935485*(1/soluteDielectric-exp(-kappa*B)/"
                        "solventDielectric)*q^2/B" + params.str();
        force->addEnergyTerm(energy, OpenMM::CustomGBForce::SingleParticle);
    } else {
        string energy = "-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*q^2/B"
                      + params.str();
        force->addEnergyTerm(energy, OpenMM::CustomGBForce::SingleParticle);
    }
    // SASA term, if applicable
    if (useSASA) {
        stringstream iss;
        iss << "28.3919551*(radius+0.14)^2*(radius/B)^6; radius=or+offset"
            << params.str();
        force->addEnergyTerm(iss.str(), OpenMM::CustomGBForce::SingleParticle);
    }
    //
    if (cutoff <= 0) {
        stringstream iss;
        iss << "-138.935485*(1/soluteDielectric-";
        if (kappa > 0)
            iss << "exp(-kappa*f)";
        else
            iss << "1";
        iss << "/solventDielectric)*q1*q2/f; f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))"
            << params.str();
        force->addEnergyTerm(iss.str(),
                             OpenMM::CustomGBForce::ParticlePairNoExclusions);
    } else {
        stringstream iss;
        iss << "-138.935485*(1/soluteDielectric-";
        if (kappa > 0)
            iss << "exp(-kappa*f)";
        else
            iss << "1";
        iss << "/solventDielectric)*q1*q2*(1/f-" << cutoff / 10
            << "; f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))" << params.str();
        force->addEnergyTerm(iss.str(),
                             OpenMM::CustomGBForce::ParticlePairNoExclusions);
    }
}

OpenMM::CustomGBForce *GB_HCT(AmberParm const& amberParm,
                              double solventDielectric,
                              double soluteDielectric,
                              bool useSASA,
                              double cutoff,
                              double kappa) {
    // Convert kappa and the cutoff from 1/A and A to 1/nm and nm, respectively
    cutoff /= 10.0;
    kappa *= 10.0;
    // Create the force
    OpenMM::CustomGBForce *force = new OpenMM::CustomGBForce();
    force->addPerParticleParameter("q");  // charge
    force->addPerParticleParameter("or"); // offset radius
    force->addPerParticleParameter("sr"); // scaled offset radius
    stringstream iss;
    iss << "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
        << "U=r+sr2;"
        << "L=max(or, D);"
        << "D=abs(r-sr2)";
    force->addComputedValue(string("I"), iss.str(),
                            OpenMM::CustomGBForce::ParticlePairNoExclusions);
    // clear the stringstream
    iss.str(string());
    iss << "1/(1/or-tanh(0.8*psi+2.909125*psi^3)/radius);"
        << "psi=I*or; radius=or+offset; offset=0.009";
    force->addComputedValue(string("B"), iss.str(),
                            OpenMM::CustomGBForce::SingleParticle);
    _createEnergyTerms(force, solventDielectric, soluteDielectric, 0.009,
                       cutoff, kappa, useSASA);
    // Now we've built our force -- populate it with the particles
    for (AmberParm::atom_iterator it = amberParm.AtomBegin();
            it != amberParm.AtomEnd(); it++) {
        vector<double> params;
        params.push_back(it->getCharge());
        params.push_back(it->getGBRadius() * 0.1);
        params.push_back(it->getGBScreen());
        force->addParticle(params);
    }
    return force;
}

OpenMM::CustomGBForce *GB_OBC1(AmberParm const& amberParm,
                               double solventDielectric,
                               double soluteDielectric,
                               bool useSASA,
                               double cutoff,
                               double kappa) {
    // Convert kappa and the cutoff from 1/A and A to 1/nm and nm, respectively
    cutoff /= 10.0;
    kappa *= 10.0;
    // Create the force
    OpenMM::CustomGBForce *force = new OpenMM::CustomGBForce();
    force->addPerParticleParameter("q");  // charge
    force->addPerParticleParameter("or"); // offset radius
    force->addPerParticleParameter("sr"); // scaled offset radius
    stringstream iss;
    iss << "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
        << "U=r+sr2;"
        << "L=max(or1, D);"
        << "D=abs(r-sr2)";
    force->addComputedValue(string("I"), iss.str(),
                            OpenMM::CustomGBForce::ParticlePairNoExclusions);
    // Clear the stringstream
    iss.str(string());
    iss << "1/(1/or-tanh(0.8*psi+2.909125*psi^3)/radius);"
        << "psi=I*or; radius=or+offset; offset=0.009";
    force->addComputedValue(string("B"), iss.str(),
                            OpenMM::CustomGBForce::SingleParticle);
    _createEnergyTerms(force, solventDielectric, soluteDielectric, 0.009,
                       cutoff, kappa, useSASA);
    // Now we've built our force -- populate it with the particles
    for (AmberParm::atom_iterator it = amberParm.AtomBegin();
            it != amberParm.AtomEnd(); it++) {
        vector<double> params;
        params.push_back(it->getCharge());
        params.push_back(it->getGBRadius() * 0.1);
        params.push_back(it->getGBScreen());
        force->addParticle(params);
    }
    return force;
}

OpenMM::CustomGBForce *GB_OBC2(AmberParm const& amberParm,
                               double solventDielectric,
                               double soluteDielectric,
                               bool useSASA,
                               double cutoff,
                               double kappa) {
    // Convert kappa and the cutoff from 1/A and A to 1/nm and nm, respectively
    cutoff /= 10.0;
    kappa *= 10.0;
    // Create the force
    OpenMM::CustomGBForce *force = new OpenMM::CustomGBForce();
    force->addPerParticleParameter("q");  // charge
    force->addPerParticleParameter("or"); // offset radius
    force->addPerParticleParameter("sr"); // scaled offset radius
    stringstream iss;
    iss << "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
        << "U=r+sr2;"
        << "L=max(or1, D);"
        << "D=abs(r-sr2)";
    force->addComputedValue(string("I"), iss.str(),
                            OpenMM::CustomGBForce::ParticlePairNoExclusions);
    // Clear the stringstream
    iss.str(string());
    iss << "1/(1/or-tanh(psi-0.8*psi^2+4.85*psi^3)/radius);"
        << "psi=I*or; radius=or+offset; offset=0.009";
    force->addComputedValue(string("B"), iss.str(),
                            OpenMM::CustomGBForce::SingleParticle);
    _createEnergyTerms(force, solventDielectric, soluteDielectric, 0.009,
                       cutoff, kappa, useSASA);
    // Now we've built our force -- populate it with the particles
    for (AmberParm::atom_iterator it = amberParm.AtomBegin();
            it != amberParm.AtomEnd(); it++) {
        vector<double> params;
        params.push_back(it->getCharge());
        params.push_back(it->getGBRadius() * 0.1);
        params.push_back(it->getGBScreen());
        force->addParticle(params);
    }
    return force;
}

OpenMM::CustomGBForce *GB_GBn(AmberParm const& amberParm,
                              double solventDielectric,
                              double soluteDielectric,
                              bool useSASA,
                              double cutoff,
                              double kappa) {
    // Convert kappa and the cutoff from 1/A and A to 1/nm and nm, respectively
    cutoff /= 10.0;
    kappa *= 10.0;
    // Create the force
    OpenMM::CustomGBForce *force = new OpenMM::CustomGBForce();
    force->addPerParticleParameter("q");
    force->addPerParticleParameter("or");
    force->addPerParticleParameter("sr");
    // Build the tabulated functions
    vector<double> d0, m0;
    vector<double>::size_type s = 21*21;
    d0.reserve(s);
    m0.reserve(s);
    for (int i = 0; i < (int)s; i++) {
        d0.push_back(D0[i]/10.0);
        m0.push_back(M0[i]*10.0);
    }
    force->addTabulatedFunction("getd0", new OpenMM::Discrete1DFunction(d0));
    force->addTabulatedFunction("getm0", new OpenMM::Discrete1DFunction(m0));
    stringstream iss;
    iss << "Ivdw+neckScale*Ineck;"
        << "Ineck=step(radius1+radius2+neckCut-r)*getm0(index)/"
        << "(1+100*(r-getd0(index))^2+0.3*1000000*(r-getd0(index))^6);"
        << "index = (radius2*200-20)*21 + (radius1*200-20);"
        << "Ivdw=step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
        << "U=r+sr2;"
        << "L=max(or1, D);"
        << "D=abs(r-sr2);"
        << "radius1=or1+offset; radius2=or2+offset;"
        << "neckScale=0.361825; neckCut=0.68; offset=0.009";
    force->addComputedValue(string("I"), iss.str(),
                            OpenMM::CustomGBForce::ParticlePairNoExclusions);
    // Clear the stringstream
    iss.str(string());
    iss << "1/(1/or-tanh(1.09511284*psi-1.907992938*psi^2+2.50798245*psi^3)/radius);"
        << "psi=I*or; radius=or+offset; offset=0.009";
    force->addComputedValue(string("B"), iss.str(),
                            OpenMM::CustomGBForce::SingleParticle);
    _createEnergyTerms(force, solventDielectric, soluteDielectric, 0.009,
                       cutoff, kappa, useSASA);
    // Now we've built our force -- populate it with the particles
    for (AmberParm::atom_iterator it = amberParm.AtomBegin();
            it != amberParm.AtomEnd(); it++) {
        vector<double> params;
        params.push_back(it->getCharge());
        params.push_back(it->getGBRadius() * 0.1);
        // Screening parameters have been replaced
        switch(it->getElement()) {
            case 1:
                params.push_back(1.09085413633);
                break;
            case 6:
                params.push_back(0.48435382330);
                break;
            case 7:
                params.push_back(0.700147318409);
                break;
            case 8:
                params.push_back(1.06557401132);
                break;
            case 16:
                params.push_back(0.602256336067);
                break;
            default:
                params.push_back(0.5); // not optimized
                break;
        }
        force->addParticle(params);
    }
    return force;
}

OpenMM::CustomGBForce *GB_GBn2(AmberParm const& amberParm,
                               double solventDielectric,
                               double soluteDielectric,
                               bool useSASA,
                               double cutoff,
                               double kappa) {
    // Convert kappa and the cutoff from 1/A and A to 1/nm and nm, respectively
    cutoff /= 10.0;
    kappa *= 10.0;
    // Create the force
    OpenMM::CustomGBForce *force = new OpenMM::CustomGBForce();
    force->addPerParticleParameter("q");
    force->addPerParticleParameter("or");
    force->addPerParticleParameter("sr");
    force->addPerParticleParameter("alpha");
    force->addPerParticleParameter("beta");
    force->addPerParticleParameter("gamma");
    // Build the tabulated functions
    vector<double> d0, m0;
    vector<double>::size_type s = 21*21;
    d0.reserve(s);
    m0.reserve(s);
    for (int i = 0; i < (int)s; i++) {
        d0.push_back(D0[i]/10.0);
        m0.push_back(M0[i]*10.0);
    }
    force->addTabulatedFunction("getd0", new OpenMM::Discrete1DFunction(d0));
    force->addTabulatedFunction("getm0", new OpenMM::Discrete1DFunction(m0));
    stringstream iss;
    iss << "Ivdw+neckScale*Ineck;"
        << "Ineck=step(radius1+radius2+neckCut-r)*getm0(index)/"
        << "(1+100*(r-getd0(index))^2+0.3*1000000*(r-getd0(index))^6);"
        << "index = (radius2*200-20)*21 + (radius1*200-20);"
        << "Ivdw=step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(r-sr2^2/r)*(1/(U^2)-1/(L^2))+0.5*log(L/U)/r);"
        << "U=r+sr2;"
        << "L=max(or1, D);"
        << "D=abs(r-sr2);"
        << "radius1=or1+offset; radius2=or2+offset;"
        << "neckScale=0.826836; neckCut=0.68; offset=0.0195141";
    force->addComputedValue(string("I"), iss.str(),
                            OpenMM::CustomGBForce::ParticlePairNoExclusions);
    // Clear the stringstream
    iss.str(string());
    iss << "1/(1/or-tanh(alpha*psi-beta*psi^2+gamma*psi^3)/radius);"
        << "psi=I*or; radius=or+offset; offset=0.0195141";
    force->addComputedValue(string("B"), iss.str(),
                            OpenMM::CustomGBForce::SingleParticle);
    _createEnergyTerms(force, solventDielectric, soluteDielectric, 0.009,
                       cutoff, kappa, useSASA);
    // Now we've built our force -- populate it with the particles
    for (AmberParm::atom_iterator it = amberParm.AtomBegin();
            it != amberParm.AtomEnd(); it++) {
        vector<double> params;
        params.push_back(it->getCharge());
        params.push_back(it->getGBRadius() * 0.1);
        // Screening parameters have been replaced
        switch(it->getElement()) {
            case 1:
                params.push_back(1.425952); // screen
                params.push_back(0.788440); // alpha
                params.push_back(0.798699); // beta
                params.push_back(0.437334); // gamma
                break;
            case 6:
                params.push_back(1.058554);
                params.push_back(0.733756);
                params.push_back(0.506378);
                params.push_back(0.205844);
                break;
            case 7:
                params.push_back(0.733599);
                params.push_back(0.503364);
                params.push_back(0.316828);
                params.push_back(0.192915);
                break;
            case 8:
                params.push_back(1.061039);
                params.push_back(0.867814);
                params.push_back(0.876635);
                params.push_back(0.387882);
                break;
            case 16:
                params.push_back(-0.703469);
                params.push_back(0.867814);
                params.push_back(0.876635);
                params.push_back(0.387882);
                break;
            default:
                params.push_back(0.5);
                params.push_back(1.0);
                params.push_back(0.8);
                params.push_back(4.85);
                break;
        }
        force->addParticle(params);
    }
    return force;
}

}; // namespace Amber
