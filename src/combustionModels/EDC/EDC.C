/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "EDC.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::EDC<ReactionThermo>::EDC
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    laminar<ReactionThermo>(modelType, thermo, turb, combustionProperties),
    version_
    (
        EDCversionNames
        [
            this->coeffs().lookupOrDefault
            (
                "version",
                word(EDCversionNames[EDCdefaultVersion])
            )
        ]
    ),
    C1_(this->coeffs().lookupOrDefault("C1", 0.05774)),
    C2_(this->coeffs().lookupOrDefault("C2", 0.5)),
    Cgamma_(this->coeffs().lookupOrDefault("Cgamma", 2.1377)),
    Ctau_(this->coeffs().lookupOrDefault("Ctau", 0.4083)),
    exp1_(this->coeffs().lookupOrDefault("exp1", EDCexp1[int(version_)])),
    exp2_(this->coeffs().lookupOrDefault("exp2", EDCexp2[int(version_)])),
    kappa_
    (
        IOobject
        (
            this->thermo().phasePropertyName(typeName + ":kappa"),
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("kappa", dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::EDC<ReactionThermo>::~EDC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Foam::combustionModels::EDC<ReactionThermo>::correct()
{
    if (this->active())
    {
        //Info << "    Nam --> in EDC: check #1" << endl; //
        tmp<volScalarField> tepsilon(this->turbulence().epsilon());
        const volScalarField& epsilon = tepsilon();

        tmp<volScalarField> tmu(this->turbulence().mu());
        const volScalarField& mu = tmu();

        tmp<volScalarField> tk(this->turbulence().k());
        const volScalarField& k = tk();

        tmp<volScalarField> trho(this->rho());
        const volScalarField& rho = trho();

        scalarField tauStar(epsilon.size(), 0);
        //Info << "    Nam --> in EDC: check #2" << endl; //

        if (version_ == EDCversions::v2016)
        {
        Info << "    Nam --> in EDC: check #3" << endl; //
            tmp<volScalarField> ttc(this->chemistryPtr_->tc());
            const volScalarField& tc = ttc();
        Info << "    Nam --> in EDC: check #4" << endl; //

            forAll(tauStar, i)
            {
        Info << "    Nam --> in EDC: check #5" << endl; //
                const scalar nu = mu[i]/(rho[i] + small);

                const scalar Da =
                    max(min(sqrt(nu/(epsilon[i] + small))/tc[i], 10), 1e-10);

                const scalar ReT = sqr(k[i])/(nu*epsilon[i] + small);
                const scalar CtauI = min(C1_/(Da*sqrt(ReT + 1)), 2.1377);

                const scalar CgammaI =
                    max(min(C2_*sqrt(Da*(ReT + 1)), 5), 0.4082);

                const scalar gammaL =
                    CgammaI*pow025(nu*epsilon[i]/(sqr(k[i]) + small));

                tauStar[i] = CtauI*sqrt(nu/(epsilon[i] + small));

        //Info << "    Nam --> in EDC: check #6" << endl; //
                if (gammaL >= 1)
                {
        //Info << "    Nam --> in EDC: check #7" << endl; //
                    kappa_[i] = 1;
                }
                else
                {
                    kappa_[i] =
                        max
                        (
                            min
                            (
                                pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                                1
                            ),
                            0
                        );
                }
            }
        }
        else
        {
        Info << "    Nam --> in EDC: check #8" << endl; //
            forAll(tauStar, i)
            {
        //Info << "    Nam --> in EDC: check #9" << endl; //
        //Info << "    after 9: mu[i] = " << mu[i]  << endl; //
        //Info << "    after 9: rho[i] = " << rho[i]  << endl; //
        //Info << "    after 9: small = " << small  << endl; //
        //Info << "    after 9: rho[i] + small = " << rho[i]+small  << endl; //

                const scalar nu = mu[i]/(rho[i] + small);
        //Info << "    Nam --> in EDC: check #9.2" << endl; //

        //Info << "    after 9.2: nu = " << nu  << endl; //
        //Info << "    after 9.2: epsilon[i] = " << epsilon[i]  << endl; //
        //Info << "    after 9.2: sqr(k[i]) = " << sqr(k[i])  << endl; //
        //Info << "    after 9.2: sqr(k[i]) + small = " << sqr(k[i])+small  << endl; //

                const scalar gammaL =
                    //Cgamma_*pow025(nu*epsilon[i]/(sqr(k[i]) + small));  //original
                    Cgamma_*pow(nu*epsilon[i]/(sqr(k[i]) + small), 0.25); //Nam
        //Info << "    Nam --> in EDC: check #10" << endl; //

                tauStar[i] = Ctau_*sqrt(nu/(epsilon[i] + small));
                if (gammaL >= 1)
                {
                    kappa_[i] = 1;
                }
                else
                {
                    kappa_[i] =
                        max
                        (
                            min
                            (
                                pow(gammaL, exp1_)/(1 - pow(gammaL, exp2_)),
                                1
                            ),
                            0
                        );
                }
            }
        }
        Info << "    Nam --> in EDC: check #11" << endl; //

        this->chemistryPtr_->solve(tauStar);
        Info << "    Nam --> in EDC: check #12" << endl; //
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::EDC<ReactionThermo>::R(volScalarField& Y) const
{
    return kappa_*laminar<ReactionThermo>::R(Y);
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::EDC<ReactionThermo>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        new volScalarField
        (
            IOobject
            (
                this->thermo().phasePropertyName(typeName + ":Qdot"),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->active())
    {
        tQdot.ref() = kappa_*this->chemistryPtr_->Qdot();
    }

    return tQdot;
}


template<class ReactionThermo>
bool Foam::combustionModels::EDC<ReactionThermo>::read()
{
    if (laminar<ReactionThermo>::read())
    {
        version_ =
        (
            EDCversionNames
            [
                this->coeffs().lookupOrDefault
                (
                    "version",
                    word(EDCversionNames[EDCdefaultVersion])
                )
            ]
        );
        C1_ = this->coeffs().lookupOrDefault("C1", 0.05774);
        C2_ = this->coeffs().lookupOrDefault("C2", 0.5);
        Cgamma_ = this->coeffs().lookupOrDefault("Cgamma", 2.1377);
        Ctau_ = this->coeffs().lookupOrDefault("Ctau", 0.4083);
        exp1_ = this->coeffs().lookupOrDefault("exp1", EDCexp1[int(version_)]);
        exp2_ = this->coeffs().lookupOrDefault("exp2", EDCexp2[int(version_)]);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
