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

Class
    Foam::thermophysicalPropertiesSelector

Description
    Wrapper class providing run-time selection of thermophysicalProperties
    for the templated thermodynamics packages.

SourceFiles
    thermophysicalPropertiesSelectorI.H
    thermophysicalPropertiesSelector.C

\*---------------------------------------------------------------------------*/

#ifndef thermophysicalPropertiesSelector_H
#define thermophysicalPropertiesSelector_H

#include "thermophysicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class thermophysicalPropertiesSelector Declaration
\*---------------------------------------------------------------------------*/

template<class ThermophysicalProperties>
class thermophysicalPropertiesSelector
{
    // Private member data

        autoPtr<ThermophysicalProperties> propertiesPtr_;


public:

    // Constructors

        //- Construct from name
        thermophysicalPropertiesSelector(const word& name);

        //- Construct from dictionary
        thermophysicalPropertiesSelector(const dictionary& dict);


    // Static data

        //- Is the equation of state is incompressible i.e. rho != f(p)
        static const bool incompressible =
            ThermophysicalProperties::incompressible;

        //- Is the equation of state is isochoric i.e. rho = const
        static const bool isochoric =
            ThermophysicalProperties::isochoric;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "thermophysicalPropertiesSelector<"
              + ThermophysicalProperties::typeName
              + '>';
        }

        //- Return reference to the selected physical properties class
        inline const ThermophysicalProperties& properties() const;


        // Physical constants which define the specie

            //- Molecular weight [kg/kmol]
            inline scalar W() const;

            //- Limit the temperature to be in the range Tlow_ to Thigh_
            inline scalar limit(const scalar T) const;


        // Fundamental equation of state properties

            //- Liquid density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Liquid compressibility rho/p [s^2/m^2]
            //  Note: currently it is assumed the liquid is incompressible
            inline scalar psi(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            //  Note: currently it is assumed the liquid is incompressible
            //  so CpMCv 0
            inline scalar CpMCv(scalar p, scalar T) const;


        // Fundamental thermodynamic properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //Nam
            //- Derivative of Sensible enthalpy [J/kgK]
            inline scalar dHsdT(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            // Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


        // Physical properties

            //- Liquid viscosity [Pa s]
            inline scalar mu(scalar p, scalar T) const;

            //- Liquid thermal conductivity  [W/(m K)]
            inline scalar kappa(scalar p, scalar T) const;

            //- Liquid thermal diffusivity of enthalpy [kg/ms]
            inline scalar alphah(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermophysicalPropertiesSelectorI.H"

#ifdef NoRepository
    #include "thermophysicalPropertiesSelector.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
