/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "parabolicVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(0),
    n_(1, 0, 0),
    y_(0, 1, 0),
    x_(0, 0, 1),
    r_(1),
    center_(0, 0, 0)
{}


parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    maxValue_(ptf.maxValue_),
    n_(ptf.n_),
    y_(ptf.y_),
    x_(ptf.x_),
    r_(ptf.r_),
    center_(ptf.center_)
{}


parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    maxValue_(readScalar(dict.lookup("maxValue"))),
    n_(dict.lookup("n")),
    y_(dict.lookup("y")),
    x_(dict.lookup("x")),
    r_(readScalar(dict.lookup("r"))),
    center_(dict.lookup("center"))
{
    /*
    Info << "mag(x) = " << mag(x_) << "\n"
         << "mag(y) = " << mag(y_) << "\n"
         << "mag(n) = " << mag(n_) << "\n";
    */
    if (mag(n_) < SMALL || mag(y_) < SMALL || mag(x_) < SMALL)
    {

        FatalErrorIn("parabolicVelocityFvPatchVectorField(dict)")
            << "n or y or x given with zero size not correct"
            << abort(FatalError);
    }

    n_ /= mag(n_);
    y_ /= mag(y_);
    x_ /= mag(x_);

    evaluate();
}


parabolicVelocityFvPatchVectorField::parabolicVelocityFvPatchVectorField
(
    const parabolicVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    maxValue_(fcvpvf.maxValue_),
    n_(fcvpvf.n_),
    y_(fcvpvf.y_),
    x_(fcvpvf.x_),
    r_(fcvpvf.r_),
    center_(fcvpvf.center_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void parabolicVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //boundBox bb(patch().patch().localPoints(), true);
    vectorField c = patch().Cf();

    //vector ctr = 0.5*(bb.min() + bb.max());

    scalarField coord1 = ((c - center_) & x_)/(r_);
    scalarField coord2 = ((c - center_) & y_)/(r_);

    scalarField coord = sqrt(sqr(coord1) + sqr(coord2));

    vectorField::operator=(n_*(1.06603 - 0.03576 * coord          - 16757.3   * pow(coord, 2) 
                                       - 70202.6 * pow(coord, 3)  - 369192648 * pow(coord, 4)));
}

/*
scalarField computeLocalCoordForBasisVector(const boundBox & bb, const vector & basisVector, const vectorField & c)
{
    vector ctr = 0.5*(bb.max() + bb.min());
    scalarField coord = 2*((c - ctr) & basisVector)/((bb.max() - bb.min()) & basisVector);

    return coord;
}
*/


// Write
void parabolicVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("maxValue")
        << maxValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("y")
        << y_ << token::END_STATEMENT << nl;
    os.writeKeyword("x")
        << x_ << token::END_STATEMENT << nl;
    os.writeKeyword("r")
        << r_ << token::END_STATEMENT << nl;
    os.writeKeyword("center")
        << center_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, parabolicVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
