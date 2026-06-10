/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::nonConformalCouplePolygons

Description
    A registered store of the intersection polygons of the non-conformal
    coupling (NCC) faces, indexed [patchi][patchFacei] as a list of point loops
    in the frame of the cell adjacent to the patch.

    These polygons are not part of the surface geometry fields (which only carry
    an area vector and centre per coupling face), but are needed by geometric
    VOF interface-reconstruction schemes (interfaceCompression PLIC/MPLIC) so
    that the coupling faces can be cut against the reconstructed interface.

    The store is populated by the fvMeshStitcher when the mesh is made
    non-conformal (unconform) and is removed when it is made conformal again
    (conform). It is only created when intersection::storePolygons_ is enabled,
    i.e., when a geometric VOF library is loaded.

SourceFiles
    nonConformalCouplePolygons.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalCouplePolygons_H
#define nonConformalCouplePolygons_H

#include "regIOobject.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                  Class nonConformalCouplePolygons Declaration
\*---------------------------------------------------------------------------*/

class nonConformalCouplePolygons
:
    public regIOobject
{
public:

    // Public Typedefs

        //- The polygons of a single coupling face; a list of point loops
        typedef List<List<point>> facePolygons;

        //- The polygons of all faces of a patch
        typedef List<facePolygons> patchPolygons;

        //- The polygons of all the boundary patches
        typedef List<patchPolygons> boundaryPolygons;


private:

    // Private Data

        //- The coupling face polygons, indexed [patchi][patchFacei]. Patches
        //  that are not non-conformal coupled, and faces that are not coupled,
        //  hold an empty list.
        boundaryPolygons boundary_;


public:

    //- Runtime type information
    TypeName("nonConformalCouplePolygons");


    // Constructors

        //- Construct empty, sized to the mesh boundary, and register
        explicit nonConformalCouplePolygons(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        nonConformalCouplePolygons(const nonConformalCouplePolygons&) = delete;


    //- Destructor
    virtual ~nonConformalCouplePolygons();


    // Member Functions

        //- Access the boundary polygons
        const boundaryPolygons& boundary() const
        {
            return boundary_;
        }

        //- Access the boundary polygons
        boundaryPolygons& boundary()
        {
            return boundary_;
        }

        //- Access the polygons of a patch
        const patchPolygons& operator[](const label patchi) const
        {
            return boundary_[patchi];
        }

        //- Dummy write. The store is transient and never written.
        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nonConformalCouplePolygons&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
