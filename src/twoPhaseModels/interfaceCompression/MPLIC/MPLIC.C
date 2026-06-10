/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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

#include "MPLIC.H"
#include "MPLICcell.H"
#include "volPointInterpolation.H"
#include "syncTools.H"
#include "slicedSurfaceFields.H"
#include "upwind.H"
#include "intersectionPatchToPatch.H"
#include "nonConformalCouplePolygons.H"
#include "nonConformalFvPatch.H"
#include "nonConformalCoupledFvPatch.H"
#include "nonConformalCyclicFvPatch.H"
#include "Map.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MPLIC, 0);

    surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<MPLIC>
        addMPLICScalarMeshFluxConstructorToTable_;

    //- Enable retention of the non-conformal coupling intersection polygons in
    //  the mesh stitcher. This is needed so that the geometric interface
    //  reconstruction can cut the coupling faces. Loading this library (i.e.,
    //  using a geometric VOF scheme) is the signal to pay the storage cost;
    //  pure conformal or non-VOF non-conformal runs are unaffected.
    static const bool enableStorePolygons_
    (
        patchToPatches::intersection::storePolygons_ = true
    );
}


// * * * * * * * * * * * * * Local Helper Functions  * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Area vector of a (near-planar, convex) polygon loop, by fan triangulation
inline vector loopAreaVector(const List<point>& ps)
{
    vector a = Zero;
    if (ps.size() < 3) return a;
    for (label i = 1; i < ps.size() - 1; ++ i)
    {
        a += 0.5*((ps[i] - ps[0])^(ps[i + 1] - ps[0]));
    }
    return a;
}


//- Submerged (liquid) area fraction of a set of polygon loops under the PLIC
//  plane through base b with interface area-vector direction n. The liquid
//  side is taken to be where (x - b) & n <= 0, i.e. n points from the liquid
//  to the gas, consistent with the cut surface area vector in MPLICcell.
//  Returns -1 if the total polygon area is negligible.
scalar submergedFraction
(
    const List<List<point>>& loops,
    const point& b,
    const vector& n
)
{
    scalar totalArea = 0;
    scalar wetArea = 0;

    forAll(loops, li)
    {
        const List<point>& loop = loops[li];
        if (loop.size() < 3) continue;

        totalArea += mag(loopAreaVector(loop));

        // Clip the loop against the half-space (x - b) & n <= 0
        DynamicList<point> wet(loop.size() + 4);
        forAll(loop, i)
        {
            const point& cur = loop[i];
            const point& nxt = loop[loop.fcIndex(i)];

            const scalar dc = (cur - b) & n;
            const scalar dn = (nxt - b) & n;

            if (dc <= 0) wet.append(cur);

            if ((dc < 0) != (dn < 0))
            {
                const scalar t = dc/(dc - dn);
                wet.append(cur + t*(nxt - cur));
            }
        }

        if (wet.size() >= 3) wetArea += mag(loopAreaVector(wet));
    }

    if (totalArea < vSmall) return -1;

    return min(max(wetArea/totalArea, scalar(0)), scalar(1));
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

void Foam::MPLIC::setCellAlphaf
(
    const label celli,
    const scalarField& phi,
    scalarField& alphaf,
    boolList& correctedFaces,
    const DynamicList<scalar>& cellAlphaf,
    const fvMesh& mesh
) const
{
    // Face owners reference
    const labelList& own = mesh.faceOwner();

    // The cell face labels
    const labelList& cFaces = mesh.cells()[celli];

    // Fill the alphaf with surface interpolation in direction of the flow
    forAll(cFaces, i)
    {
        const label facei = cFaces[i];
        const scalar phiSigni = sign(phi[facei]);

        if
        (
            (own[facei] == celli && phiSigni == 1)
         || (own[facei] != celli && phiSigni == -1)
        )
        {
            alphaf[facei] = cellAlphaf[i];
            correctedFaces[facei] = true;
        }
    }
}


void Foam::MPLIC::syncNonConformalAlphaf
(
    surfaceScalarField& alphaf,
    const surfaceScalarField& phi
) const
{
    const fvMesh& mesh = this->mesh();

    // Snapshot the owner-side reconstruction so the neighbour reads below are
    // unaffected by the in-place updates
    List<scalarField> alphaf0(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        alphaf0[patchi] = alphaf.boundaryField()[patchi];
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];

        // Same-mesh cyclic coupling. The processor-cyclic and mapped-wall
        // couplings are not yet reconciled here and retain their owner-side
        // reconstruction.
        if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

        const nonConformalCyclicFvPatch& ncc =
            refCast<const nonConformalCyclicFvPatch>(fvp);

        const label nbri = ncc.nbrPatch().index();

        const scalarField& ownAlpha = alphaf0[patchi];
        const scalarField& nbrAlpha = alphaf0[nbri];
        const scalarField& phip = phi.boundaryField()[patchi];

        scalarField& pf = alphaf.boundaryFieldRef()[patchi];

        forAll(pf, facei)
        {
            // Upwind selection: flux leaving this side (phi >= 0) keeps the own
            // donor value; flux entering (phi < 0) takes the neighbour donor
            // value. With phi_own == -phi_nbr across the couple, both sides end
            // up equal to the donor value, so the advective flux is
            // conservative. The phase fraction is a scalar, so the coupling
            // transform is the identity.
            pf[facei] = phip[facei] >= 0 ? ownAlpha[facei] : nbrAlpha[facei];
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::MPLIC::surfaceAlpha
(
    const volScalarField& alpha,
    const surfaceScalarField& phi,
    scalarField& initAlphaf,
    const bool unweighted,
    const scalar tol,
    const bool isMPLIC
) const
{
    // Finite volume mesh reference
    const fvMesh& mesh = alpha.mesh();

    // Reference to primitive mesh
    const primitiveMesh& primMesh = mesh;

    // Velocity field reference
    const volVectorField& U
    (
        mesh.lookupObject<const volVectorField>
        (
            IOobject::groupName("U", phi.group())
        )
    );

    // Interpolate alpha from volume to the points of the mesh
    const scalarField alphap
    (
        volPointInterpolation::New(mesh).interpolate(alpha)
    );

    // Interpolate U from cell centres to the points of the mesh
    vectorField Up;

    if (!unweighted)
    {
        Up = volPointInterpolation::New(mesh).interpolate(U);
    }

    // Flatten down phi flux field
    const scalarField splicedPhi
    (
        slicedSurfaceScalarField
        (
            IOobject
            (
                "splicedPhi",
                mesh.time().name(),
                mesh
            ),
            phi,
            false
        ).splice()
    );

    scalarField alphaf(mesh.nFaces(), 0);

    // Mark which faces are corrected by MPLIC
    boolList correctedFaces(mesh.nFaces(), false);

    // Look up the non-conformal coupling polygon store. This is present only
    // when the mesh is non-conformal and a geometric VOF scheme is in use. When
    // present, the coupling faces are reconstructed geometrically below.
    const nonConformalCouplePolygons* nccPolysPtr =
        mesh.foundObject<nonConformalCouplePolygons>
        (
            nonConformalCouplePolygons::typeName
        )
      ? &mesh.lookupObject<nonConformalCouplePolygons>
        (
            nonConformalCouplePolygons::typeName
        )
      : nullptr;

    // Identify the cells adjacent to non-conformal coupling faces so that their
    // reconstructed interface planes can be cached for cutting the coupling
    // faces after the main cell loop.
    labelHashSet nccOwnerCells;
    if (nccPolysPtr)
    {
        forAll(mesh.boundary(), patchi)
        {
            if (isA<nonConformalCoupledFvPatch>(mesh.boundary()[patchi]))
            {
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, i) nccOwnerCells.insert(fc[i]);
            }
        }
    }
    Map<vector> cellCutNormal;
    Map<point> cellCutBase;

    // Construct class for cell cut
    MPLICcell cutCell(unweighted, isMPLIC);

    // Loop through all the cells
    forAll(mesh.cells(), celli)
    {
        if (alpha[celli] < (1 - tol) && alpha[celli] > tol)
        {
            // Store cell information
            const MPLICcellStorage cellInfo
            (
                primMesh,
                alphap,
                Up,
                alpha[celli],
                U[celli],
                celli
            );

            // Volume ratio matching algorithm
            if (cutCell.matchAlpha(cellInfo))
            {
                // Fill cutCell.alphaf() with face values from this cell
                setCellAlphaf
                (
                    celli,
                    splicedPhi,
                    alphaf,
                    correctedFaces,
                    cutCell.alphaf(),
                    mesh
                );

                // Cache the reconstructed interface plane for cells bordering a
                // non-conformal coupling so the coupling faces can be cut
                if (nccPolysPtr && nccOwnerCells.found(celli))
                {
                    cellCutNormal.insert(celli, cutCell.cutNormal());
                    cellCutBase.insert(celli, cutCell.cutBase());
                }
            }
        }
    }

    // Synchronise across the processor and cyclic patches
    syncTools::syncFaceList(mesh, alphaf, plusEqOp<scalar>());
    syncTools::syncFaceList(mesh, correctedFaces, plusEqOp<bool>());

    // Correct selected faces
    forAll(correctedFaces, facei)
    {
        if (correctedFaces[facei])
        {
            initAlphaf[facei] = alphaf[facei];
        }
    }

    // Convert the spliced field into a surfaceScalarField. On conformal meshes
    // this is a direct slice of the poly-face-indexed initAlphaf. On
    // non-conformal meshes the poly and finite-volume face addressing differ -
    // the coupling faces have no poly start face - so the field is assembled
    // directly: the internal and conformal-boundary faces are copied from
    // initAlphaf and the coupling faces are reconstructed below.
    tmp<surfaceScalarField> tsplicedAlpha;
    if (mesh.conformal())
    {
        tsplicedAlpha =
            surfaceScalarField::New
            (
                "alphaf",
                slicedSurfaceScalarField
                (
                    IOobject
                    (
                        "alphaf",
                        mesh.time().name(),
                        mesh
                    ),
                    mesh,
                    dimless,
                    initAlphaf,
                    false
                ),
                fvsPatchField<scalar>::calculatedType()
            );
    }
    else
    {
        tsplicedAlpha =
            surfaceScalarField::New
            (
                "alphaf",
                mesh,
                dimensionedScalar(dimless, 0),
                fvsPatchField<scalar>::calculatedType()
            );
        surfaceScalarField& a = tsplicedAlpha.ref();

        // Internal faces (poly and fv internal addressing coincide)
        scalarField& ai = a.primitiveFieldRef();
        forAll(ai, facei)
        {
            ai[facei] = initAlphaf[facei];
        }

        // Conformal boundary patches map directly to their poly faces. The
        // non-conformal patches have no poly start and are left zero here; the
        // coupled ones are reconstructed below.
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];

            if (isA<nonConformalFvPatch>(fvp)) continue;

            scalarField& pf = a.boundaryFieldRef()[patchi];
            const label start = fvp.start();
            forAll(pf, facei)
            {
                pf[facei] = initAlphaf[start + facei];
            }
        }
    }
    surfaceScalarField& splicedAlpha = tsplicedAlpha.ref();

    forAll(mesh.boundary(), patchi)
    {
        if (alpha.boundaryField()[patchi].fixesValue())
        {
            splicedAlpha.boundaryFieldRef()[patchi] =
                alpha.boundaryField()[patchi];
        }
    }

    // Geometrically reconstruct the phase fraction on the non-conformal
    // coupling faces. Each coupling face is cut by the reconstructed interface
    // plane of its adjacent (owner) cell; the submerged area fraction gives the
    // face phase fraction. Faces whose owner cell is fully wet/dry, or whose
    // owner cell could not be cut, take the owner cell value. The flux-upwind
    // selection across the coupling and the matching of the two coupled sides
    // are handled subsequently by syncNonConformalAlphaf.
    if (nccPolysPtr)
    {
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];

            if (!isA<nonConformalCoupledFvPatch>(fvp)) continue;

            const labelUList& fc = fvp.faceCells();
            const nonConformalCouplePolygons::patchPolygons& patchPolys =
                (*nccPolysPtr)[patchi];

            scalarField& pf = splicedAlpha.boundaryFieldRef()[patchi];

            forAll(fc, facei)
            {
                const label celli = fc[facei];

                scalar af;
                if (alpha[celli] <= tol)
                {
                    af = 0;
                }
                else if (alpha[celli] >= 1 - tol)
                {
                    af = 1;
                }
                else if
                (
                    cellCutBase.found(celli)
                 && facei < patchPolys.size()
                )
                {
                    const scalar frac =
                        submergedFraction
                        (
                            patchPolys[facei],
                            cellCutBase[celli],
                            cellCutNormal[celli]
                        );

                    af = frac >= 0 ? frac : alpha[celli];
                }
                else
                {
                    af = alpha[celli];
                }

                pf[facei] = af;
            }
        }

        // Reconcile the two coupled sides and apply flux upwinding
        syncNonConformalAlphaf(splicedAlpha, phi);
    }

    return tsplicedAlpha;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::MPLIC::interpolate
(
    const VolField<scalar>& vf
) const
{
    tmp<surfaceScalarField> tvff(upwind<scalar>(mesh(), phi_).interpolate(vf));

    scalarField splicedTvff
    (
        slicedSurfaceScalarField
        (
            IOobject
            (
                "splicedTvff",
                mesh().time().name(),
                mesh()
            ),
            tvff,
            false
        ).splice()
    );

    return surfaceAlpha(vf, phi_, splicedTvff, true, 1e-6);
}


// ************************************************************************* //
