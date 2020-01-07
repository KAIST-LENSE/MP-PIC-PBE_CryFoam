/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "KinematicParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::KinematicParcel<ParcelType>::propertyList_ =
    Foam::KinematicParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::KinematicParcel<ParcelType>::sizeofFields_
(
    offsetof(KinematicParcel<ParcelType>, rhoc_)
  - offsetof(KinematicParcel<ParcelType>, active_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::KinematicParcel<ParcelType>::KinematicParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    active_(false),
    typeId_(0),
    nParticle_(0.0),
    d_(0.0),
    dTarget_(0.0),
    U_(Zero),
    rho_(0.0),
    age_(0.0),
    tTurb_(0.0),
    UTurb_(Zero),
    kv_(0.0),
    dr_(0.0),
    noNode_(0.0),
    minmodeTheta_(0.0), 
    totalSource_(0.0),
    rhoc_(0.0),
    Uc_(Zero),
    muc_(0.0),
    supersatc_(0.0),
    R_growthc_(0.0),
    R_nucleationc_(0.0)
{
    for(int i=0 ; i<noNode_ ; i++)
    {
        editF(i) = 0.0 ;
    }
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            active_ = readBool(is);
            typeId_ = readLabel(is);
            nParticle_ = readScalar(is);
            d_ = readScalar(is);
            dTarget_ = readScalar(is);
            is >> U_;
            rho_ = readScalar(is);
            age_ = readScalar(is);
            tTurb_ = readScalar(is);
            is >> UTurb_;
            kv_ = readScalar(is);
            dr_ = readScalar(is);
            noNode_ = readScalar(is);
            for(int i=0 ; i<noNode_ ; i++)
            {
                F_[i] = readScalar(is);
            }
            minmodeTheta_ = readScalar(is);
            totalSource_ = readScalar(is);
            rhoc_ = readScalar(is);
            is >> Uc_;
            muc_ = readScalar(is);
            supersatc_ = readScalar(is);
            R_growthc_ = readScalar(is);
            R_nucleationc_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&active_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "KinematicParcel<ParcelType>::KinematicParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<label> active
    (
        c.fieldIOobject("active", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, active);

    IOField<label> typeId
    (
        c.fieldIOobject("typeId", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, nParticle);

    IOField<scalar> d
    (
        c.fieldIOobject("d", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, d);

    IOField<scalar> dTarget
    (
        c.fieldIOobject("dTarget", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, dTarget);

    IOField<vector> U
    (
        c.fieldIOobject("U", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, U);

    IOField<scalar> rho
    (
        c.fieldIOobject("rho", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, rho);

    IOField<scalar> age
    (
        c.fieldIOobject("age", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, age);

    IOField<scalar> tTurb
    (
        c.fieldIOobject("tTurb", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, tTurb);

    IOField<vector> UTurb
    (
        c.fieldIOobject("UTurb", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, UTurb);

    IOField<scalarField> F
    (
        c.fieldIOobject("F", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, F);

    IOField<scalar> kv
    (
        c.fieldIOobject("kv", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, kv);

    IOField<scalar> dr
    (
        c.fieldIOobject("dr", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, dr);

    IOField<scalar> noNode
    (
        c.fieldIOobject("noNode", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, noNode);

    IOField<scalar> minmodeTheta
    (
        c.fieldIOobject("minmodeTheta", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, minmodeTheta);

    IOField<scalar> totalSource
    (
        c.fieldIOobject("totalSource", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, totalSource);

    IOField<scalar> rhoc
    (
        c.fieldIOobject("rhoc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, rhoc);

    IOField<vector> Uc
    (
        c.fieldIOobject("Uc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, Uc);

    IOField<scalar> muc
    (
        c.fieldIOobject("muc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, muc);

    IOField<scalar> supersatc
    (
        c.fieldIOobject("supersatc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, supersatc);


    IOField<scalar> R_growthc
    (
        c.fieldIOobject("R_growthc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, R_growthc);

    IOField<scalar> R_nucleationc
    (
        c.fieldIOobject("R_nucleationc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, R_nucleationc);


    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        KinematicParcel<ParcelType>& p = iter();

        p.active_ = active[i];
        p.typeId_ = typeId[i];
        p.nParticle_ = nParticle[i];
        p.d_ = d[i];
        p.dTarget_ = dTarget[i];
        p.U_ = U[i];
        p.rho_ = rho[i];
        p.age_ = age[i];
        p.tTurb_ = tTurb[i];
        p.UTurb_ = UTurb[i];
        p.kv_ = kv[i];
        p.dr_ = dr[i];
        p.noNode_ = noNode[i];
        for(int j=0 ; j<p.noNode() ; j++)
        {
            p.F_[j] = F[i][j];
        } 
        p.minmodeTheta_ = minmodeTheta[i];
        p.totalSource_ = totalSource[i];
        p.rhoc_ = rhoc[i];
        p.Uc_ = Uc[i];
        p.muc_ = muc[i];
        p.supersatc_ = supersatc[i];
        p.R_growthc_ = R_growthc[i];
        p.R_nucleationc_ = R_nucleationc[i];

        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::KinematicParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<label> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> nParticle
    (
        c.fieldIOobject("nParticle", IOobject::NO_READ),
        np
    );
    IOField<scalar> d(c.fieldIOobject("d", IOobject::NO_READ), np);
    IOField<scalar> dTarget(c.fieldIOobject("dTarget", IOobject::NO_READ), np);
    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> rho(c.fieldIOobject("rho", IOobject::NO_READ), np);
    IOField<scalar> age(c.fieldIOobject("age", IOobject::NO_READ), np);
    IOField<scalar> tTurb(c.fieldIOobject("tTurb", IOobject::NO_READ), np);
    IOField<vector> UTurb(c.fieldIOobject("UTurb", IOobject::NO_READ), np);
    IOField<scalarField> F(c.fieldIOobject("F", IOobject::NO_READ), np);
    IOField<scalar> kv(c.fieldIOobject("kv", IOobject::NO_READ), np);
    IOField<scalar> dr(c.fieldIOobject("dr", IOobject::NO_READ), np);
    IOField<scalar> noNode(c.fieldIOobject("noNode", IOobject::NO_READ), np);
    IOField<scalar> minmodeTheta(c.fieldIOobject("minmodeTheta", IOobject::NO_READ), np);
    IOField<scalar> totalSource(c.fieldIOobject("totalSource", IOobject::NO_READ), np);
    IOField<scalar> rhoc(c.fieldIOobject("rhoc", IOobject::NO_READ), np);
    IOField<vector> Uc(c.fieldIOobject("Uc", IOobject::NO_READ), np);
    IOField<scalar> muc(c.fieldIOobject("muc", IOobject::NO_READ), np);
    IOField<scalar> supersatc(c.fieldIOobject("supersatc", IOobject::NO_READ), np);
    IOField<scalar> R_growthc(c.fieldIOobject("R_growthc", IOobject::NO_READ), np);
    IOField<scalar> R_nucleationc(c.fieldIOobject("R_nucleationc", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const KinematicParcel<ParcelType>& p = iter();
        scalarField dataF(p.noNode());
        for(int j=0 ; j<p.noNode_ ; j++)
        {
            dataF[j] = p.F(j) ;
        }
        active[i] = p.active();
        typeId[i] = p.typeId();
        nParticle[i] = p.nParticle();
        d[i] = p.d();
        dTarget[i] = p.dTarget();
        U[i] = p.U();
        rho[i] = p.rho();
        age[i] = p.age();
        tTurb[i] = p.tTurb();
        UTurb[i] = p.UTurb();
        F[i] = dataF;
        kv[i] = p.kv();
        dr[i] = p.dr();
        noNode[i] = p.noNode();
        minmodeTheta[i] = p.minmodeTheta();
        totalSource[i] = p.totalSource();
        rhoc[i] = p.rhoc();
        Uc[i] = p.Uc();
        muc[i] = p.muc();
        supersatc[i] = p.supersatc();
        R_growthc[i] = p.R_growthc();
        R_nucleationc[i] = p.R_nucleationc();

        i++;
    }

    const bool valid = np > 0;

    active.write(valid);
    typeId.write(valid);
    nParticle.write(valid);
    d.write(valid);
    dTarget.write(valid);
    U.write(valid);
    rho.write(valid);
    age.write(valid);
    tTurb.write(valid);
    UTurb.write(valid);
    F.write(valid);
    kv.write(valid);
    dr.write(valid);
    noNode.write(valid);
    minmodeTheta.write(valid);
    totalSource.write(valid);
    rhoc.write(valid);
    Uc.write(valid);
    muc.write(valid);
    supersatc.write(valid);
    R_growthc.write(valid);
    R_nucleationc.write(valid);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const KinematicParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.active()
            << token::SPACE << p.typeId()
            << token::SPACE << p.nParticle()
            << token::SPACE << p.d()
            << token::SPACE << p.dTarget()
            << token::SPACE << p.U()
            << token::SPACE << p.rho()
            << token::SPACE << p.age()
            << token::SPACE << p.tTurb()
            << token::SPACE << p.UTurb()
            << token::SPACE << p.kv()
            << token::SPACE << p.dr()
            << token::SPACE << p.noNode()
            << token::SPACE << p.minmodeTheta()
            << token::SPACE << p.totalSource()
            << token::SPACE << p.rhoc()
            << token::SPACE << p.Uc()
            << token::SPACE << p.muc()
            << token::SPACE << p.supersatc()
            << token::SPACE << p.R_growthc()
            << token::SPACE << p.R_nucleationc();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.active_),
            KinematicParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const KinematicParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
