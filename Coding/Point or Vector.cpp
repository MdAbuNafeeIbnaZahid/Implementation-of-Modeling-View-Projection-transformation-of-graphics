#include <bits/stdc++.h>
using namespace std;

#define SIZE 4
#define eps 1e-6


//struct PointOrPointOrVector
//{
//    double ar[SIZE];
//
//    PointOrPointOrVector()
//    {
//        FOR(a,0,SIZE-1)
//        {
//            ar[a] = 0;
//        }
//        ar[SIZE-1] = 1;
//    }
//
//    double getX() const
//    {
//        return ar[0];
//    }
//
//    double getY() const
//    {
//        return ar[1];
//    }
//
//    double getZ() const
//    {
//        return ar[2];
//    }
//
//    double setX( double x )
//    {
//        ar[0] = x;
//    }
//
//    double setY(double y)
//    {
//        ar[1] = y;
//    }
//
//    double setZ( double z )
//    {
//        ar[2] = z;
//    }
//
//    PointOrPointOrVector(double x, double y, double z)
//    {
//        ar[0] = x;
//        ar[1] = y;
//        ar[2] = z;
//
//        ar[3] = 1;
//    }
//
//    PointOrVector operator - (const PointOrPointOrVector &B) const
//    {
//        double x = this->getX() - B.getX();
//        double y = this->getY() - B.getY();
//        double z = this->getZ() - B.getZ();
//
//        PointOrVector ret(x,y,z);
//        return ret;
//    }
//};
//




struct PointOrVector
{
    int X_POS = 0;
    int Y_POS = 1;
    int Z_POS = 2;
    double ar[SIZE] = {0, 0, 0, 1};

    PointOrVector()
    {
        for (int a = 0; a < SIZE; a++)
        {
            ar[a] = 0;
        }
        ar[SIZE-1] = 1;
    }

    double getX() const
    {
        return ar[X_POS];
    }

    double getY() const
    {
        return ar[Y_POS];
    }

    double getZ() const
    {
        return ar[Z_POS];
    }

    double setX( double x )
    {
        ar[X_POS] = x;
    }

    double setY(double y)
    {
        ar[Y_POS] = y;
    }

    double setZ( double z )
    {
        ar[Z_POS] = z;
    }

    PointOrVector(double x, double y, double z)
    {
        setX(x);
        setY(y);
        setZ(z);
    }

    PointOrVector operator - (const PointOrVector &B) const
    {
        double x = this->getX() - B.getX();
        double y = this->getY() - B.getY();
        double z = this->getZ() - B.getZ();

        PointOrVector ret(x,y,z);
        return ret;
    }


    PointOrVector operator+(const PointOrVector &p) const
    {
        return PointOrVector(getX()+p.getX(), getY()+p.getY(), getZ()+p.getZ());
    }

    PointOrVector operator*(double i) const
    {
        return PointOrVector( getX()*i, getY()*i, getZ()*i );
    }

    PointOrVector operator+=(const PointOrVector &p)
    {
        this->setX( this->getX() + p.getX() );
        this->setY( this->getY() + p.getY() );
        this->setZ( this->getZ() + p.getZ() );


        return *this;
    }

    PointOrVector operator-=(const PointOrVector &p)
    {
        (*this) += p * (-1);

        return *this;
    }

    PointOrVector getCrossProduct(const PointOrVector &p) const
    {
        double x = this->getY() * p.getZ() - this->getZ() * p.getY();
        double y = this->getZ() * p.getX() - this->getX() * p.getZ();
        double z = this->getX() * p.getY() - this->getY() * p.getX();

        PointOrVector ret(x,y,z);
        return ret;
    }

    double getDisFromMain() const
    {
        return sqrt(getX()*getX() + getY()*getY() + getZ()*getZ() );
    }

    double getDotProduct(const PointOrVector &p) const
    {
        return (this->getX()) * p.getX() + (this->getY()) * p.getY() + (this->getZ()) * p.getZ();
    }

    bool isParallel( const PointOrVector &p ) const
    {
        PointOrVector crossVec = getCrossProduct(p);
        if ( crossVec.getDisFromMain() <= eps )
        {
            return true;
        }

        return false;
    }

    bool isUnit()
    {
        return ( abs(getDisFromMain() - 1) <= eps );
    }

    bool isPerpendicular( const PointOrVector &p ) const
    {
        return ( abs( getDotProduct(p) ) <= eps );
    }

    PointOrVector getNormalizedVec()
    {
        double disFromMain = getDisFromMain();
        PointOrVector copyVec(this->getX(), this->getY(), this->getZ());
        PointOrVector unitPointOrVector = copyVec * ( 1 / disFromMain );

        return unitPointOrVector;
    }

//    void rotateAroundPointOrVector( const PointOrVector &vec, double radAngle )
//    {
//        assert(false);
//
//        if ( isParallel(vec) )
//        {
//            return;
//        }
//
////        cout << vec.x << " " << vec.y << " " << vec.z << endl;
//
//        PointOrVector unitVec = this->getNormalizedVec();
//
//
////        cout << unitVec.x << " " << unitVec.y << " " << unitVec.z << endl;
//
////        cout << unitVec.getDisFromMain() << endl;
//        assert( ( "unit PointOrVector is not unit", unitVec.isUnit() ) ) ;
//
//        PointOrVector perpendicularUnitPointOrVector = unitVec;
//        assert( ( " Given PointOrVector is not perpendicular ", this->isPerpendicular(perpendicularUnitPointOrVector) )   );
//
//
//        PointOrVector crossPointOrVector = perpendicularUnitPointOrVector.getCrossProduct(*this);
//
//        this->x = this->x * cos(radAngle) + crossPointOrVector.x * sin(radAngle);
//        this->y = this->y * cos(radAngle) + crossPointOrVector.y * sin(radAngle);
//        this->z = this->z * cos(radAngle) + crossPointOrVector.z * sin(radAngle);
//
//    }



};
