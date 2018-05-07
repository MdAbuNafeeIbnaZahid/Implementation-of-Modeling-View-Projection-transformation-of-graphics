struct Vector
{
	double x,y,z;
    Vector(){}
	Vector(double x, double y, double z)
	{
        this->x = x;
        this->y = y;
        this->z = z;
	}

    Vector operator+(const Vector &p) const
    {
        return Vector(x+p.x, y+p.y, z+p.z);
    }

    Vector operator*(double i) const
    {
        return Vector( x*i, y*i, z*i );
    }

    Vector operator+=(const Vector &p)
    {
        this->x += p.x;
        this->y += p.y;
        this->z += p.z;

        return *this;
    }

    Vector operator-=(const Vector &p)
    {
        (*this) += p * (-1);

        return *this;
    }

    Vector getCrossProduct(const Vector &p) const
    {
        Vector ret;

        ret.x = this->y * p.z - this->z * p.y;
        ret.y = this->z * p.x - this->x * p.z;
        ret.z = this->x * p.y - this->y * p.x;

        return ret;
    }

    double getDisFromMain() const
    {
        return sqrt(x*x + y*y + z*z);
    }

    double getDotProduct(const Vector &p) const
    {
        return (this->x) * p.x + (this->y) * p.y + (this->z) * p.z;
    }

    bool isParallel( const Vector &p ) const
    {
        Vector crossVec = getCrossProduct(p);
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

    bool isPerpendicular( const Vector &p ) const
    {
        return ( abs( getDotProduct(p) ) <= eps );
    }

    Vector getNormalizedVec()
    {
        double disFromMain = getDisFromMain();
        Vector copyVec(this.x, this.y, this.z);
        Vector unitVector = copyVec * ( 1 / disFromMain );

        return unitVector;
    }

    void rotateAroundVector( const Vector &vec, double radAngle )
    {

        if ( isParallel(vec) )
        {
            return;
        }

//        cout << vec.x << " " << vec.y << " " << vec.z << endl;

        unitVec = this.getNormalizedVec();


//        cout << unitVec.x << " " << unitVec.y << " " << unitVec.z << endl;

//        cout << unitVec.getDisFromMain() << endl;
        assert( ( "unit vector is not unit", unitVec.isUnit() ) ) ;

        Vector perpendicularUnitVector = unitVec;
        assert( ( " Given vector is not perpendicular ", this->isPerpendicular(perpendicularUnitVector) )   );


        Vector crossVector = perpendicularUnitVector.getCrossProduct(*this);

        this->x = this->x * cos(radAngle) + crossVector.x * sin(radAngle);
        this->y = this->y * cos(radAngle) + crossVector.y * sin(radAngle);
        this->z = this->z * cos(radAngle) + crossVector.z * sin(radAngle);

    }



};


const Vector unitXVec(1,0,0);
const Vector unitYVec(0,1,0);
const Vector unitZVec(0,0,1);
