#define DEBUG




#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp> // Common file
#include <ext/pb_ds/tree_policy.hpp> // Including tree_order_statistics_node_update
#include <ext/pb_ds/detail/standard_policies.hpp>
using namespace std;
using namespace __gnu_pbds;
using namespace __gnu_cxx;
// Order Statistic Tree
/* Special functions:

		find_by_order(k) --> returns iterator to the kth largest element counting from 0
		order_of_key(val) --> returns the number of items in a set that are strictly smaller than our item
*/
typedef tree<
int,
null_type,
less<int>,
rb_tree_tag,
tree_order_statistics_node_update>
ordered_set;



// Generic Triple with stream inserter
template<class Typ1, class Typ2, class Typ3> class Triple
{
public:
    Typ1 first;
    Typ2 second;
    Typ3 third;
    Triple(){}
    Triple(Typ1 first, Typ2 second, Typ3 third)
    {
        this->first = first;
        this->second = second;
        this->third = third;
    }
    bool operator < ( const Triple B ) const
    {
        if ( first == B.first )
        {
            if ( second == B.second )
            {
                return third < B.third;
            }
            return second < B.second;
        }
        return first < B.first;
    }

};
template<class Typ1, class Typ2, class Typ3>
ostream &operator << (ostream &stream, const Triple<Typ1, Typ2, Typ3> tr)
{
    stream << " (" << tr.first << ", " << tr.second << ", " << tr.third << ") ";
    return stream;
}




/// ********* debug template by Bidhan Roy *********

template < typename F, typename S >
ostream& operator << ( ostream& os, const pair< F, S > & p ) {
    return os << "(" << p.first << ", " << p.second << ")";
}

template < typename T >
ostream &operator << ( ostream & os, const vector< T > &v ) {
    os << "{";
    typename vector< T > :: const_iterator it;
    for( it = v.begin(); it != v.end(); it++ ) {
        if( it != v.begin() ) os << ", ";
        os << *it;
    }
    return os << "}";
}

template < typename T >
ostream &operator << ( ostream & os, const set< T > &v ) {
    os << "[";
    typename set< T > :: const_iterator it;
    for ( it = v.begin(); it != v.end(); it++ ) {
        if( it != v.begin() ) os << ", ";
        os << *it;
    }
    return os << "]";
}

template < typename F, typename S >
ostream &operator << ( ostream & os, const map< F, S > &v ) {
    os << "[";
    typename map< F , S >::const_iterator it;
    for( it = v.begin(); it != v.end(); it++ ) {
        if( it != v.begin() ) os << ", ";
        os << it -> first << " = " << it -> second ;
    }
    return os << "]";
}


// This debugger for Multiset is written by Nafee
template < typename T >
ostream &operator << ( ostream & os, const multiset< T > &v ) {
    os << "[";
    typename multiset< T > :: const_iterator it;
    for ( it = v.begin(); it != v.end(); it++ ) {
        if( it != v.begin() ) os << ", ";
        os << *it;
    }
    return os << "]";
}


// This debugger for queue is written by Nafee
template < typename T >
ostream &operator << ( ostream & os, queue< T > v ) {
    os << "[";
    bool isFirst = true;
    while(v.size())
    {
        T cur = v.front();
        v.pop();
        if ( !isFirst )
        {
            os << ", ";
        }
        os << cur;
        isFirst = false;
    }
    return os << "]";
}



// This debugger for Deque is written by Nafee
template < typename T >
ostream &operator << ( ostream & os, const deque< T > &v ) {
    os << "[";
    typename deque< T > :: const_iterator it;
    for ( it = v.begin(); it != v.end(); it++ ) {
        if( it != v.begin() ) os << ", ";
        os << *it;
    }
    return os << "]";
}


// This debugger for array is written by Nafee
template<typename T1>
void prAr(T1 *ar, long long si, long long ei)
{
    long long a, b, c;
    for (a = si; a <= ei; a++)
    {
        cout << "ar [ " << a << " ] = " << ar[a] << endl;
    }
}


// This debugger for unordered_set is written by Nafee
template < typename T >
ostream &operator << ( ostream & os, const unordered_set< T > &v ) {
    os << "[";
    typename unordered_set< T > :: const_iterator it;
    for ( it = v.begin(); it != v.end(); it++ ) {
        if( it != v.begin() ) os << ", ";
        os << *it;
    }
    return os << "]";
}


// This debugger for stack is written by Nafee
template < typename T >
ostream &operator << ( ostream & os, stack< T > stak ) {
    os << "[";
    while(stak.size())
    {
        os << stak.top() ;
        stak.pop();
        if ( stak.size() )
        {
            os << ",";
        }
    }
    return os << "]" << endl;
}



// This debugger for printing 2d Array is written by Nafee
//typedef array< array<LL, SIZE> , SIZE > Matrix;
//ostream &operator << ( ostream & os, const Matrix &v ) {
//    os << endl;
//    FOR(a,0,v.size())
//    {
//        FOR(b,0,v[a].size())
//        {
//            os << v[a][b] << " ";
//        }
//        os << endl;
//    }
//    return os << endl;
//}







#define LL 			long long
#define PairLL		pair<long long, long long>
#define TripleLL    Triple<long long, long long, long long>
typedef vector<long long> VL;



#define RESET(a) 		memset(a,0,sizeof(a))
#define SET(a) 			memset(a,-1,sizeof(a))
#define SORT(v)         sort(v.begin(),v.end())
#define REV(v)          reverse(v.begin(),v.end())
#define MP              make_pair
#define PUB				push_back
#define POB				pop_back
#define PUF             push_front
#define POF             pop_front


#define FOR(i,a,b) for(long long i=(a);i<(b);i++)
#define ROF(i,a,b) for(long long i=(a);i>(b);i--)




#ifdef DEBUG
    #define dbg(x) cout<<#x<<" : "<<x<<endl
    #define dbg2(x, y) cout<<#x<<" : "<<(x)<<", "<<#y<<" : "<<(y)<<endl
    #define dbg3(x, y, z) cout<<#x<<" : "<<(x)<<", "<<#y<<" : "<<(y) << ", " << #z << " : " << (z) <<endl
    #define dbg4(a,b,c,d) cout<<#a<<" : "<<(a)<<", "<<#b<<" : "<<(b)<<", "<<#c<<" : "<<(c)<<", "<<#d<<": "<<(d)<<endl
    #define dbgAr(ar, si, ei) cout<<#ar<<" : "<<endl; prAr(ar, si, ei)
#else
    #define dbg(x)
    #define dbg2(x, y)
    #define dbg3(x, y, z)
    #define dbg4(a,b,c,d)
    #define dbgAr(ar, si, ei)
#endif // DEBUG




#define PI (2.0*acos(0.0)) //#define PI acos(-1.0)
#define SPRIME 10007
#define BPRIME 1000000007
#define LLMAX ( (unsigned long long) -1LL >> 1LL )





inline void iLL(LL &u)
{
    LL a, b, c, d;
    scanf("%lld", &u);
}

inline void iLL2(LL &u, LL &v)
{
    LL a, b, c, d, e;
    scanf("%lld %lld", &u, &v);
}

inline void iLL3(LL &u, LL &v, LL &w)
{
    LL a, b, c, d, e;
    scanf("%lld %lld %lld", &u, &v, &w);
}

inline void iLL4(LL &u, LL &v, LL &w, LL &x)
{
    LL a, b, c, d, e;
    scanf("%lld %lld %lld %lld", &u, &v, &w, &x);
}


string iStr(LL maxSize)
{
    int a, b, c, d, e;
    char str[maxSize];
    scanf("%s", str);
    long long len = strlen(str);
    string ret;
    for (a = 0; a < len; a++)
    {
        ret.push_back( str[a] );
    }
    return ret;
}



LL dirR[] = {1, -1, 0, 0};
LL dirC[] = {0, 0, 1, -1};


LL bigMod(LL a, LL b, LL M)
{
    if ( b == 0 ) return 1%M;
    LL x = bigMod(a, b/2, M);
    x = (x*x)%M;
    if ( b&1LL )
    {
        x = (x*a)%M;
    }
    return x;
}

LL getGcd(LL u, LL v)
{
    if ( u > v )
    {
        swap(u, v);
    }
    if ( u == 0 )
    {
        return v;
    }
    return getGcd(v%u, u);
}

LL getLcm(LL u, LL v)
{
    LL g = getGcd(u, v);
    return (u/g) * v;
}

template<typename T1, typename T2>
T1 aMin(T1 &u, T2 v)
{
    if ( v < u )
    {
        u = v;
    }
    return u;
}


template<typename T1, typename T2>
T1 aMax(T1 &u, T2 v)
{
    if ( v > u )
    {
        u = v;
    }
    return u;
}


template<typename T1>
vector<T1> getCumSumVec( vector<T1> inpVec )
{
    long long a, b, c, d, e;
    long long len = inpVec.size();
    vector<T1> cumSumVec(len, 0);
    if (len > 0)
    {
        cumSumVec[0] = inpVec[0];
    }
    for(a = 1; a <= len; a++)
    {
        cumSumVec[a] = cumSumVec[a-1] + inpVec[a];
    }
    return cumSumVec;
}




struct DSU
{
    LL *parentAr, *countAr;
    LL siz;
    DSU(LL siz)
    {
        LL a, b, c,d ;
        this->siz = siz;
        parentAr = new LL[siz];
        countAr = new LL[siz];
        FOR(a,0,siz)
        {
            parentAr[a] = a;
            countAr[a] = 1;
        }
    }
    LL getParent( LL u )
    {
        LL ret;
        if ( parentAr[u] != u )
        {
            parentAr[u] = getParent( parentAr[u] );
        }
        return parentAr[u];
    }
    void makePair(LL u, LL v)
    {
        if ( getParent(u) == getParent(v) )
        {
            return ;
        }
        countAr[ getParent(v) ] += countAr[ getParent(u) ];
        parentAr[ getParent(u) ] = getParent( v );
    }
    LL getCount( LL u )
    {
        return countAr[ getParent(u) ];
    }
};


struct BIT
{
    LL *treeAr;
    LL siz;
    BIT(){}
    BIT(LL siz)
    {
        this->siz = siz;
        treeAr = new long long[siz+9];
        RESET( treeAr);
    }
    long long read(long long idx)
    {
        long long sum = 0;
        while(idx > 0)
        {
            sum += treeAr[ idx ];
            idx -= (idx & -idx);
        }
        return sum;
    }
    void update(long long idx, long long val)
    {
        while(idx<=siz)
        {
            treeAr[idx] += val;
            idx += (idx & -idx);
        }
    }
};


// __builtin_popcountll = long long
// void * memcpy ( void * destination, const void * source, size_t num );

LL checkBit(LL num, LL pos)
{
    return (num>>pos) & 1LL;
}

LL setBit(LL num , LL pos)
{
    return num | (1LL<<pos);
}

LL resetBit(LL num, LL pos)
{
    return num & ( ~( 1LL<<pos ) );
}








/******   END OF HEADER *********/
#define SIZE 4
#define eps 1e-6


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
        Vector copyVec(this->x, this->y, this->z);
        Vector unitVector = copyVec * ( 1 / disFromMain );

        return unitVector;
    }

    void rotateAroundVector( const Vector &vec, double radAngle )
    {
        assert(false);

        if ( isParallel(vec) )
        {
            return;
        }

//        cout << vec.x << " " << vec.y << " " << vec.z << endl;

        Vector unitVec = this->getNormalizedVec();


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




struct Point
{
    double ar[SIZE];

    Point()
    {
        FOR(a,0,SIZE-1)
        {
            ar[a] = 0;
        }
        ar[SIZE-1] = 1;
    }

    double getX() const
    {
        return ar[0];
    }

    double getY() const
    {
        return ar[1];
    }

    double getZ() const
    {
        return ar[2];
    }

    double setX( double x )
    {
        ar[0] = x;
    }

    double setY(double y)
    {
        ar[1] = y;
    }

    double setZ( double z )
    {
        ar[2] = z;
    }

    Point(double x, double y, double z)
    {
        ar[0] = x;
        ar[1] = y;
        ar[2] = z;

        ar[3] = 1;
    }

    Vector operator - (const Point &B) const
    {
        double x = this->getX() - B.getX();
        double y = this->getY() - B.getY();
        double z = this->getZ() - B.getZ();

        Vector ret(x,y,z);
        return ret;
    }
};


Point takePointInput( ifstream &in )
{
    double x, y, z;
    in >> x >> y >> z;
    Point ret(x,y,z);

    return ret;
}


struct gluLookAtParam
{

};



struct Matrix
{
    double ar[SIZE][SIZE];

    Matrix()
    {
        FOR(a,0,SIZE)
        {
            FOR(b,0,SIZE)
            {
                ar[a][b] = 0;
            }
        }
    }

    double getVal(LL row, LL col) const
    {
        assert( row >= 0 && row < SIZE );
        assert( col >= 0 && col < SIZE );

        return ar[row][col];
    }

    Matrix operator * ( const Matrix &B ) const
    {
        Matrix ret;
        FOR(a,0,SIZE)
        {
            FOR(b,0,SIZE)
            {
                ret.ar[a][b] = 0;
                FOR(c,0,SIZE)
                {
                    ret.ar[a][b] += (this->getVal(a,c) * B.getVal(c,b) );
                }
            }
        }

        return ret;
    }

    Point operator * ( const Point &P ) const
    {

        Point ret;
        FOR(a,0,SIZE)
        {
            ret.ar[a] = 0;
            FOR(b,0,SIZE)
            {
                ret.ar[a] += (this->ar)[a][b] * P.ar[b];
            }
        }

        return ret;
    }
};


stack<Matrix> matStak;


int main()
{
    ifstream fin("scene.txt");

    ofstream stage1("stage1.txt");
    ofstream stage2("stage2.txt");
    ofstream stage3("stage3.txt");




    fin.close();

    stage1.close();
    stage2.close();
    stage3.close();
}


