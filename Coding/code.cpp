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


struct PointOrVector
{
    int X_POS = 0;
    int Y_POS = 1;
    int Z_POS = 2;
    int W_POS = 3;
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

    double getW() const
    {
        return ar[W_POS];
    }


    void setX( double x )
    {
        ar[X_POS] = x;
    }

    void setY(double y)
    {
        ar[Y_POS] = y;
    }

    void setZ( double z )
    {
        ar[Z_POS] = z;
    }

    void setW( double w )
    {
        ar[W_POS] = w;
    }


    PointOrVector(double x, double y, double z)
    {
        setX(x);
        setY(y);
        setZ(z);
    }

    PointOrVector(double x, double y, double z, double w)
    {
        setX(x);
        setY(y);
        setZ(z);
        setW(w);
    }


    PointOrVector operator - (const PointOrVector &B) const
    {
        PointOrVector ret;
        for (int i = 0; i < SIZE; i++)
        {
            ret.ar[i] = this->ar[i] - B.ar[i];
        }
        return ret;
    }


    PointOrVector operator+(const PointOrVector &p) const
    {
        PointOrVector ret;
        for (int i = 0; i < SIZE; i++)
        {
            ret.ar[i] = this->ar[i] + p.ar[i];
        }
        return ret;
    }

    PointOrVector operator*(double i) const
    {
        return PointOrVector( getX()*i, getY()*i, getZ()*i );
    }

    PointOrVector operator+=(const PointOrVector &p)
    {
        for (int i = 0;  i < SIZE; i++)
        {
            this->ar[i] += p.ar[i];
        }

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





const PointOrVector unitXVec(1,0,0);
const PointOrVector unitYVec(0,1,0);
const PointOrVector unitZVec(0,0,1);






PointOrVector takePointInput( ifstream &in )
{
    double x, y, z;
    in >> x >> y >> z;
    PointOrVector ret(x,y,z);

    return ret;
}


struct Triangle
{
    PointOrVector ar[3];

    Triangle()
    { }

};


Triangle takeTriangleInput( ifstream in )
{
    Triangle ret;
    for (int a = 0; a < 3; a++)
    {
        ret.ar[a] = takePointInput(in);
    }
    return ret;
}

struct GluLookAtParam
{

    PointOrVector eyePosition;
    PointOrVector lookPosition;
    PointOrVector upDirection;

    double perspectiveAr[4];

    GluLookAtParam( ifstream &in )
    {
        eyePosition = takePointInput(in);
        lookPosition = takePointInput(in);
        upDirection = takePointInput(in);

        FOR(a,0,4)
        {
            in >> perspectiveAr[a];
        }
    }


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

    PointOrVector operator * ( const PointOrVector &P ) const
    {

        PointOrVector ret;
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

    GluLookAtParam gluLookAtParam(fin);

    string command;
    fin >> command;

    cout << command;


    fin.close();

    stage1.close();
    stage2.close();
    stage3.close();
}


