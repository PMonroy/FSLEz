#include "vectorIJK.hpp"

vectorIJK::vectorIJK(void) {//Void constructor
	i=0;
	j=0;
	k=0;
}
vectorIJK::vectorIJK(int ii,int jj,int kk) {//Constructor
	i=ii;
	j=jj;
	k=kk;
}
vectorIJK::vectorIJK(const vectorIJK &a) {//Copy vector constructor
	i=a.i;
	j=a.j;
	k=a.k;
}
void vectorIJK::operator=(vectorIJK a) {//Copy vector
	i=a.i;
	j=a.j;
	k=a.k;
}

vectorIJK &vectorIJK::operator+=(vectorIJK a) {//Addition increment
    i += a.i;
    j += a.j;
    k += a.k;
    return *this;
}
vectorIJK &vectorIJK::operator-=(vectorIJK a) {//Subtraction increment
    i -= a.i;
    j -= a.j;
    k -= a.k;
    return *this;
}
vectorIJK &vectorIJK::operator*=(vectorIJK a) {//Product increment
    i *= a.i;
    j *= a.j;
    k *= a.k;
    return *this;
}
vectorIJK &vectorIJK::operator/=(vectorIJK a) {//Division increment
    i /= a.i;
    j /= a.j;
    k /= a.k;
    return *this;
}
vectorIJK &vectorIJK::operator*=(int scalar) {//Scalar product increment  
    i *= scalar;
    j *= scalar;
    k *= scalar;
    return *this;
}

vectorIJK  vectorIJK::operator+(vectorIJK a) {//Addition
  return vectorIJK(*this)+=a;
}
vectorIJK  vectorIJK::operator-(vectorIJK a) {//Subtraction
  return vectorIJK(*this)-=a;
}
vectorIJK  vectorIJK::operator*(vectorIJK a) {//Product
  return vectorIJK(*this)*=a;
}
vectorIJK  vectorIJK::operator/(vectorIJK a) {//Division
  return vectorIJK(*this)/=a;
}
vectorIJK operator*( int scalar,vectorIJK a) {//Scalar left product
  return vectorIJK(a) *= scalar;
} 
vectorIJK operator*(vectorIJK a, int scalar) {//Scalar right product
  return vectorIJK(a) *= scalar;
} 

vectorIJK::~vectorIJK() {//Destructor 

}    

ostream &operator<<(ostream &out, vectorIJK a) { //Format output
        out<<a.i<<" "<<a.j<<" "<<a.k;
        return out;
}
istream &operator>>(istream &in, vectorIJK &a) {//Format input
  in>>a.i;
  in>>a.j;
  in>>a.k;  
  return in;
}

bool operator==(vectorIJK a,vectorIJK b) {//Boolean operator ==
  return (a.i==b.i && a.j==b.j && a.k==b.k);
}
