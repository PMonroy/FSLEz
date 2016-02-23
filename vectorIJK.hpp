#ifndef VECTORIJK
#define VECTORIJK

#include <iostream>

using namespace std;

class  vectorIJK {// Class int 3D vector

public:
  int i,j,k; // Order list of 3 elements |i|j|k|

  vectorIJK(void);                 //Void Constructor
  vectorIJK(int ii,int jj,int kk); //Constructor
  vectorIJK(const vectorIJK &a);   //Copy vector constructor
  void operator=(vectorIJK a);     //Copy vector

  vectorIJK &operator+=(vectorIJK a); //Addition increment
  vectorIJK &operator-=(vectorIJK a); //Subtraction increment
  vectorIJK &operator*=(vectorIJK a); //Product increment
  vectorIJK &operator/=(vectorIJK a); //Division increment
  vectorIJK &operator*=(int scalar);  //Scalar product increment

  vectorIJK operator+( vectorIJK a);  // Addition
  vectorIJK operator-( vectorIJK a);  // Subtraction
  vectorIJK operator*( vectorIJK a);  // Product
  vectorIJK operator/( vectorIJK a);  // Division

  friend ostream &operator<<(ostream &out, vectorIJK a); //Format output
  friend istream &operator>>(istream &in, vectorIJK &a); //Format input

  ~vectorIJK();	                            // Destructor
};

vectorIJK  operator*(int scalar, vectorIJK a); // Scalar left product
vectorIJK  operator*(vectorIJK a,int scalar);  // Scalar right product

bool operator==(vectorIJK a,vectorIJK b); //Boolean operator ==

#endif 
