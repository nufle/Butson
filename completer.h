/* matrice.h

Completer class, tries to complete a butson matrix

Arnaud Mallen mallen.arnaud@gmail.com
Vincent Popie vincent.popie@gmail.com

 */

#ifndef COMPLETER_H_
#define COMPLETER_H_
#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>

#include "ObjectPool.h"
#include "matrix.h"

using std::cout;
using std::endl;
using std::ostringstream;


template <typename T, int n, int l>
class Completer
{
public:

	Completer():nbSol(0){};

	void complete(Matrix<T, l*n , l*n > *M, int i, int j, std::vector<T>* L);

	void checkSolutions();

private:

	// Matrix object pool
	ObjectPool< Matrix<int,l*n,l*n> > opMat;

	// Vector object pool
	ObjectPool< std::vector<T> > opVec;

	int nbSol;
};

template <int l>
static int diffModuloL(int p,int q){
	int r=p-q;
	return ((p-q) < 0) ? l+r : r;
}

inline static string int2str(int i){
	ostringstream oss;
	oss << i;
    return oss.str();
}

// Vérifie que les solutions sont toutes différentes
template <typename T, int n, int l>
void Completer<T,n,l>::checkSolutions(){
	// deux matrices temporaires contenant les solutions à comparer
	const int size = l*n;
	Matrix<int,size,size> SOL_I(n*l);
	Matrix<int,size,size> SOL_J(n*l);
	bool differentes=true;

	// Loop over the solutions
	for(int i=1;i<nbSol;i++){
		for (int j=i+1; j<nbSol; j++) {
			// Read solutions
			SOL_I.readFromASCII("butson_"+int2str(n)+"_"+int2str(l)+"_"+int2str(i));
			SOL_J.readFromASCII("butson_"+int2str(n)+"_"+int2str(l)+"_"+int2str(j));

			// Are solutions different?
			if (SOL_I==SOL_J) {
				cout<<endl<<"Solutions "<<i<<" and "<<j<<" are equal!";
				differentes=false;
			}

		}
	}
	if(differentes&&nbSol>0)
	{
		cout<<endl<<"All Solutions are different"<<std::endl;
	}
	else
	{
		cout<<"No solutions"<<std::endl;
	}
}

template <typename T, int n, int l>
void Completer<T,n,l>::complete(Matrix<T, l*n , l*n > *M, int i, int j, std::vector<T>* L)
{
	// Obtain
	Matrix<T,l*n,l*n> * m = opMat.New();
	*m=*M;
	std::vector<T> * K = opVec.New();

	if (j!=0){
		bool equal=true;
		for(int ii=0;ii<i;ii++){
			if (m->operator()(ii,j-1)!=m->operator()(ii,j)){
				equal=false;
				break;
			}
		}
		int minimum=0;
		if (equal){
			minimum=m->operator()(i,j-1);
		}
		for (T curValue=minimum;curValue<l;curValue++){
			*K=*L;
			bool ok=true;
			for (T k=0;k<i;k++){
				T ind=diffModuloL<l>(curValue,m->operator()(k,j));
				K->operator[](k*l+ind)++;
				if (K->operator[](k*l+ind)==(n)+1){
					ok=false;
				}
			}
			if (ok){
				m->operator()(i,j)=curValue;
				if ( !(i==(n)*(l)-1) || !(j==(n)*(l)-1)){
					if (!(j==(n)*(l)-1)) {
						Completer<T,n,l>::complete(m,i,(j+1),K);
					}else {
						Completer<T,n,l>::complete(m,(i+1),0,K);
					}
				}else {
					if (m->minimal())
					{
						m->disp();
						m->write2ASCII("butson_"+int2str(n)+"_"+int2str(l)+"_"+int2str(nbSol));
						nbSol++;
					}
				}
			}
		}
	}else {
		if(m->minimal()){
			m->addZeroValuedRow();
			std::vector<T> TMP(i*l,0);
			*L=TMP;
			for (int row=0;row<i;row++)
				L->operator[](row*j)=1;
			Completer<T,n,l>::complete(m,i,1,L);
		}
	}
	opMat.Delete(m);
	opVec.Delete(K);
}

#endif /* COMPLETER_H_ */
