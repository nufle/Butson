/* matrice.h

Simple matrix class
Template over dimensions to obtain fast creation/computations

Arnaud Mallen mallen.arnaud@gmail.com
Vincent Popie vincent.popie@gmail.com

 */

#ifndef _CLASSE_Matrice_
#define _CLASSE_Matrice_

#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>

#define PI 3.1415926535897932384626433832795
using std::cout;
using std::endl;
using std::ostringstream;
using std::ofstream;
using std::ifstream;
using std::cin;
using std::vector;
using std::string;

/* Templated matrix class
 */
template <class T,int l,int c>
class Matrix
{
    private:
    	T * mat;
    	int curRow;

    public:

    Matrix();
    Matrix(const Matrix<T,l,c> & m);
    Matrix(T * val);
    Matrix(vector<T> vals);
    Matrix(const int& size);
	~Matrix();


	T operator()(const int& i, const int& j) const;
	T& operator()(const int& i, const int& j);
	T getVal(int x,int y) const ;
    void setVal(int x, int y, T val);
    int getNbRows() const;
    int getNbCols()const;
    bool operator==(const Matrix & matrice) const;
    Matrix& operator=(const Matrix &rhs);


	void addZeroValuedRow();
	void orderColumnwise(vector< vector < T > > & M);

    // minimality check function and utility functions
    bool minimal();

	void disp();
    void write2ASCII (string filename);
    void readFromASCII (string filename);

private:
	void rowPermutation(int row1,int row2);
	vector <vector <T> > matrix2VecVec(const Matrix& M);
	bool isOrthogonal(int rootPower);
	double scalarProduct(int ligne1, int ligne2, int rootPower);


	
	
};

    template <typename T,int l,int c>
    Matrix<T,l,c>::Matrix(){
    	mat = (T*)malloc(l*c*sizeof(T));
    	curRow=0;
    }

    template <typename T,int l,int c>
    Matrix<T,l,c>::Matrix(const int& size){
    	mat = (T*)malloc(l*c*sizeof(T));
    	curRow=size;
    	for (int i=0;i<l*c;i++)
    		mat[i]=0;
    }

    template <typename T,int l,int c>
    Matrix<T,l,c>::Matrix(T*vals){
    	mat = (T*)malloc(l*c*sizeof(T));
    	curRow=l;
    	for (int i=0;i<l*c;i++)
    		mat[i]=vals[i];
    }

    template <typename T,int l,int c>
    Matrix<T,l,c>::Matrix(vector<T>vals){
    	mat = (T*)malloc(l*c*sizeof(T));
    	curRow=l;
    	for (int i=0;i<l*c;i++)
    		mat[i]=vals[i];
    }

    template <typename T,int l,int c>
    Matrix<T,l,c>::Matrix(const Matrix & m){
    	mat = (T*)malloc(l*c*sizeof(T));
    	curRow=m.curRow;
    	for (int i = 0 ; i < l*c;i++)
    		mat[i]=m.mat[i];
    }

    template <typename T,int l,int c>
    Matrix<T,l,c>::~Matrix(){
    	free(mat);
    }

    template <typename T,int l,int c>
    T Matrix<T,l,c>::operator()(const int& i, const int& j) const {
    	return mat[i*c+j];
    }

    template <typename T,int l,int c>
    T& Matrix<T,l,c>::operator()(const int& i, const int& j) {
    	return mat[i*c+j];
    }

    template <typename T,int l,int c>
    T Matrix<T,l,c>::getVal(int x, int y) const {
    		return mat[x*c+y];
    }

    template <typename T,int l,int c>
    int Matrix<T,l,c>::getNbRows() const{
    	return l;
    }

    template <typename T,int l,int c>
    int Matrix<T,l,c>::getNbCols()const{
    	return c;
    }

    template <typename T,int l,int c>
    void Matrix<T,l,c>::setVal(int x,int y, T val) {
    		mat[x*c+y] = val;
    }

    template <typename T,int l,int c>
    Matrix<T,l,c>& Matrix<T,l,c>::operator=(const Matrix &rhs) {
    	// Do not check self-assignment on purpose : faster?
    	curRow=rhs.curRow;
    	for (int i = 0; i< l*c;i++)
    	{
    		mat[i]=rhs.mat[i];
    	}
        return *this;
    }

    template <typename T,int l,int c>
    bool Matrix<T,l,c>::operator==(const Matrix &matrice) const
    {
		bool equal = true;
		for (int i = 0; equal && (i < l*c); i++)
		{
			if (mat[i]!=matrice.mat[i])
			{
				equal = false;
			}
		}
		return equal;
    }

    // Matrix to vector of vector of T. No return by ref optimization, no use.
    template <typename T,int l,int c>
    vector <vector <T> > Matrix<T,l,c>::matrix2VecVec(const Matrix& M){
    	vector< vector<T> > res(curRow,vector<T>(c));
    	for (int i=0; i<curRow; i++) {
    		for (int j=0; j<c; j++) {
    			res[i][j]=M.mat[i*c+j];
    		}
    	}
    	return res;
    }

    /*
     * Minimality check
     * 	- Converts input matrix to vector of vector of double
     * 	- Makes use of the lexicographic order natively implemented in STL
     */
    template <typename T,int l,int c>
    bool Matrix<T,l,c>::minimal(){

        bool result=true;
    	int currentRow=0;
    	vector< vector<T> > matAsVecVec=matrix2VecVec(*this);

    	// Permutation test on each row
    	for(typename vector< vector<T> >::iterator iter_lignes=matAsVecVec.begin();iter_lignes!=matAsVecVec.end();iter_lignes++){

    		// Temporary mat, i-th row permuted with last row
    		vector <vector <T> > permutLign=matrix2VecVec(*this);
    		(*(permutLign.begin()+currentRow)).swap(*(permutLign.end()-1));

    		orderColumnwise(permutLign);

    		// mat is not minimal if temp is smaller
    		bool permutEstPlusPetite=*(permutLign.begin()+currentRow)<*(matAsVecVec.begin()+currentRow);
    		if (permutEstPlusPetite) {
    			result=false;
    			break;
    		}
    		else{
    			// If not, check that i-th row bigger than last row
    			if (*(permutLign.begin()+currentRow)==*(matAsVecVec.begin()+currentRow)){
    				if ((matAsVecVec[currentRow])>(matAsVecVec[curRow-1])){
    					result=false;
    					break;
    				}
    			}
    		}
    		currentRow++;
    	}
    	return result;
    }

    // Transposee d'un vecteur de vecteur d'T
    template <typename T,int l,int c>
    vector < vector<T> > transposeVecVec(const vector< vector < T > > & M){
    	int sizeM=M[1].size();
    	vector<T> tmp(M.size(),0);
    	vector< vector<T> > mprime(sizeM,tmp);
     	//#pragma omp parallel for
        for(unsigned int i=0;i<M[1].size();i++)
            for(unsigned int j=0;j<M.size();j++)
                mprime[i][j]=M[j][i];
        return mprime;
    }


    // Orders using lexicographic order, using STL sort
    template <typename T,int l,int c>
    void Matrix<T,l,c>::orderColumnwise(vector< vector < T > > & M){
        vector< vector<T> > tmp=transposeVecVec<T,l,c>(M);
        sort(tmp.begin(),tmp.end());
        M=transposeVecVec<T,l,c>(tmp);
    }

    // Permutation des lignes d'un vecteur
    template <typename T,int l,int c>
    void Matrix<T,l,c>::rowPermutation(int row1, int row2){
    	vector<T> row1copy(c,0);
    	for(int i=0;i<c;i++)
    		row1copy[i]=mat[c*row1+i];
    	for (int i=0; i<c; i++) {
    		this->mat[row1*c+i]=this->mat[row2*c+i];
    		this->mat[row2*c+i]=row1copy[i];
    	}
    }

    // Ajouter une ligne de zero a une matrice
    // Cree un vecteur de zeros, l'ajoute a la fin du vecteur, incremente le nombre de ligne
    template <typename T,int l,int c>
    void Matrix<T,l,c>::addZeroValuedRow(){
    	curRow++;
    }

    // Produit scalaire de deux vecteurs de nombres complexes. On rappele que pour z1 et z2, <z1,z2>=a1a2+b1b2.
    template <typename T,int l,int c>
    double Matrix<T,l,c>::scalarProduct(int ligne1, int ligne2,int rootPower){
    	double resultat=0.0;
    	for (int i=0;i<c;i++)
    		resultat+=cos(2*mat[ligne1*c+i]*PI/rootPower)*cos(2*mat[ligne2*c+i]*PI/rootPower)+sin(2*mat[ligne1*c+i]*PI/rootPower)*sin(2*mat[ligne2*c+i]*PI/rootPower);
    	return resultat;
    }

    // Check orthogonality of a matrix
    template <typename T,int l,int c>
    bool Matrix<T,l,c>::isOrthogonal(int rootPower){
    	for (int i=0;i<l;i++){
			for (int j=i+1; j<l; j++) {
				if (!(Matrix<T,l,c>::scalarProduct(i,j,rootPower)<10e-15)){
					return false;
				}
			}
    	}
    	return true;
    }

    // Basic Cout display
    template <typename T,int l,int c>
    void Matrix<T,l,c>::disp(){
        for (int ligne=0;ligne<curRow;ligne++){
            cout << "|";
            for (int colonne=0;colonne<c;colonne++){
                cout <<" "<<mat[ligne*c+colonne];
            }
            cout << "|"<<endl;
        }
    	cout<<endl;
    }


    // ASCII file writing
    template <typename T,int l,int c>
    void Matrix<T,l,c>::write2ASCII (string filename)
    {
        ofstream outfile;

        outfile.open(filename.c_str());

        if(outfile)
        {
            outfile << l << " " << c << " " << endl;
            for (int i=0; i < l; i++)
            {
                for (int j=0; j<c; j++)
                    outfile << getVal(i,j) << " ";

                outfile << endl;
            }
        }
        else
            cout << "ERR : "<<filename<<" cannot be open (WRITE)"<<std::endl;

        outfile.close();
    }

    // Lecture dans un fichier d'une matrice
    template <typename T,int l,int c>
    void Matrix<T,l,c>::readFromASCII (string filename)
    {
        ifstream infile;
        int k,m;
        infile.open(filename.c_str());

        if(infile)
        {
            infile >> k >> m;

            for (int i=0; i < l; i++)
            {
                for (int j=0; j<c; j++)
                {
                    infile >> k;
                    setVal(i,j,k);
                }
            }
        }
        else
        	cout << "ERR : "<<filename<<" cannot be open (READ)"<<std::endl;

        infile.close();

    }


#endif
