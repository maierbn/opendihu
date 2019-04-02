#include <Python.h>

#include <lapacke.h>

#include <stdlib.h>
#include <stdio.h>
#include <string>

class DmdUtility
{
public:
	static int getSpatComp(double input[], int rows, int cols, double leftSingVec[], double sigma[], double rightSingVec[], double epsilon);

	static double getL2norm(double input[], int order, int range);

	static double getL2norm(double _Complex input[], int order, int range);

	static void getReducedSnapshotsMatrix(double sigma[], double rightSingVec[], double hatT[], int min, int cols, int spatComp);

	static void resizeMatrix(double input[], double output[], int oldRows, int newRows, int firstCol, int lastCol);

	static void resizeMatrix(double _Complex input[], double _Complex output[], int oldRows, int newRows, int firstCol, int lastCol);

	static void getMatrixMult(double inputA[], double inputB[], double output[], int rowsA, int colsA_rowsB, int colsB);

	static void getMatrixMult(double _Complex inputA[], double _Complex inputB[], double _Complex output[], int rowsA, int colsA_rowsB, int colsB);

	static void getMatrixMult(double inputA[], double _Complex inputB[], double _Complex output[], int rowsA, int colsA_rowsB, int colsB);

	static void getMatrixMult(double _Complex inputA[], double _Complex inputB[], double output[], int rowsA, int colsA_rowsB, int colsB);

	static void printMatrix(std::string name, double input[], int rows, int cols);

	static void printMatrix(std::string name, double _Complex input[], int rows, int cols);

	static void getKoopmanOperator(double input[], double output[], int rows, int cols);

	static void transposeMatrix(double input[], double output[], int rows, int cols);

	static void transposeMatrix(double _Complex input[], double _Complex output[], int rows, int cols);

	static void getMatrixInverse(double a[], int order);

	static void getEigen(double input[], int order, double _Complex eigenvalues[], double _Complex eigenvectors[]);

	static void getDeltaOmega(double _Complex eigenvalues[], double growthRates[], double frequencies[], int size, double deltat);

	static void getAmplitudes(double snapshots[], double _Complex eigenvalues[], double _Complex eigenvectors[], int rows, int cols, double _Complex amplitudes[]);

	static void setZero(double input[], int rows, int cols);

	static void setZero(double _Complex input[], int rows, int cols);

	static void getMatrixPower(double _Complex input[], double _Complex output[], int order, int exponent);

	static void getMatrixLeftDivision(double inputA[], double _Complex inputB_output[], int n, int nrhs);

	static void getDmdModes(double _Complex dmdModes[], double _Complex a[], double _Complex eigenvectors[], int rows, int cols, double leftSingVec[], double leftSingVecReduced[], double amplitudes[]);

	static int getSpecComp(double _Complex dmdModes[], double growthRates[], double frequencies[], double amplitudes[], int order, double epsilon0);

	static void concatenateVector(double _Complex inputA[], double inputB[], double _Complex output[], int rows, int cols);

	static void sortMatrix(double _Complex input[], double _Complex output[], int rows, int cols);

	static void getDmdModesGrowthRatesFrequencies(double _Complex dmdModes[], double growthRates[], double frequencies[], double _Complex dmdModesReduced[], double growthRatesReduced[], double frequenciesReduced[], int rows, int cols, int compSpec);

	static void reconstructSnapshots(double _Complex dmdModes[], double growthRates[], double frequencies[], double leftSingVec[], double snapshotsReconst[], int rows, int cols, int specComp, int spatComp, double deltat);

	static void contReconst(double t, double t0, double _Complex dmdModes[], double growthRates[], double frequencies[], int rows, int cols, double output[]);

	static int getCSVRowCount(std::string filename);
};