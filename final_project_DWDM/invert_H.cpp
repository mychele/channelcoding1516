#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/GF2X.h>
#include <NTL/GF2.h>
#include <NTL/mat_GF2.h>
#include <iostream>
#include <fstream>


using namespace NTL;

void printMatrix(mat_GF2& K, const char *matrixName)
{
    std::ofstream output_file (matrixName, std::ios::out);
    if(output_file.is_open()) {
	    for(int row_index = 0; row_index < K.NumRows(); row_index++) {
	    	for(int col_index = 0; col_index < K.NumCols(); col_index++) {
	    		output_file << K[row_index][col_index] << " "; 
	    	}
	    	output_file << "\n";
	    }
	}
	output_file.close();
}


int main(int argc, char const *argv[])
{
	// define H
	int info_r = 105;
	int red_r = 7;
	int col_r = 293;
	int check_bit = red_r*col_r - 6;
	mat_GF2 H;
	H.SetDims(red_r*col_r, (info_r + red_r)*col_r);
	std::cout << "H is " << H.NumRows() << " x " << H.NumCols() << "\n";

	// define the slopes
	int slopes[7] = {1,2,3,4,5,6,7};

	// fill H
	int row_index = 0;
	// cycle on the slopes
	for(int slope_index = 0; slope_index < red_r; slope_index++) {
		// cycle on the offsets
		for(int c = 0; c < col_r; c++) {
			// cycle on a from 0 to 111
			for(int a = 0; a < (info_r + red_r); a++) {
				int j = a*col_r + (a*slopes[slope_index] + c)%col_r;
				H[row_index][j] = 1;
			}
			row_index++;
		}
	}  
	std::cout << "H filled\n";

	// this check controls if each row has 112 ones, and each column has 7 ones. It is disabled by default
	bool checkNumOnes = 0;
	if (checkNumOnes) {
		// as a check, count the number of 1 in the rows
		for (int row_index = 0; row_index < H.NumRows(); row_index++) {
			int numOnes = 0;
			for(int col_index = 0; col_index < H.NumCols(); col_index ++) {
				if (H[row_index][col_index] == 1) {numOnes++;}
			}
			if( numOnes != 112 ) {std::cout << "Row " << row_index << " has " << numOnes << " ones\n";}
		}

		// as a check, count the number of 1 in the columns
		for(int col_index = 0; col_index < H.NumCols(); col_index ++) {
			int numOnes = 0;
			for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				if (H[row_index][col_index] == 1) {numOnes++;}
			}
			if( numOnes != 7 ) {std::cout << "Column " << col_index << " has " << numOnes << " ones\n";}
		}
	}

	// copy the matrix into a matrix without six columns
	mat_GF2 H_fr;
	H_fr.SetDims(red_r*col_r, (info_r + red_r)*col_r - 6);

	// copy the first 105*293+291 columns
	for(int col_index = 0; col_index < (105*col_r+col_r-1); col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index] = H[row_index][col_index];
		}
	}
	// skip column 105*293+292, copy 292 columns
	// row 106 in data_matrix
	int offset_H = (105*col_r+col_r-1) + 1;
	int offset_H_fr = (105*col_r+col_r-1);
	for(int col_index = 0; col_index < col_r - 1; col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index + offset_H_fr] = H[row_index][col_index + offset_H];
		}
	}
	// skip column 106*293+292, copy 292 columns
	// row 107 in data_matrix
	offset_H = (106*col_r+col_r-1) + 1;
	offset_H_fr = (106*col_r+col_r-1) - 1; // 2 columns were deleted
	for(int col_index = 0; col_index < col_r - 1; col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index + offset_H_fr] = H[row_index][col_index + offset_H];
		}
	}
	// skip column 107*293+292, copy 292 columns
	// row 108 in data_matrix
	offset_H = (107*col_r+col_r-1) + 1;
	offset_H_fr = (107*col_r+col_r-1) - 2; // 3 columns were deleted
	for(int col_index = 0; col_index < col_r - 1; col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index + offset_H_fr] = H[row_index][col_index + offset_H];
		}
	}
	// skip column 108*293+292, copy 292 columns
	// row 109 in data_matrix
	offset_H = (108*col_r+col_r-1) + 1;
	offset_H_fr = (108*col_r+col_r-1) - 3; // 4 columns were deleted
	for(int col_index = 0; col_index < col_r - 1; col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index + offset_H_fr] = H[row_index][col_index + offset_H];
		}
	}
	// skip column 109*293+292, copy 292 columns
	// row 110 in data_matrix
	offset_H = (109*col_r+col_r-1) + 1;
	offset_H_fr = (109*col_r+col_r-1) - 4; // 5 columns were deleted
	for(int col_index = 0; col_index < col_r - 1; col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index + offset_H_fr] = H[row_index][col_index + offset_H];
		}
	}
	// skip column 110*293+292, copy 293 columns
	// row 111
	offset_H = (110*col_r+col_r-1) + 1;
	offset_H_fr = (110*col_r+col_r-1) - 5; // 6 columns were deleted
	for(int col_index = 0; col_index < col_r; col_index++) {
		for(int row_index = 0; row_index < H.NumRows(); row_index++) {
				H_fr[row_index][col_index + offset_H_fr] = H[row_index][col_index + offset_H];
		}
	}

	mat_GF2 H_fr_2;
	H_fr_2.SetDims(H_fr.NumRows() - 6, H_fr.NumCols());
	// copy the first 292 rows (from 0 to 291)
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index][col_index] = H_fr[row_index][col_index];
		}
	}
	// skip row 292, copy 292 rows
	offset_H_fr = col_r;
	int offset_H_fr_2 = offset_H_fr - 1; // 1 row was deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		}
	}
	// skip row 585, copy 292 rows
	offset_H_fr = 2*col_r;
	offset_H_fr_2 = offset_H_fr - 2; // 2 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		}
	}
	// skip row 878, copy 292 rows
	offset_H_fr = 3*col_r;
	offset_H_fr_2 = offset_H_fr - 3; // 3 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		}
	}
	// skip row 1171, copy 292 rows
	offset_H_fr = 4*col_r;
	offset_H_fr_2 = offset_H_fr - 4; // 4 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		}
	}
	// skip row 1464, copy 292 rows
	offset_H_fr = 5*col_r;
	offset_H_fr_2 = offset_H_fr - 5; // 5 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		}
	}
	// skip row 1757, copy 293 rows
	offset_H_fr = 6*col_r;
	offset_H_fr_2 = offset_H_fr - 6; // 5 rows were deleted
	for (int row_index = 0; row_index < col_r; row_index++) {
		for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
			H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		}
	}
	// remove useless H_fr_2 and keep just H_fr
	H_fr.kill();
	H_fr = H_fr_2;
	H_fr_2.kill();

	std::cout << "H_fr filled\n";
	std::cout << "H_fr is " << H_fr.NumRows() << " x " << H_fr.NumCols() << "\n";
	std::cout << "H_fr = [M | N]\n";

	mat_GF2 N;
	N.SetDims(check_bit, check_bit);
				//	293 x 105
	int N_offest = col_r*info_r;
	for(int row_index = 0; row_index < H_fr.NumRows(); row_index++) {
		for(int col_index =0; col_index < N.NumCols(); col_index++) {
			N[row_index][col_index] = H_fr[row_index][col_index + N_offest];
		}
	}

	mat_GF2 M;
	M.SetDims(check_bit, N_offest);
				//	293 x 105
	for(int row_index = 0; row_index < H_fr.NumRows(); row_index++) {
		for(int col_index = 0; col_index < M.NumCols(); col_index++) {
			M[row_index][col_index] = H_fr[row_index][col_index];
		}
	}

	std::cout << "M and N filled\n";

	// we need to permute the columns of H_fr so that the last 2045 become the first 2045
	mat_GF2 permutation_matrix;
	permutation_matrix.SetDims(H_fr.NumCols(), H_fr.NumCols());
	long perm_offset = col_r*info_r;
					//	293  * 105					32810
	for(int row_index = perm_offset; row_index < H_fr.NumCols(); row_index++) {
		permutation_matrix[row_index][row_index - perm_offset] = 1;
	}
	perm_offset = 2045;
	for(int col_index = perm_offset; col_index < H_fr.NumCols(); col_index++) {
		permutation_matrix[col_index - perm_offset][col_index] = 1;
	}
	std::cout << "permutation_matrix built\n";

	mat_GF2 H_toinv = H_fr*permutation_matrix;

	// append identity matrix 
	mat_GF2 H_eye;
	H_eye.SetDims(H_toinv.NumRows(), H_toinv.NumCols() + H_toinv.NumRows());
	// copy H_toinv into H_eye
	for(int col_index = 0; col_index < H_toinv.NumCols(); col_index++) {
		for(int row_index = 0; row_index < H_toinv.NumRows(); row_index++) {
				H_eye[row_index][col_index] = H_toinv[row_index][col_index];
		}
	}
	// create eye(H_toinv.NumRows())
	int offset = H_toinv.NumCols();
	for(int col_index = H_toinv.NumCols(); col_index < H_eye.NumCols(); col_index++) {
		H_eye[col_index - offset][col_index] = 1;
	}

	std::cout << "H_eye filled\n";

	gauss(H_eye);

	std::cout << "Gaussian elimination performed\n";

	//index of last nonzero row
    int last_nonzero = H_eye.NumRows() - 1;
    bool flag_nonzero = 0;
    do
    {
        flag_nonzero=0;
        for (int col_index = 0; col_index < H_fr.NumCols(); col_index++) {
            if (H_eye[last_nonzero][col_index] != 0)
                flag_nonzero = 1;
        }
        if (flag_nonzero == 0) {last_nonzero--;} // check the row above this one
        else break; // exit the cycle since a non all zero row was found
    }
    while (last_nonzero>=0);

    std::cout << "All zero rows " << H_eye.NumRows() - 1 - last_nonzero << "\n";
    std::cout << "Non zero rows " << last_nonzero + 1 << "\n";

    bool find_pivot = 0;
    if (find_pivot) { // let's find the pivot column for each row, it should be in column=row
	    long pivot[H_eye.NumRows()];
	    // initialize the array
	    for (int ar_index = 0; ar_index < H_eye.NumRows(); ar_index++) {
	    	pivot[ar_index] = H_fr.NumCols();
	    }
	    for (int row_index = 0; row_index < H_eye.NumRows(); row_index++) {
	    	for (int col_index = 0; col_index < H_fr.NumCols(); col_index++) {
	    		if (H_eye[row_index][col_index] == 1) {
	    			pivot[row_index] = col_index;
	    			break;
	    		}
	    	}
	    }
	    for (int ar_index = 0; ar_index < H_eye.NumRows(); ar_index++) {
	    	if (ar_index != pivot[ar_index]) {
	    		std::cout << ar_index << " " << pivot[ar_index] << "\n";
	    	}
	    }
	}

	std::cout << "Performing Jordan\n";
    // Jordan procedure
    for (int pivot_index = last_nonzero; pivot_index >= 0; pivot_index--) {
    	// check if the (rows, pivot) above the pivot of this row are all 0. If not,
    	// sum row=pivot to the interested row
    	for (int row_index = pivot_index - 1; row_index >= 0; row_index--) {
    		if(H_eye[row_index][pivot_index] == 1) { // sum
    			for (int col_index = pivot_index; col_index < H_eye.NumCols(); col_index++) { // just from pivot_index since it is 
    				H_eye[row_index][col_index] += H_eye[pivot_index][col_index];			  // the first non zero entry of H[pivot_index][:]
    			}
    		}
    	}
    }

    std::cout << "Jordan performed\n";

    // check if there is an identity matrix on the left
    mat_GF2 maybeEye;
    maybeEye.SetDims(2045, 2045);
    for (int row_index = 0; row_index < maybeEye.NumRows(); row_index++) {
    	for (int col_index = 0; col_index < maybeEye.NumCols(); col_index++) {
    		maybeEye[row_index][col_index] = H_eye[row_index][col_index];
    	}
    }

    std::cout << "Is the submatrix on the left (2045x2045) an identity? " << ((IsIdent(maybeEye, maybeEye.NumCols())==1)?"yes\n":"no\n");
    
    // H_eye is [ I | 	K 	| N_1 ]
    //			2045 105*293 2045

    // let's isolate N^{-1}
    mat_GF2 N_1;
    N_1.SetDims(N.NumRows(), N.NumCols());
	offset = H_toinv.NumCols();
	for(int col_index = offset; col_index < H_eye.NumCols(); col_index++) {
		for(int row_index = 0; row_index < H_eye.NumRows(); row_index++) {
			N_1[row_index][col_index - offset] = H_eye[row_index][col_index];
		}
	}
	std::cout << "N_1 created, is N_1*N an identity? " << ((IsIdent(N_1*N, N.NumCols())==1)?"yes\n":"no\n");

	// let's create K = N_1*M
	mat_GF2 K;
	K.SetDims(H_eye.NumRows(), col_r*info_r);
	offset = N.NumRows();
	for(int col_index = offset; col_index < H_toinv.NumCols(); col_index++) {
		for(int row_index = 0; row_index < K.NumRows(); row_index++) {
			K[row_index][col_index-offset] = H_eye[row_index][col_index];
		}
	}

	// check if K = N_1*M
	std::cout << "K created, is K=N_1*M? " << (IsZero(N_1*M + K) ? "yes\n":"no\n");
	printMatrix(K, "K.txt");


	// create Hprime
	mat_GF2 H_prime;
	H_prime.SetDims(H_toinv.NumRows(), H_toinv.NumCols());
	std::cout << "H prime is " << H_prime.NumRows() << " x " << H_prime.NumCols() << "\n";
	for(int col_index = 0; col_index < K.NumCols(); col_index++) {
		for(int row_index = 0; row_index < K.NumRows(); row_index++) {
			H_prime[row_index][col_index] = K[row_index][col_index];
		}
	}
	offset = K.NumCols();
	std::cout << "K copied in H_prime, offset " << offset << "\n";
	for(int col_index = offset; col_index < H_prime.NumCols(); col_index++) {
		H_prime[col_index-offset][col_index] = 1;
	}

	std::cout << "H_prime create\n";

	H_toinv.kill();
	maybeEye.kill();

	// create G_prime
	mat_GF2 G_prime;
	G_prime.SetDims(H_prime.NumCols(), K.NumCols());
	// identity on top
	for(int row_index = 0; row_index < K.NumCols(); row_index++) {
		G_prime[row_index][row_index] = 1;
	}
	// K on bottom
	offset = K.NumCols();
	for(int row_index = offset; row_index < G_prime.NumRows(); row_index++) {
		for(int col_index = 0; col_index < G_prime.NumCols(); col_index++) {
			G_prime[row_index][col_index] = K[row_index - offset][col_index];
		}
	}

	// let's check...
	std::cout << "H_prime x G_prime = 0? " << (IsZero(H_prime*G_prime) ? "yes\n":"no\n");
	std::cout << "H_fr x G_prime = 0? " << (IsZero(H_fr*G_prime) ? "yes\n":"no\n");

	return 0;
}


