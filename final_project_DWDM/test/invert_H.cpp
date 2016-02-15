#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/GF2X.h>
#include <NTL/GF2.h>
#include <NTL/mat_GF2.h>
#include <iostream>
#include <fstream>
#include <random>
#include <bitset>
#include <array>
#include "bit_io.h"


using namespace NTL;

#define INFO_ROWS 105
#define PC_ROWS 7
#define ALL_COLUMNS 293
#define INFO_BIT 30592
#define INIT_ZERO_BIT 173
#define ALL_INFO_BIT 30765
#define IND_EQ 2045

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

int mapJtoK(int j) {
	int init_zero_bit = INIT_ZERO_BIT;
   	int col_r = ALL_COLUMNS;
	int q = j + init_zero_bit; // j + 172, with j from 1 to 30592
	int r = q/col_r; // floor(q/293) since both q and 293 are positive int, the result is the floor
	int col_index = col_r*r + col_r - 1 - q; // 293*r + 292 - q
	return col_r*r + col_index;
}


int main(int argc, char const *argv[])
{
	// rng stuff
	std::random_device rd; // device entropy
    std::mt19937 m_rng(rd()); // initialize our mersenne twister with a random seed
    std::uniform_int_distribution<int> int_uni_gen(0,1);

    const int info_r = INFO_ROWS;
    const int red_r = PC_ROWS;
    const int col_r = ALL_COLUMNS;
    const int info_bit = INFO_BIT;
    const int init_zero_bit = INIT_ZERO_BIT;

	// define H
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

	// save H to file
	std::ofstream bin_H_out("H.bin", std::ios::out | std::ios::binary);
	BitIo<(info_r + red_r)*col_r> bio_H;
	for(int row_index = 0; row_index < H.NumRows(); row_index++) {
		std::bitset<(info_r + red_r)*col_r> row;
		for(int col_index = 0; col_index < H.NumCols(); col_index++) {
			row[col_index] = (H[row_index][col_index]==1); // if the value is 1, it returns true, if the value is 0, it returns false
		}
		bio_H.push_back(row);
	}
	bin_H_out << bio_H;
	bio_H.clear();
	bin_H_out.close();

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
		H_fr_2[row_index] = H_fr[row_index];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index][col_index] = H_fr[row_index][col_index];
		// }
	}
	// skip row 292, copy 292 rows
	offset_H_fr = col_r;
	int offset_H_fr_2 = offset_H_fr - 1; // 1 row was deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		H_fr_2[row_index + offset_H_fr_2] = H_fr[row_index + offset_H_fr];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		// }
	}
	// skip row 585, copy 292 rows
	offset_H_fr = 2*col_r;
	offset_H_fr_2 = offset_H_fr - 2; // 2 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		H_fr_2[row_index + offset_H_fr_2] = H_fr[row_index + offset_H_fr];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		// }
	}
	// skip row 878, copy 292 rows
	offset_H_fr = 3*col_r;
	offset_H_fr_2 = offset_H_fr - 3; // 3 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		H_fr_2[row_index + offset_H_fr_2] = H_fr[row_index + offset_H_fr];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		// }
	}
	// skip row 1171, copy 292 rows
	offset_H_fr = 4*col_r;
	offset_H_fr_2 = offset_H_fr - 4; // 4 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		H_fr_2[row_index + offset_H_fr_2] = H_fr[row_index + offset_H_fr];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		// }
	}
	// skip row 1464, copy 292 rows
	offset_H_fr = 5*col_r;
	offset_H_fr_2 = offset_H_fr - 5; // 5 rows were deleted
	for (int row_index = 0; row_index < col_r - 1; row_index++) {
		H_fr_2[row_index + offset_H_fr_2] = H_fr[row_index + offset_H_fr];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		// }
	}
	// skip row 1757, copy 293 rows
	offset_H_fr = 6*col_r;
	offset_H_fr_2 = offset_H_fr - 6; // 5 rows were deleted
	for (int row_index = 0; row_index < col_r; row_index++) {
		H_fr_2[row_index + offset_H_fr_2] = H_fr[row_index + offset_H_fr];
		// for (int col_index = 0; col_index < H_fr_2.NumCols(); col_index++) {
		// 	H_fr_2[row_index + offset_H_fr_2][col_index] = H_fr[row_index + offset_H_fr][col_index];
		// }
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
    			H_eye[row_index] += H_eye[pivot_index];
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

	// transform K into a vector<bitset> and (later on..) store it to a file. Then it will be possible to read it, & each row 
	// with the info word and sum modulo 2
	if (K.NumCols() != (int)ALL_INFO_BIT || K.NumRows() != (int)IND_EQ) {std::cout << "error, K has not the expected size"; return 1;}
	std::array< std::bitset<ALL_INFO_BIT>, IND_EQ > K_rows; // 2045

	std::cout << "Size of K_rows " << sizeof(K_rows) << "\n";
	for (int row_index = 0; row_index < K.NumRows(); row_index++) {
		std::bitset<ALL_INFO_BIT> row;
		for(int col_index = 0; col_index < K.NumCols(); col_index++) {
			row[col_index] = (K[row_index][col_index]==1); // if the value is 1, it returns true, if the value is 0, it returns false
		}
		K_rows[row_index] = row;
	}

	std::ofstream bin_out("K.bin", std::ios::out | std::ios::binary);
	BitIo<ALL_INFO_BIT> bio;
	for(int row_index = 0; row_index < K.NumRows(); row_index++) {
		bio.push_back(K_rows[row_index]);
	}
	bin_out << bio;
	bio.clear();
	bin_out.close(); 

	// once K is stored to file, check if it can be read correctly
	bool correct = 1;
	std::ifstream bin_in("K.bin", std::ios::binary);
	BitIo<ALL_INFO_BIT> bio2;
	bin_in >> bio2;
	for(int row_index = 0; row_index < K.NumRows(); row_index++) {
		std::bitset<ALL_INFO_BIT> row = bio2.pop_front();
		for (int bit_index = 0; bit_index < (int)ALL_INFO_BIT; bit_index++) {
			if(K_rows[row_index][bit_index] != row[bit_index]){
				std::cout << "Row " << row_index << " row bit " << K_rows[row_index][bit_index]
																	<< " bitset bit " << row[bit_index] << "\n";
				correct = 0;}
		}
	}
	std::cout << "Matrix K " << ((correct) ? "is":"is not") << " read correctly!\n";

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

	std::cout << "H_prime created\n";

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

	std::cout << "Perform some test: create a random vector, encode it, test on the parity check matrix and on the matrix defined in the std\n";
	for(int attempt = 0; attempt < 10; attempt ++) {

		std::cout << "attempt " << attempt << "\n";

		// create a random uncoded word
		mat_GF2 info_word;
		mat_GF2 smart_info_word;
		info_word.SetDims(info_bit, 1); 
		smart_info_word.SetDims(info_bit + init_zero_bit, 1);
		clear(smart_info_word); // all 0 for sure
		// create a bitset word too, in order to test encoding with bitset
		std::bitset<ALL_INFO_BIT> info_word_bitset;
		for(int bit_index = 0; bit_index < info_bit; bit_index++) { 
			info_word[bit_index][0] = int_uni_gen(m_rng);
			smart_info_word[mapJtoK(bit_index)][0] = info_word[bit_index][0];
			info_word_bitset[mapJtoK(bit_index)] = (info_word[bit_index][0]==1);
		}

		// fill a matrix as specified in the standard
		mat_GF2 std_M;
		std_M.SetDims(info_r+red_r, col_r);
		clear(std_M); // all 0 for sure
		int q = 0;
		int r = 0;
		for(int bit_index = 0; bit_index < info_bit; bit_index++) {
			q = bit_index + init_zero_bit; // j + 172, with j from 1 to 30592
			r = q/col_r; // floor(q/293) since both q and 293 are positive int, the result is the floor
			int col_index = col_r*r + col_r - 1 - q; // 293*r + 292 - q
			if(r == 0 && col_index == 120) {std::cout << "q = " << q << " r = " << r << " bit_index " << bit_index << "\n";}
			std_M[r][col_index] = info_word[bit_index][0];
		}

		// read the matrix by row, and encode it
		mat_GF2 inter_info_word;
		inter_info_word.SetDims(info_bit + init_zero_bit, 1);
		int read_bit_index = 0;
		for(int row_index = 0; row_index < info_r; row_index++) {
			for(int col_index = 0; col_index < std_M.NumCols(); col_index++) {
				inter_info_word[read_bit_index++][0] = std_M[row_index][col_index];
			}
		} 
		// check if bit (292-172, 292) are 0
		for(int bit_index = 292 - 172; bit_index < 292; bit_index++) {
			if(inter_info_word[bit_index][0] != 0) {std::cout << "element " << bit_index << " is not 0\n";}
		}

		// check if inter_info_word and smart_info_word are the same
		std::cout << "inter_info_word == smart_info_word ? " << (IsZero(inter_info_word + smart_info_word) ? "yes\n":"no\n");

		// encode it
		mat_GF2 code_word = G_prime*inter_info_word;

		// encode using bitset
		std::bitset<IND_EQ> parity_check_bitset;
		for(int ar_index = 0; ar_index < (int)IND_EQ; ar_index++) {
			parity_check_bitset[ar_index] = ((K_rows[ar_index]&info_word_bitset).count() & 1);
		}

		// check if this bitset is equal to the last 2045 bit of the codeword in mat_GF2
		offset = (int)ALL_INFO_BIT;
		for(int bit_index = offset; bit_index < code_word.NumRows(); bit_index++) {
			long gf2value = rep(code_word[bit_index][0]);
			long bitvalue = parity_check_bitset[bit_index-offset];
			if(gf2value != bitvalue) {std::cout << "bit " << bit_index << " is gf2value=" << gf2value << " bitvalue=" << bitvalue << "\n";}
		}

		// check it with H and H_prime
		std::cout << "H_prime x c = 0? " << (IsZero(H_prime*code_word) ? "yes\n":"no\n");
		std::cout << "H_fr x c = 0? " << (IsZero(H_fr*code_word) ? "yes\n":"no\n");

		// put the redundancy bit into std_M
		read_bit_index = info_r*col_r;
		for(int row_index = info_r; row_index < info_r + red_r - 1; row_index++) { // 292 bit in each row from 105 to 110
			for(int col_index = 0; col_index < col_r - 1; col_index++) {
				std_M[row_index][col_index] = code_word[read_bit_index++][0];
			}
		}
		// 293 bit in row 111
		for(int col_index = 0; col_index < col_r; col_index++) {
			std_M[std_M.NumRows()-1][col_index] = code_word[read_bit_index++][0];
		}

		// check if the 2051 lines sum to 0
		for(int slope_index = 0; slope_index < red_r; slope_index++) {
			// cycle on the offsets
			for(int c = 0; c < col_r; c++) {
				// cycle on a from 0 to 111
				int lineValue = 0;
				for(int a = 0; a < (info_r + red_r); a++) {
					int col_index = (a*slopes[slope_index] + c)%col_r;
					if(std_M[a][col_index] == 1) {lineValue++;}
				}
				if(lineValue%2 != 0) {std::cout << "s " << slopes[slope_index] << " c " << c << "\n";}
			}
		}

		// now read the complete codeword (with also the always 0 redundancy bit) and multiply it by the original H
		// read the matrix by row, and encode it
		mat_GF2 complete_code_word;
		complete_code_word.SetDims(std_M.NumRows()*std_M.NumCols(), 1);
		int write_bit_index = 0;
		for(int row_index = 0; row_index < info_r + red_r; row_index++) {
			for(int col_index = 0; col_index < std_M.NumCols(); col_index++) {
				complete_code_word[write_bit_index++][0] = std_M[row_index][col_index];
			}
		} 

		std::cout << "H x c (complete) = 0? " << (IsZero(H*complete_code_word) ? "yes\n":"no\n");

	}

	







	return 0;
}


