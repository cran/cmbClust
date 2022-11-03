#include <stdio.h>
#include <stdlib.h>

#include "array.h"
#include "matrix.h"





void AllPerms(int (*size1), int *perms1){

    int **perms, size;
    size = (*size1);
    int length = Factorial(size);

    MAKE_2ARRAY(perms, length, size);
    array1to2i(length,size,perms1, perms);

	int sch, i, j, v, w, finish, flag, ind;
	double **pat;
	int *cn;

	sch = 0;
	i = 0;
	j = -1;
	flag = 0;
	finish = 0;
	ind = 0;

	MAKE_MATRIX(pat, size, size);
	for (v=0; v<size; v++){
		for (w=0; w<size; w++){
			pat[v][w] = 0;
		}
	}

	MAKE_VECTOR(cn, size);
	for (v=0; v<size; v++){
		cn[v] = 0;
	}


	while (finish == 0){

		if (j != (size-1)){
			j = j+1;
		} else {
			if (flag == 1){
				j = 0;
				i = i+1;
				flag = 0;
			}
		}

		if (pat[i][j] == 0){
			for (v=0; v<size; v++){
				pat[i][v]=1;
				pat[v][j]=1;
			}

			sch = sch + 1;
			cn[sch-1] = j;
			flag = 1;
		}

		if ((sch == size) & (flag == 1)){

			for (v=0; v<size; v++){
				perms[ind][v] = cn[v];
			}

			ind++;
			flag = 0;
			sch = sch - 1;
			i = i - 1;
			j = cn[sch-1];
			sch = sch-1;

			for (v=0; v<size; v++){
				for (w=0; w<size; w++){
					pat[v][w] = 0;
				}
			}

			for (v=0; v<sch; v++){
				for (w=0; w<size; w++){
					pat[v][w] = 1;
					pat[w][cn[v]] = 1;
				}
			}

		}


		if ((j == size - 1) & (flag == 0)){
			i = i - 1;

			sch = sch-1;

			if (sch >= 0){

				j = cn[sch];

				for (v=0; v<size; v++){
					for (w=0; w<size; w++){
						pat[v][w] = 0;
					}
				}

				if (sch > 0){
					for (v=0; v<sch; v++){
						for (w=0; w<size; w++){
							pat[v][w] = 1;
							pat[w][cn[v]] = 1;
						}
					}

				}

			}

			if (i >= 0){
				pat[i][j] = 1;
			}
		}

		if (sch == -1){
			finish = 1;
		}

	}

	FREE_MATRIX(pat);
	FREE_VECTOR(cn);



    array2to1i(length, size, perms1, perms);
    FREE_2ARRAY(perms);
}
