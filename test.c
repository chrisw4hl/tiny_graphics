#include <stdio.h>

const int space_size=100;

void get_intensities(size_t r, size_t c, int cube[][100]){

	for(int j=0; j < space_size; j++){
		for(int i=0; i < space_size; i++){
			printf("%d", cube[j][i]);
		}
		printf("\n");
	}
}

void get_normals(int shape[][100]){

	const int view_angle = 45; //set camera perspective
}

int main(){

	printf("hello world\n");
	
	int cube[100][100] = {0};
	
//	int obj[][] = get_shape();
	
	size_t total = sizeof(cube)/sizeof(cube[0][0]);
	size_t c = sizeof(cube[0])/sizeof(cube[0][0]);
	size_t r = sizeof(cube)/sizeof(cube[0]);
	printf("%zu,%zu,%zu,%d\n",r,c,total,cube[r-1][c-1]);
	get_intensities(r, c, cube);

//	project_view(view);



	return 0;

}
