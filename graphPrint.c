#include <stdio.h>

int main(int argc, char ** argv){
	FILE * fp = fopen(argv[1], "r");
	FILE * fout = fopen(argv[2], "w");
	char p = 'a';
	while ((p = fgetc(fp)) != EOF){

		if (p == 'P'){
			fputc('\n', fout);
		}

		fputc(p, fout);

	}
	fclose(fp);
	fclose(fout);

}
