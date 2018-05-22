
#include <iostream>
#include <stdio.h>
#include <vector>
#define k 8808

int main(){
	

	if (k == 0){
		printf("Change the macros for the state numbers\n");
		return -1;
	}

	//intake variables
	int a, b, c;
	std::vector<int> vec;
	for (int i = 0; i < k; i++){
		std::cin >> a >> b >> c;
		vec.push_back(a);
		vec.push_back(b);
		vec.push_back(c);
	}
	
	//find max
	int max = -1;
	for (int i : vec){
		if (i > max)
		max = i;
	}
	
	//count nonrecurrent states
	int * arr = new int[max+1];
	for (int i = 0; i < max+1; i++){
		arr[i] = 0;
	}

	for (int i : vec){
		arr[i] = 1;
	}

	int num = 0;
	for (int i = 0; i < max+1; i++){
		if (arr[i] == 0){
			num++;
			printf("%d\n", i);
		}
	}
	printf("Num: %d\n", num);
	printf("Max: %d\n", max);


}
