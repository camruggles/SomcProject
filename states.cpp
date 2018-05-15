
#include <iostream>
#define n 298182
#define k 1652

int main(){
	
	int count = 0;	
	int *arr = new int[k];
	for (int i = 0; i < k; i++){
		arr[i] = 0;
	}

	int a, b, c;

	for (int i = 0; i < n; i++){
		std::cin >> a >> b >> c;
		arr[a] = arr[b] = arr[c] = 1;
	}

	for (int i = 0; i < k; i++){
		if (arr[i] != 1){
			std::cout << i << "\n";
			count++;
		}	
	}

	std::cout << "Total: " << count << "\n";

}
