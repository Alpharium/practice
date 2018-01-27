//will code for exp
/*Кубик, грани которого помечены цифрами от 1 до 6, бросают N раз. Найти вероятность того, что сумма выпавших чисел будет равна Q.
Ограничения: 1 <= N <= 500, 1 <= Q <= 3000.
Входные данные
В первой строке находятся числа N и Q через пробел.
Выходные данные
Вероятность того, что сумма выпавших чисел будет равна Q.

сумма вероятностей выпадения значений на кубике умноженых на вероятность того что на предыдущем шаге значение суммы отличалось именно на выпавшее значение на кубике на текущем шаге

*/
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iomanip>


double funct(int n, int q) {
	std::vector<double> mass(q + 1, 0);
	mass[0] = 1;
	std::vector<double> temp(q + 1, 0);
	if (q < n) return 0;
	if ((6 * n) < q) return 0;
	for (int x = 1; x <= n;x++) {
		mass.swap(temp);
		for (int y = 0;y <= q;y++) {
			mass[y] = 0;
			for (int value = 1;value <= 6;value++) {
				if (y - value >= 0) {
					mass[y] = temp[y - value] + mass[y];
				}
			}
			mass[y] /= 6;
		}
	}
	return mass[q];
};
int main()
{
	std::cout << "\n 1<=N<=500; 1<=Q<=3000\n\n\n";
	int a, b;
	std::cin >> a >> b;
	std::cout << std::setprecision(40) << std::scientific << funct(a, b);
	
	
	
	getchar();
	getchar();
	return 0;
}
