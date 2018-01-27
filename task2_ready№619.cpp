//will code for food
/*�����, ����� �������� �������� ������� �� 1 �� 6, ������� N ���. ����� ����������� ����, ��� ����� �������� ����� ����� ����� Q.
�����������: 1 <= N <= 500, 1 <= Q <= 3000.
������� ������
� ������ ������ ��������� ����� N � Q ����� ������.
�������� ������
����������� ����, ��� ����� �������� ����� ����� ����� Q.

����� ������������ ��������� �������� �� ������ ��������� �� ����������� ���� ��� �� ���������� ���� �������� ����� ���������� ������ �� �������� �������� �� ������ �� ������� ����

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