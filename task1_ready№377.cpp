/*На числовой прямой окрасили N отрезков. Известны координаты левого и правого концов каждого отрезка (Li и Ri). Найти длину окрашенной части числовой прямой.
Входные данные
В первой строке находится число N, в следующих N строках - пары Li и Ri. Li и Ri - целые, -1 000 000 000 <= Li <= Ri <= 1 000 000 000, 1 <= N <= 15 000 
Выходные данные
Вывести одно число - длину окрашенной части прямой. 



объединяю перекрывающиеся отрезки,
затем считаю сумму всех длин неперекрывающихся отрезков
*/
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
bool pairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>&secondElem) {
	return firstElem.first < secondElem.first;
}
int main()
{
	std::vector <std::pair <int, int>> mass;
	std::vector <int> lenghts;
	int n_mas;
	std::cin >> n_mas;
	mass.reserve(n_mas);
	for (int i = 0;i < n_mas;i++) {
		int left, right;
		std::cin >> left >> right;
		mass.push_back({ left,right });
	}
	int left, right;
	std::sort(mass.begin(), mass.end(), pairCompare);
	int lenght = 0;
	int cycler = 0;
	for (int n = 0;n < n_mas ;) {
		left = mass[n].first;
		right = mass[n].second;
		for (int k = n;k < n_mas;k++) {
			if (mass[k].first > right) {
				break;
			}
			cycler++;
			if (mass[k].first<=right) {
				if (mass[k].second > right) {
					right = mass[k].second;
				}
			}

		}
		lenght += right - left;
		lenghts.push_back(lenght);
		lenght = 0;
		n = n + cycler;
		cycler = 0;

	}
	int result = 0;

	for (int i = 0;i < lenghts.size();i++) {
		result += lenghts[i];
	}
	std::cout << "\n\n" << result << "\n";
	getchar();
	getchar();
	return 0;
}
