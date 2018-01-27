
/*Íà ÷èñëîâîé ïðÿìîé îêðàñèëè N îòðåçêîâ. Èçâåñòíû êîîðäèíàòû ëåâîãî è ïðàâîãî êîíöîâ êàæäîãî îòðåçêà (Li è Ri). Íàéòè äëèíó îêðàøåííîé ÷àñòè ÷èñëîâîé ïðÿìîé.
Âõîäíûå äàííûå
Â ïåðâîé ñòðîêå íàõîäèòñÿ ÷èñëî N, â ñëåäóþùèõ N ñòðîêàõ - ïàðû Li è Ri. Li è Ri - öåëûå, -1 000 000 000 <= Li <= Ri <= 1 000 000 000, 1 <= N <= 15 000 
Âûõîäíûå äàííûå
Âûâåñòè îäíî ÷èñëî - äëèíó îêðàøåííîé ÷àñòè ïðÿìîé. 



îáúåäèíÿþ ïåðåêðûâàþùèåñÿ îòðåçêè,
çàòåì ñ÷èòàþ ñóììó âñåõ äëèí íåïåðåêðûâàþùèõñÿ îòðåçêîâ
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
