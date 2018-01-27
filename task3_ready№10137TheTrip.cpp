//please give me a job
/* The Trip - "Ïóòåøåñòâèÿ" Input
Standard input will contain the information for several trips. The information for each trip consists of
a line containing a positive integer, n, the number of students on the trip, followed by n lines of input,
each containing the amount, in dollars and cents, spent by a student. There are no more than 1000
students and no student spent more than $10,000.00. A single line containing 0 follows the information
for the last trip.
Output
For each trip, output a line stating the total amount of money, in dollars and cents, that must be
exchanged to equalize the students’ costs.
*/


#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iomanip>
int main()
{
	int n;//students
	double costs;
	std::vector <double> dcosts;
	std::vector<double> exchange;
	std::vector<int> students;
	while (std::cin >> n) {
		if (n == 0) {
			std::cout << "\n you quit the program\n";
			break;
		};
		students.push_back(n);
		if (n > 0) {
			for (int i = 0;i < n;i++) {
				double in;
				std::cin >> in;
					dcosts.push_back(in);
			}
		}
	}
	int step=0;
	for (int i = 0;i < students.size();i++) {
		double totalsum = 0;
		double avg = 0;
		for (int k = 0;k < students[i];k++) {
			totalsum += dcosts[step+k];
		}
		avg = totalsum / students[i];
		double tmp=0.0;
		for (int k = 0;k < students[i];k++) {
			if (dcosts[step+k] > avg) {
				double sub=((long)((dcosts[step+k] - avg)*100.0) / 100.0);
				tmp += sub;
			}
		}
		step += students[i];
		exchange.push_back(tmp);
	}
	for (int i = 0;i < exchange.size();i++) {
		std::cout << "\n" << exchange[i];

	}
	getchar();
	getchar();
	return 0;
}
