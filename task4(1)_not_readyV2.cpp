#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
const int N = 1000;
const int K = 1000;
using std::vector;
using std::pair;
using std::ofstream;
using std::ifstream;

class all_reader
{
	
private:
	vector<pair<double, double>> mas_of_points;
	vector<double> mas_of_arguments;
	ofstream fout;
	ifstream fin;
	int n_calc;
	int k_arg;
	std::string char_mas;
	std::string	char_tmp;
	std::string file_name;

public:
	vector<pair<double,double>> fftovpoints() {
		std::cin >> file_name;
		fin.open(file_name);
		fin >> n_calc;
		for (int i = 0;i < n_calc;i++) {
			fin >> char_mas;
			fin >> char_tmp;
			mas_of_points.push_back(std::make_pair(atof(char_mas.c_str()), atof(char_tmp.c_str())));
		}
		fin.close();
		return mas_of_points;
	}
	vector<double> fftovarguments() {
		std::cin >> file_name;
		fin.open(file_name);
		fin >> k_arg;
		for (int i = 0;i < k_arg;i++) {
			fin >> char_mas;
			mas_of_arguments.push_back(atof(char_mas.c_str()));
		}
		fin.close();
		return mas_of_arguments;
	}
};

int main()
{
	vector<pair<double, double>> mas_of_points;
	vector<pair<double, double>> mas_of_arguments;
	ofstream fout;
	ifstream fin;
	int n_calc, k_arg;
	int n_mas[N];
	int k_mas[K];
	std::string char_mas,char_tmp;
	std::cin >> char_mas;
	fin.open(char_mas);
	fin >> n_calc;
	for (int i = 0;i < n_calc;i++){
		fin >> char_mas;
		fin >> char_tmp;
		mas_of_points.push_back(std::make_pair(atof(char_mas.c_str()), atof(char_tmp.c_str())));
	}
	fin.close();
	fin.open("arguments.txt");
	fin>>k_arg;
	for (int i = 0;i < k_arg;i++) {
		fin >> char_mas;
		mas_of_arguments.push_back(std::make_pair(atof(char_mas.c_str()),0));
	}
	fin.close();
	fout.open("outputed.txt");
	for (int i = 0;i < n_calc;i++) {
		fout << mas_of_points[i].first << " " << mas_of_points[i].second << "\n";
	}
	fout.close();
	//-------------------------------------------------------- linear extrapolation-------------
	double a, b,z,sx2,x_avg,sum,sum2,sum3,sum4; //a-coefficient,z-
	x_avg=sum= sum2=sum3=sum4= z = 0.0;
	for (int i = 0;i < n_calc;i++) {
		x_avg += mas_of_points[i].first;
	}
	x_avg = x_avg / n_calc;
	for (int i = 0;i < n_calc;i++) {
		sum += mas_of_points[i].second*(mas_of_points[i].first - x_avg);//sum(yi(xi-x_avg)
		sum2 += pow((mas_of_points[i].first - x_avg), 2);//sum((xi-x_avg)^2)
		sum3 += pow(mas_of_points[i].first, 2);//sum((xi)^2)
	}
	sx2 = (1.0 / n_calc)*sum2;
	z = (1.0 / n_calc)*sum3;
	a = (1.0 / (n_calc*sx2))*sum;
	for (int i = 0;i < n_calc;i++) {
		sum4 += mas_of_points[i].second*(z - mas_of_points[i].first*x_avg);//sum (yi(z-xi*x_avg))
	}
	b = (1.0 / (n_calc*sx2))*sum4;


	//-------------------------------------------------------- square extrapolation-------------
	double q, w, e, r, t, u, l, a2, b2, c;// in comments "E"-means summ from i=1 to N (where N=n_calc)
	//Q=Ex; W=Ex^2; E=Ex^3; R=Ex^4;T=E((x^2)*y); U=E(x*y); L=Ey;n=n_calc
	q = w = e = r = t = u = l = a2 = b2 = c = 0.0;
	for (int i = 0;i < n_calc;i++) {
		q += mas_of_points[i].first;
		w += pow(mas_of_points[i].first, 2);
		e += pow(mas_of_points[i].first, 3);
		r += pow(mas_of_points[i].first, 4);
		t += (pow(mas_of_points[i].first, 2))*mas_of_points[i].second;
		u += mas_of_points[i].first*mas_of_points[i].second;
		l += mas_of_points[i].second;
	}
	a2 = ((-1.0 * e*l*q + pow(q, 2)*t+e*n_calc*u-1.0*n_calc*t*w-1.0*q*u*w+l*pow(w,2)) / (pow(e,2)*n_calc+pow(q,2)-2.0*e*q*w-1.0*n_calc*r*w+pow(w,3)));
	b2 = ((l*q*r+e*n_calc*t-1.0*n_calc*r*u-1.0*e*l*w-1.0*q*t*w+u*pow(w,2)) / (pow(e,2)*n_calc+pow(q,2)*r-2.0*e*q*w-1.0*n_calc*r*w+pow(w,3)));
	c = ((pow(e,2)*l+q*r*u+w*(-1.0*l*r+t*w)-1.0*e*(q*t+u*w)) / (pow(e,2)*n_calc+pow(q,2)*r-2.0*e*q*w-1.0*n_calc*r*w+pow(w,3)));

	system("pause");
	return 0;
}
