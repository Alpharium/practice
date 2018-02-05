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
	//double a, b,z,sx2,x_avg,sum,sum2,sum3,sum4; //a-coefficient,z-
	//x_avg=sum= sum2=sum3=sum4= z = 0.0;
	//for (int i = 0;i < n_calc;i++) {
	//	x_avg += mas_of_points[i].first;
	//}
	//x_avg = x_avg / n_calc;
	//for (int i = 0;i < n_calc;i++) {
	//	sum += mas_of_points[i].second*(mas_of_points[i].first - x_avg);//sum(yi(xi-x_avg)
	//	sum2 += pow((mas_of_points[i].first - x_avg), 2);//sum((xi-x_avg)^2)
	//	sum3 += pow(mas_of_points[i].first, 2);//sum((xi)^2)
	//}
	//sx2 = (1.0 / n_calc)*sum2;
	//z = (1.0 / n_calc)*sum3;
	//a = (1.0 / (n_calc*sx2))*sum;
	//for (int i = 0;i < n_calc;i++) {
	//	sum4 += mas_of_points[i].second*(z - mas_of_points[i].first*x_avg);//sum (yi(z-xi*x_avg))
	//}
	//b = (1.0 / (n_calc*sx2))*sum4;


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

	//-------------------------------------------------------- cube extrapolation-------------
	double q, w, e, r, t, u, p, g, h, m, a3, b3, c3, d3;
	q = w = e = r = t = u = p = g = h = m = a3 = b3 = c3 = d3=0.0;
		for (int i = 0;i < n_calc;i++) {
		q += (pow(mas_of_points[i].first, 3)*mas_of_points[i].second);
		w += (pow(mas_of_points[i].first, 2)*mas_of_points[i].second);
		e += mas_of_points[i].first*mas_of_points[i].second;
		r += mas_of_points[i].second;//E(x*y)
		t += pow(mas_of_points[i].first, 6);
		u += pow(mas_of_points[i].first, 5);
		p += pow(mas_of_points[i].first, 4);
		g += pow(mas_of_points[i].first, 3);
		h += pow(mas_of_points[i].first, 2);
		m += mas_of_points[i].first;
	}
		a3 = ((-1.0*pow(h,3)*q-1.0*pow(g,2)*n_calc*q-1.0*pow(m,2)*p*q+pow(g,3)*r+m*pow(p,2)*r-1.0*pow(g,2)*m*u+g*n_calc*p*u-1.0*g*m*r*u+pow(h,2)*(g+r)*u+pow(m,2)*pow(u,2)+e*(-1.0*pow(g,2)*h+pow(h,2)*p+g*m*p-1.0*n_calc*pow(p,2)-1.0*h*m*u+g*n_calc*u)+h*(2.0*g*m*q+n_calc*p*q-2.0*g*p*r-1.0*m*p*u-1.0*n_calc*pow(u,2))) / (pow(g,4)+pow(h,2)*pow(p,2)-1.0*n_calc*pow(p,3)-1.0*pow(h,3)*t-1.0*pow(m,2)*p*t+pow(m,2)*pow(u,2)-1.0*pow(g,2)*(3.0*h*p+n_calc*t+2.0*m*u)+2.0*g*(m*(pow(p,2)+h*t)+(pow(h,2)+n_calc*p)*u)+h*(n_calc*p*t-2.0*m*p*u-1.0*n_calc*pow(u,2))));////////////////////////////////////
		b3 = -1.0*((h*m*p*q-1.0*h*pow(p,2)*r+pow(h,2)*r*t+n_calc*pow(p,2)*u-1.0*pow(m,2)*q*u+h*n_calc*q*u+m*p*r*u+pow(m,2)*t*u-1.0*h*n_calc*t*u+pow(g,2)*(m*q+p*r+h*u)-1.0*g*(pow(h,2)*q+n_calc*p*q+m*r*t+2.0*m*p*u+h*r*u)-1.0*e*(pow(g,3)+h*m*t+n_calc*p*u-1.0*g*(h*p+n_calc*t+m*u))) / (pow(g,4)+pow(h,2)*pow(p,2)-1.0*n_calc*pow(p,3)-1.0*pow(h,3)*t-1.0*pow(m,2)*p*t+pow(m,2)*pow(u,2)-1.0*pow(g,2)*(3.0*h*p+n_calc*t+2.0*m*u)+2.0*g*(m*(pow(p,2)+h*t)+(pow(h,2)+n_calc*p)*u)+h*(n_calc*p*t-2.0*m*p*u-1.0*n_calc*pow(u,2))));
		c3 = ((pow(h,2)*p*q-1.0*n_calc*pow(p,2)*q-1.0*m*r*t+pow(g,3)*u-1.0*h*m*q*u-1.0*h*p*r*u+h*m*t*u+n_calc*p*pow(u,2)+m*r*pow(u,2)-1.0*pow(g,2)*(h*q+r*u)+g*(m*p*q+pow(p,2)*r+h*r*t-1.0*h*p*u+n_calc*q*u-1.0*n_calc*t*u-1.0*m*pow(u,2))-1.0*e*(pow(g,2)*p+pow(h,2)*t-1.0*n_calc*p*t-2.0*g*h*u+n_calc*pow(u,2))) / (pow(g,4)+pow(h,2)*pow(p,2)-1.0*n_calc*pow(p,3)-1.0*pow(h,3)*t-1.0*pow(m,2)*p*t+pow(m,2)*pow(u,2)-1.0*pow(g,2)*(3.0*h*p+n_calc*t+2.0*m*u)+2.0*g*(m*(pow(p,2)+h*t)+(pow(h,2)+n_calc*p)*u)+h*(n_calc*p*t-2.0*m*p*u-1.0*n_calc*pow(u,2))));
		d3 = ((pow(g,3)*q+m*pow(p,2)*q-1.0*pow(p,3)*r+h*p*r*t+h*pow(p,2)*u+pow(h,2)*q*u-1.0*pow(h,2)*t*u-1.0*m*p*pow(u,2)-1.0*h*r*pow(u,2)-1.0*pow(g,2)*(r*t+p*u)+g*(2.0*p*r*u+m*(-1.0*q+t)*u+h*(-2.0*p*q+pow(u,2)))+e*(g*(pow(p,2)+h*t)-1.0*pow(g,2)*u-1.0*h*p*u+m*(-1.0*p*t+pow(u,2)))) / (pow(g,4)+pow(h,2)*pow(p,2)-1.0*n_calc*pow(p,3)-1.0*pow(h,3)*t-1.0*pow(m,2)*p*t+pow(m,2)*pow(u,2)-1.0*pow(g,2)*(3.0*h*p+n_calc*t+2.0*m*u)+2.0*g*(m*(pow(p,2)+h*t)+(pow(h,2)+n_calc*p)*u)+h*(n_calc*p*t-2.0*m*p*u-n_calc*pow(u,2))));


	system("pause");
	return 0;
}
