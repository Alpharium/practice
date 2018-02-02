//this code isnt finished, so be patient please
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

using namespace std;
int main()
{
	vector<pair<double, double>> mas_of_points;
	vector<pair<double, double>> mas_of_arguments;
	ofstream fout;
	ifstream fin;
	int n_calc, k_arg;
	int n_mas[N];
	int k_mas[K];
	std::string char_mas,char_tmp;/*
	fout.open("solution.txt");
	fout << "4\nHello World\n\KAPPA PRIDE\n";
	fout.close();*/
	fin.open("calculated.txt");
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
		//fin >> char_tmp;
		mas_of_arguments.push_back(std::make_pair(atof(char_mas.c_str()),0));
	}
	fin.close();

	//cout << char_mas << endl;
	/*for (int i = 0;char_mas != '\0';i++) {
		cout << char_mas[i];
	}*/
	for (int i = 0;i < n_calc;i++) {
		cout << mas_of_points[i].first << " " << mas_of_points[i].second << "\n";
	}

	system("pause");
/*
	getchar();
	getchar();*/
	return 0;
}
