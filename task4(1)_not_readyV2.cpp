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
	std::string char_mas,char_tmp;/*
	fout.open("solution.txt");
	fout << "4\nHello World\n\KAPPA PRIDE\n";
	fout.close();*/
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
		//fin >> char_tmp;
		mas_of_arguments.push_back(std::make_pair(atof(char_mas.c_str()),0));
	}
	fin.close();
/*
	cout << char_mas << endl;*/
	/*for (int i = 0;char_mas != '\0';i++) {
		cout << char_mas[i];
	}*/
	fout.open("outputed.txt");
	for (int i = 0;i < n_calc;i++) {
		fout << mas_of_points[i].first << " " << mas_of_points[i].second << "\n";
	}
	fout.open("outputed.txt");
	system("pause");
/*
	getchar();
	getchar();*/
	return 0;
}
